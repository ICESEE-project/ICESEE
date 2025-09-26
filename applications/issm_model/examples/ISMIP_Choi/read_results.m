% -----------------------------------------------------------
% @author: 	Brian Kyanjo
% @date: 		2025-04-30
% @brief: 	Reads and plot results from both ISSM and ICESEE
% ------------------------------------------------------------

%% Make the $ISSM_DIR environment variable available
% issm_dir = getenv('ISSM_DIR');  % Retrieve the ISSM_DIR environment variable
% % addpath(genpath(issm_dir));     % Add the ISSM directory and its subdirectories to the MATLAB path
% 
% % path to the results
% results_dir = fullfile(issm_dir, 'examples', 'ISMIP_Choi', 'Models','ens_id_0');  % Path to the results directory
% forecast_dir = fullfile(issm_dir, 'examples', 'ISMIP_Choi', 'Models','ens_id_0')
% 
% %% plot surface velocities
% % Load the ISSM results
% md_true = loadmodel(fullfile(results_dir, 'true_state.mat'));  % Load the ISSM results from a .mat file
% md_nurged = loadmodel(fullfile(results_dir, 'enkf_state.mat'));  % Load the ISSM results from a .mat file 
% plotmodel(md_true, 'data', md_true.results.TransientSolution.Vel, 'layer', 5, 'figure', 5);
% plotmodel(md_nurged, 'data', md_nurged.results.TransientSolution.Vel, 'layer', 5, 'figure', 6);

%% ICESEE results
% Get the Python version
% pyversion = py.sys.version;
% 
% % Add the configuration directory to the Python path
% py.sys.path().append('../../config');
% 
% % Import the Python module _utility_imports
% utility_imports = py.importlib.import_module('_utility_imports');

% Load the essential data
results_dir = 'results';
filter_type = 'true-wrong';
file_path   = fullfile(results_dir, sprintf('%s-issm.h5', filter_type));
t           = h5read(file_path,'/t');
ind_m       = h5read(file_path,'/obs_index');
tm_m        = h5read(file_path,'/obs_max_time');
run_mode    = h5read(file_path,'/run_mode');

% load the true and nurged states
file_path            = 'data/new_data/_modelrun_datasets/true_nurged_states.h5';
model_true_state     = h5read(file_path,'/true_state')';
model_nurged_state   = h5read(file_path, '/nurged_state')';

% load observation data
file_path  = 'data/new_data/_modelrun_datasets/synthetic_obs.h5';
w          = h5read(file_path, '/hu_obs')'; 

% load the ensemble data
file_path         = 'data/new_data/_modelrun_datasets/icesee_ensemble_data.h5';
ensemble_vec_full = h5read(file_path, '/ensemble'); 
ensemble_vec_mean = h5read(file_path, '/ensemble_mean')';

% Process and plot
[ndim, nt] = size(model_true_state);
num_steps = nt - 1;
hdim = floor(ndim / 3);

file_path   = fullfile("data", "ISMIP_initial_data.mat");
md = loadmodel(file_path);
md_true = md; md_nurged = md;
md_mean = md; md_ens = md;

x = md.mesh.x;
y = md.mesh.y;
k =  nt-1;
True_fcoeff = model_true_state(2*hdim+1:3*hdim, k);
True_bed = model_true_state(hdim+1:2*hdim,  k);
True_thickness= model_true_state(1:hdim,  k);
md_true.geometry.bed=True_bed;
md_true.geometry.thickness=True_thickness;
md_true.friction.coefficient=True_fcoeff;
% figure;
% plotmodel(md_true, 'data', md_true.geometry.thickness, 'title', 'Ice Thickness'); hold off;
% figure; 
% plotmodel(md_true, 'data', md_true.geometry.bed, 'title', 'Bed Topography');
% plotmodel(md_true, 'data', md_true.friction.coefficient, 'title', 'Friction Coefficient');

file_path_true= fullfile("issm_data","true_state.mat");
file_path_nurged= fullfile("issm_data","nurged_state.mat");
md_true_= loadmodel(file_path_true);
md_nurged_= loadmodel(file_path_nurged);

% load grounding lines 

for i=1:500
    gl_data_true{i} = md_true_.results.TransientSolution(i).MaskOceanLevelset;
    gl_data_nurged{i} = md_nurged_.results.TransientSolution(i).MaskOceanLevelset;

    di =  md.materials.rho_ice / md.materials.rho_water;
    Hcritical = -di * model_true_state(hdim+1:2*hdim,  i);
    ice_thickness = model_true_state(1:hdim,  i);
    % get the position at which H<=Hcritical
    pos = find((ice_thickness <= Hcritical) & (ice_thickness > 0));
    % gl_data{i} = md_true_.mesh.x(pos)
    gl_data{i} = model_true_state(1:hdim,  i) + model_true_state(hdim+1:2*hdim,  i)/di;
    % plotmodel(md_true_, 'data', gl_data{i}, 'title', sprintf('True Grounding Line Position at time %d', i)); hold off;

    gl_data_ens{i}=ensemble_vec_mean(1:hdim, i) + ensemble_vec_mean(hdim+1:2*hdim, i)/di;
    plotmodel(md_true_, 'data', gl_data_ens{i}, 'title', sprintf('Ensemble Mean Grounding Line Position at time %d', i)); hold off;
    

    % plotmodel(md_nurged_, 'data', gl, 'title', sprintf('Nurged Ice Thickness + Bed at time %d', i)); hold off;
    % subplots(2,1,1);
    % plotmodel(md_nurged_, 'data', gl_data_nurged{i}, 'title', sprintf('Nurged Grounding Line Position at time %d', i)); hold off;
    % subplots(2,1,2);
    % plotmodel(md_true, 'data', gl_data_true{i}, 'title', sprintf('Grounding Line Position at time %d', i)); hold off;
end
% figure;
% plotmodel(md_true, 'data', gl_data_true{500}, 'title', 'Grounding Line Position at Final Time'); hold off;

