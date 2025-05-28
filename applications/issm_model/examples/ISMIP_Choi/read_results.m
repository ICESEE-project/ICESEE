% -----------------------------------------------------------
% @author: 	Brian Kyanjo
% @date: 		2025-04-30
% @brief: 	Reads and plot results from both ISSM and ICESEE
% ------------------------------------------------------------

%% Make the $ISSM_DIR environment variable available
issm_dir = getenv('ISSM_DIR');  % Retrieve the ISSM_DIR environment variable
% addpath(genpath(issm_dir));     % Add the ISSM directory and its subdirectories to the MATLAB path

% path to the results
results_dir = fullfile(issm_dir, 'examples', 'ISMIP_Choi', 'Models','ens_id_0');  % Path to the results directory
forecast_dir = fullfile(issm_dir, 'examples', 'ISMIP_Choi', 'Models','ens_id_0')

%% plot surface velocities
% Load the ISSM results
md_true = loadmodel(fullfile(results_dir, 'true_state.mat'));  % Load the ISSM results from a .mat file
md_nurged = loadmodel(fullfile(results_dir, 'enkf_state.mat'));  % Load the ISSM results from a .mat file 
plotmodel(md_true, 'data', md_true.results.TransientSolution.Vel, 'layer', 5, 'figure', 5);
plotmodel(md_nurged, 'data', md_nurged.results.TransientSolution.Vel, 'layer', 5, 'figure', 6);

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
file_path            = '_modelrun_datasets/true_nurged_states.h5';
model_true_state     = h5read(file_path,'/true_state')';
model_nurged_state   = h5read(file_path, '/nurged_state')';

% load observation data
file_path  = '_modelrun_datasets/synthetic_obs.h5';
w          = h5read(file_path, '/hu_obs')'; 

% load the ensemble data
file_path         = '_modelrun_datasets/icesee_ensemble_data.h5';
ensemble_vec_full = h5read(file_path, '/ensemble'); 
ensemble_vec_mean = h5read(file_path, '/ensemble_mean')';

% Process and plot
[ndim, nt] = size(model_true_state);
num_steps = nt - 1;
hdim = floor(ndim / 4);

% velocity componets
utrue  = model_true_state(1:hdim, :);
unurge = model_nurged_state(1:hdim, :);
vtrue  = model_true_state(hdim+1:2*hdim, :);
vnurge = model_nurged_state(hdim+1:2*hdim, :);
wtrue  = model_true_state(2*hdim+1:3*hdim, :);
wnurge = model_nurged_state(2*hdim+1:3*hdim, :);
% magnitude of the velocity
vtrue_magnitude  = sqrt(utrue.^2 + vtrue.^2 + wtrue.^2);
vnurge_magnitude = sqrt(unurge.^2 + vnurge.^2 + wnurge.^2);

% pressure
ptrue  = model_true_state(3*hdim+1:4*hdim, :);
pnurge = model_nurged_state(3*hdim+1:4*hdim, :);

% plot flag
plot_flag = 'pressure';
if strcmp(plot_flag, 'velocity')
    htrue = vtrue_magnitude;
    hnurge = vnurge_magnitude;
elseif strcmp(plot_flag, 'pressure')
    htrue = ptrue;
    hnurge = pnurge;
end

profile_flag = 'middle';
if strcmp(profile_flag, 'beginning')
    h_indx = 1;
elseif strcmp(profile_flag, 'middle')
    h_indx = floor(size(htrue, 2) / 2);
elseif strcmp(profile_flag, 'end')
    h_indx = size(htrue, 2);
elseif strcmp(profile_flag, 'random')
    h_indx = randi(size(htrue, 2));
else
    error('Invalid profile_flag. Use: beginning, middle, end, random');
end

fprintf('At h_indx = %d profile\n', h_indx);

h_true = htrue(h_indx, :);
h_nurged = hnurge(h_indx, :);

figure('Position', [100, 100, 800, 400]);
plot(t, h_true, 'r', 'DisplayName', 'True'); hold on;
plot(t, h_nurged, 'g', 'DisplayName', 'Wrong');
obs = w(h_indx, :);
fprintf('obs shape: [%d, %d]\n', size(obs));
plot(t(ind_m), obs, 'kx', 'DisplayName', 'Observations');
xlabel('Time (years)');
% ylabel('Pressure (Pa)');
legend('show');
title(sprintf('%s profile %s at h_indx = %d', profile_flag, plot_flag, h_indx));
hold off;

% Plot the ensemble mean
fig = figure('Position', [100, 100, 1000, 600]);

% get the observation data
obs = w(h_indx, :);

hens = ensemble_vec_mean(1:hdim, :); % Equivalent to ensemble_vec_mean[:hdim,:]
h_true = htrue(h_indx, :);
h_nurged = hnurge(h_indx, :);
h_ens_mean = hens(h_indx, :);
h_ens_mem = ensemble_vec_full(h_indx, :, :); % 3D array
h_ens_mem = permute(h_ens_mem, [3, 2, 1]);