%% -----------------------------------------------------------
% @author: 	Brian Kyanjo
% @date: 		2025-04-30
% @brief: 	Reads and plot results from both ISSM and ICESEE
% ------------------------------------------------------------

close all; clear all

% data_file_paths='data_0/_modelrun_datasets';
% data_file_paths='data_1';
% data_file_paths='_modelrun_datasets';
global data_file_paths 
data_file_paths = '_modelrun_datasets';


% time steps
% k_array = [10, 30, 50, 80,120];  % multiple time steps
% k_array=[25, 45, 60, 80, 101, 137, 162];
% k_array=[30, 40, 80,160, 180];
k_array=[20, 80,120, 160, 240, 350, 450];
% k_array=[1];
dt = 0.2;
global nvar;
nvar = 6;

make_plots = 0;
make_multi_plots = 1;
k = 65;
% k=28;

% Load the essential data
results_dir = 'results';
filter_type = 'true-wrong';
file_path   = fullfile(results_dir, sprintf('%s-issm.h5', filter_type));
t           = h5read(file_path,'/t');
ind_m       = h5read(file_path,'/obs_index');
tm_m        = h5read(file_path,'/obs_max_time');
run_mode    = h5read(file_path,'/run_mode');

% load the true and nurged states
file_path            = fullfile(data_file_paths, 'true_nurged_states.h5');
model_true_state     = h5read(file_path,'/true_state')';
model_nurged_state   = h5read(file_path, '/nurged_state')';

% load observation data
file_path  = fullfile(data_file_paths, 'synthetic_obs.h5');
w          = h5read(file_path, '/hu_obs')'; 

% load the ensemble data
file_path         = fullfile(data_file_paths, 'icesee_ensemble_data.h5');
ensemble_vec_full = h5read(file_path, '/ensemble'); 
ensemble_vec_mean = h5read(file_path, '/ensemble_mean')';

% Or read from .mat files (from the ISSM side)
% file_path_true= fullfile("issm_data","true_state.mat");
% file_path_nurged= fullfile("issm_data","nurged_state.mat");
% md_true_= loadmodel(file_path_true);
% md_nurged_= loadmodel(file_path_nurged);

% Process and plot
[ndim, nt] = size(model_true_state);
num_steps = nt - 1;
var_inputs = ['thickness', 'Vx', 'Vy', 'friction_coefficient', 'bed_topography'];
hdim = floor(ndim / nvar);  % dimension of one variable

% file_path   = fullfile("data", "ISMIP_initial_data.mat");
% file_path   = fullfile("data", "ISMIP.Parameterization1.mat");
file_path   = fullfile("data", "ISMIP.Parameterization1.mat");
md = loadmodel(file_path);
md_true = md; md_nurged = md;
md_mean = md; md_ens = md;
% dt = t(2) - t(1);
% dt = 0.25;

k_gl = 150;             % choose one of your k_array entries
[md_true, md_nurged, md_ens] = setup_model_states(k_gl, dt, ...
    model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md);

plot_gl_on_bed(md_true, md_nurged, md_ens, round((k_gl-1)*dt));

% for k_gl = 1:499
%     % k_gl = 20;             % choose one of your k_array entries
%     [md_true, md_nurged, md_ens] = setup_model_states(k_gl, dt, ...
%         model_true_state, model_nurged_state, ensemble_vec_mean, ...
%         md_true, md_nurged, md_ens, md);
% 
%     plot_gl_on_bed(md_true, md_nurged, md_ens, round((k_gl-1)*dt));
% end

frames_plot=0;

make_gl_movie = 1;           % set true if you want an mp4
frame_stride  = 10;              % plot every 20th step (adjust)
if frames_plot
    if make_gl_movie
        v = VideoWriter('velocity_gl_evolution.mp4','MPEG-4');
        v.FrameRate = 8;
        open(v);
    end
    
    for k_gl = 1:frame_stride:nt
        % Build true / nurged / ensemble models at this time
        [md_true_k, md_nurged_k, md_ens_k] = setup_model_states( ...
            k_gl, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
    
        % Just in case: cheap diagnostic to see if GLs differ at all
        phi_true   = md_true_k.mask.ocean_levelset;
        phi_nurged = md_nurged_k.mask.ocean_levelset;
        phi_ens    = md_ens_k.mask.ocean_levelset;
        fprintf('k=%4d (t = %6.2f yr): max|true-nurged| = %.3e, max|true-ens| = %.3e\n', ...
                k_gl, (k_gl-1)*dt, ...
                max(abs(phi_true - phi_nurged)), ...
                max(abs(phi_true - phi_ens)));
    
        % Plot this timestep
        plot_velocity_with_gl(md_true_k, md_nurged_k, md_ens_k, (k_gl-1)*dt);
        drawnow;
    
        if make_gl_movie
            frame = getframe(gcf);
            writeVideo(v, frame);
        else
            pause(0.15);   % slow down when just inspecting
        end
    end
    
    if make_gl_movie
        close(v);
    end
end

if make_multi_plots 
    % %{ % Thickness difference plots
    plot_var_diff(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'geometry.thickness', 'Thicknes', 'm');
    % 
    %  % thickness plot evolution
    plot_var_evolution(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'geometry.thickness', 'Thickness', 'm');
    % 
    % % base  
    % plot_var_diff(k_array, dt, ...
    %     model_true_state, model_nurged_state, ensemble_vec_mean, ...
    %     md_true, md_nurged, md_ens, md, ...
    %     'geometry.base', 'Base', 'm');
    % 
    %  % thickness plot evolution
    % plot_var_evolution(k_array, dt, ...
    %     model_true_state, model_nurged_state, ensemble_vec_mean, ...
    %     md_true, md_nurged, md_ens, md, ...
    %     'geometry.base', 'Base', 'm');
    % 
    plot_var_diff(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'geometry.surface', 'Surface', 'm');
    % 
    %  % thickness plot evolution
    plot_var_evolution(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'geometry.surface', 'Surface', 'm');
    % 
        % velocity evolution plots
    plot_var_evolution(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'initialization.vel', 'Velocity', 'm/s');
    % 
    plot_var_diff(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'initialization.vel', 'Velocity', 'm/s'); %}
        % bed topography evolution plots        
    plot_var_evolution(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'geometry.bed', 'Bed Elevation', 'm');
    plot_var_diff(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'geometry.bed', 'Bed', 'm');
    
        % friction coefficient evolution plots  
    plot_var_evolution(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'friction.coefficient', 'Friction Coefficient', '');
    plot_var_diff(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'friction.coefficient', 'Friction', 'm');
end 

if make_plots
    % thickness
    [md_true, md_nurged, md_ens] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
    plot_triptych(md_true, md_nurged, md_ens, ...
                'geometry.thickness', sprintf('Ice Thickness after %d years', round((k-1)*dt)), parula, 'm');  
    plot_triptych(md_true, md_nurged, md_ens, ...
                'geometry.surface', sprintf('Ice Surface after %d years', round((k-1)*dt)), parula, 'm');  
    % velocity              
    plot_triptych(md_true, md_nurged, md_ens, ...
                'initialization.vel', sprintf('Ice Velocity after %d years', round((k-1)*dt)), parula, 'm/s');   
    % bed topography
    plot_triptych(md_true, md_nurged, md_ens, ...
                'geometry.bed', sprintf('Bed Elevation after %d years', round((k-1)*dt)), parula, 'm');
    %
    % Frcition coefficient
    plot_triptych(md_true, md_nurged, md_ens, ...
                'friction.coefficient', sprintf('Friction Coefficient after %d years', round((k-1)*dt)), parula, '');
    %
    % % grounding line
    % plot_triptych(md_true, md_nurged, md_ens, ...
    %               'results.TransientSolution(499).MaskOceanLevelset', ...
    %               'Grounding Line', gray, '');
    % plot_triptych(md_true, md_nurged, md_ens, ...
    %             'mask.ocean_levelset', ...
    %             sprintf('Grounding Line after %d years', round((k-1)*dt)), parula, '');
end
% create a movie for the groundingline for every 10 yrs
% Setup video writer (optional)
make_movie = false;
if make_movie
    if make_movie
        v = VideoWriter('groundingline_triptych.mp4','MPEG-4');
        v.FrameRate = 10;  % frames per second
        open(v);
    end

    % Precompute density ratio
    di = md.materials.rho_ice / md.materials.rho_water;

    for k = 1:20:500
        
        % Update ensemble model from state vector
        md_ens.geometry.bed        = ensemble_vec_mean(hdim+1:2*hdim, k);
        md_ens.geometry.thickness  = ensemble_vec_mean(1:hdim, k);
        md_ens.friction.coefficient= ensemble_vec_mean(2*hdim+1:3*hdim, k);

        % Compute flotation mask for ensemble
        md_ens.mask.ocean_levelset = md_ens.geometry.thickness + md_ens.geometry.bed/di;

        % Fetch true & nudged grounding lines
        md_true.mask.ocean_levelset   = md_true_.results.TransientSolution(k).MaskOceanLevelset;
        md_nurged.mask.ocean_levelset = md_nurged_.results.TransientSolution(k).MaskOceanLevelset;

        % Plot triptych
        plot_triptych(md_true, md_nurged, md_ens, ...
                    'mask.ocean_levelset', ...
                    sprintf('Grounding Line after %d years', round((k-1)*0.5)), ...
                    parula, '');

        drawnow;

        % Capture frame if making movie
        if make_movie
            frame = getframe(gcf);
            writeVideo(v, frame);
        else
            pause(dt); % interactive view
        end
    end

    if make_movie
        close(v);
    end
end


%% -- helper functions -- %%
function plot_var_diff(k_array, dt, ...
    model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md, field, field_title, units)
    % =========================================================================
    % plot_var_diff
    % Automatically adapts the number of subplots to the length of k_array.
    % =========================================================================

    if nargin < 12, units = ''; end
    units_str = iff(~isempty(units), [' (' units ')'], '');
    nk = length(k_array);
    nrows = 2 + nk;  % (a) True + (b) No assimilation + (c...e) Assim steps

    figure('Position',[100 100 1000 150 + 150*nrows]); clf;

    % ---- Compute global colour limits across all requested steps ----
    all_data = [];

    for k = [1, k_array]          % include the true (initial) state
        [md_true_tmp, md_nurged_tmp, md_ens_tmp] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
        data_tmp = get_nested_field(md_ens_tmp, field);
        all_data = [all_data; data_tmp(:)];
        % Also include true state for diff limits
        data_true_tmp = get_nested_field(md_true_tmp, field);
        all_data = [all_data; data_true_tmp(:)];
        % nureged state for diff limits
        % data_nurged_tmp = get_nested_field(md_nurged_tmp, field);
        % all_data = [all_data; data_nurged_tmp(:)];
    end

    cmin = min(all_data);
    cmax = max(all_data);
   
    clear all_data

    % (a) True field
    [md_true, md_nurged, md_ens] = setup_model_states(1, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);
    data_true = get_nested_field(md_true, field);
    data_ens  = get_nested_field(md_ens, field);
    data_nurged = get_nested_field(md_nurged, field);
    % cmin = min([data_true(:); data_ens(:)]);
    % cmax = max([data_true(:); data_ens(:)]);
    plotmodel(md_true, 'data', data_true, ...
        'title', sprintf('(a) True %s', field_title), ...
        'subplot', [nrows, 1, 1], 'caxis', [cmin cmax], 'colorbar', 'off');

    % (b) No assimilation − True
    % diff_noassim = data_ens - data_true;
    
    eps0 = 0.001 * max(abs(data_true(:)));
    % diff_noassim = (data_ens - data_true)./(abs(data_true) + eps0);
    [md_true, md_nurged, md_ens] = setup_model_states(1, dt, ...
    model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md);
    data_true = get_nested_field(md_true, field);
    data_ens  = get_nested_field(md_ens, field);
    data_nurged = get_nested_field(md_nurged, field);
    diff_noassim = (data_ens - data_true); 
    % diff_noassim = abs(data_ens - data_true)./(abs(data_true));
    maxAbs_noassim = max(abs(diff_noassim(:)));
    maxAbs= maxAbs_noassim;
    % maxAbs_noassim = 2;
    plotmodel(md_ens, 'data', diff_noassim, ...
        'title', '(b) No assimilation − True', ...
        'subplot', [nrows, 1, 2], 'caxis', [-maxAbs_noassim maxAbs_noassim], 'colorbar', 'off');

    % (c...): Assimilation differences at each k
    for idx = 1:nk
        label = sprintf('(%c)', 'b' + idx);
        k = k_array(idx);
        [md_true, md_nurged, md_ens] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
        true_field = get_nested_field(md_true, field);
        ens_field  = get_nested_field(md_ens, field);
        
        eps0 = 0.001 * max(abs(true_field(:)));    % 1% of max for stabilization
        % diff_data = (ens_field - true_field) ./ (abs(true_field) + eps0);

        % diff_data = abs(get_nested_field(md_ens, field) - get_nested_field(md_true, field))./abs(get_nested_field(md_true, field));
        diff_data = (get_nested_field(md_ens, field) - get_nested_field(md_true, field));
        % diff_data = (get_nested_field(md_ens, field) - get_nested_field(md_true, field));
        maxAbs = max(abs(diff_data(:)));
        % maxAbs=2;
        plotmodel(md_ens, 'data', diff_data, ...
            'title', sprintf('%s Assimilated − True (after %.1f years)', label, round((k-1)*dt)), ...
            'subplot', [nrows, 1, idx + 2], 'caxis', [-maxAbs maxAbs], 'colorbar', 'off');
    end

    % Adjust layout
    axs = flipud(findall(gcf,'Type','axes'));
    % --- Adaptive layout scaling ---
    gap = 0.02;              % small positive gap
    top = 0.95; bottom = 0.08;
    available_height = top - bottom - (nrows-1)*gap;
    height = available_height / nrows;

    % if too small (many rows), expand figure height automatically
    if height < 0.05
        fig = gcf;
        scale_factor = max(1, ceil(0.05 / height));  % ensure visible spacing
        fig.Position(4) = fig.Position(4) * scale_factor;  % increase figure height
        height = 0.05;  % set to minimum safe height
    end


    for i = 1:nrows
        pos = [0.10, bottom+(nrows-i)*(height+gap), 0.70, height];
        set(axs(i), 'Position', pos, 'FontWeight', 'bold', ...
            'LineWidth', 1.2, 'Box', 'on', 'TickDir', 'out', ...
            'Layer', 'top', 'FontSize', 11, 'TickLength',[0.005 0.005]);
        ylabel(axs(i),'y (km)','FontWeight','bold');
        if i < nrows
            set(axs(i),'XTickLabel',[]);
        else
            xlabel(axs(i),'x (km)','FontWeight','bold');
        end
    end

    % Colorbars
    cb1 = colorbar(axs(1), 'Position',[0.83 0.68 0.025 0.16]);
    ylabel(cb1,[field_title units_str],'FontSize',12,'FontWeight','bold');
    colormap(axs(1), parula);
    for i = 2:nrows, colormap(axs(i), redblue(256)); end
    cb2 = colorbar(axs(end), 'Position',[0.83 0.25 0.025 0.40]);
    % ylabel(cb2,['Δ' field_title units_str],'FontSize',12,'FontWeight','bold');
    ylabel(cb2,['Relative Error'],'FontSize',12,'FontWeight','bold');
    set(gcf,'Color','w');
end


function plot_var_evolution(k_array, dt, ...
    model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md, field, field_title, units)
    % =========================================================================
    % plot_var_evolution
    % Automatically adapts subplot layout to number of k_array elements.
    % =========================================================================

    if nargin < 13, field_title = field; end
    if nargin < 14, units = ''; end
    units_str = iff(~isempty(units), [' (' units ')'], '');
    nk = length(k_array);
    nrows = 2 + nk;

    figure('Position',[100 100 1000 150 + 150*nrows]); clf;

    % ---- Compute global colour limits across all requested steps ----
    all_data = [];

    % for k = [1, k_array]          % include the true (initial) state
    for k = [k_array(end), k_array]
        [md_true_tmp, md_nurged_tmp, md_ens_tmp] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
        data_tmp = get_nested_field(md_ens_tmp, field);
        all_data = [all_data; data_tmp(:)];
        % Also include true state for diff limits
        data_true_tmp = get_nested_field(md_true_tmp, field);
        all_data = [all_data; data_true_tmp(:)];
        % nureged state for diff limits
        % data_nurged_tmp = get_nested_field(md_nurged_tmp, field);
        % all_data = [all_data; data_nurged_tmp(:)];
    end

    cmin = min(all_data);
    cmax = max(all_data);
    clear all_data

    % (a) True
    [md_true, md_nurged, md_ens] = setup_model_states(k_array(end), dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);
    data_true = get_nested_field(md_true, field);
    data_ens  = get_nested_field(md_ens, field);
    data_nurged = get_nested_field(md_nurged, field);
    % cmin = min([data_true(:); data_ens(:)]);
    % cmax = max([data_true(:); data_ens(:)]);
    plotmodel(md_true, 'data', data_true, ...
        'title', sprintf('(a) True %s (after %.1f years)', field_title, round(k_array(end)-1)*dt), ...
        'subplot', [nrows, 1, 1], 'caxis', [cmin cmax], 'colorbar', 'off');

    % (b) No assimilation
    [md_true, md_nurged, md_ens] = setup_model_states(1, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);
    data_true = get_nested_field(md_true, field);
    data_ens  = get_nested_field(md_ens, field);
    data_nurged = get_nested_field(md_nurged, field);
    plotmodel(md_ens, 'data', data_ens, ...
        'title', sprintf('(b) No assimilation %s', field_title), ...
        'subplot', [nrows, 1, 2], 'caxis', [cmin cmax], 'colorbar', 'off');

    % (c...): Assimilated snapshots
    for idx = 1:nk
        k = k_array(idx);
        label = sprintf('(%c)', 'b' + idx);
        [md_true, md_nurged, md_ens] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
        data_ens = get_nested_field(md_ens, field);
        plotmodel(md_ens, 'data', data_ens, ...
            'title', sprintf('%s Assimilated %s (after %.1f years)', label, field_title, round((k-1)*dt)), ...
            'subplot', [nrows, 1, idx + 2], 'caxis', [cmin cmax], 'colorbar', 'off');
    end

    % Layout
    axs = flipud(findall(gcf,'Type','axes'));
    % --- Adaptive layout scaling ---
    gap = 0.02;              % small positive gap
    top = 0.95; bottom = 0.08;
    available_height = top - bottom - (nrows-1)*gap;
    height = available_height / nrows;

    % if too small (many rows), expand figure height automatically
    if height < 0.05
        fig = gcf;
        scale_factor = max(1, ceil(0.05 / height));  % ensure visible spacing
        fig.Position(4) = fig.Position(4) * scale_factor;  % increase figure height
        height = 0.05;  % set to minimum safe height
    end


    for i = 1:nrows
        pos = [0.10, bottom+(nrows-i)*(height+gap), 0.70, height];
        set(axs(i),'Position',pos, ...
            'FontWeight','bold','LineWidth',1.2,'Box','on', ...
            'TickDir','out','Layer','top','FontSize',11, ...
            'TickLength',[0.005 0.005]);
        ylabel(axs(i),'y (km)','FontWeight','bold');
        if i < nrows
            set(axs(i),'XTickLabel',[]);
        else
            xlabel(axs(i),'x (km)','FontWeight','bold');
        end
    end

    % Colorbars
    % cb1 = colorbar(axs(1), 'Position',[0.83 0.71 0.025 0.16]);
    % ylabel(cb1,[field_title units_str],'FontSize',12,'FontWeight','bold');
    % for i = 2:nrows, colormap(axs(i), parula); end
    % cb2 = colorbar(axs(end), 'Position',[0.83 0.24 0.025 0.45]);
    % ylabel(cb2,[field_title units_str],'FontSize',12,'FontWeight','bold');
    % set(gcf,'Color','w');

    for i = 1:nrows, colormap(axs(i), parula); end
    cb = colorbar(axs(end), 'Position',[0.83 0.25 0.025 0.45]);
    ylabel(cb,[field_title units_str],'FontSize',12,'FontWeight','bold');

end


function plot_triptych(md_true, md_nurged, md_ens, field, field_title, cmap, units)
% Compare true, nudged, assimilated, and difference with two separate colorbars
    global k_array;
    global dt;
    if nargin < 6 || isempty(cmap), cmap = parula; end
    if nargin < 7, units = ''; end
    units_str = iff(~isempty(units), [' (' units ')'], '');

    % --- Data ---
    data_true   = get_nested_field(md_true, field);
    data_nurged = get_nested_field(md_nurged, field);
    data_ens    = get_nested_field(md_ens, field);
    diff_data   = data_ens - data_true;
    % 
    %     % ---- Compute global colour limits across all requested steps ----
    % all_data = [];
    % 
    % for k = [1, k_array]          % include the true (initial) state
    %     [md_true_tmp, ~, md_ens_tmp] = setup_model_states(k, dt, ...
    %         model_true_state, model_nurged_state, ensemble_vec_mean, ...
    %         md_true, md_nurged, md_ens, md);
    %     data_tmp = get_nested_field(md_ens_tmp, field);
    %     all_data = [all_data; data_tmp(:)];
    %     % Also include true state for diff limits
    %     data_true_tmp = get_nested_field(md_true_tmp, field);
    %     all_data = [all_data; data_true_tmp(:)];
    % end
    % 
    % cmin = min(all_data);
    % cmax = max(all_data);
    % clear all_data

    % --- Limits ---
    cmin   = min([data_true(:); data_nurged(:); data_ens(:)]);
    cmax   = max([data_true(:); data_nurged(:); data_ens(:)]);
    maxAbs = max(abs(diff_data(:)));

    figure('Position',[100 100 1000 800]); clf;

    % 1) True
    plotmodel(md_true,'data',data_true,'title',['True ' field_title], ...
        'subplot',[4,1,1],'caxis',[cmin cmax],'colorbar','off');

    % 2) Wrong
    plotmodel(md_nurged,'data',data_nurged,'title',['Wrong ' field_title], ...
        'subplot',[4,1,2],'caxis',[cmin cmax],'colorbar','off');

    % 3) Assimilated
    plotmodel(md_ens,'data',data_ens,'title',['Assimilated ' field_title], ...
        'subplot',[4,1,3],'caxis',[cmin cmax],'colorbar','off');

    % 4) Difference
    plotmodel(md_ens,'data',diff_data, ...
        'title',['(Assmilated - True) ' field_title], ...
        'subplot',[4,1,4],'caxis',[-maxAbs maxAbs],'colorbar','off');

    % --- Axes layout ---
    axs = flipud(findall(gcf,'Type','axes'));   % 1..4 top->bottom
    gap = -0.255; top = 0.94; bottom = 0.08;    % tightened spacing
    height = (top-bottom - 3*gap)/4;

    for i = 1:4
        pos = [0.10, bottom+(4-i)*(height+gap), 0.70, height];
        set(axs(i),'Position',pos, ...
            'FontWeight','bold','LineWidth',1.5,'Box','on', ...
            'TickDir','out','TickLength',[0.005 0.005], ...
            'Layer','top');
        ylabel(axs(i),'Y (km)','FontSize',12,'FontWeight','bold');
        if i < 4
            set(axs(i),'XTickLabel',[]);  % only bottom plot shows X
        else
            xlabel(axs(i),'X (km)','FontSize',12,'FontWeight','bold');
        end
    end


    % --- First colorbar (shortened for top 3) ---
    for i = 1:3, colormap(axs(i), cmap); caxis(axs(i), [cmin cmax]); end
    cb1 = colorbar(axs(2),'Position',[0.83 0.415 0.025 0.35]); % shorter
    static_field = regexprep(field_title,'\s+after.*','');
    ylabel(cb1,[static_field units_str],'FontSize',13,'FontWeight','bold');
    cb1.FontSize = 11;
    set(cb1,'Box','on','LineWidth',1.2);

    % --- Second colorbar (shortened for difference) ---
    ax_diff = axs(4);
    colormap(ax_diff, redblue(256));
    caxis(ax_diff,[-maxAbs maxAbs]);
    pos_diff = get(ax_diff,'Position');
    cb2 = colorbar(ax_diff,'Position',[0.83 pos_diff(2)+0.14 0.025 pos_diff(4)-0.28]);
    % ylabel(cb2, ['$\mathbf{\Delta}$' static_field units_str], ...
    %    'Interpreter','latex', ...
    %    'FontSize',13);
    ylabel(cb2, ['Difference' units_str], ...
       'FontSize',13);
    cb2.FontSize = 11;
    set(cb2,'Box','on','LineWidth',1.2);
end

% ---- helpers ----
function out = get_nested_field(s, field)
    parts = strsplit(field,'.'); out = s;
    for i = 1:numel(parts)
        tok = parts{i};
        t = regexp(tok,'(.+)\((\d+)\)$','tokens');
        if ~isempty(t), out = out.(t{1}{1})(str2double(t{1}{2}));
        else, out = out.(tok);
        end
    end
end

function y = iff(c,a,b), if c, y=a; else, y=b; end, end

function cmap = redblue(n)
    if nargin<1, n=256; end
    m = n/2;
    r=[linspace(0,1,m) ones(1,m)];
    g=[linspace(0,1,m) linspace(1,0,m)];
    b=[ones(1,m) linspace(1,0,m)];
    cmap=[r(:) g(:) b(:)];
end

function [md_true, md_nurged, md_ens] = setup_model_states(k, dt, model_true_state, model_nurged_state, ensemble_vec_mean, md_true, md_nurged, md_ens, md)
% =========================================================================
% setup_model_states
%
% Purpose:
%   Initialize the true, nurged, and ensemble models at time index k.
%
% Inputs:
%   k                  - Time index (integer)
%   dt                 - Time step (scalar, e.g., 0.15)
%   model_true_state   - Matrix of true model states [n_state_vars x nt]
%   model_nurged_state - Matrix of nurged model states [n_state_vars x nt]
%   ensemble_vec_mean  - Matrix of ensemble mean states [n_state_vars x nt]
%   md_true, md_nurged, md_ens - Model structures (initialized ISSM models)
%   md                 - Reference model (for materials, constants)
%
% Outputs:
%   md_true, md_nurged, md_ens - Updated model structures at step k
%
% Author:  Brian Kyanjo
% Date:    2025-11-03
% =========================================================================

    % --- Basic setup ---
    global nvar;
    % nvar =6;
    hdim = length(model_true_state(:,1)) / nvar;  % Assuming 5 state components
    di = md.materials.rho_ice / md.materials.rho_water;
    % hdim = 3347;

    %% === TRUE STATE ===
    True_thickness = model_true_state(1:hdim, k);
    % True_base = model_true_state(hdim+1:2*hdim, k);
    True_surface = model_true_state(hdim+1:2*hdim, k);
    True_base = True_surface - True_thickness;
    Vx = model_true_state(2*hdim+1:3*hdim, k);
    Vy = model_true_state(3*hdim+1:4*hdim, k);
    Vel = sqrt(Vx.^2 + Vy.^2);
    True_bed = model_true_state(4*hdim+1:5*hdim, k);
    True_fcoeff = model_true_state(5*hdim+1:6*hdim, k);

    md_true.geometry.thickness = True_thickness;
    md_true.geometry.base = True_base;
    md_true.geometry.surface = True_surface;
    md_true.initialization.vx  = Vx;
    md_true.initialization.vy  = Vy;
    md_true.initialization.vel = Vel;
    md_true.geometry.bed       = True_bed;
    md_true.friction.coefficient = True_fcoeff;
    md_true.mask.ocean_levelset = True_thickness + True_bed / di;
    % md_true.mask.ocean_levelset = di * md_true.geometry.thickness + md_true.geometry.bed;


    %% === NURGED STATE ===
    nurged_thickness = model_nurged_state(1:hdim, k);
    % nurged_base = model_nurged_state(hdim+1:2*hdim, k);
    nurged_surface = model_nurged_state(hdim+1:2*hdim, k);
    nurged_base = nurged_surface - nurged_thickness;
    Vx = model_nurged_state(2*hdim+1:3*hdim, k);
    Vy = model_nurged_state(3*hdim+1:4*hdim, k);
    Vel = sqrt(Vx.^2 + Vy.^2);
    nurged_bed = model_nurged_state(4*hdim+1:5*hdim, k);
    nurged_fcoeff = model_nurged_state(5*hdim+1:6*hdim, k);
    % nurged_thickness = nurged_surface - nurged_bed;

    md_nurged.geometry.thickness = nurged_thickness;
    md_nurged.geometry.surface = nurged_surface;
    md_nurged.initialization.vx  = Vx;
    md_nurged.initialization.vy  = Vy;
    md_nurged.initialization.vel = Vel;
    md_nurged.geometry.bed       = nurged_bed;
    md_nurged.friction.coefficient = nurged_fcoeff;
    md_nurged.mask.ocean_levelset = nurged_thickness + nurged_bed / di;
    % md_nurged.mask.ocean_levelset = di * md_nurged.geometry.thickness + md_nurged.geometry.bed;

    %% === ENSEMBLE MEAN STATE ===
    ens_thickness = ensemble_vec_mean(1:hdim, k);
    % ens_base = ensemble_vec_mean(hdim+1:2*hdim, k);
    ens_surface = ensemble_vec_mean(hdim+1:2*hdim, k);
    ens_base = ens_surface - ens_thickness;
    Vx = ensemble_vec_mean(2*hdim+1:3*hdim, k);
    Vy = ensemble_vec_mean(3*hdim+1:4*hdim, k);
    Vel = sqrt(Vx.^2 + Vy.^2);
    ens_bed = ensemble_vec_mean(4*hdim+1:5*hdim, k);
    % ens_bed = ensemble_vec_mean(2*hdim+1:3*hdim, k);
    % ens_fcoeff = ensemble_vec_mean(5*hdim+1:6*hdim, k);
    % read friction from temp file
    icesee_path='/Users/bkyanjo3/da_project/ICESEE/applications/issm_model/examples/ISMIP_Choi';
    % data_path= '_modelrun_datasets';
    global data_file_paths;
    h5file = fullfile(icesee_path, data_file_paths, sprintf('temp_coefficient_%d.h5', 0));
    ens_fcoeff = h5read(h5file, '/coefficient'); ens_fcoeff = ens_fcoeff(k,:)';
    % Vx = h5read(h5file, '/Vx');
    % Vy = h5read(h5file, '/Vy');
    % Vx = Vx(k,:)'; Vx = Vx(k,:)';
    % Vel = sqrt(Vx.^2 + Vy.^2);


    % ens_thickness = ens_surface - ens_bed;

    md_ens.geometry.thickness = ens_thickness;
    md_ens.geometry.base = ens_base;
    md_ens.geometry.surface = ens_surface;
    md_ens.initialization.vx  = Vx;
    md_ens.initialization.vy  = Vy;
    md_ens.initialization.vel = Vel;
    md_ens.geometry.bed       = ens_bed;
    md_ens.friction.coefficient = ens_fcoeff;
    md_ens.mask.ocean_levelset = ens_thickness + ens_bed / di;
    % md_ens.mask.ocean_levelset = di * md_ens.geometry.thickness + md_ens.geometry.bed;

end


function plot_gl_on_bed(md_true, md_nurged, md_ens, t_years)
% Plot bed + smooth grounding lines (true / wrong / ens) as colored contours.

    % --- Background: bed topography ---
    % bed_true = md_true.geometry.bed;
    bed_true = md_true.initialization.vel;
    cmin = min(bed_true(:));
    cmax = max(bed_true(:));

    x = md_true.mesh.x(:);
    y = md_true.mesh.y(:);

    figure('Position',[100 100 1000 280]); clf;
    plotmodel(md_true,'data',bed_true,...
        'title',sprintf('Velocity + Grounding lines (t = %.1f years)', t_years), ...
        'caxis',[cmin cmax],'colorbar','on');
    hold on;

    % -----------------------------------------------------------------
    % 1) Build a regular grid for contouring
    % -----------------------------------------------------------------
    nx = 400;           % # of points along x (increase for smoother lines)
    ny = 60;            % # of points along y

    xg = linspace(min(x), max(x), nx);
    yg = linspace(min(y), max(y), ny);
    [Xg, Yg] = meshgrid(xg, yg);

    % Ocean level set (φ = 0 is grounding line)
    phi_true   = md_true.mask.ocean_levelset(:);
    phi_wrong  = md_nurged.mask.ocean_levelset(:);
    phi_ens    = md_ens.mask.ocean_levelset(:);

    % Interpolate to grid
    Phi_true  = griddata(x, y, phi_true,  Xg, Yg, 'linear');
    Phi_wrong = griddata(x, y, phi_wrong, Xg, Yg, 'linear');
    Phi_ens   = griddata(x, y, phi_ens,   Xg, Yg, 'linear');

    % -----------------------------------------------------------------
    % 2) Contour the zero level (grounding line) in different colours
    % -----------------------------------------------------------------
    % TRUE GL – thick black
    [~, h1] = contour(Xg, Yg, Phi_true,  [0 0], 'k-', 'LineWidth', 2.0);

    % NO-ASSIMILATION GL – dashed red
    [~, h2] = contour(Xg, Yg, Phi_wrong, [0 0], 'r--', 'LineWidth', 1.8);

    % ASSIMILATED GL – solid cyan
    [~, h3] = contour(Xg, Yg, Phi_ens,   [0 0], 'c-', 'LineWidth', 1.8);

    % -----------------------------------------------------------------
    % 3) Cosmetics for poster
    % -----------------------------------------------------------------
    axis equal tight
    xlabel('x (m)','FontWeight','bold');
    ylabel('y (m)','FontWeight','bold');

    set(gca,'FontWeight','bold','LineWidth',1.4,'Box','on', ...
        'TickDir','out','Layer','top','TickLength',[0.006 0.006]);

    legend([h1(1), h2(1), h3(1)], ...
           {'True GL','No assimilation GL','Assimilated GL'}, ...
           'Location','southwest','FontSize',11,'Box','on');

    set(gcf,'Color','w');
end

function plot_velocity_with_gl(md_true, md_nurged, md_ens, time_yrs)

    figure(10); clf
    plotmodel(md_ens, 'data', md_ens.initialization.vel, ...
        'title', sprintf('Velocity + Grounding lines (t = %.1f years)', time_yrs), ...
        'colorbar','on');
    hold on

    ok_true   = plot_gl_line(md_true,   'k', '-', 2.5);
    ok_nurg   = plot_gl_line(md_nurged, 'r', '--',2.5);
    ok_assim  = plot_gl_line(md_ens,    'c', '-', 2.5);

    % legend only for the ones that actually exist
    ax = gca;
    h = []; labels = {};
    if ok_true
        h(end+1) = plot(ax, NaN, NaN, 'k-',  'LineWidth', 2.5);
        labels{end+1} = 'True GL';
    end
    if ok_nurg
        h(end+1) = plot(ax, NaN, NaN, 'r--', 'LineWidth', 2.5);
        labels{end+1} = 'No assimilation GL';
    end
    if ok_assim
        h(end+1) = plot(ax, NaN, NaN, 'c-',  'LineWidth', 2.5);
        labels{end+1} = 'Assimilated GL';
    end
    if ~isempty(h)
        legend(h, labels, 'Location','southwest');
    end

    xlabel('x (m)','FontWeight','bold');
    ylabel('y (m)','FontWeight','bold');
    set(ax,'FontWeight','bold','LineWidth',1.2,'TickDir','out');
end

function ok = plot_gl_line(md, line_color, line_style, lw)
%PLOT_GL_LINE Plot grounding line from mask.ocean_levelset = 0
%   ok = plot_gl_line(md, color, style, lw)
%   returns ok = true if a GL was actually plotted, false otherwise.

    if nargin < 4, lw = 2.0; end
    if nargin < 3, line_style = '-'; end
    if nargin < 2, line_color = 'k'; end

    x   = md.mesh.x(:);
    y   = md.mesh.y(:);
    phi = md.mask.ocean_levelset(:);

    % quick sanity check: if phi has no sign change, there is no GL
    if all(phi >= 0) || all(phi <= 0) || all(~isfinite(phi))
        fprintf('[plot_gl_line] no sign change in ocean_levelset (no GL); skipping.\n');
        ok = false;
        return;
    end

    % --- interpolate to a regular grid for a smooth contour ---
    Nx = 400; Ny = 60;
    xg = linspace(min(x), max(x), Nx);
    yg = linspace(min(y), max(y), Ny);
    [Xg, Yg] = meshgrid(xg, yg);

    F = scatteredInterpolant(x, y, phi, 'natural', 'nearest');
    Phig = F(Xg, Yg);

    % check again after interpolation (can become “almost constant”)
    if max(Phig(:)) * min(Phig(:)) > 0
        fprintf('[plot_gl_line] interpolated mask has no sign change; skipping.\n');
        ok = false;
        return;
    end

    hold on
    [C, h] = contour(Xg, Yg, Phig, [0 0], ...
                     'Color', line_color, ...
                     'LineStyle', line_style, ...
                     'LineWidth', lw);

    if isempty(h) || ~isvalid(h)
        fprintf('[plot_gl_line] contour returned empty handle; skipping.\n');
        ok = false;
    else
        ok = true;
    end
end