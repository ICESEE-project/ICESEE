%% -----------------------------------------------------------
% @author: 	Brian Kyanjo
% @date: 		2025-04-30
% @brief: 	Reads and plot results from both ISSM and ICESEE
% ------------------------------------------------------------

close all; clear all

% data_file_paths='data_0/_modelrun_datasets';
% data_file_paths='data_1';
% data_file_paths='_modelrun_datasets';
data_file_paths = 'test_data_working_0'
% data_file_paths='cluster_data/_modelrun_datasets';
% data_file_paths='test_data_al1';
% data_file_paths='test_data_5';

% time steps
% k_array = [3, 10, 20, 30, 50];  % multiple time steps
k_array=[5,15, 30, 45, 62, 80];
dt = 0.2;

make_plots = false;
make_multi_plots = true;
k = 499;
% k=1;

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
hdim = floor(ndim / 6);  % dimension of one variable

file_path   = fullfile("data", "ISMIP_initial_data.mat");
md = loadmodel(file_path);
md_true = md; md_nurged = md;
md_mean = md; md_ens = md;
% dt = t(2) - t(1);
% dt = 0.25;

if make_multi_plots 
    (* % Thickness difference plots
    plot_var_diff(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'geometry.thickness', 'Thicknes', 'm');
    
     % thickness plot evolution
    plot_var_evolution(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'geometry.thickness', 'Thickness', 'm');
    
    % base  
    plot_var_diff(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'geometry.base', 'Base', 'm');
    
     % thickness plot evolution
    plot_var_evolution(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'geometry.base', 'Base', 'm');
    
    plot_var_diff(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'geometry.surface', 'Surface', 'm');
    
     % thickness plot evolution
    plot_var_evolution(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'geometry.surface', 'Surface', 'm');
    
        % velocity evolution plots
    plot_var_evolution(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'initialization.vel', 'Velocity', 'm/s');
    
    plot_var_diff(k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'initialization.vel', 'Velocity', 'm'); *)
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
    plot_triptych(md_true, md_nurged, md_ens, ...
                'geometry.thickness', sprintf('Ice Thickness after %d years', round((k-1)*dt)), parula, 'm');  
    % velocity              
    plot_triptych(md_true, md_nurged, md_ens, ...
                'initialization.vel', sprintf('Ice Velocity after %d years', round((k-1)*dt)), parula, 'm/s');   
    % bed topography
    plot_triptych(md_true, md_nurged, md_ens, ...
                'geometry.bed', sprintf('Bed Elevation after %d years', round((k-1)*dt)), parula, 'm');
    %
    % % Frcition coefficient
    plot_triptych(md_true, md_nurged, md_ens, ...
                'friction.coefficient', sprintf('Friction Coefficient after %d years', round((k-1)*dt)), parula, '');
    %
    % % grounding line
    % plot_triptych(md_true, md_nurged, md_ens, ...
    %               'results.TransientSolution(499).MaskOceanLevelset', ...
    %               'Grounding Line', gray, '');
    plot_triptych(md_true, md_nurged, md_ens, ...
                'mask.ocean_levelset', ...
                sprintf('Grounding Line after %d years', round((k-1)*dt)), parula, '');
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
        data_nurged_tmp = get_nested_field(md_nurged_tmp, field);
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
    diff_noassim = (data_nurged - data_true);
    % eps0 = 0.01 * max(abs(data_true(:)));
    % diff_noassim = (data_ens - data_true)./(abs(data_true) + eps0);
    maxAbs_noassim = max(abs(diff_noassim(:)));
    plotmodel(md_nurged, 'data', diff_noassim, ...
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
        
        eps0 = 0.01 * max(abs(true_field(:)));    % 1% of max for stabilization
        % diff_data = (ens_field - true_field) ./ (abs(true_field) + eps0);

        % diff_data = (get_nested_field(md_ens, field) - get_nested_field(md_true, field))./get_nested_field(md_true, field);
        diff_data = (get_nested_field(md_ens, field) - get_nested_field(md_true, field));
        maxAbs = max(abs(diff_data(:)));
        plotmodel(md_ens, 'data', diff_data, ...
            'title', sprintf('%s Assimilated − True (after %.1f years)', label, (k-1)*dt), ...
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
    ylabel(cb2,['Δ' field_title units_str],'FontSize',12,'FontWeight','bold');
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
        data_nurged_tmp = get_nested_field(md_nurged_tmp, field);
        % all_data = [all_data; data_nurged_tmp(:)];
    end

    cmin = min(all_data);
    cmax = max(all_data);
    clear all_data

    % (a) True
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

    % (b) No assimilation
    plotmodel(md_nurged, 'data', data_nurged, ...
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
            'title', sprintf('%s Assimilated %s (after %.1f years)', label, field_title, (k-1)*dt), ...
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

        % ---- Compute global colour limits across all requested steps ----
    all_data = [];

    for k = [1, k_array]          % include the true (initial) state
        [md_true_tmp, ~, md_ens_tmp] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
        data_tmp = get_nested_field(md_ens_tmp, field);
        all_data = [all_data; data_tmp(:)];
        % Also include true state for diff limits
        data_true_tmp = get_nested_field(md_true_tmp, field);
        all_data = [all_data; data_true_tmp(:)];
    end

    cmin = min(all_data);
    cmax = max(all_data);
    clear all_data

    % --- Limits ---
    % cmin   = min([data_true(:); data_nurged(:); data_ens(:)]);
    % cmax   = max([data_true(:); data_nurged(:); data_ens(:)]);
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
    hdim = length(model_true_state(:,1)) / 6;  % Assuming 5 state components
    di = md.materials.rho_ice / md.materials.rho_water;

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
    md_true.mask.ocean_levelset = True_surface + True_bed / di;


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

    %% === ENSEMBLE MEAN STATE ===
    ens_thickness = ensemble_vec_mean(1:hdim, k);
    % ens_base = ensemble_vec_mean(hdim+1:2*hdim, k);
    ens_surface = ensemble_vec_mean(hdim+1:2*hdim, k);
    ens_base = ens_surface - ens_thickness;
    Vx = ensemble_vec_mean(2*hdim+1:3*hdim, k);
    Vy = ensemble_vec_mean(3*hdim+1:4*hdim, k);
    Vel = sqrt(Vx.^2 + Vy.^2);
    ens_bed = ensemble_vec_mean(4*hdim+1:5*hdim, k);
    ens_fcoeff = ensemble_vec_mean(5*hdim+1:6*hdim, k);
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

end

