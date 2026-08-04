%% WBF / EBF / IBF comparison figure
% Produces the requested seven-panel comparison:
%   (a-d) true, fixed-wrong-friction WBF, EnKF-only EBF, and hybrid IBF
%         basal friction
%   (e-g) velocity, surface-elevation, and green-point GL RMSE vs time
%
% Put this script beside the three dataset folders, or change root_dir below.
% The original sample script is not modified.

close all; clearvars; clc

%% ---------------- User settings ----------------------------------------
% Resolve data relative to this script, so it works regardless of MATLAB's
% current working directory.
script_dir = fileparts(mfilename('fullpath'));
root_dir = script_dir;

% run_def = struct( ...
%     'key',    {'WBF','EBF','IBF'}, ...
%     'folder', {'_modelrun_datasets_method_comparison_40yr_wbf', ...
%                '_modelrun_datasets_method_comparison_40yr_ebf', ...
%                '_modelrun_datasets_method_comparison_40yr_ibf'}, ...
%     'title',  {'WBF: EnKF state/bed update with fixed wrong friction', ...
%                'EBF: EnKF-only friction recovery', ...
%                'IBF: EnKF state/bed update plus friction inversion'});

run_def = struct( ...
    'key',    {'WBF','EBF','IBF'}, ...
    'folder', {'_modelrun_datasets_p3q0_comparison_wbf_40yr', ...
               '_modelrun_datasets_p3q0_comparison_ebf_40yr', ...
               '_modelrun_datasets_p3q0_comparison_ibf_40yr'}, ...
    'title',  {'WBF: EnKF state/bed update with fixed wrong friction', ...
               'EBF: EnKF-only friction recovery', ...
               'IBF: EnKF state/bed update plus friction inversion'});

dt_fallback = 0.1;                % years, used only if no /t dataset is found
assimilation_times = 2:1:24;      % annual state observations; [] hides guides
rho_ice_fallback   = 917;         % kg m^-3
rho_water_fallback = 1028;        % kg m^-3
% Keep rendering a clearly marked diagnostic preview when an experiment
% fails the method-separation check. Set true to stop before plotting.
fail_on_uncoupled_ebf = false;
grounded_rmse_zoom_column = false; % show only the main velocity/surface RMSE panels
gl_rmse_zoom_panel = false;       % the full GL panel already resolves the method spread
show_post_da_metric_blocks = false; % keep zoom curves unobscured; report values in LaTeX table
velocity_rmse_zoom_limits = [25 120];
velocity_rmse_zoom_ticks = [25 50 75 100];
surface_rmse_zoom_limits = [40 85]; %[25 55];
surface_rmse_zoom_ticks = [40 55 70 85]; %[25 35 45 55];
gl_rmse_zoom_limits = [0 8];
gl_rmse_zoom_ticks = 0:8;
surface_rmse_domain = 'grounded_excluding_gl';
grounded_excluding_gl_xmax = 300e3; % upstream grounded ice, as in read_results_0.m

comparison_output_dir = fullfile(root_dir, '_modelrun_datasets_method_comparison_40yr');
figure_dir = fullfile(comparison_output_dir, 'figures');
if ~exist(figure_dir,'dir'), mkdir(figure_dir); end

%% ---------------- Locate the mesh/model -------------------------------
model_candidates = { ...
    fullfile(root_dir,'data','ISMIP.Parameterization1.mat'), ...
    fullfile(root_dir,'ISMIP.Parameterization1.mat'), ...
    fullfile(fileparts(root_dir),'data','ISMIP.Parameterization1.mat')};
model_file = first_existing_file(model_candidates);
if isempty(model_file)
    error(['Could not find ISMIP.Parameterization1.mat. Put it in data/ ' ...
           'under root_dir, or add its location to model_candidates.']);
end
md = read_issm_model(model_file);

x = md.mesh.x(:);
y = md.mesh.y(:);
elements = double(md.mesh.elements);
n_nodes = numel(x);
nvar = 6;

rho_ice = rho_ice_fallback;
rho_water = rho_water_fallback;
if (isstruct(md) && isfield(md,'materials')) || ...
   (isobject(md) && isprop(md,'materials'))
    try, rho_ice = md.materials.rho_ice; catch, end
    try, rho_water = md.materials.rho_water; catch, end
end
density_ratio = rho_ice/rho_water;

%% ---------------- Read runs and calculate diagnostics ------------------
n_runs = numel(run_def);
runs = repmat(struct(), n_runs, 1);

for r = 1:n_runs
    data_dir = fullfile(root_dir, run_def(r).folder);
    required = { ...
        fullfile(data_dir,'true_nurged_states.h5'), ...
        fullfile(data_dir,'icesee_ensemble_data.h5')};
    missing = required(~cellfun(@isfile,required));
    if ~isempty(missing)
        error('Missing required file for %s:\n  %s', ...
            run_def(r).key, strjoin(missing,'\n  '));
    end

    truth = read_state_matrix(required{1}, '/true_state', n_nodes, nvar);
    background = read_state_matrix(required{1}, '/nurged_state', n_nodes, nvar);
    ensemble = read_state_matrix(required{2}, '/ensemble_mean', n_nodes, nvar);

    % All experiments assimilate the same surface, velocity, and grounded-bed
    % observations. They differ only in friction treatment: WBF freezes it,
    % EBF updates it through EnKF cross-covariances, and IBF uses inversion.
    estimate = ensemble;

    nt = min([size(truth,2),size(estimate,2),size(background,2)]);
    truth = truth(:,1:nt);
    estimate = estimate(:,1:nt);

    runs(r).key = run_def(r).key;
    runs(r).title = run_def(r).title;
    runs(r).time = read_time_vector(root_dir, data_dir, nt, dt_fallback);
    runs(r).time = runs(r).time(1:nt);
    runs(r).truth = truth;
    runs(r).estimate = estimate;
    runs(r).background = background(:,1:nt);
    runs(r).valid_columns = valid_state_columns(truth,n_nodes) & ...
                            valid_state_columns(background,n_nodes) & ...
                            valid_state_columns(estimate,n_nodes);

    last_valid = find(runs(r).valid_columns,1,'last');
    if isempty(last_valid)
        error('%s has no time step with a sufficiently complete state.',run_def(r).key);
    end
    if ~all(runs(r).valid_columns(1:last_valid))
        error('%s contains an incomplete state inside its populated time range.', ...
              run_def(r).key);
    end
    runs(r).last_valid_time = runs(r).time(last_valid);
    fprintf('Loaded %s from %s: valid through year %.3f\n', ...
        runs(r).key,data_dir,runs(r).last_valid_time);
end

% Compare all methods at the same last time for which every run has a valid
% state. This avoids plotting a trailing, partially written HDF5 column.
final_time = min([runs.last_valid_time]);
display_final_time = round(final_time); % present the 39.9-yr output as the 40-yr endpoint
for r = 1:n_runs
    candidates = find(runs(r).valid_columns(:) & ...
                      runs(r).time(:) <= final_time+1e-10);
    final_k = candidates(end);
    runs(r).time = runs(r).time(1:final_k);
    runs(r).truth = runs(r).truth(:,1:final_k);
    runs(r).estimate = runs(r).estimate(:,1:final_k);
    runs(r).background = runs(r).background(:,1:final_k);
    runs(r).valid_columns = runs(r).valid_columns(1:final_k);
    runs(r).final_k = final_k;
    runs(r).final_time = runs(r).time(final_k);
    runs(r).friction_final = runs(r).estimate( ...
        5*n_nodes+1:6*n_nodes,final_k);
end

% Evaluate every method on exactly the same physical domain.  Basal friction
% affects grounded ice, so speed and surface RMSE use the dynamic
% true-grounded mask rather than a whole-ice or method-dependent mask.
true_grounded_mask = build_true_grounded_mask( ...
    runs(1).truth,n_nodes,density_ratio);
for r = 1:n_runs
    [runs(r).rmse_velocity, ...
     runs(r).rmse_surface_grounded, ...
     runs(r).rmse_surface_grounded_excluding_gl, ...
     runs(r).rmse_surface_whole_true_ice, ...
     runs(r).rmse_gl_km, runs(r).gl_green_x, runs(r).gl_green_y] = ...
        diagnostics(runs(r).truth, runs(r).estimate, ...
        x, y, density_ratio, n_nodes,true_grounded_mask, ...
        grounded_excluding_gl_xmax);
    switch lower(surface_rmse_domain)
        case 'grounded'
            runs(r).rmse_surface = runs(r).rmse_surface_grounded;
        case 'grounded_excluding_gl'
            runs(r).rmse_surface = ...
                runs(r).rmse_surface_grounded_excluding_gl;
        case 'whole_true_ice'
            runs(r).rmse_surface = runs(r).rmse_surface_whole_true_ice;
        otherwise
            error('Unknown surface_rmse_domain: %s',surface_rmse_domain);
    end
end

% Compact post-assimilation summaries used inside the dedicated zoom
% panels.  These expose the differences that can look deceptively small on
% the full-range RMSE axes.
velocity_zoom_summary = post_da_summary( ...
    runs,'rmse_velocity',assimilation_times, ...
    'Post-DA mean / final (m/yr)','%.1f');
surface_zoom_summary = post_da_summary( ...
    runs,'rmse_surface',assimilation_times, ...
    'Post-DA mean / final (m)','%.1f');
gl_zoom_summary = post_da_summary( ...
    runs,'rmse_gl_km',assimilation_times, ...
    'Post-DA mean / final (km)','%.2f');

for r = 1:n_runs
    coefficient_rows = 5*n_nodes+1:6*n_nodes;
    final_mask = true_grounded_mask(:,runs(r).final_k);
    runs(r).rmse_friction_final_grounded = vector_rmse( ...
        runs(r).estimate(coefficient_rows,runs(r).final_k), ...
        runs(r).truth(coefficient_rows,runs(r).final_k),final_mask);
end
velocity_zoom_summary = sprintf( ...
    '%s\nFinal grounded C RMSE\nW %.0f | E %.0f | I %.0f', ...
    velocity_zoom_summary, ...
    runs(1).rmse_friction_final_grounded, ...
    runs(2).rmse_friction_final_grounded, ...
    runs(3).rmse_friction_final_grounded);
if ~show_post_da_metric_blocks
    velocity_zoom_summary = '';
    surface_zoom_summary = '';
    gl_zoom_summary = '';
end

surface_spread_grounded = method_spread_score( ...
    runs,'rmse_surface_grounded',assimilation_times);
surface_spread_upstream = method_spread_score( ...
    runs,'rmse_surface_grounded_excluding_gl',assimilation_times);
surface_spread_whole = method_spread_score( ...
    runs,'rmse_surface_whole_true_ice',assimilation_times);
fprintf(['Median surface-RMSE method spread during assimilation (m): ' ...
         'grounded %.3f; ' ...
         'grounded excluding GL %.3f; whole true ice %.3f. ' ...
         'Plotting %s.\n'], ...
        surface_spread_grounded,surface_spread_upstream, ...
        surface_spread_whole,surface_rmse_domain);

% A filter comparison is valid only if all three experiments use the same
% truth, background realization, mesh, and synchronized time grid.
reference_nt = size(runs(1).truth,2);
for r = 2:n_runs
    if size(runs(r).truth,2) ~= reference_nt || ...
            numel(runs(r).time) ~= numel(runs(1).time) || ...
            max(abs(runs(r).time(:)-runs(1).time(:))) > 1e-10
        error('%s does not share the WBF synchronized time grid.',runs(r).key);
    end
    truth_difference = max(abs(runs(r).truth(:)-runs(1).truth(:)));
    background_difference = max(abs( ...
        runs(r).background(:)-runs(1).background(:)));
    if truth_difference > 1e-8 || background_difference > 1e-8
        error(['%s does not share the same truth/background realization ' ...
               '(max differences %.3g and %.3g).'], ...
              runs(r).key,truth_difference,background_difference);
    end
end

% WBF must retain its deliberately wrong initial ensemble coefficient field.
wbf_initial_friction = runs(1).estimate(5*n_nodes+1:6*n_nodes,1);
wbf_friction_drift = max(abs(runs(1).friction_final-wbf_initial_friction));
if wbf_friction_drift > 1e-8
    error('WBF friction was not fixed (maximum drift %.6g).',wbf_friction_drift);
end

% EBF is meaningful only if its changed friction is propagated through the
% forecast model. A changed coefficient history paired with bitwise-identical
% WBF/EBF geometry and velocity means that the EnKF coefficient update was
% saved, but never affected the subsequent ISSM forecast.
dynamic_rows = 1:4*n_nodes;  % H, S, Vx, Vy
ebf_friction_separation = max(abs( ...
    runs(2).estimate(5*n_nodes+1:6*n_nodes,:) - ...
    runs(1).estimate(5*n_nodes+1:6*n_nodes,:)),[],'all');
ebf_dynamic_separation = max(abs( ...
    runs(2).estimate(dynamic_rows,:) - ...
    runs(1).estimate(dynamic_rows,:)),[],'all');
ebf_response_tolerance = 1e-10;
method_comparison_valid = ~(ebf_friction_separation > 1e-6 && ...
                            ebf_dynamic_separation <= ebf_response_tolerance);
fprintf(['WBF/EBF method-separation check: max |dC| = %.6g, ' ...
         'max |d(H,S,Vx,Vy)| = %.6g\n'], ...
        ebf_friction_separation,ebf_dynamic_separation);
if ~method_comparison_valid
    validation_message = sprintf([ ...
        'INVALID EBF COUPLING: friction differs from WBF (max |dC| = %.3g), ' ...
        'but H, S, Vx and Vy are identical.'],ebf_friction_separation);
    if fail_on_uncoupled_ebf
        error('%s Rerun EBF after repairing coefficient-to-forecast propagation.', ...
              validation_message);
    else
        warning('%s Rendering a diagnostic preview only.',validation_message);
    end
else
    validation_message = '';
end

%% ---------------- Draw the requested stacked figure --------------------
colors = [0.0000 0.4470 0.7410; ...
          0.8500 0.3250 0.0980; ...
          0.4660 0.6740 0.1880];
styles = {'-','--','-.'};

true_friction = runs(1).truth(5*n_nodes+1:6*n_nodes,runs(1).final_k);
% map_data = {true_friction, runs(1).friction_final, ...
%             runs(2).friction_final, runs(3).friction_final};

map_data = {true_friction, runs(1).background(5*n_nodes+1:6*n_nodes,1), ...
    runs(2).friction_final, runs(3).friction_final};
true_H_final = runs(1).truth(1:n_nodes,runs(1).final_k);
true_floating_final = isfinite(true_H_final) & true_H_final > 0 & ...
    ~true_grounded_mask(:,runs(1).final_k);

% Plot the stored coefficient across the complete mesh. Do not replace
% floating-ice values with zero/NaN or imprint a method-dependent GL mask.
map_titles = { ...
    sprintf('True basal friction at t = %g years',display_final_time), ...
    'WBF: fixed wrong basal friction', ...
    sprintf('EBF: EnKF-only recovery at t = %g years',display_final_time), ...
    sprintf('IBF: inversion recovery at t = %g years',display_final_time)};
panel_labels = {'(a)','(b)','(c)','(d)'};

finite_friction = map_data{1}(isfinite(map_data{1}));
if isempty(finite_friction)
    error('All final true basal-friction values are non-finite.');
end
% Use the stored true-friction range for every panel. Recovered outliers
% saturate instead of stretching the scale and hiding the true structure.
friction_limits = [min(finite_friction), max(finite_friction)];
if friction_limits(1) == friction_limits(2)
    friction_limits = friction_limits + [-1 1];
end

% Use comparable physical heights for all seven panels. The map stack is
% compact, while the main RMSE panels remain wider than their zoom columns.
fig = figure('Color','w','Visible','off','Units','pixels', ...
    'Position',[80 40 2000 2100]);

map_left = 0.09; map_width = 0.76; map_height = 0.095;
map_bottom = [0.875 0.765 0.655 0.545];
map_axes = gobjects(4,1);
x_limits = [min(x) max(x)]/1000;
y_limits = [min(y) max(y)]/1000;

for p = 1:4
    ax = axes(fig,'Position',[map_left map_bottom(p) map_width map_height]);
    map_axes(p) = ax;
    plot_basal_friction(ax,elements,x,y,map_data{p}, ...
        x_limits,y_limits,friction_limits);
    % if p == 4 
    if p ~= 1
        % Basal friction is not constrained beneath floating ice. Mask that
        % part of the IBF estimate so its bound-limited values are not
        % interpreted as recovered basal properties.
        overlay_masked_region(ax,elements,x,y,true_floating_final,[0.72 0.72 0.72]);
    end
    title(ax,map_titles{p},'FontWeight','bold','FontSize',12);
    add_panel_label(ax,panel_labels{p},false);
    ylabel(ax,'y (km)','FontWeight','bold','FontSize',12);
    if p < 4
        ax.XTickLabel = [];
    else
        xlabel(ax,'x (km)','FontWeight','bold','FontSize',12);
    end
end

% One colorbar for all four maps, as in the shared evolution figures.
middle_map_position = [map_left map_bottom(2) map_width map_height];
cb = colorbar(map_axes(2));
cb.Position = [0.875 0.545 0.017 0.425];
map_axes(2).Position = middle_map_position;
cb.FontSize = 11;
cb.FontWeight = 'bold';
cb.LineWidth = 1.2;
ylabel(cb,'Basal friction (Pa m^{-1/3} yr^{-1/3})', ...
    'FontSize',13,'FontWeight','bold');

if grounded_rmse_zoom_column
    metric_main_width = 0.67;
else
    metric_main_width = 0.89;
end

% Compact the vertical stack with small, uniform gutters between RMSE panels.
metric_vel_bottom = 0.390;
metric_surf_bottom = 0.260;
metric_gl_bottom = 0.130;
metric_height = 0.105;

ax_vel = axes(fig,'Position',[0.09 metric_vel_bottom metric_main_width metric_height]);
plot_metric(ax_vel,runs,'rmse_velocity',colors,styles, ...
    assimilation_times, ...
    'Velocity', ...
    'RMSE (m yr^{-1})',false);
% ylim(ax_vel,[0 1400]);
ylim(ax_vel,[0 140]);
add_panel_label(ax_vel,'(e)',false);

ax_surf = axes(fig,'Position',[0.09 metric_surf_bottom metric_main_width metric_height]);
plot_metric(ax_surf,runs,'rmse_surface',colors,styles, ...
    assimilation_times, ...
    'Surface elevation', ...
    'RMSE (m)',false);
% ylim(ax_surf,[0 300]);
ylim(ax_surf,[0 100]);
add_panel_label(ax_surf,'(f)',false);

zoom_axes = gobjects(0);
if grounded_rmse_zoom_column
    % Retain the publication-friendly gutter between the full-range and zoom
    % axes; vertical rather than horizontal whitespace is compacted.
    ax_vel_zoom = axes(fig,'Position',[0.820 metric_vel_bottom 0.140 metric_height]);
    plot_metric_zoom(ax_vel_zoom,runs,'rmse_velocity',colors,styles, ...
        assimilation_times,velocity_rmse_zoom_limits, ...
        velocity_rmse_zoom_ticks,'Zoomed RMSE velocity',false, ...
        velocity_zoom_summary);

    ax_surf_zoom = axes(fig,'Position',[0.820 metric_surf_bottom 0.140 metric_height]);
    plot_metric_zoom(ax_surf_zoom,runs,'rmse_surface',colors,styles, ...
        assimilation_times,surface_rmse_zoom_limits, ...
        surface_rmse_zoom_ticks,'Zoomed RMSE Surface',true, ...
        surface_zoom_summary);
end

if grounded_rmse_zoom_column && gl_rmse_zoom_panel
    gl_main_width = metric_main_width;
else
    gl_main_width = 0.89;
end
ax_gl = axes(fig,'Position',[0.09 metric_gl_bottom gl_main_width metric_height]);
plot_metric(ax_gl,runs,'rmse_gl_km',colors,styles, ...
    assimilation_times,'Absolute centerline grounding-line error', ...
    '|\Delta x_{GL}| (km)',true);
ylim(ax_gl,[0 11]);
add_panel_label(ax_gl,'(g)',false);

if grounded_rmse_zoom_column && gl_rmse_zoom_panel
    ax_gl_zoom = axes(fig,'Position',[0.820 metric_gl_bottom 0.140 metric_height]);
    plot_metric_zoom(ax_gl_zoom,runs,'rmse_gl_km',colors,styles, ...
        assimilation_times,gl_rmse_zoom_limits,gl_rmse_zoom_ticks, ...
        '1-4 km zoomed',true,gl_zoom_summary);
    zoom_axes = [ax_vel_zoom; ax_surf_zoom; ax_gl_zoom];
elseif grounded_rmse_zoom_column
    zoom_axes = [ax_vel_zoom; ax_surf_zoom];
end

legend(ax_vel,{runs.key},'Location','best','Orientation','horizontal', ...
    'Box','off','FontWeight','bold','FontSize',11);

set(map_axes, ...
    'FontName','Helvetica','FontSize',11,'FontWeight','bold', ...
    'LineWidth',1.5,'Box','on','TickDir','out','Layer','top', ...
    'TickLength',[0.004 0.004]);
set([ax_vel; ax_surf; ax_gl], ...
    'FontName','Helvetica','FontSize',13,'FontWeight','bold', ...
    'LineWidth',1.7,'Box','on','TickDir','out','Layer','top', ...
    'TickLength',[0.004 0.004]);
if ~isempty(zoom_axes)
    set(zoom_axes, ...
        'FontName','Helvetica','FontSize',10,'FontWeight','bold', ...
        'LineWidth',1.4,'Box','on','TickDir','out','Layer','top', ...
        'TickLength',[0.006 0.006]);
end

if ~method_comparison_valid
    annotation(fig,'textbox',[0.09 0.965 0.84 0.022], ...
        'String',validation_message, ...
        'Color',[0.75 0 0],'EdgeColor',[0.75 0 0], ...
        'BackgroundColor',[1.0 0.94 0.94], ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Helvetica','FontSize',11,'FontWeight','bold', ...
        'FitBoxToText','off');
end

png_file = fullfile(figure_dir,'wbf_ebf_ibf_comparison.png');
pdf_file = fullfile(figure_dir,'wbf_ebf_ibf_comparison.pdf');

% Save diagnostics before the comparatively expensive figure export so the
% numerical audit remains available even if a graphics backend is interrupted.
metrics_file = fullfile(figure_dir,'wbf_ebf_ibf_metrics.mat');
metrics_runs = rmfield(runs,{'truth','estimate','background','valid_columns'});
save(metrics_file,'metrics_runs','run_def','assimilation_times', ...
    'dt_fallback','final_time','method_comparison_valid', ...
    'ebf_friction_separation','ebf_dynamic_separation', ...
    'surface_rmse_domain','surface_spread_grounded', ...
    'surface_spread_upstream','surface_spread_whole');

% Use one explicit portrait page for both formats.  On macOS, exportgraphics
% may otherwise crop a large invisible figure to a landscape screen-sized
% bounding box, producing a PDF with a different aspect ratio from the
% intended publication layout.
set(fig,'PaperUnits','inches','PaperSize',[10 12], ...
    'PaperPosition',[0 0 10 12],'PaperPositionMode','manual');
print(fig,png_file,'-dpng','-r300');
print(fig,pdf_file,'-dpdf','-painters');

fprintf('Saved:\n  %s\n  %s\n  %s\n',png_file,pdf_file,metrics_file);

%% =======================================================================
%% Local functions

function path = first_existing_file(candidates)
    path = '';
    for i = 1:numel(candidates)
        if isfile(candidates{i})
            path = candidates{i};
            return
        end
    end
end

function md = read_issm_model(model_file)
    if exist('loadmodel','file') == 2
        md = loadmodel(model_file);
        return
    end
    loaded = load(model_file);
    if isfield(loaded,'md')
        md = loaded.md;
        return
    end
    names = fieldnames(loaded);
    for i = 1:numel(names)
        candidate = loaded.(names{i});
        if (isstruct(candidate) && isfield(candidate,'mesh')) || ...
           (isobject(candidate) && isprop(candidate,'mesh'))
            md = candidate;
            return
        end
    end
    error('No ISSM model object with a mesh was found in %s.',model_file);
end

function state = read_state_matrix(file_name,dataset,n_nodes,nvar)
    state = squeeze(h5read(file_name,dataset));
    expected = n_nodes*nvar;
    if size(state,1) == expected
        % already state-by-time
    elseif size(state,2) == expected
        state = state.';
    else
        error(['%s:%s has size %s; one dimension must equal %d ' ...
               '(6 variables x %d mesh nodes).'], ...
            file_name,dataset,mat2str(size(state)),expected,n_nodes);
    end
    state = double(state);
end

function valid = valid_state_columns(state,n)
% Reject trailing HDF5 columns that are incomplete or left at zero.
    finite_fraction = sum(isfinite(state),1)/size(state,1);
    H = state(1:n,:);
    C = state(5*n+1:6*n,:);
    ice_fraction = sum(isfinite(H) & H > 0,1)/n;
    friction_fraction = sum(isfinite(C) & C > 0,1)/n;
    valid = finite_fraction >= 0.98 & ice_fraction >= 0.05 & ...
            friction_fraction >= 0.05;
    valid = valid(:);
end

function time = read_time_vector(root_dir,data_dir,nt,dt)
    candidates = {fullfile(data_dir,'true-wrong-issm.h5')};
    time = [];
    for i = 1:numel(candidates)
        if ~isfile(candidates{i}), continue, end
        try
            time = double(h5read(candidates{i},'/t'));
            time = time(:);
            break
        catch
        end
    end
    if numel(time) < nt
        time = (0:nt-1)'*dt;
    else
        time = time(1:nt);
    end
end

function grounded = build_true_grounded_mask(truth,n,density_ratio)
% One common, time-dependent evaluation mask for WBF, EBF, and IBF.
    H = truth(1:n,:);
    bed = truth(4*n+1:5*n,:);
    flotation = H + bed/density_ratio;
    grounded = isfinite(H) & isfinite(bed) & H > 0 & flotation > 0;
end

function [rmse_velocity,rmse_surface_grounded, ...
          rmse_surface_grounded_excluding_gl, ...
          rmse_surface_whole_true_ice, ...
          rmse_gl_km,gl_green_x,gl_green_y] = diagnostics( ...
        truth,estimate,x,y,density_ratio,n,true_grounded_mask, ...
        grounded_excluding_gl_xmax)
    nt = size(truth,2);
    I_H  = 1:n;
    I_S  = n+1:2*n;
    I_Vx = 2*n+1:3*n;
    I_Vy = 3*n+1:4*n;
    I_b  = 4*n+1:5*n;

    rmse_velocity = nan(nt,1);
    rmse_surface_grounded = nan(nt,1);
    rmse_surface_grounded_excluding_gl = nan(nt,1);
    rmse_surface_whole_true_ice = nan(nt,1);
    rmse_gl_km = nan(nt,1);
    gl_green_x = nan(nt,1);
    gl_green_y = nan(nt,1);

    y_center = 0.5*(min(y)+max(y));
    xg = linspace(min(x),max(x),420);
    yg = linspace(min(y),max(y),70);
    [Xg,Yg] = meshgrid(xg,yg);
    x_previous_true = NaN;
    x_previous_estimate = NaN;

    if ~isequal(size(true_grounded_mask),[n nt])
        error('The common true-grounded mask must have size %d-by-%d.',n,nt);
    end

    for k = 1:nt
        Ht = truth(I_H,k);
        St = truth(I_S,k);
        Ut = hypot(truth(I_Vx,k),truth(I_Vy,k));
        Ue = hypot(estimate(I_Vx,k),estimate(I_Vy,k));
        Se = estimate(I_S,k);

        phi_t = Ht + truth(I_b,k)/density_ratio;
        phi_e = estimate(I_H,k) + estimate(I_b,k)/density_ratio;

        % Use the same dynamic true-grounded vertices for every method.
        mask = true_grounded_mask(:,k) & isfinite(Ut) & isfinite(Ue);
        rmse_velocity(k) = vector_rmse(Ue,Ut,mask);

        mask_surface_grounded = true_grounded_mask(:,k) & ...
            isfinite(St) & isfinite(Se);
        mask_surface_upstream = mask_surface_grounded & ...
            x <= grounded_excluding_gl_xmax;
        mask_surface_whole_ice = Ht > 0 & isfinite(St) & isfinite(Se);
        rmse_surface_grounded(k) = ...
            vector_rmse(Se,St,mask_surface_grounded);
        rmse_surface_grounded_excluding_gl(k) = ...
            vector_rmse(Se,St,mask_surface_upstream);
        rmse_surface_whole_true_ice(k) = ...
            vector_rmse(Se,St,mask_surface_whole_ice);

        % Measure the same quantity represented by the centerline GL marker:
        % the absolute x-distance between the true and estimated zero-level
        % crossings at y = y_center. Track each branch through time so that
        % a changing contour topology cannot silently switch GL branches.
        Phi_t = interpolate_levelset(x,y,phi_t,Xg,Yg);
        Phi_e = interpolate_levelset(x,y,phi_e,Xg,Yg);
        [x_true,y_true] = tracked_centerline_gl( ...
            Xg,Yg,Phi_t,y_center,x_previous_true);
        [x_estimate,~] = tracked_centerline_gl( ...
            Xg,Yg,Phi_e,y_center,x_previous_estimate);
        gl_green_x(k) = x_true;
        gl_green_y(k) = y_true;
        if isfinite(x_true)
            x_previous_true = x_true;
        end
        if isfinite(x_estimate)
            x_previous_estimate = x_estimate;
        end
        if isfinite(x_true) && isfinite(x_estimate)
            rmse_gl_km(k) = abs(x_estimate-x_true)/1000;
        end
    end
end

function value = vector_rmse(a,b,mask)
    mask = logical(mask(:)) & isfinite(a(:)) & isfinite(b(:));
    if ~any(mask)
        value = NaN;
    else
        difference = a(mask)-b(mask);
        value = sqrt(mean(difference.^2));
    end
end

function score = method_spread_score(runs,field_name,assim_times)
% Median cross-method RMSE spread over the assimilation interval.
    nt = numel(runs(1).(field_name));
    curves = nan(numel(runs),nt);
    for r = 1:numel(runs)
        values = runs(r).(field_name);
        curves(r,:) = values(:).';
    end
    spread = max(curves,[],1,'omitnan') - min(curves,[],1,'omitnan');
    time = runs(1).time(:).';
    if isempty(assim_times)
        evaluation_window = true(size(time));
    else
        evaluation_window = time >= min(assim_times) & ...
                            time <= max(assim_times);
    end
    finite_spread = spread(evaluation_window & isfinite(spread));
    if isempty(finite_spread)
        score = NaN;
    else
        score = median(finite_spread);
    end
end

function Phi = interpolate_levelset(x,y,phi,Xg,Yg)
    good = isfinite(x) & isfinite(y) & isfinite(phi);
    if nnz(good) < 3
        Phi = nan(size(Xg));
        return
    end
    F = scatteredInterpolant(x(good),y(good),phi(good),'linear','nearest');
    Phi = F(Xg,Yg);
end

function [xc,yc] = tracked_centerline_gl(Xg,Yg,Phi,y_center,x_previous)
% Centerline crossing of phi=0, selected continuously through time.
    xc = NaN;
    yc = NaN;
    if ~any(isfinite(Phi(:))), return, end
    Phi(~isfinite(Phi)) = 1;
    if min(Phi(:))*max(Phi(:)) > 0, return, end
    C = contourc(Xg(1,:),Yg(:,1),Phi,[0 0]);
    [segments,lengths] = unpack_contours(C,0);
    if isempty(segments), return, end

    x_candidates = [];
    for j = 1:numel(segments)
        points = segments{j};
        gx = points(1,:);
        gy = points(2,:);
        crossing = find((gy(1:end-1)-y_center).* ...
                        (gy(2:end)-y_center) <= 0);
        for i = crossing(:).'
            if abs(gy(i+1)-gy(i)) < eps
                alpha = 0.5;
            else
                alpha = (y_center-gy(i))/(gy(i+1)-gy(i));
            end
            candidate = gx(i)+alpha*(gx(i+1)-gx(i));
            if isfinite(candidate)
                x_candidates(end+1,1) = candidate; %#ok<AGROW>
            end
        end
    end
    if isempty(x_candidates), return, end
    if isfinite(x_previous)
        x_target = x_previous;
    else
        % Seed from the main (longest) contour, exactly as the spatial
        % grounding-line plots do, before switching to temporal tracking.
        [~,main_index] = max(lengths);
        main_points = segments{main_index};
        gx = main_points(1,:);
        gy = main_points(2,:);
        crossing = find((gy(1:end-1)-y_center).* ...
                        (gy(2:end)-y_center) <= 0,1);
        if isempty(crossing)
            [~,nearest] = min(abs(gy-y_center));
            x_target = gx(nearest);
        elseif abs(gy(crossing+1)-gy(crossing)) < eps
            x_target = 0.5*(gx(crossing)+gx(crossing+1));
        else
            alpha = (y_center-gy(crossing))/ ...
                    (gy(crossing+1)-gy(crossing));
            x_target = gx(crossing)+alpha* ...
                       (gx(crossing+1)-gx(crossing));
        end
    end
    [~,idx] = min(abs(x_candidates-x_target));
    xc = x_candidates(idx);
    yc = y_center;
end

function [segments,lengths] = unpack_contours(C,min_length)
    segments = {};
    lengths = [];
    k = 1;
    while k < size(C,2)
        npts = C(2,k);
        points = C(:,k+1:k+npts);
        k = k+npts+1;
        contour_length = sum(hypot(diff(points(1,:)),diff(points(2,:))));
        if contour_length >= min_length
            segments{end+1} = points; %#ok<AGROW>
            lengths(end+1) = contour_length; %#ok<AGROW>
        end
    end
end

function plot_basal_friction(ax,elements,x,y,friction,x_limits,y_limits,c_limits)
% Match the smooth nodal rendering used by ISSM plotmodel.
    patch(ax,'Faces',elements,'Vertices',[x(:)/1000 y(:)/1000], ...
        'FaceVertexCData',friction(:),'FaceColor','interp', ...
        'EdgeColor','none');
    view(ax,2);
    xlim(ax,x_limits);
    ylim(ax,y_limits);
    daspect(ax,[1 1 1]);
    colormap(ax,parula(256));
    caxis(ax,c_limits);

    xt = ceil(x_limits(1)/100)*100 : 100 : floor(x_limits(2)/100)*100;
    yt = ceil(y_limits(1)/40)*40 : 40 : floor(y_limits(2)/40)*40;
    if ~isempty(xt), ax.XTick = xt; end
    if ~isempty(yt), ax.YTick = yt; end
end

function overlay_masked_region(ax,elements,x,y,node_mask,face_color)
% Cover triangles lying predominantly in an unconstrained domain. A
% majority-node rule keeps the gray mask aligned with the mesh while
% avoiding a visible incursion into the grounded side of the GL.
    node_mask = logical(node_mask(:));
    masked_faces = sum(node_mask(elements),2) >= 2;
    if any(masked_faces)
        hold(ax,'on');
        patch(ax,'Faces',elements(masked_faces,:), ...
            'Vertices',[x(:)/1000 y(:)/1000], ...
            'FaceColor',face_color,'EdgeColor','none');
        hold(ax,'off');
    end
end

function add_panel_label(ax,label_text,use_white_backing)
% Put panel identifiers inside the axes, consistently with the other ISSM
% evolution figures.
    args = {'Units','normalized', ...
            'VerticalAlignment','top', ...
            'HorizontalAlignment','left', ...
            'FontName','Helvetica', ...
            'FontSize',13, ...
            'FontWeight','bold', ...
            'Color',[0 0 0], ...
            'Clipping','on'};
    if use_white_backing
        args = [args, {'BackgroundColor',[1 1 1], ...
                       'EdgeColor','none','Margin',1}];
    end
    text(ax,0.018,0.94,label_text,args{:});
end

function plot_metric(ax,runs,field_name,colors,styles,assim_times, ...
        title_text,ylabel_text,show_xlabel)
    hold(ax,'on');
    for r = 1:numel(runs)
        plot(ax,runs(r).time,runs(r).(field_name), ...
            'Color',colors(r,:),'LineStyle',styles{r},'LineWidth',2.2);
    end
    if ~isempty(assim_times)
        line_handle = xline(ax,max(assim_times),'--', ...
            'Color',[0.20 0.20 0.20],'LineWidth',1.4);
        line_handle.HandleVisibility = 'off';
    end
    hold(ax,'off');
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    ax.GridAlpha = 0.18;
    title(ax,title_text,'FontWeight','bold','FontSize',13);
    ylabel(ax,ylabel_text,'FontWeight','bold','FontSize',14);
    if show_xlabel
        xlabel(ax,'Time (years)','FontWeight','bold','FontSize',14);
    else
        ax.XTickLabel = [];
    end
end

function plot_metric_zoom(ax,runs,field_name,colors,styles,assim_times, ...
        y_limits,y_ticks,title_text,show_xlabel,summary_text)
% Dedicated zoom column: expose low RMSE without obscuring the main panel.
    hold(ax,'on');
    for r = 1:numel(runs)
        plot(ax,runs(r).time,runs(r).(field_name), ...
            'Color',colors(r,:),'LineStyle',styles{r},'LineWidth',1.8);
    end
    if ~isempty(assim_times)
        line_handle = xline(ax,max(assim_times),'--', ...
            'Color',[0.20 0.20 0.20],'LineWidth',1.0);
        line_handle.HandleVisibility = 'off';
    end
    hold(ax,'off');

    x_start = floor(runs(1).time(1));
    x_end = ceil(runs(1).time(end));
    xlim(ax,[x_start x_end]);
    tick_start = ceil(x_start/10)*10;
    tick_end = floor(x_end/10)*10;
    ax.XTick = tick_start:10:tick_end;
    ylim(ax,y_limits);
    ax.YTick = y_ticks;
    if show_xlabel
        xlabel(ax,'Time (years)','FontWeight','bold','FontSize',10);
    end
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    ax.GridAlpha = 0.18;
    title(ax,title_text,'FontWeight','bold','FontSize',10);
    if nargin >= 11 && ~isempty(summary_text)
        text(ax,0.025,0.965,summary_text, ...
            'Units','normalized','VerticalAlignment','top', ...
            'HorizontalAlignment','left','Interpreter','none', ...
            'FontName','Helvetica','FontSize',7.2,'FontWeight','bold', ...
            'Color',[0.08 0.08 0.08], ...
            'BackgroundColor',[1 1 1], ...
            'EdgeColor',[0.72 0.72 0.72],'Margin',2);
    end
end

function summary_text = post_da_summary( ...
        runs,field_name,assim_times,heading,value_format)
% A compact WBF/EBF/IBF mean/final summary for a zoom panel.
    if isempty(assim_times)
        start_time = runs(1).time(1);
    else
        start_time = max(assim_times);
    end
    rows = cell(numel(runs)+1,1);
    rows{1} = heading;
    for r = 1:numel(runs)
        values = runs(r).(field_name)(:);
        time = runs(r).time(:);
        window = time >= start_time & isfinite(values);
        if any(window)
            mean_value = mean(values(window),'omitnan');
        else
            mean_value = NaN;
        end
        final_index = find(isfinite(values),1,'last');
        if isempty(final_index)
            final_value = NaN;
        else
            final_value = values(final_index);
        end
        pair_format = ['%s ' value_format ' / ' value_format];
        rows{r+1} = sprintf(pair_format,runs(r).key,mean_value,final_value);
    end
    summary_text = strjoin(rows,newline);
end
