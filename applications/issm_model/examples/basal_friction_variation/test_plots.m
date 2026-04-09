close all; clearvars; clear all

global data_file_paths nvar ensemble_vec_full ...
        label_t t nt colorbar_gap bed_obs_xy t_change ...
        t_jump k_len

nvar = 6;
colorbar_gap = 0.92;
t_change =1000;
yr_jump  = 2;

k_array = [30, 70,100, 120, 180, (t_change-1)]+1;
k_len = length(k_array);
% k_array = [20, 80, 120, 240, 480, 560, 700] +1;
dt      = 0.2;
t_jump = (yr_jump+1)/dt;

% ---------------- Load essentials --------------
results_dir = 'results';
filter_type = 'true-wrong';
file_path   = fullfile(results_dir, sprintf('%s-issm.h5', filter_type));
t        = h5read(file_path,'/t'); 
ind_m    = h5read(file_path,'/obs_index'); 
tm_m     = h5read(file_path,'/obs_max_time'); 
run_mode = h5read(file_path,'/run_mode'); 

% --------- true / wrong (nurged)
% data_file_paths = '_modelrun_datasets_0';
data_file_paths = '_modelrun_datasets';
file_path          = fullfile(data_file_paths, 'true_nurged_states.h5');
model_true_state   = h5read(file_path,'/true_state')';
model_nurged_state = h5read(file_path,'/nurged_state')';
[nd, nt] = size(model_true_state );

% obs (kept)
file_path = fullfile(data_file_paths, 'synthetic_obs.h5');
w = h5read(file_path, '/hu_obs')'; 

% ----- ensemble mean
file_path         = fullfile(data_file_paths, 'icesee_ensemble_data.h5');
ensemble_vec_mean = h5read(file_path, '/ensemble_mean')';
ensemble_vec_full = h5read(file_path, '/ensemble'); 

% ISSM model template
md = loadmodel(fullfile("data","ISMIP.Parameterization1.mat"));
md_true   = md;
md_nurged = md;
md_ens    = md;

% ----- Bed observation XY (for overlay on True Bed) 
hdim = nd / nvar;

wbed = w(4*hdim + 1 : 5*hdim, 1);

obs_idx = find(~isnan(wbed));   % <-- observations are the non-NaNs
bed_obs_xy = [md_true.mesh.x(obs_idx), md_true.mesh.y(obs_idx)];

% ---------------- GL midpoint points  ----------------
[gl_mid] = compute_gl_midpoints( ...
    k_array, dt, ...
    model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md);

% Example requested times
times_to_plot = [1.0, 6.2, 12.2, 18.2, 22.2, 30.0];

k_array = zeros(size(times_to_plot));
for i = 1:length(times_to_plot)
    [~, k_array(i)] = min(abs(t - times_to_plot(i)));
end

plot_true_friction_snapshots(k_array, model_true_state, md_true, md, '');
plot_true_friction_anomaly_snapshots(k_array, model_true_state, md_true, md, '');

scalar_means_file_true = fullfile(data_file_paths,'ensemble_true_state_scalar_0.h5');
scalar_means_file_nurged = fullfile(data_file_paths,'ensemble_nurged_state_scalar_0.h5');
scalar_means_file_ens = fullfile(data_file_paths,'ensemble_scalar_output.h5');
% scalar_means_file_ens = fullfile('_modelrun_datasets','ensemble_out_scalar_10.h5');

plot_gl_on_bed_evolution( ...
    k_array, dt, ...
    model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md, ...
    'friction.coefficient', 'Basal Friction', 'Pa m^{-1/3} yr^{-1/3}', gl_mid, ...
    scalar_means_file_true, scalar_means_file_nurged, scalar_means_file_ens );

% ============================================================
% Raw friction snapshots
% ============================================================
function plot_true_friction_snapshots(k_array, model_true_state, md_true, md, units)

    global t colorbar_gap

    if nargin < 5
        units = '';
    end

    if ~isempty(units)
        units_str = [' (' units ')'];
    else
        units_str = '';
    end

    field_title = 'Friction coefficient';

    nk = length(k_array);
    nrows = nk;

    figure('Position',[150 150 1800 1000]); clf;
    axs = gobjects(nrows,1);

    % -------------------------------------------------------------
    % Global color limits across all requested snapshots
    % -------------------------------------------------------------
    all_data = [];
    for i = 1:nk
        k = k_array(i);
        md_true_k = setup_true_state_only(k, model_true_state, md_true, md);
        tmp = md_true_k.friction.coefficient;
        all_data = [all_data; tmp(:)];
    end
    cmin = min(all_data);
    cmax = max(all_data);
    clear all_data

    % -------------------------------------------------------------
    % Plot each requested snapshot
    % -------------------------------------------------------------
    for i = 1:nk
        k = k_array(i);

        md_true_k = setup_true_state_only(k, model_true_state, md_true, md);
        fc = md_true_k.friction.coefficient;

        label_t = t(k);

        plotmodel(md_true_k, ...
            'data', fc, ...
            'title', sprintf('True friction at t = %.1f yr', label_t), ...
            'subplot', [nrows, 1, i], ...
            'caxis', [cmin cmax], ...
            'colorbar', 'off');

        ax = gca;
        axs(i) = ax;

        ttl = ax.Title;
        ttl.FontSize = 12;
        ttl.FontWeight = 'bold';
        ttl.Interpreter = 'tex';

        % km tick labels
        xt = get(ax,'XTick');
        yt = get(ax,'YTick');
        set(ax,'XTickLabel', arrayfun(@(v) sprintf('%g', v./1000), xt, 'UniformOutput', false));
        set(ax,'YTickLabel', arrayfun(@(v) sprintf('%g', v./1000), yt, 'UniformOutput', false));

        ylabel(ax,'y (km)','FontWeight','bold','FontSize',14);

        if i < nrows
            set(ax,'XTickLabel',[]);
        else
            xlabel(ax,'x (km)','FontWeight','bold','FontSize',14);
        end

        panel = sprintf('(%c)', 'a'+(i-1));
        text(ax, 0.02, 0.95, panel, 'Units','normalized', ...
            'FontWeight','bold', 'FontSize', 14, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color','k');
    end

    % -------------------------------------------------------------
    % Layout cleanup
    % -------------------------------------------------------------
    axs = flipud(findall(gcf,'Type','axes'));
    gap = 0.03; top = 0.96; bottom = 0.08;
    avail = top-bottom - (nrows-1)*gap;
    height = avail/nrows;

    for i = 1:nrows
        ax = axs(i);
        pos = [0.10, bottom+(nrows-i)*(height+gap), 0.80, height];
        set(ax, 'Position', pos, ...
            'FontWeight','bold', ...
            'FontSize',14, ...
            'Box','on', ...
            'LineWidth',1.8, ...
            'TickDir','out', ...
            'TickLength',[0.004 0.004], ...
            'XGrid','off', ...
            'YGrid','off', ...
            'YMinorGrid','off', ...
            'XTickMode','manual', ...
            'YTickMode','manual');

        yl = ax.YLim / 1000;
        step = 40;
        yt_km = ceil(yl(1)/step)*step : step : floor(yl(2)/step)*step;
        set(ax,'YTick',yt_km*1000,'YTickLabel',string(yt_km),'YTickMode','manual');

        xl = ax.XLim / 1000;
        xt_km = floor(xl(1)/100)*100 : 100 : ceil(xl(2)/100)*100;
        set(ax,'XTick',xt_km*1000, ...
               'XTickLabel',arrayfun(@num2str, xt_km, 'UniformOutput', false), ...
               'XTickMode','manual');

        ylabel(ax,'y (km)','FontWeight','bold','FontSize',14);

        if i < nrows
            set(ax,'XTickLabel',[]);
        else
            xlabel(ax,'x (km)','FontWeight','bold','FontSize',14);
        end
    end

    % shared colorbar
    for i = 1:nrows
        colormap(axs(i), parula);
    end
    cb = colorbar(axs(end), 'Position',[colorbar_gap 0.25 0.015 0.45]);
    ylabel(cb, [field_title units_str], 'FontSize', 14, 'FontWeight', 'bold');

    set(gcf,'Color','w');

    % save_figure_300dpi('true_friction_snapshots');
end

% ============================================================
% Friction anomaly snapshots: Delta C = C(t) - C(t0)
% ============================================================
function plot_true_friction_anomaly_snapshots(k_array, model_true_state, md_true, md, units)

    global t colorbar_gap

    if nargin < 5
        units = '';
    end

    if ~isempty(units)
        units_str = [' (' units ')'];
    else
        units_str = '';
    end

    field_title = 'Friction anomaly';

    nk = length(k_array);
    nrows = nk;

    figure('Position',[150 150 1800 1000]); clf;
    axs = gobjects(nrows,1);

    % ---------------------------
    % Baseline friction at t0
    % ---------------------------
    md_true_0 = setup_true_state_only(1, model_true_state, md_true, md);
    fc0 = md_true_0.friction.coefficient;

    % ---------------------------
    % Global symmetric color limits
    % ---------------------------
    dfc_all = [];
    for i = 1:nk
        k = k_array(i);
        md_true_k = setup_true_state_only(k, model_true_state, md_true, md);
        dfc = md_true_k.friction.coefficient - fc0;
        dfc_all = [dfc_all; dfc(:)];
    end
    amax = max(abs(dfc_all));
    cmin = -amax;
    cmax =  amax;

    % ---------------------------
    % Plot anomaly snapshots
    % ---------------------------
    for i = 1:nk
        k = k_array(i);

        md_true_k = setup_true_state_only(k, model_true_state, md_true, md);
        dfc = md_true_k.friction.coefficient - fc0;

        label_t = t(k);

        plotmodel(md_true_k, ...
            'data', dfc, ...
            'title', sprintf('\\Delta friction at t = %.1f yr', label_t), ...
            'subplot', [nrows, 1, i], ...
            'caxis', [cmin cmax], ...
            'colorbar', 'off');

        ax = gca;
        axs(i) = ax;

        ttl = ax.Title;
        ttl.FontSize = 12;
        ttl.FontWeight = 'bold';
        ttl.Interpreter = 'tex';

        xt = get(ax,'XTick');
        yt = get(ax,'YTick');
        set(ax,'XTickLabel', arrayfun(@(v) sprintf('%g', v./1000), xt, 'UniformOutput', false));
        set(ax,'YTickLabel', arrayfun(@(v) sprintf('%g', v./1000), yt, 'UniformOutput', false));

        ylabel(ax,'y (km)','FontWeight','bold','FontSize',14);

        if i < nrows
            set(ax,'XTickLabel',[]);
        else
            xlabel(ax,'x (km)','FontWeight','bold','FontSize',14);
        end

        panel = sprintf('(%c)', 'a'+(i-1));
        text(ax, 0.02, 0.95, panel, 'Units','normalized', ...
            'FontWeight','bold', 'FontSize', 14, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color','k');
    end

    % ---------------------------
    % Layout cleanup
    % ---------------------------
    axs = flipud(findall(gcf,'Type','axes'));
    gap = 0.03; top = 0.96; bottom = 0.08;
    avail = top-bottom - (nrows-1)*gap;
    height = avail/nrows;

    for i = 1:nrows
        ax = axs(i);
        pos = [0.10, bottom+(nrows-i)*(height+gap), 0.80, height];
        set(ax, 'Position', pos, ...
            'FontWeight','bold', ...
            'FontSize',14, ...
            'Box','on', ...
            'LineWidth',1.8, ...
            'TickDir','out', ...
            'TickLength',[0.004 0.004], ...
            'XGrid','off', ...
            'YGrid','off', ...
            'YMinorGrid','off', ...
            'XTickMode','manual', ...
            'YTickMode','manual');

        yl = ax.YLim / 1000;
        ystep = 40;
        yt_km = ceil(yl(1)/ystep)*ystep : ystep : floor(yl(2)/ystep)*ystep;
        set(ax,'YTick',yt_km*1000,'YTickLabel',string(yt_km),'YTickMode','manual');

        xl = ax.XLim / 1000;
        xt_km = floor(xl(1)/100)*100 : 100 : ceil(xl(2)/100)*100;
        set(ax,'XTick',xt_km*1000, ...
               'XTickLabel',arrayfun(@num2str, xt_km, 'UniformOutput', false), ...
               'XTickMode','manual');

        ylabel(ax,'y (km)','FontWeight','bold','FontSize',14);

        if i < nrows
            set(ax,'XTickLabel',[]);
        else
            xlabel(ax,'x (km)','FontWeight','bold','FontSize',14);
        end
    end

    % ---------------------------
    % Shared colorbar and colormap
    % ---------------------------
    for i = 1:nrows
        colormap(axs(i), redblue(256));
        caxis(axs(i), [cmin cmax]);
    end

    cb = colorbar(axs(end), 'Position',[colorbar_gap 0.25 0.015 0.45]);
    ylabel(cb, ['\Delta friction' units_str], 'FontSize', 14, 'FontWeight', 'bold');

    set(gcf,'Color','w');

    % save_figure_300dpi('true_friction_anomaly_snapshots');
end

% ============================================================
% Reconstruct ISSM true state at time index k
% ============================================================
function md_true_k = setup_true_state_only(k, model_true_state, md_true, md)

    nvar = 6;
    ndim = size(model_true_state, 1);
    hdim = ndim / nvar;
    di   = md.materials.rho_ice / md.materials.rho_water;

    H   = model_true_state(1:hdim, k);
    S   = model_true_state(hdim+1:2*hdim, k);
    B   = S - H;
    Vx  = model_true_state(2*hdim+1:3*hdim, k);
    Vy  = model_true_state(3*hdim+1:4*hdim, k);
    Vel = hypot(Vx, Vy);
    bed = model_true_state(4*hdim+1:5*hdim, k);
    fc  = model_true_state(5*hdim+1:6*hdim, k);

    md_true_k = md_true;
    md_true_k.geometry.thickness   = H;
    md_true_k.geometry.surface     = S;
    md_true_k.geometry.base        = B;
    md_true_k.geometry.bed         = bed;
    md_true_k.initialization.vx    = Vx;
    md_true_k.initialization.vy    = Vy;
    md_true_k.initialization.vel   = Vel;
    md_true_k.friction.coefficient = fc;
    md_true_k.mask.ocean_levelset  = H + bed/di;
end

% ============================================================
% Simple red-blue diverging colormap
% ============================================================
function cmap = redblue(n)
    if nargin < 1
        n = 256;
    end
    m = floor(n/2);

    r = [linspace(0,1,m), ones(1,n-m)];
    g = [linspace(0,1,m), linspace(1,0,n-m)];
    b = [ones(1,m), linspace(1,0,n-m)];

    cmap = [r(:) g(:) b(:)];
end

%% ---------------- GL evolution plot (5-panel robust) ----------------
function plot_gl_on_bed_evolution( ...
    k_array, dt, ...
    model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md, ...
    bg_field, bg_title, units, gl_mid, ...
    scalar_means_file_true, scalar_means_file_nurged, scalar_means_file_ens)

    global colorbar_gap t label_t nt t_change t_jump k_len

    if nargin < 10 || isempty(bg_field), bg_field = 'initialization.vel'; end
    if nargin < 11 || isempty(bg_title), bg_title = 'Background'; end
    if nargin < 12 || isempty(units), units = ''; end

    if ~isempty(units)
        units_str = [' (' units ')'];
    else
        units_str = '';
    end

    % nrows = 5;
    nrows = 6;
    % figure('Position',[150 150 1800 1100]); clf;
    figure('Position',[100 100 1100 800]); clf;
    axs = gobjects(nrows,1);

    % ------------------------------------------------------------
    % Read scalar time series for panel 5
    % ------------------------------------------------------------
    [ivaf_true,   t_true]   = read_ivaf_series(scalar_means_file_true,   dt);
    [ivaf_nurged, t_nurged] = read_ivaf_series(scalar_means_file_nurged, dt);
    [ivaf_ens,    t_ens]    = read_ivaf_series(scalar_means_file_ens,    dt);

    % ------------------------------------------------------------
    % Initial state: needed for mesh/grid and reference field
    % ------------------------------------------------------------
    [md_true_0, md_nurged_0, md_ens_0] = setup_model_states(1, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);

    x = md_true_0.mesh.x(:);
    y = md_true_0.mesh.y(:);

    Nx = 420; Ny = 70;
    xg = linspace(min(x), max(x), Nx);
    yg = linspace(min(y), max(y), Ny);
    [Xg, Yg] = meshgrid(xg, yg);

    % ------------------------------------------------------------
    % Global color limits from TRUE background over requested states
    % ------------------------------------------------------------
    all_data = [];
    kk_list = unique([1, k_array(:)', t_change]);
    kk_list = kk_list(kk_list >= 1);

    for kk = kk_list
        [md_true_k, ~, ~] = setup_model_states(kk, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
        bg = get_nested_field(md_true_k, bg_field);
        all_data = [all_data; bg(:)];
    end

    cmin = min(all_data, [], 'omitnan');
    cmax = max(all_data, [], 'omitnan');
    if ~isfinite(cmin) || ~isfinite(cmax) || cmin == cmax
        cmin = 0; cmax = 1;
    end
    clear all_data

    % ------------------------------------------------------------
    % GL filtering knobs
    % ------------------------------------------------------------
    minLen_true  = 3e4;
    minLen_wrong = 5e4;
    minLen_ens   = 5e4;
    minArea      = 1;

    keepLargestOnly_true  = true;
    keepLargestOnly_wrong = false;
    keepLargestOnly_ens   = false;
    keepTopK_true  = 4;
    keepTopK_wrong = 4;
    keepTopK_ens   = 4;

    panel_idx = '';

    % ------------------------------------------------------------
    % Reference field
    % ------------------------------------------------------------
    data_true = get_nested_field(md_true_0, bg_field);

    % ============================================================
    % (a) Initial state
    % ============================================================
    plotmodel(md_true_0, 'data', data_true, ...
        'title', sprintf('Initial %s', lower(bg_title)), ...
        'subplot', [nrows,1,1], ...
        'caxis', [cmin cmax], ...
        'colorbar', 'off');

    ax = gca; axs(1) = ax;
    format_map_panel(ax);
    text(ax, 0.02, 0.95, sprintf('(%c_{%d})','a',panel_idx), ...
        'Units','normalized', 'FontWeight','bold', 'FontSize',16, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', 'Color','k');

    overlay_gl_panel(ax, md_true_0, md_nurged_0, md_ens_0, ...
        x, y, Xg, Yg, ...
        minLen_true, minLen_wrong, minLen_ens, minArea, ...
        keepLargestOnly_true, keepLargestOnly_wrong, keepLargestOnly_ens, ...
        keepTopK_true, keepTopK_wrong, keepTopK_ens, ...
        [gl_mid.x(1), gl_mid.y(1)]);

    % ============================================================
    % (b) Unperturbed / no assimilation
    % ============================================================
    [md_true_1, md_nurged_1, md_ens_1] = setup_model_states(t_change, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);

    nurged_field = get_nested_field(md_nurged_1, bg_field);
    diff_no = (nurged_field - data_true)./data_true;
    maxAbs_global = safe_absmax(diff_no);

    t_change_plot = t_change;
    if exist('t','var') && ~isempty(t) && t_change <= numel(t)
        t_change_plot = t(t_change);
    end

    plotmodel(md_nurged_1, 'data', diff_no, ...
        'title', sprintf('unperturbed basal friction after %.0f years', round(t_change_plot)), ...
        'subplot', [nrows,1,2], ...
        'caxis', [-maxAbs_global maxAbs_global], ...
        'colorbar', 'off');

    ax = gca; axs(2) = ax;
    format_map_panel(ax);
    text(ax, 0.02, 0.95, sprintf('(%c_{%d})','b',panel_idx), ...
        'Units','normalized', 'FontWeight','bold', 'FontSize',16, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', 'Color','k');

    overlay_gl_panel(ax, md_true_1, md_nurged_1, md_ens_1, ...
        x, y, Xg, Yg, ...
        minLen_true, minLen_wrong, minLen_ens, minArea, ...
        keepLargestOnly_true, keepLargestOnly_wrong, keepLargestOnly_ens, ...
        keepTopK_true, keepTopK_wrong, keepTopK_ens, ...
        [gl_mid.x(k_len), gl_mid.y(k_len)]);

    % ============================================================
    % (c) Perturbed truth
    % ============================================================
    [md_true_change, md_nurged_change, md_ens_change] = setup_model_states(t_change, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);

    true_field_change = get_nested_field(md_true_change, bg_field);
    diff_true = (true_field_change - data_true)./data_true;
    maxAbs_global = max(maxAbs_global, safe_absmax(diff_true));

    plotmodel(md_true_change, 'data', diff_true, ...
        'title', sprintf('perturbed %s after %.0f years', lower(bg_title), round(t_change_plot)), ...
        'subplot', [nrows,1,3], ...
        'caxis', [-maxAbs_global maxAbs_global], ...
        'colorbar', 'off');

    ax = gca; axs(3) = ax;
    format_map_panel(ax);
    text(ax, 0.02, 0.95, sprintf('(%c_{%d})','c',panel_idx), ...
        'Units','normalized', 'FontWeight','bold', 'FontSize',16, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', 'Color','k');

    overlay_gl_panel(ax, md_true_change, md_nurged_change, md_ens_change, ...
        x, y, Xg, Yg, ...
        minLen_true, minLen_wrong, minLen_ens, minArea, ...
        keepLargestOnly_true, keepLargestOnly_wrong, keepLargestOnly_ens, ...
        keepTopK_true, keepTopK_wrong, keepTopK_ens, ...
        [gl_mid.x(k_len), gl_mid.y(k_len)]);

    % ============================================================
    % (d) Assimilated
    % ============================================================
    [md_true_a, md_nurged_a, md_ens_a] = setup_model_states(t_change, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);
    
    assim_field = get_nested_field(md_ens_a, bg_field);
    
    % zero out friction on floating ice in assimilated state
    floating = md_ens_a.mask.ocean_levelset < 0;
    % assim_field(floating) = 0;
    assim_field(floating) = data_true(floating);
    
    diff_assim = (assim_field - data_true)./data_true;
    maxAbs_global = max(maxAbs_global, safe_absmax(diff_assim));

    plotmodel(md_ens_a, 'data', diff_assim, ...
        'title', sprintf('assimilated basal friction after %.0f years', round(t_change_plot)), ...
        'subplot', [nrows,1,4], ...
        'caxis', [-maxAbs_global maxAbs_global], ...
        'colorbar', 'off');

    ax = gca; axs(4) = ax;
    format_map_panel(ax);
    text(ax, 0.02, 0.95, sprintf('(%c_{%d})','d',panel_idx), ...
        'Units','normalized', 'FontWeight','bold', 'FontSize',16, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', 'Color','k');

    overlay_gl_panel(ax, md_true_a, md_nurged_a, md_ens_a, ...
        x, y, Xg, Yg, ...
        minLen_true, minLen_wrong, minLen_ens, minArea, ...
        keepLargestOnly_true, keepLargestOnly_wrong, keepLargestOnly_ens, ...
        keepTopK_true, keepTopK_wrong, keepTopK_ens, ...
        [gl_mid.x(k_len), gl_mid.y(k_len)]);

    % --->%
    % pos = (y==40e3&& for all x);
    % plot(x, diff_assim(pos)); hold on
    % plot(x, diff_true(pos))


    % <-----%


    % ------------------------------------------------------------
    % Legend for GL overlays
    % ------------------------------------------------------------
    ax0 = axs(1);
    lg = legend(ax0);
    if ~isempty(lg) && isvalid(lg), delete(lg); end

    hold(ax0,'on');
    p1 = plot(ax0, NaN, NaN, 'k-',  'LineWidth', 3.0);
    p2 = plot(ax0, NaN, NaN, 'm--',  'LineWidth', 3.0);
    p3 = plot(ax0, NaN, NaN, 'c:', 'LineWidth', 3.0);

    lgd = legend(ax0, [p1 p2 p3], ...
        {'Perturbed GL','Unperturbed GL','Assimilated GL'}, ...
        'Orientation','horizontal', ...
        'FontSize',14, ...
        'Box','off');
    hold(ax0,'off');

    lgd.Units = 'normalized';
    lgd.Location = 'none';
    % lgd.Position = [0.52, 0.205, 0.60, 0.04];
    lgd.Position = [0.279999999999999 -0.15 0.6 0.950000000000001];


    % ============================================================
    % (e) Delta IVAF time series
    % ============================================================
    % ax = subplot(nrows,1,5);
    ax = subplot(nrows,1,[5,6]) % bottom panel spans two rows
    axs(5) = ax;
    hold(ax,'on');

    has_any = false;

    if ~isempty(ivaf_true)
        plot(ax, t_true, ivaf_true - ivaf_true(1), 'k-', 'LineWidth', 2.5);
        has_any = true;
    end

    if ~isempty(ivaf_nurged)
        plot(ax, t_true, ivaf_nurged - ivaf_nurged(1), 'm-.', 'LineWidth', 2.5);
        has_any = true;
    end

    if ~isempty(ivaf_ens)
        plot(ax, t_true, ivaf_ens - ivaf_ens(1), 'c:', 'LineWidth', 2.5);
        has_any = true;

        % ivaf_nurged_atlast=abs(ivaf_nurged(end))
        % ivaf_ens_atlast=abs(ivaf_ens(end))
        ivaf_diff_unperturbed_perturbed_200   = abs(ivaf_true(t_change)) - abs(ivaf_nurged(t_change));
        ivaf_diff_perturbed_assimilated_200   = abs(ivaf_true(t_change)) - abs(ivaf_ens(t_change));
        ivaf_diff_unperturbed_assimilated_200 = abs(ivaf_ens(t_change))  - abs(ivaf_nurged(t_change));
        fprintf('\n=== IVAF Differences at t = %.0f yr ===\n', t_true(t_change));

        fprintf('Unperturbed - Perturbed   : %+ .3e m^3 (%.2f%%)\n', ...
            ivaf_diff_unperturbed_perturbed_200, ...
            100 * ivaf_diff_unperturbed_perturbed_200 / abs(ivaf_true(t_change)));
        
        fprintf('Perturbed   - Assimilated : %+ .3e m^3 (%.2f%%)\n', ...
            ivaf_diff_perturbed_assimilated_200, ...
            100 * ivaf_diff_perturbed_assimilated_200 / abs(ivaf_true(t_change)));
        
        fprintf('Unperturbed - Assimilated : %+ .3e m^3 (%.2f%%)\n', ...
            ivaf_diff_unperturbed_assimilated_200, ...
            100 * ivaf_diff_unperturbed_assimilated_200 / abs(ivaf_true(t_change)));
        
        fprintf('=======================================\n');

    end

    if exist('t','var') && ~isempty(t) && t_jump <= numel(t)
        xline(ax, t(t_jump), 'g:', 'LineWidth', 3.0);
    else
        xline(ax, t_jump, 'g:', 'LineWidth', 3.0);
    end

    if has_any
        ylabel(ax, '\Delta IVAF', 'FontWeight','bold','FontSize',15);
        xlabel(ax, 'time (yr)', 'FontWeight','bold','FontSize',15);
        title(ax, '\Delta Ice Volume Above Floatation (IVAF)', 'FontWeight','bold', 'FontSize',15);
        grid(ax,'on');

        legend(ax, {'Perturbed','Unperturbed','Assimilated'}, ...
            'Location','best', 'FontSize',12);

        text(ax, 0.02, 0.90, sprintf('(%c_{%d})','e',panel_idx), ...
            'Units','normalized', 'FontWeight','bold', 'FontSize',16, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top', 'Color','k');
        
        ylim([-9.5e12, 0.6e12]); 
        xlim([-1.5,200]);
        % ylim([-12e12, 1e12]); 
        % xlim([-1.5,160]);
    else
        text(ax, 0.5, 0.5, 'IVAF time series not available', ...
            'Units','normalized', 'HorizontalAlignment','center', ...
            'FontWeight','bold');
        axis(ax,'off');
    end
    hold(ax,'off');

    % ------------------------------------------------------------
    % Layout / styling
    % ------------------------------------------------------------
    gap = 0.03; top = 0.96; bottom = 0.08;
    avail  = top - bottom - (nrows-1)*gap;
    height = avail / nrows;

    if height < 0.05
        fig = gcf;
        fig.Position(4) = fig.Position(4) * max(1, ceil(0.05 / height));
        height = 0.05;
    end

    % map panels
    for i = 1:4
        ax = axs(i);
        pos = [0.10, bottom+(nrows-i)*(height+gap), 0.80, height];
        set(ax, 'Position', pos, ...
            'FontWeight','bold', 'FontSize',15, 'Box','on', ...
            'LineWidth',2.0, 'TickDir','out', 'TickLength',[0.004 0.004], ...
            'XGrid','off', 'YGrid','off', 'YMinorGrid','off', ...
            'XTickMode','manual', 'YTickMode','manual');

        ttl = ax.Title;
        ttl.FontSize   = 15;
        ttl.FontWeight = 'bold';
        ttl.Interpreter = 'tex';

        yl = ax.YLim / 1000;
        step = 40;
        yt_km = ceil(yl(1)/step)*step : step : floor(yl(2)/step)*step;
        set(ax,'YTick',yt_km*1000,'YTickLabel',string(yt_km),'YTickMode','manual');

        xl = ax.XLim;
        xt_km = floor(xl(1)/1000/100)*100 : 100 : ceil(xl(2)/1000/100)*100;
        set(ax, 'XTick', xt_km*1000, ...
            'XTickLabel', arrayfun(@num2str, xt_km, 'UniformOutput', false), ...
            'XTickMode','manual');

        ylabel(ax,'y (km)','FontWeight','bold','FontSize',16);

        if i < 4
            set(ax,'XTickLabel',[]);
        else
            xlabel(ax,'x (km)','FontWeight','bold','FontSize',16);
            ax.XLabel.Units = 'normalized';
            ax.XLabel.Position(2) = -0.2;
        end
    end

    % time-series panel
    % ax = axs(5);
    % pos = [0.10, bottom-0.049, 0.80, height];
    % set(ax, 'Position', pos, ...
    %     'FontWeight','bold', 'FontSize',15, 'Box','on', ...
    %     'LineWidth',2.0, 'TickDir','out', 'TickLength',[0.004 0.004], ...
    %     'XGrid','on', 'YGrid','on');

    % time-series panel (make it taller)
    ax = axs(5);
    pos = [0.14, 0.07, 0.72, 0.21];   % [left bottom width height]
    set(ax, 'Position', pos, ...
        'FontWeight','bold', 'FontSize',15, 'Box','on', ...
        'LineWidth',2.0, 'TickDir','out', 'TickLength',[0.004 0.004], ...
        'XGrid','on', 'YGrid','on');

    % ------------------------------------------------------------
    % Colorbars
    % ------------------------------------------------------------
    colorbar_gap0 = 0.87;

    cb1 = colorbar(axs(1), 'Position',[colorbar_gap0 0.78 0.012 0.17]);
    ylabel(cb1,[bg_title units_str],'FontSize',14,'FontWeight','bold');
    colormap(axs(1), parula);

    for i = 2:4
        colormap(axs(i), redblue(256));
        caxis(axs(i), [-maxAbs_global maxAbs_global]);
    end

    cb2 = colorbar(axs(4), 'Position',[colorbar_gap0 0.38 0.012 0.38]);
    % ylabel(cb2, ['\Delta ' bg_title units_str], 'FontSize',15, 'FontWeight','bold');
    ylabel(cb2, ['\Delta ' bg_title ' / Initial Friction'], 'FontSize', 14, 'FontWeight', 'bold');

    set(gcf,'Color','w');

    % ------------------------------------------------------------
    % Save figure
    % ------------------------------------------------------------
    scriptdir = fileparts(mfilename('fullpath'));
    outdir = fullfile(scriptdir, 'figures');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    fname = fullfile(outdir, sprintf('GL_%s.png', regexprep(bg_field,'\.','_')));
    exportgraphics(gcf, fname, 'Resolution', 300);

    % ============================================================
    %  friction perturbation along centerline y = 40 km
    % ============================================================
    y0 = 40e3;

    % define full x-range
    x_line = linspace(min(x), max(x), 500);
    y_line = y0 * ones(size(x_line));

    % interpolate
    diff_true_line  = griddata(x, y, diff_true,  x_line, y_line, 'linear');
    diff_assim_line = griddata(x, y, diff_assim, x_line, y_line, 'linear');

    % handle NaNs
    diff_true_line  = fillmissing(diff_true_line,  'nearest');
    diff_assim_line = fillmissing(diff_assim_line, 'nearest');

    % plot
    figure('Position',[150 150 900 450]); clf;
    plot(x_line/1e3, diff_true_line,  'k-', 'LineWidth', 2.5); hold on;
    plot(x_line/1e3, diff_assim_line, 'c--', 'LineWidth', 2.5);

    xlabel('x (km)', 'FontWeight','bold', 'FontSize',14);
    ylabel('\Delta friction / initial friction', 'FontWeight','bold', 'FontSize',14);
    title(sprintf('Friction perturbation along centerline y = %.0f km at t = %.0f yr', ...
        y0/1e3, round(t_change_plot)), ...
        'FontWeight','bold', 'FontSize',15);

    legend({'Perturbed truth','Assimilated'}, 'Location','best', 'FontSize',12);
    grid on; box on;

    xlim([min(x) max(x)]/1e3);
    set(gca, 'FontWeight','bold', 'FontSize',15, 'LineWidth',1.5);

    % optional: save the line figure too
    fname_line = fullfile(outdir, sprintf('GL_%s_centerline_y40km.png', regexprep(bg_field,'\.','_')));
    exportgraphics(gcf, fname_line, 'Resolution', 300);

       
    % % define full x-range
    % x_line = linspace(min(x), max(x), 500);
    % y_line = y0 * ones(size(x_line));
    % 
    % % grounded masks
    % grounded_true  = md_true_change.mask.ocean_levelset >= 0;
    % grounded_assim = md_ens_a.mask.ocean_levelset      >= 0;
    % 
    % % build smooth interpolants from grounded nodes only
    % F_true = scatteredInterpolant( ...
    %     x(grounded_true), ...
    %     y(grounded_true), ...
    %     diff_true(grounded_true), ...
    %     'natural', 'none');
    % 
    % F_assim = scatteredInterpolant( ...
    %     x(grounded_assim), ...
    %     y(grounded_assim), ...
    %     diff_assim(grounded_assim), ...
    %     'natural', 'none');
    % 
    % % evaluate along centerline
    % diff_true_line  = F_true(x_line, y_line);
    % diff_assim_line = F_assim(x_line, y_line);
    % 
    % % exact imposed Gaussian along the centerline (optional diagnostic)
    % % only exact if y0 here matches the perturbation center used in truth
    % x0_gauss = ((0.5 * (min(x) + max(x)))./2.0);
    % y0_gauss = 0.5 * (min(y) + max(y));
    % sigma_gauss = 15e3;
    % gauss_exact = exp(-((x_line - x0_gauss).^2 + (y_line - y0_gauss).^2) ./ (2*sigma_gauss^2));
    % 
    % % plot
    % figure('Position',[150 150 900 450]); clf;
    % plot(x_line/1e3, diff_true_line,  'k-',  'LineWidth', 2.5); hold on;
    % plot(x_line/1e3, diff_assim_line, 'c--', 'LineWidth', 2.5);
    % % plot(x_line/1e3, gauss_exact, 'r:', 'LineWidth', 2.0);   % optional
    % 
    % xlabel('x (km)', 'FontWeight','bold', 'FontSize',14);
    % ylabel('\Delta friction / initial friction', 'FontWeight','bold', 'FontSize',14);
    % title(sprintf('Grounded friction perturbation along centerline y = %.0f km at t = %.0f yr', ...
    %     y0/1e3, round(t_change_plot)), ...
    %     'FontWeight','bold', 'FontSize',15);
    % 
    % legend({'Perturbed truth','Assimilated'}, 'Location','best', 'FontSize',12);
    % % legend({'Perturbed truth','Assimilated','Exact Gaussian'}, 'Location','best', 'FontSize',12);
    % 
    % grid on;
    % box on;
    % xlim([min(x) max(x)]/1e3);
    % set(gca, 'FontWeight','bold', 'FontSize',15, 'LineWidth',1.5);
    % 
    % % optional: save the line figure too
    % fname_line = fullfile(outdir, sprintf('GL_%s_centerline_y40km_grounded.png', ...
    %     regexprep(bg_field,'\.','_')));
    % exportgraphics(gcf, fname_line, 'Resolution', 300);

   
end

% ============================================================
% Helper: read IVAF time series from HDF5
% ============================================================
function [ivaf, t_scalar] = read_ivaf_series(fname, dt)
    ivaf = [];
    t_scalar = [];

    if ~exist('fname','var') || isempty(fname) || ~isfile(fname)
        return;
    end

    info = h5info(fname);
    dnames = string({info.Datasets.Name});

    if any(dnames == "IceVolumeAboveFloatation")
        ivaf = h5read(fname, '/IceVolumeAboveFloatation');
        ivaf = ivaf(:);
    end

    % if any(dnames == "IceVolume")
    %     ivaf = h5read(fname, '/IceVolume');
    %     ivaf = ivaf(:);
    % end

    if any(dnames == "time")
        t_scalar = h5read(fname, '/time');
        t_scalar = t_scalar(:);
    elseif ~isempty(ivaf)
        t_scalar = (0:length(ivaf)-1)' * dt;
    end
end

% ============================================================
% Helper: safe max abs for caxis
% ============================================================
function v = safe_absmax(A)
    v = max(abs(A(:)), [], 'omitnan');
    if ~isfinite(v) || v == 0
        v = 1;
    end
end

% ============================================================
% Helper: format map axes quickly
% ============================================================
function format_map_panel(ax)
    ttl = ax.Title;
    ttl.FontSize   = 10;
    ttl.FontWeight = 'bold';
    ttl.Interpreter = 'tex';

    xt = get(ax,'XTick');
    yt = get(ax,'YTick');
    set(ax,'XTickLabel', arrayfun(@(v) sprintf('%g', v./1000), xt, 'UniformOutput', false));
    set(ax,'YTickLabel', arrayfun(@(v) sprintf('%g', v./1000), yt, 'UniformOutput', false));
    ylabel(ax,'y (km)','FontWeight','bold','FontSize',16);
end

% ============================================================
% Helper: overlay GL contours
% ============================================================
function overlay_gl_panel(ax, md_true_k, md_nurged_k, md_ens_k, ...
    x, y, Xg, Yg, ...
    minLen_true, minLen_wrong, minLen_ens, minArea, ...
    keepLargestOnly_true, keepLargestOnly_wrong, keepLargestOnly_ens, ...
    keepTopK_true, keepTopK_wrong, keepTopK_ens, gl_xy)

    phi_true  = md_true_k.mask.ocean_levelset(:);
    phi_wrong = md_nurged_k.mask.ocean_levelset(:);
    phi_ens   = md_ens_k.mask.ocean_levelset(:);

    F1 = scatteredInterpolant(x, y, phi_true,  'linear','nearest');
    F2 = scatteredInterpolant(x, y, phi_wrong, 'linear','nearest');
    F3 = scatteredInterpolant(x, y, phi_ens,   'linear','nearest');

    Phi_true  = F1(Xg, Yg);
    Phi_wrong = F2(Xg, Yg);
    Phi_ens   = F3(Xg, Yg);

    hold(ax,'on');
    plot_gl_contour_filtered(Xg, Yg, Phi_true,  'k','-', 3.0, minLen_true,  minArea, keepLargestOnly_true,  keepTopK_true);
    plot_gl_contour_filtered(Xg, Yg, Phi_wrong, 'm','-.', 3.0, minLen_wrong, minArea, keepLargestOnly_wrong, keepTopK_wrong);
    plot_gl_contour_filtered(Xg, Yg, Phi_ens,   'c',':', 3.0, minLen_ens,   minArea, keepLargestOnly_ens,   keepTopK_ens);

    overlay_gl_window_points(ax, md_true_k, md_nurged_k, md_ens_k, gl_xy, ...
        'x_halfwidth',30e3, 'y_halfwidth',20e3);
    hold(ax,'off');
end

function h = plot_gl_contour_filtered( ...
    Xg, Yg, Phig, line_color, line_style, lw, minLen, minArea, keepLargestOnly, keepTopK)

    if nargin < 6 || isempty(lw), lw = 2.0; end
    if nargin < 7 || isempty(minLen), minLen = 2e4; end
    if nargin < 8 || isempty(minArea), minArea = 0; end
    if nargin < 9 || isempty(keepLargestOnly), keepLargestOnly = false; end
    if nargin < 10 || isempty(keepTopK), keepTopK = 1; end

    h = gobjects(0);

    if ~all(isfinite(Phig(:)))
        Phig(~isfinite(Phig)) = 1;
    end

    if min(Phig(:)) * max(Phig(:)) > 0
        return;
    end

    C = contourc(Xg(1,:), Yg(:,1), Phig, [0 0]);

    segs  = {};
    lens  = [];
    areas = [];

    kk = 1;
    while kk < size(C,2)
        npts = C(2,kk);
        pts  = C(:, kk+1:kk+npts);
        kk   = kk + npts + 1;

        xx = pts(1,:); yy = pts(2,:);
        L = sum(hypot(diff(xx), diff(yy)));

        isClosed = hypot(xx(1)-xx(end), yy(1)-yy(end)) < 1e-6 + 0.02*max(1, mean(hypot(diff(xx),diff(yy))));
        A = 0;
        if isClosed
            A = abs(polyarea(xx, yy));
        end

        segs{end+1}  = pts; %#ok<AGROW>
        lens(end+1)  = L;   %#ok<AGROW>
        areas(end+1) = A;   %#ok<AGROW>
    end

    if isempty(segs), return; end

    keep = true(1,numel(segs));
    keep = keep & (lens >= minLen);

    if minArea > 0
        keep = keep & (areas >= minArea | areas == 0);
    end

    idx = find(keep);
    if isempty(idx), return; end

    [~, order] = sort(lens(idx), 'descend');
    idx = idx(order);

    if keepLargestOnly
        idx = idx(1);
    else
        idx = idx(1:min(keepTopK, numel(idx)));
    end

    hold on
    for i = 1:numel(idx)
        pts = segs{idx(i)};
        h(end+1) = plot(pts(1,:), pts(2,:), ...
            'Color', line_color, 'LineStyle', line_style, 'LineWidth', lw); %#ok<AGROW>
    end
end

function [md_true, md_nurged, md_ens] = setup_model_states( ...
    k, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md)
% Initialize the true, wrong (nurged), and ensemble-mean models at index k.
    global nvar ensemble_vec_full
    ndim = size(model_true_state,1);
    hdim = ndim / nvar;
    di   = md.materials.rho_ice / md.materials.rho_water;

    % --- TRUE ---
    H  = model_true_state(1:hdim, k);
    S  = model_true_state(hdim+1:2*hdim, k);
    B  = S - H;
    Vx = model_true_state(2*hdim+1:3*hdim, k);
    Vy = model_true_state(3*hdim+1:4*hdim, k);
    Vel= hypot(Vx, Vy);
    bed= model_true_state(4*hdim+1:5*hdim, k);
    fc = model_true_state(5*hdim+1:6*hdim, k);

    md_true.geometry.thickness      = H;
    md_true.geometry.surface        = S;
    md_true.geometry.base           = B;
    md_true.geometry.bed            = bed;
    md_true.initialization.vx       = Vx;
    md_true.initialization.vy       = Vy;
    md_true.initialization.vel      = Vel;
    md_true.friction.coefficient    = fc;
    md_true.mask.ocean_levelset     = H + bed/di;

    % --- WRONG (nurged) ---
    H  = model_nurged_state(1:hdim, k);
    S  = model_nurged_state(hdim+1:2*hdim, k);
    B  = S - H;
    Vx = model_nurged_state(2*hdim+1:3*hdim, k);
    Vy = model_nurged_state(3*hdim+1:4*hdim, k);
    Vel= hypot(Vx, Vy);
    bed= model_nurged_state(4*hdim+1:5*hdim, k);
    fc = model_nurged_state(5*hdim+1:6*hdim, k);

    md_nurged.geometry.thickness    = H;
    md_nurged.geometry.surface      = S;
    md_nurged.geometry.base         = B;
    md_nurged.geometry.bed          = bed;
    md_nurged.initialization.vx     = Vx;
    md_nurged.initialization.vy     = Vy;
    md_nurged.initialization.vel    = Vel;
    md_nurged.friction.coefficient  = fc;
    md_nurged.mask.ocean_levelset   = H + bed/di;

    % --- ENSEMBLE MEAN ---
    ensemble_vec_mean(:,1) = model_nurged_state(:, 1);
    H  = ensemble_vec_mean(1:hdim, k);
    S  = ensemble_vec_mean(hdim+1:2*hdim, k);
    B  = S - H;
    Vx = ensemble_vec_mean(2*hdim+1:3*hdim, k);
    Vy = ensemble_vec_mean(3*hdim+1:4*hdim, k);
    Vel= hypot(Vx, Vy);
    bed= ensemble_vec_mean(4*hdim+1:5*hdim, k);
    fc = ensemble_vec_mean(5*hdim+1:6*hdim, k);

    md_ens.geometry.thickness       = H;
    md_ens.geometry.surface         = S;
    md_ens.geometry.base            = B;
    md_ens.geometry.bed             = bed;
    md_ens.initialization.vx        = Vx;
    md_ens.initialization.vy        = Vy;
    md_ens.initialization.vel       = Vel;
    md_ens.friction.coefficient     = fc;
    md_ens.mask.ocean_levelset      = H + bed/di;
end

function gl_mid = compute_gl_midpoints( ...
    k_array, dt, ...
    model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md)
% Compute TRUE grounding line midpoint (arc-length midpoint of longest 0-contour)
% for each k in k_array. Also stores a k->index map.

    nk = numel(k_array);
    gl_mid.k_array = k_array(:);
    gl_mid.x = nan(nk,1);
    gl_mid.y = nan(nk,1);
    gl_mid.k_to_idx = containers.Map('KeyType','double','ValueType','double');

    % contour grid resolution (match your GL plot)
    Nx = 420; Ny = 70;

    for idx = 1:nk
        k = k_array(idx);
        gl_mid.k_to_idx(k) = idx;

        [md_true_k, ~, ~] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);

        x = md_true_k.mesh.x(:);
        y = md_true_k.mesh.y(:);
        phi_true = md_true_k.mask.ocean_levelset(:);

        xg = linspace(min(x), max(x), Nx);
        yg = linspace(min(y), max(y), Ny);
        [Xg, Yg] = meshgrid(xg, yg);

        F = scatteredInterpolant(x, y, phi_true, 'linear', 'nearest');
        Phi = F(Xg, Yg);

        y_center = 0.5*(min(y) + max(y));  % channel centerline
        [xc, yc] = gl_centerline_point_from_levelset_grid(Xg, Yg, Phi, y_center);
        gl_mid.x(idx) = xc;
        gl_mid.y(idx) = yc;
    end
end

function [xc, yc] = gl_centerline_point_from_levelset_grid(Xg, Yg, Phig, y_center)
% Return GL point closest to centerline y=y_center, from longest 0-contour segment.

    xc = NaN; yc = NaN;

    if ~all(isfinite(Phig(:)))
        Phig(~isfinite(Phig)) = 1;
    end

    if min(Phig(:)) * max(Phig(:)) > 0
        return; % no GL
    end

    C = contourc(Xg(1,:), Yg(:,1), Phig, [0 0]);

    segs = {};
    lens = [];

    k = 1;
    while k < size(C,2)
        npts = C(2,k);
        pts  = C(:, k+1:k+npts);
        k    = k + npts + 1;

        gx = pts(1,:); gy = pts(2,:);
        L  = sum(hypot(diff(gx), diff(gy)));

        segs{end+1} = pts; %#ok<AGROW>
        lens(end+1) = L;   %#ok<AGROW>
    end

    if isempty(segs), return; end

    % longest segment = main GL
    [~, imax] = max(lens);
    pts = segs{imax};
    gx = pts(1,:); gy = pts(2,:);

    % pick point closest to centerline y_center
    [~, j] = min(abs(gy - y_center));

    xc = gx(j);
    yc = gy(j);

    % OPTIONAL: refine by linearly interpolating a crossing if there is one
    % (keeps point exactly on y_center when possible)
    sgn = sign(gy - y_center);
    idx = find(sgn(1:end-1).*sgn(2:end) <= 0); % sign change or hit
    if ~isempty(idx)
        i = idx(1);
        y1 = gy(i); y2 = gy(i+1);
        x1 = gx(i); x2 = gx(i+1);
        if abs(y2-y1) > eps
            t = (y_center - y1) / (y2 - y1);
            xc = x1 + t*(x2-x1);
            yc = y_center;
        end
    end
end

function y = iff(c,a,b)
    if c, y=a; else, y=b; end
end

function out = get_nested_field(s, field)
% Access nested fields with dot notation: e.g. 'geometry.base'
    parts = strsplit(field,'.');
    out = s;
    for i = 1:numel(parts)
        tok = parts{i};
        t = regexp(tok,'(.+)\((\d+)\)$','tokens');
        if ~isempty(t)
            out = out.(t{1}{1})(str2double(t{1}{2}));
        else
            out = out.(tok);
        end
    end
end

function overlay_gl_window_points(ax, md_true_k, md_nurged_k, md_ens_k, gl_mid_k, varargin)
% Overlay GL points used for windowed RMSE on an existing GL plot.
%
% Usage:
%   overlay_gl_window_points(gca, md_true_k, md_nurged_k, md_ens_k, gl_mid_k, ...
%       'x_halfwidth',30e3,'y_halfwidth',20e3,'Nx',420,'Ny',70,'minLen',3e4,'topK',4);

    % ---- parse name/value options (compatible with older MATLAB) ----
    p = inputParser;
    p.addParameter('x_halfwidth', 30e3, @(v) isnumeric(v) && isscalar(v));
    p.addParameter('y_halfwidth', 20e3, @(v) isnumeric(v) && isscalar(v));
    p.addParameter('Nx', 420, @(v) isnumeric(v) && isscalar(v));
    p.addParameter('Ny', 70, @(v) isnumeric(v) && isscalar(v));
    p.addParameter('minLen', 3e4, @(v) isnumeric(v) && isscalar(v));
    p.addParameter('topK', 4, @(v) isnumeric(v) && isscalar(v));
    p.parse(varargin{:});
    opt = p.Results;

    % gl_mid_k can be struct with fields x,y OR a 1x2 vector
    if isstruct(gl_mid_k)
        xc = gl_mid_k.x;  yc = gl_mid_k.y;
    else
        xc = gl_mid_k(1); yc = gl_mid_k(2);
    end

    if ~isfinite(xc) || ~isfinite(yc)
        return
    end

    win = [xc-opt.x_halfwidth, xc+opt.x_halfwidth, ...
           yc-opt.y_halfwidth, yc+opt.y_halfwidth];

    % Extract GL points used in RMSE (must exist in your file)
    Ptrue = extract_gl_points_in_window(md_true_k,   opt.Nx, opt.Ny, opt.minLen, opt.topK, win);
    Pno   = extract_gl_points_in_window(md_nurged_k, opt.Nx, opt.Ny, opt.minLen, opt.topK, win);
    Pas   = extract_gl_points_in_window(md_ens_k,    opt.Nx, opt.Ny, opt.minLen, opt.topK, win);

    hold(ax,'on')

    % center marker
    % plot(ax, xc, yc, 'go', 'MarkerFaceColor','g', ...
    %     'MarkerSize',8,'LineWidth',1.5,'HandleVisibility','on');

    plot(ax, xc, yc,  'Marker', 'o','Color', '#664222', 'MarkerFaceColor','#664222'  , ...
        'MarkerSize',10,'LineWidth',1.5,'HandleVisibility','on');

    hold(ax,'off')
end

function P = extract_gl_points_in_window(md_k, Nx, Ny, minLen, topK, win)
% Returns Nx2 array of [x y] points on ocean_levelset=0 within window.
% win = [xmin xmax ymin ymax].

    x = md_k.mesh.x(:);
    y = md_k.mesh.y(:);
    phi = md_k.mask.ocean_levelset(:);

    xg = linspace(min(x), max(x), Nx);
    yg = linspace(min(y), max(y), Ny);
    [Xg, Yg] = meshgrid(xg, yg);

    F = scatteredInterpolant(x, y, phi, 'linear', 'nearest');
    Phi = F(Xg, Yg);

    % no sign change => no contour
    if min(Phi(:)) * max(Phi(:)) > 0
        P = zeros(0,2);
        return;
    end

    C = contourc(Xg(1,:), Yg(:,1), Phi, [0 0]);

    segs = {};
    lens = [];

    kk = 1;
    while kk < size(C,2)
        npts = C(2,kk);
        pts  = C(:, kk+1:kk+npts);
        kk   = kk + npts + 1;

        gx = pts(1,:); gy = pts(2,:);
        L  = sum(hypot(diff(gx), diff(gy)));

        if L >= minLen
            segs{end+1} = pts; %#ok<AGROW>
            lens(end+1) = L;   %#ok<AGROW>
        end
    end

    if isempty(segs)
        P = zeros(0,2);
        return;
    end

    % keep topK longest segments (helps with messy multi-loops)
    [~, ord] = sort(lens, 'descend');
    ord = ord(1:min(topK, numel(ord)));

    % concatenate points and clip to window
    P = [];
    xmin = win(1); xmax = win(2); ymin = win(3); ymax = win(4);

    for i = 1:numel(ord)
        pts = segs{ord(i)};
        gx = pts(1,:); gy = pts(2,:);

        in = (gx >= xmin & gx <= xmax & gy >= ymin & gy <= ymax);
        P = [P; [gx(in)' gy(in)']]; %#ok<AGROW>
    end

    % optional: thin duplicates a bit
    if size(P,1) > 1
        P = unique(round(P,3),'rows'); % mm-level rounding, helps stability
    end
end

