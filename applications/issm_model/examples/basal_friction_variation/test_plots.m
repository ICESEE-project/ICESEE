close all; clearvars; clear all

global data_file_paths nvar ensemble_vec_full ...
        label_t t nt colorbar_gap bed_obs_xy

% ============================================================
% Load time vector
% ============================================================
results_dir = 'results';
filter_type = 'true-wrong';
file_path   = fullfile(results_dir, sprintf('%s-issm.h5', filter_type));
t = h5read(file_path,'/t');

% ============================================================
% Load true/nurged states
% ============================================================
data_file_paths = '_modelrun_datasets';
file_path = fullfile(data_file_paths, 'true_nurged_states.h5');

model_true_state   = h5read(file_path,'/true_state')';
model_nurged_state = h5read(file_path,'/nurged_state')';

% ============================================================
% Load ISSM model template
% ============================================================
md = loadmodel(fullfile('data','ISMIP.Parameterization1.mat'));
md_true = md;

% Optional colorbar position used by your functions
global colorbar_gap
colorbar_gap = 0.92;

% Example requested times
times_to_plot = [1.0, 6.2, 12.2, 18.2, 22.2, 30.0];

k_array = zeros(size(times_to_plot));
for i = 1:length(times_to_plot)
    [~, k_array(i)] = min(abs(t - times_to_plot(i)));
end

plot_true_friction_snapshots(k_array, model_true_state, md_true, md, '');
plot_true_friction_anomaly_snapshots(k_array, model_true_state, md_true, md, '');

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