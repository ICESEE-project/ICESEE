%% -----------------------------------------------------------
% @author:  Brian Kyanjo
% @date:    2025-06-30
% @brief:   Reads and plots results from both ISSM and ICESEE
% ------------------------------------------------------------

close all; clearvars; clear all

% shg;

global data_file_paths nvar ensemble_vec_full ...
        label_t t nt colorbar_gap bed_obs_xy
data_file_paths = '_modelrun_datasets';
% data_file_paths = '_goodgrounding';
% data_file_paths ='_modelrun_working_0';
nvar = 6;
colorbar_gap=0.92;

% ---------------- user toggles ----------------
make_plots       = 0;
make_multi_plots = 1;   % <-- ON (restored)
frames_plot      = 0;
compute_rmse     = 0;
plotgl           = 0;

% ---------------- time steps ------------------
% k_array = [0, 20,  60, 80, 89, 130, 330, 499]+1;
% k_array= [ 0, 20,80, 120, 160, 220, 250, 320, 450]+1;
% k_array = [0, 20, 80, 120, 160, 240, 360, 499] +1;
k_array = [30, 70,100, 120, 180, 249]+1;
dt      = 0.2;

% ---------------- Load essentials --------------
results_dir = 'results';
filter_type = 'true-wrong';
file_path   = fullfile(results_dir, sprintf('%s-issm.h5', filter_type));
t        = h5read(file_path,'/t'); 
ind_m    = h5read(file_path,'/obs_index'); 
tm_m     = h5read(file_path,'/obs_max_time'); 
run_mode = h5read(file_path,'/run_mode'); 

% --------- true / wrong (nurged)
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

% ---------------- GL midpoint points + pointwise RMSE ----------------
[gl_mid] = compute_gl_midpoints( ...
    k_array, dt, ...
    model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md);

%% ------ RMSE -----
if compute_rmse
    compute_rmse_timeseries(k_array, dt, t, model_true_state, model_nurged_state, ensemble_vec_mean,md_true, md_nurged, md_ens, md, 'geometry.thickness');
end

% ---------------- GL evolution plot ------------
if plotgl
    plot_gl_on_bed_evolution( ...
        k_array, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, ...
        'geometry.thickness', 'Thickness', 'm', gl_mid);
end
% ---------------- multi-plots restored ----------
if make_multi_plots
    % thickness
    plot_var_diff(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'geometry.thickness', 'Thickness', 'm');
    plot_var_evolution(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'geometry.thickness', 'Thickness', 'm');

    % surface
    plot_var_diff(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'geometry.surface', 'Surface', 'm');
    plot_var_evolution(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'geometry.surface', 'Surface', 'm');

    % velocity
    plot_var_evolution(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'initialization.vel', 'Velocity', 'm/yr');
    plot_var_diff(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'initialization.vel', 'Velocity', 'm/yr');

    % bed
    plot_var_evolution(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'geometry.bed', 'Bed Elevation', 'm');
    plot_var_diff(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'geometry.bed', 'Bed', 'm');

    % friction
    plot_var_evolution(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'friction.coefficient', 'Friction Coefficient', 'Pa m^{-1/3} yr^{-1/3}');
    plot_var_diff(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'friction.coefficient', 'Friction', 'Pa m^{-1/3} yr^{-1/3}');
end

% ---------------- optional single triptych -------
if make_plots
    global t, nt, label_t
    k = k_array(end);
    label_t = iff(k == nt-1, t(nt), t(k));
    [md_true_k, md_nurged_k, md_ens_k] = setup_model_states(k, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);

    plot_triptych(md_true_k, md_nurged_k, md_ens_k, ...
        'geometry.thickness', sprintf('\\bfIce Thickness after %s years', fmt_years(label_t)), parula, 'm');
    plot_triptych(md_true_k, md_nurged_k, md_ens_k, ...
        'geometry.surface', sprintf('\\bfIce Surface after %s years', fmt_years(label_t)), parula, 'm');
    plot_triptych(md_true_k, md_nurged_k, md_ens_k, ...
        'geometry.bed', sprintf('\\bfBed after %s years', fmt_years(label_t)), parula, 'm');
    plot_triptych(md_true_k, md_nurged_k, md_ens_k, ...
        'mask.ocean_levelset', sprintf('\\bfGrounding Line after %s years', fmt_years(label_t)), parula, 'm');
    
end

%% ========================================================================
%%  Helper functions
%% ========================================================================

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

function y = iff(c,a,b)
    if c, y=a; else, y=b; end
end

function cmap = redblue(n)
    if nargin<1, n=256; end
    m = n/2;
    r=[linspace(0,1,m) ones(1,m)];
    g=[linspace(0,1,m) linspace(1,0,m)];
    b=[ones(1,m) linspace(1,0,m)];
    cmap=[r(:) g(:) b(:)];
end

%% ---------------- GL evolution plot (robust) ----------------
function plot_gl_on_bed_evolution( ...
    k_array, dt, ...
    model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md, ...
    bg_field, bg_title, units, gl_mid)

    global colorbar_gap

    if nargin < 11 || isempty(bg_field), bg_field = 'initialization.vel'; end
    if nargin < 12 || isempty(bg_title), bg_title = 'Background'; end
    if nargin < 13, units = ''; end
    units_str = iff(~isempty(units), [' (' units ')'], '');

    nk    = numel(k_array);
    nrows = 2+nk;

    % figure('Position',[400 400 1100 (180 + 150*nrows)]); clf;
     figure('Position',[150 150 1800 1000]); clf;

    % ---- global color limits from TRUE background ----
    all_data = [];
    for kk = k_array
        [md_true_k, ~, ~] = setup_model_states(kk, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
        bg = get_nested_field(md_true_k, bg_field);
        all_data = [all_data; bg(:)]; %
    end
    cmin = min(all_data);
    cmax = max(all_data);
    clear all_data

    % ---- GL grid resolution ----
    Nx = 420; Ny = 70;

    % ---- filtering knobs ----
    minLen_true  = 3e4;
    minLen_wrong = 5e4;
    minLen_ens   = 5e4;
    minArea      = 1;

    keepLargestOnly_true  = true;
    keepLargestOnly_wrong = false;
    keepLargestOnly_ens   = false;
    keepTopK_ens          = 4;   % allow 2 longest for ensemble (prevents “loss”)
    keepTopK_true         = 4;
    keepTopK_wrong        = 4;
    global t label_t nt
    nt = 251;
    axs = gobjects(nrows,1);   % <-- store ONLY the real panel axes

    % (a) True
    [md_true_k, ~, ~] = setup_model_states(1, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);
    data_true = get_nested_field(md_true_k, bg_field);
    plotmodel(md_true_k,'data',data_true,'title',sprintf('True %s',lower(bg_title)), ...
        'subplot',[nrows,1,1],'caxis',[cmin cmax],'colorbar','off');
    ax = gca; axs(1) = ax;
    ttl = ax.Title;
    ttl.FontSize   = 10;
    ttl.FontWeight = 'bold';   % or 'normal'
    ttl.Interpreter = 'tex';
    
    % ---- km axes (ticks shown in km) ----
    xt = get(ax,'XTick'); yt = get(ax,'YTick');
    set(ax,'XTickLabel', arrayfun(@(v) sprintf('%g', v./1000), xt, 'UniformOutput', false));
    set(ax,'YTickLabel', arrayfun(@(v) sprintf('%g', v./1000), yt, 'UniformOutput', false));
    % xlabel(ax,'x (km)','FontWeight','bold');
    ylabel(ax,'y (km)','FontWeight','bold','FontSize',16);
    
    % ---- panel letter inside upper-left ----
    % panel = sprintf('(%c)', 'a'+(idx-1));
    panel_idx = 2;   % change as needed
    panel = sprintf('(%c_{%d})','a', panel_idx);
    text(ax, 0.02, 0.95, panel, 'Units','normalized', ...
        'FontWeight','bold', 'FontSize', 16, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'Color','k');

    % (b) No assimilation - True
    [md_true_1, md_nurged_1, md_ens_1] = setup_model_states(1, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);
    % diff_no = get_nested_field(md_ens_1, field) - get_nested_field(md_true_1, field);
    % ens_field = get_nested_field(md_nurged_1, bg_field);
    ens_field = get_nested_field(md_ens_1, bg_field);
    true_field = get_nested_field(md_true_1, bg_field);
   
    diff_no = ens_field - true_field;
    % diff_no = relative_error(ens_field, true_field);
    % diff_no = signed_log_relerr(ens_field, true_field);
    maxAbs_no = max(abs(diff_no(:)));
    % maxAbs_no=1;
    maxAbs_global = max(abs(diff_no(:)));
  
    % maxAbs_no = prctile(abs(diff_no(:)), 99);
    plotmodel(md_ens_1,'data',diff_no,'title',sprintf('no assimilation'), ...
        'subplot',[nrows,1,2],'caxis',[-maxAbs_no maxAbs_no],'colorbar','off');

    ax = gca; axs(2) = ax;
    ttl = ax.Title;
    ttl.FontSize   = 10;
    ttl.FontWeight = 'bold';   % or 'normal'
    ttl.Interpreter = 'tex';
    
    % ---- km axes (ticks shown in km) ----
    xt = get(ax,'XTick'); yt = get(ax,'YTick');
    set(ax,'XTickLabel', arrayfun(@(v) sprintf('%g', v./1000), xt, 'UniformOutput', false));
    set(ax,'YTickLabel', arrayfun(@(v) sprintf('%g', v./1000), yt, 'UniformOutput', false));
    % xlabel(ax,'x (km)','FontWeight','bold');
    ylabel(ax,'y (km)','FontWeight','bold','FontSize',16);
    
    % ---- panel letter inside upper-left ----
    % panel = sprintf('(%c)', 'a'+(idx-1));
    panel = sprintf('(%c_{%d})','b', panel_idx);
    text(ax, 0.02, 0.95, panel, 'Units','normalized', ...
        'FontWeight','bold', 'FontSize', 16, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'Color','k');
    % prevent plotmodel objects appearing in legend
    set(ax.Children, 'HandleVisibility','off');

    % grid for contouring
    x = md_true_1.mesh.x(:);
    y = md_true_1.mesh.y(:);
    xg = linspace(min(x), max(x), Nx);
    yg = linspace(min(y), max(y), Ny);
    [Xg, Yg] = meshgrid(xg, yg);

    phi_true  = md_true_1.mask.ocean_levelset(:);
    phi_wrong = md_nurged_1.mask.ocean_levelset(:);
    phi_ens   = md_ens_1.mask.ocean_levelset(:);

    % linear + nearest extrap => NO NaN holes that break contour
    F1 = scatteredInterpolant(x, y, phi_true,  'linear','nearest');
    F2 = scatteredInterpolant(x, y, phi_wrong, 'linear','nearest');
    % pos = find(x<=670);
    F3 = scatteredInterpolant(x, y, phi_ens,   'linear','nearest');

    Phi_true  = F1(Xg, Yg);
    Phi_wrong = F2(Xg, Yg);
    Phi_ens   = F3(Xg, Yg);

    hold(ax,'on');
    plot_gl_contour_filtered(Xg, Yg, Phi_true,  'k','-',  3.0, minLen_true,  minArea, keepLargestOnly_true,  keepTopK_true);
    plot_gl_contour_filtered(Xg, Yg, Phi_wrong, 'm','-', 3.0, minLen_wrong, minArea, keepLargestOnly_wrong,  keepTopK_wrong);
    plot_gl_contour_filtered(Xg, Yg, Phi_ens,   'c',':',  3.0, minLen_ens,   minArea, keepLargestOnly_ens,  keepTopK_ens);

    overlay_gl_window_points(gca, md_true_1, md_nurged_1, md_ens_1, ...
    [gl_mid.x(1), gl_mid.y(1)], ...
    'x_halfwidth',30e3, 'y_halfwidth',20e3);
    hold(ax,'off');


    for idx = 1:nk
        k = k_array(idx);
        label_t = iff(k == nt-1, t(nt), t(k));
        [md_true_k, md_nurged_k, md_ens_k] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);

        bg_true = get_nested_field(md_true_k, bg_field);
        bg_ens  = get_nested_field(md_ens_k, bg_field);
        % diff_k = relative_error(bg_ens, bg_true);
        diff_k = bg_ens - bg_true;
        % diff_k = signed_log_relerr(bg_ens, bg_true);
        % maxAbs = max(abs(diff_k(:)));
        % maxAbs = maxAbs_no;
        maxAbs_global = max(maxAbs_global, max(abs(diff_k(:))));
        maxAbs = maxAbs_global;

        plotmodel(md_ens_k, 'data', diff_k, ...
            'title', sprintf('after %s years of assimilation', fmt_years(label_t)), ...
            'subplot', [nrows, 1, idx+2], ...
            'caxis', [-maxAbs maxAbs], ...
            'colorbar', 'off');
        
        ax = gca;
        ttl = ax.Title;
        ttl.FontSize   = 10;
        ttl.FontWeight = 'bold';   % or 'normal'
        ttl.Interpreter = 'tex';
        
        % ---- km axes (ticks shown in km) ----
        % axs(idx) = ax;
        axs(idx+2) = ax; 
        xt = get(ax,'XTick'); yt = get(ax,'YTick');
        set(ax,'XTickLabel', arrayfun(@(v) sprintf('%g', v./1000), xt, 'UniformOutput', false));
        set(ax,'YTickLabel', arrayfun(@(v) sprintf('%g', v./1000), yt, 'UniformOutput', false));
        % xlabel(ax,'x (km)','FontWeight','bold');
        ylabel(ax,'y (km)','FontWeight','bold','FontSize',16);
        
        % ---- panel letter inside upper-left ----
        % panel = sprintf('(%c)', 'a'+(idx-1));
        % panel_idx = 2;   % change as needed
        panel = sprintf('(%c_{%d})','c'+(idx-1), panel_idx);
        text(ax, 0.02, 0.95, panel, 'Units','normalized', ...
            'FontWeight','bold', 'FontSize', 16, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color','k');

     
        % prevent plotmodel objects appearing in legend
        set(ax.Children, 'HandleVisibility','off');

        phi_true  = md_true_k.mask.ocean_levelset(:);
        phi_wrong = md_nurged_k.mask.ocean_levelset(:);
        phi_ens   = md_ens_k.mask.ocean_levelset(:);

        % print diagnostic so you can confirm if ens truly loses sign change
        fprintf('k=%d t=%.1f | ens[min,max]=[%.3e,%.3e]\n', ...
            k, t(k), min(phi_ens), max(phi_ens));

        % linear + nearest extrap => NO NaN holes that break contour
        F1 = scatteredInterpolant(x, y, phi_true,  'linear','nearest');
        F2 = scatteredInterpolant(x, y, phi_wrong, 'linear','nearest');
        % pos = find(x<=670);
        F3 = scatteredInterpolant(x, y, phi_ens,   'linear','nearest');

        Phi_true  = F1(Xg, Yg);
        Phi_wrong = F2(Xg, Yg);
        Phi_ens   = F3(Xg, Yg);

        hold(ax,'on');
        plot_gl_contour_filtered(Xg, Yg, Phi_true,  'k','-',  3.0, minLen_true,  minArea, keepLargestOnly_true,  keepTopK_true);
        plot_gl_contour_filtered(Xg, Yg, Phi_wrong, 'm','-', 3.0, minLen_wrong, minArea, keepLargestOnly_wrong,  keepTopK_wrong);
        plot_gl_contour_filtered(Xg, Yg, Phi_ens,   'c',':',  3.0, minLen_ens,   minArea, keepLargestOnly_ens,  keepTopK_ens);

        overlay_gl_window_points(gca, md_true_k, md_nurged_k, md_ens_k, ...
        [gl_mid.x(idx), gl_mid.y(idx)], ...
        'x_halfwidth',30e3, 'y_halfwidth',20e3);
        hold(ax,'off');
    end

    xlabel(ax,'x (km)','FontWeight','bold','FontSize',18);

    % ---- Layout: adaptive spacing ----
    % axs = flipud(findall(gcf,'Type','axes'));
    % gap = 0.02; top = 0.95; bottom = 0.08;
    % avail = top - bottom - (nrows-1)*gap;
    % height = avail / nrows;
    % 
    % if height < 0.05
    %     fig = gcf;
    %     scale_factor = max(1, ceil(0.05 / height));
    %     fig.Position(4) = fig.Position(4) * scale_factor;
    %     height = 0.05;
    % end
    % 
    % for i = 1:nrows
    %     pos = [0.10, bottom+(nrows-i)*(height+gap), 0.70, height];
    %     set(axs(i), 'Position', pos, ...
    %         'FontWeight','bold','LineWidth',1.2,'Box','on', ...
    %         'TickDir','out','Layer','top','FontSize',11, ...
    %         'TickLength',[0.005 0.005]);
    %     ylabel(axs(i),'y (km)','FontWeight','bold');
    %     if i < nrows
    %         set(axs(i),'XTickLabel',[]);
    %     else
    %         xlabel(axs(i),'x (km)','FontWeight','bold');
    %     end
    % end

    % ---- Layout: adaptive spacing (ONLY panel axes) ----
    % gap = 0.03; top = 0.95; bottom = 0.08;
    gap = 0.03; top = 0.96; bottom = 0.1;
    avail  = top - bottom - (nrows-1)*gap;
    height = avail / nrows;
    
    if height < 0.05
        fig = gcf;
        scale_factor = max(1, ceil(0.05 / height));
        fig.Position(4) = fig.Position(4) * scale_factor;
        height = 0.05;
    end
    titlePad = 0.025; 
    
    for i = 1:nrows
        ax = axs(i);
    
        % i=1 should be the TOP panel (subplot does this already)
        % pos = [0.10, bottom+(nrows-i)*(height+gap), 0.70, height];
        pos = [0.10, bottom+(nrows-i)*(height+gap), 0.9, height];
        set(ax, 'Position', pos, ...
        'FontWeight','bold', ...
        'FontSize',15, ...
        'Box','on', ...
        'LineWidth',2.0, ...
        'TickDir','out', ...
        'TickLength',[0.004 0.004], ...
        'XGrid','off', ...
        'YGrid','off', ...
        'YMinorGrid','off', ...
        'XTickMode','manual', ...
        'YTickMode','manual');
        ttl = ax.Title;
        ttl.FontSize   = 15;
        ttl.FontWeight = 'bold';   % or 'normal'
        ttl.Interpreter = 'tex';
        % ---- FORCE tick locations (meters) ----
        % yt_km = 0:20:80;              % desired ticks in km
        % yt_m  = yt_km * 1000;         % convert to meters
        % 
        % set(ax, ...
        %     'YTick', yt_m, ...
        %     'YTickLabel', arrayfun(@num2str, yt_km, 'UniformOutput', false), ...
        %     'YTickMode','manual');
        yl = ax.YLim / 1000;  % km
        step = 40;            % km
        yt_km = ceil(yl(1)/step)*step : step : floor(yl(2)/step)*step;
        set(ax,'YTick',yt_km*1000,'YTickLabel',string(yt_km),'YTickMode','manual');

        % same idea for x if needed
        xt = get(ax,'XLim');
        xt_km = floor(xt(1)/1000/100)*100 : 100 : ceil(xt(2)/1000/100)*100;
        set(ax, ...
            'XTick', xt_km*1000, ...
            'XTickLabel', arrayfun(@num2str, xt_km, 'UniformOutput', false), ...
            'XTickMode','manual');
    
        ylabel(ax,'y (km)','FontWeight','bold','FontSize',16);
    
        if i < nrows
            set(ax,'XTickLabel',[]);
        else
            xlabel(ax,'x (km)','FontWeight','bold','FontSize',16);
            ax.XLabel.Units = 'normalized';
            ax.XLabel.Position(2) = -0.33;
        end

    end


    % % ---- Shared colorbar ----
    % colorbar_gap0=0.8;
    % for i = 1:nrows, colormap(axs(i), parula); end
    % cb = colorbar(axs(end), 'Position',[colorbar_gap0 0.25 0.025 0.45]);
    % ylabel(cb, [bg_title units_str], 'FontSize',12,'FontWeight','bold');
    colorbar_gap0=0.91;
    cb1 = colorbar(axs(1), 'Position',[colorbar_gap0 0.75 0.015 0.16]);
    ylabel(cb1,[bg_title units_str],'FontSize',15,'FontWeight','bold');
    colormap(axs(1), parula);

    % for i = 2:nrows, colormap(axs(i), redblue(256)); end
    % cb = colorbar(axs(end), 'Position',[colorbar_gap0 0.25 0.015 0.45]);
    % % ylabel(cb,'Relative Error','FontSize',15,'FontWeight','bold');
    % ylabel(cb,['\Delta Thickness' units_str],'FontSize',15,'FontWeight','bold');
    % % set_colorbar_ticks(cb, [-maxAbs_global maxAbs_global], 6);
    % % caxis(axs(end), [-maxAbs_global maxAbs_global]);
    % for i = 2:nrows
    %     caxis(axs(i), [-maxAbs_global maxAbs_global]);
    % end
    % --- Shared DIFF colormap + shared clim across diff panels (rows 2..nrows) ---
    for i = 2:nrows
        colormap(axs(i), redblue(256));
        caxis(axs(i), [-maxAbs_global maxAbs_global]);
    end
    
    % Create the shared colorbar *after* clims are enforced
    cb = colorbar(axs(end), 'Position',[colorbar_gap0 0.25 0.015 0.45]);
    ylabel(cb,['\Delta Thickness' units_str],'FontSize',15,'FontWeight','bold');
    
    % 5–6 reasonable ticks, symmetric
    % cb.Ticks = linspace(-maxAbs_global, maxAbs_global, 6);
    % cb.TickLabels = arrayfun(@(v) sprintf('%.0f', v), cb.Ticks, 'UniformOutput', false);
   

    % ---- Clean legend (proxy only) ----
    ax0 = axs(1);
    lg = legend(ax0);
    if ~isempty(lg) && isvalid(lg), delete(lg); end

    % hold(ax0,'on');
    % p1 = plot(ax0, NaN, NaN, 'k-',  'LineWidth', 2.0);
    % p2 = plot(ax0, NaN, NaN, 'r-', 'LineWidth', 2.0);
    % p3 = plot(ax0, NaN, NaN, 'c-.',  'LineWidth', 2.0);
    % lgd = legend(ax0, [p1 p2 p3], ...
    %     {'True GL','No assimilation GL','Assimilated GL'}, ...
    %     'Location','northwest','FontSize',10,'Box','on');
    % lgd.AutoUpdate = 'on';
    % % legend(ax0,'manual');
    % hold(ax0,'off');
    hold(ax0,'on');

    p1 = plot(ax0, NaN, NaN, 'k-',  'LineWidth', 3.0);
    p2 = plot(ax0, NaN, NaN, 'm-',  'LineWidth', 3.0);
    p3 = plot(ax0, NaN, NaN, 'c-.', 'LineWidth', 3.0);
    
    lgd = legend(ax0, [p1 p2 p3], ...
        {'True GL','No assimilation GL','Assimilated GL'}, ...
        'Orientation','horizontal', ...
        'FontSize',14, ...
        'Box','off');
    
    hold(ax0,'off');
    
    % --- Force legend OUTSIDE the axes ---
    lgd.Units = 'normalized';
    lgd.Location = 'none';      % disable MATLAB auto-placement
    
    axPos = ax0.Position;       % axes position [x y w h]
    
    % place legend just outside top-right
    % lgd.Position = [ ...
    %     axPos(1) + axPos(3) - 0.185, ...   % to the right of axes
    %     axPos(2) + axPos(4) - lgd.Position(4), ... % aligned to top
    %     lgd.Position(3), ...
        % lgd.Position(4)];

    lgd.Position = [ ...
        0.525, ...   % left (centered)
        0.015, ...  % vertical position BELOW xlabel
        0.60, ...   % width
        0.04  ...   % height
    ];

    set(gcf,'Color','w');

    % ---- Save figure (300 dpi) ----
    % Use folder relative to THIS script (not MATLAB's current folder)
    scriptdir = fileparts(mfilename('fullpath'));
    outdir    = fullfile(scriptdir, 'figures');
    
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    
    fname = fullfile(outdir, sprintf('GL_%s.png', regexprep(bg_field,'\.','_')));
    
    exportgraphics(gcf, fname, 'Resolution', 300);
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
        return; % no sign change => no GL
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

        x = pts(1,:); y = pts(2,:);
        L = sum(hypot(diff(x), diff(y)));

        isClosed = hypot(x(1)-x(end), y(1)-y(end)) < 1e-6 + 0.02*max(1, mean(hypot(diff(x),diff(y))));
        A = 0;
        if isClosed
            A = abs(polyarea(x, y));
        end

        segs{end+1}  = pts; %
        lens(end+1)  = L;   %
        areas(end+1) = A;   %
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
            'Color', line_color, 'LineStyle', line_style, 'LineWidth', lw); %
    end
end

%% ---------------- multi-plot helpers (your style) ----------------
function plot_var_diff(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md, field, field_title, units)

    global t label_t nt colorbar_gap bed_obs_xy
    % lightGray = [0.85 0.85 0.85];
    % lightGray = [0.93 0.69 0.13];
    % lightGray = 'm';
    lightGray = 'k';
    if nargin < 12, units = ''; end
    units_str = iff(~isempty(units), [' (' units ')'], '');

    nk    = length(k_array);
    nrows = 2 + nk;
    % field_title=lower(field_title);

    % figure('Position',[400 400 1100 (180 + 150*nrows)]); clf;
    figure('Position',[150 150 1800 1000]); clf;
    axs = gobjects(nrows,1);   % <-- store ONLY the real panel axes

    % global limits for absolute field (panel a)
    all_data = [];
    for k = [1, k_array]
        [md_true_tmp, ~, md_ens_tmp] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
    
        tmpT = get_nested_field(md_true_tmp, field);
        tmpE = get_nested_field(md_ens_tmp,  field);
    
        all_data = [all_data; tmpT(:); tmpE(:)];
    end
    cmin = min(all_data);
    cmax = max(all_data);
    clear all_data

    % (a) True
    [md_true_k, ~, ~] = setup_model_states(1, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);
    data_true = get_nested_field(md_true_k, field);
    plotmodel(md_true_k,'data',data_true,'title',sprintf('True %s',lower(field_title)), ...
        'subplot',[nrows,1,1],'caxis',[cmin cmax],'colorbar','off');
    
    ax = gca; axs(1) = ax;
    % --- Overlay bed observation locations on panel (a1) only ---
    if strcmp(field,'geometry.bed') && exist('bed_obs_xy','var') && ~isempty(bed_obs_xy)
    
        hold(ax,'on');
    
        hObs = scatter(ax, bed_obs_xy(:,1), bed_obs_xy(:,2), 16, ...   % 22 = marker size
            'o', ...
            'MarkerFaceColor','none', ...
            'MarkerEdgeColor',[1 1 1]*0.98, ...
            'LineWidth',0.9);
    
        % transparency (scatter supports this in modern MATLAB)
        if isprop(hObs,'MarkerEdgeAlpha'), hObs.MarkerEdgeAlpha = 0.35; end
        if isprop(hObs,'MarkerFaceAlpha'), hObs.MarkerFaceAlpha = 0; end
    
        set(hObs,'HandleVisibility','off');  % don’t pollute legend
        hold(ax,'off');
    
        try uistack(hObs,'top'); catch, end
    end

    ttl = ax.Title;
    ttl.FontSize   = 10;
    ttl.FontWeight = 'bold';   % or 'normal'
    ttl.Interpreter = 'tex';
    
    % ---- km axes (ticks shown in km) ----
    xt = get(ax,'XTick'); yt = get(ax,'YTick');
    set(ax,'XTickLabel', arrayfun(@(v) sprintf('%g', v./1000), xt, 'UniformOutput', false));
    set(ax,'YTickLabel', arrayfun(@(v) sprintf('%g', v./1000), yt, 'UniformOutput', false));
    % xlabel(ax,'x (km)','FontWeight','bold');
    ylabel(ax,'y (km)','FontWeight','bold','FontSize',16);
    
    % ---- panel letter inside upper-left ----
    % panel = sprintf('(%c)', 'a'+(idx-1));
    panel_idx = 2;   % change as needed
    panel = sprintf('(%c_{%d})','a', panel_idx);
    text(ax, 0.02, 0.95, panel, 'Units','normalized', ...
        'FontWeight','bold', 'FontSize', 16, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'Color',lightGray);

    % (b) No assimilation - True
    [md_true_1, md_nurged_1, md_ens_1] = setup_model_states(1, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);
    % diff_no = get_nested_field(md_ens_1, field) - get_nested_field(md_true_1, field);
    % ens_field = get_nested_field(md_nurged_1, field);
    ens_field = get_nested_field(md_ens_1, field);
    true_field = get_nested_field(md_true_1, field);
    % if contains(field,'geometry.bed')
        % diff_no = relative_error(ens_field, true_field);
        diff_no = signed_log_relerr(ens_field, true_field);
        % diff_no = relerr_percent_clipped(ens_field, true_field);
    % else
        % diff_no = ens_field - true_field;
    % end
    maxAbs_no = max(abs(diff_no(:)));
    maxAbs_global = maxAbs_no;
  
    % maxAbs_no = prctile(abs(diff_no(:)), 99);
    plotmodel(md_ens_1,'data',diff_no,'title',sprintf('no assimilation'), ...
        'subplot',[nrows,1,2],'caxis',[-maxAbs_no maxAbs_no],'colorbar','off');

    ax = gca; axs(2) = ax;
    ttl = ax.Title;
    ttl.FontSize   = 10;
    ttl.FontWeight = 'bold';   % or 'normal'
    ttl.Interpreter = 'tex';
    
    % ---- km axes (ticks shown in km) ----
    xt = get(ax,'XTick'); yt = get(ax,'YTick');
    set(ax,'XTickLabel', arrayfun(@(v) sprintf('%g', v./1000), xt, 'UniformOutput', false));
    set(ax,'YTickLabel', arrayfun(@(v) sprintf('%g', v./1000), yt, 'UniformOutput', false));
    % xlabel(ax,'x (km)','FontWeight','bold');
    ylabel(ax,'y (km)','FontWeight','bold','FontSize',16);
    
    % ---- panel letter inside upper-left ----
    % panel = sprintf('(%c)', 'a'+(idx-1));
    panel = sprintf('(%c_{%d})','b', panel_idx);
    text(ax, 0.02, 0.95, panel, 'Units','normalized', ...
        'FontWeight','bold', 'FontSize', 16, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'Color',lightGray);
    % prevent plotmodel objects appearing in legend
    set(ax.Children, 'HandleVisibility','off');

    % (c..): Assim - True
    nt=251;
    for idx = 1:nk
        k = k_array(idx);
        label_t = iff(k == nt-1, t(nt), t(k));
        [md_true_k, ~, md_ens_k] = setup_model_states(k, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
        ens_field = get_nested_field(md_ens_k, field);
        true_field = get_nested_field(md_true_k, field);
        % if contains(field,'geometry.bed')
            % diff_k = relative_error(ens_field, true_field);
            diff_k = signed_log_relerr(ens_field, true_field);
            % diff_k = relerr_percent_clipped(ens_field, true_field);
        % else
            % diff_k = ens_field - true_field;
        % end
        % diff_k = get_nested_field(md_ens_k, field) - get_nested_field(md_true_k, field);

        % maxAbs = max(abs(diff_k(:)));
        maxAbs_global = max(maxAbs_global, max(abs(diff_k(:))));
        maxAbs=maxAbs_global;
      
        % maxAbs = prctile(abs(diff_k(:)), 99);
        % label  = sprintf('\\bf(%c)', 'b'+idx);
        plotmodel(md_ens_k,'data',diff_k, ...
            'title',sprintf('after %s years of assimilation', fmt_years(label_t)), ...
            'subplot',[nrows,1,idx+2],'caxis',[-maxAbs maxAbs],'colorbar','off');

        ax = gca;
        ttl = ax.Title;
        ttl.FontSize   = 10;
        ttl.FontWeight = 'bold';   % or 'normal'
        ttl.Interpreter = 'tex';
        
        % ---- km axes (ticks shown in km) ----
        % axs(idx) = ax;
        axs(idx+2) = ax; 
        xt = get(ax,'XTick'); yt = get(ax,'YTick');
        set(ax,'XTickLabel', arrayfun(@(v) sprintf('%g', v./1000), xt, 'UniformOutput', false));
        set(ax,'YTickLabel', arrayfun(@(v) sprintf('%g', v./1000), yt, 'UniformOutput', false));
        % xlabel(ax,'x (km)','FontWeight','bold');
        ylabel(ax,'y (km)','FontWeight','bold','FontSize',16);
        
        % ---- panel letter inside upper-left ----
        % panel = sprintf('(%c)', 'a'+(idx-1));
        % panel_idx = 2;   % change as needed
        panel = sprintf('(%c_{%d})','c'+(idx-1), panel_idx);
        text(ax, 0.02, 0.95, panel, 'Units','normalized', ...
            'FontWeight','bold', 'FontSize', 16, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color',lightGray);

     
        % prevent plotmodel objects appearing in legend
        set(ax.Children, 'HandleVisibility','off');
    end

    % layout 
    axs = flipud(findall(gcf,'Type','axes'));
    gap = 0.03; top = 0.96; bottom = 0.08;
    avail = top-bottom - (nrows-1)*gap;
    height = avail/nrows;

    if height < 0.05
        fig = gcf;
        fig.Position(4) = fig.Position(4) * max(1, ceil(0.05/height));
        height = 0.05;
    end

    % for i = 1:nrows
    %     pos = [0.10, bottom+(nrows-i)*(height+gap), 0.70, height];
    %     set(axs(i),'Position',pos,'FontWeight','bold','LineWidth',1.2,'Box','on', ...
    %         'TickDir','out','Layer','top','FontSize',11,'TickLength',[0.005 0.005]);
    %     ylabel(axs(i),'y (km)','FontWeight','bold');
    %     if i < nrows, set(axs(i),'XTickLabel',[]);
    %     else, xlabel(axs(i),'x (km)','FontWeight','bold'); end
    % end

    for i = 1:nrows
        ax = axs(i);
    
        % i=1 should be the TOP panel (subplot does this already)
        % pos = [0.10, bottom+(nrows-i)*(height+gap), 0.70, height];
        pos = [0.10, bottom+(nrows-i)*(height+gap), 0.9, height];
        set(ax, 'Position', pos, ...
        'FontWeight','bold', ...
        'FontSize',15, ...
        'Box','on', ...
        'LineWidth',2.0, ...
        'TickDir','out', ...
        'TickLength',[0.004 0.004], ...
        'XGrid','off', ...
        'YGrid','off', ...
        'YMinorGrid','off', ...
        'XTickMode','manual', ...
        'YTickMode','manual');
        ttl = ax.Title;
        ttl.FontSize   = 15;
        ttl.FontWeight = 'bold';   % or 'normal'
        ttl.Interpreter = 'tex';
        % ---- FORCE tick locations (meters) ----
        % yt_km = 0:20:80;              % desired ticks in km
        % yt_m  = yt_km * 1000;         % convert to meters
        % 
        % set(ax, ...
        %     'YTick', yt_m, ...
        %     'YTickLabel', arrayfun(@num2str, yt_km, 'UniformOutput', false), ...
        %     'YTickMode','manual');
        yl = ax.YLim / 1000;  % km
        step = 40;            % km
        yt_km = ceil(yl(1)/step)*step : step : floor(yl(2)/step)*step;
        set(ax,'YTick',yt_km*1000,'YTickLabel',string(yt_km),'YTickMode','manual');

        % same idea for x if needed
        xt = get(ax,'XLim');
        xt_km = floor(xt(1)/1000/100)*100 : 100 : ceil(xt(2)/1000/100)*100;
        set(ax, ...
            'XTick', xt_km*1000, ...
            'XTickLabel', arrayfun(@num2str, xt_km, 'UniformOutput', false), ...
            'XTickMode','manual');
    
        ylabel(ax,'y (km)','FontWeight','bold','FontSize',16);
    
        if i < nrows
            set(ax,'XTickLabel',[]);
        else
            xlabel(ax,'x (km)','FontWeight','bold','FontSize',16);
            ax.XLabel.Units = 'normalized';
            ax.XLabel.Position(2) = -0.33;
        end

    end

    cb1 = colorbar(axs(1), 'Position',[colorbar_gap 0.68 0.015 0.16]);
    ylabel(cb1,[field_title units_str],'FontSize',15,'FontWeight','bold');
    colormap(axs(1), parula);
    % 
    % for i = 2:nrows, colormap(axs(i), redblue(256)); end
    % cb2 = colorbar(axs(end), 'Position',[colorbar_gap 0.25 0.015 0.40]);
    % % ylabel(cb2,'Difference','FontSize',12,'FontWeight','bold');
    % % if contains(field,'geometry.bed')
    %     ylabel(cb2,'Relative Error','FontSize',15,'FontWeight','bold');
    % else
    %     ylabel(cb2,['Difference' units_str],'FontSize',12,'FontWeight','bold');
    % end

    for i = 2:nrows
        colormap(axs(i), redblue(256));
        caxis(axs(i), [-maxAbs_global maxAbs_global]);
    end
    cb2 = colorbar(axs(end), 'Position',[colorbar_gap 0.25 0.015 0.40]);
    ylabel(cb2,'Relative Error','FontSize',15,'FontWeight','bold');
    % if contains(field, 'initialization.vel')
    %     ylabel(cb2,['\Delta |u|' units_str],'FontSize',15,'FontWeight','bold');
    % else
    %     ylabel(cb2,['\Delta ' field_title units_str],'FontSize',15,'FontWeight','bold');
    % end

    set(gcf,'Color','w');

    % ---- Save figure ----
    fname_base = sprintf('diff_%s', regexprep(field,'\.','_'));
    save_figure_300dpi(fname_base);
end

function plot_var_evolution(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md, field, field_title, units)
    
    global t nt label_t colorbar_gap
    % lightGray = [0.85 0.85 0.85];
    % lightGray = [0.93 0.69 0.13];
    lightGray='m';

    % if nargin < 12, field_title = field; end
    % if nargin < 13, units = ''; end
    units_str = iff(~isempty(units), [' (' units ')'], '');

    nk    = length(k_array);
    nrows = 2 + nk;

    figure('Position',[400 400 1100 (180 + 150*nrows)]); clf;
    axs = gobjects(nrows,1);   % <-- store ONLY the real panel axes

    % global color limits across (true+ens) at requested steps
    all_data = [];
    for k = [k_array(end), k_array]
        [md_true_tmp, ~, md_ens_tmp] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
    
        tmpT = get_nested_field(md_true_tmp, field);
        tmpE = get_nested_field(md_ens_tmp,  field);
    
        all_data = [all_data; tmpT(:); tmpE(:)];
    end
    cmin = min(all_data);
    cmax = max(all_data);
    clear all_data

    % (a) True at last snapshot
    [md_true_last, ~, ~] = setup_model_states(1, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);
    data_true = get_nested_field(md_true_last, field);
    k1 = 1;
    nt = 251;
    label_t = iff(k1 == nt-1, t(nt), t(k1));
    plotmodel(md_true_last,'data',data_true, ...
        'title',sprintf('True %s ', field_title), ...
        'subplot',[nrows,1,1],'caxis',[cmin cmax],'colorbar','off');

    ax = gca; axs(1) = ax;
    ttl = ax.Title;
    ttl.FontSize   = 10;
    ttl.FontWeight = 'bold';   % or 'normal'
    ttl.Interpreter = 'tex';
    
    % ---- km axes (ticks shown in km) ----
    xt = get(ax,'XTick'); yt = get(ax,'YTick');
    set(ax,'XTickLabel', arrayfun(@(v) sprintf('%g', v./1000), xt, 'UniformOutput', false));
    set(ax,'YTickLabel', arrayfun(@(v) sprintf('%g', v./1000), yt, 'UniformOutput', false));
    % xlabel(ax,'x (km)','FontWeight','bold');
    ylabel(ax,'y (km)','FontWeight','bold','FontSize',16);
    
    % ---- panel letter inside upper-left ----
    % panel = sprintf('(%c)', 'a'+(idx-1));
    panel_idx = 2;   % change as needed
    panel = sprintf('(%c_{%d})','a', panel_idx);
    text(ax, 0.02, 0.95, panel, 'Units','normalized', ...
        'FontWeight','bold', 'FontSize', 16, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'Color',lightGray);

    % (b) No assimilation (k=1) ensemble
    [~, ~, md_ens_1] = setup_model_states(1, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);
    data_ens = get_nested_field(md_ens_1, field);
    plotmodel(md_ens_1,'data',data_ens, ...
        'title',sprintf(' No assimilation %s', field_title), ...
        'subplot',[nrows,1,2],'caxis',[cmin cmax],'colorbar','off');

    ax = gca; axs(2) = ax;
    ttl = ax.Title;
    ttl.FontSize   = 10;
    ttl.FontWeight = 'bold';   % or 'normal'
    ttl.Interpreter = 'tex';
    
    % ---- km axes (ticks shown in km) ----
    xt = get(ax,'XTick'); yt = get(ax,'YTick');
    set(ax,'XTickLabel', arrayfun(@(v) sprintf('%g', v./1000), xt, 'UniformOutput', false));
    set(ax,'YTickLabel', arrayfun(@(v) sprintf('%g', v./1000), yt, 'UniformOutput', false));
    % xlabel(ax,'x (km)','FontWeight','bold');
    ylabel(ax,'y (km)','FontWeight','bold','FontSize',16);
    
    % ---- panel letter inside upper-left ----
    % panel = sprintf('(%c)', 'a'+(idx-1));
    panel = sprintf('(%c_{%d})','b', panel_idx);
    text(ax, 0.02, 0.95, panel, 'Units','normalized', ...
        'FontWeight','bold', 'FontSize', 16, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'Color',lightGray);
    % prevent plotmodel objects appearing in legend
    set(ax.Children, 'HandleVisibility','off');

    % (c..) Assim snapshots
    for idx = 1:nk
        k = k_array(idx);
        label_t = iff(k == nt-1, t(nt), t(k));
        [~, ~, md_ens_k] = setup_model_states(k, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
        data_ens = get_nested_field(md_ens_k, field);
        % label = sprintf('\\bf(%c)', 'b'+idx);
        plotmodel(md_ens_k,'data',data_ens, ...
            'title',sprintf('after %s years of assimilation', fmt_years(label_t)), ...
            'subplot',[nrows,1,idx+2],'caxis',[cmin cmax],'colorbar','off');
        ax = gca;
        ttl = ax.Title;
        ttl.FontSize   = 10;
        ttl.FontWeight = 'bold';   % or 'normal'
        ttl.Interpreter = 'tex';
        
        % ---- km axes (ticks shown in km) ----
        % axs(idx) = ax;
        axs(idx+2) = ax; 
        xt = get(ax,'XTick'); yt = get(ax,'YTick');
        set(ax,'XTickLabel', arrayfun(@(v) sprintf('%g', v./1000), xt, 'UniformOutput', false));
        set(ax,'YTickLabel', arrayfun(@(v) sprintf('%g', v./1000), yt, 'UniformOutput', false));
        % xlabel(ax,'x (km)','FontWeight','bold');
        ylabel(ax,'y (km)','FontWeight','bold','FontSize',16);
        
        % ---- panel letter inside upper-left ----
        % panel = sprintf('(%c)', 'a'+(idx-1));
        % panel_idx = 2;   % change as needed
        panel = sprintf('(%c_{%d})','c'+(idx-1), panel_idx);
        text(ax, 0.02, 0.95, panel, 'Units','normalized', ...
            'FontWeight','bold', 'FontSize', 16, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color',lightGray);

     
        % prevent plotmodel objects appearing in legend
        set(ax.Children, 'HandleVisibility','off');
    end

    % layout (your adaptive spacing)
    axs = flipud(findall(gcf,'Type','axes'));
    % gap = 0.02; top = 0.95; bottom = 0.08;
    gap = 0.03; top = 0.96; bottom = 0.08;
    avail = top-bottom - (nrows-1)*gap;
    height = avail/nrows;

    if height < 0.05
        fig = gcf;
        fig.Position(4) = fig.Position(4) * max(1, ceil(0.05/height));
        height = 0.05;
    end

    % for i = 1:nrows
    %     pos = [0.10, bottom+(nrows-i)*(height+gap), 0.70, height];
    %     set(axs(i),'Position',pos,'FontWeight','bold','LineWidth',1.2,'Box','on', ...
    %         'TickDir','out','Layer','top','FontSize',11,'TickLength',[0.005 0.005]);
    %     ylabel(axs(i),'y (km)','FontWeight','bold');
    %     if i < nrows, set(axs(i),'XTickLabel',[]);
    %     else, xlabel(axs(i),'x (km)','FontWeight','bold'); end
    % end

    for i = 1:nrows
        ax = axs(i);
    
        % i=1 should be the TOP panel (subplot does this already)
        % pos = [0.10, bottom+(nrows-i)*(height+gap), 0.70, height];
        pos = [0.10, bottom+(nrows-i)*(height+gap), 0.9, height];
        set(ax, 'Position', pos, ...
        'FontWeight','bold', ...
        'FontSize',15, ...
        'Box','on', ...
        'LineWidth',2.0, ...
        'TickDir','out', ...
        'TickLength',[0.004 0.004], ...
        'XGrid','off', ...
        'YGrid','off', ...
        'YMinorGrid','off', ...
        'XTickMode','manual', ...
        'YTickMode','manual');
        ttl = ax.Title;
        ttl.FontSize   = 15;
        ttl.FontWeight = 'bold';   % or 'normal'
        ttl.Interpreter = 'tex';
        % ---- FORCE tick locations (meters) ----
        % yt_km = 0:20:80;              % desired ticks in km
        % yt_m  = yt_km * 1000;         % convert to meters
        % 
        % set(ax, ...
        %     'YTick', yt_m, ...
        %     'YTickLabel', arrayfun(@num2str, yt_km, 'UniformOutput', false), ...
        %     'YTickMode','manual');
        yl = ax.YLim / 1000;  % km
        step = 40;            % km
        yt_km = ceil(yl(1)/step)*step : step : floor(yl(2)/step)*step;
        set(ax,'YTick',yt_km*1000,'YTickLabel',string(yt_km),'YTickMode','manual');

        % same idea for x if needed
        xt = get(ax,'XLim');
        xt_km = floor(xt(1)/1000/100)*100 : 100 : ceil(xt(2)/1000/100)*100;
        set(ax, ...
            'XTick', xt_km*1000, ...
            'XTickLabel', arrayfun(@num2str, xt_km, 'UniformOutput', false), ...
            'XTickMode','manual');
    
        ylabel(ax,'y (km)','FontWeight','bold','FontSize',16);
    
        if i < nrows
            set(ax,'XTickLabel',[]);
        else
            xlabel(ax,'x (km)','FontWeight','bold','FontSize',16);
            ax.XLabel.Units = 'normalized';
            ax.XLabel.Position(2) = -0.33;
        end

    end

    for i = 1:nrows, colormap(axs(i), parula); end
    cb = colorbar(axs(end), 'Position',[colorbar_gap 0.25 0.015 0.45]);
    ylabel(cb,[field_title units_str],'FontSize',15,'FontWeight','bold');

    set(gcf,'Color','w');

    % ---- Save figure ----
    fname_base = sprintf('evol_%s', regexprep(field,'\.','_'));
    save_figure_300dpi(fname_base);
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

    % =========================
    % RELATIVE ERROR
    % =========================
    % Stabilizer avoids blowing up where true≈0
    eps0 = 0.01 * max(abs(data_true(:)));   % 1% of max(true) is a good default
    if eps0 == 0
        eps0 = 1; % fallback if true field is exactly zero everywhere
    end

    diff_noassim = (data_nurged - data_true) ./ (abs(data_true) + eps0);
    diff_assim   = (data_ens    - data_true) ./ (abs(data_true) + eps0);

    % --- Limits for the field panels (True/Wrong/Assim) ---
    cmin = min([data_true(:); data_nurged(:); data_ens(:)]);
    cmax = max([data_true(:); data_nurged(:); data_ens(:)]);

    % --- Robust limits for relative-error color axis (avoid 1–2 spikes) ---
    allerr = [diff_noassim(:); diff_assim(:)];
    allerr = allerr(isfinite(allerr));
    if isempty(allerr)
        maxAbs = 1;
    else
        maxAbs = prctile(abs(allerr), 99.5);   % robust: shows structure
        if maxAbs == 0, maxAbs = 1; end
    end

    figure('Position',[100 100 1000 820]); clf;

    % 1) True
    plotmodel(md_true,'data',data_true,'title',['(a) True ' field_title], ...
        'subplot',[4,1,1],'caxis',[cmin cmax],'colorbar','off');

    % 2) Wrong
    plotmodel(md_nurged,'data',data_nurged,'title',['(b) No assimilation ' field_title], ...
        'subplot',[4,1,2],'caxis',[cmin cmax],'colorbar','off');

    % 3) Assimilated
    plotmodel(md_ens,'data',data_ens,'title',['(c) Assimilated ' field_title], ...
        'subplot',[4,1,3],'caxis',[cmin cmax],'colorbar','off');

    % 4) Relative error (Assim − True)/(|True|+eps0)
    plotmodel(md_ens,'data',diff_assim, ...
        'title',['(d) Relative error: (Assim − True) / (|True| + \epsilon)'], ...
        'subplot',[4,1,4],'caxis',[-maxAbs maxAbs],'colorbar','off');

    % --- Axes layout ---
    axs = flipud(findall(gcf,'Type','axes'));   % 1..4 top->bottom
    gap = -0.255; top = 0.94; bottom = 0.08;
    height = (top-bottom - 3*gap)/4;

    for i = 1:4
        pos = [0.10, bottom+(4-i)*(height+gap), 0.70, height];
        set(axs(i),'Position',pos, ...
            'FontWeight','bold','LineWidth',1.5,'Box','on', ...
            'TickDir','out','TickLength',[0.005 0.005], ...
            'Layer','top','FontSize',11);
        ylabel(axs(i),'y (km)','FontSize',12,'FontWeight','bold');
        if i < 4
            set(axs(i),'XTickLabel',[]);
        else
            xlabel(axs(i),'x (km)','FontSize',12,'FontWeight','bold');
        end
    end

    % % --- First colorbar (top 3) ---
    % for i = 1:3
    %     colormap(axs(i), cmap);
    %     caxis(axs(i), [cmin cmax]);
    % end
    % cb1 = colorbar(axs(2),'Position',[0.83 0.415 0.025 0.35]);
    % static_field = regexprep(field_title,'\s+after.*','');
    % ylabel(cb1,[static_field units_str],'FontSize',13,'FontWeight','bold');
    % cb1.FontSize = 11;
    % set(cb1,'Box','on','LineWidth',1.2);
    % 
    % % --- Second colorbar (relative error) ---
    % ax_diff = axs(4);
    % colormap(ax_diff, redblue(256));
    % caxis(ax_diff,[-maxAbs maxAbs]);
    % pos_diff = get(ax_diff,'Position');
    % cb2 = colorbar(ax_diff,'Position',[0.83 pos_diff(2)+0.14 0.025 pos_diff(4)-0.28]);
    % ylabel(cb2,'Relative error','FontSize',13,'FontWeight','bold');
    % cb2.FontSize = 11;
    % set(cb2,'Box','on','LineWidth',1.2);
    % 
    % set(gcf,'Color','w');

    % --- Adaptive colorbars for triptych ---
    % Top 3 panels share cb1
    for i = 1:3
        colormap(axs(i), cmap);
        caxis(axs(i), [cmin cmax]);
    end
    
    % Bottom panel uses cb2
    colormap(axs(4), redblue(256));
    caxis(axs(4), [-maxAbs maxAbs]);
    
    pos_top3_top    = axs(1).Position;   % top of panel (a)
    pos_top3_bottom = axs(3).Position;   % bottom of panel (c)
    pos_diff        = axs(4).Position;   % panel (d)
    
    gapx = 0.02;
    cb_w = 0.022;
    
    cb_x = pos_top3_top(1) + pos_top3_top(3) + gapx;
    
    % cb1 spans panels (a)-(c)
    cb1_y = pos_top3_bottom(2);
    cb1_h = (pos_top3_top(2) + pos_top3_top(4)) - pos_top3_bottom(2);
    cb1 = colorbar(axs(2), 'Position', [cb_x cb1_y cb_w cb1_h]);
    static_field = regexprep(field_title,'\s+after.*','');
    ylabel(cb1, [static_field units_str], 'FontSize',13,'FontWeight','bold');
    cb1.FontSize = 11;
    set(cb1,'Box','on','LineWidth',1.2);
    
    % cb2 spans only panel (d)
    cb2_y = pos_diff(2);
    cb2_h = pos_diff(4);
    cb2 = colorbar(axs(4), 'Position', [cb_x cb2_y cb_w cb2_h]);
    ylabel(cb2,'Relative error','FontSize',13,'FontWeight','bold');
    cb2.FontSize = 11;
    set(cb2,'Box','on','LineWidth',1.2);

    set(gcf,'Color','w');

    % ---- Save figure ----
    fname_base = sprintf('triptych_%s', regexprep(field,'\.','_'));
    save_figure_300dpi(fname_base);
end

function rel = relative_error(a, b)
% Compute (a-b)/max(|b|, eps) safely
    eps0 = 1e-6 * max(abs(b(:)));   % scale-aware stabilization
    rel  = (a - b) ./ max(abs(b), eps0);
end

function e_log = signed_log_relerr(x, xtrue)
    eps0  = 1e-3 * prctile(abs(xtrue(:)), 95);   % robust epsilon (NOT max)
    e     = (x - xtrue) ./ (abs(xtrue) + eps0);
    e_log = sign(e) .* log10(1 + abs(e));
end

function e_pct = relerr_percent_clipped(x, xtrue)
    eps0  = 1e-3 * prctile(abs(xtrue(:)), 95);
    e_pct = 100 * (x - xtrue) ./ (abs(xtrue) + eps0);

    % robust clip so a few points don’t dominate
    cap = prctile(abs(e_pct(:)), 99);
    e_pct = max(min(e_pct, cap), -cap);
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

function [gl_mid, dist] = compute_gl_midpoints_a( ...
    k_array, dt, ...
    model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md)
% Compute TRUE, wrong, assimilated grounding line midpoint (arc-length midpoint of longest 0-contour)
% for each k in k_array. Also stores a k->index map.

    nk = numel(k_array);
    gl_mid.k_array = k_array(:);
    gl_mid.x_true = nan(nk,1); gl_mid.y_true = nan(nk,1);
    gl_mid.x_nurged = nan(nk,1); gl_mid.y_nurged = nan(nk,1);
    gl_mid.x_ens = nan(nk,1); gl_mid.y_ens = nan(nk,1);

    dist.no = nan(nk,1); dist.as = nan(nk,1);
    
    gl_mid.k_to_idx = containers.Map('KeyType','double','ValueType','double');

    % contour grid resolution (match your GL plot)
    Nx = 420; Ny = 70;

    for idx = 1:nk
        k = k_array(idx);
        gl_mid.k_to_idx(k) = idx;

        [md_true_k, md_nurged_k, md_ens_k] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);

        x = md_true_k.mesh.x(:);
        y = md_true_k.mesh.y(:);
        phi_true = md_true_k.mask.ocean_levelset(:);
        phi_nurged = md_nurged_k.mask.ocean_levelset(:);
        phi_ens = md_ens_k.mask.ocean_levelset(:);

        xg = linspace(min(x), max(x), Nx);
        yg = linspace(min(y), max(y), Ny);
        [Xg, Yg] = meshgrid(xg, yg);

        F_true = scatteredInterpolant(x, y, phi_true, 'linear', 'nearest');
        F_nurged = scatteredInterpolant(x, y, phi_nurged, 'linear', 'nearest');
        F_ens = scatteredInterpolant(x, y, phi_ens, 'linear', 'nearest');
        Phi_t = F_true(Xg, Yg);
        Phi_n = F_nurged(Xg, Yg);
        Phi_e = F_ens(Xg, Yg);

        y_center = 0.5*(min(y) + max(y));  % channel centerline
        [xct, yct] = gl_centerline_point_from_levelset_grid(Xg, Yg, Phi_t, y_center);
        [xcn, ycn] = gl_centerline_point_from_levelset_grid(Xg, Yg, Phi_n, y_center);
        [xce, yce] = gl_centerline_point_from_levelset_grid(Xg, Yg, Phi_e, y_center);
        gl_mid.x_true(idx) = xct; gl_mid.y_true(idx) = yct;
        gl_mid.x_nurged(idx) = xcn; gl_mid.y_nurged(idx) = ycn;
        gl_mid.x_ens(idx) = xce; gl_mid.y_ens(idx) = yce;

        % distance
        dist.no(idx) = abs(xct-xcn);
        dist.as(idx) = abs(xct - xce);
        
    end

    figure('Position',[200 200 950 360]); clf;
    plot(k_array, dist.no, 'r-', 'LineWidth',1.8); hold on
    plot(k_array, dist.as, 'c-', 'LineWidth',1.8); hold off
    grid on
    xlabel('Time (years)','FontWeight','bold')
    ylabel('Absolute GL distance from True GL','FontWeight','bold')
    % title('\bfWindowed grounding-line position error (bend region)','Interpreter','tex')
    legend({'No assimilation vs True','Assimilated vs True'}, 'Location','best')
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

function out = compute_point_rmse_at_gl_mid( ...
    k_array, dt, fields, gl_mid, ...
    model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md)
% RMSE over k_array of point-sampled errors at TRUE GL midpoint.

    nk = numel(k_array);
    nf = numel(fields);

    rmse_no = nan(nf,1);
    rmse_as = nan(nf,1);

    % store optional timeseries
    err_no = nan(nk,nf);
    err_as = nan(nk,nf);

    for idx = 1:nk
        k = k_array(idx);

        [md_true_k, md_nurged_k, md_ens_k] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);

        xq = gl_mid.x(idx);
        yq = gl_mid.y(idx);

        for j = 1:nf
            f = fields{j};

            vt = sample_nodal_field_at_point(md_true_k, f, xq, yq);
            vn = sample_nodal_field_at_point(md_nurged_k, f, xq, yq);
            va = sample_nodal_field_at_point(md_ens_k, f, xq, yq);

            err_no(idx,j) = vn - vt;
            err_as(idx,j) = va - vt;
        end
    end

    for j = 1:nf
        e1 = err_no(:,j); e1 = e1(isfinite(e1));
        e2 = err_as(:,j); e2 = e2(isfinite(e2));

        rmse_no(j) = sqrt(mean(e1.^2));
        rmse_as(j) = sqrt(mean(e2.^2));
    end

    T = table(fields(:), rmse_no, rmse_as, ...
        'VariableNames', {'Field','RMSE_NoAssim','RMSE_Assimilated'});

    out.table = T;
    out.err_noassim = err_no;
    out.err_assim   = err_as;
end

function val = sample_nodal_field_at_point(md, field, xq, yq)
% Sample a nodal field at (xq,yq) using scatteredInterpolant.

    data = get_nested_field(md, field);
    data = data(:);

    x = md.mesh.x(:);
    y = md.mesh.y(:);

    F = scatteredInterpolant(x, y, data, 'linear', 'nearest');
    val = F(xq, yq);
end

function plot_gl_mid_marker(ax, gl_mid, k)
% Plot stored marker for given k if available. Safe no-op if missing.

    if isempty(gl_mid) || ~isfield(gl_mid,'k_to_idx'), return; end
    if ~isKey(gl_mid.k_to_idx, k), return; end

    idx = gl_mid.k_to_idx(k);
    if ~isfinite(gl_mid.x(idx)) || ~isfinite(gl_mid.y(idx)), return; end

    hold(ax,'on');
    plot(ax, gl_mid.x(idx), gl_mid.y(idx), 'wo', ...
        'MarkerFaceColor','w','MarkerSize',6,'LineWidth',1.5, ...
        'HandleVisibility','off');
    hold(ax,'off');
end

function s = fmt_years(t)
    if abs(t - round(t)) < 1e-10
        s = sprintf('\\bf%d', round(t));   % integer, no decimal
    else
        s = sprintf('\\bf%.1f', t);        % one decimal
    end
end

function outdir = ensure_figdir()
% Create a figures folder.
    scriptdir = fileparts(mfilename('fullpath'));
    if isempty(scriptdir)
        scriptdir = pwd; % fallback
    end
    outdir = fullfile(scriptdir, 'figures');
    if ~exist(outdir,'dir'), mkdir(outdir); end
end

function save_figure_300dpi(fname_base)
% Save current figure as 300-dpi PNG in figures.
    outdir = ensure_figdir();
    fname = fullfile(outdir, [fname_base '.png']);
    drawnow;
    exportgraphics(gcf, fname, 'Resolution', 300);
end

function out = compute_rmse_timeseries_allt(dt, t, model_true_state, model_nurged_state, ensemble_vec_mean, md, nvar, varargin)

    p = inputParser;
    p.addParameter('mask_freeze_year', 0, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('useHpos', true, @(x) islogical(x) && isscalar(x));
    p.parse(varargin{:});
    mask_freeze_year = p.Results.mask_freeze_year;
    useHpos          = p.Results.useHpos;

    [nd, nt] = size(model_true_state);
    hdim = nd / nvar;

    I_h  = 1:hdim;
    I_vx = 2*hdim+1:3*hdim;
    I_vy = 3*hdim+1:4*hdim;
    I_fc = 5*hdim+1:6*hdim;

    % Fixed grounded mask from TRUTH at freeze time
    kfreeze = round(mask_freeze_year/dt) + 1;
    kfreeze = max(1, min(nt, kfreeze));
    mask_grounded = grounded_mask_from_state(model_true_state, md, hdim, kfreeze, useHpos);

    rmse_h_no  = nan(nt,1);  rmse_h_as  = nan(nt,1);
    rmse_u_no  = nan(nt,1);  rmse_u_as  = nan(nt,1);
    rmse_c_no  = nan(nt,1);  rmse_c_as  = nan(nt,1);

    for k = 1:nt-1
        % --- Thickness ---
        ht = model_true_state(I_h, k);
        hn = model_nurged_state(I_h, k);
        ha = ensemble_vec_mean(I_h, k);

        rmse_h_no(k) = rmse_vec(hn, ht);
        rmse_h_as(k) = rmse_vec(ha, ht);

        % --- Velocity magnitude ---
        vxt = model_true_state(I_vx, k); vyt = model_true_state(I_vy, k);
        vxn = model_nurged_state(I_vx, k); vyn = model_nurged_state(I_vy, k);
        vxa = ensemble_vec_mean(I_vx, k); vya = ensemble_vec_mean(I_vy, k);

        ut = hypot(vxt, vyt);
        un = hypot(vxn, vyn);
        ua = hypot(vxa, vya);

        rmse_u_no(k) = rmse_vec(un, ut);
        rmse_u_as(k) = rmse_vec(ua, ut);

        % --- Friction (grounded-only, fixed mask) ---
        ct = model_true_state(I_fc, k);
        cn = model_nurged_state(I_fc, k);
        ca = ensemble_vec_mean(I_fc, k);

        rmse_c_no(k) = rmse_vec(cn, ct, mask_grounded);
        rmse_c_as(k) = rmse_vec(ca, ct, mask_grounded);
    end

    out.time_years = (0:nt-1)' * dt;
    out.rmse_h_no  = rmse_h_no;  out.rmse_h_as = rmse_h_as;
    out.rmse_u_no  = rmse_u_no;  out.rmse_u_as = rmse_u_as;
    out.rmse_c_no  = rmse_c_no;  out.rmse_c_as = rmse_c_as;
    out.mask_grounded = mask_grounded;
end

function r = rmse_masked(a, b, mask_na, mask_true)
    % b is true
    a = a(:); b = b(:);
    mask_na = logical(mask_na(:));
    mask_true = logical(mask_true(:));
    % good_na = mask_na & isfinite(a) & isfinite(b);
    % good_true = mask_true & isfinite(a) & isfinite(b);
    % good = mask_na & isfinite(a) & isfinite(b);
    good = mask_na & isfinite(a) & isfinite(b);
    if ~any(good)
        r = NaN;
        return;
    end
    % if ~any(good_true)
    %     r = NaN;
    %     return;
    % end
    % d = a(good_na) - b(good_true);
    d = a(good) - b(good);
    r = sqrt(mean(d.^2));
end


function r = rmse_w(a, b)
    % b is true
    a = a(:); b = b(:);
    d = a - b;
    r = sqrt(mean(d.^2));
end

function r = rmse_masked_pair(a, b, mask_a, mask_b, domain)
%RMSE_MASKED_PAIR RMSE between a and b on a chosen mask domain.
% domain: 'a' | 'b' | 'intersection'

    if nargin < 5 || isempty(domain)
        domain = 'c';
    end

    a = a(:); b = b(:);
    mask_a = logical(mask_a(:));
    mask_b = logical(mask_b(:));

    switch lower(domain)
        case 'a'
            % mask_size = size(mask_a);
            good = mask_a & isfinite(a) & isfinite(b);
            % good = mask_a;
        case 'b'
            good = mask_b & isfinite(a) & isfinite(b);
            % good = mask_b;
        case 'c'
            % good = (mask_a & mask_b) & isfinite(a) & isfinite(b);
            good = mask_a & mask_b;
        otherwise
            error('domain must be ''a'', ''b'', or ''c''.');
    end

    if ~any(good)
        r = NaN;
        return;
    end

    d = a(good) - b(good);
    r = sqrt(mean(d.^2));
end


function out = compute_rmse_timeseries(k_array, dt, t, model_true_state, model_nurged_state, ensemble_vec_mean, md_true, md_nurged, md_ens, md, field)
% compute_rmse_timeseries
% Clean, plot-consistent RMSE time series for:
%   (1) Thickness RMSE on TRUE grounded ice (ocean_levelset>0 & H>0)
%   (2) Velocity  RMSE on TRUE floating  ice (ocean_levelset<0 & H>0)
%   (3) Friction  RMSE on TRUE grounded ice (ocean_levelset>0 & H>0)
%   (4) Grounding-line "star distance": |x_GL_true - x_GL_other| along centerline
%
% Notes:
% - Uses TRUE mask as the reference domain each time step (matches “RMSE at the stars” idea).
% - Uses the same centerline GL point you use for plotting (via levelset=0 crossing).
% - `field` kept for API compatibility (not used here, since you plot 4 panels anyway).

    % -------------------------------
    % Time / step indexing
    % -------------------------------
    % if isempty(k_array)
        % nt = size(model_true_state, 2);
        nt = 251;
        kvec = 1:nt-1;
    % else
    %     kvec = k_array(:)';              % enforce row
    % end

    % Map each k to a time value (robust if t is length nt)
    if numel(t) >= max(kvec)
        time_vec = t(kvec);
    else
        % fallback: assume t is already aligned with kvec or is scalar dt-based
        time_vec = (kvec-1) * dt;
    end

    nk = numel(kvec);

    % -------------------------------
    % Preallocate outputs
    % -------------------------------
    rmse_h_no   = nan(nk,1);  rmse_h_as   = nan(nk,1);
    rmse_vel_no = nan(nk,1);  rmse_vel_as = nan(nk,1);
    rmse_c_no   = nan(nk,1);  rmse_c_as   = nan(nk,1);

    % grounded ice only exluding grounding line (x<=300km)
    rmse_h_no_g   = nan(nk,1);  rmse_h_as_g   = nan(nk,1);
    rmse_vel_no_g = nan(nk,1);  rmse_vel_as_g = nan(nk,1);
    rmse_c_no_g   = nan(nk,1);  rmse_c_as_g   = nan(nk,1);

    % whole domain
    rmse_h_no_w   = nan(nk,1);  rmse_h_as_w   = nan(nk,1);
    rmse_vel_no_w = nan(nk,1);  rmse_vel_as_w = nan(nk,1);
    rmse_c_no_w   = nan(nk,1);  rmse_c_as_w   = nan(nk,1);

    gl_mid.k         = kvec(:);
    gl_mid.x_true    = nan(nk,1); gl_mid.y_true    = nan(nk,1);
    gl_mid.x_nurged  = nan(nk,1); gl_mid.y_nurged  = nan(nk,1);
    gl_mid.x_ens     = nan(nk,1); gl_mid.y_ens     = nan(nk,1);

    dist_no = nan(nk,1);
    dist_as = nan(nk,1);

    % -------------------------------
    % Contour grid (match GL plotting)
    % -------------------------------
    Nx = 420; Ny = 70;

    % Build grid ONCE (mesh is fixed)
    % We'll grab x/y bounds from the first step we load.
    k0 = kvec(1);
    [md_true_0, md_nurged_0, md_ens_0] = setup_model_states(k0, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);

    x = md_true_0.mesh.x(:);
    y = md_true_0.mesh.y(:);

    xg = linspace(min(x), max(x), Nx);
    yg = linspace(min(y), max(y), Ny);
    [Xg, Yg] = meshgrid(xg, yg);

    y_center = 0.5*(min(y) + max(y));  % channel centerline

    % freeze the reference mask:
    % pick the index corresponding to year 24
    % k_obs_end = find(abs(t - 24.0) == min(abs(t - 24.0)), 1);
     k_end =k_array(end);

     [md_true_end, md_nurged_end, md_ens_end] = setup_model_states(k_end, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);

    x = md_true_end.mesh.x;
    pos_g = find(x<=300e3);
    
    % --- Per-state masks (TRUE / NURGED / ENS) ---
    H_true = md_true_end.geometry.thickness(:);
    H_n    = md_nurged_end.geometry.thickness(:);
    H_e    = md_ens_end.geometry.thickness(:);
    
    %
    g_true = (md_true_end.mask.ocean_levelset(:)   > 0);
    g_n    = (md_nurged_end.mask.ocean_levelset(:) > 0);
    g_e    = (md_ens_end.mask.ocean_levelset(:)    > 0);

    % 
    ice_true = (H_true > 0);
    ice_n    = (H_n    > 0);
    ice_e    = (H_e    > 0);
    
    % grounded without GL
    g_trueg = (md_true_end.mask.ocean_levelset(pos_g)   > 0);
    g_ng    = (md_nurged_end.mask.ocean_levelset(pos_g) > 0);
    g_eg    = (md_ens_end.mask.ocean_levelset(pos_g)    > 0);

    H_trueg = md_true_end.geometry.thickness(pos_g);
    H_ng    = md_nurged_end.geometry.thickness(pos_g);
    H_eg    = md_ens_end.geometry.thickness(pos_g);
    
    ice_trueg = (H_trueg > 0);
    ice_ng    = (H_ng    > 0);
    ice_eg    = (H_eg    > 0);
        
    % -------------------------------
    % Main loop
    % -------------------------------
    for ii = 1:nk
        k = kvec(ii);

        [md_true_k, md_nurged_k, md_ens_k] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);

        % --- Per-state masks (TRUE / NURGED / ENS) ---
        % H_true = md_true_k.geometry.thickness(:);
        % H_n    = md_nurged_k.geometry.thickness(:);
        % H_e    = md_ens_k.geometry.thickness(:);
        % % 
        % g_true = (md_true_k.mask.ocean_levelset(:)   > 0);
        % g_n    = (md_nurged_k.mask.ocean_levelset(:) > 0);
        % g_e    = (md_ens_k.mask.ocean_levelset(:)    > 0);
        % 
        % f_true = (md_true_k.mask.ocean_levelset(:)   < 0);
        % f_n    = (md_nurged_k.mask.ocean_levelset(:) < 0);
        % f_e    = (md_ens_k.mask.ocean_levelset(:)    < 0);
        % % 
        % ice_true = (H_true > 0);
        % ice_n    = (H_n    > 0);
        % ice_e    = (H_e    > 0);
        
        % --- Evaluation domains (INTERSECTION) ---
        % Thickness + friction: grounded ice intersection
        % mask_grounded_no = g_true & g_n & ice_true & ice_n;   % true ∩ nurged
        % mask_grounded_as = g_true & g_e & ice_true & ice_e;   % true ∩ ens
        mask_grounded_no = g_n & ice_n; 
        mask_grounded_as = g_e & ice_e;
        mask_grounded_true = g_true & ice_true;

        mask_grounded_no_c = g_n & ice_n; 
        mask_grounded_as_c = g_e & ice_e;
        mask_grounded_true_c = g_true & ice_true;

        % mask grounded exluding GL
        mask_grounded_nog = g_ng & ice_ng; 
        mask_grounded_asg = g_eg & ice_eg;
        mask_grounded_trueg = g_trueg & ice_trueg;

        % % Velocity: floating ice intersection
        % mask_floating_no = f_true & f_n & ice_true & ice_n;
        % mask_floating_as = f_true & f_e & ice_true & ice_e;

        % ============================================================
        % (1) Thickness RMSE on TRUE grounded ice
        % ============================================================
        H_true = md_true_k.geometry.thickness(:);
        H_nurged = md_nurged_k.geometry.thickness(:);
        H_ens    = md_ens_k.geometry.thickness(:);

        % rmse_h_no(ii) = rmse_masked(H_nurged, H_true, mask_grounded_no, mask_grounded_true);
        % rmse_h_as(ii) = rmse_masked(H_ens,    H_true, mask_grounded_as, mask_grounded_true);
        rmse_h_no(ii) = rmse_masked_pair(H_nurged, H_true, mask_grounded_no, mask_grounded_true,'b');
        rmse_h_as(ii) = rmse_masked_pair(H_ens,    H_true, mask_grounded_as, mask_grounded_true,'b');
        
        % H_nurgedg = md_nurged_k.geometry.thickness(pos_g); H_ensg    = md_ens_k.geometry.thickness(pos_g); 
        % H_trueg = md_true_k.geometry.thickness(pos_g);
        rmse_h_no_g(ii) = rmse_masked_pair(H_nurged(pos_g), H_true(pos_g), mask_grounded_no(pos_g), mask_grounded_true(pos_g),'b');
        rmse_h_as_g(ii) = rmse_masked_pair(H_ens(pos_g),    H_true(pos_g), mask_grounded_as(pos_g), mask_grounded_true(pos_g),'b');

        rmse_h_no_w(ii) = rmse_w(H_nurged, H_true);
        rmse_h_as_w(ii) = rmse_w(H_ens,    H_true);

        % ============================================================
        % (2) Velocity RMSE on TRUE floating ice
        % ============================================================
        V_true   = md_true_k.initialization.vel(:);
        V_nurged = md_nurged_k.initialization.vel(:);
        V_ens    = md_ens_k.initialization.vel(:);
        
        % rmse_vel_no(ii) = rmse_masked(v_nurged, v_true, mask_floating_no);
        % rmse_vel_as(ii) = rmse_masked(v_ens,    v_true, mask_floating_as);
        % rmse_vel_no(ii) = rmse_masked(v_nurged, v_true, mask_grounded_no, mask_grounded_true);
        % rmse_vel_as(ii) = rmse_masked(v_ens,    v_true, mask_grounded_as, mask_grounded_true);
        rmse_vel_no(ii) = rmse_masked_pair(V_nurged, V_true, mask_grounded_no, mask_grounded_true,'b');
        rmse_vel_as(ii) = rmse_masked_pair(V_ens,    V_true, mask_grounded_as, mask_grounded_true,'b');
        
        % v_nurgedg = md_nurged_k.initialization.vel(pos_g); v_ensg    = md_ens_k.initialization.vel(pos_g);
        % V_trueg   = md_true_k.initialization.vel(pos_g);
        rmse_vel_no_g(ii) = rmse_masked_pair(V_nurged(pos_g), V_true(pos_g), mask_grounded_no(pos_g), mask_grounded_true(pos_g),'b');
        rmse_vel_as_g(ii) = rmse_masked_pair(V_ens(pos_g),    V_true(pos_g), mask_grounded_as(pos_g), mask_grounded_true(pos_g),'b');

        rmse_vel_no_w(ii) = rmse_w(V_nurged, V_true);
        rmse_vel_as_w(ii) = rmse_w(V_ens,    V_true);


        % ============================================================
        % (3) Friction coefficient RMSE on TRUE grounded ice
        % ============================================================
        C_true   = md_true_k.friction.coefficient(:);
        C_nurged = md_nurged_k.friction.coefficient(:);
        C_ens    = md_ens_k.friction.coefficient(:);

        % rmse_c_no(ii) = rmse_masked(C_nurged, C_true, mask_grounded_no, mask_grounded_true);
        % rmse_c_as(ii) = rmse_masked(C_ens,    C_true, mask_grounded_as, mask_grounded_true);
        rmse_c_no(ii) = rmse_masked_pair(C_nurged, C_true,  mask_grounded_no_c, mask_grounded_true_c,'b');
        rmse_c_as(ii) = rmse_masked_pair(C_ens,    C_true, mask_grounded_as_c, mask_grounded_true_c,'b');
        
        % C_nurgedg = md_nurged_k.friction.coefficient(pos_g); C_ensg    = md_ens_k.friction.coefficient(pos_g);
        % C_trueg   = md_true_k.friction.coefficient(pos_g);
        rmse_c_no_g(ii) = rmse_masked_pair(C_nurged(pos_g), C_true(pos_g),  mask_grounded_no(pos_g), mask_grounded_true(pos_g),'b');
        rmse_c_as_g(ii) = rmse_masked_pair(C_ens(pos_g),    C_true(pos_g), mask_grounded_as(pos_g), mask_grounded_true(pos_g),'b');

        rmse_c_no_w(ii) = rmse_w(C_nurged, C_true);
        rmse_c_as_w(ii) = rmse_w(C_ens,    C_true);

        % ============================================================
        % (4) Grounding line centerline “star distance”
        %     Use levelset grid crossing to get a single (x,y) per state.
        % ============================================================
        phi_true   = md_true_k.mask.ocean_levelset(:);
        phi_nurged = md_nurged_k.mask.ocean_levelset(:);
        phi_ens    = md_ens_k.mask.ocean_levelset(:);

        % Interpolate to plotting grid
        F_true   = scatteredInterpolant(x, y, phi_true,   'linear', 'nearest');
        F_nurged = scatteredInterpolant(x, y, phi_nurged, 'linear', 'nearest');
        F_ens    = scatteredInterpolant(x, y, phi_ens,    'linear', 'nearest');

        Phi_t = F_true(Xg, Yg);
        Phi_n = F_nurged(Xg, Yg);
        Phi_e = F_ens(Xg, Yg);

        xprev_true = NaN;
        xprev_no   = NaN;
        xprev_as   = NaN;

        % [xct, yct] = gl_centerline_point_from_levelset_grid(Xg, Yg, Phi_t, y_center);
        % [xcn, ycn] = gl_centerline_point_from_levelset_grid(Xg, Yg, Phi_n, y_center);
        % [xce, yce] = gl_centerline_point_from_levelset_grid(Xg, Yg, Phi_e, y_center);
        [xct, yct] = gl_centerline_point_from_levelset_grid_track(Xg, Yg, Phi_t, y_center, xprev_true);
        [xcn, ycn] = gl_centerline_point_from_levelset_grid_track(Xg, Yg, Phi_n, y_center, xprev_no);
        [xce, yce] = gl_centerline_point_from_levelset_grid_track(Xg, Yg, Phi_e, y_center, xprev_as);
        
        if isfinite(xct), xprev_true = xct; end
        if isfinite(xcn), xprev_no   = xcn; end
        if isfinite(xce), xprev_as   = xce; end

        gl_mid.x_true(ii)   = xct; gl_mid.y_true(ii)   = yct;
        gl_mid.x_nurged(ii) = xcn; gl_mid.y_nurged(ii) = ycn;
        gl_mid.x_ens(ii)    = xce; gl_mid.y_ens(ii)    = yce;

        % Absolute distance between stars along centerline (x-distance)
        if ~isnan(xct) && ~isnan(xcn), dist_no(ii) = abs(xct - xcn); end
        if ~isnan(xct) && ~isnan(xce), dist_as(ii) = abs(xct - xce); end

        % Hard sanity guard: if it jumps unrealistically, drop it
        if ii>1 && isfinite(dist_as(ii)) && isfinite(dist_as(ii-1))
            if dist_as(ii) > dist_as(ii-1) + 5e4   % e.g. jump > 50 km
                dist_as(ii) = NaN;
            end
        end
    end

    % -------------------------------
    % Pack outputs
    % -------------------------------
    out.k          = kvec(:);
    out.time       = time_vec(:);

    out.rmse_h_no   = rmse_h_no; out.rmse_h_as   = rmse_h_as;
    out.rmse_h_no_g   = rmse_h_no_g; out.rmse_h_as_g   = rmse_h_as_g;
    out.rmse_h_no_w   = rmse_h_no_w; out.rmse_h_as_w   = rmse_h_as_w;


    out.rmse_vel_no = rmse_vel_no; out.rmse_vel_as = rmse_vel_as;
    out.rmse_vel_no_g = rmse_vel_no_g; out.rmse_vel_as_g = rmse_vel_as_g;
    out.rmse_vel_no_w = rmse_vel_no_w; out.rmse_vel_as_w = rmse_vel_as_w;

    out.rmse_c_no   = rmse_c_no; out.rmse_c_as   = rmse_c_as;
    out.rmse_c_no_g   = rmse_c_no_g; out.rmse_c_as_g   = rmse_c_as_g;
    out.rmse_c_no_w   = rmse_c_no_w; out.rmse_c_as_w   = rmse_c_as_w;

    out.gl_mid      = gl_mid;
    out.dist_no     = dist_no;
    out.dist_as     = dist_as;

    % -------------------------------
    % figure('Position',[200 200 1100 900]); clf;
    figure('Position',[150 150 1800 1200]); clf;

    fs_axis   = 15;
    fs_label  = 17;
    fs_title  = 16;
    fs_legend = 14;
    
    % Assimilation times: year 2 to 24 every 2 years
    assim_times = 2:2:24;
    
    orange = [0.85 0.325 0.098];
    
    % ---- Create tiledlayout with extra padding for legend space ----
    tl = tiledlayout(4,1, 'TileSpacing','compact', 'Padding','loose');
    % tl = tiledlayout(2,2, 'TileSpacing','compact', 'Padding','loose');
    
    % ---------------- (a) Thickness ----------------
    ax1 = nexttile(tl,1);
    h1 = plot(out.time, out.rmse_h_no,   'r-', 'LineWidth',2.5); hold on
    h2 = plot(out.time, out.rmse_h_as,   'r:', 'LineWidth',2.5); 
    h3 = plot(out.time, out.rmse_h_no_g, 'b-', 'LineWidth',2.5); 
    h4 = plot(out.time, out.rmse_h_as_g, 'b:', 'LineWidth',2.5); 
    h5 = plot(out.time, out.rmse_h_no_w, 'c-', 'LineWidth',2.5); 
    h6 = plot(out.time, out.rmse_h_as_w, 'c:', 'LineWidth',2.5); hold off

    % try semilogy
    % h1 = semilogy(out.time, rmse_h_no,   'r-', 'LineWidth',2.5); hold on
    % h2 = semilogy(out.time, rmse_h_as,   'r:', 'LineWidth',2.5);
    % h3 = semilogy(out.time, rmse_h_no_g, 'b-', 'LineWidth',2.5);
    % h4 = semilogy(out.time, rmse_h_as_g, 'b:', 'LineWidth',2.5);
    % h5 = semilogy(out.time, rmse_h_no_w, 'c-', 'LineWidth',2.5);
    % h6 = semilogy(out.time, rmse_h_as_w, 'c:', 'LineWidth',2.5); hold off
    
    ylabel('RMSE (m)','FontWeight','bold','FontSize',fs_label);
    title('Thickness','FontWeight','bold','FontSize',fs_title);
    % ylim([-0.5,410]); xlim([-1.5,50])
    ylim([-10 440]); xlim([-1.5,50]);
    % yticks([10 30 80 200 420]);
    yticks([30 100 200 300 420]);
    yticklabels({'30', '100', '200', '300', '420'});
    set(gca,'YMinorTick','off');
    
    text(ax1,0.01,0.93,'(a)','Units','normalized', ...
        'FontWeight','bold','FontSize',fs_title, ...
        'HorizontalAlignment','left','VerticalAlignment','top');
    
    % ---------------- (b) Velocity ----------------
    ax2 = nexttile(tl,2);
    plot(out.time, out.rmse_vel_no, 'r-', 'LineWidth',2.5); hold on
    plot(out.time, out.rmse_vel_as, 'r:', 'LineWidth',2.5);
    plot(out.time, out.rmse_vel_no_g, 'b-', 'LineWidth',2.5); 
    plot(out.time, out.rmse_vel_as_g, 'b:', 'LineWidth',2.5); 
    plot(out.time, out.rmse_vel_no_w, 'c-', 'LineWidth',2.5); 
    plot(out.time, out.rmse_vel_as_w, 'c:', 'LineWidth',2.5); hold off

    % semilogy(out.time, rmse_vel_no,   'r-', 'LineWidth',2.5); hold on
    % semilogy(out.time, rmse_vel_as,   'r:', 'LineWidth',2.5);
    % semilogy(out.time, rmse_vel_no_g, 'b-', 'LineWidth',2.5);
    % semilogy(out.time, rmse_vel_as_g, 'b:', 'LineWidth',2.5);
    % semilogy(out.time, rmse_vel_no_w, 'c-', 'LineWidth',2.5);
    % semilogy(out.time, rmse_vel_as_w, 'c:', 'LineWidth',2.5); hold off
    
    ylabel('RMSE (m/yr)','FontWeight','bold','FontSize',fs_label);
    title('Velocity','FontWeight','bold','FontSize',fs_title);
    % ylim([-20,950]); xlim([-1.5,50])
    ylim([-20 960]); xlim([-1.5,50]);
    % yticks([0, 50 200 400  600 800])
    yticks([20 200 400  600  800])
    yticklabels({'20', '200', '400', '600' ,'800'})
    set(gca,'YMinorTick','off');
    
    text(ax2,0.01,0.93,'(b)','Units','normalized', ...
        'FontWeight','bold','FontSize',fs_title, ...
        'HorizontalAlignment','left','VerticalAlignment','top');
    
    % ---------------- (c) Friction ----------------
    ax3 = nexttile(tl,3);
    plot(out.time, out.rmse_c_no, 'r-', 'LineWidth',2.5); hold on
    plot(out.time, out.rmse_c_as, 'r:', 'LineWidth',2.5); 
    plot(out.time, out.rmse_c_no_g, 'b-', 'LineWidth',2.5); 
    plot(out.time, out.rmse_c_as_g, 'b:', 'LineWidth',2.5); hold off
    % semilogy(out.time, rmse_c_no,   'r-', 'LineWidth',2.5); hold on
    % semilogy(out.time, rmse_c_as,   'r:', 'LineWidth',2.5);
    % semilogy(out.time, rmse_c_no_g, 'b-', 'LineWidth',2.5);
    % semilogy(out.time, rmse_c_as_g, 'b:', 'LineWidth',2.5); hold off
    
    
    ylabel('RMSE (Pa m^{-1/3} yr^{-1/3})','FontWeight','bold','FontSize',fs_label);
    title('Friction coefficient','FontWeight','bold','FontSize',fs_title);
    ylim([200,900]); xlim([-1.5,50])
    
    text(ax3,0.01,0.93,'(c)','Units','normalized', ...
        'FontWeight','bold','FontSize',fs_title, ...
        'HorizontalAlignment','left','VerticalAlignment','top');
    
    % ---------------- (d) GL distance ----------------
    ax4 = nexttile(tl,4);
    plot(out.time, out.dist_no./1000, 'r-', 'LineWidth',2.5); hold on
    plot(out.time, out.dist_as./1000, 'r:', 'LineWidth',2.5); hold off
    % semilogy(out.time, out.dist_no./1000, 'r-', 'LineWidth',2.5); hold on
    % semilogy(out.time, out.dist_as./1000, 'r:', 'LineWidth',2.5); hold off
    
    % xlabel('Time (years)','FontWeight','bold','FontSize',fs_label);
    ylabel('|Δx| (km)','FontWeight','bold','FontSize',fs_label);
    title('Absolute distance between GL positions along the centerline', ...
          'FontWeight','bold','FontSize',fs_title);
    ylim([-0.5e4/1000,4e4/1000]); xlim([-1.5,50])
    % yticks([0 1e4/1000 2e4/1000  3e4/1000 4e4/1000])
    % yticklabels({'0','10','20','30','40'})
    set(gca,'YMinorTick','off');
    
    text(ax4,0.01,0.93,'(d)','Units','normalized', ...
        'FontWeight','bold','FontSize',fs_title, ...
        'HorizontalAlignment','left','VerticalAlignment','top');
    
    % ---------------- Global axes style ----------------
    axs = [ax1 ax2 ax3 ax4];
    set(axs, ...
        'FontWeight','bold', ...
        'FontSize',fs_axis, ...
        'Box','on', ...
        'LineWidth',2.0, ...
        'TickDir','out', ...
        'TickLength',[0.004 0.004], ...
        'XGrid','off', ...
        'YGrid','on', ...
        'YMinorGrid','off');
    
    % ---------------- Vertical assimilation lines (thin dotted) -------------
    for ax = axs
        hold(ax,'on')
        for t = assim_times
            xl = xline(ax, t, ':', ...
                'Color',[0.1 0.1 0.1], ...
                'LineWidth',1.8);
            xl.HandleVisibility = 'off';
            xl.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        hold(ax,'off')
    end
    
    % ---------------- Shared legend (reserve space + place it) --------------
    % Create legend from the thickness handles only (6 entries)
    lgd = legend(ax1, [h1 h2 h3 h4 h5 h6], ...
        {'No assimilation (grounded)', ...
         'Assimilated (grounded)', ...
         'No assimilation (grounded excluding GL)', ...
         'Assimilated (grounded excluding GL)', ...
         'No assimilation (whole domain)', ...
         'Assimilated (whole domain)'}, ...
        'FontSize',fs_legend, ...
        'Box','off', ...
        'Orientation','horizontal', ...
        'NumColumns', 3);     % <-- 6 items -> 3 columns => 2 rows
    
    lgd.ItemTokenSize = [16 8];
    lgd.Units = 'normalized';
    
    % --- Reserve a bottom strip for the legend by shrinking tile area ---
    tl.Units = 'normalized';
    tlPos = tl.OuterPosition;            % [x y w h] normalized
    legendStrip = 0.003;                  % <-- bigger strip for 2 rows
    tl.OuterPosition = [tlPos(1) tlPos(2)+legendStrip tlPos(3) tlPos(4)-legendStrip];
    
    % Make sure layout is finalized before positioning legend
    drawnow;
    
    % --- Place legend centered inside the reserved bottom strip ---
    lgdPos = lgd.Position;               % [x y w h]
    lgdPos(1) = 0.5 - lgdPos(3)/2;       % center horizontally
    lgdPos(2) = 0.02;                    % inside bottom margin (raised a bit)
    lgd.Position = lgdPos;
    
    % % Optional safety: if legend is too wide, shrink a little
    % if lgd.Position(1) < 0.01
    %     lgd.Position(1) = 0.01;
    % end
    % if lgd.Position(1) + lgd.Position(3) > 0.99
    %     lgd.Position(1) = 0.99 - lgd.Position(3);
    % end
    xlabel(tl,'Time (years)','FontWeight','bold','FontSize',fs_label);

    % ---- Save figure (300 dpi) ----
    % Use folder relative to THIS script (not MATLAB's current folder)
    scriptdir = fileparts(mfilename('fullpath'));
    outdir    = fullfile(scriptdir, 'figures');
    
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    
    fname = fullfile(outdir, sprintf('RMSE_%s.png', regexprep('hucgl','\.','_')));
    
    exportgraphics(gcf, fname, 'Resolution', 300);
end

function [xc, yc, info] = gl_centerline_point_from_levelset_grid_track( ...
    Xg, Yg, Phig, y_center, x_prev)
% Robust GL "star" on centerline y=y_center:
%   - Find ALL 0-contour segments.
%   - For each segment, find intersections with y=y_center (line crossing).
%   - Pick the intersection closest to x_prev (continuity).
%   - If x_prev is NaN, pick the intersection closest to mid-domain.

    xc = NaN; yc = NaN;

    info.hasGL          = false;
    info.nSegments      = 0;
    info.nIntersections = 0;
    info.reason         = '';

    if ~all(isfinite(Phig(:)))
        Phig(~isfinite(Phig)) = 1;
    end

    % If no sign change anywhere -> no contour
    if min(Phig(:)) * max(Phig(:)) > 0
        info.reason = 'no sign change in grid';
        return;
    end
    info.hasGL = true;

    C = contourc(Xg(1,:), Yg(:,1), Phig, [0 0]);

    % gather all candidate intersections
    xCand = [];
    yCand = [];

    kk = 1; segCount = 0;
    while kk < size(C,2)
        segCount = segCount + 1;
        npts = C(2,kk);
        pts  = C(:, kk+1:kk+npts);
        kk   = kk + npts + 1;

        gx = pts(1,:); gy = pts(2,:);

        % find indices where the polyline crosses y_center
        s = gy - y_center;
        crossIdx = find(s(1:end-1).*s(2:end) <= 0);  % includes hits

        for i = crossIdx(:)'
            y1 = gy(i);   y2 = gy(i+1);
            x1 = gx(i);   x2 = gx(i+1);

            if abs(y2-y1) < eps
                % segment is (nearly) horizontal at y_center; just take midpoint
                xt = 0.5*(x1+x2);
                yt = y_center;
            else
                t  = (y_center - y1) / (y2 - y1);
                xt = x1 + t*(x2-x1);
                yt = y_center;
            end

            if isfinite(xt)
                xCand(end+1,1) = xt; %#ok<AGROW>
                yCand(end+1,1) = yt; %#ok<AGROW>
            end
        end
    end

    info.nSegments = segCount;
    info.nIntersections = numel(xCand);

    if isempty(xCand)
        info.reason = 'no intersection with centerline';
        return;
    end

    % choose target x
    if nargin < 5 || ~isfinite(x_prev)
        x_target = 0.5*(min(Xg(1,:)) + max(Xg(1,:)));  % mid-domain
    else
        x_target = x_prev;
    end

    [~, j] = min(abs(xCand - x_target));
    xc = xCand(j);
    yc = yCand(j);
end

function mask = grounded_mask_from_state(model_state, md, hdim, k, useHpos)
% grounded if phi = H + bed/di > 0 (and optionally H>0)

    if nargin < 5, useHpos = true; end
    di = md.materials.rho_ice / md.materials.rho_water;

    I_h   = 1:hdim;
    I_bed = 4*hdim+1:5*hdim;

    H   = model_state(I_h, k);
    bed = model_state(I_bed, k);
    phi = H + bed/di;

    mask = (phi > 0);
    if useHpos
        mask = mask & (H > 0);
    end
end

function r = rmse_vec(a, b, mask)
% r = sqrt(mean((a-b).^2)) over mask (if provided)

    a = a(:); b = b(:);

    if nargin < 3 || isempty(mask)
        good = isfinite(a) & isfinite(b);
    else
        mask = logical(mask(:));
        good = mask & isfinite(a) & isfinite(b);
    end

    if ~any(good), r = NaN; return; end
    d = a(good) - b(good);
    r = sqrt(mean(d.^2));
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

    % window box
    % plot(ax, ...
    %     [win(1) win(2) win(2) win(1) win(1)], ...
    %     [win(3) win(3) win(4) win(4) win(3)], ...
    %     'w--','LineWidth',1.2,'HandleVisibility','off');

    % GL points (exactly what RMSE uses)
    % if ~isempty(Ptrue)
    %     plot(ax, Ptrue(:,1), Ptrue(:,2), 'k.', 'MarkerSize',10, 'HandleVisibility','on');
    % end
    % if ~isempty(Pno)
    %     plot(ax, Pno(:,1), Pno(:,2), 'r.', 'MarkerSize',10, 'HandleVisibility','on');
    % end
    % if ~isempty(Pas)
    %     plot(ax, Pas(:,1), Pas(:,2), 'c.', 'MarkerSize',10, 'HandleVisibility','on');
    % end

    % center marker
    plot(ax, xc, yc, 'go', 'MarkerFaceColor','g', ...
        'MarkerSize',8,'LineWidth',1.5,'HandleVisibility','on');

    hold(ax,'off')
end

function out = compute_gl_position_rmse_windowed( ...
    k_array, dt, t, ...
    model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md, ...
    gl_mid, varargin)

% Windowed grounding line position RMSE between:
%   - Wrong vs True
%   - Assim vs True
%
% Uses distance between 0-contour polylines of ocean_levelset, restricted to a window.
%
% Required:
%   gl_mid.x, gl_mid.y from your compute_gl_midpoints (TRUE bend-center estimate)
%
% Options:
%   'x_halfwidth' (m): half width of window in x around gl_mid.x (default 30e3)
%   'y_halfwidth' (m): half width of window in y around gl_mid.y (default 20e3)
%   'Nx','Ny': grid for contouring (default 420,70)
%   'minLen': min contour length to keep (default 3e4)
%   'topK': number of longest segments to keep (default 4)

    p = inputParser;
    p.addParameter('x_halfwidth', 30e3, @(v) isnumeric(v) && isscalar(v));
    p.addParameter('y_halfwidth', 20e3, @(v) isnumeric(v) && isscalar(v));
    p.addParameter('Nx', 420, @(v) isnumeric(v) && isscalar(v));
    p.addParameter('Ny', 70,  @(v) isnumeric(v) && isscalar(v));
    p.addParameter('minLen', 3e4, @(v) isnumeric(v) && isscalar(v));
    p.addParameter('topK', 4, @(v) isnumeric(v) && isscalar(v));
    p.parse(varargin{:});
    opt = p.Results;

    nk = numel(k_array);
    rmse_no = nan(nk,1);
    rmse_as = nan(nk,1);

    npts_true = nan(nk,1);
    npts_no   = nan(nk,1);
    npts_as   = nan(nk,1);
    bias_no = nan(nk,1);
    bias_as = nan(nk,1);
    nmatch_true = nan(nk,1);
    nmatch_as   = nan(nk,1);

    for idx = 1:nk
        k = k_array(idx);

        [md_true_k, md_nurged_k, md_ens_k] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);

        % --- window center (use your TRUE bend-center estimate) ---
        xc = gl_mid.x(idx);
        yc = gl_mid.y(idx);
        if ~isfinite(xc) || ~isfinite(yc)
            continue;
        end
        win = [xc-opt.x_halfwidth, xc+opt.x_halfwidth, ...
               yc-opt.y_halfwidth, yc+opt.y_halfwidth];

        % --- extract GL points in window for each state ---
        Ptrue = extract_gl_points_in_window(md_true_k,   opt.Nx, opt.Ny, opt.minLen, opt.topK, win);
        Pno   = extract_gl_points_in_window(md_nurged_k, opt.Nx, opt.Ny, opt.minLen, opt.topK, win);
        Pas   = extract_gl_points_in_window(md_ens_k,    opt.Nx, opt.Ny, opt.minLen, opt.topK, win);

        
        npts_as(idx)   = size(Pas,1);

        if size(Ptrue,1) < 5 || size(Pno,1) < 5 || size(Pas,1) < 5
            continue;
        end

        % --- symmetric RMS distance between curves (point clouds) ---
        % rmse_no(idx) = symmetric_rms_distance(Ptrue, Pno);
        % rmse_as(idx) = symmetric_rms_distance(Ptrue, Pas);
        % rmse_no(idx) = gl_x_offset_rmse(Ptrue, Pno);
        % rmse_as(idx) = gl_x_offset_rmse(Ptrue, Pas);

        s_no = gl_dx_misplacement_stats(Ptrue, Pno);
        s_as = gl_dx_misplacement_stats(Ptrue, Pas);
        
        rmse_no(idx) = s_no.rmse_dx;
        rmse_as(idx) = s_as.rmse_dx;
        
        bias_no(idx) = s_no.bias_dx;
        bias_as(idx) = s_as.bias_dx;
        
        nmatch_true(idx) = s_no.n;
        nmatch_as(idx)   = s_as.n;

        npts_true(idx) = size(Ptrue,1);
        npts_no(idx)   = size(Pno,1);
    end

    out.k = k_array(:);
    out.time = t(k_array(:));
    out.rmse_gl_no = rmse_no;
    out.rmse_gl_as = rmse_as;
    out.npts_true  = npts_true;
    out.npts_no    = npts_no;
    out.npts_as    = npts_as;
    out.bias_dx_no = bias_no;
    out.bias_dx_as = bias_as;
    out.nmatch_no  = nmatch_true;
    out.nmatch_as  = nmatch_as;
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
function rmsd = symmetric_rms_distance(A, B)
% Symmetric RMS nearest-neighbor distance between point sets A and B.
% A,B: (n x 2) arrays of [x y]

    dAB = nn_distances(A, B); % for each point in A: dist to nearest in B
    dBA = nn_distances(B, A); % for each point in B: dist to nearest in A

    d = [dAB; dBA];
    d = d(isfinite(d));

    if isempty(d)
        rmsd = NaN;
    else
        rmsd = sqrt(mean(d.^2));
    end
end

function d = nn_distances(A, B)
% For each point in A, compute distance to nearest point in B.

    nA = size(A,1);
    d  = nan(nA,1);

    % Vectorized block approach to avoid huge memory if needed
    % (but your GL windows are small so this is fine)
    for i = 1:nA
        dx = B(:,1) - A(i,1);
        dy = B(:,2) - A(i,2);
        d(i) = sqrt(min(dx.^2 + dy.^2));
    end
end

function rmse_dx = gl_x_offset_rmse(Ptrue, Pmodel)
% RMSE of x-offsets between model GL and true GL, matched by y
% Ptrue, Pmodel: [x y] arrays within the same window

    if size(Ptrue,1) < 5 || size(Pmodel,1) < 5
        rmse_dx = NaN;
        return
    end

    % sort by y for stable matching
    Ptrue  = sortrows(Ptrue,  2);
    Pmodel = sortrows(Pmodel, 2);

    dx = nan(size(Ptrue,1),1);

    for i = 1:size(Ptrue,1)
        y0 = Ptrue(i,2);
        [~, j] = min(abs(Pmodel(:,2) - y0));  % nearest in y
        dx(i) = Pmodel(j,1) - Ptrue(i,1);     % x-offset ONLY
    end

    dx = dx(isfinite(dx));
    rmse_dx = sqrt(mean(dx.^2));
end

function stats = gl_dx_misplacement_stats(Ptrue, Pmodel)
% Positional misplacement in (approx) normal direction ~ x
% Match model points to true points by nearest y.
%
% Returns:
%   stats.rmse_dx  : sqrt(mean(dx^2))
%   stats.bias_dx  : mean(dx) (signed)
%   stats.med_dx   : median(dx)
%   stats.n        : number of matched points used

    stats.rmse_dx = NaN;
    stats.bias_dx = NaN;
    stats.med_dx  = NaN;
    stats.n       = 0;

    if size(Ptrue,1) < 5 || size(Pmodel,1) < 5
        return
    end

    Ptrue  = sortrows(Ptrue,  2); % sort by y
    Pmodel = sortrows(Pmodel, 2);

    dx = nan(size(Ptrue,1),1);

    for i = 1:size(Ptrue,1)
        y0 = Ptrue(i,2);
        [~, j] = min(abs(Pmodel(:,2) - y0));
        dx(i) = Pmodel(j,1) - Ptrue(i,1);  % +dx means model GL is to the right of true
    end

    dx = dx(isfinite(dx));
    if isempty(dx), return; end

    stats.rmse_dx = sqrt(mean(dx.^2));
    stats.bias_dx = mean(dx);
    stats.med_dx  = median(dx);
    stats.n       = numel(dx);
end
