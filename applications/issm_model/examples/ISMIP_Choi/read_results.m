%% -----------------------------------------------------------
% @author:  Brian Kyanjo
% @date:    2025-06-30
% @brief:   Reads and plots results from both ISSM and ICESEE
% ------------------------------------------------------------

close all; clearvars; clear all

global data_file_paths nvar
data_file_paths = '_modelrun_datasets';
nvar = 6;

% ---------------- user toggles ----------------
make_plots       = 0;
make_multi_plots = 0;   % <-- ON (restored)
frames_plot      = 0;

% ---------------- time steps ------------------
% k_array = [1, 20,  50, 80, 120, 170, 220];
k_array= [20,80, 120, 160, 190]+1;
% k_array = [20, 50, 80, 120, 160, 200, 240, 320] +1;
% k_array = [20, 80, 160, 200, 250, 300]+1;
dt      = 0.2;

% ---------------- Load essentials --------------
results_dir = 'results';
filter_type = 'true-wrong';
file_path   = fullfile(results_dir, sprintf('%s-issm.h5', filter_type));
t        = h5read(file_path,'/t'); 
ind_m    = h5read(file_path,'/obs_index'); 
tm_m     = h5read(file_path,'/obs_max_time'); 
run_mode = h5read(file_path,'/run_mode'); 

% true / wrong (nurged)
file_path          = fullfile(data_file_paths, 'true_nurged_states.h5');
model_true_state   = h5read(file_path,'/true_state')';
model_nurged_state = h5read(file_path,'/nurged_state')';

% obs (kept)
file_path = fullfile(data_file_paths, 'synthetic_obs.h5');
w = h5read(file_path, '/hu_obs')'; 

% ensemble mean
file_path         = fullfile(data_file_paths, 'icesee_ensemble_data.h5');
ensemble_vec_mean = h5read(file_path, '/ensemble_mean')';

% ISSM model template
md = loadmodel(fullfile("data","ISMIP.Parameterization1.mat"));
md_true   = md;
md_nurged = md;
md_ens    = md;

% ---------------- GL evolution plot ------------
plot_gl_on_bed_evolution( ...
    k_array, dt, ...
    model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md, ...
    'geometry.thickness', 'Thickness', 'm');

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
        md_true, md_nurged, md_ens, md, 'initialization.vel', 'Velocity', 'm/s');
    plot_var_diff(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'initialization.vel', 'Velocity', 'm/s');

    % bed
    plot_var_evolution(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'geometry.bed', 'Bed Elevation', 'm');
    plot_var_diff(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'geometry.bed', 'Bed', 'm');

    % friction
    plot_var_evolution(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'friction.coefficient', 'Friction Coefficient', '');
    plot_var_diff(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md, 'friction.coefficient', 'Friction', '');
end

% ---------------- optional single triptych -------
if make_plots
    k = k_array(end);
    [md_true_k, md_nurged_k, md_ens_k] = setup_model_states(k, dt, ...
        model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);

    plot_triptych(md_true_k, md_nurged_k, md_ens_k, ...
        'geometry.thickness', sprintf('Ice Thickness after %.1f years', (k-1)*dt), parula, 'm');
    plot_triptych(md_true_k, md_nurged_k, md_ens_k, ...
        'geometry.surface', sprintf('Ice Surface after %.1f years', (k-1)*dt), parula, 'm');
    plot_triptych(md_true_k, md_nurged_k, md_ens_k, ...
        'geometry.bed', sprintf('Bed after %.1f years', (k-1)*dt), parula, 'm');
    plot_triptych(md_true_k, md_nurged_k, md_ens_k, ...
        'mask.ocean_levelset', sprintf('Grounding Line after %.1f years', (k-1)*dt), parula, 'm');
    
end

%% ========================================================================
%%  Helper functions
%% ========================================================================

function [md_true, md_nurged, md_ens] = setup_model_states( ...
    k, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md)
% Initialize the true, wrong (nurged), and ensemble-mean models at index k.
    global nvar
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
    bg_field, bg_title, units)

    if nargin < 11 || isempty(bg_field), bg_field = 'initialization.vel'; end
    if nargin < 12 || isempty(bg_title), bg_title = 'Background'; end
    if nargin < 13, units = ''; end
    units_str = iff(~isempty(units), [' (' units ')'], '');

    nk    = numel(k_array);
    nrows = nk;

    figure('Position',[100 100 1000 (180 + 150*nrows)]); clf;

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
    minLen_wrong = 3e4;
    minLen_ens   = 3e4;
    minArea      = 2;

    keepLargestOnly_true  = true;
    keepLargestOnly_wrong = true;
    keepTopK_ens          = 4;   % allow 2 longest for ensemble (prevents “loss”)

    for idx = 1:nk
        k = k_array(idx);

        [md_true_k, md_nurged_k, md_ens_k] = setup_model_states(k, dt, ...
            model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);

        bg = get_nested_field(md_true_k, bg_field);

        plotmodel(md_true_k, 'data', bg, ...
            'title', sprintf('(%c) %s + GLs (t = %.1f years)', ...
                'a'+(idx-1), bg_title, (k-1)*dt), ...
            'subplot', [nrows, 1, idx], ...
            'caxis', [cmin cmax], ...
            'colorbar', 'off');

        ax = gca;

        % prevent plotmodel objects appearing in legend
        set(ax.Children, 'HandleVisibility','off');

        % grid for contouring
        x = md_true_k.mesh.x(:);
        y = md_true_k.mesh.y(:);
        xg = linspace(min(x), max(x), Nx);
        yg = linspace(min(y), max(y), Ny);
        [Xg, Yg] = meshgrid(xg, yg);

        phi_true  = md_true_k.mask.ocean_levelset(:);
        phi_wrong = md_nurged_k.mask.ocean_levelset(:);
        phi_ens   = md_ens_k.mask.ocean_levelset(:);

        % print diagnostic so you can confirm if ens truly loses sign change
        fprintf('k=%d t=%.1f | ens[min,max]=[%.3e,%.3e]\n', ...
            k, (k-1)*dt, min(phi_ens), max(phi_ens));

        % linear + nearest extrap => NO NaN holes that break contour
        F1 = scatteredInterpolant(x, y, phi_true,  'linear','nearest');
        F2 = scatteredInterpolant(x, y, phi_wrong, 'linear','nearest');
        F3 = scatteredInterpolant(x, y, phi_ens,   'linear','nearest');

        Phi_true  = F1(Xg, Yg);
        Phi_wrong = F2(Xg, Yg);
        Phi_ens   = F3(Xg, Yg);

        hold(ax,'on');
        plot_gl_contour_filtered(Xg, Yg, Phi_true,  'k','-',  2.0, minLen_true,  minArea, keepLargestOnly_true,  1);
        plot_gl_contour_filtered(Xg, Yg, Phi_wrong, 'r','-', 2.0, minLen_wrong, minArea, keepLargestOnly_wrong, 1);
        plot_gl_contour_filtered(Xg, Yg, Phi_ens,   'c','-',  2.0, minLen_ens,   minArea, false,               keepTopK_ens);
        hold(ax,'off');
    end

    % ---- Layout: adaptive spacing ----
    axs = flipud(findall(gcf,'Type','axes'));
    gap = 0.02; top = 0.95; bottom = 0.08;
    avail = top - bottom - (nrows-1)*gap;
    height = avail / nrows;

    if height < 0.05
        fig = gcf;
        scale_factor = max(1, ceil(0.05 / height));
        fig.Position(4) = fig.Position(4) * scale_factor;
        height = 0.05;
    end

    for i = 1:nrows
        pos = [0.10, bottom+(nrows-i)*(height+gap), 0.70, height];
        set(axs(i), 'Position', pos, ...
            'FontWeight','bold','LineWidth',1.2,'Box','on', ...
            'TickDir','out','Layer','top','FontSize',11, ...
            'TickLength',[0.005 0.005]);
        ylabel(axs(i),'y (m)','FontWeight','bold');
        if i < nrows
            set(axs(i),'XTickLabel',[]);
        else
            xlabel(axs(i),'x (m)','FontWeight','bold');
        end
    end

    % ---- Shared colorbar ----
    for i = 1:nrows, colormap(axs(i), parula); end
    cb = colorbar(axs(end), 'Position',[0.83 0.25 0.025 0.45]);
    ylabel(cb, [bg_title units_str], 'FontSize',12,'FontWeight','bold');

    % ---- Clean legend (proxy only) ----
    ax0 = axs(1);
    lg = legend(ax0);
    if ~isempty(lg) && isvalid(lg), delete(lg); end

    hold(ax0,'on');
    p1 = plot(ax0, NaN, NaN, 'k-',  'LineWidth', 2.0);
    p2 = plot(ax0, NaN, NaN, 'r-', 'LineWidth', 2.0);
    p3 = plot(ax0, NaN, NaN, 'c-',  'LineWidth', 2.0);
    lgd = legend(ax0, [p1 p2 p3], ...
        {'True GL','No assimilation GL','Assimilated GL'}, ...
        'Location','northwest','FontSize',10,'Box','on');
    lgd.AutoUpdate = 'on';
    % legend(ax0,'manual');
    hold(ax0,'off');

    set(gcf,'Color','w');
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

    if nargin < 12, units = ''; end
    units_str = iff(~isempty(units), [' (' units ')'], '');

    nk    = length(k_array);
    nrows = 2 + nk;

    figure('Position',[100 100 1000 150 + 150*nrows]); clf;

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
    plotmodel(md_true_k,'data',data_true,'title',sprintf('(a) True %s',field_title), ...
        'subplot',[nrows,1,1],'caxis',[cmin cmax],'colorbar','off');

    % (b) No assimilation - True
    [md_true_1, ~, md_ens_1] = setup_model_states(1, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);
    % diff_no = get_nested_field(md_ens_1, field) - get_nested_field(md_true_1, field);
    ens_field = get_nested_field(md_ens_1, field);
    true_field = get_nested_field(md_true_1, field);
    if contains(field,'geometry.bed')
        % diff_no = relative_error(ens_field, true_field);
        diff_no = signed_log_relerr(ens_field, true_field);
        % diff_no = relerr_percent_clipped(ens_field, true_field);
    else
        diff_no = ens_field - true_field;
    end
    maxAbs_no = max(abs(diff_no(:)));
  
    % maxAbs_no = prctile(abs(diff_no(:)), 99);
    plotmodel(md_ens_1,'data',diff_no,'title','(b) No assimilation − True', ...
        'subplot',[nrows,1,2],'caxis',[-maxAbs_no maxAbs_no],'colorbar','off');

    % (c..): Assim - True
    for idx = 1:nk
        k = k_array(idx);
        [md_true_k, ~, md_ens_k] = setup_model_states(k, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
        ens_field = get_nested_field(md_ens_k, field);
        true_field = get_nested_field(md_true_k, field);
        if contains(field,'geometry.bed')
            diff_k = relative_error(ens_field, true_field);
            % diff_k = signed_log_relerr(ens_field, true_field);
            % diff_k = relerr_percent_clipped(ens_field, true_field);
        else
            diff_k = ens_field - true_field;
        end
        % diff_k = get_nested_field(md_ens_k, field) - get_nested_field(md_true_k, field);

        maxAbs = max(abs(diff_k(:)));
      
        % maxAbs = prctile(abs(diff_k(:)), 99);
        label  = sprintf('(%c)', 'b'+idx);
        plotmodel(md_ens_k,'data',diff_k, ...
            'title',sprintf('%s Assimilated − True (after %.1f years)', label, (k-1)*dt), ...
            'subplot',[nrows,1,idx+2],'caxis',[-maxAbs maxAbs],'colorbar','off');
    end

    % layout (your adaptive spacing)
    axs = flipud(findall(gcf,'Type','axes'));
    gap = 0.02; top = 0.95; bottom = 0.08;
    avail = top-bottom - (nrows-1)*gap;
    height = avail/nrows;

    if height < 0.05
        fig = gcf;
        fig.Position(4) = fig.Position(4) * max(1, ceil(0.05/height));
        height = 0.05;
    end

    for i = 1:nrows
        pos = [0.10, bottom+(nrows-i)*(height+gap), 0.70, height];
        set(axs(i),'Position',pos,'FontWeight','bold','LineWidth',1.2,'Box','on', ...
            'TickDir','out','Layer','top','FontSize',11,'TickLength',[0.005 0.005]);
        ylabel(axs(i),'y (km)','FontWeight','bold');
        if i < nrows, set(axs(i),'XTickLabel',[]);
        else, xlabel(axs(i),'x (km)','FontWeight','bold'); end
    end

    cb1 = colorbar(axs(1), 'Position',[0.83 0.68 0.025 0.16]);
    ylabel(cb1,[field_title units_str],'FontSize',12,'FontWeight','bold');
    colormap(axs(1), parula);

    for i = 2:nrows, colormap(axs(i), redblue(256)); end
    cb2 = colorbar(axs(end), 'Position',[0.83 0.25 0.025 0.40]);
    % ylabel(cb2,'Difference','FontSize',12,'FontWeight','bold');
    if contains(field,'geometry.bed')
        ylabel(cb2,'Relative Error','FontSize',12,'FontWeight','bold');
    else
        ylabel(cb2,['Difference' units_str],'FontSize',12,'FontWeight','bold');
    end

    set(gcf,'Color','w');
end

function plot_var_evolution(k_array, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
    md_true, md_nurged, md_ens, md, field, field_title, units)

    if nargin < 12, field_title = field; end
    if nargin < 13, units = ''; end
    units_str = iff(~isempty(units), [' (' units ')'], '');

    nk    = length(k_array);
    nrows = 2 + nk;

    figure('Position',[100 100 1000 150 + 150*nrows]); clf;

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
    [md_true_last, ~, ~] = setup_model_states(k_array(end), dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);
    data_true = get_nested_field(md_true_last, field);
    plotmodel(md_true_last,'data',data_true, ...
        'title',sprintf('(a) True %s (after %.1f years)', field_title, (k_array(end)-1)*dt), ...
        'subplot',[nrows,1,1],'caxis',[cmin cmax],'colorbar','off');

    % (b) No assimilation (k=1) ensemble
    [~, ~, md_ens_1] = setup_model_states(1, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
        md_true, md_nurged, md_ens, md);
    data_ens = get_nested_field(md_ens_1, field);
    plotmodel(md_ens_1,'data',data_ens, ...
        'title',sprintf('(b) No assimilation %s', field_title), ...
        'subplot',[nrows,1,2],'caxis',[cmin cmax],'colorbar','off');

    % (c..) Assim snapshots
    for idx = 1:nk
        k = k_array(idx);
        [~, ~, md_ens_k] = setup_model_states(k, dt, model_true_state, model_nurged_state, ensemble_vec_mean, ...
            md_true, md_nurged, md_ens, md);
        data_ens = get_nested_field(md_ens_k, field);
        label = sprintf('(%c)', 'b'+idx);
        plotmodel(md_ens_k,'data',data_ens, ...
            'title',sprintf('%s Assimilated %s (after %.1f years)', label, field_title, (k-1)*dt), ...
            'subplot',[nrows,1,idx+2],'caxis',[cmin cmax],'colorbar','off');
    end

    % layout (your adaptive spacing)
    axs = flipud(findall(gcf,'Type','axes'));
    gap = 0.02; top = 0.95; bottom = 0.08;
    avail = top-bottom - (nrows-1)*gap;
    height = avail/nrows;

    if height < 0.05
        fig = gcf;
        fig.Position(4) = fig.Position(4) * max(1, ceil(0.05/height));
        height = 0.05;
    end

    for i = 1:nrows
        pos = [0.10, bottom+(nrows-i)*(height+gap), 0.70, height];
        set(axs(i),'Position',pos,'FontWeight','bold','LineWidth',1.2,'Box','on', ...
            'TickDir','out','Layer','top','FontSize',11,'TickLength',[0.005 0.005]);
        ylabel(axs(i),'y (km)','FontWeight','bold');
        if i < nrows, set(axs(i),'XTickLabel',[]);
        else, xlabel(axs(i),'x (km)','FontWeight','bold'); end
    end

    for i = 1:nrows, colormap(axs(i), parula); end
    cb = colorbar(axs(end), 'Position',[0.83 0.25 0.025 0.45]);
    ylabel(cb,[field_title units_str],'FontSize',12,'FontWeight','bold');

    set(gcf,'Color','w');
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

    % --- First colorbar (top 3) ---
    for i = 1:3
        colormap(axs(i), cmap);
        caxis(axs(i), [cmin cmax]);
    end
    cb1 = colorbar(axs(2),'Position',[0.83 0.415 0.025 0.35]);
    static_field = regexprep(field_title,'\s+after.*','');
    ylabel(cb1,[static_field units_str],'FontSize',13,'FontWeight','bold');
    cb1.FontSize = 11;
    set(cb1,'Box','on','LineWidth',1.2);

    % --- Second colorbar (relative error) ---
    ax_diff = axs(4);
    colormap(ax_diff, redblue(256));
    caxis(ax_diff,[-maxAbs maxAbs]);
    pos_diff = get(ax_diff,'Position');
    cb2 = colorbar(ax_diff,'Position',[0.83 pos_diff(2)+0.14 0.025 pos_diff(4)-0.28]);
    ylabel(cb2,'Relative error','FontSize',13,'FontWeight','bold');
    cb2.FontSize = 11;
    set(cb2,'Box','on','LineWidth',1.2);

    set(gcf,'Color','w');
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

