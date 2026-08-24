%% -----------------------------------------------------------
% @author:  Brian Kyanjo
% @date:    2025-04-30  (revised: contour overlays)
% @brief:   Reads and plots ISSM/ICESEE results with triptych and contours
% ------------------------------------------------------------

close all; clearvars;

data_file_paths = 'data3/_modelrun_datasets';
results_dir     = data_file_paths;
filter_type     = 'true-wrong';

% --- Load metadata and states ---
file_path   = fullfile(results_dir, sprintf('%s-issm.h5', filter_type));
t           = h5read(file_path,'/t');
run_mode    = h5read(file_path,'/run_mode');

file_path         = fullfile(data_file_paths, 'synthetic_obs.h5');
[ind_m, tm_m]     = read_observation_metadata(file_path);
w                 = h5read(file_path, '/hu_obs')'; %#ok<NASGU>

file_path            = fullfile(data_file_paths, 'true_nurged_states.h5');
model_true_state     = h5read(file_path,'/true_state')';     % (nd, nt)
model_nurged_state   = h5read(file_path,'/nurged_state')';   % (nd, nt)

file_path         = fullfile(data_file_paths, 'icesee_ensemble_data.h5');
ensemble_vec_full = h5read(file_path, '/ensemble');          %#ok<NASGU>
ensemble_vec_mean = h5read(file_path, '/ensemble_mean')';    % (nd, nt)

% --- Geometry / base model ---
md = loadmodel(fullfile("data","ISMIP_initial_data.mat"));
md_true   = md; md_nurged = md; md_ens = md;

[ndim, nt] = size(model_true_state);
hdim       = floor(ndim/3);
k          = nt-1;                             % frame to show
years      = (k-1)*0.5;                        % your time conversion

% True state
md_true.geometry.bed        = model_true_state(hdim+1:2*hdim,  k);
md_true.geometry.thickness  = model_true_state(1:hdim,        k);
md_true.friction.coefficient= model_true_state(2*hdim+1:3*hdim,k);

% Nudged state
md_nurged.geometry.bed        = model_nurged_state(hdim+1:2*hdim,  k);
md_nurged.geometry.thickness  = model_nurged_state(1:hdim,        k);
md_nurged.friction.coefficient= model_nurged_state(2*hdim+1:3*hdim,k);

% Ensemble mean
md_ens.geometry.bed         = ensemble_vec_mean(hdim+1:2*hdim, k);
md_ens.geometry.thickness   = ensemble_vec_mean(1:hdim,        k);
md_ens.friction.coefficient = ensemble_vec_mean(2*hdim+1:3*hdim,k);

% Grounding-line mask (flotation approx here; use TransientSolution if available)
di = md.materials.rho_ice / md.materials.rho_water;
md_true. ( 'mask').( 'ocean_levelset')   = md_true.geometry.thickness   + md_true.geometry.bed   / di;
md_nurged.('mask').('ocean_levelset')    = md_nurged.geometry.thickness + md_nurged.geometry.bed / di;
md_ens.   ('mask').('ocean_levelset')    = md_ens.geometry.thickness    + md_ens.geometry.bed    / di;

%% ---------------- Triptych with contour overlays ----------------
% Example 1: Bed + GL=0 contour (orange)
ov1.field   = 'mask.ocean_levelset';
ov1.levels  = 0;                      % can be scalar or vector [0 100 200]
ov1.color   = [0.85 0.33 0.10];
ov1.width   = 1.5;

plot_triptych(md_true, md_nurged, md_ens, ...
    'geometry.bed', sprintf('Bed Elevation after %d years', years), ...
    parula, 'm', ov1);

% Example 2: Thickness + optional thickness contours every 100 m
% ov2.field  = 'geometry.thickness';
% ov2.levels = 0:100:1000; ov2.color=[1 1 1]; ov2.width=0.8;
% plot_triptych(md_true, md_nurged, md_ens, ...
%     'geometry.thickness', sprintf('Ice Thickness after %d years', years), ...
%     parula, 'm', ov2);

% Example 3: Friction (no overlay)
plot_triptych(md_true, md_nurged, md_ens, ...
    'friction.coefficient', sprintf('Friction Coefficient after %d years', years), ...
    parula, '');

% Example 4: GL mask background + GL=0 contour (redundant but okay)
plot_triptych(md_true, md_nurged, md_ens, ...
    'mask.ocean_levelset', sprintf('Grounding Line after %d years', years), ...
    parula, '', ov1);

%% ---------------- Movie of GL every 10 years ----------------
make_movie = true;
if make_movie
    figure_dir = fullfile(data_file_paths, 'figures');
    if ~exist(figure_dir,'dir'), mkdir(figure_dir); end
    v = VideoWriter(fullfile(figure_dir,'groundingline_triptych.mp4'),'MPEG-4');
    v.FrameRate = 10; open(v);
    hfig = figure('Position',[100 100 1000 700]);

    for k = 1:20:501 % ends at 250 years
        % update ensemble mean
        md_ens.geometry.bed         = ensemble_vec_mean(hdim+1:2*hdim, k);
        md_ens.geometry.thickness   = ensemble_vec_mean(1:hdim,        k);
        md_ens.friction.coefficient = ensemble_vec_mean(2*hdim+1:3*hdim,k);
        md_ens.mask.ocean_levelset  = md_ens.geometry.thickness + md_ens.geometry.bed/di;

        % update true/nudged GL (use TransientSolution if available)
        md_true.mask.ocean_levelset   = md_true.geometry.thickness   + md_true.geometry.bed/di;
        md_nurged.mask.ocean_levelset = md_nurged.geometry.thickness + md_nurged.geometry.bed/di;

        clf(hfig);
        plot_triptych(md_true, md_nurged, md_ens, ...
            'mask.ocean_levelset', ...
            sprintf('Grounding Line after %d years', (k-1)*0.5), ...
            parula, '', ov1);

        drawnow; writeVideo(v, getframe(hfig));
    end
    close(v);
end

%% ====================== Helper functions ======================

function plot_triptych(md_true, md_nurged, md_ens, field, field_title, cmap, units, overlay)
% PLOT_TRIPTYCH plots 3 panels (true/nudged/assimilated) with optional overlays.
    if nargin < 6 || isempty(cmap), cmap = parula; end
    if nargin < 7, units = ''; end
    if nargin < 8, overlay = []; end

    units_str   = iff(~isempty(units), [' (' units ')'], '');
    data_true   = get_nested_field(md_true,   field);
    data_nurged = get_nested_field(md_nurged, field);
    data_ens    = get_nested_field(md_ens,    field);

    cmin = min([min(data_true(:)), min(data_nurged(:)), min(data_ens(:))]);
    cmax = max([max(data_true(:)), max(data_nurged(:)), max(data_ens(:))]);

    figure('Position',[100 100 900 620]); clf;

    plotmodel(md_true,  'data',data_true,  'title',['True '        field_title], 'subplot',[3,1,1], 'caxis',[cmin cmax], 'colorbar','off');
    plotmodel(md_nurged,'data',data_nurged,'title',['Nudged '      field_title], 'subplot',[3,1,2], 'caxis',[cmin cmax], 'colorbar','off');
    plotmodel(md_ens,   'data',data_ens,   'title',['Assimilated ' field_title], 'subplot',[3,1,3], 'caxis',[cmin cmax], 'colorbar','off');

    axs = flipud(findall(gcf,'Type','axes'));
    n = numel(axs); gap = -0.125; top = 0.95; bottom = 0.08;
    height = (top-bottom - (n-1)*gap)/n;
    cb_height = 0.75*(top-bottom); cb_bottom = bottom + (top-bottom-cb_height)/2;

    cb = colorbar('Position',[0.88 cb_bottom 0.03 cb_height]); colormap(cmap);
    static_field = regexprep(field_title, '\s+after.*', '');
    ylabel(cb,[static_field units_str],'FontSize',14,'FontWeight','bold'); cb.FontSize = 12;

    % Layout + km ticks
    for i = 1:n
        pos = [0.10, bottom+(n-i)*(height+gap), 0.75, height];
        set(axs(i),'Position',pos);
        xt = get(axs(i),'XTick'); yt = get(axs(i),'YTick');
        set(axs(i),'XTickLabel',xt/1000,'FontSize',12);
        set(axs(i),'YTickLabel',yt/1000,'FontSize',12);
        xlabel(axs(i),'X (km)','FontSize',12);
        ylabel(axs(i),'Y (km)','FontSize',12);
    end

    % ---------- Robust contour overlay (scatteredInterpolant + contour) ----------
    if ~isempty(overlay)
        if ~isfield(overlay,'field'),  error('overlay.field is required'); end
        lvls  = getfielddef(overlay,'levels',0);
        col   = getfielddef(overlay,'color',[0.85 0.33 0.10]);
        lw    = getfielddef(overlay,'width',1.5);

        mds = {md_true, md_nurged, md_ens};
        for pi = 1:3
            vals = get_nested_field(mds{pi}, overlay.field);      % node values
            draw_contour_scattered(mds{pi}, vals, lvls, axs(pi), col, lw);
        end
    end
    % ---------------------------------------------------------------------------

    sgtitle(['Comparison of ' field_title ' States'],'FontSize',16);
end

function draw_contour_scattered(md, node_vals, levels, ax, color, lw, do_labels)
% Interpolate scattered node data to a grid and draw contours at given levels.
% levels: scalar or vector of contour values (e.g., 0 or [0 100 200])
% do_labels (optional): true/false to place labels. Default: false.

    if nargin < 7, do_labels = false; end
    assert(numel(node_vals)==numel(md.mesh.x), ...
        'node_vals must be length(md.mesh.x)');

    % Grid resolution (tune as needed)
    nx = 600; ny = 140;
    xi = linspace(min(md.mesh.x), max(md.mesh.x), nx);
    yi = linspace(min(md.mesh.y), max(md.mesh.y), ny);
    [XI,YI] = meshgrid(xi, yi);

    % Interpolate scattered node data onto grid
    F = scatteredInterpolant(md.mesh.x, md.mesh.y, node_vals, 'linear', 'nearest');
    ZI = F(XI, YI);

    axes(ax); hold(ax, 'on');

    % Ensure levels is a row vector
    levels = levels(:).';

    % Draw contours (x,y in km)
    [C,h] = contour(XI/1000, YI/1000, ZI, levels, ...
                    'LineColor', color, 'LineWidth', lw);

    % Only try to label if there are contours and labeling was requested
    if do_labels && ~isempty(C) && ~isempty(h) && size(C,1) == 2
        % Use automatic placement; remove 'manual' to avoid interactive blocking
        clabel(C, h, 'Color', color, 'FontSize', 9, 'LabelSpacing', 400);
    end

    hold(ax, 'off');
end

function out = get_nested_field(s, field)
% Access nested fields: 'geometry.thickness' or 'results.TransientSolution(10).MaskOceanLevelset'
    parts = strsplit(field,'.'); out = s;
    for i = 1:numel(parts)
        tok = parts{i};
        t = regexp(tok,'(.+)\((\d+)\)$','tokens');
        if ~isempty(t)
            base = t{1}{1}; k = str2double(t{1}{2});
            out = out.(base)(k);
        else
            out = out.(tok);
        end
    end
end

function y = getfielddef(s, name, default)
% get field or default
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), y = s.(name); else, y = default; end
end

function y = iff(cond, a, b)
% inline ternary
    if cond, y = a; else, y = b; end
end
