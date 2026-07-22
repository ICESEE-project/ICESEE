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

run_def = struct( ...
    'key',    {'WBF','EBF','IBF'}, ...
    'folder', {'_reviewer_wrong_friction_fixed_v5', ...
               '_reviewer_friction_enkf_only_v5', ...
               '_reviewer_friction_inversion_hybrid_v5'}, ...
    'title',  {'WBF: EnKF state/bed update with fixed wrong friction', ...
               'EBF: EnKF-only friction recovery', ...
               'IBF: EnKF state/bed update plus friction inversion'});

dt_fallback = 0.2;                % years, used only if no /t dataset is found
assimilation_times = 2:2:24;      % years; set to [] to suppress guide lines
rho_ice_fallback   = 917;         % kg m^-3
rho_water_fallback = 1028;        % kg m^-3

figure_dir = fullfile(script_dir, 'figures');
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

    [runs(r).rmse_velocity, runs(r).rmse_surface, ...
     runs(r).rmse_gl_km, runs(r).gl_green_x, runs(r).gl_green_y] = ...
        diagnostics(runs(r).truth, runs(r).estimate, ...
        x, y, density_ratio, n_nodes);
end

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

%% ---------------- Draw the requested stacked figure --------------------
colors = [0.0000 0.4470 0.7410; ...
          0.8500 0.3250 0.0980; ...
          0.4660 0.6740 0.1880];
styles = {'-','--','-.'};

true_friction = runs(1).truth(5*n_nodes+1:6*n_nodes,runs(1).final_k);
map_data = {true_friction, runs(1).friction_final, ...
            runs(2).friction_final, runs(3).friction_final};

% Basal friction is meaningful only beneath grounded ice. Apply the same
% final true-grounded mask to every map for a like-for-like comparison.
true_H_final = runs(1).truth(1:n_nodes,runs(1).final_k);
true_bed_final = runs(1).truth(4*n_nodes+1:5*n_nodes,runs(1).final_k);
true_grounded_final = true_H_final > 0 & ...
    true_H_final + true_bed_final/density_ratio > 0;
for p = 1:numel(map_data)
    map_data{p}(~true_grounded_final) = NaN;
end
map_titles = { ...
    sprintf('(a) True basal friction at t = %g years',final_time), ...
    '(b) WBF: fixed wrong basal friction', ...
    sprintf('(c) EBF: EnKF-only recovery at t = %g years',runs(2).final_time), ...
    sprintf('(d) IBF: inversion recovery at t = %g years',runs(3).final_time)};

finite_friction = map_data{1}(isfinite(map_data{1}));
if isempty(finite_friction)
    error('All final true basal-friction values are non-finite.');
end
% Use the true grounded-friction range for every panel. Recovered outliers
% saturate instead of stretching the scale and hiding the true structure.
friction_limits = [min(finite_friction), max(finite_friction)];
if friction_limits(1) == friction_limits(2)
    friction_limits = friction_limits + [-1 1];
end

% Manual panel positions follow the shared script: broad, shallow map panels,
% bold boxed axes, one shared map colorbar, and three aligned RMSE panels.
fig = figure('Color','w','Units','pixels','Position',[80 40 1800 1500]);

map_left = 0.09; map_width = 0.78; map_height = 0.075;
map_bottom = [0.890 0.790 0.690 0.590];
map_axes = gobjects(4,1);
x_limits = [min(x) max(x)]/1000;
y_limits = [min(y) max(y)]/1000;

for p = 1:4
    ax = axes(fig,'Position',[map_left map_bottom(p) map_width map_height]);
    map_axes(p) = ax;
    plot_basal_friction(ax,elements,x,y,map_data{p}, ...
        x_limits,y_limits,friction_limits);
    title(ax,map_titles{p},'FontWeight','bold','FontSize',14);
    ylabel(ax,'y (km)','FontWeight','bold','FontSize',14);
    if p < 4
        ax.XTickLabel = [];
    else
        xlabel(ax,'x (km)','FontWeight','bold','FontSize',14);
    end
end

% One colorbar for all four maps, as in the shared evolution figures.
middle_map_position = [map_left map_bottom(2) map_width map_height];
cb = colorbar(map_axes(2));
cb.Position = [0.89 0.590 0.017 0.375];
map_axes(2).Position = middle_map_position;
cb.FontSize = 13;
cb.FontWeight = 'bold';
cb.LineWidth = 1.2;
ylabel(cb,'Basal friction (Pa m^{-1/3} yr^{-1/3})', ...
    'FontSize',15,'FontWeight','bold');

ax_vel = axes(fig,'Position',[0.09 0.395 0.84 0.115]);
plot_metric(ax_vel,runs,'rmse_velocity',colors,styles, ...
    assimilation_times,'(e) Grounded-ice speed','RMSE (m yr^{-1})',false);

ax_surf = axes(fig,'Position',[0.09 0.225 0.84 0.115]);
plot_metric(ax_surf,runs,'rmse_surface',colors,styles, ...
    assimilation_times,'(f) Grounded-ice surface elevation','RMSE (m)',false);

ax_gl = axes(fig,'Position',[0.09 0.065 0.84 0.115]);
plot_metric(ax_gl,runs,'rmse_gl_km',colors,styles, ...
    assimilation_times,'(g) Grounding-line RMSE near green midpoint', ...
    'RMSE \Delta x (km)',true);

legend(ax_vel,{runs.key},'Location','northwest','Orientation','horizontal', ...
    'Box','off','FontWeight','bold','FontSize',13);

set([map_axes; ax_vel; ax_surf; ax_gl], ...
    'FontName','Helvetica','FontSize',15,'FontWeight','bold', ...
    'LineWidth',2.0,'Box','on','TickDir','out','Layer','top', ...
    'TickLength',[0.004 0.004]);

png_file = fullfile(figure_dir,'wbf_ebf_ibf_comparison.png');
pdf_file = fullfile(figure_dir,'wbf_ebf_ibf_comparison.pdf');
exportgraphics(fig,png_file,'Resolution',300);
exportgraphics(fig,pdf_file,'ContentType','vector');

metrics_file = fullfile(figure_dir,'wbf_ebf_ibf_metrics.mat');
metrics_runs = rmfield(runs,{'truth','estimate','background','valid_columns'});
save(metrics_file,'metrics_runs','run_def','assimilation_times', ...
    'dt_fallback','final_time');

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
    candidates = { ...
        fullfile(data_dir,'true-wrong-issm.h5'), ...
        fullfile(data_dir,'results','true-wrong-issm.h5'), ...
        fullfile(root_dir,'results','true-wrong-issm.h5')};
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

function [rmse_velocity,rmse_surface,rmse_gl_km,gl_green_x,gl_green_y] = diagnostics( ...
        truth,estimate,x,y,density_ratio,n)
    nt = size(truth,2);
    I_H  = 1:n;
    I_S  = n+1:2*n;
    I_Vx = 2*n+1:3*n;
    I_Vy = 3*n+1:4*n;
    I_b  = 4*n+1:5*n;

    rmse_velocity = nan(nt,1);
    rmse_surface = nan(nt,1);
    rmse_gl_km = nan(nt,1);
    gl_green_x = nan(nt,1);
    gl_green_y = nan(nt,1);

    y_center = 0.5*(min(y)+max(y));
    xg = linspace(min(x),max(x),420);
    yg = linspace(min(y),max(y),70);
    [Xg,Yg] = meshgrid(xg,yg);

    for k = 1:nt
        Ht = truth(I_H,k);
        St = truth(I_S,k);
        Ut = hypot(truth(I_Vx,k),truth(I_Vy,k));
        Ue = hypot(estimate(I_Vx,k),estimate(I_Vy,k));
        Se = estimate(I_S,k);

        phi_t = Ht + truth(I_b,k)/density_ratio;
        phi_e = estimate(I_H,k) + estimate(I_b,k)/density_ratio;

        % The true grounded-ice domain is the common evaluation mask.
        mask = Ht > 0 & phi_t > 0 & isfinite(Ut) & isfinite(Ue);
        rmse_velocity(k) = vector_rmse(Ue,Ut,mask);
        mask_s = Ht > 0 & phi_t > 0 & isfinite(St) & isfinite(Se);
        rmse_surface(k) = vector_rmse(Se,St,mask_s);

        % Match the green-point diagnostic in the shared script: find the
        % true GL midpoint on the longest zero contour, then calculate the
        % x-position RMSE inside a 60 km x 40 km window around that point.
        Phi_t = interpolate_levelset(x,y,phi_t,Xg,Yg);
        Phi_e = interpolate_levelset(x,y,phi_e,Xg,Yg);
        [xc,yc] = green_gl_midpoint(Xg,Yg,Phi_t,y_center);
        gl_green_x(k) = xc;
        gl_green_y(k) = yc;
        if isfinite(xc) && isfinite(yc)
            win = [xc-30e3 xc+30e3 yc-20e3 yc+20e3];
            Ptrue = gl_points_in_window(Xg,Yg,Phi_t,win,30e3,4);
            Pest = gl_points_in_window(Xg,Yg,Phi_e,win,30e3,4);
            rmse_gl_km(k) = gl_x_rmse(Ptrue,Pest)/1000;
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

function Phi = interpolate_levelset(x,y,phi,Xg,Yg)
    good = isfinite(x) & isfinite(y) & isfinite(phi);
    if nnz(good) < 3
        Phi = nan(size(Xg));
        return
    end
    F = scatteredInterpolant(x(good),y(good),phi(good),'linear','nearest');
    Phi = F(Xg,Yg);
end

function [xc,yc] = green_gl_midpoint(Xg,Yg,Phi,y_center)
% Green point: centerline crossing of the longest true GL contour.
    xc = NaN; yc = NaN;
    if ~any(isfinite(Phi(:))), return, end
    Phi(~isfinite(Phi)) = 1;
    if min(Phi(:))*max(Phi(:)) > 0, return, end
    C = contourc(Xg(1,:),Yg(:,1),Phi,[0 0]);
    [segments,lengths] = unpack_contours(C,0);
    if isempty(segments), return, end
    [~,idx] = max(lengths);
    points = segments{idx};
    gx = points(1,:); gy = points(2,:);

    crossing = find((gy(1:end-1)-y_center).*(gy(2:end)-y_center) <= 0,1);
    if isempty(crossing)
        [~,j] = min(abs(gy-y_center));
        xc = gx(j); yc = gy(j);
    else
        i = crossing;
        if abs(gy(i+1)-gy(i)) < eps
            alpha = 0.5;
        else
            alpha = (y_center-gy(i))/(gy(i+1)-gy(i));
        end
        xc = gx(i)+alpha*(gx(i+1)-gx(i));
        yc = y_center;
    end
end

function P = gl_points_in_window(Xg,Yg,Phi,win,min_length,top_k)
    P = zeros(0,2);
    if ~any(isfinite(Phi(:))), return, end
    Phi(~isfinite(Phi)) = 1;
    if min(Phi(:))*max(Phi(:)) > 0, return, end
    C = contourc(Xg(1,:),Yg(:,1),Phi,[0 0]);
    [segments,lengths] = unpack_contours(C,min_length);
    if isempty(segments), return, end
    [~,order] = sort(lengths,'descend');
    order = order(1:min(top_k,numel(order)));
    for j = order(:).'
        points = segments{j};
        inside = points(1,:) >= win(1) & points(1,:) <= win(2) & ...
                 points(2,:) >= win(3) & points(2,:) <= win(4);
        P = [P; points(:,inside).']; %#ok<AGROW>
    end
    if size(P,1) > 1, P = unique(round(P,3),'rows'); end
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

function value = gl_x_rmse(Ptrue,Pestimate)
    value = NaN;
    if size(Ptrue,1) < 5 || size(Pestimate,1) < 5, return, end
    Ptrue = sortrows(Ptrue,2);
    Pestimate = sortrows(Pestimate,2);
    dx = nan(size(Ptrue,1),1);
    for i = 1:size(Ptrue,1)
        [~,j] = min(abs(Pestimate(:,2)-Ptrue(i,2)));
        dx(i) = Pestimate(j,1)-Ptrue(i,1);
    end
    dx = dx(isfinite(dx));
    if ~isempty(dx), value = sqrt(mean(dx.^2)); end
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

function plot_metric(ax,runs,field_name,colors,styles,assim_times, ...
        title_text,ylabel_text,show_xlabel)
    hold(ax,'on');
    for r = 1:numel(runs)
        plot(ax,runs(r).time,runs(r).(field_name), ...
            'Color',colors(r,:),'LineStyle',styles{r},'LineWidth',2.2);
    end
    for ta = assim_times
        if ta == max(assim_times)
            line_handle = xline(ax,ta,'--','Color',[0.20 0.20 0.20], ...
                'LineWidth',1.4);
        else
            line_handle = xline(ax,ta,':','Color',[0.65 0.65 0.65], ...
                'LineWidth',0.8);
        end
        line_handle.HandleVisibility = 'off';
    end
    hold(ax,'off');
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    ax.GridAlpha = 0.18;
    title(ax,title_text,'FontWeight','bold','FontSize',15);
    ylabel(ax,ylabel_text,'FontWeight','bold','FontSize',16);
    if show_xlabel
        xlabel(ax,'Time (years)','FontWeight','bold','FontSize',16);
    else
        ax.XTickLabel = [];
    end
end
