clear all; close all;

data_path='/Users/bkyanjo3/da_project/ISSM-matlab/examples/ISMIP_Choi/Models/ens_id_0';
file_path_true= fullfile(data_path,"true_state.mat");
file_path_nurged= fullfile(data_path,"nurged_state.mat");
md_true= loadmodel(file_path_true);
md_nurged= loadmodel(file_path_nurged);

k = 80;
gl_true=md_true.results.TransientSolution(k).MaskOceanLevelset;
gl_nurged=md_nurged.results.TransientSolution(k).MaskOceanLevelset;

bed_true=md_true.geometry.bed;
bed_nurged=md_nurged.geometry.bed;

friction_true=md_true.friction.coefficient;
friction_nurged=md_nurged.friction.coefficient;

thickness_true=md_true.results.TransientSolution(k).Thickness;
thickness_nurged=md_nurged.results.TransientSolution(k).Thickness;
% figure(1);
% plotmodel(md_true, 'data', thickness_true, 'title', 'True Ice Thickness'); 
% hold off;

% --- simulate datasets
% load observation data
data_file_paths='_modelrun_datasets';
file_path  = fullfile(data_file_paths, 'synthetic_obs.h5');
w          = h5read(file_path, '/hu_obs')'; 

% load the ensemble data
file_path         = fullfile(data_file_paths, 'icesee_ensemble_data.h5');
ensemble_vec_full = h5read(file_path, '/ensemble'); 
ensemble_vec_mean = h5read(file_path, '/ensemble_mean')';
data_file_paths='_modelrun_datasets';
file_path   = fullfile("data", "ISMIP_initial_data.mat");
md = loadmodel(file_path);
md_ens = md;

% update ens loads
[ndim, nt] = size(ensemble_vec_mean);
num_steps = nt - 1;
hdim = floor(ndim / 3);
md_ens.geometry.bed=ensemble_vec_mean(hdim+1:2*hdim, k);
md_ens.geometry.thickness=ensemble_vec_mean(1:hdim,  k);
md_ens.friction.coefficient=ensemble_vec_mean(2*hdim+1:3*hdim, k);

% figure;
% plotmodel(md_ens, 'data', ensemble_vec_mean(1:hdim,  k), 'title', 'Assimilated Ice Thickness'); hold off;


% thickness
plot_triptych(md_true, md_nurged, md_ens, ...
              'geometry.thickness', sprintf('Ice Thickness after %d years', round((k-1)*0.5)), parula, 'm');   
% bed topography
plot_triptych(md_true, md_nurged, md_ens, ...
              'geometry.bed', sprintf('Bed Elevation after %d years', round((k-1)*0.5)), parula, 'm');
%
% % Frcition coefficient
plot_triptych(md_true, md_nurged, md_ens, ...
              'friction.coefficient', sprintf('Friction Coefficient after %d years', round((k-1)*0.5)), parula, '');
%
% % grounding line
% plot_triptych(md_true, md_nurged, md_ens, ...
%               'results.TransientSolution(499).MaskOceanLevelset', ...
%               'Grounding Line', gray, '');
plot_triptych(md_true, md_nurged, md_ens, ...
              'mask.ocean_levelset', ...
              sprintf('Grounding Line after %d years', round((k-1)*0.5)), parula, '');

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
            pause(0.1); % interactive view
        end
    end

    if make_movie
        close(v);
    end
end


%% -- helper functions -- %%
function plot_triptych(md_true, md_nurged, md_ens, field, field_title, cmap, units)
% Compare true, nudged, assimilated, and difference with two separate colorbars

    if nargin < 6 || isempty(cmap), cmap = parula; end
    if nargin < 7, units = ''; end
    units_str = iff(~isempty(units), [' (' units ')'], '');

    % --- Data ---
    data_true   = get_nested_field(md_true, field);
    data_nurged = get_nested_field(md_nurged, field);
    data_ens    = get_nested_field(md_ens, field);
    diff_data   = data_ens - data_true;

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
    plotmodel(md_true,'data',diff_data, ...
        'title',['Assimilated - True ' field_title], ...
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

