% if any(steps == 8)
    % Load models
    initial_state = fullfile("Data", "ISMIP.reference_simulation_0.mat");
    % true_state = fullfile("data", "true_data.mat");
    % nurged_state = fullfile("data", "nurged_data.mat");
    true_state = fullfile("issm_data","true_state.mat");
    nurged_state = fullfile("issm_data","nurged_state.mat");
    md_initial = loadmodel(initial_state);
    md_true = loadmodel(true_state);
    md_wrong = loadmodel(nurged_state);

    % Update initial state fields
    % md_initial.geometry.thickness = md_initial.results.TransientSolution(end).Thickness;
    % md_initial.geometry.surface   = md_initial.results.TransientSolution(end).Surface;
    % md_initial.geometry.base      = md_initial.results.TransientSolution(end).Base;
    % md_initial.initialization.vx  = md_initial.results.TransientSolution(end).Vx;
    % md_initial.initialization.vy  = md_initial.results.TransientSolution(end).Vy;
    % md_initial.initialization.vel = md_initial.results.TransientSolution(end).Vel;
    % md_initial.initialization.pressure = md_initial.results.TransientSolution(end).Pressure;
    % md_initial.smb.mass_balance   = md_initial.results.TransientSolution(end).SmbMassBalance;
    % md_initial.mask.ocean_levelset = md_initial.results.TransientSolution(end).MaskOceanLevelset;

    % Update true state fields
    md_true.geometry.thickness = md_true.results.TransientSolution(end).Thickness;
    md_true.geometry.surface   = md_true.results.TransientSolution(end).Surface;
    md_true.geometry.base      = md_true.results.TransientSolution(end).Base;
    md_true.initialization.vx  = md_true.results.TransientSolution(end).Vx;
    md_true.initialization.vy  = md_true.results.TransientSolution(end).Vy;
    md_true.initialization.vel = md_true.results.TransientSolution(end).Vel;
    md_true.initialization.pressure = md_true.results.TransientSolution(end).Pressure;
    md_true.smb.mass_balance   = md_true.results.TransientSolution(end).SmbMassBalance;
    md_true.mask.ocean_levelset = md_true.results.TransientSolution(end).MaskOceanLevelset;

    % Update nurged state fields
    md_wrong.geometry.thickness = md_wrong.results.TransientSolution(end).Thickness;
    md_wrong.geometry.surface   = md_wrong.results.TransientSolution(end).Surface;
    md_wrong.geometry.base      = md_wrong.results.TransientSolution(end).Base;
    md_wrong.initialization.vx  = md_wrong.results.TransientSolution(end).Vx;
    md_wrong.initialization.vy  = md_wrong.results.TransientSolution(end).Vy;
    md_wrong.initialization.vel = md_wrong.results.TransientSolution(end).Vel;
    md_wrong.initialization.pressure = md_wrong.results.TransientSolution(end).Pressure;
    md_wrong.smb.mass_balance   = md_wrong.results.TransientSolution(end).SmbMassBalance;
    md_wrong.mask.ocean_levelset = md_wrong.results.TransientSolution(end).MaskOceanLevelset;
% Define interpolation grid
    [X, Y] = meshgrid(1:1000:700000, 1:500:80000); % Same grid for all
    X_km = X / 1000;
    Y_km = Y / 1000;

    % Interpolate model fields for initial state
    bed_grid_initial = griddata(md_initial.mesh.x, md_initial.mesh.y, md_initial.geometry.bed, X, Y, 'linear');
    base_grid_initial = griddata(md_initial.mesh.x, md_initial.mesh.y, md_initial.geometry.base, X, Y, 'linear');
    surface_grid_initial = griddata(md_initial.mesh.x, md_initial.mesh.y, md_initial.geometry.surface, X, Y, 'linear');

    % Interpolate model fields for true state
    bed_grid_true = griddata(md_true.mesh.x, md_true.mesh.y, md_true.geometry.bed, X, Y, 'linear');
    base_grid_true = griddata(md_true.mesh.x, md_true.mesh.y, md_true.geometry.base, X, Y, 'linear');
    surface_grid_true = griddata(md_true.mesh.x, md_true.mesh.y, md_true.geometry.surface, X, Y, 'linear');

    % Interpolate model fields for nurged state
    bed_grid_wrong = griddata(md_wrong.mesh.x, md_wrong.mesh.y, md_wrong.geometry.bed, X, Y, 'linear');
    base_grid_wrong = griddata(md_wrong.mesh.x, md_wrong.mesh.y, md_wrong.geometry.base, X, Y, 'linear');
    surface_grid_wrong = griddata(md_wrong.mesh.x, md_wrong.mesh.y, md_wrong.geometry.surface, X, Y, 'linear');

    % Create figure with three subplots stacked vertically
    fig = figure('Color', 'w', 'Position', [100, 100, 2400, 3600]); % Adjusted for vertical layout

    % Subplot 1: Initial State (top)
    ax1 = subplot(3, 1, 1);
    plot_subplot(ax1, X_km, Y_km, bed_grid_initial, base_grid_initial, surface_grid_initial, 'Initial State');

    % Subplot 2: True State (middle)
    ax2 = subplot(3, 1, 2);
    plot_subplot(ax2, X_km, Y_km, bed_grid_true, base_grid_true, surface_grid_true, 'True State');

    % Subplot 3: Nurged State (bottom)
    ax3 = subplot(3, 1, 3);
    plot_subplot(ax3, X_km, Y_km, bed_grid_wrong, base_grid_wrong, surface_grid_wrong, 'Nurged State');

    % Add a single colorbar for all subplots
    cbar = colorbar(ax3, 'Location', 'eastoutside');
    cbar.Label.String = 'Elevation (m)';
    cbar.Label.FontSize = 16;
    cbar.Position = [0.92, 0.15, 0.02, 0.7]; % Spans the height of all subplots
    colormap(fig, jet); % Apply colormap to the entire figure
    caxis([-1500, 2500]); % Apply caxis to all subplots

    % Adjust layout and add a super title for the presentation
    sgtitle('Comparison of Ice Sheet Model States', 'FontSize', 24);
    set(gcf, 'Position', [100, 100, 2400, 3600]);
% end

% Helper function to plot a subplot
function plot_subplot(ax, X_km, Y_km, bed_grid, base_grid, surface_grid, title_str)
    hold on;
    surf(ax, X_km, Y_km, bed_grid, 'EdgeColor', 'none', 'FaceColor', 'interp');
    surf(ax, X_km, Y_km, base_grid, 'EdgeColor', 'none', 'FaceColor', 'interp');
    surf_handle = surf(ax, X_km, Y_km, surface_grid, ...
        'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 0.6);
    view(ax, 7, 10);
    grid on;
    daspect([16, 3.5, 180]);
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.ZColor = 'k';
    set(ax, 'BoxStyle', 'full', 'Box', 'off');
    xlabel(ax, 'x (km)', 'FontSize', 16);
    ylabel(ax, 'y (km)', 'FontSize', 16);
    zlabel(ax, 'z (m)', 'FontSize', 16);
    title(ax, title_str, 'FontSize', 18);
    set(ax, 'FontSize', 16);
    camlight headlight;
    lighting gouraud;
    hold off;
end
%     % Define interpolation grid
%     [X, Y] = meshgrid(1:1000:700000, 1:500:80000); % Same grid for all
%     X_km = X / 1000;
%     Y_km = Y / 1000;
% 
%     % Interpolate model fields for initial state
%     bed_grid_initial = griddata(md_initial.mesh.x, md_initial.mesh.y, md_initial.geometry.bed, X, Y, 'linear');
%     base_grid_initial = griddata(md_initial.mesh.x, md_initial.mesh.y, md_initial.geometry.base, X, Y, 'linear');
%     surface_grid_initial = griddata(md_initial.mesh.x, md_initial.mesh.y, md_initial.geometry.surface, X, Y, 'linear');
% 
%     % Interpolate model fields for true state
%     bed_grid_true = griddata(md_true.mesh.x, md_true.mesh.y, md_true.geometry.bed, X, Y, 'linear');
%     base_grid_true = griddata(md_true.mesh.x, md_true.mesh.y, md_true.geometry.base, X, Y, 'linear');
%     surface_grid_true = griddata(md_true.mesh.x, md_true.mesh.y, md_true.geometry.surface, X, Y, 'linear');
% 
%     % Interpolate model fields for nurged state
%     bed_grid_wrong = griddata(md_wrong.mesh.x, md_wrong.mesh.y, md_wrong.geometry.bed, X, Y, 'linear');
%     base_grid_wrong = griddata(md_wrong.mesh.x, md_wrong.mesh.y, md_wrong.geometry.base, X, Y, 'linear');
%     surface_grid_wrong = griddata(md_wrong.mesh.x, md_wrong.mesh.y, md_wrong.geometry.surface, X, Y, 'linear');
% 
%     % Create figure with three subplots stacked vertically
%     fig = figure('Color', 'w', 'Position', [100, 100, 2400, 3600]); % Adjusted for vertical layout
% 
%     % Subplot 1: Initial State
%     ax1 = subplot(3, 1, 1);
%     plot_subplot(ax1, X_km, Y_km, bed_grid_initial, base_grid_initial, surface_grid_initial, 'Initial State');
% 
%     % Subplot 2: True State
%     ax2 = subplot(3, 1, 2);
%     plot_subplot(ax2, X_km, Y_km, bed_grid_true, base_grid_true, surface_grid_true, 'True State');
% 
%     % Subplot 3: Nurged State
%     ax3 = subplot(3, 1, 3);
%     plot_subplot(ax3, X_km, Y_km, bed_grid_wrong, base_grid_wrong, surface_grid_wrong, 'Nurged State');
% 
%     % Add a single colorbar for all subplots
%     cbar = colorbar(ax3, 'Location', 'eastoutside');
%     cbar.Label.String = 'Elevation (m)';
%     cbar.Label.FontSize = 16;
%     cbar.Position = [0.92, 0.15, 0.02, 0.7]; % Spans the height of all subplots
%     colormap(fig, jet); % Apply colormap to the entire figure
%     caxis([-1500, 2500]); % Apply caxis to all subplots
% 
%     % Adjust layout and add a super title for the presentation
%     sgtitle('Comparison of Ice Sheet Model States', 'FontSize', 24);
%     set(gcf, 'Position', [100, 100, 2400, 3600]);
% % end
% 
% % Helper function to plot a subplot
% function plot_subplot(ax, X_km, Y_km, bed_grid, base_grid, surface_grid, title_str)
%     hold on;
%     surf(ax, X_km, Y_km, bed_grid, 'EdgeColor', 'none', 'FaceColor', 'interp');
%     surf(ax, X_km, Y_km, base_grid, 'EdgeColor', 'none', 'FaceColor', 'interp');
%     surf_handle = surf(ax, X_km, Y_km, surface_grid, ...
%         'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 0.6);
%     view(ax, 7, 10);
%     grid on;
%     daspect([16, 3.5, 180]);
%     ax.XColor = 'k';
%     ax.YColor = 'k';
%     ax.ZColor = 'k';
%     set(ax, 'BoxStyle', 'full', 'Box', 'off');
%     xlabel(ax, 'x (km)', 'FontSize', 16);
%     ylabel(ax, 'y (km)', 'FontSize', 16);
%     zlabel(ax, 'z (m)', 'FontSize', 16);
%     title(ax, title_str, 'FontSize', 18);
%     set(ax, 'FontSize', 16);
%     camlight headlight;
%     lighting gouraud;
%     hold off;
% end