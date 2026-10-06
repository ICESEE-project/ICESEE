steps = [8];

if any(steps == 8)
  
    initial_state= fullfile("data", "ISMIP_initial_data.mat");
    true_state = fullfile("issm_data","true_state.mat");
    nurged_state = fullfile("issm_data","nurged_state.mat");
    md_true = loadmodel(true_state);
    md_wrong = loadmodel(nurged_state);

    % md = transientrestart(md);
    % end_ = 700;
    md_true.geometry.thickness = md_true.results.TransientSolution(end).Thickness;
    md_true.geometry.surface   = md_true.results.TransientSolution(end).Surface;
    md_true.geometry.base      = md_true.results.TransientSolution(end).Base;

    % _truemd_trueUpdate other fields
    md_true.initialization.vx        = md_true.results.TransientSolution(end).Vx;
    md_true.initialization.vy        = md_true.results.TransientSolution(end).Vy;
    md_true.initialization.vel       = md_true.results.TransientSolution(end).Vel;
    md_true.initialization.pressure  = md_true.results.TransientSolution(end).Pressure;
    md_true.smb.mass_balance         = md_true.results.TransientSolution(end).SmbMassBalance;
    md_true.mask.ocean_levelset      = md_true.results.TransientSolution(end).MaskOceanLevelset
    
    % nurged_state  upddate
    md_wrong.geometry.thickness = md_wrong.results.TransientSolution(end).Thickness;
    md_wrong.geometry.surface   = md_wrong.results.TransientSolution(end).Surface;
    md_wrong.geometry.base      = md_wrong.results.TransientSolution(end).Base;

    % _wrongmd_wrongUpdate other fields
    md_wrong.initialization.vx        = md_wrong.results.TransientSolution(end).Vx;
    md_wrong.initialization.vy        = md_wrong.results.TransientSolution(end).Vy;
    md_wrong.initialization.vel       = md_wrong.results.TransientSolution(end).Vel;
    md_wrong.initialization.pressure  = md_wrong.results.TransientSolution(end).Pressure;
    md_wrong.smb.mass_balance         = md_wrong.results.TransientSolution(end).SmbMassBalance;
    md_wrong.mask.ocean_levelset      = md_wrong.results.TransientSolution(end).MaskOceanLevelset;

    % Export the models to netCDF files
    % check if true_state.nc and nurged_state.nc don't exist
    % if exist('true_state.nc', 'file') || exist('nurged_state.nc', 'file')
    %     disp('NetCDF files already exist. Skipping export.');
    % else
    %     export_netCDF(md_true, 'true_state.nc');
    %     export_netCDF(md_wrong, 'nurged_state.nc');
    % end

    % Define interpolation grid
    [X, Y] = meshgrid(1:1000:700000, 1:500:80000); % Same grid for both
    X_km = X / 1000;
    Y_km = Y / 1000;

    % Interpolate model fields for true state
    bed_grid_true = griddata(md_true.mesh.x, md_true.mesh.y, md_true.geometry.bed, X, Y, 'linear');
    base_grid_true = griddata(md_true.mesh.x, md_true.mesh.y, md_true.geometry.base, X, Y, 'linear');
    surface_grid_true = griddata(md_true.mesh.x, md_true.mesh.y, md_true.geometry.surface, X, Y, 'linear');

    % Interpolate model fields for nurged state
    bed_grid_wrong = griddata(md_wrong.mesh.x, md_wrong.mesh.y, md_wrong.geometry.bed, X, Y, 'linear');
    base_grid_wrong = griddata(md_wrong.mesh.x, md_wrong.mesh.y, md_wrong.geometry.base, X, Y, 'linear');
    surface_grid_wrong = griddata(md_wrong.mesh.x, md_wrong.mesh.y, md_wrong.geometry.surface, X, Y, 'linear');

    % Create figure with two subplots
    fig = figure('Color', 'w', 'Position', [100, 100, 2400, 1200]);

    % Subplot 1: True State
    ax1 = subplot(1, 2, 1);
    hold on;
    surf(ax1, X_km, Y_km, bed_grid_true, 'EdgeColor', 'none', 'FaceColor', 'interp');
    surf(ax1, X_km, Y_km, base_grid_true, 'EdgeColor', 'none', 'FaceColor', 'interp');
    surf_handle1 = surf(ax1, X_km, Y_km, surface_grid_true, ...
        'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 0.6);
    colormap(ax1, jet);
    caxis([-1500, 2500]);
    cbar1 = colorbar('Location', 'eastoutside');
    cbar1.Label.String = 'Elevation (m)';
    cbar1.Label.FontSize = 14;
    view(ax1, 7, 10);
    grid on;
    daspect([16, 3.5, 180]);
    ax1.XColor = 'k';
    ax1.YColor = 'k';
    ax1.ZColor = 'k';
    set(ax1, 'BoxStyle', 'full', 'Box', 'off');
    xlabel(ax1, 'x (km)', 'FontSize', 14);
    ylabel(ax1, 'y (km)', 'FontSize', 14);
    zlabel(ax1, 'z (m)', 'FontSize', 14);
    title(ax1, 'True State', 'FontSize', 16);
    set(ax1, 'FontSize', 14);
    camlight headlight;
    lighting gouraud;
    hold off;

    % Subplot 2: Nurged State
    ax2 = subplot(1, 2, 2);
    hold on;
    surf(ax2, X_km, Y_km, bed_grid_wrong, 'EdgeColor', 'none', 'FaceColor', 'interp');
    surf(ax2, X_km, Y_km, base_grid_wrong, 'EdgeColor', 'none', 'FaceColor', 'interp');
    surf_handle2 = surf(ax2, X_km, Y_km, surface_grid_wrong, ...
        'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 0.6);
    colormap(ax2, jet);
    caxis([-1500, 2500]);
    cbar2 = colorbar('Location', 'eastoutside');
    cbar2.Label.String = 'Elevation (m)';
    cbar2.Label.FontSize = 14;
    view(ax2, 7, 10);
    grid on;
    daspect([16, 3.5, 180]);
    ax2.XColor = 'k';
    ax2.YColor = 'k';
    ax2.ZColor = 'k';
    set(ax2, 'BoxStyle', 'full', 'Box', 'off');
    xlabel(ax2, 'x (km)', 'FontSize', 14);
    ylabel(ax2, 'y (km)', 'FontSize', 14);
    zlabel(ax2, 'z (m)', 'FontSize', 14);
    title(ax2, 'Nurged State', 'FontSize', 16);
    set(ax2, 'FontSize', 14);
    camlight headlight;
    lighting gouraud;
    hold off;

    % Adjust layout to prevent overlap
    set(gcf, 'Position', [100, 100, 2400, 1200]);
end