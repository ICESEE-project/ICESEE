function variable_size = initialize_model(rank, nprocs, ens_id)
    % Initialize ISSM model for MISMIP-like experiment
    % Inputs: rank, nprocs (MPI settings), ens_id (ensemble ID)

    % Read kwargs from .mat file
    model_kwargs = sprintf('model_kwargs_%d.mat', ens_id);
    kwargs       = load(model_kwargs);
    ParamFile    = char(kwargs.ParamFile);
    Lx           = double(kwargs.Lx); % 640000 m
    Ly           = double(kwargs.Ly); % 80000 m
    cluster_name = char(kwargs.cluster_name);
    steps        = double(kwargs.steps);
    icesee_path  = char(kwargs.icesee_path);
    data_path    = char(kwargs.data_path);
    devmode      = logical(kwargs.devmode); % Development mode flag

    % get the current working directory
    cwd = pwd;
    [issmroot,~,~]=fileparts(fileparts(cwd));
    if devmode
        newpath=fullfile(issmroot,'/src/m/dev');
        addpath(newpath);
        devpath;
    end

    if iscell(kwargs.vec_inputs)
        vec_inputs = kwargs.vec_inputs;
    elseif isstring(kwargs.vec_inputs)
        vec_inputs = cellstr(kwargs.vec_inputs(:));
    elseif ischar(kwargs.vec_inputs)
        vec_inputs = cellstr(kwargs.vec_inputs);
    else
        error('Unsupported type for kwargs.vec_inputs');
    end

    vec_inputs = reshape(vec_inputs, 1, []);
    vec_inputs = cellfun(@strtrim, vec_inputs, 'UniformOutput', false);

    folder = sprintf('./Models/ens_id_%d', ens_id);
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    % disp(['[MATLAB] Initializing model with rank: ', num2str(rank), ', nprocs: ', num2str(nprocs), ', ens_id: ', num2str(ens_id)]);

	% steps = [1:4]; 
    steps = [55];
    % steps = [1:9];
    % steps=[83];

    % clear all; close all;

    % ens_id = 0; only set to zero during spin-up

    folder = sprintf('./Models');
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    % Mesh generation (Step 1)
    if any(steps == 1)
        md = model();

        % Define the domain outline (x, y coordinates, closing the contour)
        domain = [
            0, 0, 0
            0, Lx, 0;      % Point 1
            Lx, Lx, Ly;    % Point 2
            Lx, 0, Ly;     % Point 3
            0, 0, 0;       % Point 4
        ];
        
        % Define the output filename
        filename = 'Domain.exp';
        
        % Open file for writing
        fid = fopen(filename, 'w');
        if fid == -1
            error('Unable to open %s for writing', filename);
        end
        
        % Write header
        fprintf(fid, '## Name:DomainOutline\n');
        fprintf(fid, '## Icon:0\n');
        fprintf(fid, '# Points Count  Value\n');
        fprintf(fid, '%d 1\n', size(domain, 1));
        fprintf(fid, '# X pos Y pos\n');
        
        % Write points (using columns 2 and 3 for x, y)
        for i = 1:size(domain, 1)
            fprintf(fid, '%f %f\n', domain(i, 2), domain(i, 3));
        end
        
        % Close file
        fclose(fid);

        % Use bamg for variable-resolution mesh (500 m to 10 km)
        % md = bamg(md, 'domain', './Domain.exp', 'hmax', 2000, 'splitcorners', 1);
        % hvertices=[10000;500;500;10000];
        hvertices=[10000;500;5000;7500];
        md = bamg(md, 'domain', 'Domain.exp', 'hvertices',hvertices);
        % md = bamg(md, 'domain', 'Domain.exp', 'hmax', 2000, 'splitcorners', 1);
        % md = triangle(md,'Domain.exp',27000);
        
        plotmodel(md,'data','mesh');
        % Save mesh
        filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
        save(filename, 'md');
    end

    % Parameterization (Step 2)
    if any(steps == 2)
        % filename = fullfile(folder, 'ISMIP.SetMask.mat');
        filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
        md = loadmodel(filename);
        md=setmask(md,'','');
        % md=setflowequation(md,'SSA','all');
        ParamFile = 'Mismip2_matlab.par';
        md = parameterize(md, ParamFile); % Use Mismip2.par
        
        filename = fullfile(folder, 'ISMIP.Parameterization.mat');
        save(filename, 'md');

        % write_netCDF(md, 'ISMIP_Parameterization.nc');
        % export_netCDF(md,'ISMIP_Parameterization.nc');
        plotmodel(md, 'data', md.geometry.bed, 'title', 'Ice bed_t=0'); 
        % plotmodel(md, 'data', md.initialization.vel, 'title', 'Ice velocity');
        % plotmodel(md, 'data', md.friction.coefficient, 'title', 'Ice friction');
    end

    % Adding bed roughness
    if any(steps == 3)
        filename = fullfile(folder, 'ISMIP.Parameterization.mat');
        md = loadmodel(filename);
        
        % plotmodel(md, 'data', md.geometry.bed);

        % read in bed roughness data
        icesee_path='/Users/bkyanjo3/da_project/ICESEE/applications/issm_model/examples/ISMIP_Choi';
        filename = fullfile(icesee_path,'_modelrun_datasets/', 'friction_bed_0.h5');
        % filename = fullfile('friction_bed_0.h5');
        bed = h5read(filename, '/bed');
        nsize = md.mesh.numberofvertices;
        % md.geometry.bed = md.geometry.bed + bed(1:nsize) -  284 * rand(nsize,1);
        md.geometry.bed = md.geometry.bed + bed(1:nsize) -  284 * ones(nsize,1);
        % pos = (md.mesh.x > 375e3) & (md.mesh.x < 595e3) & (md.mesh.y > 65e3);
        % pos = (md.mesh.x < 640e3) & (md.mesh.y < 20e3);
        % md.geometry.bed(pos) = md.geometry.bed(pos) - 200;
        % md.geometry.base = md.geometry.base + bed(1:nsize);

        plotmodel(md, 'data', md.geometry.bed);
       
        filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
        save(filename, 'md');
    end

    % solve steady state
    if any(steps == 4)
        filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
        md = loadmodel(filename);

        md=setflowequation(md,'SSA','all');

        % Time stepping
        md.timestepping=timesteppingadaptive();
        md.timestepping.time_step_max=100;
        md.timestepping.time_step_min=0.1;
        md.timestepping.final_time=5000;

        md.settings.output_frequency=100;
        md.stressbalance.maxiter=100;
        md.stressbalance.abstol = NaN;
        md.stressbalance.restol = 1;
        md.settings.solver_residue_threshold=5e-5;
        
        % Verbose
        md.verbose = verbose('all');

        md.cluster=generic('name',oshostname(),'np',4);
        md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
        md = solve(md, 'Transient','runtimename',false);

        filename = fullfile(folder, 'Transient_steadystate_0.mat');
        save(filename, 'md');
    end

    % solve steady state
    if any(steps == 5)
        filename = fullfile(folder, 'Transient_steadystate_0.mat');
        md_ref = loadmodel(filename);

        filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
        md = loadmodel(filename);
        md=setflowequation(md,'SSA','all');
        
        md.geometry.thickness = md_ref.results.TransientSolution(end).Thickness;
        md.geometry.surface   = md_ref.results.TransientSolution(end).Surface;
        md.geometry.base      = md_ref.results.TransientSolution(end).Base;

        % Update other fields
        md.initialization.vx        = md_ref.results.TransientSolution(end).Vx;
        md.initialization.vy        = md_ref.results.TransientSolution(end).Vy;
        md.initialization.vel       = md_ref.results.TransientSolution(end).Vel;
        md.initialization.pressure  = md_ref.results.TransientSolution(end).Pressure;
        md.smb.mass_balance         = md_ref.results.TransientSolution(end).SmbMassBalance;
        md.mask.ocean_levelset      = md_ref.results.TransientSolution(end).MaskOceanLevelset;

        md.geometry.bed = md_ref.results.TransientSolution(end).Bed;
        md.friction.coefficient = md_ref.results.TransientSolution(end).FrictionCoefficient;
        
        % Time stepping
        md.timestepping=timesteppingadaptive();
        md.timestepping.time_step_max=100;
        md.timestepping.time_step_min=0.1;
        md.timestepping.final_time=5000;

        md.settings.output_frequency=100;
        md.stressbalance.maxiter=100;
        md.stressbalance.abstol = NaN;
        md.stressbalance.restol = 1;
        md.settings.solver_residue_threshold=5e-5;
        
        % Verbose
        md.verbose = verbose('all');

        md.cluster=generic('name',oshostname(),'np',4);
        md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
        md = solve(md, 'Transient','runtimename',false);

        filename = fullfile(folder, 'Transient_steadystate_00.mat');
        save(filename, 'md');
    end

    % Steady state simulation 2
    if any(steps == 6)

        filename = fullfile(folder, 'Transient_steadystate_00.mat');
        md_ref = loadmodel(filename);

        filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
        md = loadmodel(filename);
        md=setflowequation(md,'SSA','all');
        
        md.geometry.thickness = md_ref.results.TransientSolution(end).Thickness;
        md.geometry.surface   = md_ref.results.TransientSolution(end).Surface;
        md.geometry.base      = md_ref.results.TransientSolution(end).Base;

        % Update other fields
        md.initialization.vx        = md_ref.results.TransientSolution(end).Vx;
        md.initialization.vy        = md_ref.results.TransientSolution(end).Vy;
        md.initialization.vel       = md_ref.results.TransientSolution(end).Vel;
        md.initialization.pressure  = md_ref.results.TransientSolution(end).Pressure;
        md.smb.mass_balance         = md_ref.results.TransientSolution(end).SmbMassBalance;
        md.mask.ocean_levelset      = md_ref.results.TransientSolution(end).MaskOceanLevelset;

        md.geometry.bed = md_ref.results.TransientSolution(end).Bed;
        md.friction.coefficient = md_ref.results.TransientSolution(end).FrictionCoefficient;
        
        md.timestepping=timesteppingadaptive();
        md.timestepping.final_time=7500;
        md.settings.output_frequency=100;
        md.stressbalance.maxiter=100;

        md.stressbalance.abstol = NaN;
        md.stressbalance.restol = 1;

        md.verbose = verbose('all');

        md.cluster=generic('name',oshostname(),'np',4);
        md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
        md = solve(md, 'Transient','runtimename',false);

        filename = fullfile(folder, 'Transient_steadystate_1.mat');
        save(filename, 'md');
    end

    % Steady state simulation 2
    if any(steps == 7)

        filename = fullfile(folder, 'Transient_steadystate_1.mat');
        md_ref = loadmodel(filename);

        filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
        md = loadmodel(filename);
        md=setflowequation(md,'SSA','all');
        
        md.geometry.thickness = md_ref.results.TransientSolution(end).Thickness;
        md.geometry.surface   = md_ref.results.TransientSolution(end).Surface;
        md.geometry.base      = md_ref.results.TransientSolution(end).Base;

        % Update other fields
        md.initialization.vx        = md_ref.results.TransientSolution(end).Vx;
        md.initialization.vy        = md_ref.results.TransientSolution(end).Vy;
        md.initialization.vel       = md_ref.results.TransientSolution(end).Vel;
        md.initialization.pressure  = md_ref.results.TransientSolution(end).Pressure;
        md.smb.mass_balance         = md_ref.results.TransientSolution(end).SmbMassBalance;
        md.mask.ocean_levelset      = md_ref.results.TransientSolution(end).MaskOceanLevelset;

        md.geometry.bed = md_ref.results.TransientSolution(end).Bed;
        md.friction.coefficient = md_ref.results.TransientSolution(end).FrictionCoefficient;

        % filename = fullfile(folder, 'Transient_steadystate_11.mat');
        % save(filename, 'md');
        % md = loadmodel(filename);

        md.timestepping=timesteppingadaptive();
        md.timestepping.final_time=7500;
        md.settings.output_frequency=100;
        md.stressbalance.maxiter=100;

        md.stressbalance.abstol = NaN;
        md.stressbalance.restol = 1;

        md.verbose = verbose('all');

        md.cluster=generic('name',oshostname(),'np',4);
        md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
        md = solve(md, 'Transient','runtimename',false);

        filename = fullfile(folder, 'Transient_steadystate_2.mat');
        save(filename, 'md');
    end
    % Steady state simulation 2
    if any(steps == 8)

        filename = fullfile(folder, 'Transient_steadystate_2.mat');
        md_ref = loadmodel(filename);

        filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
        md = loadmodel(filename);
        md=setflowequation(md,'SSA','all');
        
        md.geometry.thickness = md_ref.results.TransientSolution(end).Thickness;
        md.geometry.surface   = md_ref.results.TransientSolution(end).Surface;
        md.geometry.base      = md_ref.results.TransientSolution(end).Base;

        % Update other fields
        md.initialization.vx        = md_ref.results.TransientSolution(end).Vx;
        md.initialization.vy        = md_ref.results.TransientSolution(end).Vy;
        md.initialization.vel       = md_ref.results.TransientSolution(end).Vel;
        md.initialization.pressure  = md_ref.results.TransientSolution(end).Pressure;
        md.smb.mass_balance         = md_ref.results.TransientSolution(end).SmbMassBalance;
        md.mask.ocean_levelset      = md_ref.results.TransientSolution(end).MaskOceanLevelset;

        md.geometry.bed = md_ref.results.TransientSolution(end).Bed;
        md.friction.coefficient = md_ref.results.TransientSolution(end).FrictionCoefficient;

        % filename = fullfile(folder, 'Transient_steadystate_21.mat');
        % save(filename, 'md');
        % md = loadmodel(filename);

        md.timestepping=timestepping();
        md.timestepping.time_step=0.1;
        md.timestepping.start_time=0;
        md.timestepping.final_time=200;
        md.settings.output_frequency=10;
        md.stressbalance.maxiter=100;

        md.stressbalance.abstol = NaN;
        md.stressbalance.restol = 1;

        md.verbose = verbose('all');

        md.cluster=generic('name',oshostname(),'np',4);
        md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
        md = solve(md, 'Transient','runtimename',false);

        filename = fullfile(folder, 'reference_simulation.mat');
        save(filename, 'md');
    end

    if any(steps == 9)
        % clear all; close all;
        ens_id = 0;
        % folder = sprintf('./Models/ens_id_%d', ens_id);
        % if ~exist(folder, 'dir')
        %     mkdir(folder);
        % end
        % filename = fullfile(folder, 'Transient_steadystate1.mat');

        filename = fullfile(folder, 'reference_simulation.mat');
        md = loadmodel(filename);
        md_ref = loadmodel(filename);

        filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
        md = loadmodel(filename);
        md=setflowequation(md,'SSA','all');
        
        md.geometry.thickness = md_ref.results.TransientSolution(end).Thickness;
        md.geometry.surface   = md_ref.results.TransientSolution(end).Surface;
        md.geometry.base      = md_ref.results.TransientSolution(end).Base;

        % Update other fields
        md.initialization.vx        = md_ref.results.TransientSolution(end).Vx;
        md.initialization.vy        = md_ref.results.TransientSolution(end).Vy;
        md.initialization.vel       = md_ref.results.TransientSolution(end).Vel;
        md.initialization.pressure  = md_ref.results.TransientSolution(end).Pressure;
        md.smb.mass_balance         = md_ref.results.TransientSolution(end).SmbMassBalance;
        md.mask.ocean_levelset      = md_ref.results.TransientSolution(end).MaskOceanLevelset;

        md.geometry.bed = md_ref.results.TransientSolution(end).Bed;
        md.friction.coefficient = md_ref.results.TransientSolution(end).FrictionCoefficient;

        filename = fullfile(folder, 'ISMIP.reference_simulation_0.mat');
        save(filename, 'md','-v7.3');
        
        md.smb.mass_balance=-0.3*ones(md.mesh.numberofvertices,1);
        % md.smb.mass_balance=-1.5*ones(md.mesh.numberofvertices,1);
        md.transient.ismovingfront=0;
        % 
        md.basalforcings=linearbasalforcings();
        md.basalforcings.deepwater_melting_rate=150;
        md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
        
        % 
        md.timestepping=timestepping();
        md.timestepping.start_time = 0;
        md.timestepping.time_step = 0.2;
        md.timestepping.final_time = 200;
        md.settings.output_frequency = 10;
        md.stressbalance.maxiter = 100;
        md.stressbalance.restol = 1;
        md.stressbalance.reltol = 0.001;
        md.stressbalance.abstol = NaN;
        md.settings.solver_residue_threshold=5e-2;
        % 
        % md.verbose = verbose('all');
        % 
        % 
        md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
        md.cluster=generic('name',oshostname(),'np',4);
        md.settings.waitonlock = 1;
        md = solve(md, 'Transient', 'runtimename',false);

        % % export_netCDF(md,'Transient_steadystate4.nc');
        % 
        filename = fullfile(folder, 'ISMIP.reference_simulation1.mat');
        save(filename, 'md', '-v7.3');

        % save file for ICESEE runs
        md_ref = loadmodel(filename);

        filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
        md = loadmodel(filename);
        md=setflowequation(md,'SSA','all');
        
        md.geometry.thickness = md_ref.results.TransientSolution(end).Thickness;
        md.geometry.surface   = md_ref.results.TransientSolution(end).Surface;
        md.geometry.base      = md_ref.results.TransientSolution(end).Base;

        % Update other fields
        md.initialization.vx        = md_ref.results.TransientSolution(end).Vx;
        md.initialization.vy        = md_ref.results.TransientSolution(end).Vy;
        md.initialization.vel       = md_ref.results.TransientSolution(end).Vel;
        md.initialization.pressure  = md_ref.results.TransientSolution(end).Pressure;
        md.smb.mass_balance         = md_ref.results.TransientSolution(end).SmbMassBalance;
        md.mask.ocean_levelset      = md_ref.results.TransientSolution(end).MaskOceanLevelset;

        md.geometry.bed = md_ref.results.TransientSolution(end).Bed;
        md.friction.coefficient = md_ref.results.TransientSolution(end).FrictionCoefficient;


        filename = fullfile(folder, 'ISMIP.reference_simulation_1.mat');
        save(filename, 'md','-v7.3');
    
    end

    if any(steps == 83)
        % clear all; close all;
        ens_id = 0;
        % folder = sprintf('./Models/ens_id_%d', ens_id);
        % if ~exist(folder, 'dir')
        %     mkdir(folder);
        % end

        filename = fullfile(folder, 'ISMIP.reference_simulation1.mat');
        % filename = fullfile(folder, 'reference_simulation.mat');
        % filename = fullfile(folder, 'Transient_steadystate_00.mat');
        % filename = fullfile(folder, 'Transient_steadystate_2.mat');
        % filename = fullfile(folder, 'Transient_steadystate_2.mat');
        % filename = fullfile(folder, 'Transient_steadystate_2.mat');
        % data_path='/Users/bkyanjo3/da_project/ICESEE/applications/issm_model/examples/ISMIP_Choi/data';
        % filename=fullfile(data_path,'ISMIP.reference_simulation2.mat');
        md = loadmodel(filename);
        end_=1;
        md.geometry.thickness = md.results.TransientSolution(end).Thickness;
        md.geometry.surface   = md.results.TransientSolution(end).Surface;
        md.geometry.base      = md.results.TransientSolution(end).Base;

        % Update other fields
        md.initialization.vx        = md.results.TransientSolution(end).Vx;
        md.initialization.vy        = md.results.TransientSolution(end).Vy;
        md.initialization.vel       = md.results.TransientSolution(end).Vel;
        md.initialization.pressure  = md.results.TransientSolution(end).Pressure;
        md.smb.mass_balance         = md.results.TransientSolution(end).SmbMassBalance;
        md.mask.ocean_levelset      = md.results.TransientSolution(end).MaskOceanLevelset;

        % md.friction.coefficient = md.results.TransientSolution(end).FrictionCoefficient;
        % md.geometry.bed      = md.results.TransientSolution(end).Bed;

        filename = fullfile(folder, 'ISMIP.initial_reference.mat');
        save(filename, 'md','-v7.3');

        md.smb.mass_balance=-0.3*ones(md.mesh.numberofvertices,1);
        % md.smb.mass_balance=-1.5*ones(md.mesh.numberofvertices,1);
        md.transient.ismovingfront=0;
        % 
        md.basalforcings=linearbasalforcings();
        md.basalforcings.deepwater_melting_rate=150;
        md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
        
        % 
        md.timestepping=timestepping();
        md.timestepping.start_time = 0;
        md.timestepping.time_step = 0.1;
        md.timestepping.final_time = 0.2;
        md.settings.output_frequency = 1;
        md.stressbalance.maxiter = 100;
        md.stressbalance.restol = 1;
        md.stressbalance.reltol = 0.001;
        md.stressbalance.abstol = NaN;
        md.settings.solver_residue_threshold=5e-2;
        % 
        % md.verbose = verbose('all');
        % 
        % 
        md.cluster=generic('name',oshostname(),'np',4);
        md.settings.waitonlock = 1;
        md = solve(md, 'Transient', 'runtimename',false);

        filename = fullfile(folder, 'ISMIP.true_simulation.mat');
        % filename = fullfile(folder, 'ISMIP.initial_reference.mat');
        save(filename, 'md','-v7.3');

        md =loadmodel(fullfile(folder, 'ISMIP.initial_reference.mat'));
        % export_netCDF(md,'Transient_steadystate4.nc')


        % Define interpolation grid
        [X, Y] = meshgrid(1:1000:700000, 1:500:80000);  % Similar to Python np.meshgrid
        % [X, Y] = meshgrid(1:1000:640000, 1:500:80000);
        
        % Interpolate model fields onto regular grid
        bed_grid = griddata(md.mesh.x, md.mesh.y, md.geometry.bed, X, Y, 'linear');
        base_grid = griddata(md.mesh.x, md.mesh.y, md.geometry.base, X, Y, 'linear');
        surface_grid = griddata(md.mesh.x, md.mesh.y, md.geometry.surface, X, Y, 'linear');
        
        % Convert to kilometers
        X_km = X / 1000;
        Y_km = Y / 1000;
        
        % Create figure
        fig = figure('Color', 'w', 'Position', [100, 100, 2000, 1200]);
        % fig = figure('Position', [100, 100, 1200, 720]);
        ax = axes('Parent', fig);
        hold on;
        
        % Plot surfaces
        surf(ax, X_km, Y_km, bed_grid, 'EdgeColor', 'none', 'FaceColor', 'interp');
        surf(ax, X_km, Y_km, base_grid, 'EdgeColor', 'none', 'FaceColor', 'interp');
        surf_handle = surf(ax, X_km, Y_km, surface_grid, ...
            'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 0.6);
        
        % Set color map and color limits
        colormap(jet);
        caxis([-1500, 2500]);
        
        % Add colorbar
        cbar = colorbar('Location', 'eastoutside');
        cbar.Label.String = 'Elevation (m)';
        cbar.Label.FontSize = 14;
        
        % Set viewing angle
        % view(ax, 30, -60);
        view(ax, 7, 10);

        grid on;

        % set 
        
        % Set aspect ratio (roughly matching Python's ax.set_box_aspect)
        % daspect([8, 0.02, 3]);
        % daspect([16, 3.5, 180]);
        
        % Remove pane fills (simulates transparency)
        ax.XColor = 'k';
        ax.YColor = 'k';
        ax.ZColor = 'k';
        set(ax, 'BoxStyle', 'full', 'Box', 'off');
        
        % Ticks
        % yticks([0, 40, 80]);
        
        % Axis labels
        xlabel(ax, 'x (km)', 'FontSize', 14);
        ylabel(ax, 'y (km)', 'FontSize', 14);
        zlabel(ax, 'z (m)', 'FontSize', 14);
        
        % Axis label padding and alignment (not identical to Python but close)
        ax.LabelFontSizeMultiplier = 1.2;
        ax.TitleFontSizeMultiplier = 1.2;
        
        % Set font size of all tick labels
        set(ax, 'FontSize', 14);
        
        % Optional: Lighting (MATLAB has lighting tools but not pane edge controls like matplotlib)
        camlight headlight;
        lighting gouraud;
        
        hold off;

    end

    %  initialize from the reference simulation (Step 5)
    use_reference_data = true; % Flag to use reference data
    if any(steps == 55)
    % if use_reference_data
        reference_data = char(kwargs.reference_data); % Path to reference data
        ens_id_init = 0;
        folder = sprintf('./Models/ens_id_%d',  ens_id_init);
        if ~exist(folder, 'dir')
            mkdir(folder);
        end

        filename = fullfile(folder,reference_data);
        md = loadmodel(filename);

        % Extract fields to save
        result_0 = md.initialization;
        result_1 = md.geometry;
        result_2 = md.friction;

        % save mesh cordinates
        % ens_id = 0;
        h5file = fullfile(icesee_path, data_path, sprintf('mesh_idxy_%d.h5', ens_id));

        % Always use column vectors for consistency
        x_param  = double(md.mesh.x(:));                      % [N x 1]
        y_param  = double(md.mesh.y(:));                      % [N x 1]
        fric_idx = double((1:md.mesh.numberofvertices).');    % [N x 1]

        % If file exists from a previous run (possibly with different sizes), remove it
        if exist(h5file, 'file')
            delete(h5file);
        end

        % Create datasets with the current sizes (N x 1)
        % size_x = size(x_param)
        % size_y = size(y_param)
        % size_fric_idx = size(fric_idx)
        h5create(h5file, '/fric_x',   size(x_param),  'Datatype', 'double');
        h5create(h5file, '/fric_y',   size(y_param),  'Datatype', 'double');
        h5create(h5file, '/fric_idx', size(fric_idx), 'Datatype', 'double');

        % Now write the data – shapes match exactly
        h5write(h5file, '/fric_x',   x_param);
        h5write(h5file, '/fric_y',   y_param);
        h5write(h5file, '/fric_idx', fric_idx);


        % 	% --- fetch and save data for ensemble use
		filename = fullfile(icesee_path, data_path, sprintf('ensemble_init_%d.h5', ens_id));
        
		% Ensure the directory exists
		[filepath, ~, ~] = fileparts(filename);
		if ~exist(filepath, 'dir')
			mkdir(filepath);
		end

        % Check if the file exists and delete it if it does
		if isfile(filename)
			delete(filename);
		end

        %  save the fields to the file
        data = { 'Thickness', result_1, 'thickness';
                % 'Base', result_0, 'Base';
                'Surface', result_1, 'surface';
                'Vx', result_0, 'vx';
                'Vy', result_0, 'vy';
                'bed', result_1, 'bed';
                'coefficient', result_2, 'coefficient'
        };
                
        writeToHDF5(filename, data);

    end

end


function writeToHDF5(filename, data)
    % WRITETOHDF5 Writes variables to an HDF5 file.
    % Inputs:
    %   filename - Name of the HDF5 file
    %   data - Cell array with columns: {var_name, source_object, field_name}

    [filepath, ~, ~] = fileparts(filename);
    if ~exist(filepath, 'dir')
        mkdir(filepath);
    end
    if isfile(filename)
        delete(filename);
    end
    
    for i = 1:size(data, 1)
        var_name = data{i, 1};
        var_value = data{i, 2}.(data{i, 3});
        h5create(filename, ['/' var_name], size(var_value));
        h5write(filename, ['/' var_name], var_value);
    end
end

