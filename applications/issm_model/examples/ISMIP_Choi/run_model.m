function run_model(data_fname, ens_id, rank, nprocs, k, dt, tinitial, tfinal)
    % Run ISSM model for transient simulation with ensemble support
    % Inputs: data_fname (output file name), ens_id (ensemble ID), rank, nprocs (MPI settings),
    %         k (time step index), dt (time step), tinitial, tfinal (time bounds)

    % Read kwargs from .mat file
    model_kwargs = sprintf('model_kwargs_%d.mat', ens_id);
    kwargs       = load(model_kwargs);
    cluster_name = char(kwargs.cluster_name);
    steps        = double(kwargs.steps);
    icesee_path  = char(kwargs.icesee_path);
    data_path    = char(kwargs.data_path);
    hpcmode      = logical(kwargs.hpcmode); % HPC mode flag
    devmode      = logical(kwargs.devmode); % Development mode flag

    reference_data = char(kwargs.reference_data);

    % Get the current working directory
    cwd = pwd;
    [issmroot,~,~] = fileparts(fileparts(cwd));
    if devmode
        newpath = fullfile(issmroot,'/src/m/dev');
        addpath(newpath);
        devpath;
    end

    % fprintf('[MATLAB] Running model with ens_id: %d, rank: %d, nprocs: %d, filename: %s\n', ens_id, rank, nprocs, data_fname);

    folder = sprintf('./Models/ens_id_%d', ens_id);
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    % Set up model for each EnKF stage
    if strcmp(data_fname, 'true_state.mat')
        % Special case for true state
        if k == 0 || isempty(k)
            % Initial run: load boundary conditions
            filename = fullfile(folder, reference_data);
            md = loadmodel(filename);

            md = setflowequation(md,'SSA','all');

            md.smb.mass_balance = -0.3 * ones(md.mesh.numberofvertices, 1); % m/yr
            
            md.basalforcings = linearbasalforcings();
            md.basalforcings.deepwater_melting_rate = 200; % m/yr
            md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices, 1);

            md.transient.ismovingfront = 0;   

            filename = fullfile(folder, data_fname);
            save(filename, 'md');
        end

    elseif strcmp(data_fname, 'nurged_state.mat')
        % Special case for ensemble input
        if k == 0 || isempty(k)
            % Initial run: load boundary conditions
            filename = fullfile(folder, reference_data);
            md = loadmodel(filename);

            md = setflowequation(md,'SSA','all');

            md.smb.mass_balance = -0.3 * ones(md.mesh.numberofvertices, 1); % m/yr
            
            md.basalforcings = linearbasalforcings();
            md.basalforcings.deepwater_melting_rate = 200; % m/yr
            md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices, 1);

            md.transient.ismovingfront = 0;   

            filename = fullfile(folder, data_fname);
            save(filename, 'md');
        end

    elseif strcmp(data_fname, 'initialize_ensemble.mat')
        % Special case for ensemble initialization
        if k == 0 || isempty(k)
            % Initial run: load boundary conditions
            filename = fullfile(folder, reference_data);
            md = loadmodel(filename);

            md.timestepping = timestepping();
            md.timestepping.start_time = tinitial;
            md.timestepping.time_step  = 0.1;
            md.timestepping.final_time = tfinal;

            % friction coefficient
            fcoeff = 2500*ones(md.mesh.numberofvertices,1);
            
            vx = md.initialization.vx;
            pos = find(vx ~= -888888);
            md.initialization.vx(pos) = vx(pos);
            vy = md.initialization.vy;
            pos = find(vy ~= -888888);
            md.initialization.vy(pos) = vy(pos);
            pos = find(fcoeff ~= -888888);
            md.friction.coefficient(pos) = fcoeff(pos); 
            
            thk_ref = md.geometry.thickness;
            pos = find(thk_ref ~= -888888);
            md.geometry.thickness(pos) = thk_ref(pos);
            bed_ref = md.geometry.bed;
            pos = find(bed_ref ~= -888888);
            md.geometry.bed(pos) = bed_ref(pos);
            base_ref = md.geometry.base;
            pos = find(base_ref ~= -888888);
            md.geometry.base(pos) = base_ref(pos);
            surf_ref = md.geometry.surface;
            pos = find(surf_ref ~= -888888);
            md.geometry.surface(pos) = surf_ref(pos);

            md.inversion.iscontrol=0;
        
            % Set minimum thickness to 1
            pos = find(md.geometry.thickness < 1);
            md.geometry.thickness(pos) = 1;

            % Calculate density ratio
            di = md.materials.rho_ice / md.materials.rho_water;

            % Calculate ocean levelset
            md.mask.ocean_levelset = md.geometry.thickness + md.geometry.bed / di;

            % Find positions where ocean_levelset < 0
            pos = find(md.mask.ocean_levelset < 0);

            % Update surface for floating ice
            md.geometry.surface(pos) = md.geometry.thickness(pos) * ...
                (md.materials.rho_water - md.materials.rho_ice) / md.materials.rho_water;

            % Update base
            md.geometry.base = md.geometry.surface - md.geometry.thickness;

            % Ensure base is not below bed
            pos = find(md.geometry.base < md.geometry.bed);
            md.geometry.base(pos) = md.geometry.bed(pos);

            % For grounded ice (ocean_levelset > 0)
            pos = find(md.mask.ocean_levelset > 0);
            md.geometry.base(pos) = md.geometry.bed(pos);

            % Final surface update
            md.geometry.surface = md.geometry.base + md.geometry.thickness;

            % Cluster setup
            if hpcmode
                md.settings.waitonlock = 0;
                md.cluster = generic('name', oshostname(), 'np', nprocs);
                md.cluster.codepath = [issmroot , '/bin'];
                md.cluster.login = 'arobel3';
                if ~exist(sprintf('%s/execution/color_%d', issmroot, ens_id), 'dir')
                    mkdir(sprintf('%s/execution/color_%d', issmroot, ens_id));
                end
                md.cluster.executionpath = sprintf('%s/execution/color_%d', issmroot, ens_id);
            else
                md.cluster = generic('name', cluster_name, 'np', nprocs);
                md.settings.waitonlock = 1;
                md.miscellaneous.name = sprintf('color_%d', ens_id);
            end

            % Verbose settings
            md.verbose = verbose('convergence', false, 'solution', true);

            % Solve transient
            md = solve(md, 'Transient');

            filename = fullfile(folder, data_fname);
            save(filename, 'md');

            % Save ensemble outputs in HDF5
            fields = {'Thickness', 'bed', 'coefficient'};
            result_0 = md.results.TransientSolution(end);
            result_1 = md.geometry;
            result_2 = md.friction;

            filename = fullfile(icesee_path, data_path, sprintf('ensemble_init_%d.h5', ens_id));

            data = {'Thickness', result_0, 'Thickness';
                    'Surface', result_0, 'Surface';
                    'bed', result_1, 'bed';
                    'coefficient', result_2, 'coefficient'};

            writeToHDF5(filename, data);

            % Break and return to avoid further processing
            return;
        end

    elseif strcmp(data_fname, 'enkf_state.mat')
        % Special case for ensemble assimilation
        if k == 0 || isempty(k)
            % Initial run: load boundary conditions
            filename = fullfile(folder, reference_data);
            md = loadmodel(filename);

            md = setflowequation(md,'SSA','all');

            thk_ref = md.geometry.thickness;
            bed_ref = md.geometry.bed;
            base_ref = md.geometry.base;

            md.smb.mass_balance = -0.3 * ones(md.mesh.numberofvertices, 1); % m/yr
            
            md.basalforcings = linearbasalforcings();
            md.basalforcings.deepwater_melting_rate = 200; % m/yr
            md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices, 1);

            md.transient.ismovingfront = 0;   
            
            md.transient.ismovingfront=0;
            md.transient.isthermal=0;
            md.transient.isstressbalance=1;
            md.transient.ismasstransport=1;
            md.transient.isgroundingline=1;
            md.groundingline.migration = 'SubelementMigration';
            md.groundingline.friction_interpolation='SubelementFriction1';
            md.groundingline.melt_interpolation='NoMeltOnPartiallyFloating';

            md.initialization.pressure = zeros(md.mesh.numberofvertices,1);
            md.masstransport.spcthickness = NaN*ones(md.mesh.numberofvertices,1);
    
            filename = fullfile(folder, data_fname);
            save(filename, 'md');

            % Save ensemble outputs in HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
            result_0 = md.results.TransientSolution(end);
            result_1 = md.geometry;
            result_2 = md.friction;

            data = {'Thickness', result_0, 'Thickness';
                    % 'Surface', result_0, 'Surface';
                    'bed', result_1, 'bed';
                    'coefficient', result_2, 'coefficient'};

            writeToHDF5(filename, data);
        end
    end

    if ismember(data_fname, {'true_state.mat', 'nurged_state.mat', 'enkf_state.mat'})
        % if k > 0
            % fprintf('[MATLAB ---] Running model for ensemble ID %d, step %d\n', ens_id, k);
            % Subsequent time steps: 
            filename = fullfile(folder, data_fname);
            md = loadmodel(filename);
            md = transientrestart(md);
            md = setflowequation(md,'SSA','all');

            % Load ensemble input from HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
            md.geometry.thickness = h5read(filename, '/Thickness');
            % md.geometry.surface   = h5read(filename, '/Surface');
            md.geometry.bed = h5read(filename, '/bed');
            md.friction.coefficient = h5read(filename, '/coefficient');

            % Update geometry
            % md.geometry.thickness = md.results.TransientSolution(end).Thickness;
            % md.geometry.surface   = md.results.TransientSolution(end).Surface;
            % md.geometry.base      = md.results.TransientSolution(end).Base;
            % md.geometry.base = md.geometry.surface  - md.geometry.thickness; % Update base from thickness and surface

            % Update other fields
            % md.initialization.vx        = md.results.TransientSolution(end).Vx;
            % md.initialization.vy        = md.results.TransientSolution(end).Vy;
            % md.initialization.vel       = md.results.TransientSolution(end).Vel;
            % % md.initialization.pressure  = md.results.TransientSolution(end).Pressure;
            % md.smb.mass_balance         = md.results.TransientSolution(end).SmbMassBalance;
            % md.mask.ocean_levelset      = md.results.TransientSolution(end).MaskOceanLevelset;

            % Calculate thickness
            % md.geometry.thickness = md.geometry.surface - md.geometry.base;

            % Set minimum thickness to 1
            pos = find(md.geometry.thickness < 1);
            md.geometry.thickness(pos) = 1;

            % Update surface
            md.geometry.surface = md.geometry.base + md.geometry.thickness;

            % disp('      -- ice shelf base based on hydrostatic equilibrium');

            % Calculate density ratio
            di = md.materials.rho_ice / md.materials.rho_water;

            % Calculate ocean levelset
            md.mask.ocean_levelset = md.geometry.thickness + md.geometry.bed / di;

            % Find positions where ocean_levelset < 0
            pos = find(md.mask.ocean_levelset < 0);

            % Update surface for floating ice
            md.geometry.surface(pos) = md.geometry.thickness(pos) * ...
                (md.materials.rho_water - md.materials.rho_ice) / md.materials.rho_water;

            % Update base
            md.geometry.base = md.geometry.surface - md.geometry.thickness;

            % Ensure base is not below bed
            pos = find(md.geometry.base < md.geometry.bed);
            md.geometry.base(pos) = md.geometry.bed(pos);

            % For grounded ice (ocean_levelset > 0)
            pos = find(md.mask.ocean_levelset > 0);
            md.geometry.base(pos) = md.geometry.bed(pos);

            % Final surface update
            md.geometry.surface = md.geometry.base + md.geometry.thickness;


            % Time stepping
            md.timestepping = timestepping();
            md.timestepping.time_step = dt;
            md.timestepping.start_time = tinitial;
            md.timestepping.final_time = tfinal;
            md.settings.output_frequency = 10;
            md.stressbalance.maxiter = 100;
            md.stressbalance.restol = 1;
            md.stressbalance.reltol = 0.001;
            md.stressbalance.abstol = NaN;

            % Cluster setup
            if hpcmode
                md.settings.waitonlock = 0;
                md.cluster = generic('name', oshostname(), 'np', nprocs);
                md.cluster.codepath = [issmroot , '/bin'];
                md.cluster.login = 'arobel3';
                if ~exist(sprintf('%s/execution/color_%d', issmroot, ens_id), 'dir')
                    mkdir(sprintf('%s/execution/color_%d', issmroot, ens_id));
                end
                md.cluster.executionpath = sprintf('%s/execution/color_%d', issmroot, ens_id);
            else
                md.cluster = generic('name', cluster_name, 'np', nprocs);
                md.settings.waitonlock = 1;
                md.miscellaneous.name = sprintf('color_%d', ens_id);
            end

            % Verbose settings
            md.verbose = verbose('convergence', false, 'solution', true);

            % Solve transient
            md = solve(md, 'Transient');

            % Save model
            filename = fullfile(folder, data_fname);
            save(filename, 'md');

            % Save ensemble outputs in HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));

            fields = {'Thickness', 'bed', 'coefficient'};
            result_0 = md.results.TransientSolution(end);
            result_1 = md.geometry;
            result_2 = md.friction;

            data = {'Thickness', result_0, 'Thickness';
                    % 'Surface', result_0, 'Surface';
                    'bed', result_1, 'bed';
                    'coefficient', result_2, 'coefficient'};

            writeToHDF5(filename, data);
        % end
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
