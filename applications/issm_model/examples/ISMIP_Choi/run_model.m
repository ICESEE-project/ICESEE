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


    % get the current working directory
    cwd = pwd;
    [issmroot,~,~]=fileparts(fileparts(cwd));
    if devmode
        newpath=fullfile(issmroot,'/src/m/dev');
        addpath(newpath);
        devpath;
    end

    % fprintf('[MATLAB] Running model with ens_id: %d, rank: %d, nprocs: %d, filename: %s\n', ens_id, rank, nprocs, data_fname);

    folder = sprintf('./Models/ens_id_%d', ens_id);
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    if any(steps == 6)
      
        if k == 0 || isempty(k)
            % Initial run: load boundary conditions
            filename = fullfile(folder, reference_data);
            md = loadmodel(filename);

            % md = transientrestart(md);
            % md.results = [];
            md=setflowequation(md,'SSA','all');
            
             % Update geometry
             md.geometry.thickness =  md.results.TransientSolution(end).Thickness;
             md.geometry.surface   =  md.results.TransientSolution(end).Surface;
             md.geometry.base      =  md.results.TransientSolution(end).Base;
 
             % update other fields
             md.initialization.vel      = md.results.TransientSolution(end).Vel;
             md.initialization.vx        = md.results.TransientSolution(end).Vx;
             md.initialization.vy        = md.results.TransientSolution(end).Vy;
             md.initialization.pressure = md.results.TransientSolution(end).Pressure;
             md.smb.mass_balance        = md.results.TransientSolution(end).SmbMassBalance;
             md.mask.ocean_levelset     = md.results.TransientSolution(end).MaskOceanLevelset;


            md.smb.mass_balance = -0.3 * ones(md.mesh.numberofvertices, 1); % m/yr
            
            md.basalforcings = linearbasalforcings();
            md.basalforcings.deepwater_melting_rate = 200; % m/yr
            md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices, 1);

            md.transient.ismovingfront=0;   


            % change from adaptive time stepping to fixed time stepping (must be done after model update)
            md.timestepping = timestepping();
            md.timestepping.time_step  = dt;
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
                md.cluster = generic('name',oshostname(),'np',nprocs);
                md.cluster.codepath = [issmroot , '/bin'];
                md.cluster.login = 'arobel3';
                % md.cluster.executionpath = [issmroot, '/execution'];
                md.cluster.executionpath = sprintf('%s/execution/color_%d', issmroot, ens_id);
            else
                md.cluster=generic('name',cluster_name,'np',nprocs);
                md.settings.waitonlock = 1;
                md.miscellaneous.name =  sprintf('color_%d', ens_id);
                % md.cluster.executionpath = sprintf('%s/execution/color_%d', issmroot, ens_id);
            end

            % Verbose settings
            md.verbose = verbose('convergence', false, 'solution', true);
           
            md = solve(md, 'Transient');

            % Save model
            filename = fullfile(folder, data_fname);
            save(filename, 'md');

            % Save ensemble outputs in HDF5
            fields = {'Thickness', 'bed', 'coefficient'};
            result_0 = md.results.TransientSolution(end);
            result_1 = md.geometry;
            result_2 = md.friction;

            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));

            data = {'Thickness', result_0, 'Thickness';
                    'bed', result_1, 'bed';
                    'coefficient', result_2, 'coefficient'};

            writeToHDF5(filename, data);

        else

            % Subsequent time steps: load previous transient state
            filename = fullfile(folder, data_fname);
            md = loadmodel(filename);
            % md = transientrestart(md);
            md=setflowequation(md,'SSA','all');

            % Load ensemble input from HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
            md.geometry.thickness = h5read(filename, '/Thickness');
            md.geometry.bed = h5read(filename, '/bed');
            md.friction.coefficient = h5read(filename, '/coefficient');

            % Update geometry
            md.geometry.surface   =  md.results.TransientSolution(end).Surface;
            md.geometry.base      =  md.results.TransientSolution(end).Base;

            % update other fields
            md.initialization.vel      = md.results.TransientSolution(end).Vel;
            md.initialization.pressure = md.results.TransientSolution(end).Pressure;
            md.smb.mass_balance        = md.results.TransientSolution(end).SmbMassBalance;
            md.mask.ocean_levelset     = md.results.TransientSolution(end).MaskOceanLevelset;

            % Time stepping
            md.timestepping.time_step = dt;
            md.timestepping.start_time = tinitial;
            md.timestepping.final_time = tfinal;
            % md.transient.requested_outputs = {'Vx', 'Vy', 'Thickness', 'Surface', 'Bed', 'FrictionCoefficient'};

            % Cluster setup
            if hpcmode
                md.settings.waitonlock = 0;
                md.cluster = generic('name',oshostname(),'np',nprocs);
                md.cluster.codepath = [issmroot , '/bin'];
                md.cluster.login = 'arobel3';
                % md.cluster.executionpath = [issmroot, '/execution'];
                if ~exist(sprintf('%s/execution/color_%d', issmroot, ens_id), 'dir')
                    mkdir(sprintf('%s/execution/color_%d', issmroot, ens_id));
                end
                md.cluster.executionpath = sprintf('%s/execution/color_%d', issmroot, ens_id);
            else
                md.cluster=generic('name',cluster_name,'np',nprocs);
                md.settings.waitonlock = 1;
                md.miscellaneous.name =  sprintf('color_%d', ens_id);
                % md=loadresultsfromcluster(md);
                % md.cluster.executionpath = sprintf('%s/execution/color_%d', issmroot, ens_id);
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
                    'bed', result_1, 'bed';
                    'coefficient', result_2, 'coefficient'};

            writeToHDF5(filename, data);
            
            % save_ensemble_hdf5(filename, result_0, {'Vx', 'Vy', 'Thickness'});
            % save_ensemble_hdf5(filename, result_1, {'bed'});
            % save_ensemble_hdf5(filename, result_2, {'coefficient'});
        end
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


function save_ensemble_hdf5(filename, result, field_names)
    % Save model outputs to HDF5 file
    [filepath, ~, ~] = fileparts(filename);
    if ~exist(filepath, 'dir')
        mkdir(filepath);
    end
    if isfile(filename)
        delete(filename);
    end
    for i = 1:length(field_names)
        field = field_names{i};
        if isfield(result, field)
            data = result.(field);
            h5create(filename, ['/' field], size(data));
            h5write(filename, ['/' field], data);
        else
            warning('Field "%s" not found in result. Skipping.', field);
        end
    end
    fprintf('[HDF5] Saved: %s\n', filename);
end