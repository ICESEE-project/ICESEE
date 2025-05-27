function run_model(data_fname, ens_id, rank, nprocs, k, dt, tinitial, tfinal)
    % Run ISSM model for transient simulation with ensemble support
    % Inputs: data_fname (output file name), ens_id (ensemble ID), rank, nprocs (MPI settings),
    %         k (time step index), dt (time step), tinitial, tfinal (time bounds)

    % Read kwargs from .mat file
    model_kwargs = sprintf('model_kwargs_%d.mat', ens_id);
    kwargs = load(model_kwargs);
    cluster_name = char(kwargs.cluster_name);
    steps = double(kwargs.steps);
    icesee_path = char(kwargs.icesee_path);
    data_path = char(kwargs.data_path);

    % fprintf('[MATLAB] Running model with ens_id: %d, rank: %d, nprocs: %d, filename: %s\n', ens_id, rank, nprocs, data_fname);

    folder = sprintf('./Models/ens_id_%d', ens_id);
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    if any(steps == 6)
      
        if k == 0 || isempty(k)
            % Initial run: load boundary conditions
            filename = fullfile(folder, 'ISMIP.BoundaryCondition.mat');
            md = loadmodel(filename);

            % Set transient parameters
            md.timestepping.time_step  = dt;
            md.timestepping.start_time = tinitial;
            md.timestepping.final_time = tfinal;

            % md.transient.requested_outputs = {'Vx', 'Vy', 'Thickness', 'Surface', 'Bed', 'FrictionCoefficient'};
            % md.transient.ismovingfront = 0;
            % md.transient.isthermal = 0;
            % md.transient.isstressbalance = 1;
            % md.transient.ismasstransport = 1;
            % md.transient.isgroundingline = 1;
            % md.groundingline.migration = 'SubelementMigration';
            % md.groundingline.friction_interpolation = 'SubelementFriction1';
            % md.groundingline.melt_interpolation = 'NoMeltOnPartiallyFloating';

            % % Set forcings (from notebook Cell 11)
            % md.smb.mass_balance = -0.3 * ones(md.mesh.numberofvertices, 1); % m/yr
            % md.basalforcings = linearbasalforcings();
            % md.basalforcings.deepwater_melting_rate = 200; % m/yr
            % md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices, 1);

            md.cluster=generic('name',cluster_name,'np',nprocs);

            % Verbose settings
            md.verbose = verbose('convergence', false, 'solution', true);
            md.transient.ismovingfront = 0;
            md.transient.isthermal = 0;
            md.miscellaneous.name =  sprintf('color_%d', ens_id);

            md = solve(md, 'Transient');

            % Save model
            filename = fullfile(folder, data_fname);
            save(filename, 'md');

            % Save ensemble outputs in HDF5
            fields = {'Thickness', 'Vx', 'Vy', 'bed', 'coefficient'};
            result_0 = md.results.TransientSolution(end);
            result_1 = md.geometry;
            result_2 = md.friction;

            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));

            data = {'Vx', result_0, 'Vx';
                    'Vy', result_0, 'Vy';
                    'Thickness', result_0, 'Thickness';
                    'bed', result_1, 'bed';
                    'coefficient', result_2, 'coefficient'};

            writeToHDF5(filename, data);

        else

            % Subsequent time steps: load previous transient state
            filename = fullfile(folder, data_fname);
            md = loadmodel(filename);

            % Load ensemble input from HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
            md.initialization.vx = h5read(filename, '/Vx');
            md.initialization.vy = h5read(filename, '/Vy');
            md.geometry.thickness = h5read(filename, '/Thickness');
            md.geometry.bed = h5read(filename, '/bed');
            md.friction.coefficient = h5read(filename, '/coefficient');

            % Update geometry
            md.geometry.thickness =  md.results.TransientSolution(end).Thickness;
            md.geometry.surface   =  md.results.TransientSolution(end).Surface;
            md.geometry.base      =  md.results.TransientSolution(end).Base;

            % update other fields
            md.initialization.vel      = md.results.TransientSolution(end).Vel;
            md.initialization.pressure = md.results.TransientSolution(end).Pressure;
            md.smb.mass_balance        = md.results.TransientSolution(end).SmbMassBalance;
            md.mask.ocean_levelset     = md.results.TransientSolution(end).MaskOceanLevelset;


            % pos = find(md.geometry.thickness < 1);
            % md.geometry.thickness(pos) = 1;
            % di = md.materials.rho_ice / md.materials.rho_water;
            % md.mask.ocean_levelset = md.geometry.thickness + md.geometry.bed / di;
            % pos = find(md.mask.ocean_levelset < 0);
            % md.geometry.surface(pos) = md.geometry.thickness(pos) * (md.materials.rho_water - md.materials.rho_ice) / md.materials.rho_water;
            % md.geometry.base = md.geometry.surface - md.geometry.thickness;
            % pos = find(md.geometry.base < md.geometry.bed);
            % md.geometry.base(pos) = md.geometry.bed(pos);
            % pos = find(md.mask.ocean_levelset > 0);
            % md.geometry.base(pos) = md.geometry.bed(pos);
            % md.geometry.surface = md.geometry.base + md.geometry.thickness;

            % Time stepping
            md.timestepping.time_step = dt;
            md.timestepping.start_time = tinitial;
            md.timestepping.final_time = tfinal;
            % md.transient.requested_outputs = {'Vx', 'Vy', 'Thickness', 'Surface', 'Bed', 'FrictionCoefficient'};

            % Cluster setup
            md.cluster = generic('name', cluster_name, 'np', nprocs);

            % Verbose settings
            md.verbose = verbose('convergence', false, 'solution', true);

            md.transient.ismovingfront = 0;
            md.transient.isthermal = 0;

            % Solve transient
            md.miscellaneous.name =  sprintf('color_%d', ens_id);
            md = solve(md, 'Transient');

            % Save model
            filename = fullfile(folder, data_fname);
            save(filename, 'md');

            % Save ensemble outputs in HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));

            fields = {'Thickness', 'Vx', 'Vy', 'bed', 'coefficient'};
            result_0 = md.results.TransientSolution(end);
            result_1 = md.geometry;
            result_2 = md.friction;

            data = {'Vx', result_0, 'Vx';
                    'Vy', result_0, 'Vy';
                    'Thickness', result_0, 'Thickness';
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