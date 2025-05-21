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

    fprintf('[MATLAB] Running model with ens_id: %d, rank: %d, nprocs: %d, filename: %s\n', ens_id, rank, nprocs, data_fname);

    folder = sprintf('./Models/ens_id_%d', ens_id);
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    if any(steps == 8)
        % Transient simulation (equivalent to notebook's Experiment 1)
        if k == 0 || isempty(k)
            % Initial run: load boundary conditions
            filename = fullfile(folder, 'ISMIP.BoundaryCondition.mat');
            md = loadmodel(filename);

            % Set transient parameters
            md.timestepping.time_step = dt; % Monthly steps (e.g., 0.0833 years)
            md.timestepping.start_time = tinitial;
            md.timestepping.final_time = tfinal;
            md.transient.requested_outputs = {'Vx', 'Vy', 'Thickness', 'Surface', 'Bed', 'FrictionCoefficient'};
            md.transient.ismovingfront = 0;
            md.transient.isthermal = 0;
            md.transient.isstressbalance = 1;
            md.transient.ismasstransport = 1;
            md.transient.isgroundingline = 1;
            md.groundingline.migration = 'SubelementMigration';
            md.groundingline.friction_interpolation = 'SubelementFriction1';
            md.groundingline.melt_interpolation = 'NoMeltOnPartiallyFloating';

            % Set forcings (from notebook Cell 11)
            md.smb.mass_balance = -0.3 * ones(md.mesh.numberofvertices, 1); % m/yr
            md.basalforcings = linearbasalforcings();
            md.basalforcings.deepwater_melting_rate = 200; % m/yr
            md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices, 1);

            % Cluster setup
            md.cluster = discover('numnodes', 1, 'cpuspernode', 40, 'time', 60*60*2, 'processor', 'sky', ...
                                 'queue', 'allnccs', 'port', 0, ...
                                 'executionpath', '/discover/nobackup/ychoi14/execution', ...
                                 'srcpath', '/discover/nobackup/ychoi14/trunk-python', ...
                                 'codepath', '/discover/nobackup/ychoi14/trunk-python/bin', ...
                                 'login', 'ychoi14', 'grouplist', 's2772');
            md.cluster.modules = {'comp/intel/2021.4.0', 'mpi/impi/2021.4.0'};
            md.settings.waitonlock = 0;

            % Verbose settings
            md.verbose = verbose('convergence', false, 'solution', true);

            % Set flow equation
            md = setflowequation(md, 'SSA', 'all');

            % Solve transient
            md.miscellaneous.name = sprintf('Transient_ENS_%d', ens_id);
            md = solve(md, 'Transient', 'runtimename', false);

            % Save model
            filename = fullfile(folder, data_fname);
            save(filename, 'md');

            % Save ensemble outputs in HDF5
            fields = {'Vx', 'Vy', 'Thickness', 'Surface', 'Bed', 'FrictionCoefficient'};
            result = md.results.TransientSolution(end);
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
            save_ensemble_hdf5(filename, result, fields);
        else
            % Subsequent time steps: load previous transient state
            filename = fullfile(folder, data_fname);
            md = loadmodel(filename);

            % Load ensemble input from HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
            md.initialization.vx = h5read(filename, '/Vx');
            md.initialization.vy = h5read(filename, '/Vy');
            md.initialization.thickness = h5read(filename, '/Thickness');
            md.geometry.surface = h5read(filename, '/Surface');
            md.geometry.bed = h5read(filename, '/Bed');
            md.friction.coefficient = h5read(filename, '/FrictionCoefficient');

            % Update geometry
            md.geometry.base = md.geometry.surface - md.geometry.thickness;
            pos = find(md.geometry.thickness < 1);
            md.geometry.thickness(pos) = 1;
            di = md.materials.rho_ice / md.materials.rho_water;
            md.mask.ocean_levelset = md.geometry.thickness + md.geometry.bed / di;
            pos = find(md.mask.ocean_levelset < 0);
            md.geometry.surface(pos) = md.geometry.thickness(pos) * (md.materials.rho_water - md.materials.rho_ice) / md.materials.rho_water;
            md.geometry.base = md.geometry.surface - md.geometry.thickness;
            pos = find(md.geometry.base < md.geometry.bed);
            md.geometry.base(pos) = md.geometry.bed(pos);
            pos = find(md.mask.ocean_levelset > 0);
            md.geometry.base(pos) = md.geometry.bed(pos);
            md.geometry.surface = md.geometry.base + md.geometry.thickness;

            % Time stepping
            md.timestepping.time_step = dt;
            md.timestepping.start_time = tinitial;
            md.timestepping.final_time = tfinal;
            md.transient.requested_outputs = {'Vx', 'Vy', 'Thickness', 'Surface', 'Bed', 'FrictionCoefficient'};

            % Cluster setup
            md.cluster = discover('numnodes', 1, 'cpuspernode', 40, 'time', 60*60*2, 'processor', 'sky', ...
                                 'queue', 'allnccs', 'port', 0, ...
                                 'executionpath', '/discover/nobackup/ychoi14/execution', ...
                                 'srcpath', '/discover/nobackup/ychoi14/trunk-python', ...
                                 'codepath', '/discover/nobackup/ychoi14/trunk-python/bin', ...
                                 'login', 'ychoi14', 'grouplist', 's2772');
            md.cluster.modules = {'comp/intel/2021.4.0', 'mpi/impi/2021.4.0'};
            md.settings.waitonlock = 0;

            % Verbose settings
            md.verbose = verbose('convergence', false, 'solution', true);

            % Solve transient
            md.miscellaneous.name = sprintf('Transient_ENS_%d', ens_id);
            md = solve(md, 'Transient', 'runtimename', false);

            % Save model
            filename = fullfile(folder, data_fname);
            save(filename, 'md');

            % Save ensemble outputs in HDF5
            fields = {'Vx', 'Vy', 'Thickness', 'Surface', 'Bed', 'FrictionCoefficient'};
            result = md.results.TransientSolution(end);
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
            save_ensemble_hdf5(filename, result, fields);
        end
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