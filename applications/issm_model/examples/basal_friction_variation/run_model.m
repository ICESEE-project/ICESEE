
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
    devmode      = logical(kwargs.devmode);
    issm_example_dir     = char(kwargs.issm_examples_dir);
    deepwater_melting_rate = double(kwargs.deepwater_melting_rate);
    smb = double(kwargs.smb);
    mean_friction  = double(kwargs.mean_friction);
    reference_data = char(kwargs.reference_data);
    nens = double(kwargs.Nens);
    wrong_reference_data = 'wrong_reference_data.mat';
    min_friction = double(kwargs.min_friction);
    max_friction = double(kwargs.max_friction);
    abs_vel_weight = double(kwargs.abs_vel_weight);
    rel_vel_weight = double(kwargs.rel_vel_weight);
    tikhonov_regularization_weight = double(kwargs.tikhonov_regularization_weight);

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

    % Get the current working directory
    cwd = pwd;
    [issmroot,~,~] = fileparts(fileparts(cwd));

    % number of variables
    nvar = 8;
    rng(1000 + ens_id, 'twister');   % ens_id = ensemble index

    % set initail ens_id
    ens_id_init = 0;
    s_perturb = double(kwargs.s_nurge);
    b_perturb = double(kwargs.b_nurge);

    output_frequency = 1; % make sure this is set to 1 for coupling with ICESEE

    % Set up model for each EnKF stage
    if strcmp(data_fname, 'true_state.mat')
        % Special case for true state
        folder = sprintf('./Models/ens_id_%d', ens_id_init);
        if ~exist(folder, 'dir')
            mkdir(folder);
        end

        % Initial run: load boundary conditions
        filename = fullfile(folder, reference_data);
        md = loadmodel(filename);

        md = setflowequation(md,'SSA','all');

        md.smb.mass_balance = smb * ones(md.mesh.numberofvertices,1);
        md.transient.ismovingfront = 0;

        md.basalforcings = linearbasalforcings();
        md.basalforcings.deepwater_melting_rate = deepwater_melting_rate;
        md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);

        % -- time stepping
        md.timestepping = timestepping();
        md.timestepping.time_step = dt;
        md.settings.output_frequency = output_frequency;   % make sure this is 1
        md.stressbalance.maxiter = 100;
        md.stressbalance.restol = 1;
        md.stressbalance.reltol = 0.001;
        md.stressbalance.abstol = NaN;
        md.settings.solver_residue_threshold = 5e-2;

        % Cluster setup
        md.cluster = generic('name', oshostname(), 'np', nprocs);
        md.settings.waitonlock = 1;
        md.miscellaneous.name = sprintf('color_%d', ens_id);

        % Verbose settings
        md.verbose = verbose('convergence', false, 'solution', true);

        md.transient.requested_outputs = {'default', 'FrictionCoefficient', ...
                                        'Thickness', 'Surface', 'Base', 'Bed', ...
                                        'IceVolume', 'IceVolumeAboveFloatation'};

        % ------------------------------------------------------------------
        % Synthetic truth perturbation settings
        % ------------------------------------------------------------------
        % Example: 4 friction changes during observed period
        % change_times = [6.0, 12.0, 18.0, 22.0];   % years
        % amps         = [0.08, 0.15, 0.22, 0.12];  % percentage relative increases
        % sigma        = 3e3;                       % 3 km Gaussian width
        % change_times = [6.0, 12.0, 18.0, 22.0];
        % amps         = [0.25, 0.50, 0.75, 0.35];
        % sigma        = 20e3;   % 20 km
        change_times = [2.0];   % years
        amps         = [0.50];
        sigma        = 20e3;   % 20 km

        % Background friction coefficient
        C0 = md.friction.coefficient;

        % Gaussian patch centered in the middle of the mesh
        x = md.mesh.x;
        y = md.mesh.y;
        x0 = (0.5 * (min(x) + max(x)))./2.0  % middle of the icesheet 
        y0 = 0.5 * (min(y) + max(y));
        gauss = exp(-((x - x0).^2 + (y - y0).^2) ./ (2 * sigma^2));

        % Time grid
        time_vec = tinitial:dt:tfinal;
        if abs(time_vec(end) - tfinal) > 1e-12
            time_vec = [time_vec, tfinal];
        end
        N = length(time_vec) - 1;

        % Preallocate results container
        transient_results = cell(N,1);

        % ------------------------------------------------------------------
        % March forward one interval at a time
        % ------------------------------------------------------------------
        for it = 1:N
            t_now  = time_vec(it);
            t_next = time_vec(it+1);

            % Determine which friction perturbation is active
            active_idx = find(t_now >= change_times, 1, 'last');

            if isempty(active_idx)
                md.friction.coefficient = C0;
                fprintf('Time [%g, %g]: using baseline friction\n', t_now, t_next);
            else
                md.friction.coefficient = C0 .* (1 + amps(active_idx) * gauss);
                fprintf('Time [%g, %g]: using friction regime %d (amp = %.3f)\n', ...
                 t_now, t_next, active_idx, amps(active_idx));
            end

            % Solve only over the current interval
            md.timestepping.start_time = t_now;
            md.timestepping.final_time = t_next;
            md.timestepping.time_step  = t_next - t_now;

            md = solve(md, 'Transient', 'runtimename', false);

            % Save final state of this interval
            sol = md.results.TransientSolution(end);
            transient_results{it} = sol;

            % Update model state for next interval
            md.geometry.thickness = sol.Thickness;
            md.geometry.surface   = sol.Surface;
            md.geometry.base      = sol.Base;

            md.initialization.vx       = sol.Vx;
            md.initialization.vy       = sol.Vy;
            md.initialization.vel      = sol.Vel;
            md.initialization.pressure = sol.Pressure;

            if isfield(sol, 'SmbMassBalance')
                md.smb.mass_balance = sol.SmbMassBalance;
            end
            if isfield(sol, 'MaskOceanLevelset')
                md.mask.ocean_levelset = sol.MaskOceanLevelset;
            end

            % Clear results before next short run
            md.results = struct();
        end

        % Save updated final model
        filename = fullfile(folder, data_fname);
        save(filename, 'md', '-v7.3');

        % ------------------------------------------------------------------
        % Pack outputs for HDF5 export
        % ------------------------------------------------------------------
        scalar_inputs = {'IceVolume', 'IceVolumeAboveFloatation'};

        data = cell(N * length(vec_inputs), 3);
        data_scalar = cell(N * length(scalar_inputs), 3);

        idx = 1;
        idxs = 1;

        for k = 1:N
            sol = transient_results{k};

            % mesh / state-like outputs
            for j = 1:length(vec_inputs)
                key = vec_inputs{j};

                data{idx, 1} = sprintf('%s_%d', key, k);
                data{idx, 2} = sol;
                data{idx, 3} = key;
                idx = idx + 1;
            end

            % scalar diagnostics at each time step
            for j = 1:length(scalar_inputs)
                key = scalar_inputs{j};

                data_scalar{idxs, 1} = sprintf('%s_%d', key, k);
                data_scalar{idxs, 2} = sol;
                data_scalar{idxs, 3} = key;
                idxs = idxs + 1;
            end
        end

        data = data(1:idx-1, :);
        data_scalar = data_scalar(1:idxs-1, :);

        filename = fullfile(icesee_path, data_path, sprintf('ensemble_true_state_%d.h5', ens_id));
        writeToHDF5(filename, data);

        filename = fullfile(icesee_path, data_path, sprintf('ensemble_true_state_scalar_%d.h5', ens_id));
        writeToHDF5(filename, data_scalar);

    elseif strcmp(data_fname, 'nurged_state.mat')

        folder = sprintf('./Models/ens_id_%d', ens_id_init);
        if ~exist(folder, 'dir')
            mkdir(folder);
        end
            
        filename = fullfile(folder, reference_data);
        md = loadmodel(filename);

        md = setflowequation(md,'SSA','all');

        % setup nugged state
        friction_ref = mean_friction*ones(md.mesh.numberofvertices,1);
        % friction_ref = md.friction.coefficient;
        % friction_ref = md_ref.friction.coefficient;
        thickness_ref = md.geometry.thickness;
        bed_ref = md.geometry.bed;
        base_ref = md.geometry.base;

        % % read the friction_bed file
        filename = fullfile(icesee_path, data_path, sprintf('friction_bed_%d.h5', ens_id));
        % bed = h5read(filename, '/bed');
        coefficient = h5read(filename, '/FrictionCoefficient');


        md.friction.coefficient = friction_ref + coefficient;
        % md.friction.coefficient = friction_ref;
        md.friction.p=ones(md.mesh.numberofelements,1);
        md.friction.q=ones(md.mesh.numberofelements,1);

        % bed_err = bed - bed_ref;
        % md.geometry.bed = (bed_ref + bed_err) - b_perturb*randn(md.mesh.numberofvertices, 1);
        % md.geometry.base = (base_ref + bed_err) - b_perturb*randn(md.mesh.numberofvertices, 1);
        % md.geometry.surface = (md.geometry.surface + bed_err) - s_perturb*randn(md.mesh.numberofvertices, 1);

        % md.geometry.thickness = md.geometry.surface - md.geometry.base; %- 50*ones(md.mesh.numberofvertices,1);
        
        % Ensure minimum ice thickness of 1 m
        pos = find(md.geometry.thickness < 1);
        md.geometry.thickness(pos) = 1;
        md.geometry.surface = md.geometry.base + md.geometry.thickness;

        % md.geometry.surface = md.geometry.surface + s_perturb*ones(md.mesh.numberofvertices,1);

        disp('      -- ice shelf base based on hydrostatic equilibrium');
        di = md.materials.rho_ice / md.materials.rho_water;

        % Compute ocean level set based on hydrostatic equilibrium
        md.mask.ocean_levelset = md.geometry.thickness + md.geometry.bed / di;

        % Floating ice (ocean_levelset < 0)
        pos = find(md.mask.ocean_levelset < 0);
        md.geometry.surface(pos) = md.geometry.thickness(pos) .* ...
        (md.materials.rho_water - md.materials.rho_ice) / md.materials.rho_water;
        md.geometry.base = md.geometry.surface - md.geometry.thickness;

        % Ensure base not below bedrock
        pos = find(md.geometry.base < md.geometry.bed);
        md.geometry.base(pos) = md.geometry.bed(pos);
        % md.geometry.base(pos) = md.geometry.base(pos);

        % Grounded ice (ocean_levelset > 0)
        pos = find(md.mask.ocean_levelset > 0);
        md.geometry.base(pos) = md.geometry.bed(pos);
        md.geometry.surface = md.geometry.base + md.geometry.thickness;

        md.smb.mass_balance=smb*ones(md.mesh.numberofvertices,1);
        md.transient.ismovingfront=0;
        % 
        md.initialization.pressure       = zeros(md.mesh.numberofvertices, 1);
        md.masstransport.spcthickness    = NaN * ones(md.mesh.numberofvertices, 1);
        md.basalforcings=linearbasalforcings();
        md.basalforcings.deepwater_melting_rate=deepwater_melting_rate;
        md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);

       % --time stepping
       md.timestepping = timestepping();
       md.timestepping.time_step = dt;
       md.timestepping.start_time = tinitial;
       md.timestepping.final_time = tfinal;
       md.settings.output_frequency = output_frequency; %make sure this is set to 1 for 
       md.stressbalance.maxiter = 100;
       md.stressbalance.restol = 1;
       md.stressbalance.reltol = 0.001;
       md.stressbalance.abstol = NaN;
       md.settings.solver_residue_threshold=5e-2;

        % Cluster setup
        md.cluster = generic('name', oshostname(), 'np', nprocs);
        md.settings.waitonlock = Inf;
        md.settings.waitonlock=1;
        md.miscellaneous.name = sprintf('color_%d', ens_id);

        % Verbose settings
        md.verbose = verbose('convergence', false, 'solution', true);

        md.transient.requested_outputs = {'default', 'FrictionCoefficient', 'Thickness', ...
                                            'Surface','Base','Bed', 'IceVolume', 'IceVolumeAboveFloatation'};

        % Solve transient
        md = solve(md, 'Transient','runtimename',false);
            
        filename = fullfile(folder, data_fname);
        save(filename, 'md', '-v7.3');

        N = length(md.results.TransientSolution);

        scalar_inputs = {'IceVolume', 'IceVolumeAboveFloatation'};

        data = cell(N * length(vec_inputs), 3);
        data_scalar = cell(N * length(scalar_inputs), 3);

        idx = 1;
        idxs = 1;

        for k = 1:N
            sol = md.results.TransientSolution(k);

            for j = 1:length(vec_inputs)
                key = vec_inputs{j};
                data{idx, 1} = sprintf('%s_%d', key, k);
                data{idx, 2} = sol;
                data{idx, 3} = key;
                idx = idx + 1;
            end

            for j = 1:length(scalar_inputs)
                key = scalar_inputs{j};
                data_scalar{idxs, 1} = sprintf('%s_%d', key, k);
                data_scalar{idxs, 2} = sol;
                data_scalar{idxs, 3} = key;
                idxs = idxs + 1;
            end
        end

        data = data(1:idx-1, :);
        data_scalar = data_scalar(1:idxs-1, :);

        writeToHDF5(fullfile(icesee_path, data_path, sprintf('ensemble_nurged_state_%d.h5', ens_id)), data);
        writeToHDF5(fullfile(icesee_path, data_path, sprintf('ensemble_nurged_state_scalar_%d.h5', ens_id)), data_scalar);

    elseif strcmp(data_fname, 'initialize_ensemble.mat')
        % Special case for ensemble initialization
        if k+1 == 0 || isempty(k)
            % Initial run: load boundary conditions
            % filename = fullfile(folder, reference_data);
            folder = sprintf('./Models/ens_id_%d', ens_id_init);
            % folder = sprintf('./Models/ens_id_%d', ens_id);
            if ~exist(folder, 'dir')
                mkdir(folder);
            end

            % seed the random number generator for reproducibility
            % rng(ens_id + 1000); % Offset seed to avoid overlap with other uses

            filename = fullfile(folder, reference_data);
            % filename = fullfile(icesee_path, 'data', wrong_reference_data);
            md = loadmodel(filename);
            md = setflowequation(md,'SSA','all');

            friction_ref = mean_friction*ones(md.mesh.numberofvertices,1);

            thickness_ref = md.geometry.thickness;
            bed_ref = md.geometry.bed;
            base_ref = md.geometry.base;

             % % read the friction_bed file
            filename = fullfile(icesee_path, data_path, sprintf('friction_bed_%d.h5', ens_id));
            % bed = h5read(filename, '/bed');
            coefficient = h5read(filename, '/FrictionCoefficient');

            %  update the friction and bed
            md.friction.coefficient = friction_ref + coefficient;
            % md.friction.coefficient = friction_ref;
            md.friction.p=ones(md.mesh.numberofelements,1);
            md.friction.q=ones(md.mesh.numberofelements,1);

 
            % bed_err = bed - bed_ref;
            % md.geometry.bed = (bed_ref + bed_err) - b_perturb*randn(md.mesh.numberofvertices, 1);
            % md.geometry.base = (base_ref + bed_err) - b_perturb*randn(md.mesh.numberofvertices, 1);
            % md.geometry.surface = (md.geometry.surface + bed_err) - s_perturb*randn(md.mesh.numberofvertices, 1);

            md.initialization.pressure       = zeros(md.mesh.numberofvertices, 1);
            md.masstransport.spcthickness    = NaN * ones(md.mesh.numberofvertices, 1);
            md.transient.ismovingfront=0;

            md.smb.mass_balance=smb*ones(md.mesh.numberofvertices,1);
            md.basalforcings=linearbasalforcings();
            md.basalforcings.deepwater_melting_rate=deepwater_melting_rate;
            md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);

            md.geometry.thickness = md.geometry.surface - md.geometry.base;

            % Ensure minimum ice thickness of 1 m
            pos = find(md.geometry.thickness < 1);
            md.geometry.thickness(pos) = 1;
            % md.geometry.thickness(pos) = max(1, min(thickness_ref));
            md.geometry.surface = md.geometry.base + md.geometry.thickness;
            % md.geometry.surface = md.geometry.surface + s_perturb*ones(md.mesh.numberofvertices,1);

            disp('      -- ice shelf base based on hydrostatic equilibrium');
            di = md.materials.rho_ice / md.materials.rho_water;

            % Compute ocean level set based on hydrostatic equilibrium
            md.mask.ocean_levelset = md.geometry.thickness + md.geometry.bed / di;

            % Floating ice (ocean_levelset < 0)
            pos = find(md.mask.ocean_levelset < 0);
            md.geometry.surface(pos) = md.geometry.thickness(pos) .* ...
                (md.materials.rho_water - md.materials.rho_ice) / md.materials.rho_water;
            md.geometry.base = md.geometry.surface - md.geometry.thickness;

            % Ensure base not below bedrock
            pos = find(md.geometry.base < md.geometry.bed);
            md.geometry.base(pos) = md.geometry.bed(pos);

            % Grounded ice (ocean_levelset > 0)
            pos = find(md.mask.ocean_levelset > 0);
            md.geometry.base(pos) = md.geometry.bed(pos);
            md.geometry.surface = md.geometry.base + md.geometry.thickness;

            % pos = find(md.mask.ocean_levelset < 0);
            % md.geometry.thickness(pos)=1/(1-di)*md.geometry.surface(pos);


            % --time stepping
            md.timestepping = timestepping();
            md.timestepping.time_step = 0.2;
            md.timestepping.start_time = 0;
            md.timestepping.final_time = 1.0;
            md.settings.output_frequency = output_frequency; %make sure this is set to 1 for 
            md.stressbalance.maxiter = 100;
            md.stressbalance.restol = 1;
            md.stressbalance.reltol = 0.001;
            md.stressbalance.abstol = NaN;
            md.settings.solver_residue_threshold=5e-2;

            % Cluster setup
            md.cluster = generic('name', oshostname(), 'np', nprocs);
            md.settings.waitonlock = Inf;
            md.settings.waitonlock=1;
            md.miscellaneous.name = sprintf('color_%d', ens_id);

            % Verbose settings
            md.verbose = verbose('convergence', false, 'solution', true);

            md.transient.requested_outputs = {'default', 'FrictionCoefficient', 'Thickness', ...
                                            'Surface','Base','Bed', 'IceVolume', 'IceVolumeAboveFloatation'};

            % Solve transient
            md = solve(md, 'Transient','runtimename',false);
             
            % save updated model to every ensemble folder
            folder = sprintf('./Models/ens_id_%d', ens_id);
            if ~exist(folder, 'dir')
                mkdir(folder);
            end
            filename = fullfile(folder, data_fname);
            save(filename, 'md', '-v7.3');

            % Save ensemble outputs in HDF5
            result_0 = md.results.TransientSolution(end);

            filename = fullfile(icesee_path, data_path, sprintf('ensemble_out_%d.h5', ens_id));

            data = cell(length(vec_inputs), 3);
            for j = 1:length(vec_inputs)
                key = vec_inputs{j};
                data{j, 1} = key;
                data{j, 2} = result_0;
                data{j, 3} = key;
            end

            writeToHDF5(filename, data);

            filename = fullfile(icesee_path, data_path, sprintf('ensemble_out_scalar_%d.h5', ens_id));
            h5write(filename, '/IceVolume', result_0.IceVolume, [k+1], [1]);
            h5write(filename, '/IceVolumeAboveFloatation', result_0.IceVolumeAboveFloatation, [k+1], [1]);

            % Break and return to avoid further processing
            return;
        end

    elseif strcmp(data_fname, 'enkf_state.mat')
        % Special case for ensemble assimilation
        folder = sprintf('./Models/ens_id_%d', ens_id);
        if ~exist(folder, 'dir')
            mkdir(folder);
        end
        
        if k == 0 || isempty(k)
            % Initial run: load boundary conditions
            % filename = fullfile(folder, reference_data);
            % filename_ens_init = fullfile(folder, 'initialize_ensemble.mat');

            folder_init = sprintf('./Models/ens_id_%d', ens_id_init);
            % folder = sprintf('./Models/ens_id_%d', ens_id);
            if ~exist(folder_init, 'dir')
                mkdir(folder_init);
            end
            
            folder_true = sprintf('./Models/ens_id_%d', 0);
            if ~exist(folder_true, 'dir')
                mkdir(folder_true);
            end
            filename = fullfile(folder_true, 'true_state.mat');
            % filename = fullfile(folder_init, reference_data);
            % filename = fullfile(icesee_path, 'data', wrong_reference_data);

            % filename = fullfile(folder, 'initialize_ensemble.mat');
            md = loadmodel(filename);
            
            md.inversion.iscontrol            = 0;
            md.transient.ismovingfront        = 0;
            md.transient.isthermal            = 0;
            md.transient.isstressbalance      = 1;
            md.transient.ismasstransport      = 1;
            md.transient.isgroundingline      = 1;

            md.groundingline.migration                = 'SubelementMigration';
            md.groundingline.friction_interpolation   = 'SubelementFriction1';
            md.groundingline.melt_interpolation       = 'NoMeltOnPartiallyFloating';

            md.initialization.pressure       = zeros(md.mesh.numberofvertices, 1);
            md.masstransport.spcthickness    = NaN * ones(md.mesh.numberofvertices, 1);

            md.verbose.solution              = 1;

            mask_all = zeros(md.mesh.numberofvertices,1);
            md.smb.mass_balance=smb*ones(md.mesh.numberofvertices,1);
            md.basalforcings=linearbasalforcings();
            md.basalforcings.deepwater_melting_rate=deepwater_melting_rate;
            md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);


             % Load ensemble input from HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
            md.geometry.surface = h5read(filename, '/Surface');
            % md.geometry.base = h5read(filename, '/Base');
            md.geometry.thickness = h5read(filename, '/Thickness');
            md.initialization.vx = h5read(filename, '/Vx');
            md.initialization.vy = h5read(filename, '/Vy');
            md.initialization.vel = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);
        
            % parameters for bed and friction
            % md.geometry.bed = h5read(filename, '/bed');
            md.friction.coefficient = h5read(filename, '/FrictionCoefficient');

            % --time stepping
            md.timestepping = timestepping();
            md.timestepping.time_step = 0.2;
            md.timestepping.start_time = 0;
            md.timestepping.final_time = 0.2;
            md.settings.output_frequency = output_frequency; %make sure this is set to 1 for
            
            % Ensure minimum ice thickness
            pos = find(md.geometry.thickness < 1);
            md.geometry.thickness(pos) = 1;

            % Compute density ratio
            di = md.materials.rho_ice / md.materials.rho_water;

            % Compute ocean level set
            md.mask.ocean_levelset = md.geometry.thickness + md.geometry.bed / di;

            % Floating ice (ocean_levelset < 0)
            pos = find(md.mask.ocean_levelset < 0);
            md.geometry.surface(pos) = md.geometry.thickness(pos) .* ...
                (md.materials.rho_water - md.materials.rho_ice) / md.materials.rho_water;

            % Update base geometry
            md.geometry.base = md.geometry.surface - md.geometry.thickness;

            % Ensure base is not below bedrock
            pos = find(md.geometry.base < md.geometry.bed);
            md.geometry.base(pos) = md.geometry.base(pos);
            %md.geometry.base(pos) = md.geometry.bed(pos);

            % Grounded ice (ocean_levelset > 0)
            pos = find(md.mask.ocean_levelset > 0);
            md.geometry.base(pos) = md.geometry.bed(pos);

            % Update surface geometry
            md.geometry.surface = md.geometry.base + md.geometry.thickness;

            % % Outputs and verbosity
            md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed','IceVolume','IceVolumeAboveFloatation'};
            % md.transient.requested_outputs = {'default','Thickness','Surface','Base','Bed'};
            md.verbose = verbose('all', false);
            md.verbose.solution = true;

            % Cluster setup
            md.cluster = generic('name', oshostname(), 'np', nprocs);
            md.settings.waitonlock = Inf;
            md.settings.waitonlock=1;
            md.miscellaneous.name = sprintf('color_%d', ens_id);

            % % Verbose settings
            md.verbose = verbose('convergence', false, 'solution', true);

            % % Solve transient
            md = solve(md, 'Transient','runtimename',false); %TODO: instead of solving just take th initial solution

            % Save model
            filename = fullfile(folder, data_fname);
            save(filename, 'md', '-v7.3');

            % Save ensemble outputs in HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
            result_0 = md.results.TransientSolution(end);
            % result_1 = md.geometry;
            result_1 = md.results.TransientSolution(end);
            % result_2 = md.friction;
            result_2 = md.results.TransientSolution(end);

            data = cell(length(vec_inputs), 3);
            for j = 1:length(vec_inputs)
                key = vec_inputs{j};
                data{j, 1} = key;
                data{j, 2} = result_0;
                data{j, 3} = key;
            end
                    
            writeToHDF5(filename, data);

            filename = fullfile(icesee_path, data_path, sprintf('ensemble_out_scalar_%d.h5', ens_id));
            h5write(filename, '/IceVolume', result_0.IceVolume, [k+1], [1]);
            h5write(filename, '/IceVolumeAboveFloatation', result_0.IceVolumeAboveFloatation, [k+1], [1]);

        else
          
            % fprintf('[MATLAB ---] Running model for ensemble ID %d, step %d\n', ens_id, k);
            
            % Subsequent time steps: 
            filename = fullfile(folder, data_fname);
            md = loadmodel(filename);

            % Load ensemble input from HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
            md.geometry.surface = h5read(filename, '/Surface');
            % md.geometry.base = h5read(filename, '/Base');
            md.geometry.thickness = h5read(filename, '/Thickness');
            md.initialization.vx = h5read(filename, '/Vx');
            md.initialization.vy = h5read(filename, '/Vy');
            md.initialization.vel = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);
            md.initialization.pressure=md.materials.rho_ice*md.constants.g*h5read(filename, '/Thickness');
        
            % parameters for bed and friction
            % md.geometry.bed = h5read(filename, '/bed');
            md.friction.coefficient = h5read(filename, '/FrictionCoefficient');

            % Ensure minimum ice thickness
            pos = find(md.geometry.thickness < 1);
            md.geometry.thickness(pos) = 1;

            % Compute density ratio
            di = md.materials.rho_ice / md.materials.rho_water;

            % Compute ocean level set
            md.mask.ocean_levelset = md.geometry.thickness + md.geometry.bed / di;

            % Floating ice (ocean_levelset < 0)
            pos = find(md.mask.ocean_levelset < 0);
            md.geometry.surface(pos) = md.geometry.thickness(pos) .* ...
                (md.materials.rho_water - md.materials.rho_ice) / md.materials.rho_water;

            % Update base geometry
            md.geometry.base = md.geometry.surface - md.geometry.thickness;

            % Ensure base is not below bedrock
            pos = find(md.geometry.base < md.geometry.bed);
            md.geometry.base(pos) = md.geometry.base(pos);
            % md.geometry.base(pos) = md.geometry.bed(pos);

            % Grounded ice (ocean_levelset > 0)
            pos = find(md.mask.ocean_levelset > 0);
            md.geometry.base(pos) = md.geometry.bed(pos);

            % Update surface geometry
            md.geometry.surface = md.geometry.base + md.geometry.thickness;
            % md.geometry.surface = md.geometry.bed + md.geometry.thickness;

            md.smb.mass_balance=smb*ones(md.mesh.numberofvertices,1);
            md.transient.ismovingfront=0;
            % 
            md.basalforcings=linearbasalforcings();
            md.basalforcings.deepwater_melting_rate=deepwater_melting_rate;
            md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);

            % Time stepping
            md.timestepping = timestepping();
            md.timestepping.time_step = dt;
            md.timestepping.start_time = tinitial;
            md.timestepping.final_time = tfinal;
            md.settings.output_frequency = output_frequency;
            md.stressbalance.maxiter = 100;
            md.stressbalance.restol = 1;
            md.stressbalance.reltol = 0.001;
            md.stressbalance.abstol = NaN;

            % Cluster setup
            md.cluster = generic('name', oshostname(), 'np', nprocs);
            md.settings.waitonlock = Inf;
            md.settings.waitonlock=1;
            md.miscellaneous.name = sprintf('color_%d', ens_id);

            % Verbose settings
            md.verbose = verbose('convergence', false, 'solution', true);
            md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness', ...
                                            'Base','Bed', 'IceVolume', 'IceVolumeAboveFloatation'};
            % md.transient.requested_outputs = {'default','Thickness','Surface','Base','Bed'};

            % Solve transient
            md = solve(md, 'Transient','runtimename',false);

            % Save model
            filename = fullfile(folder, data_fname);
            save(filename, 'md', '-v7.3');

            % Save ensemble outputs in HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));

            result_0 = md.results.TransientSolution(end);

            data = cell(length(vec_inputs), 3);
            for j = 1:length(vec_inputs)
                key = vec_inputs{j};
                data{j, 1} = key;
                data{j, 2} = result_0;
                data{j, 3} = key;
            end
            
            writeToHDF5(filename, data);

            filename = fullfile(icesee_path, data_path, sprintf('ensemble_out_scalar_%d.h5', ens_id));
            h5write(filename, '/IceVolume', result_0.IceVolume, [k+1], [1]);
            h5write(filename, '/IceVolumeAboveFloatation', result_0.IceVolumeAboveFloatation, [k+1], [1]);

        end

    elseif strcmp(data_fname, 'inverse_state.mat')
        % folder = sprintf('./Models/ens_id_%d', ens_id_init);
        % if ~exist(folder, 'dir')
        %     mkdir(folder);
        % end
        % filename = fullfile(folder, reference_data);

        folder_true = sprintf('./Models/ens_id_%d', 0);
        % folder_true = sprintf('/Users/bkyanjo3/da_project/ISSM-matlab/examples/ISMIP_Choi/Models/ens_id_%d', 0);
        if ~exist(folder_true, 'dir')
            mkdir(folder_true);
        end
        % filename = fullfile(folder_true, 'true_state.mat');
        filename = fullfile(folder_true, reference_data);
              
        % load true state model for boundary conditions and other settings
        md = loadmodel(filename);

        vel_idx = double(kwargs.vel_idx);
        % km = double(kwargs.km);
        km = k+1; % matlab indexing starts at 1

        maxsteps = 40;

        % read in bed roughness data
        % filename = fullfile(icesee_path,'data/', 'synthetic_obs_0.h5');
        filename = fullfile(icesee_path, data_path, sprintf('synthetic_obs.h5'));
        obs_u = h5read(filename, '/hu_obs');
        nsize = md.mesh.numberofvertices;  % or: nsize = size(md.initialization.vx, 1);

        disp(['--- Ensemble ID: ', num2str(ens_id), '  Inverse Assimilation step: ', num2str(k)]);
        obs_col = obs_u(km,:)';            
     
        vx_obs = obs_col(vel_idx*nsize + 1 : (vel_idx+1)*nsize); 
        vy_obs = obs_col((vel_idx+1)*nsize + 1 : (vel_idx+2)*nsize);
        vel_obs = sqrt(vx_obs.^2 + vy_obs.^2);  

        % fetch the updated, vx, vy, h, s, bed, and base
        filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
        % md.geometry.thickness = h5read(filename, '/Thickness');
        % md.geometry.surface   = h5read(filename, '/Surface');
        md.initialization.vx   = h5read(filename, '/Vx');
        md.initialization.vy   = h5read(filename, '/Vy');
        md.initialization.vel  = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);
        % md.geometry.bed       = h5read(filename, '/bed');
        % md.geometry.base      = h5read(filename, '/Surface') - h5read(filename, '/Thickness');
        md.friction.coefficient = h5read(filename, '/FrictionCoefficient');
        % md.friction.coefficient = mean_friction*ones(md.mesh.numberofvertices,1);
        md.initialization.pressure=md.materials.rho_ice*md.constants.g*h5read(filename, '/Thickness');

        % Compute density ratio
        di = md.materials.rho_ice / md.materials.rho_water;

        % Compute ocean level set
        % md.mask.ocean_levelset = md.geometry.thickness + md.geometry.bed / di;
        md.mask.ocean_levelset = h5read(filename, '/Thickness') + h5read(filename, '/Bed') / di;

        % no friction applied on floating ice
        pos = find(md.mask.ocean_levelset < 0);
        md.friction.coefficient(pos)=0; %TODO: check the impact of this
        md.groundingline.migration='SubelementMigration';

        % set boundary conditions and other parameters

        md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices,1);
        md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);
        md.thermal.spctemperature                 = md.initialization.temperature;
        md.masstransport.spcthickness             = NaN*ones(md.mesh.numberofvertices,1);


        hvertices=[10000;500;5000;7500];
        gradation=1.7;
        err=8.0;
        md = bamg(md, 'domain', 'Domain.exp', 'hvertices',hvertices,'gradation',gradation,'field',md.initialization.vel,'err',err);
        % size(md.initialization.vx)


        %results of previous run are taken as observations
        md.inversion=m1qn3inversion();

        md.inversion.vx_obs = vx_obs;
        md.inversion.vy_obs = vy_obs;
        md.inversion.vel_obs = vel_obs;

        % Control general
        md.inversion.iscontrol=1;
        md.inversion.maxiter=40;
        md.inversion.dxmin=0.1;
        md.inversion.gttol=1.0e-4;
        md.verbose=verbose('control',true);

        md.inversion.maxsteps = maxsteps;
        md.inversion.cost_functions=[101 103 501];
        md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
        md.inversion.cost_functions_coefficients(:,1)=abs_vel_weight;
        md.inversion.cost_functions_coefficients(:,2)=rel_vel_weight;
        md.inversion.cost_functions_coefficients(:,3)=tikhonov_regularization_weight;

        md.inversion.control_parameters={'FrictionCoefficient'};
        md.inversion.min_parameters=min_friction*ones(md.mesh.numberofvertices,1);
        md.inversion.max_parameters=max_friction*ones(md.mesh.numberofvertices,1);

        md.stressbalance.restol=0.01;
        md.stressbalance.reltol=0.1;
        md.stressbalance.abstol=NaN;

        md.toolkits=toolkits;
        md.cluster=generic('name',oshostname,'np',nprocs);
        md.miscellaneous.name = sprintf('inverse_%d', ens_id);
        md=solve(md,'Stressbalance','runtimename',false);

        fcoef = md.friction.coefficient;
        md.friction.coefficient = md.results.StressbalanceSolution.FrictionCoefficient;

        md.initialization.vx = md.initialization.vx;
        md.initialization.vy = md.initialization.vy;
        md.geometry.thickness = h5read(filename, '/Thickness');
        md.geometry.surface   = h5read(filename, '/Surface');
        md.geometry.bed       = h5read(filename, '/Bed');


        % Save ensemble outputs in HDF5
        filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));

        % result_0 = md.results.TransientSolution(end);
        result_0 = md.initialization;
        % result_1 = md.results.TransientSolution(end);
        result_1 = md.geometry;
        result_2 = md.friction;
        % result_2 = md.results.TransientSolution(end);

        data = cell(length(vec_inputs), 3);
        for j = 1:length(vec_inputs)
            key = vec_inputs{j};
            data{j, 1} = key;
            data{j, 2} = result_0;
            data{j, 3} = key;
        end

        
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

function field = generate_correlated_field(md, ref_field, corr_length, std_dev)
%==========================================================================
% generate_correlated_field  (toolbox-free version)
%
% Purpose:
%   Generate a spatially correlated Gaussian random field suitable for EnKF
%   ensemble initialization in ISSM/ICESEE, without requiring pdist().
%
% Author:  Brian Kyanjo (2025)
%==========================================================================

    x = md.mesh.x(:);
    y = md.mesh.y(:);
    n = md.mesh.numberofvertices;

    rng('shuffle');

    % --- Compute pairwise distance matrix manually (memory-optimized) ---
    D2 = zeros(n, n);
    for i = 1:n
        dx = x - x(i);
        dy = y - y(i);
        D2(:, i) = dx.^2 + dy.^2;
    end

    % Gaussian covariance model
    C = exp(-D2 / (2 * corr_length^2));

    % Add small diagonal regularization
    C = C + 1e-6 * eye(n);

    % Cholesky decomposition (may require cholcov for large n)
    L = chol(C, 'lower');

    % Generate correlated perturbation
    z = randn(n,1);
    perturbation = std_dev * (L * z);

    % Apply perturbation
    field = ref_field + perturbation;

    % Enforce minimum physical constraint
    pos = find(field < 1);
    field(pos) = 1;

    disp(['[generate_correlated_field] Applied correlated noise with L = ' ...
        num2str(corr_length/1e3, '%.1f') ' km, std = ' num2str(std_dev) ' m']);
end
