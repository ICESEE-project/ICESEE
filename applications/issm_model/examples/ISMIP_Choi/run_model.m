
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


    % Get the current working directory
    cwd = pwd;
    [issmroot,~,~] = fileparts(fileparts(cwd));

    % number of variables
    nvar = 6;
    rng(1000 + ens_id, 'twister');   % ens_id = ensemble index

    % set initail ens_id
    ens_id_init = 0;
    s_perturb = double(kwargs.s_nurge);
    b_perturb = double(kwargs.b_nurge);

    output_frequency = 1; % make sure this is set to 1 for coupling with ICESEE

    % Set up model for each EnKF stage
    if strcmp(data_fname, 'initial_true_state.mat')
        folder = sprintf('./Models/ens_id_%d', ens_id_init);
        md = loadmodel(fullfile(folder, reference_data));
        writeInitialStateHDF5(fullfile(icesee_path, data_path, ...
            sprintf('ensemble_true_state_%d.h5', ens_id)), md);

    elseif strcmp(data_fname, 'initial_nurged_state.mat')
        folder = sprintf('./Models/ens_id_%d', ens_id_init);
        md = loadmodel(fullfile(folder, reference_data));
        md = setflowequation(md, 'SSA', 'all');
        prior_file = fullfile(icesee_path, data_path, ...
            sprintf('friction_bed_%d.h5', ens_id));
        bed = h5read(prior_file, '/bed');
        coefficient = h5read(prior_file, '/coefficient');
        md.friction.coefficient = mean_friction .* ...
            ones(md.mesh.numberofvertices, 1) + coefficient;
        md.friction.p = ones(md.mesh.numberofelements, 1);
        md.friction.q = ones(md.mesh.numberofelements, 1);
        md = apply_configured_initial_geometry(md, bed, kwargs);

        % Diagnose the velocity implied by this geometry without advancing
        % thickness, bed, grounding line, or time.
        md.cluster = generic('name', oshostname(), 'np', nprocs);
        md.settings.waitonlock = 1;
        md.verbose = verbose('convergence', false, 'solution', false);
        md = solve(md, 'Stressbalance', 'runtimename', false);
        md.initialization.vx = md.results.StressbalanceSolution.Vx;
        md.initialization.vy = md.results.StressbalanceSolution.Vy;
        md.initialization.vel = md.results.StressbalanceSolution.Vel;
        writeInitialStateHDF5(fullfile(icesee_path, data_path, ...
            sprintf('ensemble_nurged_state_%d.h5', ens_id)), md);

    elseif strcmp(data_fname, 'true_state.mat')
        % Special case for true state
        % if k == 0 || isempty(k)
        folder = sprintf('./Models/ens_id_%d', ens_id_init);
        if ~exist(folder, 'dir')
            mkdir(folder);
        end
    
        % Initial run: load boundary conditions
        filename = fullfile(folder, reference_data);
        md = loadmodel(filename);

        md = setflowequation(md,'SSA','all');

        md.smb.mass_balance=smb*ones(md.mesh.numberofvertices,1);
        md.transient.ismovingfront=0;
        % 
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

        md.transient.requested_outputs = {'default', 'FrictionCoefficient', 'Thickness', 'Surface','Base','Bed'};
        % md.transient.requested_outputs = {'default',  'Thickness', 'Surface','Base','Bed'};

        % Solve transient
        md = solve(md, 'Transient','runtimename',false);

        % update geometry
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

        % save updated model
        filename = fullfile(folder, data_fname);
        save(filename, 'md', '-v7.3');

        N = length(md.results.TransientSolution);
        % data = cell(N * 7, 7);   % 5 variables per step
        % data = cell(N*6,6); 
        data = cell(N * nvar, nvar);   % 5 variables per step

        idx = 1;
        for k = 1:N
            %  Thickness
            data{idx, 1} = sprintf('Thickness_%d', k);
            data{idx, 2} = md.results.TransientSolution(k);
            data{idx, 3} = 'Thickness';
            idx = idx + 1;

            % Base
            % data{idx, 1} = sprintf('Base_%d', k);
            % data{idx, 2} = md.results.TransientSolution(k);
            % data{idx, 3} = 'Base';
            % idx = idx + 1;

            % Surface
            data{idx, 1} = sprintf('Surface_%d', k);
            data{idx, 2} = md.results.TransientSolution(k);
            data{idx, 3} = 'Surface';
            idx = idx + 1;

            % Vx 
            data{idx, 1} = sprintf('Vx_%d', k);
            data{idx, 2} = md.results.TransientSolution(k);
            data{idx, 3} = 'Vx';    
            idx = idx + 1;

            % Vy
            data{idx, 1} = sprintf('Vy_%d', k);
            data{idx, 2} = md.results.TransientSolution(k);
            data{idx, 3} = 'Vy';
            idx = idx + 1;

            % Bed (from results)
            data{idx, 1} = sprintf('bed_%d', k);
            data{idx, 2} = md.results.TransientSolution(k);
            data{idx, 3} = 'Bed';
            idx = idx + 1;

            % Friction Coefficient (from results)
            data{idx, 1} = sprintf('coefficient_%d', k);
            data{idx, 2} = md.results.TransientSolution(k);
            data{idx, 3} = 'FrictionCoefficient';
            idx = idx + 1;
        end

        filename = fullfile(icesee_path, data_path, sprintf('ensemble_true_state_%d.h5', ens_id));
        writeToHDF5(filename, data);

    elseif strcmp(data_fname, 'nurged_state.mat')

        folder = sprintf('./Models/ens_id_%d', ens_id_init);
        if ~exist(folder, 'dir')
            mkdir(folder);
        end
            
        filename = fullfile(folder, reference_data);
        md = loadmodel(filename);

        md = setflowequation(md,'SSA','all');

        % setup nugged state
        friction_ref = mean_friction * ones(md.mesh.numberofvertices,1);

        filename = fullfile(icesee_path, data_path, sprintf('friction_bed_%d.h5', ens_id));
        bed = h5read(filename, '/bed');
        coefficient = h5read(filename, '/coefficient');

        md.friction.coefficient = friction_ref + coefficient;
        md.friction.p = ones(md.mesh.numberofelements,1);
        md.friction.q = ones(md.mesh.numberofelements,1);

        md = apply_configured_initial_geometry(md, bed, kwargs);

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

        md.transient.requested_outputs = {'default', 'FrictionCoefficient', 'Thickness', 'Surface','Base','Bed'};
        % md.transient.requested_outputs = {'default',  'Thickness', 'Surface','Base','Bed'};

        % Solve transient
        md = solve(md, 'Transient','runtimename',false);
            
        filename = fullfile(folder, data_fname);
        save(filename, 'md', '-v7.3');

        N = length(md.results.TransientSolution);
        data = cell(N * nvar, nvar);   % 5 variables per step

        idx = 1;
        for k = 1:N
            %  Thickness
            data{idx, 1} = sprintf('Thickness_%d', k);
            data{idx, 2} = md.results.TransientSolution(k);
            data{idx, 3} = 'Thickness';
            idx = idx + 1;

            % Base
            % data{idx, 1} = sprintf('Base_%d', k);
            % data{idx, 2} = md.results.TransientSolution(k);
            % data{idx, 3} = 'Base';  
            % idx = idx + 1;

            % Surface
            data{idx, 1} = sprintf('Surface_%d', k);
            data{idx, 2} = md.results.TransientSolution(k);
            data{idx, 3} = 'Surface';
            idx = idx + 1;

            % Vx
            data{idx, 1} = sprintf('Vx_%d', k);
            data{idx, 2} = md.results.TransientSolution(k);
            data{idx, 3} = 'Vx';
            idx = idx + 1;

            % Vy
            data{idx, 1} = sprintf('Vy_%d', k);
            data{idx, 2} = md.results.TransientSolution(k);
            data{idx, 3} = 'Vy';
            idx = idx + 1;

            % Bed (from results)
            data{idx, 1} = sprintf('bed_%d', k);
            data{idx, 2} = md.results.TransientSolution(k);
            data{idx, 3} = 'Bed';
            idx = idx + 1;

            % Friction Coefficient (from results)
            data{idx, 1} = sprintf('coefficient_%d', k);
            data{idx, 2} = md.results.TransientSolution(k);
            data{idx, 3} = 'FrictionCoefficient';
            idx = idx + 1;
        end


        filename = fullfile(icesee_path, data_path, sprintf('ensemble_nurged_state_%d.h5', ens_id));
        writeToHDF5(filename, data);


    elseif strcmp(data_fname, 'initialize_ensemble.mat')
        % Special case for ensemble initialization
        if k == 0 || isempty(k)

            % % call the wrong data and only fetch the first iteration
            % filename = fullfile(icesee_path, data_path, 'true_nurged_states.h5');
            % model_nurged_state = h5read(filename,'/nurged_state')';

            % k = 1; % first time step
            % [nd, nt] = size(model_nurged_state);
            % nvar = 6; % thickness, surface, base, Vx, Vy, bed, coefficient  
            % hdim = nd / nvar; % number of vertices (assuming 6 variables: thickness, surface, base, Vx, Vy, bed, coefficient)
            % H  = model_nurged_state(1:hdim, k);
            % S  = model_nurged_state(hdim+1:2*hdim, k);
            % B  = S - H;
            % Vx = model_nurged_state(2*hdim+1:3*hdim, k);
            % Vy = model_nurged_state(3*hdim+1:4*hdim, k);
            % Vel= hypot(Vx, Vy);
            % % bed= model_nurged_state(4*hdim+1:5*hdim, k);
            % % fc = model_nurged_state(5*hdim+1:6*hdim, k);

            % filename = fullfile(icesee_path, data_path, sprintf('friction_bed_%d.h5', ens_id));
            % bed = h5read(filename, '/bed');
            % fc = h5read(filename, '/coefficient');

            %  % Ensure base not below bedrock
            % pos = find(B < bed);
            % B(pos) = bed(pos);

            % folder_true = sprintf('./Models/ens_id_%d', 0);
            % % folder_true = sprintf('/Users/bkyanjo3/da_project/ISSM-matlab/examples/ISMIP_Choi/Models/ens_id_%d', 0);
            % if ~exist(folder_true, 'dir')
            %     mkdir(folder_true);
            % end
            % % filename = fullfile(folder_true, 'true_state.mat');
            % filename = fullfile(folder_true, reference_data);
                
            % % load true state model for boundary conditions and other settings
            % md = loadmodel(filename);

            % % Grounded ice (ocean_levelset > 0)
            % ocean_levelset = H + bed/(md.materials.rho_ice/md.materials.rho_water);
            % pos = find(ocean_levelset > 0);
            % B(pos) = bed(pos);

            % md.geometry.thickness    = H;
            % md.geometry.surface      = S;
            % md.geometry.base         = md.geometry.base - 0.25*bed;
            % md.geometry.bed          = md.geometry.bed - 0.25*bed;
            % md.initialization.vx     = Vx;
            % md.initialization.vy     = Vy;
            % md.initialization.vel    = Vel;
            % md.friction.coefficient  = fc;
            % di = md.materials.rho_ice / md.materials.rho_water;
            % md.mask.ocean_levelset   = H + bed/di;
      
            % filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
            % data = {'Thickness', md.geometry, 'thickness';
            %         'Surface', md.geometry, 'surface';
            %         'Base', md.geometry, 'base';
            %         'bed', md.geometry, 'bed';
            %         'Vx', md.initialization, 'vx';
            %         'Vy', md.initialization, 'vy';
            %         'Vel', md.initialization, 'vel';
            %         'coefficient', md.friction, 'coefficient'
            % };
            % writeToHDF5(filename, data);

            % filename = fullfile(icesee_path, data_path, sprintf('ensemble_out_%d.h5', ens_id));
            % writeToHDF5(filename, data);


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

             % % read the friction_bed file
            filename = fullfile(icesee_path, data_path, sprintf('friction_bed_%d.h5', ens_id));
            bed = h5read(filename, '/bed');
            coefficient = h5read(filename, '/coefficient');

            %  update the friction and bed
            md.friction.coefficient = friction_ref + coefficient;
            % md.friction.coefficient = friction_ref;
            md.friction.p=ones(md.mesh.numberofelements,1);
            md.friction.q=ones(md.mesh.numberofelements,1);

 
            md = apply_configured_initial_geometry(md, bed, kwargs);

            % pos = find(md.mask.ocean_levelset < 0);
            % md.geometry.thickness(pos)=1/(1-di)*md.geometry.surface(pos);


            % --time stepping
            % dt/tinitial/tfinal are supplied by initialize_ensemble().  Their
            % historical defaults still give one 0.2-year step, while an
            % experiment may request a longer pre-DA dynamic spin-up.
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

            md.transient.requested_outputs = {'default', 'FrictionCoefficient', 'Thickness', 'Surface','Base','Bed'};
            % md.transient.requested_outputs = {'default',  'Thickness', 'Surface','Base','Bed'};

            % Solve transient
            md = solve(md, 'Transient','runtimename',false);
            spinup_speed = hypot(md.results.TransientSolution(end).Vx, ...
                                 md.results.TransientSolution(end).Vy);
            fprintf(['[ICESEE] Ensemble %d initialization spin-up: ', ...
                     'duration=%.6g yr, dt=%.6g yr, max(speed)=%.6g m/yr\n'], ...
                    ens_id, tfinal - tinitial, dt, max(spinup_speed));
             
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

            data = {'Thickness', result_0, 'Thickness';
                    % 'Base', result_0, 'Base';
                    'Surface', result_0, 'Surface';
                    'Vx', result_0, 'Vx';
                    'Vy', result_0, 'Vy';
                    'bed', result_0, 'Bed';
                    'coefficient', result_0, 'FrictionCoefficient'
            };

            writeToHDF5(filename, data);

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
            
            % md.inversion.iscontrol            = 0;
            % md.transient.ismovingfront        = 0;
            % md.transient.isthermal            = 0;
            % md.transient.isstressbalance      = 1;
            % md.transient.ismasstransport      = 1;
            % md.transient.isgroundingline      = 1;

            % md.groundingline.migration                = 'SubelementMigration';
            % md.groundingline.friction_interpolation   = 'SubelementFriction1';
            % md.groundingline.melt_interpolation       = 'NoMeltOnPartiallyFloating';

            % md.initialization.pressure       = zeros(md.mesh.numberofvertices, 1);
            % md.masstransport.spcthickness    = NaN * ones(md.mesh.numberofvertices, 1);

            % md.verbose.solution              = 1;

            % mask_all = zeros(md.mesh.numberofvertices,1);
            % md.smb.mass_balance=smb*ones(md.mesh.numberofvertices,1);
            % md.basalforcings=linearbasalforcings();
            % md.basalforcings.deepwater_melting_rate=deepwater_melting_rate;
            % md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);


             % Load ensemble input from HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
            % Keep the last state accepted by ISSM.  The ensemble HDF5 file is
            % written by the filter and can contain a one-cycle member outlier;
            % the model file remains a safe persistence fallback.
            fallback_vx = md.initialization.vx;
            fallback_vy = md.initialization.vy;
            md.geometry.surface = h5read(filename, '/Surface');
            % md.geometry.base = h5read(filename, '/Base');
            md.geometry.thickness = h5read(filename, '/Thickness');
            md.initialization.vx = h5read(filename, '/Vx');
            md.initialization.vy = h5read(filename, '/Vy');
            md.initialization.vel = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);
            % md.initialization.pressure=md.materials.rho_ice*md.constants.g*h5read(filename, '/Thickness');
        
            % parameters for bed and friction
            md.geometry.bed = h5read(filename, '/bed');
            md.friction.coefficient = h5read(filename, '/coefficient');

            % Reject a divergent filter member before ISSM's less-informative
            % geometry-consistency check. These limits are deliberately far
            % outside the ISMIP-Choi solution range.
            input_speed = hypot(md.initialization.vx, md.initialization.vy);
            geometry_is_bad = any(~isfinite(md.geometry.thickness)) || ...
                    any(~isfinite(md.geometry.surface)) || ...
                    max(md.geometry.thickness) > 1.0e4 || ...
                    max(abs(md.geometry.surface)) > 2.0e4;
            velocity_is_bad = any(~isfinite(input_speed)) || ...
                    max(input_speed) > 5.0e3;
            fallback_speed = hypot(fallback_vx, fallback_vy);
            if ~geometry_is_bad && velocity_is_bad && ...
                    all(isfinite(fallback_speed)) && max(fallback_speed) <= 5.0e3
                warning(['[ICESEE] Rejecting catastrophic member velocity ', ...
                         '(max %.6g m/yr); retaining the last ISSM-accepted ', ...
                         'velocity for this forecast.'], max(input_speed));
                md.initialization.vx = fallback_vx;
                md.initialization.vy = fallback_vy;
                md.initialization.vel = fallback_speed;
                input_speed = fallback_speed;
                velocity_is_bad = false;
            end
            if geometry_is_bad || velocity_is_bad
                error(['[ICESEE] Catastrophic EnKF member detected before ISSM: ', ...
                       'max(H)=%.6g m, max(abs(S))=%.6g m, max(speed)=%.6g m/yr'], ...
                      max(md.geometry.thickness), max(abs(md.geometry.surface)), ...
                      max(input_speed));
            end

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
            % md.geometry.base(pos) = md.geometry.bed(pos);

            % Grounded ice (ocean_levelset > 0)
            pos = find(md.mask.ocean_levelset > 0);
            md.geometry.base(pos) = md.geometry.bed(pos);

            % Update surface geometry
            md.geometry.surface = md.geometry.base + md.geometry.thickness;
            % Make ISSM's checked identity exact in its own evaluation order.
            md.geometry.thickness = md.geometry.surface - md.geometry.base;
            md.initialization.pressure = md.materials.rho_ice * ...
                md.constants.g * md.geometry.thickness;

            % % Outputs and verbosity
            md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
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

            % Retain a physically accepted state so that a rare nonlinear
            % transient burst cannot poison the next EnKF cycle.
            accepted_thickness = md.geometry.thickness;
            accepted_surface = md.geometry.surface;
            accepted_vx = md.initialization.vx;
            accepted_vy = md.initialization.vy;

            % % Solve transient
            md = solve(md, 'Transient','runtimename',false); %TODO: instead of solving just take th initial solution

            result_speed = hypot(md.results.TransientSolution(end).Vx, ...
                                 md.results.TransientSolution(end).Vy);
            if any(~isfinite(result_speed)) || max(result_speed) > 5.0e3
                warning(['[ICESEE] Rejecting catastrophic ISSM forecast velocity ', ...
                         '(max %.6g m/yr) for member %d; using persistence ', ...
                         'for this member and cycle.'], max(result_speed), ens_id);
                md.results.TransientSolution(end).Thickness = accepted_thickness;
                md.results.TransientSolution(end).Surface = accepted_surface;
                md.results.TransientSolution(end).Base = accepted_surface - accepted_thickness;
                md.results.TransientSolution(end).Vx = accepted_vx;
                md.results.TransientSolution(end).Vy = accepted_vy;
                md.results.TransientSolution(end).Vel = hypot(accepted_vx, accepted_vy);
            end

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

            data = {'Thickness', result_0, 'Thickness';
                    % 'Base', result_0, 'Base';
                    'Surface', result_0, 'Surface';
                    'Vx', result_0, 'Vx';
                    'Vy', result_0, 'Vy';
                    'bed', result_1, 'Bed';
                    % 'coefficient', result_2, 'coefficient'};
                    'coefficient', result_2, 'FrictionCoefficient'};
                    

            writeToHDF5(filename, data);
        % ;

        else
          
            % fprintf('[MATLAB ---] Running model for ensemble ID %d, step %d\n', ens_id, k);
            
            % Subsequent time steps: 
            filename = fullfile(folder, data_fname);
            if exist(filename, 'file')
                md = loadmodel(filename);
            else
                % A partial-run resume deliberately skips ensemble
                % initialization, so its per-member MAT cache may not exist.
                % Rebuild only the model container from the verified truth
                % model; all six evolving member fields are replaced from the
                % resumed ensemble HDF5 immediately below.
                bootstrap_file = fullfile('./Models/ens_id_0', 'true_state.mat');
                if ~exist(bootstrap_file, 'file')
                    error(['[ICESEE] Cannot bootstrap resumed member: ', ...
                           '%s does not exist'], bootstrap_file);
                end
                warning(['[ICESEE] Rebuilding missing member model cache ', ...
                         '%s from %s'], filename, bootstrap_file);
                md = loadmodel(bootstrap_file);
            end

            % Load ensemble input from HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
            % Keep the last state accepted by ISSM as a member-level fallback.
            fallback_vx = md.initialization.vx;
            fallback_vy = md.initialization.vy;
            md.geometry.surface = h5read(filename, '/Surface');
            % md.geometry.base = h5read(filename, '/Base');
            md.geometry.thickness = h5read(filename, '/Thickness');
            md.initialization.vx = h5read(filename, '/Vx');
            md.initialization.vy = h5read(filename, '/Vy');
            md.initialization.vel = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);
            md.initialization.pressure=md.materials.rho_ice*md.constants.g*h5read(filename, '/Thickness');
        
            % parameters for bed and friction
            md.geometry.bed = h5read(filename, '/bed');
            md.friction.coefficient = h5read(filename, '/coefficient');

            % Reject a divergent filter member before ISSM's less-informative
            % geometry-consistency check. These limits are deliberately far
            % outside the ISMIP-Choi solution range.
            input_speed = hypot(md.initialization.vx, md.initialization.vy);
            geometry_is_bad = any(~isfinite(md.geometry.thickness)) || ...
                    any(~isfinite(md.geometry.surface)) || ...
                    max(md.geometry.thickness) > 1.0e4 || ...
                    max(abs(md.geometry.surface)) > 2.0e4;
            velocity_is_bad = any(~isfinite(input_speed)) || ...
                    max(input_speed) > 5.0e3;
            fallback_speed = hypot(fallback_vx, fallback_vy);
            if ~geometry_is_bad && velocity_is_bad && ...
                    all(isfinite(fallback_speed)) && max(fallback_speed) <= 5.0e3
                warning(['[ICESEE] Rejecting catastrophic member velocity ', ...
                         '(max %.6g m/yr); retaining the last ISSM-accepted ', ...
                         'velocity for this forecast.'], max(input_speed));
                md.initialization.vx = fallback_vx;
                md.initialization.vy = fallback_vy;
                md.initialization.vel = fallback_speed;
                input_speed = fallback_speed;
                velocity_is_bad = false;
            end
            if geometry_is_bad || velocity_is_bad
                error(['[ICESEE] Catastrophic EnKF member detected before ISSM: ', ...
                       'max(H)=%.6g m, max(abs(S))=%.6g m, max(speed)=%.6g m/yr'], ...
                      max(md.geometry.thickness), max(abs(md.geometry.surface)), ...
                      max(input_speed));
            end

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
            % Make ISSM's checked identity exact in its own evaluation order.
            md.geometry.thickness = md.geometry.surface - md.geometry.base;
            md.initialization.pressure = md.materials.rho_ice * ...
                md.constants.g * md.geometry.thickness;

            md.smb.mass_balance=smb*ones(md.mesh.numberofvertices,1);
            % md.transient.ismovingfront=0;
            % 
            % md.basalforcings=linearbasalforcings();
            md.basalforcings.deepwater_melting_rate=deepwater_melting_rate;
            % md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);

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
            md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
            % md.transient.requested_outputs = {'default','Thickness','Surface','Base','Bed'};

            % Retain a physically accepted state so that a rare nonlinear
            % transient burst cannot poison the next EnKF cycle.
            accepted_thickness = md.geometry.thickness;
            accepted_surface = md.geometry.surface;
            accepted_vx = md.initialization.vx;
            accepted_vy = md.initialization.vy;

            % Solve transient
            md = solve(md, 'Transient','runtimename',false);

            result_speed = hypot(md.results.TransientSolution(end).Vx, ...
                                 md.results.TransientSolution(end).Vy);
            if any(~isfinite(result_speed)) || max(result_speed) > 5.0e3
                warning(['[ICESEE] Rejecting catastrophic ISSM forecast velocity ', ...
                         '(max %.6g m/yr) for member %d; using persistence ', ...
                         'for this member and cycle.'], max(result_speed), ens_id);
                md.results.TransientSolution(end).Thickness = accepted_thickness;
                md.results.TransientSolution(end).Surface = accepted_surface;
                md.results.TransientSolution(end).Base = accepted_surface - accepted_thickness;
                md.results.TransientSolution(end).Vx = accepted_vx;
                md.results.TransientSolution(end).Vy = accepted_vy;
                md.results.TransientSolution(end).Vel = hypot(accepted_vx, accepted_vy);
            end

            % Save model
            filename = fullfile(folder, data_fname);
            save(filename, 'md', '-v7.3');

            % md = transientrestart(md);
            md.geometry.thickness = md.results.TransientSolution(end).Thickness;
            md.geometry.surface   = md.results.TransientSolution(end).Surface;
            md.geometry.base      = md.results.TransientSolution(end).Base;

            % Update other fields
            md.initialization.vx        = md.results.TransientSolution(end).Vx;
            md.initialization.vy        = md.results.TransientSolution(end).Vy;
            md.initialization.vel       = md.results.TransientSolution(end).Vel;
            % md.initialization.pressure  = md.results.TransientSolution(end).Pressure;
            % md.smb.mass_balance         = md.results.TransientSolution(end).SmbMassBalance;
            md.mask.ocean_levelset      = md.results.TransientSolution(end).MaskOceanLevelset;

            md.geometry.bed = md.results.TransientSolution(end).Bed;
            md.friction.coefficient = md.results.TransientSolution(end).FrictionCoefficient;

            % *--
            % Ensure minimum ice thickness of 1 m
            % pos = find(md.geometry.thickness < 1);
            % md.geometry.thickness(pos) = 1;

            % % Density ratio
            % di = md.materials.rho_ice / md.materials.rho_water;

            % % Compute ocean level set based on hydrostatic equilibrium
            % md.mask.ocean_levelset = md.geometry.thickness + md.geometry.bed / di;

            % % Floating ice (ocean_levelset < 0)
            % pos = find(md.mask.ocean_levelset < 0);
            % md.geometry.surface(pos) = md.geometry.thickness(pos) .* ...
            %     (md.materials.rho_water - md.materials.rho_ice) / md.materials.rho_water;

            % % Update base geometry
            % md.geometry.base = md.geometry.surface - md.geometry.thickness;

            % % Ensure base not below bedrock
            % pos = find(md.geometry.base < md.geometry.bed);
            % % md.geometry.base(pos) = md.geometry.base(pos);
            % md.geometry.base(pos) = md.geometry.bed(pos);

            % % Grounded ice (ocean_levelset > 0)
            % pos = find(md.mask.ocean_levelset > 0);
            % md.geometry.base(pos) = md.geometry.bed(pos);

            % % Update surface geometry
            % md.geometry.surface = md.geometry.base + md.geometry.thickness;

            % Save ensemble outputs in HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));

            result_0 = md.results.TransientSolution(end);
            % result_0 = md.initialization;
            result_1 = md.results.TransientSolution(end);
            % result_1 = md.geometry;
            % result_2 = md.friction;
            result_2 = md.results.TransientSolution(end);

            data = {'Thickness', result_1, 'Thickness';
                    % 'Base', result_1, 'base';
                    'Surface', result_1, 'Surface';
                    'Vx', result_0, 'Vx';
                    'Vy', result_0, 'Vy';
                    'bed', result_1, 'Bed';
                    'coefficient', result_2, 'FrictionCoefficient'};
            

            writeToHDF5(filename, data);

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
        md.friction.coefficient = h5read(filename, '/coefficient');
        % md.friction.coefficient = mean_friction*ones(md.mesh.numberofvertices,1);
        md.initialization.pressure=md.materials.rho_ice*md.constants.g*h5read(filename, '/Thickness');

        % Compute density ratio
        di = md.materials.rho_ice / md.materials.rho_water;

        % Compute ocean level set
        % md.mask.ocean_levelset = md.geometry.thickness + md.geometry.bed / di;
        md.mask.ocean_levelset = h5read(filename, '/Thickness') + h5read(filename, '/bed') / di;

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
        md.geometry.bed       = h5read(filename, '/bed');


        % Save ensemble outputs in HDF5
        filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));

        % result_0 = md.results.TransientSolution(end);
        result_0 = md.initialization;
        % result_1 = md.results.TransientSolution(end);
        result_1 = md.geometry;
        result_2 = md.friction;
        % result_2 = md.results.TransientSolution(end);

        data = {'Thickness', result_1, 'thickness';
                % 'Base', result_1, 'base';
                'Surface', result_1, 'surface';
                'Vx', result_0, 'vx';
                'Vy', result_0, 'vy';
                'bed', result_1, 'bed';
                'coefficient', result_2, 'coefficient'};
        
        writeToHDF5(filename, data);
    end
end

function md = apply_configured_initial_geometry(md, bed_candidate, kwargs)
%APPLY_CONFIGURED_INITIAL_GEOMETRY Build one physically consistent prior.
% The same construction is used for the no-assimilation trajectory and each
% ensemble member.  Surface is never biased independently: it is recovered
% from bed/base and thickness, and floating ice is put in hydrostatic
% equilibrium.  Defaults reproduce the historical behaviour.

    thickness_scale = 1.0;
    bed_offset_m = 0.0;
    bed_domain = 'all';
    bed_gl_buffer_m = 0.0;
    floating_bed_anomaly_factor = 0.0;
    floating_bed_max_error_m = 100.0;
    floating_bed_transition_m = 25000.0;
    floating_bed_flotation_margin_m = 5.0;
    bed_smoothing_iterations = 35;
    bed_smoothing_strength = 0.65;
    bed_seed_max_x_m = 300000.0;
    bed_downstream_anomaly_factor = 0.60;
    thickness_anomaly_fraction = 0.0;
    thickness_anomaly_m = 0.0;
    thickness_delta_min_m = -500.0;
    thickness_delta_max_m = 500.0;
    floating_thickness_anomaly_factor = 1.0;
    gl_seaward_thickness_m = 0.0;
    gl_seaward_width_m = 50000.0;
    bed_anomaly_m = 0.0;
    bed_delta_min_m = -500.0;
    bed_delta_max_m = 500.0;
    pattern_length_x_m = 120000.0;
    pattern_length_y_m = 40000.0;
    pattern_phase = 0.0;
    thickness_factor_min = 0.60;
    thickness_factor_max = 1.25;
    if isfield(kwargs, 'initial_thickness_scale')
        thickness_scale = double(kwargs.initial_thickness_scale);
    end
    if isfield(kwargs, 'initial_bed_offset_m')
        bed_offset_m = double(kwargs.initial_bed_offset_m);
    end
    if isfield(kwargs, 'initial_bed_background_domain')
        bed_domain = lower(strtrim(char(kwargs.initial_bed_background_domain)));
    end
    if isfield(kwargs, 'initial_bed_gl_buffer_m')
        bed_gl_buffer_m = double(kwargs.initial_bed_gl_buffer_m);
    end
    if isfield(kwargs, 'initial_floating_bed_anomaly_factor')
        floating_bed_anomaly_factor = double( ...
            kwargs.initial_floating_bed_anomaly_factor);
    end
    if isfield(kwargs, 'initial_floating_bed_max_error_m')
        floating_bed_max_error_m = double( ...
            kwargs.initial_floating_bed_max_error_m);
    end
    if isfield(kwargs, 'initial_floating_bed_transition_m')
        floating_bed_transition_m = double( ...
            kwargs.initial_floating_bed_transition_m);
    end
    if isfield(kwargs, 'initial_floating_bed_flotation_margin_m')
        floating_bed_flotation_margin_m = double( ...
            kwargs.initial_floating_bed_flotation_margin_m);
    end
    if isfield(kwargs, 'initial_bed_smoothing_iterations')
        bed_smoothing_iterations = double( ...
            kwargs.initial_bed_smoothing_iterations);
    end
    if isfield(kwargs, 'initial_bed_smoothing_strength')
        bed_smoothing_strength = double( ...
            kwargs.initial_bed_smoothing_strength);
    end
    if isfield(kwargs, 'initial_bed_seed_max_x_m')
        bed_seed_max_x_m = double(kwargs.initial_bed_seed_max_x_m);
    end
    if isfield(kwargs, 'initial_bed_downstream_anomaly_factor')
        bed_downstream_anomaly_factor = double( ...
            kwargs.initial_bed_downstream_anomaly_factor);
    end
    if isfield(kwargs, 'initial_thickness_anomaly_fraction')
        thickness_anomaly_fraction = double(kwargs.initial_thickness_anomaly_fraction);
    end
    if isfield(kwargs, 'initial_thickness_anomaly_m')
        thickness_anomaly_m = double(kwargs.initial_thickness_anomaly_m);
    end
    if isfield(kwargs, 'initial_thickness_delta_min_m')
        thickness_delta_min_m = double(kwargs.initial_thickness_delta_min_m);
    end
    if isfield(kwargs, 'initial_thickness_delta_max_m')
        thickness_delta_max_m = double(kwargs.initial_thickness_delta_max_m);
    end
    if isfield(kwargs, 'initial_floating_thickness_anomaly_factor')
        floating_thickness_anomaly_factor = double( ...
            kwargs.initial_floating_thickness_anomaly_factor);
    end
    if isfield(kwargs, 'initial_gl_seaward_thickness_m')
        gl_seaward_thickness_m = double( ...
            kwargs.initial_gl_seaward_thickness_m);
    end
    if isfield(kwargs, 'initial_gl_seaward_width_m')
        gl_seaward_width_m = double(kwargs.initial_gl_seaward_width_m);
    end
    if isfield(kwargs, 'initial_bed_anomaly_m')
        bed_anomaly_m = double(kwargs.initial_bed_anomaly_m);
    end
    if isfield(kwargs, 'initial_bed_delta_min_m')
        bed_delta_min_m = double(kwargs.initial_bed_delta_min_m);
    end
    if isfield(kwargs, 'initial_bed_delta_max_m')
        bed_delta_max_m = double(kwargs.initial_bed_delta_max_m);
    end
    if isfield(kwargs, 'initial_prior_length_x_m')
        pattern_length_x_m = double(kwargs.initial_prior_length_x_m);
    end
    if isfield(kwargs, 'initial_prior_length_y_m')
        pattern_length_y_m = double(kwargs.initial_prior_length_y_m);
    end
    if isfield(kwargs, 'initial_prior_pattern_phase')
        pattern_phase = double(kwargs.initial_prior_pattern_phase);
    end
    if isfield(kwargs, 'initial_thickness_factor_min')
        thickness_factor_min = double(kwargs.initial_thickness_factor_min);
    end
    if isfield(kwargs, 'initial_thickness_factor_max')
        thickness_factor_max = double(kwargs.initial_thickness_factor_max);
    end
    if ~isfinite(thickness_scale) || thickness_scale <= 0 || thickness_scale > 2
        error('[ICESEE] initial_thickness_scale must be in (0, 2].');
    end
    if ~isfinite(bed_offset_m) || abs(bed_offset_m) > 2000
        error('[ICESEE] initial_bed_offset_m must be finite and within 2000 m.');
    end
    if ~isfinite(bed_gl_buffer_m) || bed_gl_buffer_m < 0 || ...
            bed_gl_buffer_m > 200000
        error('[ICESEE] initial_bed_gl_buffer_m must be in [0, 200000] m.');
    end
    if ~isfinite(floating_bed_anomaly_factor) || ...
            floating_bed_anomaly_factor < 0 || ...
            floating_bed_anomaly_factor > 1
        error(['[ICESEE] initial_floating_bed_anomaly_factor must be ', ...
               'in [0, 1].']);
    end
    if ~isfinite(floating_bed_max_error_m) || ...
            floating_bed_max_error_m <= 0 || ...
            floating_bed_max_error_m > 1000
        error(['[ICESEE] initial_floating_bed_max_error_m must be ', ...
               'in (0, 1000] m.']);
    end
    if ~isfinite(floating_bed_transition_m) || ...
            floating_bed_transition_m <= 0 || ...
            floating_bed_transition_m > 200000
        error(['[ICESEE] initial_floating_bed_transition_m must be ', ...
               'in (0, 200000] m.']);
    end
    if ~isfinite(floating_bed_flotation_margin_m) || ...
            floating_bed_flotation_margin_m < 0 || ...
            floating_bed_flotation_margin_m > 100
        error(['[ICESEE] initial_floating_bed_flotation_margin_m must be ', ...
               'in [0, 100] m.']);
    end
    if ~isfinite(bed_smoothing_iterations) || ...
            bed_smoothing_iterations < 0 || ...
            bed_smoothing_iterations > 200 || ...
            bed_smoothing_iterations ~= floor(bed_smoothing_iterations)
        error(['[ICESEE] initial_bed_smoothing_iterations must be an ', ...
               'integer in [0, 200].']);
    end
    if ~isfinite(bed_smoothing_strength) || ...
            bed_smoothing_strength < 0 || bed_smoothing_strength > 1
        error(['[ICESEE] initial_bed_smoothing_strength must be in ', ...
               '[0, 1].']);
    end
    if ~isfinite(bed_seed_max_x_m) || bed_seed_max_x_m <= 0
        error('[ICESEE] initial_bed_seed_max_x_m must be positive.');
    end
    if ~isfinite(bed_downstream_anomaly_factor) || ...
            bed_downstream_anomaly_factor < 0 || ...
            bed_downstream_anomaly_factor > 1
        error(['[ICESEE] initial_bed_downstream_anomaly_factor must be ', ...
               'in [0, 1].']);
    end
    if ~isfinite(thickness_anomaly_fraction) || ...
            thickness_anomaly_fraction < 0 || thickness_anomaly_fraction > 0.5
        error('[ICESEE] initial_thickness_anomaly_fraction must be in [0, 0.5].');
    end
    if ~isfinite(thickness_anomaly_m) || thickness_anomaly_m < 0 || ...
            thickness_anomaly_m > 1000
        error('[ICESEE] initial_thickness_anomaly_m must be in [0, 1000] m.');
    end
    if ~isfinite(thickness_delta_min_m) || ...
            ~isfinite(thickness_delta_max_m) || ...
            thickness_delta_max_m <= thickness_delta_min_m
        error('[ICESEE] Invalid initial thickness-delta bounds.');
    end
    if ~isfinite(floating_thickness_anomaly_factor) || ...
            floating_thickness_anomaly_factor < 0 || ...
            floating_thickness_anomaly_factor > 1
        error('[ICESEE] initial_floating_thickness_anomaly_factor must be in [0, 1].');
    end
    if ~isfinite(gl_seaward_thickness_m) || ...
            gl_seaward_thickness_m < 0 || gl_seaward_thickness_m > 1000
        error('[ICESEE] initial_gl_seaward_thickness_m must be in [0, 1000] m.');
    end
    if ~isfinite(gl_seaward_width_m) || ...
            gl_seaward_width_m <= 0 || gl_seaward_width_m > 200000
        error('[ICESEE] initial_gl_seaward_width_m must be in (0, 200000] m.');
    end
    if ~isfinite(bed_anomaly_m) || bed_anomaly_m < 0 || bed_anomaly_m > 1000
        error('[ICESEE] initial_bed_anomaly_m must be in [0, 1000] m.');
    end
    if ~isfinite(bed_delta_min_m) || ~isfinite(bed_delta_max_m) || ...
            bed_delta_max_m <= bed_delta_min_m
        error('[ICESEE] Invalid initial bed-delta bounds.');
    end
    if ~isfinite(pattern_length_x_m) || pattern_length_x_m <= 0 || ...
            ~isfinite(pattern_length_y_m) || pattern_length_y_m <= 0
        error('[ICESEE] Initial-prior pattern lengths must be positive.');
    end
    if ~isfinite(thickness_factor_min) || ~isfinite(thickness_factor_max) || ...
            thickness_factor_min <= 0 || ...
            thickness_factor_max <= thickness_factor_min
        error('[ICESEE] Invalid initial thickness-factor bounds.');
    end

    bed_background = md.geometry.bed(:);
    bed_candidate = double(bed_candidate(:));
    if numel(bed_candidate) ~= md.mesh.numberofvertices
        error('[ICESEE] Initial bed candidate has an incompatible size.');
    end
    initial_grounded = md.geometry.thickness(:) > 0 & ...
                       md.mask.ocean_levelset(:) > 0;
    tapered_floating_bed = false;
    switch bed_domain
        case 'all'
            apply_bed = true(md.mesh.numberofvertices, 1);
        case 'grounded_only'
            % The survey-derived kriging field has no observational support
            % beneath the initial shelf. Retain the background there.
            apply_bed = initial_grounded;
        case 'grounded_plus_tapered_floating'
            % Retain the validated grounded kriging prior. Beneath floating
            % ice, perturb the background with a smaller independent field
            % that tapers to zero at the GL and cannot change flotation.
            apply_bed = initial_grounded;
            tapered_floating_bed = true;
        otherwise
            error(['[ICESEE] initial_bed_background_domain must be all, ', ...
                   'grounded_only, or grounded_plus_tapered_floating.']);
    end

    x = double(md.mesh.x(:));
    y = double(md.mesh.y(:));
    ice_mask = md.geometry.thickness(:) > 1.0;

    % Keep the initial bed challenge away from the grounding transition. This
    % prevents the bed prior itself from changing the mask before the first
    % analysis while retaining heterogeneous upstream bed errors.
    distance_to_gl = inf(size(x));
    distance_to_upstream_gl = inf(size(x));
    if bed_gl_buffer_m > 0 || tapered_floating_bed || ...
            gl_seaward_thickness_m > 0
        elements = double(md.mesh.elements);
        element_grounded = initial_grounded(elements);
        transition_elements = any(element_grounded, 2) & ...
                              any(~element_grounded, 2);
        gl_nodes = unique(elements(transition_elements, :));
        if ~isempty(gl_nodes)
            distance_to_gl = inf(size(x));
            for j = 1:numel(gl_nodes)
                node = gl_nodes(j);
                distance_to_gl = min(distance_to_gl, ...
                    hypot(x - x(node), y - y(node)));
            end
            if gl_seaward_thickness_m > 0
                % The truth GL is U-shaped. Euclidean distance to all of it
                % would make nearly the entire shelf close to a lateral arm.
                % Instead, find the leftmost (upstream) transition locally in
                % y and measure only the along-flow distance to that front.
                y_band = max(5000.0, 0.05 .* (max(y) - min(y)));
                for i = 1:numel(x)
                    local_gl = gl_nodes(abs(y(gl_nodes) - y(i)) <= y_band);
                    if isempty(local_gl)
                        [~, nearest] = min(abs(y(gl_nodes) - y(i)));
                        local_gl = gl_nodes(nearest);
                    end
                    upstream_gl_x = min(x(local_gl));
                    distance_to_upstream_gl(i) = abs(x(i) - upstream_gl_x);
                end
            end
            if bed_gl_buffer_m > 0 && ~tapered_floating_bed
                apply_bed = apply_bed & distance_to_gl >= bed_gl_buffer_m;
            end
        end
    end

    % Use independent, broad deterministic modes for bed and thickness.  The
    % modes are normalized on their application domains, so their configured
    % amplitudes are standard deviations and do not change the requested mean
    % biases.  They depend only on mesh coordinates and explicit configuration,
    % making the experiment reproducible without consulting the hidden truth.
    thickness_pattern = broad_thickness_pattern(x, y, initial_grounded);
    bed_pattern = spatial_prior_pattern( ...
        x, y, 1.25 .* pattern_length_x_m, 0.85 .* pattern_length_y_m, ...
        pattern_phase + 2.10, apply_bed);

    bed_delta = bed_offset_m + bed_anomaly_m .* bed_pattern;
    bed_delta = min(max(bed_delta, bed_delta_min_m), bed_delta_max_m);
    bed_new = bed_background;
    if tapered_floating_bed
        % Use the survey/kriging-supported upstream correction as the spatial
        % template. Stretch that same cross-flow structure continuously across
        % the domain rather than stitching or inventing a second realization.
        % Only its anomaly amplitude decreases smoothly downstream.
        seed_error = bed_candidate - bed_background;
        coherent_pattern = stretched_seed_bed_pattern( ...
            x, y, seed_error, ice_mask, initial_grounded, ...
            bed_seed_max_x_m);
        xn = (x - min(x(ice_mask))) ./ ...
             max(max(x(ice_mask)) - min(x(ice_mask)), eps);
        downstream_taper = min(max(xn, 0), 1);
        downstream_taper = downstream_taper .^ 2 .* ...
                           (3 - 2 .* downstream_taper);
        anomaly_envelope = 1 - ...
            (1 - bed_downstream_anomaly_factor) .* downstream_taper;
        coherent_delta = bed_offset_m + bed_anomaly_m .* ...
            anomaly_envelope .* coherent_pattern;
        % Smoothly approach the configured bounds rather than hard-clipping.
        % Hard clipping creates large constant-color plateaus in the prior.
        bound_center = 0.5 .* (bed_delta_min_m + bed_delta_max_m);
        bound_half_width = 0.5 .* ...
            (bed_delta_max_m - bed_delta_min_m);
        coherent_delta = bound_center + bound_half_width .* tanh( ...
            (coherent_delta - bound_center) ./ bound_half_width);

        % Use the same field on both sides of the GL. Its amplitude approaches
        % the configured floating factor continuously through the grounded
        % transition zone and stays at that moderate level beneath the shelf.
        domain_factor = floating_bed_anomaly_factor .* ones(size(x));
        grounded_apply = ice_mask & initial_grounded;
        if bed_gl_buffer_m > 0
            grounded_taper = min(max(distance_to_gl ./ bed_gl_buffer_m, 0), 1);
            grounded_taper(~isfinite(grounded_taper)) = 1;
            grounded_taper = grounded_taper .^ 2 .* ...
                              (3 - 2 .* grounded_taper);
        else
            grounded_taper = ones(size(x));
        end
        domain_factor(grounded_apply) = floating_bed_anomaly_factor + ...
            (1 - floating_bed_anomaly_factor) .* ...
            grounded_taper(grounded_apply);
        applied_delta = domain_factor .* coherent_delta;
        floating_apply = ice_mask & ~initial_grounded;
        applied_delta(floating_apply) = min(max( ...
            applied_delta(floating_apply), -floating_bed_max_error_m), ...
            floating_bed_max_error_m);
        bed_new(ice_mask) = bed_background(ice_mask) + ...
                            applied_delta(ice_mask);
    else
        bed_new(apply_bed) = bed_candidate(apply_bed) + bed_delta(apply_bed);
    end
    thickness_factor = thickness_scale + ...
                       thickness_anomaly_fraction .* thickness_pattern;
    thickness_factor = min(max(thickness_factor, thickness_factor_min), ...
                           thickness_factor_max);
    thickness_delta = thickness_anomaly_m .* thickness_pattern;
    thickness_delta = min(max(thickness_delta, thickness_delta_min_m), ...
                          thickness_delta_max_m);
    thickness_delta(~initial_grounded) = floating_thickness_anomaly_factor .* ...
                                         thickness_delta(~initial_grounded);
    thickness_new = max(thickness_factor .* md.geometry.thickness(:) + ...
                        thickness_delta, 1.0);

    di = md.materials.rho_ice / md.materials.rho_water;
    if tapered_floating_bed
        % Preserve the original floating topology without copying the hidden
        % bed for the ordinary prior. The margin is expressed in equivalent
        % ice-thickness metres. This safeguard intentionally precedes the
        % optional seaward-GL bump so that the robustness experiment can
        % physically ground a controlled strip of the original shelf.
        original_floating = ice_mask & ~initial_grounded;
        maximum_floating_bed = -di .* ...
            (thickness_new + floating_bed_flotation_margin_m);
        bed_new(original_floating) = min(bed_new(original_floating), ...
                                         maximum_floating_bed(original_floating));
    end
    if gl_seaward_thickness_m > 0
        % Controlled robustness experiment: start with a grounding line on
        % the seaward side of truth. Add a smooth, compact thickness bump
        % across both sides of the grounding zone. Including the grounded
        % side first repairs any local retreat caused by the heterogeneous
        % background prior; the same continuous bump can then ground a
        % controlled strip of shelf. This modifies the flotation function
        % physically; it does not shift a plotted contour or copy hidden bed.
        gl_taper = max(1 - distance_to_upstream_gl ./ ...
                       gl_seaward_width_m, 0);
        gl_taper(~isfinite(gl_taper)) = 0;
        gl_taper = gl_taper .^ 2 .* (3 - 2 .* gl_taper);
        thickness_new(ice_mask) = thickness_new(ice_mask) + ...
            gl_seaward_thickness_m .* gl_taper(ice_mask);
    end
    ocean_levelset = thickness_new + bed_new ./ di;
    grounded = ocean_levelset >= 0;
    floating = ~grounded;

    base_new = bed_new;
    base_new(floating) = -di .* thickness_new(floating);
    surface_new = base_new + thickness_new;

    md.geometry.bed = bed_new;
    md.geometry.base = base_new;
    md.geometry.thickness = thickness_new;
    md.geometry.surface = surface_new;
    md.mask.ocean_levelset = ocean_levelset;

    consistency_error = max(abs(md.geometry.surface - md.geometry.base - ...
                                md.geometry.thickness));
    if consistency_error > 1e-8
        error('[ICESEE] Failed to construct a consistent initial geometry.');
    end
    fprintf(['[ICESEE] Initial prior: mean H scale=%.3f, H fractional anomaly SD=%.3f, ' ...
             'H additive anomaly SD=%.1f m, floating factor=%.2f, ' ...
             'H delta range=[%.1f, %.1f] m, ' ...
             'H factor range=[%.3f, %.3f], seaward GL bump=%.1f m/%.1f km, ' ...
             'bed offset=%+.1f m, ' ...
             'bed anomaly SD=%.1f m, bed delta range=[%.1f, %.1f] m, ' ...
             'bed domain=%s, GL buffer=%.1f km, grounded=%d/%d\n'], ...
            thickness_scale, thickness_anomaly_fraction, thickness_anomaly_m, ...
            floating_thickness_anomaly_factor, ...
            min(thickness_delta(ice_mask)), max(thickness_delta(ice_mask)), ...
            min(thickness_factor(ice_mask)), max(thickness_factor(ice_mask)), ...
            gl_seaward_thickness_m, gl_seaward_width_m ./ 1000, ...
            bed_offset_m, bed_anomaly_m, min(bed_delta(apply_bed)), ...
            max(bed_delta(apply_bed)), bed_domain, bed_gl_buffer_m ./ 1000, ...
            nnz(grounded), ...
            md.mesh.numberofvertices);
end

function pattern = broad_thickness_pattern(x, y, mask)
%BROAD_THICKNESS_PATTERN Low-order mixed-sign geometry error.
% This field contains no repeating short-wavelength lobes.  It represents a
% domain-scale thickness tilt plus a gentle cross-flow component and is
% normalized over existing ice.
    xn = (x - min(x)) ./ max(max(x) - min(x), eps);
    yn = (y - min(y)) ./ max(max(y) - min(y), eps);
    raw = 0.70 .* cos(2 .* pi .* yn) + ...
          0.25 .* sin(2 .* pi .* xn) + ...
          0.15 .* cos(pi .* xn) .* sin(2 .* pi .* yn);
    mask = logical(mask(:));
    pattern = zeros(size(raw));
    if ~any(mask)
        return;
    end
    sigma = std(raw(mask));
    if ~isfinite(sigma) || sigma < 1.0e-12
        return;
    end
    pattern = (raw - mean(raw(mask))) ./ sigma;
end

function pattern = spatial_prior_pattern(x, y, length_x, length_y, phase, mask)
%SPATIAL_PRIOR_PATTERN Smooth reproducible sign-changing field on a mesh.
% A small spectral mixture avoids a single artificial stripe direction.  The
% output has zero mean and unit standard deviation over MASK.

    raw = sin(2 .* pi .* x ./ length_x + phase) .* ...
          cos(2 .* pi .* y ./ length_y - 0.61) + ...
          0.55 .* cos(2 .* pi .* x ./ (2.30 .* length_x) + ...
                     2 .* pi .* y ./ (1.70 .* length_y) + 1.37 + phase) + ...
          0.30 .* sin(2 .* pi .* x ./ (0.72 .* length_x) - ...
                     2 .* pi .* y ./ (2.40 .* length_y) + 0.83 - phase);
    mask = logical(mask(:));
    pattern = zeros(size(raw));
    if ~any(mask)
        return;
    end
    mu = mean(raw(mask));
    sigma = std(raw(mask));
    if ~isfinite(sigma) || sigma < 1.0e-12
        return;
    end
    pattern = (raw - mu) ./ sigma;
end

function pattern = stretched_seed_bed_pattern(x, y, seed_error, mask, ...
                                              norm_mask, seed_max_x)
%STRETCHED_SEED_BED_PATTERN Extend the supported upstream error structure.
% The upstream kriged-minus-background field is treated as a structural
% template, not as truth. Its x-coordinate is stretched continuously over the
% full ice domain, so there are no repeated tiles or stitched realizations.
    x = double(x(:));
    y = double(y(:));
    seed_error = double(seed_error(:));
    mask = logical(mask(:));
    norm_mask = logical(norm_mask(:)) & mask;
    pattern = zeros(size(x));
    if ~any(mask) || ~any(norm_mask)
        return;
    end

    source = mask & isfinite(seed_error) & x <= seed_max_x;
    if nnz(source) < 10
        error(['[ICESEE] Too few upstream vertices to construct the ', ...
               'stretched bed-error template.']);
    end

    source_x_min = min(x(source));
    source_x_max = max(x(source));
    domain_x_min = min(x(mask));
    domain_x_max = max(x(mask));
    x_query = source_x_min + (x - domain_x_min) .* ...
        (source_x_max - source_x_min) ./ ...
        max(domain_x_max - domain_x_min, eps);

    interpolant = scatteredInterpolant( ...
        x(source), y(source), seed_error(source), 'natural', 'nearest');
    raw = interpolant(x_query, y);

    mu = mean(raw(norm_mask));
    sigma = std(raw(norm_mask));
    if ~isfinite(sigma) || sigma < 1.0e-12
        return;
    end
    pattern(mask) = (raw(mask) - mu) ./ sigma;
end

function writeInitialStateHDF5(filename, md)
%WRITEINITIALSTATEHDF5 Save an unadvanced state using the transient schema.
    if isfile(filename)
        delete(filename);
    end
    names = {'Thickness_1', 'Surface_1', 'Vx_1', 'Vy_1', ...
             'bed_1', 'coefficient_1'};
    values = {md.geometry.thickness(:), md.geometry.surface(:), ...
              md.initialization.vx(:), md.initialization.vy(:), ...
              md.geometry.bed(:), md.friction.coefficient(:)};
    for i = 1:numel(names)
        h5create(filename, ['/' names{i}], size(values{i}));
        h5write(filename, ['/' names{i}], values{i});
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
