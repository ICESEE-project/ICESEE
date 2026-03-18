
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

    deepwater_melting_rate = double(kwargs.deepwater_melting_rate);
    smb = double(kwargs.smb);
    mean_friction  = double(kwargs.mean_friction);

    reference_data = char(kwargs.reference_data);
    
    wrong_reference_data = 'wrong_reference_data.mat';

    % Get the current working directory
    cwd = pwd;
    [issmroot,~,~] = fileparts(fileparts(cwd));

    % set initail ens_id
    ens_id_init = 0;
    h_perturb = 250;
    s_perturb = 200;
    b_perturb = 100;
    nurged_entries_percentage = 0.25;

    output_frequency = 1; % make sure this is set to 1 for coupling with ICESEE

    % Set up model for each EnKF stage
    if strcmp(data_fname, 'true_state.mat')
        % Special case for true state
        % if k == 0 || isempty(k)
        folder = sprintf('./Models/ens_id_%d', ens_id_init);
        if ~exist(folder, 'dir')
            mkdir(folder);
        end
    
        % Initial run: load boundary conditions
        filename = fullfile(folder, reference_data);
        md = loadmodel(filename);

        % md = transientrestart(md);
        % update geometry
        md.geometry.thickness = md.results.TransientSolution(end).Thickness;
        md.geometry.surface   = md.results.TransientSolution(end).Surface;
        md.geometry.base      = md.results.TransientSolution(end).Base;

        % Update other fields
        md.initialization.vx        = md.results.TransientSolution(end).Vx;
        md.initialization.vy        = md.results.TransientSolution(end).Vy;
        md.initialization.vel       = md.results.TransientSolution(end).Vel;
        md.initialization.pressure  = md.results.TransientSolution(end).Pressure;
        % md.smb.mass_balance         = md.results.TransientSolution(end).SmbMassBalance;
        md.mask.ocean_levelset      = md.results.TransientSolution(end).MaskOceanLevelset;


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

        % Solve transient
        md = solve(md, 'Transient','runtimename',false);

        % update geometry
        % md.geometry.thickness = md.results.TransientSolution(end).Thickness;
        % md.geometry.surface   = md.results.TransientSolution(end).Surface;
        % md.geometry.base      = md.results.TransientSolution(end).Base;

        % % Update other fields
        % md.initialization.vx        = md.results.TransientSolution(end).Vx;
        % md.initialization.vy        = md.results.TransientSolution(end).Vy;
        % md.initialization.vel       = md.results.TransientSolution(end).Vel;
        % % md.initialization.pressure  = md.results.TransientSolution(end).Pressure;
        % % md.smb.mass_balance         = md.results.TransientSolution(end).SmbMassBalance;
        % md.mask.ocean_levelset      = md.results.TransientSolution(end).MaskOceanLevelset;

        % save updated model
        filename = fullfile(folder, data_fname);
        save(filename, 'md', '-v7.3');

        N = length(md.results.TransientSolution);
        % data = cell(N * 7, 7);   % 5 variables per step
        data = cell(N*6,6); 

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
        % filename_ref = fullfile(folder, 'true_state.mat');
        % md_ref = loadmodel(filename_ref);
        % filename = fullfile(folder, 'true_state.mat');
        % filename = fullfile(icesee_path, 'data', wrong_reference_data);
        md = loadmodel(filename);
        % md = transientrestart(md);

    %     % md = setflowequation(md,'SSA','all');

    %     % first_soln = 1; % start with the inital wrong state
    %    % update geometry
       md.geometry.thickness = md.results.TransientSolution(end).Thickness;
       md.geometry.surface   = md.results.TransientSolution(end).Surface;
       md.geometry.base      = md.results.TransientSolution(end).Base;

       % Update other fields
       md.initialization.vx        = md.results.TransientSolution(end).Vx;
       md.initialization.vy        = md.results.TransientSolution(end).Vy;
       md.initialization.vel       = md.results.TransientSolution(end).Vel;
       md.initialization.pressure  = md.results.TransientSolution(end).Pressure;
    %    md.smb.mass_balance         = md.results.TransientSolution(end).SmbMassBalance;
       md.mask.ocean_levelset      = md.results.TransientSolution(end).MaskOceanLevelset;
       md.geometry.bed              = md.results.TransientSolution(end).Bed;
       md.friction.coefficient      = md.results.TransientSolution(end).FrictionCoefficient;

        % setup nugged state
        friction_ref = mean_friction*ones(md.mesh.numberofvertices,1);
        % friction_ref = md.friction.coefficient;
        % friction_ref = md_ref.friction.coefficient;
        thickness_ref = md.geometry.thickness;
        bed_ref = md.geometry.bed;
        base_ref = md.geometry.base;

        % read the friction_bed file
        filename = fullfile(icesee_path, data_path, sprintf('friction_bed_%d.h5', ens_id));
        bed = h5read(filename, '/bed');
        coefficient = h5read(filename, '/coefficient');

        %  update the friction and bed
        md.friction.coefficient = friction_ref + coefficient;
        % md.friction.coefficient = coefficient;

        bed_err = bed - bed_ref;
        % md.geometry.bed = bed_ref + bed_err;
        % md.geometry.base = base_ref + bed_err;
        md.geometry.bed = bed_ref + bed - b_perturb*ones(md.mesh.numberofvertices,1);
        md.geometry.base = base_ref + bed - b_perturb*ones(md.mesh.numberofvertices,1);

        % md.geometry.bed = md.geometry.bed - b_perturb*ones(md.mesh.numberofvertices,1);
        % md.geometry.base = md.geometry.base - b_perturb*ones(md.mesh.numberofvertices,1);
        
        % Compute ice thickness
        % md.geometry.thickness = generate_correlated_field(md, thickness_ref, 10e3, 300);
        rand_field = h_perturb * (randn(md.mesh.numberofvertices,1) - min(randn(md.mesh.numberofvertices,1))) ./ (max(randn(md.mesh.numberofvertices,1)) - min(randn(md.mesh.numberofvertices,1)));
        % md.geometry.thickness = thickness_ref + rand_field;
        % md.geometry.thickness = thickness_ref + linspace(-h_perturb, 0, md.mesh.numberofvertices)';
        hdim = md.mesh.numberofvertices;
        h_indx = ceil(nurged_entries_percentage * hdim + 1);
        h_bump = linspace(-h_perturb, 0, h_indx)';
        h_with_bump = thickness_ref(1:h_indx) + h_bump;
        md.geometry.thickness = md.geometry.surface - md.geometry.base; 
        thickness_ref = md.geometry.thickness;
        md.geometry.thickness = [h_with_bump; thickness_ref(h_indx+1:end)];
        md.geometry.thickness = thickness_ref - h_perturb*ones(md.mesh.numberofvertices,1);
        
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

        % Solve transient
        md = solve(md, 'Transient','runtimename',false);
            
        filename = fullfile(folder, data_fname);
        save(filename, 'md', '-v7.3');

        % update geometry
        % md.geometry.thickness = md.results.TransientSolution(end).Thickness;
        % md.geometry.surface   = md.results.TransientSolution(end).Surface;
        % md.geometry.base      = md.results.TransientSolution(end).Base;

        % % Update other fields
        % md.initialization.vx        = md.results.TransientSolution(end).Vx;
        % md.initialization.vy        = md.results.TransientSolution(end).Vy;
        % md.initialization.vel       = md.results.TransientSolution(end).Vel;
        % % md.initialization.pressure  = md.results.TransientSolution(end).Pressure;
        % % md.smb.mass_balance         = md.results.TransientSolution(end).SmbMassBalance;
        % md.mask.ocean_levelset      = md.results.TransientSolution(end).MaskOceanLevelset;

        % save updated model
        % filename = fullfile(folder, data_fname);
        % save(filename, 'md', '-v7.3');

        N = length(md.results.TransientSolution);
        data = cell(N * 6, 6);   % 5 variables per step

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
            % Initial run: load boundary conditions
            % filename = fullfile(folder, reference_data);
            folder = sprintf('./Models/ens_id_%d', ens_id_init);
            % folder = sprintf('./Models/ens_id_%d', ens_id);
            if ~exist(folder, 'dir')
                mkdir(folder);
            end

            % filename_ref = fullfile(folder, reference_data);
            % filename_ref = fullfile(folder, 'true_state.mat');
            % md_ref = loadmodel(filename_ref);
            %*-----------------------
            % filename = fullfile(folder, 'true_state.mat');
            filename = fullfile(folder, reference_data);
            % filename = fullfile(icesee_path, 'data', wrong_reference_data);
            md = loadmodel(filename);
          

            md.geometry.thickness = md.results.TransientSolution(end).Thickness;
            md.geometry.surface   = md.results.TransientSolution(end).Surface;
            md.geometry.base      = md.results.TransientSolution(end).Base;

            % Update other fields
            md.initialization.vx        = md.results.TransientSolution(end).Vx;
            md.initialization.vy        = md.results.TransientSolution(end).Vy;
            md.initialization.vel       = md.results.TransientSolution(end).Vel;
            md.initialization.pressure  = md.results.TransientSolution(end).Pressure;
            % md.smb.mass_balance         = md.results.TransientSolution(end).SmbMassBalance;
            md.mask.ocean_levelset      = md.results.TransientSolution(end).MaskOceanLevelset;
            md.geometry.bed             = md.results.TransientSolution(end).Bed;
            md.friction.coefficient     = md.results.TransientSolution(end).FrictionCoefficient;

            friction_ref = mean_friction*ones(md.mesh.numberofvertices,1);
            % friction_ref = md.friction.coefficient;
            % friction_ref = md_ref.friction.coefficient;
            thickness_ref = md.geometry.thickness;
            bed_ref = md.geometry.bed;
            base_ref = md.geometry.base;

            % % read the friction_bed file
            filename = fullfile(icesee_path, data_path, sprintf('friction_bed_%d.h5', ens_id));
            bed = h5read(filename, '/bed');
            coefficient = h5read(filename, '/coefficient');

            %  update the friction and bed
            md.friction.coefficient = friction_ref + coefficient;
            % md.friction.coefficient = coefficient;
            bed_err = bed - bed_ref;
            % md.geometry.bed = bed_ref + bed_err;
            % md.geometry.base = base_ref + bed_err;
            md.geometry.bed = bed_ref + bed - b_perturb*ones(md.mesh.numberofvertices,1);
            md.geometry.base = base_ref + bed - b_perturb*ones(md.mesh.numberofvertices,1);
    
            md.smb.mass_balance=smb*ones(md.mesh.numberofvertices,1);

            % Compute ice thickness
            % md.geometry.thickness = md.geometry.surface - md.geometry.base;
            % md.geometry.thickness = generate_correlated_field(md, thickness_ref, 10e3, 300);
            % rand_field = h_perturb * (randn(md.mesh.numberofvertices,1) - min(randn(md.mesh.numberofvertices,1))) ./ (max(randn(md.mesh.numberofvertices,1)) - min(randn(md.mesh.numberofvertices,1)));
        % md.geometry.thickness = thickness_ref + rand_field;
            hdim = md.mesh.numberofvertices;
            h_indx = ceil(nurged_entries_percentage * hdim + 1);
            h_bump = linspace(-h_perturb, 0, h_indx)';
            h_with_bump = thickness_ref(1:h_indx) + h_bump;
            md.geometry.thickness = md.geometry.surface - md.geometry.base; 
            thickness_ref = md.geometry.thickness;
            md.geometry.thickness = [h_with_bump; thickness_ref(h_indx+1:end)];
            md.geometry.thickness = thickness_ref - h_perturb*ones(md.mesh.numberofvertices,1);
            
            % md.geometry.bed = md.geometry.bed - b_perturb*ones(md.mesh.numberofvertices,1);
            % md.geometry.base = md.geometry.base - b_perturb*ones(md.mesh.numberofvertices,1);

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

            % Grounded ice (ocean_levelset > 0)
            pos = find(md.mask.ocean_levelset > 0);
            md.geometry.base(pos) = md.geometry.bed(pos);
            md.geometry.surface = md.geometry.base + md.geometry.thickness;


            % --time stepping
            md.timestepping = timestepping();
            md.timestepping.time_step = 0.1;
            md.timestepping.start_time = 0;
            md.timestepping.final_time = 2.0;
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

            data = {'Thickness', result_0, 'Thickness';
                    % 'Base', result_0, 'Base';
                    'Surface', result_0, 'Surface';
                    'Vx', result_0, 'Vx';
                    'Vy', result_0, 'Vy';
                    'bed', result_0, 'Bed';
                    'coefficient', result_0, 'FrictionCoefficient'};

            writeToHDF5(filename, data);

            % write vx,vy, fcoeff, thickness, surface, bed, oceanlevelset, maskground to h5 file for each ensemble member
            filename_ens_init = fullfile(icesee_path, data_path, sprintf('ens_init%d.h5', ens_id));
            data = {'vx', result_0, 'Vx';
                    'vy', result_0, 'Vy';
                    'vel', result_0, 'Vel';
                    'fcoeff', result_0, 'FrictionCoefficient';
                    'thickness', result_0, 'Thickness';
                    'surface', result_0, 'Surface';
                    'bed', result_0, 'Bed';
                    'oceanlevelset', result_0, 'MaskOceanLevelset';
                    'base', result_0, 'Base';
                    %'maskground', result_0, 'MaskGroundedice'
                    };
            writeToHDF5(filename_ens_init, data);

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

            % update from initial ensemble file
            filename = fullfile(icesee_path, data_path, sprintf('ens_init%d.h5', ens_id));
            bed = h5read(filename, '/bed');
            coefficient = h5read(filename, '/fcoeff');
            thickness = h5read(filename, '/thickness');
            surface = h5read(filename, '/surface');
            vx = h5read(filename, '/vx');
            vy = h5read(filename, '/vy');
            vel = h5read(filename, '/vel');
            % base = h5read(filename, '/base');
            oceanlevelset = h5read(filename, '/oceanlevelset');
            md.initialization.vx   = vx;
            md.initialization.vy   = vy;
            md.initialization.vel  = sqrt(vx.^2 + vy.^2);
            md.geometry.thickness  = thickness;
            md.geometry.surface    = surface;
            md.geometry.bed        = bed;
            md.mask.ocean_levelset = oceanlevelset;
            md.friction.coefficient = coefficient;
            % md.geometry.base       = base;

            % --time stepping
            md.timestepping = timestepping();
            md.timestepping.time_step = 0.1;
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
            md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
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

            data = {'Thickness', result_0, 'Thickness';
                    % 'Base', result_0, 'Base';
                    'Surface', result_0, 'Surface';
                    'Vx', result_0, 'Vx';
                    'Vy', result_0, 'Vy';
                    'bed', result_1, 'Bed';
                    % 'coefficient', result_2, 'coefficient'};
                    'coefficient', result_2, 'FrictionCoefficient'};

            writeToHDF5(filename, data);

        else
          
            % fprintf('[MATLAB ---] Running model for ensemble ID %d, step %d\n', ens_id, k);
            
            % Subsequent time steps: 
            filename = fullfile(folder, data_fname);
            md = loadmodel(filename);
            
            % update geometry
            % md.geometry.thickness = md.results.TransientSolution(end).Thickness;
            % md.geometry.surface   = md.results.TransientSolution(end).Surface;
            md.geometry.base      = md.results.TransientSolution(end).Base;

            % Update other fields
            % md.initialization.vx        = md.results.TransientSolution(end).Vx;
            % md.initialization.vy        = md.results.TransientSolution(end).Vy;
            % md.initialization.vel       = md.results.TransientSolution(end).Vel;
            % md.initialization.pressure  = md.results.TransientSolution(end).Pressure;
            % md.smb.mass_balance         = md.results.TransientSolution(end).SmbMassBalance;
            md.mask.ocean_levelset      = md.results.TransientSolution(end).MaskOceanLevelset;

            % md.geometry.bed = md.results.TransientSolution(end).Bed;
            % md.friction.coefficient = md.results.TransientSolution(end).FrictionCoefficient;

            % Load ensemble input from HDF5
            filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
            md.geometry.surface = h5read(filename, '/Surface');
            % md.geometry.base = h5read(filename, '/Base');
            md.geometry.thickness = h5read(filename, '/Thickness');
            md.initialization.vx = h5read(filename, '/Vx');
            md.initialization.vy = h5read(filename, '/Vy');
            md.initialization.vel = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);
        
            % parameters for bed and friction
            md.geometry.bed = h5read(filename, '/bed');
            md.friction.coefficient = h5read(filename, '/coefficient');

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
            md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};

            % Solve transient
            md = solve(md, 'Transient','runtimename',false);

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

            % base data
            % data_base = {'Base', result_1, 'base'};
            % % write base to h5 file
            % filename = fullfile(icesee_path, data_path, sprintf('ensemble_base_%d.h5', ens_id));
            % writeToHDF5(filename, data_base);
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
