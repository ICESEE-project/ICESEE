function variable_size = initialize_model(rank, nprocs, ens_id)
    % Initialize ISSM model for MISMIP-like experiment
    % Inputs: rank, nprocs (MPI settings), ens_id (ensemble ID)
    % Output: variable_size (number of vertices for DART state vector)

    % Read kwargs from .mat file
    model_kwargs = sprintf('model_kwargs_%d.mat', ens_id);
    kwargs = load(model_kwargs);
    ParamFile = char(kwargs.ParamFile);
    Lx = double(kwargs.Lx); % 640000 m
    Ly = double(kwargs.Ly); % 80000 m
    cluster_name = char(kwargs.cluster_name);
    steps = double(kwargs.steps);
    icesee_path = char(kwargs.icesee_path);
    data_path = char(kwargs.data_path);

    folder = sprintf('./Models/ens_id_%d', ens_id);
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    disp(['[MATLAB] Initializing model with rank: ', num2str(rank), ', nprocs: ', num2str(nprocs), ', ens_id: ', num2str(ens_id)]);

	steps = [1:6]; 
	
    % Mesh generation (Step 1)
    if any(steps == 1)
        md = model();
        % Use bamg for variable-resolution mesh (500 m to 10 km)
        domain = [0, Lx, 0; Lx, Lx, Ly; Lx, 0, Ly; 0, 0, 0]; % Domain outline
        fid = fopen('Domain.exp', 'w');
        fprintf(fid, '0\n4\n');
        for i = 1:4
            fprintf(fid, '%f %f\n', domain(i, 2), domain(i, 3));
        end
        fclose(fid);
        md = bamg(md, 'domain', 'Domain.exp', 'hmax', 10000, 'hmin', 500, 'splitcorners', 0);
        % Save mesh
        filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
        save(filename, 'md');
    end

    % Masks (Step 2)
    if any(steps == 2)
        filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
        md = loadmodel(filename);
        md = setmask(md, '', ''); % All grounded, no ice shelves
        filename = fullfile(folder, 'ISMIP.SetMask.mat');
        save(filename, 'md');
    end

    % Parameterization (Step 3)
    if any(steps == 3)
        filename = fullfile(folder, 'ISMIP.SetMask.mat');
        md = loadmodel(filename);
        md = parameterize(md, ParamFile); % Use Mismip2.par
        % Initialize friction coefficient (Weertman law, reference value)
        md.friction = frictionweertman();
        md.friction.coefficient = 2500 * ones(md.mesh.numberofvertices, 1);
        md.friction.m = 3; % Weertman exponent
        filename = fullfile(folder, 'ISMIP.Parameterization.mat');
        save(filename, 'md');
    end

    % Set flow equation (Step 4, no extrusion)
    if any(steps == 4)
        filename = fullfile(folder, 'ISMIP.Parameterization.mat');
        md = loadmodel(filename);
        md = setflowequation(md, 'SSA', 'all'); % Shelfy-stream approximation
        filename = fullfile(folder, 'ISMIP.SetFlow.mat');
        save(filename, 'md');
    end

    % Set boundary conditions (Step 5)
    if any(steps == 5)
        filename = fullfile(folder, 'ISMIP.SetFlow.mat');
        md = loadmodel(filename);
        % Dirichlet boundary conditions (no-slip at base)
        md.stressbalance.spcvx = NaN * ones(md.mesh.numberofvertices, 1);
        md.stressbalance.spcvy = NaN * ones(md.mesh.numberofvertices, 1);
        basalnodes = find(md.mesh.vertexonbase);
        md.stressbalance.spcvx(basalnodes) = 0;
        md.stressbalance.spcvy(basalnodes) = 0;
        % Periodic boundaries on lateral sides
        maxX = find(md.mesh.x == max(md.mesh.x));
        minX = find(md.mesh.x == min(md.mesh.x));
        maxY = find(md.mesh.y == max(md.mesh.y) & md.mesh.x ~= max(md.mesh.x) & md.mesh.x ~= min(md.mesh.x));
        minY = find(md.mesh.y == min(md.mesh.y) & md.mesh.x ~= max(md.mesh.x) & md.mesh.x ~= min(md.mesh.x));
        md.stressbalance.vertex_pairing = [minX, maxX; minY, maxY];
        md.masstransport.vertex_pairing = md.stressbalance.vertex_pairing;
        filename = fullfile(folder, 'ISMIP.BoundaryCondition.mat');
        save(filename, 'md');
    end

    % Initialize ensemble fields (Step 6, equivalent to notebook's InitEnsemble)
    if any(steps == 6)
        filename = fullfile(folder, 'ISMIP.BoundaryCondition.mat');
        md = loadmodel(filename);
        % Load ensemble perturbations (from notebook Cell 12)
        ncfile = fullfile(icesee_path, data_path, sprintf('uncondition_fcoeff_err_ens1000.nc'));
        fcoeff = h5read(ncfile, '/fcoeff', [ens_id, 1], [1, md.mesh.numberofvertices]);
        md.friction.coefficient = 2500 + fcoeff';
        ncfile = fullfile(icesee_path, data_path, sprintf('condition_bed_err_30km_ens1000.nc'));
        bed_err = h5read(ncfile, '/bed_err', [ens_id, 1], [1, md.mesh.numberofvertices]);
        md.geometry.bed = md.geometry.bed + bed_err';
        md.geometry.base = md.geometry.base + bed_err';
        md.geometry.thickness = md.geometry.surface - md.geometry.base;
        pos = find(md.geometry.thickness < 1);
        md.geometry.thickness(pos) = 1;
        md.geometry.surface = md.geometry.base + md.geometry.thickness;
        % Hydrostatic adjustment (from notebook Cell 12)
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
        % Save ensemble initialization
        fields = {'vx', 'vy', 'thickness', 'surface', 'bed', 'friction.coefficient'};
        result = struct('vx', md.initialization.vx, 'vy', md.initialization.vy, ...
                        'thickness', md.geometry.thickness, 'surface', md.geometry.surface, ...
                        'bed', md.geometry.bed, 'friction_coefficient', md.friction.coefficient);
        filename = fullfile(icesee_path, data_path, sprintf('ensemble_init_%d.h5', ens_id));
        save_ensemble_hdf5(filename, result, fields);
        filename = fullfile(folder, 'ISMIP.EnsembleInit.mat');
        save(filename, 'md');
    end

    variable_size = md.mesh.numberofvertices;
end

function save_ensemble_hdf5(filename, result, field_names)
    % Save ensemble initialization to HDF5
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