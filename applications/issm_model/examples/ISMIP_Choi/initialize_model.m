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

    folder = sprintf('./Models/ens_id_%d', ens_id);
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    % disp(['[MATLAB] Initializing model with rank: ', num2str(rank), ', nprocs: ', num2str(nprocs), ', ens_id: ', num2str(ens_id)]);

	steps = [1:5]; 
	
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
        hvertices=[10000;500;500;10000];
        md = bamg(md, 'domain', 'Domain.exp', 'hvertices',hvertices);
        % md = bamg(md, 'domain', 'Domain.exp', 'hmax', 10000, 'hmin', 500, 'splitcorners', 0);
        % md = triangle(md,'Domain.exp',27000);
        
        % plotmodel(md,'data','mesh');
        % Save mesh
        filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
        save(filename, 'md');
    end

    % Parameterization (Step 2)
    if any(steps == 2)
        filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
        md = loadmodel(filename);
        md = setflowequation(md, 'SSA', 'all'); % Shelfy-stream approximation
        ParamFile = 'Mismip2.par'
        md = parameterize(md, ParamFile); % Use Mismip2.par
        
        filename = fullfile(folder, 'ISMIP.Parameterization.mat');
        save(filename, 'md');

        % write_netCDF(md, 'ISMIP_Parameterization.nc');
        filename = fullfile(folder, 'ISMIP_Parameterization.nc');
        if ens_id == 0
            disp('[MATLAB] Exporting netCDF file for parameterization...');
            export_netCDF(md, filename);
        end
       
    end

    % Set boundary conditions (Step 5)
    if any(steps == 3)
        filename = fullfile(folder, 'ISMIP.Parameterization.mat');
        md = loadmodel(filename);
    
        % % Stressbalance referential
        md.stressbalance.referential = NaN(md.mesh.numberofvertices, 6);
       
        md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices, 1);
        md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices, 1);
      
        md.masstransport.spcthickness = NaN(md.mesh.numberofvertices, 1);
        
        filename = fullfile(folder, 'ISMIP.BoundaryCondition.mat');
        save(filename, 'md');
    end

    % Initialize ensemble fields (Step 6, equivalent to notebook's InitEnsemble)
    if any(steps == 4)
        filename = fullfile(folder, 'ISMIP.BoundaryCondition.mat');
        md = loadmodel(filename);

        fields = {'thickness', 'bed', 'coefficient'};

        result_0 = md.initialization(end);
        result_1 = md.geometry;
        result_2 = md.friction;

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
        data = {'Thickness', result_1, 'thickness';
                'bed', result_1, 'bed';
                'coefficient', result_2, 'coefficient'};
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

