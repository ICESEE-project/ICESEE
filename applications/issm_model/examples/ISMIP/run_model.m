function run_model(data_fname,ens_id,rank,nprocs,k,dt,tinitial,tfinal)
	% function run_model
	
		%  read kwargs from a .mat file
		model_kwargs = sprintf('model_kwargs_%d.mat', ens_id);
		kwargs 			= load(model_kwargs);
		cluster_name    = char(kwargs.cluster_name);
		steps 			= double(kwargs.steps);
		icesee_path     = char(kwargs.icesee_path);
		data_path       = char(kwargs.data_path);

		
		fprintf('[MATLAB] Running model with rank: %d, nprocs: %d filename: %s\n', rank, nprocs, data_fname);

		
		
		if any(steps==7)
			% load the preceding step #help loadmodel
			% path is given by the organizer with the name of the given step
			%->
			md = loadmodel('./Models/ISMIP.BoundaryCondition');
			% Set cluster #md.cluster
			% generic parameters #help generic
			% set only the name and number of process
			%->
			% cluster_name = 'cos2a16204.local'
			md.cluster=generic('name',cluster_name,'np',nprocs);
			% md.cluster=generic('name','cos2a16204.local','np',nprocs );
			% md.cluster=generic('name',oshostname(),'np',nprocs );
			% Set which control message you want to see #help verbose
			%->
			md.verbose=verbose('convergence',true);
			% Solve #help solve
			% we are solving a StressBalanc
	
			%  add time stepping parameters
			md.timestepping.time_step  = dt;
			md.timestepping.start_time = tinitial;
			md.timestepping.final_time = tfinal;
			%->
			md=solve(md,'Stressbalance');
			% save the given model
			%->
			save ./Models/ISMIP.StressBalance md;
			% plot the surface velocities #plotdoc
			%->
			% plotmodel(md,'data',md.results.StressbalanceSolution.Vel)
		end

		folder = sprintf('./Models/ens_id_%d', ens_id);
		% Only create if it doesn't exist
		if ~exist(folder, 'dir')
			mkdir(folder);
		end
		if any(steps==8)
			% load the preceding step #help loadmodel
			% path is given by the organizer with the name of the given step
			%->
			if k == 0 || isempty(k)
				% load Boundary conditions from the inital conditions
				% md = loadmodel('./Models/ISMIP.BoundaryCondition');
				% filename = sprintf('./Models/ISMIP.BoundaryCondition_%d.mat', rank);
				filename = fullfile(folder,'ISMIP.BoundaryCondition.mat');
				md = loadmodel(filename);

				% time stepping parameters
				md.timestepping.time_step=dt;
				md.timestepping.start_time=tinitial;
				md.timestepping.final_time=tfinal;

				md.cluster=generic('name',cluster_name,'np',nprocs);
				% Set which control message you want to see #help verbose
				%->
				md.verbose=verbose('convergence',true);
				% set the transient model to ignore the thermal model
				% #md.transient
				%->
				md.transient.isthermal=0;
				md.miscellaneous.name =  sprintf('color_%d', ens_id);
				md=solve(md,'Transient');
				% save the given model
				%->
				% save ./Models/ISMIP.Transient md;
				% filename_transient = sprintf('./Models/ISMIP.Transient_%d.mat', rank);
				% filename = fullfile(folder,'ISMIP.Transient.mat');
				filename = fullfile(folder, data_fname);	
				save(filename, 'md');

				% fprintf('[MATLAB] Running model at k: %d, with rank: %d, nprocs: %d filename: %s\n', k,rank, nprocs, data_fname);

				% save these fields to a file for ensemble use
				fields = {'Vx', 'Vy', 'Vz', 'Pressure'};
				result = md.results.TransientSolution(end);
				filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
				save_ensemble_hdf5(filename, result, fields);
				% disp['skipping the first step'];
			else
				
				% Load previous model
				% md = loadmodel('./Models/ISMIP.Transient');
				% filename = sprintf('./Models/ISMIP.Transient_%d.mat', rank);
				% filename = fullfile(folder,'ISMIP.Transient.mat');
				filename = fullfile(folder,data_fname);
				md = loadmodel(filename);

				% load from an ensemble_input file
				filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
				md.initialization.vx       = h5read(filename, '/Vx');
				md.initialization.vy       = h5read(filename, '/Vy');
				md.initialization.vz       = h5read(filename, '/Vz');
				md.initialization.pressure = h5read(filename, '/Pressure');

				% Time stepping parameters
				md.timestepping.time_step  = dt;
				md.timestepping.start_time = tinitial;
				md.timestepping.final_time = tfinal;

				% Cluster setup
				md.cluster = generic('name', cluster_name, 'np', nprocs);
				md.verbose = verbose('convergence', true);
				md.transient.isthermal = 0;

				% save ens_id as color to .mat file to be read inside solve function for each rank
				% color_name = sprintf('color_%d', ens_id);
				% fid = fopen(color_name, 'w');
				% fprintf(fid, '%d', ens_id);
				% fclose(fid);

				% now read the color name inside the solve function
				% md = loadmodel(color_name);

				md.miscellaneous.name =  sprintf('color_%d', ens_id);
				% Solve
				md = solve(md, 'Transient');

				% Save model
				% save('./Models/ISMIP.Transient', 'md');
				% filename_transient = sprintf('./Models/ISMIP.Transient_%d.mat', rank);
				% filename = fullfile(folder,'ISMIP.Transient.mat');
				filename_transient = fullfile(folder, data_fname);
				save(filename_transient, 'md');

				% save these fields to a file for ensemble use
				fields = {'Vx', 'Vy', 'Vz', 'Pressure'};
				result = md.results.TransientSolution(end);
				filename = fullfile(icesee_path, data_path, sprintf('ensemble_output_%d.h5', ens_id));
				save_ensemble_hdf5(filename, result, fields);
			end
	
		end
	end

	function save_ensemble_hdf5(filename, result, field_names)

		% Ensure the directory exists
		[filepath, ~, ~] = fileparts(filename);
		if ~exist(filepath, 'dir')
			mkdir(filepath);
		end
		
		% Remove file if it already exists
		if isfile(filename)
			delete(filename);
		end
	
		% Iterate over each requested field
		for i = 1:length(field_names)
			field = field_names{i};
	
			% Check field exists in result
			if isfield(result, field)
				data = result.(field);
				h5create(filename, ['/' field], size(data));
				h5write(filename, ['/' field], data);
			else
				warning('Field "%s" not found in result. Skipping.', field);
			end
		end
	
		% fprintf('[HDF5] Saved: %s\n', filename);
	end