function initialize_model(rank, nprocs)

    %  read kwargs from a .mat file
	kwargs = load('model_kwargs.mat');

	%  access the values of the dictionary
	ParamFile 			 = char(kwargs.ParamFile);
	Lx 					 = double(kwargs.Lx); % length of the domain in x direction
	Ly 					 = double(kwargs.Ly); % length of the domain in y direction
	nx 					 = double(kwargs.nx); % number of nodes in x direction
	ny 					 = double(kwargs.ny); % number of nodes in y direction
    extrusion_layers     = double(kwargs.extrusion_layers); % number of layers for extrusion
    extrusion_exponent	 = double(kwargs.extrusion_exponent); % exponent for extrusion
	flow_model			 = char(kwargs.flow_model); % flow model to use
	sliding_vx			 = double(kwargs.sliding_vx); % sliding velocity in x
	sliding_vy			 = double(kwargs.sliding_vy); % sliding velocity in y
	cluster_name 		 = char(kwargs.cluster_name); % cluster name
	step_ens  		     = double(kwargs.steps); % step for ensemble 
	icesee_path		     = char(kwargs.icesee_path); % path to icesee
	data_path		     = char(kwargs.data_path); % path to data

	folder = sprintf('./Models/ens_id_%04d', rank);
	% Only create if it doesn't exist
	if ~exist(folder, 'dir')
		mkdir(folder);
	end

    steps = [1:6]; 

	disp(['[MATLAB] Running model with rank: ', num2str(rank), ', nprocs: ', num2str(nprocs)]);

    % Mesh generation #1
    if any(steps==1)
        %initialize md as a new model #help model
	    %->
	    md=model();
	    % generate a squaremesh #help squaremesh
		md = squaremesh(md,Lx,Ly,nx,ny);
	   
	    % save the given model
	    %->
		% filename = sprintf('./Models/ISMIP.Mesh_generation_%d.mat', rank);
		filename = fullfile(folder,'ISMIP.Mesh_generation.mat');
		save(filename, 'md');
	    % save ./Models/ISMIP.Mesh_generation md;
    end

    % Masks #2
    if any(steps==2)
	    % load the preceding step #help loadmodel
	    % path is given by the organizer with the name of the given step
	    %->
	    % md = loadmodel('./Models/ISMIP.Mesh_generation');
		% filename = sprintf('./Models/ISMIP.Mesh_generation_%d.mat', rank);
		filename = fullfile(folder,'ISMIP.Mesh_generation.mat');
		md = loadmodel(filename);

	    % set the mask #help setmask
	    % all MISMIP nodes are grounded
	    %->
	    md=setmask(md,'','');
	    % plot the given mask #md.mask to locate the field
	
	    % save the given model
	    %->
	    % save ./Models/ISMIP.SetMask md;
		% filename = sprintf('./Models/ISMIP.SetMask_%d.mat', rank);
		filename = fullfile(folder,'ISMIP.SetMask.mat');
		save(filename, 'md');
	  
    end
    
    %Parameterization #3
    if any(steps==3)
	    % load the preceding step #help loadmodel
	    % path is given by the organizer with the name of the given step
	    %->
	    % md = loadmodel('./Models/ISMIP.SetMask');
		% filename = sprintf('./Models/ISMIP.SetMask_%d.mat', rank);
		filename = fullfile(folder,'ISMIP.SetMask.mat');
		md = loadmodel(filename);

	    % parametrize the model # help parameterize
	    % you will need to fil-up the parameter file defined by the
	    % ParamFile variable
	    %->
	    md=parameterize(md,ParamFile);
	    % save the given model
	    %->
	    % save ./Models/ISMIP.Parameterization md;
		% filename = sprintf('./Models/ISMIP.Parameterization_%d.mat', rank);
		filename = fullfile(folder,'ISMIP.Parameterization.mat');
		save(filename, 'md');
	    % save ./Models/ISMIP.Parameterization md;
    end

    %Extrusion #4
    if any(steps==4)
	    
	    % load the preceding step #help loadmodel
	    % path is given by the organizer with the name of the given step
	    %->
	    % md = loadmodel('./Models/ISMIP.Parameterization');
		% filename = sprintf('./Models/ISMIP.Parameterization_%d.mat', rank);
		filename = fullfile(folder,'ISMIP.Parameterization.mat');
		md = loadmodel(filename);

	    % vertically extrude the preceding mesh #help extrude
	    % only 5 layers exponent 1
	    %->
	    md=extrude(md,extrusion_layers,extrusion_exponent);
	    % plot the 3D geometry #plotdoc
	    %->
	    % plotmodel(md,'data',md.geometry.base)
	    % save the given model
	    %->
	    % save ./Models/ISMIP.Extrusion md;
		% filename = sprintf('./Models/ISMIP.Extrusion_%d.mat', rank);
		filename = fullfile(folder,'ISMIP.Extrusion.mat');
		save(filename, 'md');
	    % save ./Models/ISMIP.Extrusion md;
    end

    %Set the flow computing method #5
    if any(steps==5)
    
	    % load the preceding step #help loadmodel
	    % path is given by the organizer with the name of the given step
	    %->
	    % md = loadmodel('./Models/ISMIP.Extrusion');
		% filename = sprintf('./Models/ISMIP.Extrusion_%d.mat', rank);
		filename = fullfile(folder,'ISMIP.Extrusion.mat');
		md = loadmodel(filename);
	    % set the approximation for the flow computation #help setflowequation
	    % We will be using the Higher Order Model (HO)
	    %->
	    md=setflowequation(md,flow_model,'all');
	    % save the given model
	    %->
	    % save ./Models/ISMIP.SetFlow md;
		% filename = sprintf('./Models/ISMIP.SetFlow_%d.mat', rank);
		filename = fullfile(folder,'ISMIP.SetFlow.mat');
		save(filename, 'md');
    end
    
    %Set Boundary Conditions #6
    if any(steps==6)
    
	    % load the preceding step #help loadmodel
	    % path is given by the organizer with the name of the given step
	    %->
	    % md = loadmodel('./Models/ISMIP.SetFlow');
		% filename = sprintf('./Models/ISMIP.SetFlow_%d.mat', rank);
		filename = fullfile(folder,'ISMIP.SetFlow.mat');
		md = loadmodel(filename);

	    % dirichlet boundary condition are known as SPCs
	    % ice frozen to the base, no velocity	#md.stressbalance
	    % SPCs are initialized at NaN one value per vertex
	    %->
	    md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
	    %->
	    md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
	    %->
	    md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
	    % extract the nodenumbers at the base #md.mesh.vertexonbase
	    %->
	    basalnodes=find(md.mesh.vertexonbase);
	    % set the sliding to zero on the bed
	    %->
	    md.stressbalance.spcvx(basalnodes)=sliding_vx;
	    %->
	    md.stressbalance.spcvy(basalnodes)=sliding_vy;
	    % periodic boundaries have to be fixed on the sides
	    % Find the indices of the sides of the domain, for x and then for y
	    % for x
	    % create maxX, list of indices where x is equal to max of x (use >> help find)
	    %->
	    maxX=find(md.mesh.x==max(md.mesh.x));
	    % create minX, list of indices where x is equal to min of x
	    %->
	    minX=find(md.mesh.x==min(md.mesh.x));
	    % for y
	    % create maxY, list of indices where y is equal to max of y
	    %  but not where x is equal to max or min of x
	    % (i.e, indices in maxX and minX should be excluded from maxY and minY)
	    %->
	    maxY=find(md.mesh.y==max(md.mesh.y) & md.mesh.x~=max(md.mesh.x) & md.mesh.x~=min(md.mesh.x));
	    % create minY, list of indices where y is equal to max of y
	    % but not where x is equal to max or min of x
	    %->
	    minY=find(md.mesh.y==min(md.mesh.y) & md.mesh.x~=max(md.mesh.x) & md.mesh.x~=min(md.mesh.x));
	    % set the node that should be paired together, minX with maxX and minY with maxY
	    % #md.stressbalance.vertex_pairing
	    %->
	    md.stressbalance.vertex_pairing=[minX,maxX;minY,maxY];
	    if (ParamFile=='IsmipF.par')
		    % if we are dealing with IsmipF the solution is in
		    % masstransport
		    md.masstransport.vertex_pairing=md.stressbalance.vertex_pairing;
	    end
	    % save the given model
	    %->
	    % save ./Models/ISMIP.BoundaryCondition md;
		% filename = sprintf('./Models/ISMIP.BoundaryCondition_%d.mat', rank);
		filename = fullfile(folder,'ISMIP.BoundaryCondition.mat');
		save(filename, 'md');
    end

	if step_ens == 8

	% 	% load Boundary conditions from the inital conditions
		% md = loadmodel('./Models/ISMIP.BoundaryCondition');
		% filename = sprintf('./Models/ISMIP.BoundaryCondition_%d.mat', rank);
		filename = fullfile(folder,'ISMIP.BoundaryCondition.mat');
		md = loadmodel(filename);
	
	% 	% save these fields to a file for ensemble use
		fields = {'vx', 'vy', 'vz', 'pressure'};
	% 	result = md.results.TransientSolution(end);
		result = md.initialization(end);
		
		% 	% --- fetch and save data for ensemble use
		filename = fullfile(icesee_path, data_path, sprintf('ensemble_init_%d.h5', rank));
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
		vx = result.vx;
		vy = result.vy;
		vz = result.vz;
		pressure = result.pressure;
		h5create(filename, '/Vx', size(vx));
		h5write(filename, '/Vx', vx);
		h5create(filename, '/Vy', size(vy));
		h5write(filename, '/Vy', vy);
		h5create(filename, '/Vz', size(vz));
		h5create(filename, '/Pressure', size(pressure));
		h5write(filename, '/Vz', vz);
		h5write(filename, '/Pressure', pressure);
	end
    
end

