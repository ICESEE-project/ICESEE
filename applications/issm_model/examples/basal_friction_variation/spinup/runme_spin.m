Lx = 640000;
Ly = 80000;

% debug steps
steps=[1:3];

% run steps
% steps=[1:4]; 
% steps=[5:8];
% steps=[83];
% steps=[9]; 
% steps=[80:82];
% steps=[85];

% clear all; close all;

ens_id = 0;

folder = sprintf('./Models');
if ~exist(folder, 'dir')
    mkdir(folder);
end
% clear all; close all;

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
    % hvertices=[10000;500;500;10000];
    hvertices=[10000;500;5000;7500];
    md = bamg(md, 'domain', 'Domain.exp', 'hvertices',hvertices);
    % md = bamg(md, 'domain', 'Domain.exp', 'hmax', 2000, 'splitcorners', 1);
    % md = triangle(md,'Domain.exp',27000);
    
    plotmodel(md,'data','mesh');
    % Save mesh
    filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
    save(filename, 'md');
end

% Parameterization (Step 2)
if any(steps == 2)
    % filename = fullfile(folder, 'ISMIP.SetMask.mat');
    filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
    md = loadmodel(filename);
    md=setmask(md,'','');
    % md=setflowequation(md,'SSA','all');
    ParamFile = 'Mismip2_matlab.par';
    md = parameterize(md, ParamFile); % Use Mismip2.par
    
    filename = fullfile(folder, 'ISMIP.Parameterization.mat');
    save(filename, 'md');

    % write_netCDF(md, 'ISMIP_Parameterization.nc');
    % export_netCDF(md,'ISMIP_Parameterization.nc');
    plotmodel(md, 'data', md.geometry.bed, 'title', 'Ice bed_t=0'); 
    % plotmodel(md, 'data', md.initialization.vel, 'title', 'Ice velocity');
    % plotmodel(md, 'data', md.friction.coefficient, 'title', 'Ice friction');
end

% Adding bed roughness
if any(steps == 3)
    filename = fullfile(folder, 'ISMIP.Parameterization.mat');
    md = loadmodel(filename);
    
    % plotmodel(md, 'data', md.geometry.bed);

    % read in bed roughness data
    icesee_path='/Users/bkyanjo3/da_project/ICESEE/applications/issm_model/examples/ISMIP_Choi';
    filename = fullfile(icesee_path,'_modelrun_datasets/', 'friction_bed_0.h5');
    % filename = fullfile('friction_bed_0.h5');
    bed = h5read(filename, '/bed');
    nsize = md.mesh.numberofvertices;
    % md.geometry.bed = md.geometry.bed + bed(1:nsize) -  284 * rand(nsize,1);
    md.geometry.bed = md.geometry.bed + bed(1:nsize) -  284 * ones(nsize,1);
    % pos = (md.mesh.x > 375e3) & (md.mesh.x < 595e3) & (md.mesh.y > 65e3);
    % pos = (md.mesh.x < 640e3) & (md.mesh.y < 20e3);
    % md.geometry.bed(pos) = md.geometry.bed(pos) - 200;
    % md.geometry.base = md.geometry.base + bed(1:nsize);

    plotmodel(md, 'data', md.geometry.bed);
    
    filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
    save(filename, 'md');
end

% solve steady state
if any(steps == 4)
    filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
    md = loadmodel(filename);

    md=setflowequation(md,'SSA','all');

    % Time stepping
    md.timestepping=timesteppingadaptive();
    md.timestepping.time_step_max=100;
    md.timestepping.time_step_min=0.1;
    md.timestepping.final_time=5000;

    md.settings.output_frequency=100;
    md.stressbalance.maxiter=100;
    md.stressbalance.abstol = NaN;
    md.stressbalance.restol = 1;
    md.settings.solver_residue_threshold=5e-5;
    
    % Verbose
    md.verbose = verbose('all');

    md.cluster=generic('name',oshostname(),'np',4);
    md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
    md = solve(md, 'Transient','runtimename',false);

    filename = fullfile(folder, 'Transient_steadystate_0.mat');
    save(filename, 'md');
end

% solve steady state
if any(steps == 5)
    filename = fullfile(folder, 'Transient_steadystate_0.mat');
    md_ref = loadmodel(filename);

    filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
    md = loadmodel(filename);
    md=setflowequation(md,'SSA','all');
    
    md.geometry.thickness = md_ref.results.TransientSolution(end).Thickness;
    md.geometry.surface   = md_ref.results.TransientSolution(end).Surface;
    md.geometry.base      = md_ref.results.TransientSolution(end).Base;

    % Update other fields
    md.initialization.vx        = md_ref.results.TransientSolution(end).Vx;
    md.initialization.vy        = md_ref.results.TransientSolution(end).Vy;
    md.initialization.vel       = md_ref.results.TransientSolution(end).Vel;
    md.initialization.pressure  = md_ref.results.TransientSolution(end).Pressure;
    md.smb.mass_balance         = md_ref.results.TransientSolution(end).SmbMassBalance;
    md.mask.ocean_levelset      = md_ref.results.TransientSolution(end).MaskOceanLevelset;

    md.geometry.bed = md_ref.results.TransientSolution(end).Bed;
    md.friction.coefficient = md_ref.results.TransientSolution(end).FrictionCoefficient;
    
    % Time stepping
    md.timestepping=timesteppingadaptive();
    md.timestepping.time_step_max=100;
    md.timestepping.time_step_min=0.1;
    md.timestepping.final_time=5000;

    md.settings.output_frequency=100;
    md.stressbalance.maxiter=100;
    md.stressbalance.abstol = NaN;
    md.stressbalance.restol = 1;
    md.settings.solver_residue_threshold=5e-5;
    
    % Verbose
    md.verbose = verbose('all');

    md.cluster=generic('name',oshostname(),'np',4);
    md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
    md = solve(md, 'Transient','runtimename',false);

    filename = fullfile(folder, 'Transient_steadystate_00.mat');
    save(filename, 'md');
end

% Steady state simulation 2
if any(steps == 6)

    filename = fullfile(folder, 'Transient_steadystate_00.mat');
    md_ref = loadmodel(filename);

    filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
    md = loadmodel(filename);
    md=setflowequation(md,'SSA','all');
    
    md.geometry.thickness = md_ref.results.TransientSolution(end).Thickness;
    md.geometry.surface   = md_ref.results.TransientSolution(end).Surface;
    md.geometry.base      = md_ref.results.TransientSolution(end).Base;

    % Update other fields
    md.initialization.vx        = md_ref.results.TransientSolution(end).Vx;
    md.initialization.vy        = md_ref.results.TransientSolution(end).Vy;
    md.initialization.vel       = md_ref.results.TransientSolution(end).Vel;
    md.initialization.pressure  = md_ref.results.TransientSolution(end).Pressure;
    md.smb.mass_balance         = md_ref.results.TransientSolution(end).SmbMassBalance;
    md.mask.ocean_levelset      = md_ref.results.TransientSolution(end).MaskOceanLevelset;

    md.geometry.bed = md_ref.results.TransientSolution(end).Bed;
    md.friction.coefficient = md_ref.results.TransientSolution(end).FrictionCoefficient;
    
    md.timestepping=timesteppingadaptive();
    md.timestepping.final_time=7500;
    md.settings.output_frequency=100;
    md.stressbalance.maxiter=100;

    md.stressbalance.abstol = NaN;
    md.stressbalance.restol = 1;

    md.verbose = verbose('all');

    md.cluster=generic('name',oshostname(),'np',4);
    md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
    md = solve(md, 'Transient','runtimename',false);

    filename = fullfile(folder, 'Transient_steadystate_1.mat');
    save(filename, 'md');
end

% Steady state simulation 2
if any(steps == 7)

    filename = fullfile(folder, 'Transient_steadystate_1.mat');
    md_ref = loadmodel(filename);

    filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
    md = loadmodel(filename);
    md=setflowequation(md,'SSA','all');
    
    md.geometry.thickness = md_ref.results.TransientSolution(end).Thickness;
    md.geometry.surface   = md_ref.results.TransientSolution(end).Surface;
    md.geometry.base      = md_ref.results.TransientSolution(end).Base;

    % Update other fields
    md.initialization.vx        = md_ref.results.TransientSolution(end).Vx;
    md.initialization.vy        = md_ref.results.TransientSolution(end).Vy;
    md.initialization.vel       = md_ref.results.TransientSolution(end).Vel;
    md.initialization.pressure  = md_ref.results.TransientSolution(end).Pressure;
    md.smb.mass_balance         = md_ref.results.TransientSolution(end).SmbMassBalance;
    md.mask.ocean_levelset      = md_ref.results.TransientSolution(end).MaskOceanLevelset;

    md.geometry.bed = md_ref.results.TransientSolution(end).Bed;
    md.friction.coefficient = md_ref.results.TransientSolution(end).FrictionCoefficient;

    % filename = fullfile(folder, 'Transient_steadystate_11.mat');
    % save(filename, 'md');
    % md = loadmodel(filename);

    md.timestepping=timesteppingadaptive();
    md.timestepping.final_time=7500;
    md.settings.output_frequency=100;
    md.stressbalance.maxiter=100;

    md.stressbalance.abstol = NaN;
    md.stressbalance.restol = 1;

    md.verbose = verbose('all');

    md.cluster=generic('name',oshostname(),'np',4);
    md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
    md = solve(md, 'Transient','runtimename',false);

    filename = fullfile(folder, 'Transient_steadystate_2.mat');
    save(filename, 'md');
end
% Steady state simulation 2
if any(steps == 8)

    filename = fullfile(folder, 'Transient_steadystate_2.mat');
    md_ref = loadmodel(filename);

    filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
    md = loadmodel(filename);
    md=setflowequation(md,'SSA','all');
    
    md.geometry.thickness = md_ref.results.TransientSolution(end).Thickness;
    md.geometry.surface   = md_ref.results.TransientSolution(end).Surface;
    md.geometry.base      = md_ref.results.TransientSolution(end).Base;

    % Update other fields
    md.initialization.vx        = md_ref.results.TransientSolution(end).Vx;
    md.initialization.vy        = md_ref.results.TransientSolution(end).Vy;
    md.initialization.vel       = md_ref.results.TransientSolution(end).Vel;
    md.initialization.pressure  = md_ref.results.TransientSolution(end).Pressure;
    md.smb.mass_balance         = md_ref.results.TransientSolution(end).SmbMassBalance;
    md.mask.ocean_levelset      = md_ref.results.TransientSolution(end).MaskOceanLevelset;

    md.geometry.bed = md_ref.results.TransientSolution(end).Bed;
    md.friction.coefficient = md_ref.results.TransientSolution(end).FrictionCoefficient;

    % filename = fullfile(folder, 'Transient_steadystate_21.mat');
    % save(filename, 'md');
    % md = loadmodel(filename);

    md.timestepping=timestepping();
    md.timestepping.time_step=0.1;
    md.timestepping.start_time=0;
    md.timestepping.final_time=200;
    md.settings.output_frequency=10;
    md.stressbalance.maxiter=100;

    md.stressbalance.abstol = NaN;
    md.stressbalance.restol = 1;

    md.verbose = verbose('all');

    md.cluster=generic('name',oshostname(),'np',4);
    md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
    md = solve(md, 'Transient','runtimename',false);

    filename = fullfile(folder, 'reference_simulation.mat');
    save(filename, 'md');
end

if any(steps == 9)
    % clear all; close all;
    ens_id = 0;
    % folder = sprintf('./Models/ens_id_%d', ens_id);
    % if ~exist(folder, 'dir')
    %     mkdir(folder);
    % end
    % filename = fullfile(folder, 'Transient_steadystate1.mat');

    filename = fullfile(folder, 'reference_simulation.mat');
    md = loadmodel(filename);
    md_ref = loadmodel(filename);

    filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
    md = loadmodel(filename);
    md=setflowequation(md,'SSA','all');
    
    md.geometry.thickness = md_ref.results.TransientSolution(end).Thickness;
    md.geometry.surface   = md_ref.results.TransientSolution(end).Surface;
    md.geometry.base      = md_ref.results.TransientSolution(end).Base;

    % Update other fields
    md.initialization.vx        = md_ref.results.TransientSolution(end).Vx;
    md.initialization.vy        = md_ref.results.TransientSolution(end).Vy;
    md.initialization.vel       = md_ref.results.TransientSolution(end).Vel;
    md.initialization.pressure  = md_ref.results.TransientSolution(end).Pressure;
    md.smb.mass_balance         = md_ref.results.TransientSolution(end).SmbMassBalance;
    md.mask.ocean_levelset      = md_ref.results.TransientSolution(end).MaskOceanLevelset;

    md.geometry.bed = md_ref.results.TransientSolution(end).Bed;
    md.friction.coefficient = md_ref.results.TransientSolution(end).FrictionCoefficient;

    filename = fullfile(folder, 'ISMIP.reference_simulation_0.mat');
    save(filename, 'md','-v7.3');
    
    md.smb.mass_balance=-0.3*ones(md.mesh.numberofvertices,1);
    % md.smb.mass_balance=-1.5*ones(md.mesh.numberofvertices,1);
    md.transient.ismovingfront=0;
    % 
    md.basalforcings=linearbasalforcings();
    md.basalforcings.deepwater_melting_rate=150;
    md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
    
    % 
    md.timestepping=timestepping();
    md.timestepping.start_time = 0;
    md.timestepping.time_step = 0.2;
    md.timestepping.final_time = 200;
    md.settings.output_frequency = 10;
    md.stressbalance.maxiter = 100;
    md.stressbalance.restol = 1;
    md.stressbalance.reltol = 0.001;
    md.stressbalance.abstol = NaN;
    md.settings.solver_residue_threshold=5e-2;
    % 
    % md.verbose = verbose('all');
    % 
    % 
    md.transient.requested_outputs = {'default','FrictionCoefficient','Thickness','Base','Bed'};
    md.cluster=generic('name',oshostname(),'np',4);
    md.settings.waitonlock = 1;
    md = solve(md, 'Transient', 'runtimename',false);

    % % export_netCDF(md,'Transient_steadystate4.nc');
    % 
    filename = fullfile(folder, 'ISMIP.reference_simulation1.mat');
    save(filename, 'md', '-v7.3');

    % save file for ICESEE runs
    md_ref = loadmodel(filename);

    filename = fullfile(folder, 'ISMIP.Parameterization1.mat');
    md = loadmodel(filename);
    md=setflowequation(md,'SSA','all');
    
    md.geometry.thickness = md_ref.results.TransientSolution(end).Thickness;
    md.geometry.surface   = md_ref.results.TransientSolution(end).Surface;
    md.geometry.base      = md_ref.results.TransientSolution(end).Base;

    % Update other fields
    md.initialization.vx        = md_ref.results.TransientSolution(end).Vx;
    md.initialization.vy        = md_ref.results.TransientSolution(end).Vy;
    md.initialization.vel       = md_ref.results.TransientSolution(end).Vel;
    md.initialization.pressure  = md_ref.results.TransientSolution(end).Pressure;
    md.smb.mass_balance         = md_ref.results.TransientSolution(end).SmbMassBalance;
    md.mask.ocean_levelset      = md_ref.results.TransientSolution(end).MaskOceanLevelset;

    md.geometry.bed = md_ref.results.TransientSolution(end).Bed;
    md.friction.coefficient = md_ref.results.TransientSolution(end).FrictionCoefficient;


    filename = fullfile(folder, 'ISMIP.reference_simulation_1.mat');
    save(filename, 'md','-v7.3');

end