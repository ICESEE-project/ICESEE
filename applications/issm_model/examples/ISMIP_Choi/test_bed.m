% ============================================================
% Verify ICESEE-ISSM initialized geometry after interface writes
% ============================================================

clear; close all; clc;

data_path = '_modelrun_datasets';
ens_id = 0;
k = 1;

md = loadmodel(fullfile('data','ISMIP.Parameterization1.mat'));

% ---------------- True/nurged state ----------------
state_file = fullfile(data_path, 'true_nurged_states.h5');
true_state   = h5read(state_file,'/true_state')';
nurged_state = h5read(state_file,'/nurged_state')';

nvar = 6;
[nd, nt] = size(true_state);
hdim = nd / nvar;

Ih   = 1:hdim;
Is   = hdim+1:2*hdim;
Ibed = 4*hdim+1:5*hdim;

H_true    = true_state(Ih,k);
S_true    = true_state(Is,k);
bed_true  = true_state(Ibed,k);
base_true = S_true - H_true;

% ---------------- What Python wrote before ISSM ----------------
fb_file = fullfile(data_path, sprintf('friction_bed_%d.h5', ens_id));
bed_input = h5read(fb_file, '/bed');
coef_input = h5read(fb_file, '/coefficient');

bed_input = bed_input(:);

% ---------------- What ISSM interface wrote after initialization ----------------
ens_file = fullfile(data_path, sprintf('ensemble_out_%d.h5', ens_id));

H_init   = h5read(ens_file, '/Thickness');
S_init   = h5read(ens_file, '/Surface');
bed_init = h5read(ens_file, '/bed');

H_init   = H_init(:);
S_init   = S_init(:);
bed_init = bed_init(:);

base_raw = S_init - H_init;

% ---------------- Raw consistency checks ----------------
di = md.materials.rho_ice / md.materials.rho_water;
ocean_raw = H_init + bed_init / di;

pos_float_raw = find(ocean_raw < 0);
pos_grounded_raw = find(ocean_raw > 0);

raw_violations = sum(base_raw < bed_init);

% ---------------- Final projected geometry check ----------------
H_proj = H_init;
S_proj = S_init;
bed_proj = bed_init;
base_proj = S_proj - H_proj;

% Minimum thickness
pos = find(H_proj < 1);
H_proj(pos) = 1;
S_proj = base_proj + H_proj;

% Recompute ocean levelset
ocean_proj = H_proj + bed_proj / di;

% Floating ice: hydrostatic base
pos_float = find(ocean_proj < 0);
S_proj(pos_float) = H_proj(pos_float) .* ...
    (md.materials.rho_water - md.materials.rho_ice) / md.materials.rho_water;

base_proj = S_proj - H_proj;

% Base cannot be below bed
pos = find(base_proj < bed_proj);
base_proj(pos) = bed_proj(pos);

% Grounded ice: base equals bed
pos_grounded = find(ocean_proj > 0);
base_proj(pos_grounded) = bed_proj(pos_grounded);

% Final surface and final ocean levelset
S_proj = base_proj + H_proj;
ocean_proj = H_proj + bed_proj / di;

proj_violations = sum(base_proj < bed_proj);

% ---------------- Diagnostics ----------------
fprintf('\nICESEE-ISSM interface diagnostics\n');
fprintf('----------------------------------------\n');
fprintf('hdim                              = %d\n', hdim);
fprintf('bed_input length                  = %d\n', length(bed_input));
fprintf('bed_init length                   = %d\n', length(bed_init));

fprintf('\nRaw interface geometry:\n');
fprintf('RMSE input bed vs true            = %.3f m\n', sqrt(mean((bed_input-bed_true).^2)));
fprintf('RMSE initialized bed vs true      = %.3f m\n', sqrt(mean((bed_init-bed_true).^2)));
fprintf('RMSE raw base vs true             = %.3f m\n', sqrt(mean((base_raw-base_true).^2)));
fprintf('RMSE raw surface vs true          = %.3f m\n', sqrt(mean((S_init-S_true).^2)));
fprintf('Min raw thickness                 = %.3f m\n', min(H_init));
fprintf('Raw floating points               = %d\n', numel(pos_float_raw));
fprintf('Raw grounded points               = %d\n', numel(pos_grounded_raw));
fprintf('Raw base below bed violations     = %d\n', raw_violations);

fprintf('\nProjected geometry:\n');
fprintf('RMSE projected base vs true       = %.3f m\n', sqrt(mean((base_proj-base_true).^2)));
fprintf('RMSE projected surface vs true    = %.3f m\n', sqrt(mean((S_proj-S_true).^2)));
fprintf('Min projected thickness           = %.3f m\n', min(H_proj));
fprintf('Projected floating points         = %d\n', numel(pos_float));
fprintf('Projected grounded points         = %d\n', numel(pos_grounded));
fprintf('Projected base below bed          = %d\n', proj_violations);

fprintf('\nBed insertion check:\n');
fprintf('Max |bed_init - bed_input|        = %.6e m\n', max(abs(bed_init - bed_input)));

% ---------------- Plot ----------------
x = md.mesh.x / 1000;
y = md.mesh.y / 1000;
tri = md.mesh.elements;

figure('Color','w','Position',[100 100 1300 1200]);
tiledlayout(8,1,'TileSpacing','compact','Padding','compact');

nexttile
trisurf(tri,x,y,bed_true,'EdgeColor','none');
view(2); axis equal tight; colorbar; colormap(turbo);
title('True Bed');

nexttile
trisurf(tri,x,y,bed_input,'EdgeColor','none');
view(2); axis equal tight; colorbar;
title('Bed Written by Python: friction\_bed file');

nexttile
trisurf(tri,x,y,bed_init,'EdgeColor','none');
view(2); axis equal tight; colorbar;
title('Bed Returned by ISSM Interface: ensemble\_out');

nexttile
trisurf(tri,x,y,bed_init-bed_true,'EdgeColor','none');
view(2); axis equal tight; colorbar;
title('Initialized Bed Error: bed\_init - true');

nexttile
trisurf(tri,x,y,base_raw,'EdgeColor','none');
view(2); axis equal tight; colorbar;
title(sprintf('Raw Base, violations = %d', raw_violations));

nexttile
trisurf(tri,x,y,base_proj,'EdgeColor','none');
view(2); axis equal tight; colorbar;
title(sprintf('Projected Base, violations = %d', proj_violations));

nexttile
trisurf(tri,x,y,S_proj,'EdgeColor','none');
view(2); axis equal tight; colorbar;
title('Projected Surface');

nexttile
trisurf(tri,x,y,ocean_proj,'EdgeColor','none');
view(2); axis equal tight; colorbar;
title('Projected Ocean Levelset: >0 Grounded, <0 Floating');
xlabel('x (km)');
ylabel('y (km)');