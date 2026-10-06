% ============================================================
% Test kriged ensemble-mean bed before EnKF initialization
% ============================================================

clear; close all; clc;

data_path = '_modelrun_datasets';

% ---------------- Load ISSM model template ----------------
md = loadmodel(fullfile('data','ISMIP.Parameterization1.mat'));

% ---------------- Load true and nurged states ----------------
state_file = fullfile(data_path, 'true_nurged_states.h5');

true_state   = h5read(state_file,'/true_state')';
nurged_state = h5read(state_file,'/nurged_state')';

nvar = 6;
[nd, nt] = size(true_state);
hdim = nd / nvar;

k = 1;

% State indexing
Ih    = 1:hdim;
Is    = hdim+1:2*hdim;
Ivx   = 2*hdim+1:3*hdim;
Ivy   = 3*hdim+1:4*hdim;
Ibed  = 4*hdim+1:5*hdim;
Ifric = 5*hdim+1:6*hdim;

% ---------------- Reference/nurged state ----------------
H_ref    = nurged_state(Ih,k);
S_ref    = nurged_state(Is,k);
Vx_ref   = nurged_state(Ivx,k);
Vy_ref   = nurged_state(Ivy,k);
bed_ref  = nurged_state(Ibed,k);
fc_ref   = nurged_state(Ifric,k);

base_ref = S_ref - H_ref;

% ---------------- True state ----------------
H_true    = true_state(Ih,k);
S_true    = true_state(Is,k);
bed_true  = true_state(Ibed,k);
base_true = S_true - H_true;

% ---------------- Load kriged ensemble ----------------
krig_file = fullfile(data_path, 'bed_kriging_results.h5');
bed_ens = h5read(krig_file,'/bed_ens');

% Fix orientation if needed
if size(bed_ens,1) ~= hdim && size(bed_ens,2) == hdim
    bed_ens = bed_ens';
end

fprintf('bed_ens shape = %d x %d\n', size(bed_ens,1), size(bed_ens,2));
fprintf('hdim          = %d\n', hdim);

if size(bed_ens,1) ~= hdim
    error('bed_ens first dimension does not match hdim.');
end

% ============================================================
% USE ENSEMBLE MEAN BED
% ============================================================
bed_new = mean(bed_ens,2);

% ============================================================
% Apply Youngmin-style bed/base correction
% ============================================================

bed_err = bed_new - bed_ref;

% Move bed and base consistently
bed_test  = bed_ref  + bed_err;
base_test = base_ref + bed_err;

% Keep surface initially fixed
S_test = S_ref;

% Recompute thickness
H_test = S_test - base_test;

% Enforce minimum thickness
pos = find(H_test < 1);
H_test(pos) = 1;

% Update surface
S_test = base_test + H_test;

% ============================================================
% Hydrostatic consistency
% ============================================================

disp('Applying hydrostatic consistency ...');

di = md.materials.rho_ice / md.materials.rho_water;

% Ocean levelset
ocean_levelset = H_test + bed_test / di;

% Floating ice
pos_float = find(ocean_levelset < 0);

S_test(pos_float) = H_test(pos_float) .* ...
    (md.materials.rho_water - md.materials.rho_ice) ...
    / md.materials.rho_water;

base_test = S_test - H_test;

% Prevent base below bed
pos = find(base_test < bed_test);
base_test(pos) = bed_test(pos);

% Grounded ice
pos_grounded = find(ocean_levelset > 0);
base_test(pos_grounded) = bed_test(pos_grounded);

% Final surface update
S_test = base_test + H_test;

% ============================================================
% Diagnostics
% ============================================================

rmse_bed = sqrt(mean((bed_new - bed_true).^2));
rmse_base = sqrt(mean((base_test - base_true).^2));
rmse_surface = sqrt(mean((S_test - S_true).^2));

fprintf('\nDiagnostics:\n');
fprintf('-----------------------------------------\n');
fprintf('RMSE bed (ensemble mean vs true)   = %.3f m\n', rmse_bed);
fprintf('RMSE base (corrected vs true)      = %.3f m\n', rmse_base);
fprintf('RMSE surface (corrected vs true)   = %.3f m\n', rmse_surface);
fprintf('Min thickness                      = %.3f m\n', min(H_test));
fprintf('Max thickness                      = %.3f m\n', max(H_test));
fprintf('Floating points                    = %d\n', numel(pos_float));
fprintf('Grounded points                    = %d\n', numel(pos_grounded));
fprintf('Base below bed violations          = %d\n', sum(base_test < bed_test));

% ============================================================
% Coordinates
% ============================================================

x = md.mesh.x / 1000;
y = md.mesh.y / 1000;
tri = md.mesh.elements;

% ============================================================
% Plot
% ============================================================

figure('Color','w','Position',[100 100 1300 1100]);

tiledlayout(7,1,'TileSpacing','compact','Padding','compact');

% True bed
nexttile
trisurf(tri,x,y,bed_true,'EdgeColor','none');
view(2); axis equal tight;
colorbar; colormap(turbo);
title('True Bed');

% Nurged/reference bed
nexttile
trisurf(tri,x,y,bed_ref,'EdgeColor','none');
view(2); axis equal tight;
colorbar;
title('Reference/Nurged Bed');

% Ensemble mean bed
nexttile
trisurf(tri,x,y,bed_new,'EdgeColor','none');
view(2); axis equal tight;
colorbar;
title('Kriged Ensemble Mean Bed');

% Bed perturbation
nexttile
trisurf(tri,x,y,bed_err,'EdgeColor','none');
view(2); axis equal tight;
colorbar;
title('Bed Perturbation (bed_{mean} - bed_{ref})');

% Corrected base
nexttile
trisurf(tri,x,y,base_test,'EdgeColor','none');
view(2); axis equal tight;
colorbar;
title('Corrected Base');

% Corrected surface
nexttile
trisurf(tri,x,y,S_test,'EdgeColor','none');
view(2); axis equal tight;
colorbar;
title('Corrected Surface');

% Ocean levelset
nexttile
trisurf(tri,x,y,ocean_levelset,'EdgeColor','none');
view(2); axis equal tight;
colorbar;
title('Ocean Levelset (>0 grounded, <0 floating)');

xlabel('x (km)');
ylabel('y (km)');