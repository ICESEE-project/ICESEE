%% Extended spatial-diagnostic entry point
% Retains the full-grounded friction diagnostic and uses signed differences
% in physical units for bed, thickness, speed, and surface. A bright
% assimilated grounding line is reserved for estimated-friction panels; it
% is intentionally absent from truth and geometry/state panels. The main
% read_results.m publication view remains unchanged.

setenv('ICESEE_PLOT_ALL_GROUNDED_FRICTION','true');
setenv('ICESEE_RELATIVE_ERROR_MAPS','false');
setenv('ICESEE_OVERLAY_ASSIMILATED_GL','true');
read_results
