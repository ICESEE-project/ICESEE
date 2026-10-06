close all; clearvars; clear all

dt =0.2;

% file_paths='_modelrun_datasets_1';
file_paths='_modelrun_datasets';

scalar_means_file_true = fullfile(file_paths,'ensemble_true_state_scalar_0.h5');
scalar_means_file_nurged = fullfile(file_paths,'ensemble_nurged_state_scalar_0.h5');
scalar_means_file_ens = fullfile(file_paths,'ensemble_scalar_output.h5');
% scalar_means_file_ens = fullfile(file_paths,'ensemble_out_scalar_0.h5');

[ivaf_true,   t_true]   = read_ivaf_series(scalar_means_file_true,   dt);
[ivaf_nurged, t_nurged] = read_ivaf_series(scalar_means_file_nurged, dt);
[ivaf_ens,    t_ens]    = read_ivaf_series(scalar_means_file_ens,    dt);

plot(t_true, ivaf_true - ivaf_nurged(1), 'k-', 'LineWidth', 2.5); hold on
plot(t_nurged, ivaf_nurged - ivaf_nurged(1), 'm--', 'LineWidth', 2.5); hold on;
plot(t_ens, ivaf_ens - ivaf_nurged(1), 'c:', 'LineWidth', 2.5);
xlim([-1.5,188]);


function [ivaf, t_scalar] = read_ivaf_series(fname, dt)
    ivaf = [];
    t_scalar = [];

    if ~exist('fname','var') || isempty(fname) || ~isfile(fname)
        return;
    end

    info = h5info(fname);
    dnames = string({info.Datasets.Name});

    if any(dnames == "IceVolumeAboveFloatation")
        ivaf = h5read(fname, '/IceVolumeAboveFloatation');
        ivaf = ivaf(:);
    end

    % if any(dnames == "IceVolume")
    %     ivaf = h5read(fname, '/IceVolume');
    %     ivaf = ivaf(:);
    % end

    if any(dnames == "time")
        t_scalar = h5read(fname, '/time');
        t_scalar = t_scalar(:);
    elseif ~isempty(ivaf)
        t_scalar = (0:length(ivaf)-1)' * dt;
    end
end

