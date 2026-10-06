function [ind_m, tm_m] = read_observation_metadata(file_path)
%READ_OBSERVATION_METADATA Read observation timing from synthetic_obs.h5.
% New files expose /obs_index and /obs_max_time.  Files generated before
% these aliases were introduced contain the equivalent /ind_m and /obs_t.

info = h5info(file_path);
dataset_names = {info.Datasets.Name};

if any(strcmp(dataset_names, 'obs_index'))
    ind_m = h5read(file_path, '/obs_index');
elseif any(strcmp(dataset_names, 'ind_m'))
    ind_m = h5read(file_path, '/ind_m');
else
    error('ICESEE:MissingObservationIndex', ...
        '%s contains neither /obs_index nor /ind_m.', file_path);
end

if any(strcmp(dataset_names, 'obs_max_time'))
    tm_m = h5read(file_path, '/obs_max_time');
elseif any(strcmp(dataset_names, 'obs_t'))
    obs_t = h5read(file_path, '/obs_t');
    if isempty(obs_t)
        error('ICESEE:EmptyObservationTimes', ...
            '%s contains an empty /obs_t dataset.', file_path);
    end
    tm_m = max(obs_t(:));
else
    error('ICESEE:MissingObservationTime', ...
        '%s contains neither /obs_max_time nor /obs_t.', file_path);
end
end
