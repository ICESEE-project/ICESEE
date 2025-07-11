issm_dir = getenv('ISSM_DIR');
if ~isempty(issm_dir)
    bin_path = fullfile(issm_dir, 'bin');
    lib_path = fullfile(issm_dir, 'lib');
    if ~contains(path, bin_path) || ~contains(path, lib_path)
        addpath(bin_path, lib_path);
        disp(['Added ISSM bin and lib directories from: ', issm_dir]);
    else
        disp(['ISSM bin and lib directories already in path: ', issm_dir]);
    end
else
    error('ISSM_DIR is not set.');
end



% issm_dir = getenv('ISSM_DIR');
% if ~isempty(issm_dir)
%     addpath(genpath(issm_dir));
%     disp(['Added ISSM directory and subdirectories from path: ', issm_dir]);
% else
%     error('ISSM_DIR is not set.');
% end