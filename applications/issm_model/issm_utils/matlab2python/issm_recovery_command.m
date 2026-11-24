function issm_recovery_command(command, errorfile)
% Preserve the MATLAB exception identifier across the file-based bridge.
    try
        evalin('base', command);
    catch problem
        detail = struct('identifier', problem.identifier, 'message', problem.message);
        temporary = [errorfile, '.tmp'];
        fid = fopen(temporary, 'w');
        if fid < 0, rethrow(problem); end
        fprintf(fid, '%s', jsonencode(detail));
        fclose(fid);
        movefile(temporary, errorfile, 'f');
        rethrow(problem);
    end
end
