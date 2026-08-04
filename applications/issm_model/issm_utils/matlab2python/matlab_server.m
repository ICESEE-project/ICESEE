% ==================================================================
% @author: Brian Kyanjo
% @description: MATLAB server for executing commands from Python.
% @date: 2025-04-16
% ==================================================================

function matlab_server(cmdfile, statusfile)
    % cmdfile: File where commands are written by Python
    % statusfile: File to signal server status to Python
    
    disp('Server starting...');
    disp(['Command file: ', cmdfile]);
    disp(['Status file: ', statusfile]);
    
    % Publish statuses atomically.  Opening statusfile with 'w' briefly made
    % it visible as an empty file, so Python could read '' and kill a healthy
    % server before fprintf/fclose completed.
    write_status_atomic(statusfile, 'ready');
    disp(' Server initialized and ready.');
    
    % Main loop
    while true
        pause(0.1);  % Reduced pause for responsiveness
        if isfile(cmdfile)
            disp('Server Detected command file.');
            try
                % Read command
                fid = fopen(cmdfile, 'r');
                command = strtrim(fgetl(fid));
                fclose(fid);
                
                if isempty(command)
                    disp('Empty command, skipping.');
                    delete(cmdfile);
                    continue;
                end
                
                if strcmp(command, 'exit')
                    disp(' Received exit command.');
                    delete(cmdfile);
                    delete(statusfile);
                    disp('Server shutting down.');
                    break;
                end
                
                disp(['Executing: ', command]);
                evalin('base', command);  % Execute in base workspace
                disp('Command completed.');
                write_status_atomic(statusfile, 'done');
                drawnow;  % Ensure GUI updates if needed
                pause(0.1);  % Allow time for GUI updates
                disp('Status updated to done.');
                delete(cmdfile);  % Clean up
                
            catch ME
                disp(['[MATLAB ERROR] ', getReport(ME)]);
                write_status_atomic(statusfile, 'error');
                delete(cmdfile);  % Clean up even on error
            end
        end
    end
end

function write_status_atomic(statusfile, value)
    tmpfile = [statusfile, '.tmp'];
    fid = fopen(tmpfile, 'w');
    if fid < 0
        error('Could not open temporary status file: %s', tmpfile);
    end
    fprintf(fid, '%s', value);
    fclose(fid);
    [ok, message] = movefile(tmpfile, statusfile, 'f');
    if ~ok
        error('Could not publish MATLAB server status: %s', message);
    end
end
