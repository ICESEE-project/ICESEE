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
    
    % Write "ready" to statusfile to signal Python that server is up
    fid = fopen(statusfile, 'w');
    fprintf(fid, 'ready');
    fclose(fid);
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
                delete(cmdfile);  % Clean up
                
            catch ME
                disp(['[MATLAB ERROR] ', getReport(ME)]);
                delete(cmdfile);  % Clean up even on error
            end
        end
    end
end