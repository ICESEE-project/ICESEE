function issm_recovery_rng(action, filename)
% A function scope avoids leaving transaction variables in the model workspace.
    if strcmp(action, 'save')
        state = rng;
        save(filename, 'state');
    elseif strcmp(action, 'restore')
        saved = load(filename, 'state');
        rng(saved.state);
    else
        error('ICESEE:ISSM:RecoveryProtocol', 'Unknown RNG transaction action');
    end
end
