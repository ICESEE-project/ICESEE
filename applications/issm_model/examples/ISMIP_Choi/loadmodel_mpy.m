function md = loadmodel_mpy(ncpath, varargin)
% LOADMODEL_MPY  Load ISSM .nc via Python 'loadmodel' and copy into MATLAB md=model.
% Usage:
%   md = loadmodel_mpy('MISMIP_2km_Experiment0.nc');
%   md = loadmodel_mpy('MISMIP_2km_Experiment0.nc','singletime',-1,'singleres','Vel');
%
% Notes:
%   - Requires the Python 'loadmodel' module available on sys.path (your working one).
%   - Copies common sub-objects and fields; robust to numpy arrays and masked arrays.
%   - Results are NOT reconstructed as ISSM classes here (can be added if needed).

    % ---------- 1) Load Python md ----------
    py_md = py_loadmodel(ncpath, varargin{:});

    % ---------- 2) Create MATLAB ISSM model ----------
    md = model;

    % ---------- 3) Sub-objects to copy (expand as needed) ----------
    subobjs = {'mesh','geometry','materials','friction','mask','constants', ...
               'timestepping','transient','flowequation','solver','cluster'};

    for i = 1:numel(subobjs)
        s = subobjs{i};
        if has_py_attr(py_md, s) && has_mat_prop(md, s)
            try
                copy_matching_props(md.(s), py_md.(s));
            catch ME
                fprintf(2,'Skip subobject %s: %s\n', s, ME.message);
            end
        end
    end

    % ---------- 4) A few top-level scalars if present ----------
    tops = {'deltat','dimension','numberofelements','numberofvertices'};
    for i = 1:numel(tops)
        p = tops{i};
        if has_py_attr(py_md,p)
            try, md.(p) = py2mat(py_md.(p)); end %#ok<TRYNC>
        end
    end

    % ---------- 5) Common arrays (explicit) ----------
    if has_py_attr(py_md,'mesh')
        try_assign(md.mesh, 'x',        py_md.mesh);
        try_assign(md.mesh, 'y',        py_md.mesh);
        try_assign(md.mesh, 'z',        py_md.mesh);
        try_assign(md.mesh, 'elements', py_md.mesh);
    end
    if has_py_attr(py_md,'geometry')
        try_assign(md.geometry, 'thickness', py_md.geometry);
        try_assign(md.geometry, 'bed',       py_md.geometry);
        try_assign(md.geometry, 'surface',   py_md.geometry);
        try_assign(md.geometry, 'base',      py_md.geometry);
    end

    % ---------- 6) Done ----------
end

% ======================== Helpers (inline) ========================

function py_md = py_loadmodel(ncpath, varargin)
    ensure_python_ready();
    lm = py.importlib.import_module('loadmodel');
    % forward only supported kwargs
    kv = {};
    valid = ["singletime","singleres"];
    for k = 1:2:numel(varargin)
        key = string(varargin{k});
        val = varargin{k+1};
        if any(valid == lower(key))
            kv(end+1:end+2) = {char(lower(key)), to_py_scalar(val)}; %#ok<AGROW>
        end
    end
    if isempty(kv)
        py_md = lm.loadmodel(ncpath);
    else
        py_md = lm.loadmodel(ncpath, pyargs(kv{:}));
    end
end

function ensure_python_ready()
    % Use your environment; OutOfProcess is safest for large objects
    try
        pe = pyenv;
        if pe.Status == "NotLoaded"
            pyenv("ExecutionMode","OutOfProcess");
        elseif pe.ExecutionMode ~= "OutOfProcess"
            pyenv("ExecutionMode","OutOfProcess");
        end
    catch
        % ok on older MATLAB
    end
    % Optionally add ISSM_DIR/{bin,lib}
    d = getenv('ISSM_DIR');
    if ~isempty(d)
        sys = py.importlib.import_module('sys');
        pbin = fullfile(d,'bin'); plib = fullfile(d,'lib');
        if isfolder(pbin) && ~any(strcmp(string(sys.path), string(pbin))), sys.path.append(pbin); end
        if isfolder(plib) && ~any(strcmp(string(sys.path), string(plib))), sys.path.append(plib); end
    end
    % Sanity check: something known
    py.importlib.import_module('issmversion');
end

function tf = has_py_attr(obj, name)
    try, tf = logical(py.hasattr(obj, name)); catch, tf = false; end
end

function tf = has_mat_prop(obj, name)
    try, tf = isprop(obj, name); catch, tf = false; end
end

function try_assign(tgt, field, src_pyobj)
    if has_py_attr(src_pyobj, field) && has_mat_prop(tgt, field)
        val = src_pyobj.(field);
        % Fill masked arrays with NaN if needed
        val = fill_masked(val);
        mval = py2mat(val);
        if ~isempty(mval)
            try, tgt.(field) = mval; catch, end
        end
    end
end

function copy_matching_props(tgt, src)
    props = string(properties(tgt));
    for i = 1:numel(props)
        p = char(props(i));
        if has_py_attr(src, p)
            v = src.(p);
            v = fill_masked(v);
            mv = py2mat(v);
            if ~isempty(mv)
                try, tgt.(p) = mv; catch, end
            end
        end
    end
end

function v = fill_masked(v)
% If v is numpy.ma.MaskedArray, fill masked entries with NaN
    try
        ma = py.importlib.import_module('numpy.ma');
        np = py.importlib.import_module('numpy');
        if py.isinstance(v, ma.MaskedArray)
            v = v.filled(np.nan);
        end
    catch
        % ignore if numpy not available
    end
end

function m = py2mat(x)
% Robust converter covering numpy ndarray, numpy.ma.MaskedArray, 0-d arrays,
% Python scalars/strings, list/tuple, dict, and None.

    % ---- NumPy / masked arrays / array-likes ----
    if has_py_attr(x,'dtype') || has_py_attr(x,'shape') || has_py_attr(x,'mask')
        m = np2mat(x);
        if ~isempty(m), return; end
    end

    % ---- scalars / strings ----
    if isa(x,'py.int') || isa(x,'py.long') || isa(x,'py.float')
        m = double(x); return
    end
    if isa(x,'py.bool'), m = logical(x); return; end
    if isa(x,'py.str'),  m = char(x);    return; end

    % ---- list / tuple ----
    if isa(x,'py.list') || isa(x,'py.tuple')
        C = cell(x);
        if ~isempty(C) && (isa(C{1},'py.int') || isa(C{1},'py.float') || has_py_attr(C{1},'dtype'))
            try, m = cellfun(@double,C); m = m(:); return; end
        end
        m = cell(size(C));
        for i = 1:numel(C), m{i} = py2mat(C{i}); end
        return
    end

    % ---- dict -> struct (shallow) ----
    if isa(x,'py.dict')
        keys = cell(x.keys());
        m = struct();
        for i = 1:numel(keys)
            k = keys{i};
            m.(matlab.lang.makeValidName(char(k))) = py2mat(x.get(k));
        end
        return
    end

    % ---- None / fallback ----
    if isa(x,'py.NoneType'), m = []; return; end
    try, m = char(py.str(x)); catch, m = []; end
end

function A = np2mat(x)
% Convert numpy ndarray / masked array to MATLAB double, regardless of ndim.
    A = [];
    try
        np = py.importlib.import_module('numpy');
        ma = py.importlib.import_module('numpy.ma');

        % Fill masks
        if py.isinstance(x, ma.MaskedArray)
            x = x.filled(np.nan);
        end

        % Ensure ndarray
        if ~py.isinstance(x, np.ndarray)
            x = np.asarray(x);
        end

        % 0-d arrays (numpy scalars)
        shp = cellfun(@double, cell(x.shape));
        if isempty(shp) || prod(shp)==0
            y = np.asarray(x).reshape(py.tuple({int32(1)}));
            A = double(y);
            A = A(1);
            return
        end

        % Convert to float64 for reliability
        try
            y = np.asarray(x, dtype=np.float64);
        catch
            y = np.asarray(x);
        end

        % Use tolist() to cross the MATLAB/Python boundary cleanly
        T = y.tolist();
        A = nested_pylist_to_double(T);

    catch
        % last resort direct conversion
        try, A = double(x); catch, A = []; end
    end
end

function M = nested_pylist_to_double(P)
% Recursively convert nested Python lists of numbers to MATLAB double array.
    if isa(P,'py.list') || isa(P,'py.tuple')
        % Determine shape by walking down the first branch
        shp = [];
        Q = P;
        while isa(Q,'py.list') || isa(Q,'py.tuple')
            Qc = cell(Q);
            shp(end+1) = numel(Qc); %#ok<AGROW>
            if isempty(Qc), break; end
            Q = Qc{1};
        end
        flat = flatten_num_list(P);
        M = reshape(double(cell2mat(flat)), fliplr(shp));
        M = permute(M, numel(shp):-1:1);  % row->col major
    else
        M = double(P);
    end
end

function out = flatten_num_list(P)
% Flatten nested py.list/py.tuple of numbers to a 1D cell array of doubles.
    out = {};
    if isa(P,'py.list') || isa(P,'py.tuple')
        C = cell(P);
        for i=1:numel(C)
            out = [out, flatten_num_list(C{i})]; %#ok<AGROW>
        end
    else
        out = {double(P)};
    end
end

function p = to_py_scalar(v)
    if isstring(v) || ischar(v)
        p = py.str(string(v));
    elseif islogical(v)
        p = py.bool(v);
    elseif isnumeric(v) && isscalar(v) && isfinite(v)
        p = py.int(int64(v)); % keep negative indices semantics
    else
        p = v;
    end
end