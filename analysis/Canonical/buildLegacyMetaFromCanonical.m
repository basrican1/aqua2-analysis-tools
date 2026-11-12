function [opts, DB] = buildLegacyMetaFromCanonical(filesTbl, fts1)
    % Reconstruct legacy-style opts/DB metadata from canonical files table.

    if istable(filesTbl)
        T = filesTbl;
    else
        error('filesTbl must be a table loaded from canonical files.mat');
    end

    if nargin < 2 || ~isstruct(fts1) || ~isfield(fts1, 'fileID')
        fileIDs = canonCellstr(T.file_id);
    else
        fileIDs = unique(canonCellstr(fts1.fileID), 'stable');
    end

    T_file_id = canonCellstr(T.file_id);
    n = numel(fileIDs);

    DB.fileID = cell(1,n);
    DB.animalID = cell(1,n);
    DB.sliceNum = cell(1,n);
    DB.nEvents = nan(1,n);
    DB.nCells = nan(1,n);
    DB.nSomas = nan(1,n);

    opts.sz = nan(n,4);
    opts.frameRate = nan(1,n);
    opts.spatialRes = nan(1,n);

    for ii = 1:n
        thisID = fileIDs{ii};
        hit = find(strcmp(T_file_id, thisID), 1, 'first');
        if isempty(hit)
            continue
        end

        DB.fileID{ii} = thisID;
        DB.animalID{ii} = asCharVal(T.animal_id(hit));
        DB.sliceNum{ii} = asCharVal(T.slice_num(hit));
        DB.nEvents(ii) = asNumVal(T.n_events(hit));
        DB.nCells(ii) = asNumVal(T.n_cells(hit));
        DB.nSomas(ii) = asNumVal(T.n_somas(hit));

        opts.sz(ii,1) = getNumCol(T, 'sz_x', hit, NaN);
        opts.sz(ii,2) = getNumCol(T, 'sz_y', hit, NaN);
        opts.sz(ii,3) = getNumCol(T, 'sz_z', hit, NaN);
        opts.sz(ii,4) = getNumCol(T, 'sz_t', hit, getNumCol(T, 'n_frames', hit, NaN));
        opts.frameRate(ii) = asNumVal(T.frame_rate_s(hit));
        opts.spatialRes(ii) = asNumVal(T.spatial_res(hit));
    end
end

function c = canonCellstr(v)
    if iscell(v)
        c = cellfun(@asCharVal, v, 'UniformOutput', false);
    elseif isstring(v)
        c = cellstr(v(:));
    elseif ischar(v)
        c = {v};
    else
        c = cellstr(string(v(:)));
    end
end

function c = asCharVal(v)
    if iscell(v), v = v{1}; end
    if ischar(v), c = v;
    elseif isstring(v), c = char(v);
    elseif isempty(v), c = '';
    else, c = char(string(v));
    end
end

function x = asNumVal(v)
    if iscell(v), v = v{1}; end
    if isnumeric(v), x = double(v(1));
    elseif islogical(v), x = double(v(1));
    else
        x = str2double(asCharVal(v));
    end
    if isnan(x), x = NaN; end
end

function x = getNumCol(T, colName, idx, defaultVal)
    if ismember(colName, T.Properties.VariableNames)
        x = asNumVal(T.(colName)(idx));
    else
        x = defaultVal;
    end
end
