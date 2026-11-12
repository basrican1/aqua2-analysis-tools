function BWA_AQ2 = attachEventFieldsFromCanonical(BWA_AQ2, eventsTbl, filesTbl)
    % Attach legacy-style event indexing fields from canonical events/files tables:
    %   BWA_AQ2.fts1.loc.t0, BWA_AQ2.fts1.loc.t1, BWA_AQ2.fts1.fileID

    hasT0 = isfield(BWA_AQ2, 'fts1') && isfield(BWA_AQ2.fts1, 'loc') && isfield(BWA_AQ2.fts1.loc, 't0');
    hasT1 = isfield(BWA_AQ2, 'fts1') && isfield(BWA_AQ2.fts1, 'loc') && isfield(BWA_AQ2.fts1.loc, 't1');
    hasFID = isfield(BWA_AQ2, 'fts1') && isfield(BWA_AQ2.fts1, 'fileID');
    if hasT0 && hasT1 && hasFID
        return
    end

    if ~istable(eventsTbl)
        error('eventsTbl must be canonical events table.');
    end
    if nargin < 3
        filesTbl = table();
    end

    E = eventsTbl;
    if ~ismember('file_id', E.Properties.VariableNames)
        error('eventsTbl missing file_id column.');
    end
    if ~ismember('local_event_id', E.Properties.VariableNames)
        error('eventsTbl missing local_event_id column.');
    end
    if ~ismember('t0_frame', E.Properties.VariableNames) || ~ismember('t1_frame', E.Properties.VariableNames)
        error('eventsTbl missing t0_frame/t1_frame columns.');
    end

    % Preserve original file ordering from canonical files table when available.
    if istable(filesTbl) && ismember('file_id', filesTbl.Properties.VariableNames) && ~isempty(filesTbl)
        fileOrder = canonCellstr(filesTbl.file_id);
        eFile = canonCellstr(E.file_id);
        fileRank = nan(height(E),1);
        for ii = 1:numel(fileOrder)
            fileRank(strcmp(eFile, fileOrder{ii})) = ii;
        end
        known = fileRank(~isnan(fileRank));
        if isempty(known)
            fileRank(isnan(fileRank)) = 1;
        else
            fileRank(isnan(fileRank)) = max(known) + 1;
        end
        E.file_rank = fileRank;
        E = sortrows(E, {'file_rank', 'local_event_id'});
    else
        E = sortrows(E, {'file_id', 'local_event_id'});
    end

    BWA_AQ2.fts1.loc.t0 = rowNum(E.t0_frame);
    BWA_AQ2.fts1.loc.t1 = rowNum(E.t1_frame);
    BWA_AQ2.fts1.fileID = rowCellstr(E.file_id);
end

function x = rowNum(v)
    if iscell(v)
        x = cellfun(@(z) double(z(1)), v);
    else
        x = double(v(:));
    end
    x = x(:)';
end

function c = rowCellstr(v)
    c = canonCellstr(v);
    c = c(:)';
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
    if ischar(v), c=v;
    elseif isstring(v), c=char(v);
    elseif isempty(v), c='';
    else, c=char(string(v));
    end
end
