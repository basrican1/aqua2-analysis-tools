function qc = markQC(qc, level, ids, keep, reason, reviewer)
% Batch-update QC decisions by ID.
%
% level: 'files' | 'cells' | 'events'
% ids  : cellstr/string list of file_id/cell_uid/event_uid

if nargin < 6 || isempty(reviewer), reviewer = ""; end
if nargin < 5 || isempty(reason), reason = ""; end

tblField = "";
idField = "";
switch lower(char(string(level)))
    case 'files'
        tblField = "qc_files"; idField = "file_id";
    case 'cells'
        tblField = "qc_cells"; idField = "cell_uid";
    case 'events'
        tblField = "qc_events"; idField = "event_uid";
    otherwise
        error('markQC:BadLevel', 'level must be files/cells/events');
end

if ~isstruct(qc) || ~isfield(qc, tblField)
    error('markQC:MissingTable', 'qc struct missing %s', char(tblField));
end
T = qc.(tblField);
if ~istable(T) || ~ismember(idField, T.Properties.VariableNames)
    error('markQC:BadTable', '%s must be a table with %s', char(tblField), char(idField));
end
if ~ismember('keep', T.Properties.VariableNames)
    T.keep = true(height(T),1);
end
if ~ismember('reason', T.Properties.VariableNames)
    T.reason = strings(height(T),1);
end
if ~ismember('reviewer', T.Properties.VariableNames)
    T.reviewer = strings(height(T),1);
end
if ~ismember('timestamp', T.Properties.VariableNames)
    T.timestamp = strings(height(T),1);
end

idList = string(ids(:));
allIDs = string(T.(idField));
[tf, loc] = ismember(idList, allIDs);
if any(~tf)
    warning('markQC:UnknownIDs', '%d IDs were not found in %s', nnz(~tf), char(tblField));
end

hitRows = unique(loc(tf));
if isempty(hitRows)
    qc.(tblField) = T;
    return
end

T.keep(hitRows) = logical(keep);
T.reason(hitRows) = string(reason);
T.reviewer(hitRows) = string(reviewer);
T.timestamp(hitRows) = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

qc.(tblField) = T;
end

