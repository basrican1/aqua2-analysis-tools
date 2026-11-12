function qc = createQCTables(files, cells, events, varargin)
% Create QC decision tables for canonical data.
%
% Outputs:
%   qc.qc_files  : file_id, keep, reason, reviewer, timestamp
%   qc.qc_cells  : cell_uid, keep, reason, reviewer, timestamp
%   qc.qc_events : event_uid, keep, reason, reviewer, timestamp

p = inputParser;
p.StructExpand = false;
addRequired(p, 'files', @istable);
addRequired(p, 'cells', @istable);
addRequired(p, 'events', @istable);
addParameter(p, 'defaultKeep', true, @islogical);
addParameter(p, 'reviewer', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'reason', "", @(x) ischar(x) || isstring(x));
parse(p, files, cells, events, varargin{:});

requireCol(files, 'file_id');
requireCol(cells, 'cell_uid');
requireCol(events, 'event_uid');

defaultKeep = p.Results.defaultKeep;
reviewer = string(p.Results.reviewer);
reason = string(p.Results.reason);
ts = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

qc_files = table();
qc_files.file_id = string(files.file_id);
qc_files.keep = repmat(defaultKeep, height(files), 1);
qc_files.reason = repmat(reason, height(files), 1);
qc_files.reviewer = repmat(reviewer, height(files), 1);
qc_files.timestamp = repmat(ts, height(files), 1);

qc_cells = table();
qc_cells.cell_uid = string(cells.cell_uid);
qc_cells.keep = repmat(defaultKeep, height(cells), 1);
qc_cells.reason = repmat(reason, height(cells), 1);
qc_cells.reviewer = repmat(reviewer, height(cells), 1);
qc_cells.timestamp = repmat(ts, height(cells), 1);

qc_events = table();
qc_events.event_uid = string(events.event_uid);
qc_events.keep = repmat(defaultKeep, height(events), 1);
qc_events.reason = repmat(reason, height(events), 1);
qc_events.reviewer = repmat(reviewer, height(events), 1);
qc_events.timestamp = repmat(ts, height(events), 1);

qc = struct('qc_files', qc_files, 'qc_cells', qc_cells, 'qc_events', qc_events);
end

function requireCol(T, c)
if ~ismember(c, T.Properties.VariableNames)
    error('createQCTables:MissingColumn', 'Missing required column: %s', c);
end
end

