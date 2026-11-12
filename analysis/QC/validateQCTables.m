function report = validateQCTables(qc_files, qc_cells, qc_events, files, cells, events, varargin)
% Validate QC tables against canonical IDs.
%
% report fields:
%   ok
%   errors
%   warnings

p = inputParser;
p.StructExpand = false;
addRequired(p, 'qc_files', @(x) istable(x) || isempty(x));
addRequired(p, 'qc_cells', @(x) istable(x) || isempty(x));
addRequired(p, 'qc_events', @(x) istable(x) || isempty(x));
addRequired(p, 'files', @istable);
addRequired(p, 'cells', @istable);
addRequired(p, 'events', @istable);
addParameter(p, 'strict', false, @islogical);
parse(p, qc_files, qc_cells, qc_events, files, cells, events, varargin{:});

report = struct();
report.ok = true;
report.errors = {};
report.warnings = {};

report = checkOne(report, qc_files, files, 'file_id', 'qc_files');
report = checkOne(report, qc_cells, cells, 'cell_uid', 'qc_cells');
report = checkOne(report, qc_events, events, 'event_uid', 'qc_events');

if p.Results.strict && ~report.ok
    error('validateQCTables:Invalid', '%s', strjoin(report.errors, ' | '));
end
end

function report = checkOne(report, qc, src, key, name)
if isempty(qc)
    return
end

req = {key, 'keep'};
for i = 1:numel(req)
    if ~ismember(req{i}, qc.Properties.VariableNames)
        report.ok = false;
        report.errors{end+1} = sprintf('%s missing required column %s', name, req{i}); %#ok<AGROW>
    end
end
if ~report.ok
    return
end

id = string(qc.(key));
if numel(unique(id)) ~= numel(id)
    report.ok = false;
    report.errors{end+1} = sprintf('%s has duplicate %s values', name, key); %#ok<AGROW>
end

keep = qc.keep;
if ~(islogical(keep) || isnumeric(keep))
    report.ok = false;
    report.errors{end+1} = sprintf('%s.keep must be logical or numeric', name); %#ok<AGROW>
end

srcID = string(src.(key));
unknown = ~ismember(id, srcID);
if any(unknown)
    report.warnings{end+1} = sprintf('%s has %d unknown IDs not in canonical %s', name, nnz(unknown), key); %#ok<AGROW>
end
end

