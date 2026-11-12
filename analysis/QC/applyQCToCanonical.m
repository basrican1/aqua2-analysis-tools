function out = applyQCToCanonical(files, cells, events, event_cell_map, qc_files, qc_cells, qc_events, varargin)
% Apply file/cell/event QC decisions to canonical tables.
%
% Rules:
%   - Unspecified IDs default to keep=true.
%   - Effective keep requires both object keep flag and parent file keep.
%   - event_cell_map keeps rows where both event and cell survive.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'files', @istable);
addRequired(p, 'cells', @istable);
addRequired(p, 'events', @istable);
addRequired(p, 'event_cell_map', @istable);
addRequired(p, 'qc_files', @(x) istable(x) || isempty(x));
addRequired(p, 'qc_cells', @(x) istable(x) || isempty(x));
addRequired(p, 'qc_events', @(x) istable(x) || isempty(x));
addParameter(p, 'strictQC', false, @islogical);
parse(p, files, cells, events, event_cell_map, qc_files, qc_cells, qc_events, varargin{:});

mustHave(files, {'file_id'});
mustHave(cells, {'cell_uid','file_id'});
mustHave(events, {'event_uid','file_id'});
mustHave(event_cell_map, {'event_uid','cell_uid'});

report = validateQCTables(qc_files, qc_cells, qc_events, files, cells, events, 'strict', p.Results.strictQC);

fileKeep = idKeepMask(string(files.file_id), qc_files, 'file_id');
keepFileIDs = string(files.file_id(fileKeep));

cellKeep0 = idKeepMask(string(cells.cell_uid), qc_cells, 'cell_uid');
cellKeep = cellKeep0 & ismember(string(cells.file_id), keepFileIDs);
keepCellIDs = string(cells.cell_uid(cellKeep));

eventKeep0 = idKeepMask(string(events.event_uid), qc_events, 'event_uid');
eventKeep = eventKeep0 & ismember(string(events.file_id), keepFileIDs);
keepEventIDs = string(events.event_uid(eventKeep));

mapKeep = ismember(string(event_cell_map.event_uid), keepEventIDs) & ...
          ismember(string(event_cell_map.cell_uid), keepCellIDs);

out = struct();
out.files = files(fileKeep, :);
out.cells = cells(cellKeep, :);
out.events = events(eventKeep, :);
out.event_cell_map = event_cell_map(mapKeep, :);

out.masks = struct();
out.masks.fileKeep = fileKeep;
out.masks.cellKeep = cellKeep;
out.masks.eventKeep = eventKeep;
out.masks.mapKeep = mapKeep;

out.report = report;
end

function keep = idKeepMask(allIDs, qc, key)
keep = true(numel(allIDs), 1);
if isempty(qc)
    return
end
if ~ismember(key, qc.Properties.VariableNames) || ~ismember('keep', qc.Properties.VariableNames)
    return
end
qID = string(qc.(key));
qKeep = logical(qc.keep);
[tf, loc] = ismember(allIDs, qID);
keep(tf) = qKeep(loc(tf));
end

function mustHave(T, cols)
for i = 1:numel(cols)
    if ~ismember(cols{i}, T.Properties.VariableNames)
        error('applyQCToCanonical:MissingColumn', 'Missing required column: %s', cols{i});
    end
end
end

