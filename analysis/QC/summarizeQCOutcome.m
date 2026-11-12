function summary = summarizeQCOutcome(original, filtered, qc_files, qc_cells, qc_events)
% Build compact QC summary tables for reporting.
%
% Inputs:
%   original: struct with canonical tables (files/cells/events/event_cell_map) before QC
%   filtered: struct with canonical tables after QC
%   qc_*: QC decision tables
%
% Output fields:
%   summary.overall
%   summary.by_file
%   summary.by_animal (if animal_id exists in files table)
%   summary.reasons.files / summary.reasons.cells / summary.reasons.events

summary = struct();

% Overall counts.
entities = {'files','cells','events','event_cell_map'};
n0 = zeros(numel(entities),1);
n1 = zeros(numel(entities),1);
for i = 1:numel(entities)
    e = entities{i};
    n0(i) = tableHeightOrZero(original, e);
    n1(i) = tableHeightOrZero(filtered, e);
end
summary.overall = table(string(entities(:)), n0, n1, n0 - n1, ...
    'VariableNames', {'entity','n_before','n_after','n_removed'});

% By file counts.
summary.by_file = table();
if isfield(original, 'files') && isfield(original, 'events') && istable(original.files) && istable(original.events) ...
        && ismember('file_id', original.files.Properties.VariableNames) && ismember('file_id', original.events.Properties.VariableNames)
    fileIDs = unique(string(original.files.file_id), 'stable');
    nFiles = numel(fileIDs);
    evBefore = zeros(nFiles,1);
    evAfter = zeros(nFiles,1);
    keepFile = true(nFiles,1);
    if istable(qc_files) && ~isempty(qc_files) && ismember('file_id', qc_files.Properties.VariableNames) && ismember('keep', qc_files.Properties.VariableNames)
        [tf, loc] = ismember(fileIDs, string(qc_files.file_id));
        keepFile(tf) = logical(qc_files.keep(loc(tf)));
    end
    for i = 1:nFiles
        fid = fileIDs(i);
        evBefore(i) = nnz(string(original.events.file_id) == fid);
        if isfield(filtered, 'events') && istable(filtered.events) && ismember('file_id', filtered.events.Properties.VariableNames)
            evAfter(i) = nnz(string(filtered.events.file_id) == fid);
        end
    end
    summary.by_file = table(fileIDs, keepFile, evBefore, evAfter, evBefore-evAfter, ...
        'VariableNames', {'file_id','keep_file','events_before','events_after','events_removed'});
end

% By animal counts (if available).
summary.by_animal = table();
if isfield(original, 'files') && istable(original.files) && ismember('animal_id', original.files.Properties.VariableNames)
    A0 = original.files(:, intersect({'animal_id'}, original.files.Properties.VariableNames, 'stable'));
    A0.n_files_before = ones(height(A0),1);
    A0g = groupsumOne(A0, 'animal_id', 'n_files_before');

    if isfield(filtered, 'files') && istable(filtered.files) && ismember('animal_id', filtered.files.Properties.VariableNames)
        A1 = filtered.files(:, intersect({'animal_id'}, filtered.files.Properties.VariableNames, 'stable'));
        A1.n_files_after = ones(height(A1),1);
        A1g = groupsumOne(A1, 'animal_id', 'n_files_after');
        summary.by_animal = outerjoin(A0g, A1g, 'Keys', 'animal_id', 'MergeKeys', true, 'Type', 'full');
    else
        summary.by_animal = A0g;
        summary.by_animal.n_files_after = zeros(height(summary.by_animal),1);
    end
    if ~ismember('n_files_before', summary.by_animal.Properties.VariableNames)
        summary.by_animal.n_files_before = zeros(height(summary.by_animal),1);
    end
    if ~ismember('n_files_after', summary.by_animal.Properties.VariableNames)
        summary.by_animal.n_files_after = zeros(height(summary.by_animal),1);
    end
    summary.by_animal.n_files_removed = summary.by_animal.n_files_before - summary.by_animal.n_files_after;
end

% Reason tallies.
summary.reasons = struct();
summary.reasons.files = reasonCounts(qc_files, 'file_id');
summary.reasons.cells = reasonCounts(qc_cells, 'cell_uid');
summary.reasons.events = reasonCounts(qc_events, 'event_uid');
end

function n = tableHeightOrZero(S, fieldName)
if isfield(S, fieldName) && istable(S.(fieldName))
    n = height(S.(fieldName));
else
    n = 0;
end
end

function out = reasonCounts(qc, keyName)
out = table();
if ~istable(qc) || isempty(qc) || ~ismember('keep', qc.Properties.VariableNames) || ~ismember('reason', qc.Properties.VariableNames) || ~ismember(keyName, qc.Properties.VariableNames)
    return
end
dropMask = ~logical(qc.keep);
if ~any(dropMask)
    return
end
R = table();
R.reason = string(qc.reason(dropMask));
R.reason(R.reason=="") = "<unspecified>";
R.n = ones(nnz(dropMask),1);
out = groupsumOne(R, 'reason', 'n');
out = sortrows(out, 'n', 'descend');
end

function out = groupsumOne(T, keyCol, valCol)
[G, keys] = findgroups(T.(keyCol));
vals = T.(valCol);
if isempty(G)
    out = table();
    return
end
s = accumarray(G, vals, [max(G),1], @sum, 0);
out = table(keys, s, 'VariableNames', {keyCol, valCol});
end

