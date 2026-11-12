function out = computeNEventsPerMinCanonical(events, event_cell_map, files, epochSpec, varargin)
% Compute event frequency (events/min) from canonical tables, no plotting.
%
% Inputs:
%   events: canonical events table
%   event_cell_map: canonical mapping table
%   files: canonical files table
%   epochSpec: struct from frameTimeConversionsByEpoch()
%
% Options:
%   'areaPercentiles' [min max], default [0 100]
%   'groupBy' cellstr of grouping columns from files table (default {'file_id'})
%   'primaryOnly' logical, default true
%   'dropZeroRates' logical, default true
%
% Output fields:
%   out.ratesByCellEpoch
%   out.summaryByGroupEpoch
%   out.eventsAssigned
%   out.epochTable

p = inputParser;
p.StructExpand = false;
addRequired(p, 'events', @istable);
addRequired(p, 'event_cell_map', @istable);
addRequired(p, 'files', @istable);
addRequired(p, 'epochSpec', @isstruct);
addParameter(p, 'areaPercentiles', [0 100], @isnumeric);
addParameter(p, 'groupBy', {'file_id'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'primaryOnly', true, @islogical);
addParameter(p, 'dropZeroRates', true, @islogical);
addParameter(p, 'cells', table(), @(x) istable(x) || isempty(x));
parse(p, events, event_cell_map, files, epochSpec, varargin{:});

if ~ismember('event_uid', events.Properties.VariableNames)
    error('events must include event_uid');
end
if ~ismember('file_id', events.Properties.VariableNames)
    error('events must include file_id');
end

grpCols = cellstr(string(p.Results.groupBy(:)));
grpCols = grpCols(:)';

% 1) Filter events by canonical raw area percentile.
E = events;
nInput = height(E);
if ismember('area_px2', E.Properties.VariableNames)
    ap = p.Results.areaPercentiles;
    lo = prctile(E.area_px2, ap(1));
    hi = prctile(E.area_px2, ap(2));
    E = E(E.area_px2 >= lo & E.area_px2 <= hi, :);
end

% 2) Attach epoch.
A = assignEpochsToEvents(E, epochSpec, 'method', 'auto');
E = A.events;
nExcludedNonEpoch = sum(E.epoch_idx <= 0);
E = E(E.epoch_idx > 0, :);
nAssignedEpoch = height(E);

% 3) Select one cell per event (primary default).
M = event_cell_map;
if p.Results.primaryOnly && ismember('is_primary_cell', M.Properties.VariableNames)
    M = M(M.is_primary_cell, :);
end
if isempty(M)
    error('No event->cell mappings available after primaryOnly filter.');
end

ME = innerjoin(M, E(:, {'event_uid','file_id','epoch_idx','epoch_name'}), 'Keys', 'event_uid');
if ~ismember('cell_uid', ME.Properties.VariableNames)
    error('event_cell_map must include cell_uid');
end

% 4) Join file metadata for requested groups.
needCols = unique([{'file_id'}, grpCols], 'stable');
needCols = needCols(ismember(needCols, files.Properties.VariableNames));
if ~ismember('file_id', needCols)
    needCols = [{'file_id'}, needCols];
end
T = outerjoin(ME, files(:, needCols), 'Keys', 'file_id', 'Type', 'left', 'MergeKeys', true);

% 5) Build cell-epoch rate table.
for i = 1:numel(grpCols)
    c = grpCols{i};
    if ~ismember(c, T.Properties.VariableNames)
        T.(c) = repmat("", height(T), 1);
    end
end

keyCols = unique([{'file_id'}, grpCols, {'cell_uid','epoch_idx','epoch_name'}], 'stable');
[G1, grp1] = groupTableRows(T, keyCols);
if isempty(G1)
    nEvents = zeros(0,1);
else
    nEvents = accumarray(G1, 1, [max(G1) 1], @sum, 0);
end
counts = grp1;
counts.n_events = nEvents;

% Optionally include rows that are implicitly zero for missing cell/epoch combos.
if ~p.Results.dropZeroRates
    counts = fillMissingCellEpochCombos(counts, T, grpCols, A, p.Results.cells, files);
end
counts = attachEpochDurations(counts, A);
counts.rate_events_per_min = counts.n_events ./ counts.duration_min;

% 6) Group-level summary.
sumCols = unique([grpCols, {'epoch_idx','epoch_name'}], 'stable');
[G2, grp2] = groupTableRows(counts, sumCols);
vals = counts.rate_events_per_min;
if isempty(G2)
    meanRate = zeros(0,1);
    medianRate = zeros(0,1);
    stdRate = zeros(0,1);
    nCells = zeros(0,1);
else
    nGrp = max(G2);
    meanRate = NaN(nGrp,1);
    medianRate = NaN(nGrp,1);
    stdRate = NaN(nGrp,1);
    nCells = zeros(nGrp,1);
    for gg = 1:nGrp
        x = vals(G2 == gg);
        x = x(isfinite(x));
        nCells(gg) = numel(x);
        if isempty(x)
            continue
        end
        meanRate(gg) = mean(x);
        medianRate(gg) = median(x);
        stdRate(gg) = std(x);
    end
end
S = grp2;
S.mean_rate_events_per_min = meanRate;
S.median_rate_events_per_min = medianRate;
S.std_rate_events_per_min = stdRate;
S.n_cells = nCells;
S.sem_rate_events_per_min = S.std_rate_events_per_min ./ sqrt(max(S.n_cells,1));

out = struct();
out.ratesByCellEpoch = counts;
out.summaryByGroupEpoch = S;
out.eventsAssigned = E;
out.epochTable = A.epochTable;
out.groupBy = grpCols;
out.options = p.Results;
out.n_events_input = nInput;
out.n_events_excluded_non_epoch = nExcludedNonEpoch;
out.n_events_assigned_epoch = nAssignedEpoch;
end

function d = epochDurationsMin(epochTable)
if ismember('duration_min', epochTable.Properties.VariableNames) && all(isfinite(epochTable.duration_min))
    d = epochTable.duration_min;
elseif ismember('time_start_min', epochTable.Properties.VariableNames) && ismember('time_end_min', epochTable.Properties.VariableNames)
    d = epochTable.time_end_min - epochTable.time_start_min;
elseif ismember('frame_start', epochTable.Properties.VariableNames) && ismember('frame_end', epochTable.Properties.VariableNames)
    % Fallback: frame counts if time not provided.
    d = epochTable.frame_end - epochTable.frame_start + 1;
else
    error('epochTable missing duration fields');
end
d = d(:);
end

function out = fillMissingCellEpochCombos(counts, T, grpCols, assigned, cells, files)
cellCols = unique([{'file_id'}, grpCols, {'cell_uid'}], 'stable');
if ~isempty(cells) && all(ismember({'file_id','cell_uid'}, cells.Properties.VariableNames))
    baseCols = {'file_id','cell_uid'};
    cellKeys = cells(:, baseCols);
    metadataCols = unique([{'file_id'}, grpCols], 'stable');
    metadataCols = metadataCols(ismember(metadataCols, files.Properties.VariableNames));
    cellKeys = outerjoin(cellKeys, files(:, metadataCols), 'Keys', 'file_id', ...
        'Type', 'left', 'MergeKeys', true);
    for ii = 1:numel(cellCols)
        if ~ismember(cellCols{ii}, cellKeys.Properties.VariableNames)
            cellKeys.(cellCols{ii}) = repmat("", height(cellKeys), 1);
        end
    end
    cellKeys = unique(cellKeys(:, cellCols));
else
    cellKeys = unique(T(:, cellCols));
end
if isfield(assigned, 'assignmentTable') && istable(assigned.assignmentTable)
    allEpochs = unique(assigned.assignmentTable(:, {'file_id','epoch_idx','epoch_name'}), 'rows', 'stable');
    fullGrid = innerjoin(cellKeys, allEpochs, 'Keys', 'file_id');
else
    allEpochs = assigned.epochTable(:, {'epoch_idx','epoch_name'});
    cellKeys.tmp_join_key = ones(height(cellKeys),1);
    allEpochs.tmp_join_key = ones(height(allEpochs),1);
    fullGrid = innerjoin(cellKeys, allEpochs, 'Keys', 'tmp_join_key');
    fullGrid.tmp_join_key = [];
end

joinCols = unique([{'file_id'}, grpCols, {'cell_uid','epoch_idx','epoch_name'}], 'stable');
out = outerjoin(fullGrid, counts, 'Keys', joinCols, 'MergeKeys', true, 'Type', 'left');
if ~ismember('n_events', out.Properties.VariableNames)
    out.n_events = zeros(height(out),1);
end
out.n_events(isnan(out.n_events)) = 0;
end

function counts = attachEpochDurations(counts, assigned)
if isfield(assigned, 'assignmentTable') && istable(assigned.assignmentTable)
    D = unique(assigned.assignmentTable(:, {'file_id','epoch_idx','duration_min'}), 'rows', 'stable');
    counts = outerjoin(counts, D, 'Keys', {'file_id','epoch_idx'}, 'Type', 'left', 'MergeKeys', true);
else
    dur = epochDurationsMin(assigned.epochTable);
    counts.duration_min = NaN(height(counts),1);
    for ii = 1:height(counts)
        ei = counts.epoch_idx(ii);
        if ei >= 1 && ei <= numel(dur), counts.duration_min(ii) = dur(ei); end
    end
end
if any(~isfinite(counts.duration_min) | counts.duration_min <= 0)
    error('computeNEventsPerMinCanonical:EpochDuration', ...
        'Every cell/epoch row requires a positive duration in minutes.');
end
end

function [G, grpTbl] = groupTableRows(T, cols)
n = height(T);
if n == 0
    G = zeros(0,1);
    grpTbl = T(:, cols);
    return
end

key = strings(n,1);
for ii = 1:numel(cols)
    c = cols{ii};
    key = key + "|" + toKeyStrings(T.(c));
end
[~, ~, G] = unique(key, 'stable');
G = full(double(G(:)));

firstIdx = accumarray(G, (1:n)', [], @min);
grpTbl = T(firstIdx, cols);
end

function s = toKeyStrings(v)
if isstring(v)
    s = v(:);
elseif iscell(v)
    s = strings(numel(v),1);
    for i = 1:numel(v)
        vi = v{i};
        if ischar(vi)
            s(i) = string(vi);
        elseif isstring(vi)
            s(i) = vi;
        elseif isempty(vi)
            s(i) = "";
        else
            s(i) = string(vi);
        end
    end
elseif ischar(v)
    s = repmat(string(v), size(v,1), 1);
elseif isnumeric(v) || islogical(v)
    s = string(v(:));
elseif iscategorical(v)
    s = string(v(:));
else
    s = string(v);
    s = s(:);
end
end
