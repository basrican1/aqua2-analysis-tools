function out = subsetLocalResultEvents(result, varargin)
%SUBSETLOCALRESULTEVENTS Filter a local-analysis result to a chosen event subset.
%   OUT = SUBSETLOCALRESULTEVENTS(RESULT, 'EventUID', keepUID)
%   OUT = SUBSETLOCALRESULTEVENTS(RESULT, 'MinParticipatingCells', 2)
%
%   This helper keeps only rows whose event_uid belongs to the selected
%   subset across metrics tables and canonical subset tables when present.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'result', @isstruct);
addParameter(p, 'EventUID', strings(0,1), @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'MinParticipatingCells', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'MinRecruitedCells', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'MinAreaPx2', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'MaxAreaPx2', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'MinAreaUm2', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'MaxAreaUm2', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'MinDurationSec', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'MaxDurationSec', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'PruneForSave', true, @islogical);
addParameter(p, 'KeepCanonicalMetrics', false, @islogical);
addParameter(p, 'KeepDerived', false, @islogical);
parse(p, result, varargin{:});

out = result;
keepUID = localResolveKeepUID(result, string(p.Results.EventUID(:)), ...
    double(p.Results.MinParticipatingCells), double(p.Results.MinRecruitedCells), ...
    double(p.Results.MinAreaPx2), double(p.Results.MaxAreaPx2), ...
    double(p.Results.MinAreaUm2), double(p.Results.MaxAreaUm2), ...
    double(p.Results.MinDurationSec), double(p.Results.MaxDurationSec));
keepUID = unique(keepUID(strlength(keepUID) > 0), 'stable');
if isempty(keepUID)
    warning('subsetLocalResultEvents:EmptySubset', ...
        'No events matched the requested subset. Returning tables with zero kept events.');
end

if isfield(out, 'metrics') && isstruct(out.metrics) && ...
        isfield(out.metrics, 'tables') && isstruct(out.metrics.tables)
    tnames = fieldnames(out.metrics.tables);
    for i = 1:numel(tnames)
        T = out.metrics.tables.(tnames{i});
        if istable(T) && ismember('event_uid', T.Properties.VariableNames)
            out.metrics.tables.(tnames{i}) = T(ismember(string(T.event_uid), keepUID), :);
        end
    end
end

if isfield(out, 'canonical') && isstruct(out.canonical) && ...
        isfield(out.canonical, 'subset') && isstruct(out.canonical.subset)
    out.canonical.subset = localSubsetCanonicalSubset(out.canonical.subset, keepUID);
end

if p.Results.PruneForSave
    out = localPruneSubsetResult(out, ...
        logical(p.Results.KeepCanonicalMetrics), ...
        logical(p.Results.KeepDerived));
end
end

function keepUID = localResolveKeepUID(result, explicitUID, minCells, minRecruited, ...
    minAreaPx2, maxAreaPx2, minAreaUm2, maxAreaUm2, minDurationSec, maxDurationSec)
keepUID = unique(explicitUID(strlength(explicitUID) > 0), 'stable');
if ~isempty(keepUID)
    return
end

if ~(isfield(result, 'metrics') && isstruct(result.metrics) && ...
        isfield(result.metrics, 'tables') && isstruct(result.metrics.tables) && ...
        isfield(result.metrics.tables, 'event_partition') && istable(result.metrics.tables.event_partition))
    keepUID = strings(0,1);
    return
end

P = result.metrics.tables.event_partition;
if ~ismember('event_uid', P.Properties.VariableNames)
    keepUID = strings(0,1);
    return
end
cellCol = localFirstPresent(P.Properties.VariableNames, {'partition_cell_uid','cell_uid'});
if isempty(cellCol)
    keepUID = unique(string(P.event_uid), 'stable');
    return
end

eventUID = string(P.event_uid);
cellUID = string(P.(cellCol));
[G, ~] = findgroups(eventUID);
nCells = splitapply(@(x) numel(unique(x(x ~= ""))), cellUID, G);
baseUID = splitapply(@(x) string(x(1)), eventUID, G);

keep = true(size(baseUID));
if isfinite(minCells)
    keep = keep & nCells >= minCells;
end

if isfinite(minRecruited)
    nRecruited = zeros(size(baseUID));
    if ismember('is_recruited_partition', P.Properties.VariableNames)
        isRec = logical(double(P.is_recruited_partition));
        nRecruited = splitapply(@(cid, rec) numel(unique(cid(rec & cid ~= ""))), cellUID, isRec, G);
    end
    keep = keep & nRecruited >= minRecruited;
end

eventTbl = localBuildEventLevelTable(result);
if ~isempty(eventTbl)
    eventKeep = true(height(eventTbl), 1);
    if isfinite(minAreaPx2) && ismember('area_px2', eventTbl.Properties.VariableNames)
        eventKeep = eventKeep & double(eventTbl.area_px2) >= minAreaPx2;
    end
    if isfinite(maxAreaPx2) && ismember('area_px2', eventTbl.Properties.VariableNames)
        eventKeep = eventKeep & double(eventTbl.area_px2) <= maxAreaPx2;
    end
    if isfinite(minAreaUm2) && ismember('area_um2', eventTbl.Properties.VariableNames)
        eventKeep = eventKeep & double(eventTbl.area_um2) >= minAreaUm2;
    end
    if isfinite(maxAreaUm2) && ismember('area_um2', eventTbl.Properties.VariableNames)
        eventKeep = eventKeep & double(eventTbl.area_um2) <= maxAreaUm2;
    end
    if isfinite(minDurationSec) && ismember('duration_sec', eventTbl.Properties.VariableNames)
        eventKeep = eventKeep & double(eventTbl.duration_sec) >= minDurationSec;
    end
    if isfinite(maxDurationSec) && ismember('duration_sec', eventTbl.Properties.VariableNames)
        eventKeep = eventKeep & double(eventTbl.duration_sec) <= maxDurationSec;
    end
    eventKeepUID = unique(string(eventTbl.event_uid(eventKeep)), 'stable');
    keep = keep & ismember(baseUID, eventKeepUID);
end

keepUID = unique(baseUID(keep), 'stable');
end

function T = localBuildEventLevelTable(result)
T = table();
if isfield(result, 'metrics') && isstruct(result.metrics) && ...
        isfield(result.metrics, 'tables') && isstruct(result.metrics.tables)
    tables = result.metrics.tables;
    T = localMergeEventLevelSource(T, tables, 'event_canonical', {'event_uid','area_px2','area_um2','duration_sec'});
    T = localMergeEventLevelSource(T, tables, 'event_shape', {'event_uid','area_px2','area_um2','duration_sec'});
    T = localMergeEventLevelSource(T, tables, 'events', {'event_uid','area_px2','area_um2','duration_sec'});
end

if isempty(T)
    T = localBuildEventLevelTableFromCanonical(result);
end

if ~isempty(T)
    T = unique(T, 'rows', 'stable');
end
end

function T = localMergeEventLevelSource(T, tables, tableName, keepCols)
if ~(isfield(tables, tableName) && istable(tables.(tableName)))
    return
end

keep = intersect(keepCols, tables.(tableName).Properties.VariableNames, 'stable');
if ~ismember('event_uid', keep)
    return
end

Ts = unique(tables.(tableName)(:, keep), 'rows', 'stable');
if isempty(T)
    T = Ts;
else
    T = outerjoin(T, Ts, 'Keys', 'event_uid', 'Type', 'left', 'MergeKeys', true);
end
end

function T = localBuildEventLevelTableFromCanonical(result)
T = table();
if ~(isfield(result, 'canonical') && isstruct(result.canonical))
    return
end

if isfield(result.canonical, 'subset') && isstruct(result.canonical.subset) && ...
        isfield(result.canonical.subset, 'events') && istable(result.canonical.subset.events)
    E = result.canonical.subset.events;
    keep = intersect({'event_uid','area_px2','area_um2','duration_sec'}, E.Properties.VariableNames, 'stable');
    if ismember('event_uid', keep)
        T = unique(E(:, keep), 'rows', 'stable');
    end
end

if ~(isfield(result.canonical, 'metrics') && isstruct(result.canonical.metrics))
    return
end

T = localAttachCanonicalMetricEventValues(T, result.canonical.metrics, 'area_um2', 'area_um2');
T = localAttachCanonicalMetricEventValues(T, result.canonical.metrics, 'duration_sec', 'duration_sec');
end

function T = localAttachCanonicalMetricEventValues(T, metrics, metricField, valueName)
if ~(isfield(metrics, metricField) && isstruct(metrics.(metricField)) && ...
        isfield(metrics.(metricField), 'eventValues') && istable(metrics.(metricField).eventValues))
    return
end

V = metrics.(metricField).eventValues;
keep = intersect({'event_uid', valueName}, V.Properties.VariableNames, 'stable');
if numel(keep) < 2 || ~ismember('event_uid', keep)
    return
end

Vs = unique(V(:, keep), 'rows', 'stable');
if isempty(T)
    T = Vs;
else
    T = outerjoin(T, Vs, 'Keys', 'event_uid', 'Type', 'left', 'MergeKeys', true);
end
end

function subset = localSubsetCanonicalSubset(subset, keepUID)
if isfield(subset, 'events') && istable(subset.events) && ismember('event_uid', subset.events.Properties.VariableNames)
    subset.events = subset.events(ismember(string(subset.events.event_uid), keepUID), :);
end
if isfield(subset, 'event_cell_map') && istable(subset.event_cell_map) && ismember('event_uid', subset.event_cell_map.Properties.VariableNames)
    subset.event_cell_map = subset.event_cell_map(ismember(string(subset.event_cell_map.event_uid), keepUID), :);
end
if isfield(subset, 'cells') && istable(subset.cells)
    subset.cells = localSubsetCellsByEventMap(subset, subset.cells);
end
end

function cells = localSubsetCellsByEventMap(subset, cells)
if ~(isfield(subset, 'event_cell_map') && istable(subset.event_cell_map) && ...
        ismember('cell_uid', subset.event_cell_map.Properties.VariableNames) && ...
        ismember('cell_uid', cells.Properties.VariableNames))
    return
end
keepCell = unique(string(subset.event_cell_map.cell_uid), 'stable');
cells = cells(ismember(string(cells.cell_uid), keepCell), :);
end

function name = localFirstPresent(names, candidates)
name = '';
for i = 1:numel(candidates)
    if any(strcmp(names, candidates{i}))
        name = candidates{i};
        return
    end
end
end

function out = localPruneSubsetResult(out, keepCanonicalMetrics, keepDerived)
topDrop = {'canonical_norm','canonical_norms','metrics_by_transform', ...
    'normalization_by_transform','normalized','stats','metric_qc_report'};
for i = 1:numel(topDrop)
    if isfield(out, topDrop{i})
        out = rmfield(out, topDrop{i});
    end
end

if isfield(out, 'canonical') && isstruct(out.canonical)
    canonDrop = {'loaded','frequency','plan','planMeta','qcSummary'};
    if ~keepCanonicalMetrics
        canonDrop{end+1} = 'metrics'; %#ok<AGROW>
    end
    for i = 1:numel(canonDrop)
        if isfield(out.canonical, canonDrop{i})
            out.canonical = rmfield(out.canonical, canonDrop{i});
        end
    end
end

if ~keepDerived && isfield(out, 'derived')
    out = rmfield(out, 'derived');
end
end
