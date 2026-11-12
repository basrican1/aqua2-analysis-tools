function T = list_top_recruited_events(source, varargin)
%LIST_TOP_RECRUITED_EVENTS Rank events by recruited-cell count.
%   T = LIST_TOP_RECRUITED_EVENTS(result)
%   T = LIST_TOP_RECRUITED_EVENTS(pathIn)
%
%   Returns one row per event, ranked by number of recruited cells
%   (neighboring non-primary cells in metrics.tables.event_partition).
%
%   Name-value options:
%     'TopN'      : number of rows to return (default 10)
%     'EpochName' : optional epoch filter
%     'FileID'    : optional file_id filter
%     'EventsTbl' : optional canonical events table override for tests

ip = inputParser;
ip.addRequired('source');
ip.addParameter('TopN', 10, @(x) isnumeric(x) && isscalar(x) && x >= 1);
ip.addParameter('EpochName', "", @(x) ischar(x) || isstring(x));
ip.addParameter('FileID', "", @(x) ischar(x) || isstring(x));
ip.addParameter('EventsTbl', table(), @(x) istable(x));
ip.parse(source, varargin{:});
opts = ip.Results;

result = localResolveResult(source);
eventPartition = localRequireEventPartition(result);
eventPartition = localApplyFilters(eventPartition, string(opts.EpochName), string(opts.FileID));

if isempty(eventPartition)
    T = table();
    return
end

isNeighbor = true(height(eventPartition), 1);
if ismember('is_primary_cell', eventPartition.Properties.VariableNames)
    isPrimary = false(height(eventPartition), 1);
    try
        isPrimary = logical(eventPartition.is_primary_cell);
    catch
        isPrimary = false(height(eventPartition), 1);
    end
    isNeighbor = ~isPrimary;
elseif all(ismember({'partition_cell_uid','primary_cell_uid'}, eventPartition.Properties.VariableNames))
    isNeighbor = string(eventPartition.partition_cell_uid) ~= string(eventPartition.primary_cell_uid);
end

neighborRows = eventPartition(isNeighbor, :);
if isempty(neighborRows)
    T = table();
    return
end

groupVars = {'event_uid'};
if ismember('file_id', neighborRows.Properties.VariableNames)
    groupVars{end+1} = 'file_id'; %#ok<AGROW>
end
if ismember('epoch_name', neighborRows.Properties.VariableNames)
    groupVars{end+1} = 'epoch_name'; %#ok<AGROW>
end

cellVar = localPickCellVar(neighborRows);
G = groupsummary(neighborRows, groupVars, @(x) numel(unique(string(x))), cellVar);
G.Properties.VariableNames{end} = 'n_recruited_cells';

if isfield(result.metrics.tables, 'event_spread') && ~isempty(result.metrics.tables.event_spread)
    eventSpread = result.metrics.tables.event_spread;
    keep = intersect({'event_uid','file_id','spread_count','is_multicell','spread_threshold_pct'}, ...
        eventSpread.Properties.VariableNames, 'stable');
    if ~isempty(keep)
        G = outerjoin(G, eventSpread(:, keep), 'Keys', intersect({'event_uid','file_id'}, keep, 'stable'), ...
            'MergeKeys', true, 'Type', 'left');
    end
end

eventsTbl = opts.EventsTbl;
if isempty(eventsTbl)
    eventsTbl = localResolveEventsTable(source, result);
end
if ~isempty(eventsTbl)
    eventKeep = intersect({'event_uid','file_id','local_event_id','area_px2','area_um2','duration_sec'}, ...
        eventsTbl.Properties.VariableNames, 'stable');
    if ~isempty(eventKeep)
        G = outerjoin(G, eventsTbl(:, eventKeep), 'Keys', intersect({'event_uid','file_id'}, eventKeep, 'stable'), ...
            'MergeKeys', true, 'Type', 'left');
    end
end

sortVars = intersect({'n_recruited_cells','area_um2','area_px2','spread_count','event_uid'}, G.Properties.VariableNames, 'stable');
sortDir = repmat({'descend'}, 1, max(0, numel(sortVars) - 1));
if ~isempty(sortVars)
    sortDir{end+1} = 'ascend'; %#ok<AGROW>
    G = sortrows(G, sortVars, sortDir);
end

topN = min(height(G), floor(double(opts.TopN)));
T = G(1:topN, :);
end

function result = localResolveResult(source)
if isstruct(source) && isfield(source, 'metrics')
    result = source;
    return
end

if ischar(source) || isstring(source)
    root = char(string(source));
    if exist(root, 'dir') ~= 7
        error('list_top_recruited_events:InvalidPath', ...
            'Expected a result struct or dataset root path. Missing: %s', root);
    end
    d = dir(fullfile(root, 'derived', 'local_outputs', 'preset_*', 'analysis_result.mat'));
    if isempty(d)
        error('list_top_recruited_events:MissingResult', ...
            'No analysis_result.mat found under %s/derived/local_outputs.', root);
    end
    [~, idx] = max([d.datenum]);
    S = load(fullfile(d(idx).folder, d(idx).name), 'result');
    result = load_local_analysis_result(S.result);
    return
end

error('list_top_recruited_events:InvalidSource', ...
    'Expected a local-analysis result struct or dataset root path.');
end

function eventPartition = localRequireEventPartition(result)
eventPartition = table();
if isfield(result, 'metrics') && isfield(result.metrics, 'tables') && ...
        isfield(result.metrics.tables, 'event_partition')
    eventPartition = result.metrics.tables.event_partition;
end
if isempty(eventPartition)
    error('list_top_recruited_events:MissingEventPartition', ...
        'Result does not contain metrics.tables.event_partition.');
end
end

function T = localApplyFilters(T, epochName, fileID)
if strlength(epochName) > 0 && ismember('epoch_name', T.Properties.VariableNames)
    T = T(string(T.epoch_name) == epochName, :);
end
if strlength(fileID) > 0 && ismember('file_id', T.Properties.VariableNames)
    T = T(string(T.file_id) == fileID, :);
end
end

function cellVar = localPickCellVar(T)
choices = {'partition_cell_uid','cell_uid','secondary_cell_uid'};
idx = find(ismember(choices, T.Properties.VariableNames), 1, 'first');
if isempty(idx)
    error('list_top_recruited_events:MissingCellColumn', ...
        'Could not identify a recruited-cell column in event_partition.');
end
cellVar = choices{idx};
end

function eventsTbl = localResolveEventsTable(source, result)
eventsTbl = table();
pathIn = "";
if isstruct(source) && isfield(source, 'pathIn') && ~isempty(source.pathIn)
    pathIn = string(source.pathIn);
elseif isstruct(result) && isfield(result, 'pathIn') && ~isempty(result.pathIn)
    pathIn = string(result.pathIn);
elseif ischar(source) || isstring(source)
    pathIn = string(source);
end

if strlength(pathIn) == 0 || exist(char(pathIn), 'dir') ~= 7
    return
end

loaded = loadClusterOutputs(char(pathIn), 'attachLegacyMeta', false, 'attachEventMeta', false);
if isfield(loaded, 'events') && istable(loaded.events)
    eventsTbl = loaded.events;
elseif isfield(loaded, 'canonical') && isfield(loaded.canonical, 'subset') && ...
        isfield(loaded.canonical.subset, 'events') && istable(loaded.canonical.subset.events)
    eventsTbl = loaded.canonical.subset.events;
end
end
