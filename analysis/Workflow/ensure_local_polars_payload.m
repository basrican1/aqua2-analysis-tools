function [polData, result, meta] = ensure_local_polars_payload(source, varargin)
% Ensure a usable derived.polars payload exists for spider plotting.

p = inputParser;
p.StructExpand = false;
addParameter(p, 'requireEntityValues', false, @islogical);
addParameter(p, 'requireEntityValuesByFile', false, @islogical);
parse(p, varargin{:});

result = load_local_analysis_result(source);
meta = struct('status', 'existing', 'message', '');
polData = struct();

if isfield(result, 'derived') && isstruct(result.derived) && isfield(result.derived, 'polars')
    polData = result.derived.polars;
end
if hasRequiredPolars(polData) && ...
        (~p.Results.requireEntityValues || hasEntityPolars(polData)) && ...
        (~p.Results.requireEntityValuesByFile || hasEntityPolarsByFile(polData))
    return
end

try
    polData = rebuildPolars(result);
    if ~isfield(result, 'derived') || ~isstruct(result.derived)
        result.derived = struct();
    end
    result.derived.polars = polData;
    meta.status = 'rebuilt';
catch ME
    polData = struct();
    meta.status = 'missing';
    meta.message = ME.message;
end
end

function tf = hasEntityPolars(polData)
tf = isstruct(polData) && isfield(polData, 'labels') && ...
    isfield(polData, 'entityValuesA') && iscell(polData.entityValuesA) && ...
    isfield(polData, 'entityValuesC') && iscell(polData.entityValuesC) && ...
    numel(polData.entityValuesA) == numel(polData.labels) && ...
    numel(polData.entityValuesC) == numel(polData.labels);
end

function tf = hasEntityPolarsByFile(polData)
tf = isstruct(polData) && isfield(polData, 'labels') && ...
    isfield(polData, 'entityValuesAByFile') && iscell(polData.entityValuesAByFile) && ...
    isfield(polData, 'entityValuesCByFile') && iscell(polData.entityValuesCByFile) && ...
    numel(polData.entityValuesAByFile) == numel(polData.labels) && ...
    numel(polData.entityValuesCByFile) == numel(polData.labels);
end

function tf = hasRequiredPolars(polData)
% subplotGrid is layout metadata and can be filled later by the plotting
% wrapper, so it should not force a full rebuild here.
need = {'labels','valuesA','valuesC','uniqueFiles'};
tf = isstruct(polData);
if ~tf
    return
end
for i = 1:numel(need)
    if ~isfield(polData, need{i}) || isempty(polData.(need{i}))
        tf = false;
        return
    end
end
end

function polData = rebuildPolars(result)
if ~isfield(result, 'pathIn') || isempty(result.pathIn)
    error('ensure_local_polars_payload:NoPath', ...
        'Result does not include pathIn, so spider payload cannot be rebuilt automatically.');
end

loaded = loadClusterOutputs(char(string(result.pathIn)), 'attachLegacyMeta', false, 'attachEventMeta', false);
if ~isfield(loaded, 'BWA_AQ2') || ~isstruct(loaded.BWA_AQ2) || isempty(fieldnames(loaded.BWA_AQ2))
    error('ensure_local_polars_payload:NoBWA', ...
        'BWA_AQ2.mat is required to rebuild spider payload from a canonical-only export.');
end

tfre = [];
if isfield(result, 'derived') && isstruct(result.derived) && ...
        isfield(result.derived, 'context') && isstruct(result.derived.context) && ...
        isfield(result.derived.context, 'timeFrameRangeEpochs')
    tfre = result.derived.context.timeFrameRangeEpochs;
elseif isfield(result, 'canonical') && isstruct(result.canonical) && isfield(result.canonical, 'epochTable')
    tfre = epochSpecFromEpochTable(result.canonical.epochTable);
end
if isempty(tfre)
    error('ensure_local_polars_payload:NoEpochSpec', ...
        'Missing epoch specification for spider payload rebuild.');
end

filesTbl = table();
eventsTbl = table();
if isfield(result, 'canonical') && isstruct(result.canonical) && ...
        isfield(result.canonical, 'subset') && isstruct(result.canonical.subset)
    if isfield(result.canonical.subset, 'files') && istable(result.canonical.subset.files)
        filesTbl = result.canonical.subset.files;
    end
    if isfield(result.canonical.subset, 'events') && istable(result.canonical.subset.events)
        eventsTbl = result.canonical.subset.events;
    end
end
if isempty(filesTbl) && isfield(loaded, 'files') && istable(loaded.files)
    filesTbl = loaded.files;
end
if isempty(eventsTbl) && isfield(loaded, 'events') && istable(loaded.events)
    eventsTbl = loaded.events;
end

dataLimitRules = {};
if isfield(result, 'dataLimitRules') && iscell(result.dataLimitRules)
    dataLimitRules = result.dataLimitRules;
end
combineEpochNames = {};
combinedEpochLabel = 'CNO_post';
if isfield(result, 'spiderCombineEpochNames')
    combineEpochNames = result.spiderCombineEpochNames;
elseif isfield(result, 'derived') && isstruct(result.derived) && isfield(result.derived, 'polars') && ...
        isstruct(result.derived.polars) && isfield(result.derived.polars, 'combineEpochNames')
    combineEpochNames = result.derived.polars.combineEpochNames;
end
if isfield(result, 'spiderCombinedEpochLabel') && ~isempty(result.spiderCombinedEpochLabel)
    combinedEpochLabel = result.spiderCombinedEpochLabel;
elseif isfield(result, 'derived') && isstruct(result.derived) && isfield(result.derived, 'polars') && ...
        isstruct(result.derived.polars) && isfield(result.derived.polars, 'combinedEpochLabel') && ...
        ~isempty(result.derived.polars.combinedEpochLabel)
    combinedEpochLabel = result.derived.polars.combinedEpochLabel;
end

polData = computeBWA_Aqua2PolarsBySliceData(loaded.BWA_AQ2, ...
    'timeFrameRangeEpochs', tfre, ...
    'filesTbl', filesTbl, ...
    'eventsTbl', eventsTbl, ...
    'dataLimitRules', dataLimitRules, ...
    'combineEpochNames', combineEpochNames, ...
    'combinedEpochLabel', combinedEpochLabel, ...
    'includeStandardFeatures', true, ...
    'filePlots', false);
end

function tfre = epochSpecFromEpochTable(T)
tfre = [];
if ~istable(T) || isempty(T) || ~all(ismember({'epoch_idx','epoch_name'}, T.Properties.VariableNames))
    return
end
n = max(double(T.epoch_idx));
if ~isfinite(n) || n <= 0
    return
end
tfre = struct();
tfre.epochNames = cell(n,1);
tfre.frameRange = nan(n,2);
tfre.timeRange = nan(n,2);
for i = 1:n
    row = T(double(T.epoch_idx) == i, :);
    if isempty(row)
        continue
    end
    tfre.epochNames{i,1} = char(string(row.epoch_name(1)));
    if ismember('frame_start', row.Properties.VariableNames)
        tfre.frameRange(i,1) = double(row.frame_start(1));
    end
    if ismember('frame_end', row.Properties.VariableNames)
        tfre.frameRange(i,2) = double(row.frame_end(1));
    end
    if ismember('time_start', row.Properties.VariableNames)
        tfre.timeRange(i,1) = double(row.time_start(1));
    end
    if ismember('time_end', row.Properties.VariableNames)
        tfre.timeRange(i,2) = double(row.time_end(1));
    end
end
end
