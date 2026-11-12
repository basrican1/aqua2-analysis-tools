function prov = bwaResolveNormalizationProvenance(result, dataView, metricIDs, unitLevels)
%BWARESOLVENORMALIZATIONPROVENANCE Describe effective normalization usage.
%
% Inputs:
%   result     : local analysis result struct
%   dataView   : 'raw' | 'norm' | 'recommended'
%   metricIDs  : optional metric IDs to resolve per-metric transforms
%   unitLevels : optional unit levels aligned with metricIDs

if nargin < 3 || isempty(metricIDs)
    metricIDs = strings(0,1);
end
if nargin < 4 || isempty(unitLevels)
    unitLevels = strings(size(metricIDs));
end

metricIDs = string(metricIDs(:));
unitLevels = string(unitLevels(:));
if numel(unitLevels) ~= numel(metricIDs)
    unitLevels = repmat("", size(metricIDs));
end

view = lower(char(string(dataView)));
mode = "";
defaultTransform = "";
enabled = false;
baselineEpochIdx = NaN;
baselineEpochName = "";
if isstruct(result) && isfield(result, 'normalization') && isstruct(result.normalization)
    if isfield(result.normalization, 'mode') && ~isempty(result.normalization.mode)
        mode = string(result.normalization.mode);
        enabled = ~strcmpi(mode, "none");
    end
    if isfield(result.normalization, 'transform') && ~isempty(result.normalization.transform)
        defaultTransform = string(result.normalization.transform);
    end
    if isfield(result.normalization, 'baselineEpochIdx') && ~isempty(result.normalization.baselineEpochIdx)
        baselineEpochIdx = double(result.normalization.baselineEpochIdx);
    end
end
if defaultTransform == "" && isstruct(result) && isfield(result, 'normalized') && isstruct(result.normalized) && ...
        isfield(result.normalized, 'default_transform') && ~isempty(result.normalized.default_transform)
    defaultTransform = string(result.normalized.default_transform);
end
if mode == ""
    mode = "none";
end
if defaultTransform == ""
    defaultTransform = "delta";
end
baselineEpochName = inferBaselineEpochName(result, baselineEpochIdx);

resolvedTransforms = strings(size(metricIDs));
for i = 1:numel(metricIDs)
    mid = metricIDs(i);
    ul = unitLevels(i);
    if strlength(ul) == 0
        ul = inferUnitLevel(mid);
    end
    resolvedTransforms(i) = resolveTransformForMetric(mid, ul, view, defaultTransform);
end

uniqTransforms = unique(resolvedTransforms(strlength(resolvedTransforms) > 0), 'stable');
summaryLabel = buildSummaryLabel(view, mode, defaultTransform, enabled, uniqTransforms);
compactLabel = buildCompactLabel(view, defaultTransform, enabled, uniqTransforms);

prov = struct();
prov.requestedView = view;
prov.normalizationEnabled = logical(enabled);
prov.normalizationMode = char(mode);
prov.defaultTransform = char(defaultTransform);
prov.baselineEpochIdx = baselineEpochIdx;
prov.baselineEpochName = char(baselineEpochName);
prov.summaryLabel = char(summaryLabel);
prov.compactLabel = char(compactLabel);
prov.metricIDs = cellstr(metricIDs)';
prov.unitLevels = cellstr(unitLevels)';
prov.metricTransforms = cellstr(resolvedTransforms)';
prov.uniqueTransforms = cellstr(uniqTransforms)';
if isempty(metricIDs)
    prov.metricTable = table(strings(0,1), strings(0,1), strings(0,1), ...
        'VariableNames', {'metric_id','unit_level','effective_transform'});
else
    prov.metricTable = table(metricIDs, unitLevels, resolvedTransforms, ...
        'VariableNames', {'metric_id','unit_level','effective_transform'});
end
end

function out = resolveTransformForMetric(metricID, unitLevel, view, defaultTransform)
switch view
    case 'raw'
        out = "raw";
    case 'norm'
        out = string(defaultTransform);
    case 'recommended'
        if strlength(metricID) == 0
            out = "mixed";
        else
            out = string(bwaPreferredTransformForMetricID(metricID, unitLevel));
        end
    otherwise
        out = "";
end
end

function out = inferUnitLevel(metricID)
metricID = string(metricID);
if startsWith(metricID, "cell.")
    out = "cell";
elseif startsWith(metricID, "slice.")
    out = "slice";
elseif startsWith(metricID, "event.")
    out = "event";
else
    out = "canonical";
end
end

function out = inferBaselineEpochName(result, baselineEpochIdx)
out = "";
if ~(isfinite(baselineEpochIdx) && baselineEpochIdx > 0)
    return
end
if isstruct(result) && isfield(result, 'canonical') && isstruct(result.canonical) && ...
        isfield(result.canonical, 'frequency') && isstruct(result.canonical.frequency) && ...
        isfield(result.canonical.frequency, 'epochTable') && istable(result.canonical.frequency.epochTable)
    T = result.canonical.frequency.epochTable;
    if ismember('epoch_idx', T.Properties.VariableNames) && ismember('epoch_name', T.Properties.VariableNames)
        hit = find(double(T.epoch_idx) == baselineEpochIdx, 1, 'first');
        if ~isempty(hit)
            out = string(T.epoch_name(hit));
        end
    end
end
end

function out = buildSummaryLabel(view, mode, defaultTransform, enabled, uniqTransforms)
switch view
    case 'raw'
        out = "raw values";
    case 'norm'
        if enabled
            out = sprintf('default normalized view [%s; mode=%s]', char(defaultTransform), char(mode));
        else
            out = 'requested norm view, but normalization is disabled';
        end
    case 'recommended'
        if isempty(uniqTransforms)
            out = 'recommended metric-aware normalization';
        elseif numel(uniqTransforms) == 1
            out = sprintf('recommended normalization [%s]', char(uniqTransforms));
        else
            out = sprintf('recommended normalization [mixed: %s]', strjoin(cellstr(uniqTransforms'), ', '));
        end
    otherwise
        out = string(view);
end
end

function out = buildCompactLabel(view, defaultTransform, enabled, uniqTransforms)
switch view
    case 'raw'
        out = "norm=raw";
    case 'norm'
        if enabled
            out = sprintf('norm=default [%s]', char(defaultTransform));
        else
            out = 'norm=disabled';
        end
    case 'recommended'
        if isempty(uniqTransforms)
            out = 'norm=recommended';
        elseif numel(uniqTransforms) == 1
            out = sprintf('norm=%s', char(uniqTransforms));
        else
            out = sprintf('norm=mixed [%s]', strjoin(cellstr(uniqTransforms'), ', '));
        end
    otherwise
        out = sprintf('norm=%s', view);
end
end
