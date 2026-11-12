function [metricOut, info] = bwaSelectCanonicalMetricSummary(result, metricName, varargin)
% Select a canonical metric summary view with explicit aggregation policy.
%
% Options:
%   'dataView'             : 'raw' | 'norm' | 'recommended'
%   'canonicalAggregation' : 'cell_promoted' | 'event_native'

p = inputParser;
p.StructExpand = false;
addRequired(p, 'result', @isstruct);
addRequired(p, 'metricName', @(x) ischar(x) || isstring(x));
addParameter(p, 'dataView', 'raw', @(x) ischar(x) || isstring(x));
addParameter(p, 'canonicalAggregation', 'cell_promoted', @(x) ischar(x) || isstring(x));
parse(p, result, metricName, varargin{:});

metricName = char(string(metricName));
dataView = lower(char(string(p.Results.dataView)));
canonicalAggregation = lower(char(string(p.Results.canonicalAggregation)));

metricOut = localSelectNativeCanonicalMetric(result, metricName, dataView);
metricOut.titleName = localResolveDisplayName(metricName) + " [event-native]";
info = struct( ...
    'requestedAggregation', canonicalAggregation, ...
    'effectiveAggregation', "event_native", ...
    'canonical_metric_id', string(metricName), ...
    'promoted_metric_id', "", ...
    'usedPromotedMetric', false);

if ~strcmp(canonicalAggregation, 'cell_promoted')
    return
end

promotedMetricID = localMapCanonicalToPromoted(metricName);
if strlength(promotedMetricID) == 0
    return
end

[tablesOut, regOut, tname, ycol] = localResolvePromotedMetricView(result, promotedMetricID, dataView);
if isempty(tablesOut) || ~isfield(tablesOut, tname)
    return
end
T = tablesOut.(tname);
if ~istable(T) || isempty(T) || ~ismember(ycol, T.Properties.VariableNames)
    return
end

epochTable = localResolveEpochTable(result, dataView);
if isempty(epochTable) || ~istable(epochTable)
    return
end

metricOut = bwaBuildMetricOutFromTable(T, ycol, epochTable, 'entityKeys', {'primary_cell_uid'});
metricOut.titleName = localResolveDisplayName(metricName) + " [cell-promoted]";

info.effectiveAggregation = "cell_promoted";
info.promoted_metric_id = promotedMetricID;
info.usedPromotedMetric = true;

if istable(regOut) && ~isempty(regOut)
    hit = find(string(regOut.metric_id) == promotedMetricID, 1, 'first');
    if ~isempty(hit) && ismember('public_metric_id', regOut.Properties.VariableNames)
        val = string(regOut.public_metric_id(hit));
        if strlength(val) > 0
            metricOut.titleName = val + " [cell-promoted]";
        end
    end
end
end

function metricOut = localSelectNativeCanonicalMetric(result, metricName, dataView)
switch dataView
    case 'recommended'
        metricOut = localSelectRecommendedCanonicalMetric(result, metricName);
    case 'norm'
        if isfield(result, 'canonical_norm') && isfield(result.canonical_norm, 'metrics') && ...
                isfield(result.canonical_norm.metrics, metricName)
            metricOut = result.canonical_norm.metrics.(metricName);
        else
            metricOut = result.canonical.metrics.(metricName);
        end
    otherwise
        metricOut = result.canonical.metrics.(metricName);
end
end

function metricOut = localSelectRecommendedCanonicalMetric(result, metricName)
metricOut = result.canonical.metrics.(metricName);
tr = bwaPreferredTransformForMetricID(metricName, 'canonical');
if ~isfield(result, 'canonical_norms') || ~isstruct(result.canonical_norms)
    return
end
field = matlab.lang.makeValidName(tr);
if isfield(result.canonical_norms, field) && isfield(result.canonical_norms.(field), 'metrics') && ...
        isfield(result.canonical_norms.(field).metrics, metricName)
    metricOut = result.canonical_norms.(field).metrics.(metricName);
elseif isfield(result, 'canonical_norm') && isfield(result.canonical_norm, 'metrics') && ...
        isfield(result.canonical_norm.metrics, metricName)
    metricOut = result.canonical_norm.metrics.(metricName);
end
end

function [tablesOut, regOut, tname, ycol] = localResolvePromotedMetricView(result, promotedMetricID, dataView)
tablesOut = struct();
regOut = table();
tname = '';
ycol = '';

if ~(isfield(result, 'metrics') && isstruct(result.metrics) && ...
        isfield(result.metrics, 'registry') && istable(result.metrics.registry) && ...
        isfield(result.metrics, 'tables') && isstruct(result.metrics.tables))
    return
end

regOut = result.metrics.registry;
hit = find(string(regOut.metric_id) == promotedMetricID, 1, 'first');
if isempty(hit)
    return
end

tablesOut = result.metrics.tables;
rawCol = char(string(regOut.value_column_raw(hit)));
normCol = char(string(regOut.value_column_norm(hit)));
if strcmp(dataView, 'recommended')
    tr = bwaPreferredTransformForMetricID(char(promotedMetricID), 'cell');
    if isfield(result, 'metrics_by_transform') && isstruct(result.metrics_by_transform)
        field = matlab.lang.makeValidName(tr);
        if isfield(result.metrics_by_transform, field) && isfield(result.metrics_by_transform.(field), 'tables')
            tablesOut = result.metrics_by_transform.(field).tables;
            if isfield(result.metrics_by_transform.(field), 'registry') && ...
                    istable(result.metrics_by_transform.(field).registry)
                reg2 = result.metrics_by_transform.(field).registry;
                hit2 = find(string(reg2.metric_id) == promotedMetricID, 1, 'first');
                if ~isempty(hit2)
                    regOut = reg2;
                    rawCol = char(string(reg2.value_column_raw(hit2)));
                    normCol = char(string(reg2.value_column_norm(hit2)));
                    hit = hit2;
                end
            end
        end
    end
end

tname = char(string(regOut.source_table(hit)));
ycol = rawCol;
if (strcmp(dataView, 'norm') || strcmp(dataView, 'recommended')) && strlength(string(normCol)) > 0
    if isfield(tablesOut, tname)
        T = tablesOut.(tname);
        if istable(T) && ismember(normCol, T.Properties.VariableNames)
            ycol = normCol;
        end
    end
end
end

function epochTable = localResolveEpochTable(result, dataView)
epochTable = table();
switch dataView
    case 'norm'
        if isfield(result, 'canonical_norm') && isfield(result.canonical_norm, 'frequency') && ...
                isfield(result.canonical_norm.frequency, 'epochTable')
            epochTable = result.canonical_norm.frequency.epochTable;
            return
        end
    case 'recommended'
        if isfield(result, 'canonical_norms') && isstruct(result.canonical_norms)
            tr = bwaPreferredTransformForMetricID('frequency', 'canonical');
            field = matlab.lang.makeValidName(tr);
            if isfield(result.canonical_norms, field) && isfield(result.canonical_norms.(field), 'frequency') && ...
                    isfield(result.canonical_norms.(field).frequency, 'epochTable')
                epochTable = result.canonical_norms.(field).frequency.epochTable;
                return
            end
        end
        if isfield(result, 'canonical_norm') && isfield(result.canonical_norm, 'frequency') && ...
                isfield(result.canonical_norm.frequency, 'epochTable')
            epochTable = result.canonical_norm.frequency.epochTable;
            return
        end
end
if isfield(result, 'canonical') && isfield(result.canonical, 'frequency') && ...
        isfield(result.canonical.frequency, 'epochTable')
    epochTable = result.canonical.frequency.epochTable;
end
end

function out = localResolveDisplayName(metricID)
out = string(metricID);
vocabRow = bwaLookupMetricVocabulary(metricID);
if istable(vocabRow) && ~isempty(vocabRow) && ismember('public_metric_id', vocabRow.Properties.VariableNames)
    val = string(vocabRow.public_metric_id(1));
    if strlength(val) > 0
        out = val;
    end
end
end

function promotedMetricID = localMapCanonicalToPromoted(metricName)
promotedMetricID = bwaCanonicalPromotedMetricID(metricName);
end
