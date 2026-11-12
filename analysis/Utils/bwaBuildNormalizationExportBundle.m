function [bundle, summaryTable, metricTable] = bwaBuildNormalizationExportBundle(result)
%BWABUILDNORMALIZATIONEXPORTBUNDLE Build self-describing normalization export metadata.

[metricIDs, unitLevels] = collectMetricRegistry(result);

rawProv = bwaResolveNormalizationProvenance(result, 'raw', metricIDs, unitLevels);
normProv = bwaResolveNormalizationProvenance(result, 'norm', metricIDs, unitLevels);
recommendedProv = bwaResolveNormalizationProvenance(result, 'recommended', metricIDs, unitLevels);

summaryTable = provenanceSummaryTable(rawProv, normProv, recommendedProv);
metricTable = provenanceMetricTable(rawProv, normProv, recommendedProv);

bundle = struct();
bundle.raw = rawProv;
bundle.norm = normProv;
bundle.recommended = recommendedProv;
bundle.summaryByView = summaryTable;
bundle.metricTransforms = metricTable;
end

function [metricIDs, unitLevels] = collectMetricRegistry(result)
metricIDs = strings(0,1);
unitLevels = strings(0,1);
if ~isstruct(result) || ~isfield(result, 'metrics') || ~isstruct(result.metrics) || ...
        ~isfield(result.metrics, 'registry') || ~istable(result.metrics.registry)
    return
end

reg = result.metrics.registry;
if ~all(ismember({'metric_id','unit_level'}, reg.Properties.VariableNames))
    return
end

metricIDs = string(reg.metric_id(:));
unitLevels = string(reg.unit_level(:));
end

function T = provenanceSummaryTable(rawProv, normProv, recommendedProv)
views = ["raw"; "norm"; "recommended"];
P = {rawProv; normProv; recommendedProv};
n = numel(P);

enabled = false(n,1);
mode = strings(n,1);
defaultTransform = strings(n,1);
baselineEpochIdx = nan(n,1);
baselineEpochName = strings(n,1);
compactLabel = strings(n,1);
summaryLabel = strings(n,1);
uniqueTransforms = strings(n,1);

for i = 1:n
    p = P{i};
    enabled(i) = logical(p.normalizationEnabled);
    mode(i) = string(p.normalizationMode);
    defaultTransform(i) = string(p.defaultTransform);
    baselineEpochIdx(i) = double(p.baselineEpochIdx);
    baselineEpochName(i) = string(p.baselineEpochName);
    compactLabel(i) = string(p.compactLabel);
    summaryLabel(i) = string(p.summaryLabel);
    uniqueTransforms(i) = strjoin(string(p.uniqueTransforms), '|');
end

T = table(views, enabled, mode, defaultTransform, baselineEpochIdx, baselineEpochName, compactLabel, summaryLabel, uniqueTransforms, ...
    'VariableNames', {'requested_view','normalization_enabled','normalization_mode','default_transform', ...
    'baseline_epoch_idx','baseline_epoch_name','compact_label','summary_label','unique_transforms'});
end

function T = provenanceMetricTable(rawProv, normProv, recommendedProv)
rawT = rawProv.metricTable;
normT = normProv.metricTable;
recT = recommendedProv.metricTable;

if isempty(rawT) && isempty(normT) && isempty(recT)
    T = table(strings(0,1), strings(0,1), strings(0,1), strings(0,1), strings(0,1), ...
        'VariableNames', {'metric_id','unit_level','raw_transform','norm_transform','recommended_transform'});
    return
end

metricID = string(rawT.metric_id);
unitLevel = string(rawT.unit_level);
rawTransform = string(rawT.effective_transform);
normTransform = string(normT.effective_transform);
recommendedTransform = string(recT.effective_transform);

T = table(metricID, unitLevel, rawTransform, normTransform, recommendedTransform, ...
    'VariableNames', {'metric_id','unit_level','raw_transform','norm_transform','recommended_transform'});
end
