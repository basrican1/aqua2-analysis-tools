function [mOut, ok] = bwaBuildDynamicSliceMetricOut(result, metricID, dataView)
% Build a dynamic slice-level metric output for promoted slice summaries.

mOut = struct();
ok = false;
metricID = string(metricID);
dataView = lower(char(string(dataView)));

spec = dynamicSliceMetricSpec(metricID);
if isempty(spec)
    return
end

transform = resolveTransform(result, dataView, metricID);
T = bwaBuildDynamicSliceRecruitmentTable(result, 'transform', transform);
if ~istable(T) || isempty(T)
    return
end

ycol = spec.rawCol;
if (strcmp(dataView, 'norm') || strcmp(dataView, 'recommended')) && ismember(spec.normCol, T.Properties.VariableNames)
    ycol = spec.normCol;
end
if ~ismember(ycol, T.Properties.VariableNames)
    return
end

epochTable = table();
if isfield(result, 'canonical') && isfield(result.canonical, 'frequency') && ...
        isfield(result.canonical.frequency, 'epochTable')
    epochTable = result.canonical.frequency.epochTable;
end
if isempty(epochTable) && ismember('epoch_idx', T.Properties.VariableNames) && ismember('epoch_name', T.Properties.VariableNames)
    epochTable = unique(T(:, {'epoch_idx','epoch_name'}), 'rows', 'stable');
end
if isempty(epochTable)
    return
end

mOut = bwaBuildMetricOutFromTable(T, ycol, epochTable);
vocabRow = bwaLookupMetricVocabulary(metricID);
if istable(vocabRow) && ~isempty(vocabRow) && ismember('public_metric_id', vocabRow.Properties.VariableNames)
    val = string(vocabRow.public_metric_id(1));
    if strlength(val) > 0
        mOut.titleName = char(val);
    end
end
ok = true;
end

function spec = dynamicSliceMetricSpec(metricID)
spec = struct([]);
pairs = { ...
    'slice.frac_primary_events_with_neighbor_recruitment','frac_primary_events_with_neighbor_recruitment','frac_primary_events_with_neighbor_recruitment_norm'; ...
    'slice.mean_neighbor_entry_delay_sec','mean_neighbor_entry_delay_sec','mean_neighbor_entry_delay_sec_norm'; ...
    'slice.mean_neighbor_peak_overlap_frac_of_cell','mean_neighbor_peak_overlap_frac_of_cell','mean_neighbor_peak_overlap_frac_of_cell_norm'; ...
    'slice.mean_neighbor_partitions_per_event','mean_neighbor_partitions_per_event','mean_neighbor_partitions_per_event_norm'};
for i = 1:size(pairs,1)
    if string(pairs{i,1}) == metricID
        spec = struct('rawCol', char(pairs{i,2}), 'normCol', char(pairs{i,3}));
        return
    end
end
end

function transform = resolveTransform(result, dataView, metricID)
switch dataView
    case 'recommended'
        transform = bwaPreferredTransformForMetricID(metricID, 'slice');
    case 'norm'
        transform = 'delta';
        if isfield(result, 'normalization') && isstruct(result.normalization) && ...
                isfield(result.normalization, 'transform')
            transform = lower(char(string(result.normalization.transform)));
        end
    otherwise
        transform = 'delta';
end
end
