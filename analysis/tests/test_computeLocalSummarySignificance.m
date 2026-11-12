function tests = test_computeLocalSummarySignificance
 tests = functiontests(localfunctions);
end

function testComputeFromStructIncludesCanonicalAndCellMetrics(testCase)
result = makeSyntheticResult();
out = compute_local_summary_significance(result, 'comparisonLevel', 'by_slice');

verifyTrue(testCase, isfield(out, 'frequency'));
verifyTrue(testCase, isfield(out, 'metrics'));
verifyTrue(testCase, isfield(out.metrics, 'area_px2'));
verifyTrue(testCase, isfield(out, 'cell_metrics'));
verifyFalse(testCase, isfield(out.cell_metrics, 'cell_frequency'));
verifyTrue(testCase, istable(out.frequency.significance));
verifyTrue(testCase, istable(out.metrics.area_px2.significance));
verifyEqual(testCase, string(out.metrics.area_px2.significance.epoch_name(2)), "drug");
verifyTrue(testCase, any(out.metrics.area_px2.significance.is_sig(2:end)));
verifyTrue(testCase, isfield(out, 'normalizationProvenance'));
verifyEqual(testCase, string(out.normalizationProvenance.compactLabel), "norm=raw");
verifyEqual(testCase, string(out.canonicalAggregation), "cell_promoted");
end

function testNormViewUsesCanonicalNorm(testCase)
result = makeSyntheticResult();
out = compute_local_summary_significance(result, 'dataView', 'norm', 'comparisonLevel', 'by_slice');
verifyEqual(testCase, out.dataView, 'norm');
verifyTrue(testCase, any(out.metrics.area_px2.significance.is_sig(2:end)));
verifyEqual(testCase, string(out.normalizationProvenance.compactLabel), "norm=default [log2_fold]");
end

function testMetricListFiltersCanonicalOnly(testCase)
result = makeSyntheticResult();
out = compute_local_summary_significance(result, 'metricList', {'duration_frames'}, 'includeDerived', false);
verifyFalse(testCase, isfield(out.metrics, 'area_px2'));
verifyTrue(testCase, isfield(out.metrics, 'duration_frames'));
end

function testMetricFamiliesFilterCanonicalAndCellMetrics(testCase)
result = makeSyntheticResult();
out = compute_local_summary_significance(result, ...
    'metricFamilies', {'shape'}, ...
    'includeDerived', false, ...
    'comparisonLevel', 'by_slice');

verifyTrue(testCase, isfield(out.metrics, 'area_px2'));
verifyFalse(testCase, isfield(out.metrics, 'duration_frames'));
verifyFalse(testCase, isfield(out.cell_metrics, 'cell_mean_event_area_um2'));
verifyFalse(testCase, isfield(out, 'frequency') && isfield(out.frequency, 'significance'));
end

function testEventNativeRetainsPromotedCellMetricSurface(testCase)
result = makeSyntheticResult();
out = compute_local_summary_significance(result, ...
    'metricFamilies', {'shape'}, ...
    'includeDerived', false, ...
    'comparisonLevel', 'by_slice', ...
    'canonicalAggregation', 'event_native');

verifyTrue(testCase, isfield(out.metrics, 'area_px2'));
verifyTrue(testCase, isfield(out.cell_metrics, 'cell_mean_event_area_um2'));
end

function testCanonicalDefaultsToCellPromotedWhenAvailable(testCase)
result = makeSyntheticResult();
out = compute_local_summary_significance(result, 'metricList', {'area_um2'}, 'includeDerived', false, 'comparisonLevel', 'by_slice');

verifyTrue(testCase, isfield(out.metrics, 'area_um2'));
verifyEqual(testCase, string(out.metrics.area_um2.effectiveAggregation), "cell_promoted");
verifyEqual(testCase, string(out.metrics.area_um2.promoted_metric_id), "cell.mean_event_area_um2");
verifyEqual(testCase, string(out.metrics.area_um2.significance.epoch_name(2)), "drug");
drugVals = double(out.metrics.area_um2.traceTable.value(out.metrics.area_um2.traceTable.epoch_idx == 2));
verifyEqual(testCase, mean(drugVals), 2.5, 'AbsTol', 1e-9);
end

function testCanonicalCanUseEventNativeWhenRequested(testCase)
result = makeSyntheticResult();
out = compute_local_summary_significance(result, ...
    'metricList', {'area_um2'}, ...
    'includeDerived', false, ...
    'comparisonLevel', 'by_slice', ...
    'canonicalAggregation', 'event_native');

verifyEqual(testCase, string(out.canonicalAggregation), "event_native");
verifyEqual(testCase, string(out.metrics.area_um2.effectiveAggregation), "event_native");
verifyEqual(testCase, strlength(string(out.metrics.area_um2.promoted_metric_id)), 0);
drugVals = double(out.metrics.area_um2.traceTable.value(out.metrics.area_um2.traceTable.epoch_idx == 2));
verifyEqual(testCase, mean(drugVals), 3.5, 'AbsTol', 1e-9);
end

function testCanonicalSelectorLabelsAggregationMode(testCase)
result = makeSyntheticResult();

[mCell, infoCell] = bwaSelectCanonicalMetricSummary(result, 'area_um2');
verifyEqual(testCase, string(infoCell.effectiveAggregation), "cell_promoted");
verifyTrue(testCase, contains(string(mCell.titleName), "cell-promoted"));

[mEvent, infoEvent] = bwaSelectCanonicalMetricSummary(result, 'area_um2', 'canonicalAggregation', 'event_native');
verifyEqual(testCase, string(infoEvent.effectiveAggregation), "event_native");
verifyTrue(testCase, contains(string(mEvent.titleName), "event-native"));
end

function result = makeSyntheticResult()
epochTable = table([1;2],["baseline";"drug"], 'VariableNames', {'epoch_idx','epoch_name'});
files = table(["f1";"f2"],["a1";"a2"],[1;1], 'VariableNames', {'file_id','animal_id','slice_num'});

freqSummary = table(["a1";"a1";"a2";"a2"], [1;2;1;2], ["baseline";"drug";"baseline";"drug"], ...
    ["f1";"f1";"f2";"f2"], [1;1;1;1], [1;4;1;5], [1;4;1;5], [0;0;0;0], ...
    'VariableNames', {'animal_id','epoch_idx','epoch_name','file_id','slice_num', ...
    'mean_rate_events_per_min','median_rate_events_per_min','sem_rate_events_per_min'});
ratesByCellEpoch = table(["a1";"a1";"a2";"a2"], [1;2;1;2], ["baseline";"drug";"baseline";"drug"], ...
    ["f1";"f1";"f2";"f2"], [1;1;1;1], ["c1";"c1";"c2";"c2"], [1;4;1;5], ...
    'VariableNames', {'animal_id','epoch_idx','epoch_name','file_id','slice_num','cell_uid','rate_events_per_min'});
frequency = struct('summaryByGroupEpoch', freqSummary, 'ratesByCellEpoch', ratesByCellEpoch, 'epochTable', epochTable);

areaSummary = table(["a1";"a1";"a2";"a2"], [1;2;1;2], ["baseline";"drug";"baseline";"drug"], ...
    ["f1";"f1";"f2";"f2"], [1;1;1;1], [1;3;1;4], [1;3;1;4], [0.1;0.1;0.1;0.1], ...
    'VariableNames', {'animal_id','epoch_idx','epoch_name','file_id','slice_num', ...
    'mean_area_px2','median_area_px2','sem'});
areaEvents = table(["f1";"f1";"f2";"f2"], [1;2;1;2], [1;3;1;4], ["c1";"c1";"c2";"c2"], ["a1";"a1";"a2";"a2"], [1;1;1;1], ["baseline";"drug";"baseline";"drug"], ...
    'VariableNames', {'file_id','epoch_idx','area_px2','primary_cell_uid','animal_id','slice_num','epoch_name'});
areaMetric = struct('metricCol', 'area_px2', 'summaryByGroupEpoch', areaSummary, 'eventValues', areaEvents, 'epochTable', epochTable);

areaUmSummary = areaSummary;
areaUmSummary.Properties.VariableNames{'mean_area_px2' == string(areaUmSummary.Properties.VariableNames)} = 'mean_area_um2';
areaUmSummary.Properties.VariableNames{'median_area_px2' == string(areaUmSummary.Properties.VariableNames)} = 'median_area_um2';
areaUmEvents = areaEvents;
areaUmEvents.area_um2 = [1;3;1;4];
areaUmMetric = struct('metricCol', 'area_um2', 'summaryByGroupEpoch', areaUmSummary, 'eventValues', areaUmEvents, 'epochTable', epochTable);

durSummary = areaSummary;
durSummary.mean_area_px2 = [2;2;2;2];
durSummary.median_area_px2 = [2;2;2;2];
durSummary.Properties.VariableNames{'mean_area_px2' == string(durSummary.Properties.VariableNames)} = 'mean_duration_frames';
durSummary.Properties.VariableNames{'median_area_px2' == string(durSummary.Properties.VariableNames)} = 'median_duration_frames';
durEvents = areaEvents;
durEvents.duration_frames = [2;2;2;2];
durEvents = removevars(durEvents, 'area_px2');
durationMetric = struct('metricCol', 'duration_frames', 'summaryByGroupEpoch', durSummary, 'eventValues', durEvents, 'epochTable', epochTable);

canonical = struct();
canonical.frequency = frequency;
canonical.metrics = struct('area_px2', areaMetric, 'area_um2', areaUmMetric, 'duration_frames', durationMetric);
canonical.subset = struct('files', files);

areaSummaryN = areaSummary;
areaSummaryN.mean_area_px2 = [0;1.2;0;1.5];
areaSummaryN.median_area_px2 = areaSummaryN.mean_area_px2;
areaEventsN = areaEvents;
areaEventsN.area_px2 = [0;1.2;0;1.5];
areaMetricN = struct('metricCol', 'area_px2', 'summaryByGroupEpoch', areaSummaryN, 'eventValues', areaEventsN, 'epochTable', epochTable);

areaSummaryUmN = areaUmSummary;
areaSummaryUmN.mean_area_um2 = [0;1.2;0;1.5];
areaSummaryUmN.median_area_um2 = areaSummaryUmN.mean_area_um2;
areaEventsUmN = areaUmEvents;
areaEventsUmN.area_um2 = [0;1.2;0;1.5];
areaMetricUmN = struct('metricCol', 'area_um2', 'summaryByGroupEpoch', areaSummaryUmN, 'eventValues', areaEventsUmN, 'epochTable', epochTable);

freqSummaryN = freqSummary;
freqSummaryN.mean_rate_events_per_min = [0;1.5;0;2.0];
freqSummaryN.median_rate_events_per_min = freqSummaryN.mean_rate_events_per_min;
ratesByCellEpochN = ratesByCellEpoch;
ratesByCellEpochN.rate_events_per_min = [0;1.5;0;2.0];
frequencyN = struct('summaryByGroupEpoch', freqSummaryN, 'ratesByCellEpoch', ratesByCellEpochN, 'epochTable', epochTable);

canonical_norm = canonical;
canonical_norm.frequency = frequencyN;
canonical_norm.metrics.area_px2 = areaMetricN;
canonical_norm.metrics.area_um2 = areaMetricUmN;

cellTbl = table(["a1";"a1";"a2";"a2"], [1;2;1;2], ["baseline";"drug";"baseline";"drug"], ...
    ["f1";"f1";"f2";"f2"], [1;1;1;1], ["c1";"c1";"c2";"c2"], [1;4;1;5], [0;2;0;2], ...
    'VariableNames', {'animal_id','epoch_idx','epoch_name','file_id','slice_num','primary_cell_uid','rate_events_per_min','rate_events_per_min_norm'});
cellCanonicalTbl = table(["a1";"a1";"a2";"a2"], [1;2;1;2], ["baseline";"drug";"baseline";"drug"], ...
    ["f1";"f1";"f2";"f2"], [1;1;1;1], ["c1";"c1";"c2";"c2"], [1;2;1;3], [0;1;0;2], ...
    'VariableNames', {'animal_id','epoch_idx','epoch_name','file_id','slice_num','primary_cell_uid','mean_event_area_um2','mean_event_area_um2_norm'});
registry = table( ...
    ["cell.frequency";"cell.mean_event_area_um2"], ...
    ["cell";"cell"], ...
    ["cell_metrics";"cell_canonical"], ...
    ["rate_events_per_min";"mean_event_area_um2"], ...
    ["rate_events_per_min_norm";"mean_event_area_um2_norm"], ...
    'VariableNames', {'metric_id','unit_level','source_table','value_column_raw','value_column_norm'});
metrics = struct('registry', registry, 'tables', struct('cell_metrics', cellTbl, 'cell_canonical', cellCanonicalTbl));

metrics_by_transform = struct();
metrics_by_transform.delta = struct('registry', registry, 'tables', struct('cell_metrics', cellTbl, 'cell_canonical', cellCanonicalTbl));
metrics_by_transform.log2_fold = struct('registry', registry, 'tables', struct('cell_metrics', cellTbl, 'cell_canonical', cellCanonicalTbl));

result = struct();
result.canonical = canonical;
result.canonical_norm = canonical_norm;
result.metrics = metrics;
result.metrics_by_transform = metrics_by_transform;
result.normalization = struct('mode', 'by_cell_baseline_mean', 'transform', 'log2_fold', 'baselineEpochIdx', 1);
result.derived = struct();
end
