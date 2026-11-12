function tests = test_flattenLocalSummarySignificance
 tests = functiontests(localfunctions);
end

function testFlattenFromComputedStruct(testCase)
S = makeSyntheticSigStruct();
T = flatten_local_summary_significance(S);

verifyTrue(testCase, istable(T));
verifyTrue(testCase, height(T) >= 4);
verifyTrue(testCase, all(ismember({'family','metric','epoch_idx','epoch_name','p_value','is_sig','method_used'}, T.Properties.VariableNames)));
verifyTrue(testCase, any(T.family == "frequency"));
verifyTrue(testCase, any(T.family == "canonical_metric" & T.metric == "area_px2"));
verifyFalse(testCase, any(T.family == "cell_metric" & T.metric == "cell_frequency"));
verifyTrue(testCase, any(T.family == "slice_metric" & T.metric == "slice_frac_active_cells"));
end

function testFlattenFromResultRunsCompute(testCase)
result = makeSyntheticResult();
T = flatten_local_summary_significance(result, 'comparisonLevel', 'by_slice');
verifyTrue(testCase, any(T.metric == "area_px2"));
verifyTrue(testCase, any(T.metric == "frequency"));
end

function S = makeSyntheticSigStruct()
sig = table([1;2],["baseline";"drug"],[2;2],[NaN;0.01],[false;true],["baseline";"paired_signrank"], ...
    'VariableNames', {'epoch_idx','epoch_name','n','p_value','is_sig','method_used'});
S = struct();
S.source = 'compute_local_summary_significance';
S.frequency = struct('valueVar', 'mean_rate_events_per_min', 'significance', sig);
S.metrics = struct('area_px2', struct('valueVar', 'mean_area_px2', 'significance', sig));
S.derived = struct('derived_prop_distance', struct('valueVar', 'mean_value', 'significance', sig));
S.cell_metrics = struct('cell_frequency', struct('valueVar', 'mean_rate_events_per_min', 'significance', sig));
S.slice_metrics = struct('slice_frac_active_cells', struct('valueVar', 'mean_frac_active_cells', 'significance', sig));
end

function result = makeSyntheticResult()
epochTable = table([1;2],["baseline";"drug"], 'VariableNames', {'epoch_idx','epoch_name'});
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

canonical = struct('frequency', frequency, 'metrics', struct('area_px2', areaMetric));
cellTbl = table(["a1";"a1";"a2";"a2"], [1;2;1;2], ["baseline";"drug";"baseline";"drug"], ...
    ["f1";"f1";"f2";"f2"], [1;1;1;1], ["c1";"c1";"c2";"c2"], [1;4;1;5], [0;2;0;2], ...
    'VariableNames', {'animal_id','epoch_idx','epoch_name','file_id','slice_num','primary_cell_uid','rate_events_per_min','rate_events_per_min_norm'});
registry = table("cell.frequency", "cell", "cell_metrics", "rate_events_per_min", "rate_events_per_min_norm", ...
    'VariableNames', {'metric_id','unit_level','source_table','value_column_raw','value_column_norm'});
metrics = struct('registry', registry, 'tables', struct('cell_metrics', cellTbl));

result = struct('canonical', canonical, 'metrics', metrics, 'derived', struct());
end
