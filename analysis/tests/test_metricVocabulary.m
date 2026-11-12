function tests = test_metricVocabulary
tests = functiontests(localfunctions);
end

function testCatalogContainsContainmentRename(testCase)
T = bwaMetricVocabularyCatalog(["cell.mean_frac_outside_origin_cell"; "cell.frac_escape_events_gt_threshold"]);
verifyEqual(testCase, string(T.public_metric_id)', ["cell.mean_frac_within_origin_cell" "cell.frac_contained_events"]);
verifyEqual(testCase, string(T.family)', ["containment" "containment"]);
end

function testRegistryGetsFamilyAndLineageMetadata(testCase)
canonicalOut = struct();
canonicalNorm = struct();
ctx = struct();
epochSpec = struct('mode','frames','frameRange',[1 10; 11 20], 'epochNames', {{'baseline','CNO1'}});

E = table( ...
    string(["f1_E1";"f1_E2";"f1_E3";"f1_E4"]), ...
    string(repmat("A1",4,1)), ...
    string(repmat("S1",4,1)), ...
    string(repmat("f1",4,1)), ...
    string(repmat("f1_C1",4,1)), ...
    [1;2;12;13], [2;3;13;14], [10;12;11;13], [0.80;0.90;0.70;0.60], [2;7;1;9], [10;20;30;40], ...
    'VariableNames', {'event_uid','animal_id','slice_num','file_id','primary_cell_uid','t0_frame','t1_frame','perimeter_px','circularity','propagation_distance_px','propagation_angle_deg'});
M = table( ...
    string(["f1_E1";"f1_E2";"f1_E3";"f1_E4"]), ...
    string(repmat("f1_C1",4,1)), ...
    true(4,1), [0.05;0.20;0.15;0.40], [1;2;3;4], [0.95;0.80;0.85;0.60], ...
    'VariableNames', {'event_uid','cell_uid','is_primary_cell','frac_outside_cell','outside_cell_area_npix','cell_overlap_frac_of_event'});
canonicalOut.subset = struct('events', E, 'event_cell_map', M);

S = build_local_metrics_schema(canonicalOut, canonicalNorm, ctx, epochSpec, 'propagationThresholdPx', 5, 'escapeThresholdFrac', 0.10);
reg = S.registry;
row = reg(string(reg.metric_id) == "cell.mean_frac_outside_origin_cell", :);
verifyEqual(testCase, string(row.family), "containment");
verifyEqual(testCase, string(row.public_metric_id), "cell.mean_frac_within_origin_cell");
verifyEqual(testCase, string(row.native_unit), "event");
verifyEqual(testCase, string(row.aggregation_path), "event -> cell");
verifyEqual(testCase, string(row.storage_status), "stored");
verifyEqual(testCase, string(row.storage_level), "cell");
verifyEqual(testCase, string(row.metric_class), "promoted");
verifyEqual(testCase, string(row.default_slice_summary), "dynamic_cell_to_slice");
end

function testDescribeMetricLineageReturnsCatalogRows(testCase)
T = describe_metric_lineage("slice.mean_pairwise_cell_correlation");
verifyEqual(testCase, string(T.public_metric_id), "slice.mean_pairwise_event_rate_correlation");
verifyEqual(testCase, string(T.family), "synchrony");
verifyEqual(testCase, string(T.storage_level), "slice");
verifyEqual(testCase, string(T.default_slice_summary), "native_slice");
end

function testCatalogStorageContractDistinguishesEventAndCellMetrics(testCase)
T = bwaMetricVocabularyCatalog(["area_um2"; "cell.mean_event_area_um2"]);
verifyEqual(testCase, string(T.storage_status)', ["stored" "stored"]);
verifyEqual(testCase, string(T.storage_level)', ["event" "cell"]);
verifyEqual(testCase, string(T.metric_class)', ["native" "promoted"]);
verifyEqual(testCase, string(T.default_slice_summary)', ...
    ["dynamic_event_to_cell_to_slice" "dynamic_cell_to_slice"]);
verifyEqual(testCase, string(T.alternate_slice_summary)', ...
    ["dynamic_event_to_slice" "dynamic_event_to_slice"]);
end

function testPromotedSliceRecruitmentIsDynamicOnlyInContract(testCase)
T = bwaMetricVocabularyCatalog("slice.mean_neighbor_entry_delay_sec");
verifyEqual(testCase, string(T.storage_status), "dynamic_only");
verifyEqual(testCase, string(T.storage_level), "dynamic");
verifyEqual(testCase, string(T.metric_class), "promoted_slice_summary");
verifyEqual(testCase, string(T.default_slice_summary), "dynamic_event_to_slice");
end
