function tests = test_plotCapabilityMessaging
tests = functiontests(localfunctions);
end

function testCapabilityHelperRecognizesCanonicalOnlyResult(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));
R = makeCanonicalOnlyResult(root);
caps = bwaDescribeResultCapabilities(R);
verifyTrue(testCase, caps.hasCanonical);
verifyFalse(testCase, caps.hasLegacyPayload);
verifyFalse(testCase, caps.hasPolarsPayload);
verifyFalse(testCase, caps.canPlotSpider);
verifyFalse(testCase, caps.canPlotFileMetrics);
end

function testDescribeResultReportsCanonicalOnlyCapabilities(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));
R = makeCanonicalOnlyResult(root);
out = describe_local_analysis_result(R);
verifyTrue(testCase, istable(out.table));
verifyTrue(testCase, any(string(out.table.capability) == "canonical tables"));
verifyTrue(testCase, any(string(out.table.capability) == "can plot spider"));
verifyFalse(testCase, out.capabilities.canPlotSpider);
end

function testSpiderWrapperFailsEarlyForCanonicalOnlyResult(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));
R = makeCanonicalOnlyResult(root);
f = @() plot_local_analysis_spider(R);
verifyError(testCase, f, 'plot_local_analysis_spider:MissingPolarsPayload');
try
    f();
catch ME
    verifyTrue(testCase, contains(string(ME.message), "canonical-only"));
    verifyTrue(testCase, contains(string(ME.message), "Canonical summary plots remain available"));
end
end

function testLegacyPanelsFailEarlyForCanonicalOnlyResult(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));
R = makeCanonicalOnlyResult(root);
f = @() plot_local_analysis_legacy_panels(R);
verifyError(testCase, f, 'plot_local_analysis_legacy_panels:MissingStandard');
try
    f();
catch ME
    verifyTrue(testCase, contains(string(ME.message), "legacy-enriched compute run"));
    verifyTrue(testCase, contains(string(ME.message), "Canonical summary plots remain available"));
end
end

function testFileMetricsFailEarlyForCanonicalOnlyResult(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));
R = makeCanonicalOnlyResult(root);
f = @() plot_local_analysis_file_metrics(R, 'basic_area');
verifyError(testCase, f, 'plot_local_analysis_file_metrics:MissingLegacyPayload');
try
    f();
catch ME
    verifyTrue(testCase, contains(string(ME.message), "requires BWA_AQ2.mat"));
    verifyTrue(testCase, contains(string(ME.message), "Canonical summary plots remain available"));
end
end

function R = makeCanonicalOnlyResult(root)
makeCanonicalOnlyRun(root);
epochTable = table([1;2], ["baseline";"drug"], 'VariableNames', {'epoch_idx','epoch_name'});
summaryByGroupEpoch = table( ...
    ["all";"all"], [1;2], ["baseline";"drug"], [1.0; 1.5], ...
    'VariableNames', {'group_id','epoch_idx','epoch_name','mean_rate_events_per_min'});
R = struct();
R.pathIn = root;
R.canonical = struct('epochTable', epochTable, ...
    'frequency', struct('summaryByGroupEpoch', summaryByGroupEpoch));
R.derived = struct('status', 'skipped_missing_legacy', ...
    'error', 'BWA_AQ2.mat not available; skipped legacy-derived metrics.', ...
    'context', struct('timeFrameRangeEpochs', struct('epochNames', {{'baseline','drug'}}, 'frameRange', [1 10;11 20], 'timeRange', [0 9;10 19])));
end

function makeCanonicalOnlyRun(root)
runDir = fullfile(root, 'derived', 'ingest_runs', '20260224_000101');
mkdir(runDir);

files = table(string("f1"), string("a1"), 'VariableNames', {'file_id','animal_id'});
cells = table(string("f1_C1"), string("f1"), 'VariableNames', {'cell_uid','file_id'});
events = table(string("f1_E1"), string("f1"), 1, 10, 12, 'VariableNames', {'event_uid','file_id','local_event_id','t0_frame','t1_frame'});
event_cell_map = table(string("f1_E1"), string("f1_C1"), true, 'VariableNames', {'event_uid','cell_uid','is_primary_cell'});

save(fullfile(runDir, 'files.mat'), 'files');
save(fullfile(runDir, 'cells.mat'), 'cells');
save(fullfile(runDir, 'events.mat'), 'events');
save(fullfile(runDir, 'event_cell_map.mat'), 'event_cell_map');

fid = fopen(fullfile(root, 'derived', 'ingest_runs', 'latest_ingest.txt'), 'w');
fprintf(fid, '%s\n', runDir);
fclose(fid);
end

function safeRmdir(pathIn)
if exist(pathIn, 'dir')
    rmdir(pathIn, 's');
end
end
