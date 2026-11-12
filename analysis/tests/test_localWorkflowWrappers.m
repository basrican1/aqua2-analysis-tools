function tests = test_localWorkflowWrappers
% Wrapper-level tests for local workflow helpers.

tests = functiontests(localfunctions);
end

function testLoadResultAcceptsCanonicalStruct(testCase)
in = struct();
in.canonical = struct('frequency', struct(), 'metrics', struct());
out = load_local_analysis_result(in);
testCase.verifyTrue(isfield(out, 'canonical'));
end

function testLoadClusterOutputsAddsTreatmentFallback(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));
makeMiniRun(root);

out = loadClusterOutputs(root, 'attachLegacyMeta', false, 'attachEventMeta', false);
testCase.verifyTrue(ismember('treatment', out.files.Properties.VariableNames));
testCase.verifyEqual(string(out.files.treatment(1)), "CNO");
testCase.verifyEqual(string(out.files.treatment(2)), "VirusX");
end

function testComputeRespectsModeAndGroupValues(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));
makeMiniRun(root);
makeEpochConfig(root);

out = run_local_analysis_compute( ...
    'pathIn', root, ...
    'mode', 'genotype', ...
    'groupValues', {'WT'}, ...
    'useQC', false, ...
    'exportOutputs', false, ...
    'runStats', false);

testCase.verifyEqual(unique(string(out.canonical.subset.files.genotype)), "WT");
testCase.verifyEqual(out.runMode, 'genotype');
testCase.verifyTrue(isfield(out.derived, 'standard'));
testCase.verifyTrue(isfield(out.derived, 'polars'));
testCase.verifyEqual(string(out.normalization.mode), "by_cell_baseline_mean");
testCase.verifyEqual(string(out.normalized.default_transform), "delta");
end

function testDefaultComputeEnablesBaselineNormalization(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));
makeMiniRun(root);
makeEpochConfig(root);

out = run_local_analysis_compute( ...
    'pathIn', root, ...
    'mode', 'epoch', ...
    'useQC', false, ...
    'exportOutputs', false, ...
    'runStats', false);

testCase.verifyEqual(string(out.normalization.mode), "by_cell_baseline_mean");
testCase.verifyTrue(logical(out.normalization.enabled));
testCase.verifyEqual(string(out.normalization.transform), "delta");
testCase.verifyEqual(string(out.normalized.default_transform), "delta");
end

function testCanonicalOnlyComputeSkipsLegacyDerivedMetrics(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));
makeMiniRun(root, false);
makeEpochConfig(root);

out = run_local_analysis_compute( ...
    'pathIn', root, ...
    'mode', 'epoch', ...
    'useQC', false, ...
    'exportOutputs', false, ...
    'runStats', false);

testCase.verifyFalse(isfield(out.derived.context, 'BWA_AQ2'));
testCase.verifyEqual(string(out.derived.status), "skipped_missing_legacy");
testCase.verifyEqual(string(out.derived.error), "BWA_AQ2.mat not available; skipped legacy-derived metrics.");
testCase.verifyEqual(height(out.canonical.frequency.summaryByGroupEpoch), 2);
end

function testExportedResultCarriesNormalizationProvenance(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));
makeMiniRun(root);
makeEpochConfig(root);

exportDir = fullfile(root, 'exports');
out = run_local_analysis_compute( ...
    'pathIn', root, ...
    'mode', 'epoch', ...
    'useQC', false, ...
    'runStats', false, ...
    'exportOutputs', true, ...
    'exportSlim', true, ...
    'exportDir', exportDir);

testCase.verifyTrue(exist(fullfile(exportDir, 'analysis_result.mat'), 'file') == 2);
testCase.verifyTrue(exist(fullfile(exportDir, 'normalization_provenance.mat'), 'file') == 2);
testCase.verifyTrue(exist(fullfile(exportDir, 'normalization_provenance.csv'), 'file') == 2);
testCase.verifyTrue(exist(fullfile(exportDir, 'normalization_metric_provenance.csv'), 'file') == 2);
testCase.verifyTrue(exist(fullfile(exportDir, 'result_capabilities.mat'), 'file') == 2);
testCase.verifyTrue(exist(fullfile(exportDir, 'result_capabilities.csv'), 'file') == 2);
testCase.verifyTrue(exist(fullfile(exportDir, 'metric_vocabulary.mat'), 'file') == 2);
testCase.verifyTrue(exist(fullfile(exportDir, 'metric_vocabulary.csv'), 'file') == 2);
testCase.verifyTrue(isfield(out, 'exportMeta'));
testCase.verifyTrue(isfield(out.exportMeta, 'normalization'));
testCase.verifyTrue(isfield(out.exportMeta, 'capabilities'));
testCase.verifyEqual(string(out.exportMeta.normalization.norm.defaultTransform), "delta");
testCase.verifyTrue(logical(out.exportMeta.capabilities.hasLegacyPayload));

S = load(fullfile(exportDir, 'analysis_result.mat'), 'result');
testCase.verifyTrue(isfield(S.result, 'exportMeta'));
testCase.verifyTrue(isfield(S.result.exportMeta, 'normalization'));
testCase.verifyTrue(isfield(S.result.exportMeta, 'capabilities'));
testCase.verifyTrue(isfield(S.result.exportMeta, 'metricVocabularyRows'));
testCase.verifyTrue(logical(S.result.exportMeta.slim));
testCase.verifyEqual(string(S.result.exportMeta.normalization.recommended.compactLabel), "norm=delta");
testCase.verifyTrue(logical(S.result.exportMeta.capabilities.canPlotSpider));

T = readtable(fullfile(exportDir, 'normalization_provenance.csv'));
testCase.verifyTrue(any(string(T.requested_view) == "recommended"));

Tm = readtable(fullfile(exportDir, 'normalization_metric_provenance.csv'));
testCase.verifyTrue(any(string(Tm.metric_id) == "cell.frequency"));
testCase.verifyTrue(any(string(Tm.recommended_transform) == "delta"));

Tc = readtable(fullfile(exportDir, 'result_capabilities.csv'));
testCase.verifyTrue(any(string(Tc.capability) == "canPlotSpider" & logical(Tc.available)));
testCase.verifyTrue(any(string(Tc.capability) == "hasLegacyPayload" & logical(Tc.available)));

Tv = readtable(fullfile(exportDir, 'metric_vocabulary.csv'));
testCase.verifyTrue(any(string(Tv.family) == "activity"));
testCase.verifyTrue(any(string(Tv.public_metric_id) == "cell.event_rate_per_min"));
end

function testPresetEpochFallbackToDefaultConfigPreset(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));
makeMiniRun(root);
makeEpochConfig(root);
makeAnalysisPlanWithMissingEpochPreset(root);

out = run_local_analysis_compute( ...
    'pathIn', root, ...
    'mode', 'preset', ...
    'presetName', 'drug_epoch_compare', ...
    'useQC', false, ...
    'exportOutputs', false, ...
    'runStats', false);

testCase.verifyTrue(isfield(out.canonical, 'epochMeta'));
if isfield(out.canonical.epochMeta, 'fallbackUsed')
    testCase.verifyTrue(logical(out.canonical.epochMeta.fallbackUsed));
end
end

function makeMiniRun(root, includeLegacy)
if nargin < 2
    includeLegacy = true;
end
runDir = fullfile(root, 'derived', 'ingest_runs', '20260224_000003');
mkdir(runDir);

if includeLegacy
    % Minimal BWA payload: enough for loader; derived compute may skip.
    BWA_AQ2 = struct();
    BWA_AQ2.fts1 = struct();
    save(fullfile(root, 'BWA_AQ2.mat'), 'BWA_AQ2');
end

files = table( ...
    string(["f1";"f2"]), ...
    string(["a1";"a2"]), ...
    string(["WT";"KO"]), ...
    string(["CNO";""]), ...
    string(["";"VirusX"]), ...
    [80;80], ...
    [1;1], ...
    'VariableNames', {'file_id','animal_id','genotype','drug','virus','sz_t','frame_rate_s'});

cells = table( ...
    string(["f1_C1";"f2_C1"]), ...
    string(["f1";"f2"]), ...
    'VariableNames', {'cell_uid','file_id'});

events = table( ...
    string(["f1_E1";"f2_E1"]), ...
    string(["f1";"f2"]), ...
    [1;1], ...
    [10;10], ...
    [12;12], ...
    [10;20], ...
    [5;6], ...
    'VariableNames', {'event_uid','file_id','local_event_id','t0_sec','t1_sec','area_px2','duration_frames'});

event_cell_map = table( ...
    string(["f1_E1";"f2_E1"]), ...
    string(["f1_C1";"f2_C1"]), ...
    true(2,1), ...
    'VariableNames', {'event_uid','cell_uid','is_primary_cell'});

save(fullfile(runDir, 'files.mat'), 'files');
save(fullfile(runDir, 'cells.mat'), 'cells');
save(fullfile(runDir, 'events.mat'), 'events');
save(fullfile(runDir, 'event_cell_map.mat'), 'event_cell_map');

fid = fopen(fullfile(root, 'derived', 'ingest_runs', 'latest_ingest.txt'), 'w');
fprintf(fid, '%s\n', runDir);
fclose(fid);
end

function makeEpochConfig(root)
cfg = struct();
cfg.defaultPreset = 'all_epoch';
cfg.presets = struct();
cfg.presets.all_epoch = struct('epochs', struct('mode', 'all', 'epochNames', {{'All'}}));
fid = fopen(fullfile(root, 'analysisConfig.json'), 'w');
fwrite(fid, jsonencode(cfg, 'PrettyPrint', true), 'char');
fclose(fid);
end

function makeAnalysisPlanWithMissingEpochPreset(root)
plan = struct();
plan.defaultPreset = 'drug_epoch_compare';
plan.presets = struct();
plan.presets.drug_epoch_compare = struct( ...
    'epochPreset', 'drug_application', ... % intentionally missing from analysisConfig
    'subset', struct(), ...
    'metrics', {{'area_px2','duration_frames'}}, ...
    'freqGroupBy', {{'file_id'}}, ...
    'metricGroupBy', {{'file_id'}}, ...
    'areaPercentiles', [0 100], ...
    'makePlots', false, ...
    'strictQC', false);
fid = fopen(fullfile(root, 'analysisPlan.json'), 'w');
fwrite(fid, jsonencode(plan, 'PrettyPrint', true), 'char');
fclose(fid);
end

function safeRmdir(pathIn)
if exist(pathIn, 'dir')
    rmdir(pathIn, 's');
end
end
