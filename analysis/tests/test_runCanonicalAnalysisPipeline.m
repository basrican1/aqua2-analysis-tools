function tests = test_runCanonicalAnalysisPipeline
% Unit tests for runCanonicalAnalysisPipeline

tests = functiontests(localfunctions);
end

function testPipelineSmokeWithSubset(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));

runDir = makeMiniRun(root);
ep = struct();
ep.timeRange = [0 1; 1 2];
ep.epochNames = {'E1';'E2'};

out = runCanonicalAnalysisPipeline(root, ep, ...
    'subset', struct('genotypes', {{'WT'}}), ...
    'metrics', {'area_px2','duration_frames','does_not_exist'}, ...
    'makePlots', false);

testCase.verifyTrue(isfield(out, 'frequency'));
testCase.verifyTrue(isfield(out, 'metrics'));
testCase.verifyTrue(isfield(out.metrics, 'area_px2'));
testCase.verifyTrue(isfield(out.metrics, 'duration_frames'));
testCase.verifyFalse(isfield(out.metrics, 'does_not_exist'));
testCase.verifyEqual(unique(string(out.subset.files.genotype)), "WT");
testCase.verifyEqual(out.loaded.runDir, runDir);
end

function runDir = makeMiniRun(root)
runDir = fullfile(root, 'derived', 'ingest_runs', '20260224_000000');
mkdir(runDir);

BWA_AQ2 = struct();
BWA_AQ2.fts1 = struct();
save(fullfile(root, 'BWA_AQ2.mat'), 'BWA_AQ2');

files = table( ...
    string(["f1";"f2"]), ...
    string(["a1";"a2"]), ...
    string(["WT";"KO"]), ...
    'VariableNames', {'file_id','animal_id','genotype'});

cells = table( ...
    string(["f1_C1";"f2_C1"]), ...
    string(["f1";"f2"]), ...
    'VariableNames', {'cell_uid','file_id'});

events = table( ...
    string(["f1_E1";"f1_E2";"f2_E1"]), ...
    string(["f1";"f1";"f2"]), ...
    [1;2;1], ...
    [10;70;10], ...
    [12;72;12], ...
    [10;20;30], ...
    [5;6;7], ...
    'VariableNames', {'event_uid','file_id','local_event_id','t0_sec','t1_sec','area_px2','duration_frames'});

event_cell_map = table( ...
    string(["f1_E1";"f1_E2";"f2_E1"]), ...
    string(["f1_C1";"f1_C1";"f2_C1"]), ...
    true(3,1), ...
    'VariableNames', {'event_uid','cell_uid','is_primary_cell'});

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
