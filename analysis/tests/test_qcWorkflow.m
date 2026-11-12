function tests = test_qcWorkflow
% Unit tests for QC workflow helpers.

tests = functiontests(localfunctions);
end

function testCreateAndValidateQCTables(testCase)
[files, cells, events, ~] = miniCanonical();

qc = createQCTables(files, cells, events, 'defaultKeep', true, 'reviewer', 'tester');
testCase.verifyEqual(height(qc.qc_files), height(files));
testCase.verifyEqual(height(qc.qc_cells), height(cells));
testCase.verifyEqual(height(qc.qc_events), height(events));
testCase.verifyTrue(all(qc.qc_files.keep));

rep = validateQCTables(qc.qc_files, qc.qc_cells, qc.qc_events, files, cells, events);
testCase.verifyTrue(rep.ok);
end

function testApplyQCToCanonical(testCase)
[files, cells, events, map] = miniCanonical();
qc = createQCTables(files, cells, events);

% Drop file f2 and one event from f1.
qc.qc_files.keep(qc.qc_files.file_id == "f2") = false;
qc.qc_files.reason(qc.qc_files.file_id == "f2") = "exclude_file";
qc.qc_events.keep(qc.qc_events.event_uid == "f1_E2") = false;

out = applyQCToCanonical(files, cells, events, map, qc.qc_files, qc.qc_cells, qc.qc_events);

testCase.verifyEqual(unique(string(out.files.file_id)), "f1");
testCase.verifyEqual(unique(string(out.cells.file_id)), "f1");
testCase.verifyEqual(unique(string(out.events.file_id)), "f1");
testCase.verifyFalse(any(string(out.events.event_uid) == "f1_E2"));

% map should only contain surviving event/cell ids.
testCase.verifyTrue(all(ismember(string(out.event_cell_map.event_uid), string(out.events.event_uid))));
testCase.verifyTrue(all(ismember(string(out.event_cell_map.cell_uid), string(out.cells.cell_uid))));
end

function testRunPipelineWithQC(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));
makeMiniRun(root);

loaded = loadClusterOutputs(root, 'attachLegacyMeta', false, 'attachEventMeta', false);
qc = createQCTables(loaded.files, loaded.cells, loaded.events);
qc.qc_files.keep(qc.qc_files.file_id == "f2") = false;

ep = struct();
ep.timeRange = [0 1; 1 2];
ep.epochNames = {'E1';'E2'};

out = runCanonicalAnalysisPipelineWithQC(root, ep, ...
    'qc_files', qc.qc_files, 'qc_cells', qc.qc_cells, 'qc_events', qc.qc_events, ...
    'metrics', {'area_px2'}, 'makePlots', false);

testCase.verifyEqual(unique(string(out.subset.files.file_id)), "f1");
testCase.verifyTrue(isfield(out.metrics, 'area_px2'));
testCase.verifyTrue(isfield(out, 'qcSummary'));
testCase.verifyTrue(istable(out.qcSummary.overall));
testCase.verifyTrue(any(out.qcSummary.overall.entity == "files"));
end

function testInvalidLegacyDataLimitTypeErrors(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));
makeMiniRun(root);

ep = struct();
ep.timeRange = [0 1; 1 2];
ep.epochNames = {'E1';'E2'};

didError = false;
try
    runCanonicalAnalysisPipelineWithQC(root, ep, ...
        'dataLimitRules', {'percentiles','area',[0 100]}, ...
        'makePlots', false);
catch ME
    didError = true;
    testCase.verifyEqual(ME.identifier, 'runCanonicalAnalysisPipelineWithQC:UnknownDataLimitType');
end
testCase.verifyTrue(didError);
end

function [files, cells, events, map] = miniCanonical()
files = table( ...
    string(["f1";"f2"]), ...
    string(["a1";"a2"]), ...
    string(["WT";"KO"]), ...
    'VariableNames', {'file_id','animal_id','genotype'});

cells = table( ...
    string(["f1_C1";"f1_C2";"f2_C1"]), ...
    string(["f1";"f1";"f2"]), ...
    'VariableNames', {'cell_uid','file_id'});

events = table( ...
    string(["f1_E1";"f1_E2";"f2_E1"]), ...
    string(["f1";"f1";"f2"]), ...
    [1;2;1], ...
    [10;70;10], ...
    [11;71;11], ...
    [10;20;30], ...
    [5;6;7], ...
    'VariableNames', {'event_uid','file_id','local_event_id','t0_sec','t1_sec','area_px2','duration_frames'});

map = table( ...
    string(["f1_E1";"f1_E2";"f2_E1"]), ...
    string(["f1_C1";"f1_C2";"f2_C1"]), ...
    true(3,1), ...
    'VariableNames', {'event_uid','cell_uid','is_primary_cell'});
end

function makeMiniRun(root)
runDir = fullfile(root, 'derived', 'ingest_runs', '20260224_111111');
mkdir(runDir);
[files, cells, events, event_cell_map] = miniCanonical();

BWA_AQ2 = struct();
BWA_AQ2.fts1 = struct();
save(fullfile(root, 'BWA_AQ2.mat'), 'BWA_AQ2');
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
