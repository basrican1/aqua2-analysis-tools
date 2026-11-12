function tests = test_qcHelpers
% Tests for loadOrCreateQC / markQC / saveQC helper functions.

tests = functiontests(localfunctions);
end

function testLoadOrCreateAndSaveRoundTrip(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));
makeMiniRun(root);

[qc1, loaded, meta1] = loadOrCreateQC(root, 'reviewer', 'tester');
testCase.verifyEqual(meta1.source, 'created');
testCase.verifyEqual(height(qc1.qc_files), height(loaded.files));

qc1 = markQC(qc1, 'files', {'f2'}, false, 'motion_artifact', 'tester');
saveOut = saveQC(root, qc1, 'versioned', true);
testCase.verifyTrue(exist(saveOut.filesPath, 'file') == 2);
testCase.verifyTrue(exist(saveOut.versionDir, 'dir') == 7);
testCase.verifyTrue(exist(fullfile(root, 'qc_versions', 'latest_qc.txt'), 'file') == 2);

[qc2, ~, meta2] = loadOrCreateQC(root);
testCase.verifyEqual(meta2.source, 'loaded');
testCase.verifyFalse(qc2.qc_files.keep(qc2.qc_files.file_id=="f2"));
end

function testMarkQCBasic(testCase)
qc = struct();
qc.qc_files = table(string(["f1";"f2"]), true(2,1), strings(2,1), strings(2,1), strings(2,1), ...
    'VariableNames', {'file_id','keep','reason','reviewer','timestamp'});
qc.qc_cells = table(string(["c1"]), true, strings(1,1), strings(1,1), strings(1,1), ...
    'VariableNames', {'cell_uid','keep','reason','reviewer','timestamp'});
qc.qc_events = table(string(["e1"]), true, strings(1,1), strings(1,1), strings(1,1), ...
    'VariableNames', {'event_uid','keep','reason','reviewer','timestamp'});

qc = markQC(qc, 'files', {'f1'}, false, 'bad_segmentation', 'me');
testCase.verifyFalse(qc.qc_files.keep(qc.qc_files.file_id=="f1"));
testCase.verifyEqual(qc.qc_files.reason(qc.qc_files.file_id=="f1"), "bad_segmentation");
end

function makeMiniRun(root)
runDir = fullfile(root, 'derived', 'ingest_runs', '20260224_333333');
mkdir(runDir);

BWA_AQ2 = struct(); BWA_AQ2.fts1 = struct();
save(fullfile(root, 'BWA_AQ2.mat'), 'BWA_AQ2');

files = table(string(["f1";"f2"]), string(["a1";"a2"]), ...
    'VariableNames', {'file_id','animal_id'});
cells = table(string(["f1_C1";"f2_C1"]), string(["f1";"f2"]), ...
    'VariableNames', {'cell_uid','file_id'});
events = table(string(["f1_E1";"f2_E1"]), string(["f1";"f2"]), [1;1], [10;10], [12;12], ...
    'VariableNames', {'event_uid','file_id','local_event_id','t0_frame','t1_frame'});
event_cell_map = table(string(["f1_E1";"f2_E1"]), string(["f1_C1";"f2_C1"]), true(2,1), ...
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

