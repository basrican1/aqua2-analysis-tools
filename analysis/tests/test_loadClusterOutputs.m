function tests = test_loadClusterOutputs
% Unit tests for loadClusterOutputs

tests = functiontests(localfunctions);
end

function testLoadsFromRootAndRunDir(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));

runDir = makeMiniRun(root);

out1 = loadClusterOutputs(root, 'attachLegacyMeta', false, 'attachEventMeta', false);
out2 = loadClusterOutputs(runDir, 'attachLegacyMeta', false, 'attachEventMeta', false);

testCase.verifyEqual(out1.runDir, runDir);
testCase.verifyEqual(out2.runDir, runDir);
testCase.verifyEqual(height(out1.files), 1);
testCase.verifyEqual(height(out2.events), 1);
end

function testMissingRunDirErrors(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));

didError = false;
try
    loadClusterOutputs(root, 'attachLegacyMeta', false, 'attachEventMeta', false);
catch ME
    didError = true;
    testCase.verifyTrue(contains(ME.message, 'Could not infer canonical run dir'));
end
testCase.verifyTrue(didError);
end

function testLoadsCanonicalWithoutLegacyPayload(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));

runDir = makeMiniRun(root, false);
out = loadClusterOutputs(runDir, 'attachLegacyMeta', false, 'attachEventMeta', false);

testCase.verifyEqual(out.runDir, runDir);
testCase.verifyFalse(logical(out.hasLegacyPayload));
testCase.verifyTrue(isstruct(out.BWA_AQ2));
testCase.verifyTrue(isempty(fieldnames(out.BWA_AQ2)));
testCase.verifyEqual(height(out.files), 1);
testCase.verifyEqual(height(out.events), 1);
end

function runDir = makeMiniRun(root, includeLegacy)
if nargin < 2
    includeLegacy = true;
end
runDir = fullfile(root, 'derived', 'ingest_runs', '20260224_000001');
mkdir(runDir);

if includeLegacy
    BWA_AQ2 = struct();
    BWA_AQ2.fts1 = struct();
    save(fullfile(root, 'BWA_AQ2.mat'), 'BWA_AQ2');
end

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
