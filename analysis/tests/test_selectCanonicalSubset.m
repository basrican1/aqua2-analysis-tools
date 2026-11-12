function tests = test_selectCanonicalSubset
% Unit tests for selectCanonicalSubset

tests = functiontests(localfunctions);
end

function testFilterByGenotype(testCase)
[files, cells, events, map] = makeCanonicalFixtures();

S = selectCanonicalSubset(files, cells, events, map, 'genotypes', {'WT'});

testCase.verifyEqual(height(S.files), 1);
testCase.verifyEqual(string(S.files.file_id), "f1");
testCase.verifyEqual(unique(string(S.events.file_id)), "f1");
testCase.verifyEqual(unique(string(S.cells.file_id)), "f1");

% Map rows should only reference kept event/cell ids
keptEvents = string(S.events.event_uid);
keptCells = string(S.cells.cell_uid);
testCase.verifyTrue(all(ismember(string(S.event_cell_map.event_uid), keptEvents)));
testCase.verifyTrue(all(ismember(string(S.event_cell_map.cell_uid), keptCells)));
end

function testFilterByQCMask(testCase)
[files, cells, events, map] = makeCanonicalFixtures();

S = selectCanonicalSubset(files, cells, events, map, 'qc_mask_files', [false; true]);

testCase.verifyEqual(height(S.files), 1);
testCase.verifyEqual(string(S.files.file_id), "f2");
testCase.verifyEqual(unique(string(S.events.file_id)), "f2");
end

function [files, cells, events, map] = makeCanonicalFixtures()
files = table( ...
    string(["f1"; "f2"]), ...
    string(["a1"; "a2"]), ...
    string(["WT"; "KO"]), ...
    string(["veh"; "drug"]), ...
    'VariableNames', {'file_id','animal_id','genotype','treatment'});

cells = table( ...
    string(["f1_C1"; "f2_C1"]), ...
    string(["f1"; "f2"]), ...
    'VariableNames', {'cell_uid','file_id'});

events = table( ...
    string(["f1_E1"; "f2_E1"]), ...
    string(["f1"; "f2"]), ...
    [10; 20], ...
    'VariableNames', {'event_uid','file_id','t0_frame'});

map = table( ...
    string(["f1_E1"; "f2_E1"]), ...
    string(["f1_C1"; "f2_C1"]), ...
    'VariableNames', {'event_uid','cell_uid'});
end
