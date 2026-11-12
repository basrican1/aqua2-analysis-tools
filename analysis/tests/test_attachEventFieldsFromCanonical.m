function tests = test_attachEventFieldsFromCanonical
% Unit tests for attachEventFieldsFromCanonical

tests = functiontests(localfunctions);
end

function testReconstructsAndRespectsFileOrder(testCase)
BWA_AQ2 = struct();
BWA_AQ2.fts1 = struct();

% Events intentionally not ordered by file_id to verify filesTbl ordering is used.
events = table( ...
    string(["f2_E1";"f1_E1";"f1_E2"]), ...
    string(["f2";"f1";"f1"]), ...
    [1;1;2], ...
    [201;101;102], ...
    [211;111;112], ...
    'VariableNames', {'event_uid','file_id','local_event_id','t0_frame','t1_frame'});

files = table(string(["f1";"f2"]), 'VariableNames', {'file_id'});

out = attachEventFieldsFromCanonical(BWA_AQ2, events, files);

testCase.verifyEqual(out.fts1.loc.t0, [101 102 201]);
testCase.verifyEqual(out.fts1.loc.t1, [111 112 211]);
testCase.verifyEqual(string(out.fts1.fileID), string(["f1" "f1" "f2"]));
end

function testNoOverwriteWhenAlreadyPresent(testCase)
BWA_AQ2.fts1.loc.t0 = [1 2 3];
BWA_AQ2.fts1.loc.t1 = [4 5 6];
BWA_AQ2.fts1.fileID = {'x','x','x'};

events = table(string("f1_E1"), string("f1"), 1, 10, 12, 'VariableNames', {'event_uid','file_id','local_event_id','t0_frame','t1_frame'});
files = table(string("f1"), 'VariableNames', {'file_id'});

out = attachEventFieldsFromCanonical(BWA_AQ2, events, files);

testCase.verifyEqual(out.fts1.loc.t0, [1 2 3]);
testCase.verifyEqual(out.fts1.loc.t1, [4 5 6]);
testCase.verifyEqual(string(out.fts1.fileID), string(["x" "x" "x"]));
end
