function tests = test_bwaMakeSliceLabels
tests = functiontests(localfunctions);
end

function testUsesCanonicalDateSliceWhenAvailable(testCase)
fileIDs = ["uuid_a"; "uuid_b"];
filesTbl = table( ...
    string(["uuid_a"; "uuid_b"]), ...
    string(["2025-02-05"; "2025-02-05"]), ...
    [1; 2], ...
    'VariableNames', {'file_id','exp_date','slice_num'});

labels = bwaMakeSliceLabels(fileIDs, filesTbl);

testCase.verifyEqual(labels{1}, '250205_S1');
testCase.verifyEqual(labels{2}, '250205_S2');
end

function testFallsBackToUniqueLabels(testCase)
fileIDs = ["file_one"; "file_one"];
labels = bwaMakeSliceLabels(fileIDs, table());

testCase.verifyEqual(numel(labels), 2);
testCase.verifyNotEqual(labels{1}, labels{2});
end
