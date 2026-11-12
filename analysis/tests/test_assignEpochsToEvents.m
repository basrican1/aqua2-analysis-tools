function tests = test_assignEpochsToEvents
% Unit tests for assignEpochsToEvents

tests = functiontests(localfunctions);
end

function testFrameAssignment(testCase)
events = table([5;15;25], 'VariableNames', {'t0_frame'});
epochSpec = struct();
epochSpec.frameRange = [1 10; 11 20];
epochSpec.epochNames = {'A'; 'B'};

out = assignEpochsToEvents(events, epochSpec, 'method', 'frame');

testCase.verifyEqual(out.events.epoch_idx, [1;2;0]);
testCase.verifyEqual(string(out.events.epoch_name), string(["A";"B";""]));
testCase.verifyEqual(height(out.epochTable), 2);
end

function testTimeAssignment(testCase)
events = table([30;90], 'VariableNames', {'t0_sec'});
epochSpec = struct();
epochSpec.timeRange = [0 1; 1 2];
epochSpec.epochNames = {'Base'; 'Drug'};

out = assignEpochsToEvents(events, epochSpec, 'method', 'time');

testCase.verifyEqual(out.events.epoch_idx, [1;2]);
testCase.verifyEqual(string(out.events.epoch_name), string(["Base";"Drug"]));
end

function testMissingRangesErrors(testCase)
events = table([1;2], 'VariableNames', {'t0_frame'});
epochSpec = struct();

didError = false;
try
    assignEpochsToEvents(events, epochSpec);
catch ME
    didError = true;
    testCase.verifyTrue(contains(ME.message, 'epochSpec must include frameRange or timeRange'));
end
testCase.verifyTrue(didError);
end
