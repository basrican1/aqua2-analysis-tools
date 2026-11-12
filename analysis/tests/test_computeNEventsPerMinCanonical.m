function tests = test_computeNEventsPerMinCanonical
% Unit tests for computeNEventsPerMinCanonical

tests = functiontests(localfunctions);
end

function testBalancedRatesByEpoch(testCase)
[events, map, files, epochSpec] = makeBalancedFixture();

out = computeNEventsPerMinCanonical(events, map, files, epochSpec, ...
    'groupBy', {'file_id'}, 'dropZeroRates', true, 'primaryOnly', true);

S = out.summaryByGroupEpoch;
[~, ord] = sort(S.epoch_idx);
S = S(ord,:);

testCase.verifyEqual(height(S), 2);
testCase.verifyEqual(S.mean_rate_events_per_min, [1;1], 'AbsTol', 1e-12);
testCase.verifyEqual(S.median_rate_events_per_min, [1;1], 'AbsTol', 1e-12);
testCase.verifyEqual(S.std_rate_events_per_min, [0;0], 'AbsTol', 1e-12);
testCase.verifyEqual(S.n_cells, [2;2]);
end

function testDropZeroRatesFalseFillsMissingCombos(testCase)
[events, map, files, epochSpec] = makeSparseFixture();

out = computeNEventsPerMinCanonical(events, map, files, epochSpec, ...
    'groupBy', {'file_id'}, 'dropZeroRates', false, 'primaryOnly', true);

S = out.summaryByGroupEpoch;
[~, ord] = sort(S.epoch_idx);
S = S(ord,:);

testCase.verifyEqual(height(S), 2);
testCase.verifyEqual(S.mean_rate_events_per_min, [0.5;0.5], 'AbsTol', 1e-12);
testCase.verifyEqual(S.n_cells, [2;2]);
end

function [events, map, files, epochSpec] = makeBalancedFixture()
% Two cells, each has one event in each 1-min epoch => rate = 1 for all cells/epochs.
events = table( ...
    string(["e1";"e2";"e3";"e4"]), ...
    string(["f1";"f1";"f1";"f1"]), ...
    [10;20;70;80], ... % seconds
    [10;20;30;40], ... % area_px2
    'VariableNames', {'event_uid','file_id','t0_sec','area_px2'});

map = table( ...
    string(["e1";"e2";"e3";"e4"]), ...
    string(["f1_C1";"f1_C2";"f1_C1";"f1_C2"]), ...
    true(4,1), ...
    'VariableNames', {'event_uid','cell_uid','is_primary_cell'});

files = table(string("f1"), string("a1"), 'VariableNames', {'file_id','animal_id'});

epochSpec = struct();
epochSpec.timeRange = [0 1; 1 2];
epochSpec.epochNames = {'E1';'E2'};
end

function [events, map, files, epochSpec] = makeSparseFixture()
% Cell1 only epoch1, Cell2 only epoch2. With filled zeros, mean rate is 0.5 each epoch.
events = table( ...
    string(["e1";"e2"]), ...
    string(["f1";"f1"]), ...
    [10;70], ...
    [10;20], ...
    'VariableNames', {'event_uid','file_id','t0_sec','area_px2'});

map = table( ...
    string(["e1";"e2"]), ...
    string(["f1_C1";"f1_C2"]), ...
    true(2,1), ...
    'VariableNames', {'event_uid','cell_uid','is_primary_cell'});

files = table(string("f1"), string("a1"), 'VariableNames', {'file_id','animal_id'});

epochSpec = struct();
epochSpec.timeRange = [0 1; 1 2];
epochSpec.epochNames = {'E1';'E2'};
end
