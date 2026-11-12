function tests = test_computeMetricBySliceCanonical
% Unit tests for computeMetricBySliceCanonical

tests = functiontests(localfunctions);
end

function testComputesExpectedMeans(testCase)
[events, files, ep] = makeFixture();

out = computeMetricBySliceCanonical(events, files, 'area_px2', ep, 'groupBy', {'file_id'});
S = out.summaryByGroupEpoch;
[~, ord] = sort(S.epoch_idx);
S = S(ord,:);

testCase.verifyEqual(height(S), 2);
testCase.verifyEqual(S.mean_area_px2, [15; 35], 'AbsTol', 1e-12);
testCase.verifyEqual(S.median_area_px2, [15; 35], 'AbsTol', 1e-12);
end

function testMetricFilterDropsNaN(testCase)
[events, files, ep] = makeFixture();
events.area_px2(2) = NaN;

out = computeMetricBySliceCanonical(events, files, 'area_px2', ep, ...
    'groupBy', {'file_id'}, 'metricFilter', @(v) isfinite(v));

S = out.summaryByGroupEpoch;
[~, ord] = sort(S.epoch_idx);
S = S(ord,:);

testCase.verifyEqual(S.n_events, [1;2]);
testCase.verifyEqual(S.mean_area_px2, [10;35], 'AbsTol', 1e-12);
end

function testMissingMetricErrors(testCase)
[events, files, ep] = makeFixture();

didError = false;
try
    computeMetricBySliceCanonical(events, files, 'not_a_metric', ep);
catch ME
    didError = true;
    testCase.verifyTrue(contains(ME.message, 'Metric column not found in events'));
end
testCase.verifyTrue(didError);
end

function [events, files, ep] = makeFixture()
events = table( ...
    string(["e1";"e2";"e3";"e4"]), ...
    string(["f1";"f1";"f1";"f1"]), ...
    [10;20;70;80], ...
    [10;20;30;40], ...
    'VariableNames', {'event_uid','file_id','t0_sec','area_px2'});

files = table(string("f1"), string("a1"), string("WT"), ...
    'VariableNames', {'file_id','animal_id','genotype'});

ep = struct();
ep.timeRange = [0 1; 1 2];
ep.epochNames = {'E1';'E2'};
end
