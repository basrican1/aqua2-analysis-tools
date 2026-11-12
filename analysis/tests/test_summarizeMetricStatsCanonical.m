function tests = test_summarizeMetricStatsCanonical
% Unit tests for summarizeMetricStatsCanonical
tests = functiontests(localfunctions);
end

function testUnpairedSummary(testCase)
T = table();
T.value = [1 2 3 6 7 8]';
T.group = categorical({'A','A','A','B','B','B'})';
T.animal_id = [1 2 3 1 2 3]';

out = summarizeMetricStatsCanonical(T, 'value', 'group', ...
    'statsConfig', struct('paired', false, 'parametric', true));

verifyTrue(testCase, isfield(out, 'stats'));
verifyEqual(testCase, out.stats.method, 'ttest2');
verifyTrue(testCase, isfinite(out.stats.pValue));
end

function testPairedSummaryWithPairVar(testCase)
T = table();
T.value = [1.0 2.0 3.0 1.4 2.5 3.6]';
T.group = categorical({'A','A','A','B','B','B'})';
T.animal_id = [1 2 3 1 2 3]';

out = summarizeMetricStatsCanonical(T, 'value', 'group', ...
    'pairVar', 'animal_id', ...
    'statsConfig', struct('paired', true, 'parametric', false));

verifyEqual(testCase, out.stats.method, 'signrank');
verifyTrue(testCase, isfinite(out.stats.pValue));
end

function testMissingVarErrors(testCase)
T = table();
T.value = [1 2]';
T.group = categorical({'A','B'})';

didErr = false;
try
    summarizeMetricStatsCanonical(T, 'not_here', 'group');
catch ME
    didErr = contains(ME.identifier, 'MissingValueVar');
end
verifyTrue(testCase, didErr);
end
