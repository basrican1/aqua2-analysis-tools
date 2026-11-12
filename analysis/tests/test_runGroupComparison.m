function tests = test_runGroupComparison
% Unit tests for runGroupComparison
tests = functiontests(localfunctions);
end

function testTwoGroupUnpairedAutoParametric(testCase)
values = [1 2 3 6 7 8]';
groups = categorical({'A','A','A','B','B','B'})';
out = runGroupComparison(values, groups, 'method', 'auto', 'paired', false, 'parametric', true);
verifyEqual(testCase, out.method, 'ttest2');
verifyTrue(testCase, isfinite(out.pValue));
verifySize(testCase, out.pairwiseC, [1 6]);
end

function testTwoGroupPairedAutoSignrank(testCase)
values = [1.0 1.2 0.8 2.0 2.3 1.8]';
groups = categorical({'A','A','A','B','B','B'})';
pairID = [1 2 3 1 2 3]';
out = runGroupComparison(values, groups, 'method', 'auto', 'paired', true, 'parametric', false, 'pairID', pairID);
verifyEqual(testCase, out.method, 'signrank');
verifyTrue(testCase, isfinite(out.pValue));
end

function testMultiGroupUnpairedAutoKruskal(testCase)
values = [1 2 3 6 7 8 10 11 12]';
groups = categorical({'A','A','A','B','B','B','C','C','C'})';
out = runGroupComparison(values, groups, 'method', 'auto', 'paired', false, 'parametric', false);
verifyEqual(testCase, out.method, 'kruskal');
verifyTrue(testCase, isfinite(out.pValue));
verifyGreaterThan(testCase, size(out.pairwiseC,1), 0);
end

function testMultiGroupPairedAutoFriedman(testCase)
values = [1 3 2 2 4 2.5 3 5 4]';
groups = categorical({'A','B','C','A','B','C','A','B','C'})';
pairID = [1 1 1 2 2 2 3 3 3]';
out = runGroupComparison(values, groups, 'method', 'auto', 'paired', true, 'parametric', true, 'pairID', pairID);
verifyEqual(testCase, out.method, 'friedman');
verifyTrue(testCase, isfinite(out.pValue));
end

function testPairedRequiresPairID(testCase)
values = [1 2 3 4]';
groups = categorical({'A','A','B','B'})';
didErr = false;
try
    runGroupComparison(values, groups, 'method', 'ttest_paired');
catch ME
    didErr = contains(ME.identifier, 'PairIDRequired');
end
verifyTrue(testCase, didErr);
end
