function tests = test_plotBAviolinSubplotAnova
% Unit tests for plotBAviolinSubplotAnova
tests = functiontests(localfunctions);
end

function testUnpairedAuto(testCase)
values = [1 2 3 7 8 9]';
groups = categorical({'A','A','A','B','B','B'})';

h = figure('Visible', 'off');
cleanup = onCleanup(@() closeIfValid(h)); %#ok<NASGU>
ax1 = axes('Parent', h, 'Position', [0.08 0.12 0.40 0.78]);
ax2 = axes('Parent', h, 'Position', [0.52 0.12 0.40 0.78]);

out = plotBAviolinSubplotAnova(values, groups, ...
    'figure', h, 'axesPlot', ax1, 'axesTable', ax2, ...
    'statsConfig', struct('method', 'auto', 'paired', false, 'parametric', true));

verifyEqual(testCase, out.stats.method, 'ttest2');
verifyTrue(testCase, isfinite(out.pvalue));
verifyTrue(testCase, istable(out.stats.summaryTable));
end

function testPairedAuto(testCase)
values = [1 2 3 1.5 2.5 3.5]';
groups = categorical({'A','A','A','B','B','B'})';
pairID = [1 2 3 1 2 3]';

h = figure('Visible', 'off');
cleanup = onCleanup(@() closeIfValid(h)); %#ok<NASGU>
ax1 = axes('Parent', h, 'Position', [0.08 0.12 0.40 0.78]);
ax2 = axes('Parent', h, 'Position', [0.52 0.12 0.40 0.78]);

out = plotBAviolinSubplotAnova(values, groups, ...
    'figure', h, 'axesPlot', ax1, 'axesTable', ax2, ...
    'statsConfig', struct('method', 'auto', 'paired', true, 'parametric', false), ...
    'pairID', pairID);

verifyEqual(testCase, out.stats.method, 'signrank');
verifyTrue(testCase, isfinite(out.pvalue));
end

function closeIfValid(h)
if ~isempty(h) && isgraphics(h)
    close(h);
end
end
