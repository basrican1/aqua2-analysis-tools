function tests = test_plotCorrelationStatsHeatmap
 tests = functiontests(localfunctions);
end

function testPlotsFromCorrelationOutput(testCase)
out = struct();
out.source = 'CorrelationStats';
out.method = 'spearman';
out.unitLevel = 'slice';
out.dataView = 'raw';
out.epochName = 'baseline';
out.metricIDs = {'area_um2','frequency','slice.frac_active_cells'};
out.R = [1 0.7 -0.3; 0.7 1 0.2; -0.3 0.2 1];
out.P = zeros(3);
out.N = 5*ones(3);
fig = plot_CorrelationStats_heatmap(out, 'cluster', false, 'showValues', true);
verifyTrue(testCase, isstruct(fig));
verifyTrue(testCase, isgraphics(fig.figure));
close(fig.figure);
end
