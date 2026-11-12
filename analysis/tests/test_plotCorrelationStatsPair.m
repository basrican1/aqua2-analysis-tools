function tests = test_plotCorrelationStatsPair
 tests = functiontests(localfunctions);
end

function testPlotsSelectedPair(testCase)
out = struct();
out.source = 'CorrelationStats';
out.method = 'spearman';
out.metricIDs = {'area_um2','frequency'};
out.metricNames = {'area_um2','frequency'};
out.wideTable = table( ...
    string(["slice";"slice";"slice"]), ...
    string(["f1";"f2";"f3"]), ...
    [1;2;3], ...
    [2;4;6], ...
    'VariableNames', {'unit_level','entity_id','area_um2','frequency'});
fig = plot_CorrelationStats_pair(out, 'area_um2', 'frequency', 'showFit', true);
verifyTrue(testCase, isstruct(fig));
verifyTrue(testCase, isgraphics(fig.figure));
close(fig.figure);
end
