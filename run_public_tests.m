function results = run_public_tests()
% Selected existing tests with self-created artificial inputs, not private runs.
root = setup_aqua2_tools();
names = {'test_assignEpochsToEvents.m','test_computeNEventsPerMinCanonical.m', ...
    'test_computeMetricBySliceCanonical.m','test_qcHelpers.m','test_qcWorkflow.m', ...
    'test_runGroupComparison.m','test_selectCanonicalSubset.m', ...
    'test_CorrelationStats.m','test_metricVocabulary.m', ...
    'test_plotCorrelationStatsHeatmap.m','test_plotCorrelationStatsPair.m'};
if exist('spider_plot', 'file') == 2
    names = [names {'test_spiderMeanAggregation.m','test_spiderSliceAggregation.m'}];
else
    fprintf('Optional spider rendering suites not run: install spider_plot separately.\n');
end
paths = cellfun(@(n) fullfile(root,'analysis','tests',n), names, 'UniformOutput', false);
results = runtests(paths);
assertSuccess(results);
run_synthetic_workflow();
end
