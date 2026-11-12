function tests = test_spiderMeanAggregation
tests = functiontests(localfunctions);
end

function testDefaultSliceAggregationUsesSliceMeans(testCase)
cleanup = suppressFigures();
R = makeMiniSpiderResult();
out = plot_local_analysis_spider_mean(R, ...
    'spiderInputScale', 'log2', ...
    'branchMode', 'both', ...
    'metricFamilies', {'intensity'}, ...
    'showStatsFigure', false);
verifyEqual(testCase, out.meanValuesA, 2, 'AbsTol', 1e-9);
verifyEqual(testCase, out.meanValuesC, 5, 'AbsTol', 1e-9);
verifyEqual(testCase, out.statsSummary.n_E, 2);
verifyEqual(testCase, out.statsSummary.n_C, 2);
closeIfOpen(out);
clear cleanup
end

function testEntityAggregationUsesPooledEventsAndCells(testCase)
cleanup = suppressFigures();
R = makeMiniSpiderResult();
out = plot_local_analysis_spider_mean(R, ...
    'spiderInputScale', 'log2', ...
    'branchMode', 'both', ...
    'metricFamilies', {'intensity'}, ...
    'comparisonLevel', 'by_branch', ...
    'showStatsFigure', false);
verifyEqual(testCase, out.meanValuesA, mean([1 1 1 3]), 'AbsTol', 1e-9);
verifyEqual(testCase, out.meanValuesC, mean([10 10 0]), 'AbsTol', 1e-9);
verifyEqual(testCase, out.statsSummary.n_E, 4);
verifyEqual(testCase, out.statsSummary.n_C, 3);
verifyEqual(testCase, string(out.statsInfo.comparison_level), "by_branch");
closeIfOpen(out);
clear cleanup
end

function testEventModeRelabelsCellDerivedEventAnalogs(testCase)
cleanup = suppressFigures();
R = makeMiniSpiderResult();
R.derived.polars.labels = {'cell_mean_perimeter_um', 'cell_mean_circularity'};
R.derived.polars.valuesA = [1 2; 3 4];
R.derived.polars.valuesC = [10 20; 0 40];
R.derived.polars.entityValuesA = {[1;1;3], [2;2;4]};
R.derived.polars.entityValuesC = {[10;0], [20;40]};
out = plot_local_analysis_spider_mean(R, ...
    'spiderInputScale', 'log2', ...
    'branchMode', 'E', ...
    'comparisonLevel', 'by_event', ...
    'metricFamilies', {'shape'}, ...
    'showStatsFigure', false);
verifyEqual(testCase, string(out.statsSummary.metric)', ["circularity" "perimeter_um"]);
closeIfOpen(out);
clear cleanup
end

function testByBranchUsesGenericLabelsAndPooledCells(testCase)
cleanup = suppressFigures();
R = makeMiniSpiderResult();
R.derived.polars.labels = {'cell_mean_perimeter_um', 'cell_mean_circularity'};
R.derived.polars.valuesA = [1 2; 3 4];
R.derived.polars.valuesC = [10 20; 0 40];
R.derived.polars.entityValuesA = {[1;1;3], [2;2;4]};
R.derived.polars.entityValuesC = {[10;20;30], [20;40;60]};
out = plot_local_analysis_spider_mean(R, ...
    'spiderInputScale', 'log2', ...
    'branchMode', 'both', ...
    'comparisonLevel', 'by_branch', ...
    'metricFamilies', {'shape'}, ...
    'showStatsFigure', false);
verifyEqual(testCase, sort(string(out.statsSummary.metric))', ["circularity" "perimeter_um"]);
verifyEqual(testCase, sort(double(out.statsSummary.n_E))', [3 3]);
verifyEqual(testCase, sort(double(out.statsSummary.n_C))', [3 3]);
closeIfOpen(out);
clear cleanup
end

function testStatsSummarySortedByMetricFamily(testCase)
cleanup = suppressFigures();
R = makeMiniSpiderResult();
R.derived.polars.labels = {'area_px2', 'dffauc', 'duration_s'};
R.derived.polars.valuesA = [4 1 2; 8 3 6];
R.derived.polars.valuesC = [40 10 20; 80 30 60];
R.derived.polars.entityValuesA = {[4;8], [1;3], [2;6]};
R.derived.polars.entityValuesC = {[40;80], [10;30], [20;60]};
out = plot_local_analysis_spider_mean(R, ...
    'spiderInputScale', 'log2', ...
    'branchMode', 'E', ...
    'comparisonLevel', 'by_event', ...
    'metricFamilies', 'all', ...
    'showStatsFigure', false);
verifyEqual(testCase, string(out.statsSummary.metric)', ["dffauc" "duration_s" "area_px2"]);
verifyEqual(testCase, string(out.statsSummary.metric_family)', ["intensity" "kinetics" "shape"]);
closeIfOpen(out);
clear cleanup
end

function testEntityPointOverlayForEventMode(testCase)
cleanup = suppressFigures();
R = makeMiniSpiderResult();
R.derived.polars.labels = {'dffauc'};
R.derived.polars.valuesA = [1; 3];
R.derived.polars.valuesC = [10; 0];
R.derived.polars.entityValuesA = {[1; 1; 1; 3]};
out = plot_local_analysis_spider_mean(R, ...
    'spiderInputScale', 'log2', ...
    'branchMode', 'E', ...
    'comparisonLevel', 'by_event', ...
    'metricFamilies', {'intensity'}, ...
    'showEntityPoints', true, ...
    'showStatsFigure', false);
verifyTrue(testCase, isfield(out, 'entityPoints'));
verifyGreaterThan(testCase, numel(out.entityPoints), 0);
closeIfOpen(out);
clear cleanup
end

function cleanup = suppressFigures()
old = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
cleanup = onCleanup(@() set(groot, 'DefaultFigureVisible', old));
end

function R = makeMiniSpiderResult()
pol = struct();
pol.labels = {'dffauc'};
pol.valuesA = [1; 3];
pol.valuesC = [10; 0];
pol.entityValuesA = {[1; 1; 1; 3]};
pol.entityValuesC = {[10; 10; 0]};
pol.uniqueFiles = ["f1"; "f2"];
pol.subplotGrid = [1 2];
pol.fileLabels = {'f1','f2'};

R = struct();
R.canonical = struct();
R.derived = struct('polars', pol);
R.normalization = struct('mode', 'by_cell_baseline_mean', 'transform', 'delta');
end

function closeIfOpen(out)
if isstruct(out) && isfield(out, 'figure') && ~isempty(out.figure) && isgraphics(out.figure)
    close(out.figure);
end
if isstruct(out) && isfield(out, 'statsFigure') && ~isempty(out.statsFigure) && isgraphics(out.statsFigure)
    close(out.statsFigure);
end
end
