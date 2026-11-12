function tests = test_spiderSliceAggregation
tests = functiontests(localfunctions);
end

function testByEventUsesPerFileEventEntities(testCase)
cleanup = suppressFigures();
pol = makeMiniPolars();
out = plotBWA_Aqua2PolarsBySlice(pol, ...
    'showLegend', false, ...
    'exportCompanionTable', false, ...
    'spiderInputScale', 'log2', ...
    'branchMode', 'E', ...
    'comparisonLevel', 'by_event');
T = out.companionTable;
TA = T(T.group == "E", :);
verifyEqual(testCase, double(TA.value_log2)', [2 3], 'AbsTol', 1e-9);
verifyEqual(testCase, string(TA.aggregate_by)', ["events" "events"]);
verifyEqual(testCase, logical(TA.entity_backed)', [true true]);
verifyEqual(testCase, string(out.statsInfo.comparison_level), "by_event");
verifyEqual(testCase, string(out.statsSummary.group), "E");
verifyEqual(testCase, strlength(string(out.statsSummary.metric_family)) >= 0, true);
verifyEqual(testCase, double(out.statsSummary.f1_log2), 2, 'AbsTol', 1e-9);
verifyEqual(testCase, double(out.statsSummary.f2_log2), 3, 'AbsTol', 1e-9);
verifyEqual(testCase, double(out.statsSummary.f1_n), 2, 'AbsTol', 1e-9);
verifyEqual(testCase, double(out.statsSummary.f2_n), 2, 'AbsTol', 1e-9);
closeIfOpen(out);
clear cleanup
end

function testByBranchUsesPerFileEventAndCellEntities(testCase)
cleanup = suppressFigures();
pol = makeMiniPolars();
out = plotBWA_Aqua2PolarsBySlice(pol, ...
    'showLegend', false, ...
    'exportCompanionTable', false, ...
    'spiderInputScale', 'log2', ...
    'branchMode', 'both', ...
    'comparisonLevel', 'by_branch');
T = out.companionTable;
TA = T(T.group == "E", :);
TC = T(T.group == "C", :);
verifyEqual(testCase, double(TA.value_log2)', [2 3], 'AbsTol', 1e-9);
verifyEqual(testCase, double(TC.value_log2)', [20 30], 'AbsTol', 1e-9);
verifyEqual(testCase, string(TA.aggregate_by)', ["events" "events"]);
verifyEqual(testCase, string(TC.aggregate_by)', ["cells" "cells"]);
verifyEqual(testCase, logical(TA.entity_backed)', [true true]);
verifyEqual(testCase, logical(TC.entity_backed)', [true true]);
verifyEqual(testCase, string(out.statsInfo.comparison_level), "by_branch");
verifyEqual(testCase, sort(string(out.statsSummary.group))', ["C" "E"]);
closeIfOpen(out);
clear cleanup
end


function testSliceStatsSummarySortedByMetricFamily(testCase)
cleanup = suppressFigures();
pol = makeMiniPolars();
pol.labels = {'area_px2','dffauc','duration_s'};
pol.valuesA = [4 1 2; 8 3 6];
pol.valuesC = [40 10 20; 80 30 60];
pol.entityValuesA = {[4;8], [1;3], [2;6]};
pol.entityValuesC = {[40;80], [10;30], [20;60]};
pol.entityValuesAByFile = {{[4],[8]}, {[1],[3]}, {[2],[6]}};
pol.entityValuesCByFile = {{[40],[80]}, {[10],[30]}, {[20],[60]}};
out = plotBWA_Aqua2PolarsBySlice(pol, ...
    'showLegend', false, ...
    'exportCompanionTable', false, ...
    'spiderInputScale', 'log2', ...
    'branchMode', 'E', ...
    'comparisonLevel', 'by_event');
verifyEqual(testCase, string(out.statsSummary.metric)', ["dffauc" "duration_s" "area_px2"]);
verifyEqual(testCase, string(out.statsSummary.metric_family)', ["intensity" "kinetics" "shape"]);
closeIfOpen(out);
clear cleanup
end

function pol = makeMiniPolars()
pol = struct();
pol.labels = {'m1'};
pol.valuesA = [10; 100];
pol.valuesC = [40; 400];
pol.entityValuesA = {[2; 3]};
pol.entityValuesC = {[20; 30]};
pol.entityValuesAByFile = {{[1; 3], [2; 4]}};
pol.entityValuesCByFile = {{[10; 30], [20; 40]}};
pol.uniqueFiles = ["f1"; "f2"];
pol.fileLabels = {'f1','f2'};
pol.subplotGrid = [1 2];
end

function cleanup = suppressFigures()
old = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
cleanup = onCleanup(@() set(groot, 'DefaultFigureVisible', old));
end

function closeIfOpen(out)
if isstruct(out) && isfield(out, 'figure') && ~isempty(out.figure) && isgraphics(out.figure)
    close(out.figure);
end
end
