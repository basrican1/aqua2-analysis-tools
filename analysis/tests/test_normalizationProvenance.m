function tests = test_normalizationProvenance
tests = functiontests(localfunctions);
end

function testRawViewReportsRaw(testCase)
R = makeMiniResult();
prov = bwaResolveNormalizationProvenance(R, 'raw', ["area_um2"; "cell.mean_perimeter_um"], ["canonical"; "cell"]);
verifyEqual(testCase, string(prov.requestedView), "raw");
verifyEqual(testCase, string(prov.compactLabel), "norm=raw");
verifyEqual(testCase, string(prov.metricTable.effective_transform)', ["raw" "raw"]);
end

function testNormViewReportsDefaultTransform(testCase)
R = makeMiniResult();
prov = bwaResolveNormalizationProvenance(R, 'norm', ["area_um2"; "cell.mean_perimeter_um"], ["canonical"; "cell"]);
verifyEqual(testCase, string(prov.normalizationMode), "by_cell_baseline_mean");
verifyEqual(testCase, string(prov.defaultTransform), "delta");
verifyEqual(testCase, string(prov.compactLabel), "norm=default [delta]");
verifyEqual(testCase, string(prov.metricTable.effective_transform)', ["delta" "delta"]);
verifyEqual(testCase, string(prov.baselineEpochName), "baseline");
end

function testRecommendedViewResolvesMixedTransforms(testCase)
R = makeMiniResult();
prov = bwaResolveNormalizationProvenance(R, 'recommended', ["area_um2"; "cell.mean_perimeter_um"], ["canonical"; "cell"]);
verifyEqual(testCase, string(prov.compactLabel), "norm=mixed [log2_fold, delta]");
verifyEqual(testCase, string(prov.metricTable.effective_transform)', ["log2_fold" "delta"]);
end

function R = makeMiniResult()
epochTable = table([1;2], ["baseline";"drug"], 'VariableNames', {'epoch_idx','epoch_name'});
R = struct();
R.normalization = struct('mode', 'by_cell_baseline_mean', 'transform', 'delta', 'baselineEpochIdx', 1);
R.normalized = struct('default_transform', 'delta');
R.canonical = struct('frequency', struct('epochTable', epochTable));
end
