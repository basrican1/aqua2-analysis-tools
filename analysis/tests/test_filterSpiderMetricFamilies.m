function tests = test_filterSpiderMetricFamilies
tests = functiontests(localfunctions);
end

function testKeepShapeRecruitmentSynchrony(testCase)
pol = struct();
pol.labels = {'area_px2','footprint','mean_perimeter_um','frac_events_with_neighbor_recruitment','frac_primary_events_with_neighbor_recruitment','freq','coincidence_rate'};
pol.valuesA = reshape(1:14, 2, 7);
pol.valuesC = reshape(15:28, 2, 7);
out = filter_spider_metric_families(pol, 'metricFamilies', {'shape','recruitment','synchrony'});
verifyEqual(testCase, string(out.labels(:))', ["area_px2","footprint","mean_perimeter_um","frac_events_with_neighbor_recruitment","frac_primary_events_with_neighbor_recruitment","coincidence_rate"]);
verifyEqual(testCase, size(out.valuesA), [2 6]);
verifyEqual(testCase, size(out.valuesC), [2 6]);
end

function testKeepPropagation(testCase)
pol = struct();
pol.labels = {'area_px2','propdist','propagation_distance_um','mean_perimeter_um'};
pol.valuesA = reshape(1:8, 2, 4);
pol.valuesC = reshape(9:16, 2, 4);
out = filter_spider_metric_families(pol, 'metricFamilies', 'propagation');
verifyEqual(testCase, string(out.labels(:))', ["propagation_distance_um"]);
end

function testDropLegacyPropdistForAllFamilies(testCase)
pol = struct();
pol.labels = {'propdist','centroid','propagation_distance_um'};
pol.valuesA = reshape(1:6, 2, 3);
pol.valuesC = reshape(7:12, 2, 3);
out = filter_spider_metric_families(pol, 'metricFamilies', 'all');
verifyEqual(testCase, string(out.labels(:))', ["centroid","propagation_distance_um"]);
end
