function tests = test_analysisPlanWriter
tests = functiontests(localfunctions);
end

function testWriteAndLoadPreset(testCase)
root = tempname;
mkdir(root);

outPath = writeAnalysisPlanPreset(root, ...
    'presetName', 'my_preset', ...
    'setDefault', true, ...
    'epochPreset', 'drug_application', ...
    'dataLimitType', 'duration_frames', ...
    'dataLimitPercentiles', [10 95], ...
    'metrics', {'area_px2','duration_frames'});

testCase.verifyTrue(exist(outPath, 'file') == 2);
[plan, meta] = loadAnalysisPlan(root, 'presetName', 'my_preset');
testCase.verifyEqual(string(meta.presetName), "my_preset");
testCase.verifyEqual(string(plan.epochPreset), "drug_application");
testCase.verifyEqual(string(plan.dataLimitType), "duration_frames");
testCase.verifyEqual(double(plan.dataLimitPercentiles(:))', [10 95]);
testCase.verifyEqual(cellstr(string(plan.metrics(:)')), {'area_px2','duration_frames'});
end
