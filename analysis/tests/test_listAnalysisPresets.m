function tests = test_listAnalysisPresets
% Tests for listing analysis presets.

tests = functiontests(localfunctions);
end

function testListFromTemplate(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));

% Minimal root marker.
BWA_AQ2 = struct(); %#ok<NASGU>
save(fullfile(root, 'BWA_AQ2.mat'), 'BWA_AQ2');
writeAnalysisPlanTemplate(root);

T = listAnalysisPresets(root, 'printTable', false);
testCase.verifyTrue(istable(T));
testCase.verifyGreaterThan(height(T), 0);
testCase.verifyTrue(all(ismember({'preset','is_default','epoch_preset'}, T.Properties.VariableNames)));
end

function safeRmdir(pathIn)
if exist(pathIn, 'dir')
    rmdir(pathIn, 's');
end
end

