function tests = test_shared_manifest_preflight
tests = functiontests(localfunctions);
end

function testValidSharedFixtureRecordsProvenance(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
sharedRoot = fullfile(fileparts(fileparts(repoRoot)), 'SharedAnalysisMetadata');
fixture = fullfile(sharedRoot, 'tests', 'fixtures', 'v1', 'mixed_factor_valid');
testCase.assumeTrue(isfolder(fixture), 'SharedAnalysisMetadata fixture is unavailable locally.');

addpath(genpath(repoRoot));
preflight = preflightSharedAnalysisProject(fixture);
verifyTrue(testCase, preflight.validated);
verifyEqual(testCase, string(preflight.schema_version), "1.0.0");
verifyEqual(testCase, strlength(string(preflight.shared_metadata_revision)), 40);
verifyTrue(testCase, islogical(preflight.shared_metadata_dirty));
end

function testOptionalAbsentManifestIsNoOp(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
preflight = preflightSharedAnalysisProject('', 'Required', false);
verifyFalse(testCase, preflight.requested);
verifyFalse(testCase, preflight.validated);
end

function testRequiredAbsentManifestFails(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
verifyError(testCase, @() preflightSharedAnalysisProject('', 'Required', true), ...
    'post2:MissingProjectManifest');
end
