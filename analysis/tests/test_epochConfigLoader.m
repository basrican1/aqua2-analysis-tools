function tests = test_epochConfigLoader
% Tests for config-driven epoch loading.

tests = functiontests(localfunctions);
end

function testDefaultFallbackAllEpoch(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));

[ep, meta] = loadEpochSpecFromConfig(root, 'sz_t', 100, 'frameRate_s', 1);
testCase.verifyTrue(isstruct(ep));
testCase.verifyTrue(meta.usedDefault);
testCase.verifyEqual(size(ep.frameRange,1), 1);
testCase.verifyEqual(ep.epochNames{1}, 'All');
end

function testJsonFramesConfig(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));

cfg = struct();
cfg.epochs = struct();
cfg.epochs.mode = 'frames';
cfg.epochs.frameRange = [10 20; 21 30];
cfg.epochs.epochNames = {'B'; 'D'};

fid = fopen(fullfile(root, 'analysisConfig.json'), 'w');
fwrite(fid, jsonencode(cfg, 'PrettyPrint', true), 'char');
fclose(fid);

[ep, meta] = loadEpochSpecFromConfig(root, 'sz_t', 100, 'frameRate_s', 1);
testCase.verifyFalse(meta.usedDefault);
testCase.verifyEqual(size(ep.frameRange,1), 2);
testCase.verifyEqual(ep.epochNames{1}, 'B');
testCase.verifyEqual(ep.epochNames{2}, 'D');
end

function testJsonPresetDefaultAndExplicit(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));

cfg = struct();
cfg.defaultPreset = 'drug_application';
cfg.presets = struct();
cfg.presets.drug_application = struct('epochs', struct( ...
    'mode', 'frames', ...
    'frameRange', [10 20; 21 30], ...
    'epochNames', {{'B','D'}}));
cfg.presets.genotype_compare = struct('epochs', struct( ...
    'mode', 'all', ...
    'epochNames', {{'All'}}));

fid = fopen(fullfile(root, 'analysisConfig.json'), 'w');
fwrite(fid, jsonencode(cfg, 'PrettyPrint', true), 'char');
fclose(fid);

% Uses defaultPreset
[ep1, meta1] = loadEpochSpecFromConfig(root, 'sz_t', 100, 'frameRate_s', 1);
testCase.verifyEqual(size(ep1.frameRange,1), 2);
testCase.verifyEqual(meta1.presetName, 'drug_application');

% Explicit preset override
[ep2, meta2] = loadEpochSpecFromConfig(root, 'sz_t', 100, 'frameRate_s', 1, 'presetName', 'genotype_compare');
testCase.verifyEqual(size(ep2.frameRange,1), 1);
testCase.verifyEqual(meta2.presetName, 'genotype_compare');
end

function testJsonPresetMissingErrors(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));

cfg = struct();
cfg.defaultPreset = 'missing_name';
cfg.presets = struct();
cfg.presets.some_other = struct('epochs', struct('mode', 'all', 'epochNames', {{'All'}}));

fid = fopen(fullfile(root, 'analysisConfig.json'), 'w');
fwrite(fid, jsonencode(cfg, 'PrettyPrint', true), 'char');
fclose(fid);

didErr = false;
try
    loadEpochSpecFromConfig(root, 'sz_t', 100, 'frameRate_s', 1);
catch ME
    didErr = true;
    testCase.verifyTrue(contains(ME.identifier, 'PresetNotFound') || contains(ME.message, 'Preset'));
end
testCase.verifyTrue(didErr);
end

function testWriterParamFramesConfig(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));

writeEpochConfigTemplate(root, ...
    'overwrite', true, ...
    'includeTemplatePresets', false, ...
    'presetName', 'sample_epochs', ...
    'defaultPreset', 'sample_epochs', ...
    'mode', 'frames', ...
    'frameRange', [15 266; 281 532; 533 784], ...
    'epochNames', {'baseline','CNO1','CNO2'});

[ep, meta] = loadEpochSpecFromConfig(root, 'sz_t', 1000, 'frameRate_s', 1);
testCase.verifyFalse(meta.usedDefault);
testCase.verifyEqual(meta.presetName, 'sample_epochs');
testCase.verifyEqual(size(ep.frameRange,1), 3);
testCase.verifyEqual(ep.epochNames{1}, 'baseline');
testCase.verifyEqual(ep.epochNames{3}, 'CNO2');
end

function testFrameSplitConvenienceWriter(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));

writeEpochConfigFromFrameSplit(root, ...
    [15 266; 281 532; 533 784], ...
    {'baseline','CNO1','CNO2'}, ...
    'presetName', 'sample_epochs', ...
    'defaultPreset', 'sample_epochs', ...
    'overwrite', true);

[ep, meta] = loadEpochSpecFromConfig(root, 'sz_t', 1000, 'frameRate_s', 1);
testCase.verifyFalse(meta.usedDefault);
testCase.verifyEqual(meta.presetName, 'sample_epochs');
testCase.verifyEqual(size(ep.frameRange,1), 3);
testCase.verifyEqual(ep.epochNames{2}, 'CNO1');
end

function testWriterEnsuresPlanCompatiblePresetAliases(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() safeRmdir(root));

writeEpochConfigTemplate(root, ...
    'overwrite', true, ...
    'includeTemplatePresets', false, ...
    'ensurePlanCompatiblePresets', true, ...
    'presetName', 'sample_epochs', ...
    'defaultPreset', 'sample_epochs', ...
    'mode', 'frames', ...
    'frameRange', [10 20; 21 30], ...
    'epochNames', {'B','D'});

cfg = jsondecode(fileread(fullfile(root, 'analysisConfig.json')));
testCase.verifyTrue(isfield(cfg.presets, 'sample_epochs'));
testCase.verifyTrue(isfield(cfg.presets, 'drug_application'));
testCase.verifyTrue(isfield(cfg.presets, 'genotype_compare'));
end

function safeRmdir(pathIn)
if exist(pathIn, 'dir')
    rmdir(pathIn, 's');
end
end
