function [plan, meta] = loadAnalysisPlan(pathIn, varargin)
% Load analysis plan and select one preset.
%
% Config discovery:
%   <root>/analysisPlan.json
%   <root>/analysisPlan.mat (variable analysisPlan or first variable)
%
% Plan shape:
%   defaultPreset: "<name>"
%   presets.<name> = struct with fields used by runAnalysisPreset.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'pathIn', @ischar);
addParameter(p, 'presetName', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'jsonName', 'analysisPlan.json', @ischar);
addParameter(p, 'matName', 'analysisPlan.mat', @ischar);
parse(p, pathIn, varargin{:});

rootDir = detectRootDir(pathIn);
jsonPath = fullfile(rootDir, p.Results.jsonName);
matPath = fullfile(rootDir, p.Results.matName);

cfg = struct();
sourceFile = '';
if exist(jsonPath, 'file') == 2
    cfg = jsondecode(fileread(jsonPath));
    sourceFile = jsonPath;
elseif exist(matPath, 'file') == 2
    S = load(matPath);
    if isfield(S, 'analysisPlan')
        cfg = S.analysisPlan;
    else
        f = fieldnames(S);
        if ~isempty(f), cfg = S.(f{1}); end
    end
    sourceFile = matPath;
else
    error('loadAnalysisPlan:MissingConfig', 'No analysisPlan.json/.mat found under: %s', rootDir);
end

preset = char(string(p.Results.presetName));
if isempty(strtrim(preset))
    if isfield(cfg, 'defaultPreset') && ~isempty(cfg.defaultPreset)
        preset = char(string(cfg.defaultPreset));
    else
        error('loadAnalysisPlan:MissingPreset', 'No presetName provided and defaultPreset missing.');
    end
end

if ~isfield(cfg, 'presets') || ~isstruct(cfg.presets) || ~isfield(cfg.presets, preset)
    error('loadAnalysisPlan:PresetNotFound', 'Preset "%s" not found in analysis plan.', preset);
end

plan = cfg.presets.(preset);
plan.presetName = preset;

meta = struct();
meta.rootDir = rootDir;
meta.sourceFile = sourceFile;
meta.presetName = preset;
end

function rootDir = detectRootDir(pathIn)
if exist(fullfile(pathIn, 'BWA_AQ2.mat'), 'file') == 2
    rootDir = pathIn; return
end
if exist(fullfile(pathIn, 'files.mat'), 'file') == 2 && exist(fullfile(pathIn, 'event_cell_map.mat'), 'file') == 2
    rootDir = fileparts(fileparts(fileparts(pathIn))); return
end
rootDir = pathIn;
end

