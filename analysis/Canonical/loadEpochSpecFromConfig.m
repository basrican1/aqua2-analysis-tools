function [epochSpec, meta] = loadEpochSpecFromConfig(pathIn, varargin)
% Load epoch definition from config file near BWA_AQ2.mat.
%
% Config discovery order:
%   1) <root>/analysisConfig.json
%   2) <root>/analysisConfig.mat (variable: analysisConfig or epochs)
%
% Supported config shape:
%   top-level struct with field "epochs", or epoch fields at top-level.
%   optionally with presets:
%     presets.<name>.epochs = struct(...)
%     defaultPreset = "<name>"
%
% Epoch fields:
%   mode: 'frames' | 'times' | 'all'
%   frameRange: [N x 2]
%   timeRangeMin (or timeRange): [N x 2] in minutes
%   epochNames: cellstr/string array

p = inputParser;
p.StructExpand = false;
addRequired(p, 'pathIn', @ischar);
addParameter(p, 'sz_t', NaN, @isnumeric);
addParameter(p, 'frameRate_s', NaN, @isnumeric);
addParameter(p, 'jsonName', 'analysisConfig.json', @ischar);
addParameter(p, 'matName', 'analysisConfig.mat', @ischar);
addParameter(p, 'defaultEpochNames', {'All'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'presetName', '', @(x) ischar(x) || isstring(x));
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
    if isfield(S, 'analysisConfig')
        cfg = S.analysisConfig;
    elseif isfield(S, 'epochs')
        cfg = struct('epochs', S.epochs);
    else
        f = fieldnames(S);
        if ~isempty(f), cfg = S.(f{1}); end
    end
    sourceFile = matPath;
end

[sz_t, frameRate] = resolveSizes(p.Results.sz_t, p.Results.frameRate_s, rootDir);
if ~(isfinite(sz_t) && sz_t > 0), sz_t = 1; end
if ~(isfinite(frameRate) && frameRate > 0), frameRate = 1; end

usedDefault = false;
preset = '';
if isempty(sourceFile)
    epochSpec = frameTimeConversionsByEpoch([1 inf], [], sz_t, frameRate, p.Results.defaultEpochNames);
    mode = 'all';
    usedDefault = true;
else
    preset = resolvePresetName(cfg, p.Results.presetName);
    E = selectEpochConfig(cfg, preset);
    mode = lower(string(getFieldOr(E, 'mode', 'all')));
    names = getFieldOr(E, 'epochNames', p.Results.defaultEpochNames);
    if isstring(names), names = cellstr(names(:)); end

    switch char(mode)
        case 'frames'
            framesIn = toNumeric2(getFieldOr(E, 'frameRange', []));
            if isempty(framesIn)
                error('loadEpochSpecFromConfig:BadConfig', 'mode=frames requires frameRange.');
            end
            epochSpec = frameTimeConversionsByEpoch(framesIn, [], sz_t, frameRate, names);
        case 'times'
            timesIn = toNumeric2(getFieldOr(E, 'timeRangeMin', getFieldOr(E, 'timeRange', [])));
            if isempty(timesIn)
                error('loadEpochSpecFromConfig:BadConfig', 'mode=times requires timeRangeMin/timeRange.');
            end
            epochSpec = frameTimeConversionsByEpoch([], timesIn, sz_t, frameRate, names);
        case {'all','none'}
            epochSpec = frameTimeConversionsByEpoch([1 inf], [], sz_t, frameRate, names);
        otherwise
            error('loadEpochSpecFromConfig:BadMode', 'Unknown epoch mode: %s', char(mode));
    end
end

meta = struct();
meta.rootDir = rootDir;
meta.sourceFile = sourceFile;
meta.mode = char(mode);
meta.usedDefault = usedDefault;
meta.presetName = char(string(preset));
meta.sz_t = sz_t;
meta.frameRate_s = frameRate;
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

function [sz_t, frameRate] = resolveSizes(szIn, frIn, rootDir)
sz_t = szIn; frameRate = frIn;
if isfinite(sz_t) && sz_t > 0 && isfinite(frameRate) && frameRate > 0
    return
end
try
    L = loadClusterOutputs(rootDir, 'attachLegacyMeta', false, 'attachEventMeta', false);
    if ~(isfinite(sz_t) && sz_t > 0) && ismember('sz_t', L.files.Properties.VariableNames) && ~isempty(L.files)
        sz_t = double(L.files.sz_t(1));
    end
    if ~(isfinite(frameRate) && frameRate > 0) && ismember('frame_rate_s', L.files.Properties.VariableNames) && ~isempty(L.files)
        frameRate = double(L.files.frame_rate_s(1));
    end
catch
end
end

function E = selectEpochConfig(cfg, preset)
if ~isempty(preset)
    if ~isfield(cfg, 'presets')
        error('loadEpochSpecFromConfig:MissingPresets', 'Requested preset "%s" but config has no presets field.', preset);
    end
    if ~isfield(cfg.presets, preset)
        error('loadEpochSpecFromConfig:PresetNotFound', 'Preset "%s" not found in config.presets.', preset);
    end
    P = cfg.presets.(preset);
    if isfield(P, 'epochs'), E = P.epochs; else, E = P; end
    return
end
if isfield(cfg, 'epochs'), E = cfg.epochs; else, E = cfg; end
end

function preset = resolvePresetName(cfg, requestedPreset)
preset = char(string(requestedPreset));
if ~isempty(strtrim(preset)), return; end
preset = '';
if isstruct(cfg) && isfield(cfg, 'defaultPreset') && ~isempty(cfg.defaultPreset)
    preset = char(string(cfg.defaultPreset));
end
end

function v = getFieldOr(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function M = toNumeric2(v)
if isempty(v)
    M = [];
elseif isnumeric(v)
    M = double(v);
elseif iscell(v)
    try, M = cell2mat(v); catch, M = []; end
else
    M = [];
end
if ~isempty(M) && size(M,2) ~= 2
    error('loadEpochSpecFromConfig:BadShape', 'Epoch ranges must be N x 2.');
end
end
