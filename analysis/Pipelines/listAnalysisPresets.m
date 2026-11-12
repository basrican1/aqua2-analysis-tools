function T = listAnalysisPresets(pathIn, varargin)
% List available analysis presets from analysisPlan config.
%
% Usage:
%   T = listAnalysisPresets(pathIn)
%   listAnalysisPresets(pathIn)  % also prints table

p = inputParser;
p.StructExpand = false;
addRequired(p, 'pathIn', @ischar);
addParameter(p, 'printTable', true, @islogical);
parse(p, pathIn, varargin{:});

cfg = loadPlanConfig(pathIn);
defaultPreset = "";
if isfield(cfg, 'defaultPreset') && ~isempty(cfg.defaultPreset)
    defaultPreset = string(cfg.defaultPreset);
end

if ~isfield(cfg, 'presets') || ~isstruct(cfg.presets)
    T = table();
    if p.Results.printTable
        fprintf('[listAnalysisPresets] No presets found.\\n');
    end
    return
end

names = string(fieldnames(cfg.presets));
n = numel(names);
isDefault = false(n,1);
epochPreset = strings(n,1);
nMetrics = zeros(n,1);
freqGroupBy = strings(n,1);
metricGroupBy = strings(n,1);
dataLimitType = strings(n,1);
dataLimitPercentiles = strings(n,1);

for i = 1:n
    nm = names(i);
    P = cfg.presets.(char(nm));
    isDefault(i) = (nm == defaultPreset);
    epochPreset(i) = string(getFieldOr(P, 'epochPreset', ""));

    m = getFieldOr(P, 'metrics', {});
    if isstring(m)
        nMetrics(i) = numel(m);
    elseif iscell(m)
        nMetrics(i) = numel(m);
    else
        nMetrics(i) = 0;
    end

    fgb = getFieldOr(P, 'freqGroupBy', {});
    mgb = getFieldOr(P, 'metricGroupBy', {});
    freqGroupBy(i) = joinStrings(fgb, ',');
    metricGroupBy(i) = joinStrings(mgb, ',');
    dataLimitType(i) = string(getFieldOr(P, 'dataLimitType', ""));
    dlp = getFieldOr(P, 'dataLimitPercentiles', []);
    if isnumeric(dlp) && numel(dlp) == 2
        dataLimitPercentiles(i) = string(sprintf('[%g %g]', dlp(1), dlp(2)));
    else
        dataLimitPercentiles(i) = "";
    end
end

T = table(names, isDefault, epochPreset, nMetrics, freqGroupBy, metricGroupBy, dataLimitType, dataLimitPercentiles, ...
    'VariableNames', {'preset','is_default','epoch_preset','n_metrics','freq_group_by','metric_group_by','data_limit_type','data_limit_percentiles'});

if p.Results.printTable
    disp(T);
end
end

function cfg = loadPlanConfig(pathIn)
rootDir = detectRootDir(pathIn);
jsonPath = fullfile(rootDir, 'analysisPlan.json');
matPath = fullfile(rootDir, 'analysisPlan.mat');

if exist(jsonPath, 'file') == 2
    cfg = jsondecode(fileread(jsonPath));
    return
end
if exist(matPath, 'file') == 2
    S = load(matPath);
    if isfield(S, 'analysisPlan')
        cfg = S.analysisPlan;
    else
        f = fieldnames(S);
        if isempty(f)
            cfg = struct();
        else
            cfg = S.(f{1});
        end
    end
    return
end
error('listAnalysisPresets:MissingConfig', 'No analysisPlan.json/.mat found in %s', rootDir);
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

function v = getFieldOr(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = d;
end
end

function out = joinStrings(v, delim)
if nargin < 2, delim = ','; end
if isempty(v)
    out = "";
elseif isstring(v)
    out = strjoin(cellstr(v(:)'), delim);
elseif iscell(v)
    out = strjoin(cellstr(string(v(:)')), delim);
else
    out = string(v);
end
out = string(out);
end
