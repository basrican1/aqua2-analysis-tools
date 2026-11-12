function outPath = writeAnalysisPlanPreset(pathIn, varargin)
% Create or update one preset in analysisPlan.json without manual editing.
%
% Example:
%   writeAnalysisPlanPreset(pathIn, ...
%       'presetName', 'drug_epoch_compare', ...
%       'setDefault', true, ...
%       'epochPreset', 'drug_application', ...
%       'dataLimitType', 'duration_frames', ...
%       'dataLimitPercentiles', [10 95]);

p = inputParser;
p.StructExpand = false;
addRequired(p, 'pathIn', @ischar);
addParameter(p, 'presetName', 'drug_epoch_compare', @(x) ischar(x) || isstring(x));
addParameter(p, 'setDefault', true, @islogical);
addParameter(p, 'epochPreset', 'drug_application', @(x) ischar(x) || isstring(x));
addParameter(p, 'subset', struct(), @isstruct);
addParameter(p, 'metrics', {'area_px2','duration_frames','rise19_frames','fall91_frames','decay_tau_frames','dff_auc','dff_max2','max_prop_speed','perimeter_px','circularity','propagation_distance_px'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'freqGroupBy', {'file_id'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'metricGroupBy', {'animal_id','slice_num','file_id'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'dataLimitType', 'area_px2', @(x) ischar(x) || isstring(x));
addParameter(p, 'dataLimitPercentiles', [0 100], @isnumeric);
addParameter(p, 'dataLimitRules', {}, @(x) iscell(x) || isempty(x));
addParameter(p, 'areaPercentiles', [0 100], @isnumeric); % backward-compatible field
addParameter(p, 'makePlots', false, @islogical);
addParameter(p, 'strictQC', false, @islogical);
parse(p, pathIn, varargin{:});

rootDir = detectRootDir(pathIn);
outPath = fullfile(rootDir, 'analysisPlan.json');
cfg = struct('defaultPreset', '', 'presets', struct());
if exist(outPath, 'file') == 2
    cfg = jsondecode(fileread(outPath));
    if ~isfield(cfg, 'presets') || ~isstruct(cfg.presets)
        cfg.presets = struct();
    end
end

presetName = char(string(p.Results.presetName));
P = struct();
P.epochPreset = char(string(p.Results.epochPreset));
P.subset = p.Results.subset;
P.metrics = cellstr(string(p.Results.metrics(:)'));
P.freqGroupBy = cellstr(string(p.Results.freqGroupBy(:)'));
P.metricGroupBy = cellstr(string(p.Results.metricGroupBy(:)'));
P.dataLimitType = char(string(p.Results.dataLimitType));
P.dataLimitPercentiles = p.Results.dataLimitPercentiles;
P.dataLimitRules = p.Results.dataLimitRules;
P.areaPercentiles = p.Results.areaPercentiles;
P.makePlots = logical(p.Results.makePlots);
P.strictQC = logical(p.Results.strictQC);

cfg.presets.(presetName) = P;
if p.Results.setDefault || ~isfield(cfg, 'defaultPreset') || isempty(cfg.defaultPreset)
    cfg.defaultPreset = presetName;
end

fid = fopen(outPath, 'w');
if fid < 0
    error('writeAnalysisPlanPreset:OpenFail', 'Could not write: %s', outPath);
end
fwrite(fid, jsonencode(cfg, 'PrettyPrint', true), 'char');
fclose(fid);
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
