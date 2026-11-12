function outPath = writeEpochConfigTemplate(pathIn, varargin)
% Write analysisConfig.json near BWA_AQ2.mat.
%
% Default behavior (no extra args) writes a starter template with
% drug/genotype example presets.
%
% Param-driven behavior example:
%   writeEpochConfigTemplate(pathIn, ...
%       'overwrite', true, ...
%       'includeTemplatePresets', false, ...
%       'presetName', 'sample_epochs', ...
%       'defaultPreset', 'sample_epochs', ...
%       'mode', 'frames', ...
%       'frameRange', [15 266; 281 532; 533 784], ...
%       'epochNames', {'baseline','ExperimentalEpoch1','ExperimentalEpoch2'});

p = inputParser;
p.StructExpand = false;
addRequired(p, 'pathIn', @ischar);
addParameter(p, 'overwrite', false, @islogical);
addParameter(p, 'includeTemplatePresets', true, @islogical);
addParameter(p, 'ensurePlanCompatiblePresets', true, @islogical);
addParameter(p, 'presetName', 'drug_application', @(x) ischar(x) || isstring(x));
addParameter(p, 'defaultPreset', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'mode', 'frames', @(x) any(strcmpi(string(x), ["frames","times","all"])));
addParameter(p, 'frameRange', [1, 300; 301, 600], @isnumeric);
addParameter(p, 'timeRangeMin', [], @isnumeric);
addParameter(p, 'epochNames', {'Baseline', 'Drug'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'notes', '', @(x) ischar(x) || isstring(x));
parse(p, pathIn, varargin{:});

rootDir = pathIn;
if exist(fullfile(pathIn, 'BWA_AQ2.mat'), 'file') ~= 2
    if exist(fullfile(pathIn, 'files.mat'), 'file') == 2 && exist(fullfile(pathIn, 'event_cell_map.mat'), 'file') == 2
        rootDir = fileparts(fileparts(fileparts(pathIn)));
    end
end

outPath = fullfile(rootDir, 'analysisConfig.json');
if exist(outPath, 'file') == 2 && ~p.Results.overwrite
    error('writeEpochConfigTemplate:Exists', 'File exists: %s (use overwrite=true)', outPath);
end

cfg = struct();
cfg.presets = struct();

if p.Results.includeTemplatePresets
    cfg.presets.drug_application = struct();
    cfg.presets.drug_application.epochs = struct();
    cfg.presets.drug_application.epochs.mode = 'frames'; % 'frames' | 'times' | 'all'
    cfg.presets.drug_application.epochs.frameRange = [1, 300; 301, 600];
    cfg.presets.drug_application.epochs.epochNames = {'Baseline', 'Drug'};

    cfg.presets.genotype_compare = struct();
    cfg.presets.genotype_compare.epochs = struct();
    cfg.presets.genotype_compare.epochs.mode = 'all';
    cfg.presets.genotype_compare.epochs.epochNames = {'All'};
end

presetName = char(string(p.Results.presetName));
epochs = buildEpochStruct(p.Results.mode, p.Results.frameRange, p.Results.timeRangeMin, p.Results.epochNames);
cfg.presets.(presetName) = struct('epochs', epochs);

if p.Results.ensurePlanCompatiblePresets
    if ~isfield(cfg.presets, 'drug_application')
        % Alias to a frame/time/all epoch spec compatible with plan templates.
        cfg.presets.drug_application = struct('epochs', epochs);
    end
    if ~isfield(cfg.presets, 'genotype_compare')
        cfg.presets.genotype_compare = struct('epochs', struct('mode', 'all', 'epochNames', {{'All'}}));
    end
end

defPreset = char(string(p.Results.defaultPreset));
if isempty(strtrim(defPreset))
    if p.Results.includeTemplatePresets
        defPreset = 'drug_application';
    else
        defPreset = presetName;
    end
end
cfg.defaultPreset = defPreset;

noteTxt = char(string(p.Results.notes));
if isempty(strtrim(noteTxt))
    noteTxt = 'Set defaultPreset or pass presetName to loadEpochSpecFromConfig(pathIn, ''presetName'', ''name'').';
end
cfg.notes = noteTxt;

fid = fopen(outPath, 'w');
if fid < 0
    error('writeEpochConfigTemplate:OpenFail', 'Could not write: %s', outPath);
end
fwrite(fid, jsonencode(cfg, 'PrettyPrint', true), 'char');
fclose(fid);
end

function epochs = buildEpochStruct(modeIn, frameRangeIn, timeRangeIn, epochNamesIn)
mode = lower(char(string(modeIn)));
epochs = struct();
epochs.mode = mode;

switch mode
    case 'frames'
        frameRange = toNx2(frameRangeIn, 'frameRange');
        epochs.frameRange = frameRange;
        n = size(frameRange, 1);
    case 'times'
        timeRange = toNx2(timeRangeIn, 'timeRangeMin');
        epochs.timeRangeMin = timeRange;
        n = size(timeRange, 1);
    case 'all'
        n = 1;
    otherwise
        error('writeEpochConfigTemplate:BadMode', 'Unsupported mode: %s', mode);
end

epochNames = cellstr(string(epochNamesIn(:)));
if isempty(epochNames)
    epochNames = autoEpochNames(n);
elseif numel(epochNames) ~= n
    error('writeEpochConfigTemplate:NameCountMismatch', ...
        'epochNames count (%d) must match number of epochs (%d).', numel(epochNames), n);
end
epochs.epochNames = epochNames;
end

function M = toNx2(v, fieldName)
if isempty(v) || ~isnumeric(v)
    error('writeEpochConfigTemplate:MissingRange', '%s must be numeric N x 2.', fieldName);
end
M = double(v);
if size(M,2) ~= 2
    error('writeEpochConfigTemplate:BadShape', '%s must be numeric N x 2.', fieldName);
end
end

function names = autoEpochNames(n)
names = cell(n,1);
for ii = 1:n
    names{ii} = ['Epoch' num2str(ii)];
end
end
