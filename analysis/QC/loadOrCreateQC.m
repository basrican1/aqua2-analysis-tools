function [qc, loaded, meta] = loadOrCreateQC(pathIn, varargin)
% Load existing QC tables or create defaults from canonical data.
%
% Inputs:
%   pathIn: analysis root (with BWA_AQ2.mat) or canonical run folder
%
% Options:
%   reviewer, reason, defaultKeep (passed to createQCTables when creating)

p = inputParser;
p.StructExpand = false;
addRequired(p, 'pathIn', @ischar);
addParameter(p, 'reviewer', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'reason', "default_keep", @(x) ischar(x) || isstring(x));
addParameter(p, 'defaultKeep', true, @islogical);
parse(p, pathIn, varargin{:});

rootDir = detectRootDir(pathIn);
loaded = loadClusterOutputs(pathIn, 'attachLegacyMeta', false, 'attachEventMeta', false);

filesPath = fullfile(rootDir, 'qc_files.mat');
cellsPath = fullfile(rootDir, 'qc_cells.mat');
eventsPath = fullfile(rootDir, 'qc_events.mat');

haveAll = exist(filesPath, 'file') == 2 && exist(cellsPath, 'file') == 2 && exist(eventsPath, 'file') == 2;
if haveAll
    s1 = load(filesPath); s2 = load(cellsPath); s3 = load(eventsPath);
    qc = struct();
    qc.qc_files = getVarOrFirstTable(s1, 'qc_files');
    qc.qc_cells = getVarOrFirstTable(s2, 'qc_cells');
    qc.qc_events = getVarOrFirstTable(s3, 'qc_events');
    source = 'loaded';
else
    qc = createQCTables(loaded.files, loaded.cells, loaded.events, ...
        'reviewer', p.Results.reviewer, 'reason', p.Results.reason, 'defaultKeep', p.Results.defaultKeep);
    source = 'created';
end

meta = struct('rootDir', rootDir, 'source', source, ...
    'filesPath', filesPath, 'cellsPath', cellsPath, 'eventsPath', eventsPath);
end

function rootDir = detectRootDir(pathIn)
if exist(fullfile(pathIn, 'BWA_AQ2.mat'), 'file') == 2
    rootDir = pathIn;
elseif exist(fullfile(pathIn, 'files.mat'), 'file') == 2 && exist(fullfile(pathIn, 'event_cell_map.mat'), 'file') == 2
    rootDir = fileparts(fileparts(fileparts(pathIn)));
else
    rootDir = pathIn;
end
end

function T = getVarOrFirstTable(S, name)
if isfield(S, name)
    T = S.(name);
    return
end
f = fieldnames(S);
T = table();
for i = 1:numel(f)
    v = S.(f{i});
    if istable(v)
        T = v;
        return
    end
end
end

