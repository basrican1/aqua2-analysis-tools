function out = saveQC(pathIn, qc, varargin)
% Save QC tables to analysis root, optionally with versioned snapshot.
%
% Saves:
%   qc_files.mat, qc_cells.mat, qc_events.mat
%
% Optional versioned snapshot:
%   qc_versions/qc_YYYYMMDD_HHMMSS/{qc_files.mat,qc_cells.mat,qc_events.mat}

p = inputParser;
p.StructExpand = false;
addRequired(p, 'pathIn', @ischar);
addRequired(p, 'qc', @isstruct);
addParameter(p, 'versioned', true, @islogical);
parse(p, pathIn, qc, varargin{:});

rootDir = detectRootDir(pathIn);
if ~exist(rootDir, 'dir')
    error('saveQC:RootMissing', 'Root directory not found: %s', rootDir);
end

requireField(qc, 'qc_files');
requireField(qc, 'qc_cells');
requireField(qc, 'qc_events');

qc_files = qc.qc_files; %#ok<NASGU>
qc_cells = qc.qc_cells; %#ok<NASGU>
qc_events = qc.qc_events; %#ok<NASGU>

filesPath = fullfile(rootDir, 'qc_files.mat');
cellsPath = fullfile(rootDir, 'qc_cells.mat');
eventsPath = fullfile(rootDir, 'qc_events.mat');

save(filesPath, 'qc_files');
save(cellsPath, 'qc_cells');
save(eventsPath, 'qc_events');

out = struct();
out.rootDir = rootDir;
out.filesPath = filesPath;
out.cellsPath = cellsPath;
out.eventsPath = eventsPath;
out.versionDir = '';

if p.Results.versioned
    stamp = datestr(now, 'yyyymmdd_HHMMSS');
    verBase = fullfile(rootDir, 'qc_versions');
    verDir = fullfile(verBase, ['qc_' stamp]);
    if ~exist(verDir, 'dir')
        mkdir(verDir);
    end
    save(fullfile(verDir, 'qc_files.mat'), 'qc_files');
    save(fullfile(verDir, 'qc_cells.mat'), 'qc_cells');
    save(fullfile(verDir, 'qc_events.mat'), 'qc_events');

    fid = fopen(fullfile(verBase, 'latest_qc.txt'), 'w');
    if fid > 0
        fprintf(fid, '%s\n', verDir);
        fclose(fid);
    end
    out.versionDir = verDir;
end
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

function requireField(S, f)
if ~isfield(S, f) || ~istable(S.(f))
    error('saveQC:BadQC', 'qc.%s missing or not a table', f);
end
end

