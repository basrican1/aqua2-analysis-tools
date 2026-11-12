function out = loadClusterOutputs(pathIn, varargin)
    % Load canonical tables from a cluster run output, with optional legacy payload.
    %
    % Inputs:
    %   pathIn: either
    %     1) run folder: .../derived/ingest_runs/<timestamp>
    %     2) analysis root: folder containing derived/ingest_runs, optionally BWA_AQ2.mat
    %
    % Options:
    %   'attachLegacyMeta' (default true): rebuild BWA_AQ2.opts/DB from files table.
    %   'attachEventMeta' (default true): rebuild BWA_AQ2.fts1.loc.t0/t1 and fts1.fileID from events/files tables.
    %   'runDir' (default ''): explicit canonical run folder override.

    p = inputParser;
    p.StructExpand = false;
    addRequired(p, 'pathIn', @ischar);
    addParameter(p, 'attachLegacyMeta', true, @islogical);
    addParameter(p, 'attachEventMeta', true, @islogical);
    addParameter(p, 'runDir', '', @ischar);
    parse(p, pathIn, varargin{:});

    attachLegacyMeta = p.Results.attachLegacyMeta;
    attachEventMeta = p.Results.attachEventMeta;
    runDirOverride = p.Results.runDir;

    inputPath = pathIn;
    if startsWith(inputPath, '~')
        homeDir = char(java.lang.System.getProperty('user.home'));
        rel = char(extractAfter(string(inputPath), 1));
        if startsWith(rel, filesep)
            rel = rel(2:end);
        end
        inputPath = fullfile(homeDir, rel);
    end

    if ~isempty(runDirOverride)
        runDir = runDirOverride;
    else
        runDir = detectRunDir(inputPath);
    end

    rootDir = detectRootDir(inputPath, runDir);

    bwaPath = fullfile(rootDir, 'BWA_AQ2.mat');

    filesPath = fullfile(runDir, 'files.mat');
    cellsPath = fullfile(runDir, 'cells.mat');
    eventsPath = fullfile(runDir, 'events.mat');
    mapPath = fullfile(runDir, 'event_cell_map.mat');
    if ~exist(filesPath, 'file') || ~exist(cellsPath, 'file') || ~exist(eventsPath, 'file') || ~exist(mapPath, 'file')
        error('Run dir missing one or more canonical tables: %s', runDir);
    end

    out.BWA_AQ2 = struct();
    out.hasLegacyPayload = exist(bwaPath, 'file') == 2;
    out.legacyPayloadPath = bwaPath;
    if out.hasLegacyPayload
        sBwa = load(bwaPath);
        if ~isfield(sBwa, 'BWA_AQ2')
            error('BWA_AQ2.mat does not contain variable BWA_AQ2');
        end
        out.BWA_AQ2 = sBwa.BWA_AQ2;
    end

    sFiles = load(filesPath);
    sCells = load(cellsPath);
    sEvents = load(eventsPath);
    sMap = load(mapPath);
    out.files = sFiles.files;
    out.cells = sCells.cells;
    out.events = sEvents.events;
    out.event_cell_map = sMap.event_cell_map;
    out.files = ensureDerivedFileColumns(out.files);
    out.runDir = runDir;
    out.rootDir = rootDir;

    if attachLegacyMeta && out.hasLegacyPayload
        out.BWA_AQ2 = attachLegacyMetaFromCanonical(out.BWA_AQ2, out.files);
    end
    if attachEventMeta && out.hasLegacyPayload
        out.BWA_AQ2 = attachEventFieldsFromCanonical(out.BWA_AQ2, out.events, out.files);
    end
end

function files = ensureDerivedFileColumns(files)
    if ~istable(files)
        return
    end
    % Provide a normalized treatment column for downstream compare modes.
    if ~ismember('treatment', files.Properties.VariableNames)
        treatment = repmat("", height(files), 1);
        if ismember('drug', files.Properties.VariableNames)
            d = string(files.drug);
            treatment(d ~= "") = d(d ~= "");
        end
        if ismember('virus', files.Properties.VariableNames)
            v = string(files.virus);
            fillMask = treatment == "" & v ~= "";
            treatment(fillMask) = v(fillMask);
        end
        files.treatment = treatment;
    else
        files.treatment = string(files.treatment);
    end
end

function runDir = detectRunDir(inputPath)
    if exist(fullfile(inputPath, 'files.mat'), 'file') && exist(fullfile(inputPath, 'event_cell_map.mat'), 'file')
        runDir = inputPath;
        return
    end
    latestPtr = fullfile(inputPath, 'derived', 'ingest_runs', 'latest_ingest.txt');
    if exist(latestPtr, 'file')
        txt = strtrim(fileread(latestPtr));
        if exist(txt, 'dir')
            runDir = txt;
            return
        end
        rel = fullfile(inputPath, txt);
        if exist(rel, 'dir')
            runDir = rel;
            return
        end
    end
    error('Could not infer canonical run dir from path: %s', inputPath);
end

function rootDir = detectRootDir(inputPath, runDir)
    if exist(fullfile(inputPath, 'BWA_AQ2.mat'), 'file')
        rootDir = inputPath;
        return
    end
    % runDir expected at <root>/derived/ingest_runs/<timestamp>
    rootDir = fileparts(fileparts(fileparts(runDir)));
end
