function out = runLegacyComputeThenPlot(pathIn, cfg)
% Config-driven compute-first legacy pipeline with optional plotting.
%
% This wrapper separates workflow into:
%   1) load canonical + BWA_AQ2
%   2) optional QC + subset on canonical tables
%   3) optional BWA_AQ2 reduction to selected files
%   4) compute modules
%   5) plotting modules (optional)
%
% Usage:
%   cfg = struct();
%   cfg.modules = {'standard','polars','baviolins'};
%   cfg.doPlots = true;
%   out = runLegacyComputeThenPlot(pathIn, cfg);

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

if nargin < 2 || isempty(cfg)
    cfg = struct();
end
if ~ischar(pathIn) && ~(isstring(pathIn) && isscalar(pathIn))
    error('runLegacyComputeThenPlot:BadPath', 'pathIn must be char or scalar string.');
end
pathIn = char(string(pathIn));

d = defaultCfg();
cfg = mergeStructs(d, cfg);

loaded = loadClusterOutputs(pathIn, ...
    'attachLegacyMeta', false, ...
    'attachEventMeta', false, ...
    'runDir', char(string(cfg.runDir)));

if cfg.useQC
    qc_files = cfg.qc_files;
    qc_cells = cfg.qc_cells;
    qc_events = cfg.qc_events;
    if isempty(qc_files) && isempty(qc_cells) && isempty(qc_events)
        qc = createQCTables(loaded.files, loaded.cells, loaded.events, ...
            'defaultKeep', true, 'reason', 'default_keep');
        qc_files = qc.qc_files;
        qc_cells = qc.qc_cells;
        qc_events = qc.qc_events;
    end
    qcApplied = applyQCToCanonical(loaded.files, loaded.cells, loaded.events, loaded.event_cell_map, ...
        qc_files, qc_cells, qc_events, 'strictQC', logical(cfg.strictQC));
else
    qc_files = table();
    qc_cells = table();
    qc_events = table();
    qcApplied = struct();
    qcApplied.files = loaded.files;
    qcApplied.cells = loaded.cells;
    qcApplied.events = loaded.events;
    qcApplied.event_cell_map = loaded.event_cell_map;
    qcApplied.report = struct();
end

subsetArgs = structToNameValue(cfg.subset);
S = selectCanonicalSubset(qcApplied.files, qcApplied.cells, qcApplied.events, qcApplied.event_cell_map, subsetArgs{:});
selectedFileIDs = cellstr(string(S.files.file_id));

BWA_use = loaded.BWA_AQ2;
if cfg.reduceToSubsetFiles
    if isempty(selectedFileIDs)
        warning('runLegacyComputeThenPlot:EmptySubset', 'Subset produced 0 files; using original BWA_AQ2.');
    else
        BWA_use = ReduceAQ2ByFile(BWA_use, selectedFileIDs, 'filesTbl', qcApplied.files, 'eventsTbl', qcApplied.events);
    end
end

tfre = resolveEpochSpec(pathIn, cfg);
ctx = prepareLegacyAqua2Context(BWA_use, ...
    'timeFrameRangeEpochs', tfre, ...
    'timeRange', cfg.timeRange, ...
    'frameRange', cfg.frameRange, ...
    'EpochList', cfg.EpochList, ...
    'memberlist', cfg.memberlist, ...
    'associatedCoverageTraces', cfg.associatedCoverageTraces, ...
    'memberCoverage', cfg.memberCoverage, ...
    'filesTbl', S.files, ...
    'eventsTbl', S.events, ...
    'eventCellMapTbl', S.event_cell_map);

moduleList = cellstr(string(cfg.modules));
out = struct();
out.config = cfg;
out.loaded = loaded;
out.qc = struct('qc_files', qc_files, 'qc_cells', qc_cells, 'qc_events', qc_events, 'applied', qcApplied);
out.subset = S;
out.context = ctx;
out.modules = struct();

for i = 1:numel(moduleList)
    m = lower(strtrim(moduleList{i}));
    switch m
        case {'standard','bwa_aqua2_standard'}
            modOut = struct();
            if cfg.returnCompute
                modOut.context = ctx;
                try
                    modOut.compute = computeBWA_Aqua2StandardData(ctx);
                catch ME
                    modOut.compute = struct('status', 'error', 'error', ME.message);
                end
            end
            if cfg.doPlots
                modOut.figures = plotBWA_Aqua2StandardPanels(ctx, ...
                    'style', cfg.plotStyle, ...
                    'statsConfig', cfg.statsConfig);
            end
            out.modules.standard = modOut;

        case {'polars','bwa_polars','bwa_aqua2_polars'}
            polData = computeBWA_Aqua2PolarsBySliceData(ctx.BWA_AQ2, ...
                'timeFrameRangeEpochs', ctx.timeFrameRangeEpochs, ...
                'memberlist', ctx.memberlist, ...
                'associatedCoverageTraces', ctx.associatedCoverageTraces, ...
                'memberCoverage', ctx.memberCoverage, ...
                'filesTbl', S.files, ...
                'eventsTbl', S.events, ...
                'includeStandardFeatures', logical(cfg.includeStandardFeatures));

            modOut = struct();
            if cfg.returnCompute
                modOut.compute = polData;
            end
            if cfg.doPlots
                modOut.figures = plotBWA_Aqua2PolarsBySlice(polData, ...
                    'spiderInputScale', cfg.spiderInputScale, ...
                    'spiderSymmetricLimit', cfg.spiderSymmetricLimit, ...
                    'exportCompanionTable', logical(cfg.exportCompanionTable), ...
                    'companionTableVarName', cfg.companionTableVarName, ...
                    'style', cfg.plotStyle);
            end
            out.modules.polars = modOut;

        case {'baviolins','ba_violin_batch'}
            trees = cellstr(string(cfg.baViolinDatTrees));
            batch = struct();
            for j = 1:numel(trees)
                datTree = trees{j};
                key = matlab.lang.makeValidName(strrep(datTree, '.', '_'));
                cdat = computeBAviolinData(ctx.BWA_AQ2.fts1, ctx.opts, datTree, ctx.timeFrameRangeEpochs, ctx.memberlistRow3, ...
                    'eventsTbl', ctx.eventsTbl, 'filesTbl', ctx.filesTbl);
                entry = struct();
                if cfg.returnCompute
                    entry.compute = cdat;
                end
                if cfg.doPlots
                    [h_fig, means_alln, means_celln] = plotBAviolinPanels(cdat, ctx.timeFrameRangeEpochs, ctx.memberlistRow3, ...
                        'optFigSuffix', cfg.baOptFigSuffix, ...
                        'doByCellOnly', logical(cfg.baDoByCellOnly), ...
                        'style', cfg.plotStyle, ...
                        'statsConfig', cfg.statsConfig);
                    entry.plot = struct('figure', h_fig, 'means_alln', means_alln, 'means_celln', means_celln);
                end
                batch.(key) = entry;
            end
            out.modules.baviolins = batch;

        otherwise
            error('runLegacyComputeThenPlot:UnknownModule', 'Unknown module: %s', m);
    end
end
end

function tfre = resolveEpochSpec(pathIn, cfg)
tfre = struct();
if isstruct(cfg.epochSpec) && ~isempty(fieldnames(cfg.epochSpec))
    tfre = cfg.epochSpec;
    return
end
if ~isempty(cfg.epochPreset)
    [tfre, ~] = loadEpochSpecFromConfig(pathIn, 'presetName', char(string(cfg.epochPreset)));
end
end

function nv = structToNameValue(s)
if isempty(s)
    nv = {};
    return
end
f = fieldnames(s);
nv = cell(1, 2*numel(f));
for i = 1:numel(f)
    nv{2*i-1} = f{i};
    nv{2*i} = s.(f{i});
end
end

function out = mergeStructs(defaults, in)
out = defaults;
if isempty(in)
    return
end
f = fieldnames(in);
for i = 1:numel(f)
    out.(f{i}) = in.(f{i});
end
end

function cfg = defaultCfg()
cfg = struct();
cfg.runDir = '';
cfg.useQC = false;
cfg.qc_files = table();
cfg.qc_cells = table();
cfg.qc_events = table();
cfg.strictQC = false;

cfg.subset = struct();
cfg.reduceToSubsetFiles = true;

cfg.epochSpec = struct();
cfg.epochPreset = '';
cfg.timeRange = [];
cfg.frameRange = [];
cfg.EpochList = {};

cfg.memberlist = [];
cfg.associatedCoverageTraces = {};
cfg.memberCoverage = [];

cfg.modules = {'standard','polars'};
cfg.doPlots = true;
cfg.returnCompute = true;
cfg.plotStyle = struct();
cfg.statsConfig = struct();

cfg.includeStandardFeatures = true;
cfg.spiderInputScale = 'fold';
cfg.spiderSymmetricLimit = [];
cfg.exportCompanionTable = true;
cfg.companionTableVarName = 'BWA_spiderFoldTable';

cfg.baViolinDatTrees = { ...
    'basic.area', ...
    'curve.duration', ...
    'curve.rise19', ...
    'curve.fall91', ...
    'curve.decayTau', ...
    'curve.dffAUC', ...
    'curve.dffMax2', ...
    'propagation.maxPropSpeed'};
cfg.baDoByCellOnly = false;
cfg.baOptFigSuffix = '';
end
