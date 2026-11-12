function out = compute_local_file_metrics(source, datTree, varargin)
% Compute file-level BAviolin metrics without necessarily plotting.
result = load_local_analysis_result(source);
ctx = get_derived_context(result);
p = inputParser;
p.StructExpand = false;
addOptional(p, 'datTree', datTree, @(x) ischar(x) || isstring(x));
addParameter(p, 'fileIDs', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'percentiles', [0 100], @isnumeric);
addParameter(p, 'plotFigures', false, @islogical);
parse(p, datTree, varargin{:});

metric = char(string(p.Results.datTree));
files = cellstr(string(p.Results.fileIDs(:)));
if isempty(files)
    files = unique(ctx.BWA_AQ2.fts1.fileID, 'stable');
end

BWA_use = bwaReduceByFiles(ctx.BWA_AQ2, files, ctx.filesTbl, ctx.eventsTbl);
memberlist = bwaComputeMemberlistForBWA(BWA_use, ctx.opts, ctx.DB);
dataType = bwaDatTreeToDataType(metric);
if isempty(dataType)
    metadata = {};
else
    metadata = AddLimitingRuleset({}, 'percentiles', dataType, p.Results.percentiles);
end
[figs, means_alln, means_celln] = checkBySlice(BWA_use.fts1, ctx.opts, metric, ctx.timeFrameRangeEpochs, memberlist, metadata, 'plotFigures', p.Results.plotFigures);
out = struct();
out.metric = metric;
out.fileIDs = files;
out.figures = {figs};
out.means_alln = {means_alln};
out.means_celln = {means_celln};
out.context = ctx;
end

function ctx = get_derived_context(result)
if isfield(result, 'derived') && isstruct(result.derived) && isfield(result.derived, 'context') ...
        && isstruct(result.derived.context) && isfield(result.derived.context, 'BWA_AQ2')
    ctx = result.derived.context;
    if ~isfield(ctx, 'filesTbl'), ctx.filesTbl = table(); end
    if ~isfield(ctx, 'eventsTbl'), ctx.eventsTbl = table(); end
    return
end
ctx = rebuild_derived_context(result);
end

function ctx = rebuild_derived_context(result)
if ~isfield(result, 'pathIn') || ~exist(char(string(result.pathIn)), 'dir')
    error('compute_local_file_metrics:MissingContext', ...
        'Slim result missing in-memory context; rerun with exportSlim=false or include result.pathIn.');
end
pathIn = char(string(result.pathIn));
loaded = loadClusterOutputs(pathIn, 'attachLegacyMeta', false, 'attachEventMeta', false);

if isfield(result, 'derived') && isfield(result.derived, 'context') && ...
        isstruct(result.derived.context) && isfield(result.derived.context, 'timeFrameRangeEpochs')
    epochSpec = result.derived.context.timeFrameRangeEpochs;
else
    if isfield(result, 'presetName') && ~isempty(result.presetName)
        epochSpec = loadEpochSpecFromConfig(pathIn, 'presetName', result.presetName);
    else
        epochSpec = loadEpochSpecFromConfig(pathIn);
    end
end

ctx = prepareLegacyAqua2Context(loaded.BWA_AQ2, ...
    'timeFrameRangeEpochs', epochSpec, ...
    'filesTbl', loaded.files, ...
    'eventsTbl', loaded.events, ...
    'eventCellMapTbl', loaded.event_cell_map);
end
