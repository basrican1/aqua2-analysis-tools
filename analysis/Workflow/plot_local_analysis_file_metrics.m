function out = plot_local_analysis_file_metrics(source, datTree, varargin)
% Plot file-by-file BAviolin plots for a single metric using saved data.
result = load_local_analysis_result(source);
caps = bwaDescribeResultCapabilities(result);
ctx = get_derived_context(result, caps);
p = inputParser;
p.StructExpand = false;
addOptional(p, 'datTree', datTree, @(x) ischar(x) || isstring(x));
addParameter(p, 'fileIDs', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'percentiles', [0 100], @isnumeric);
addParameter(p, 'plotArgs', {}, @iscell);
parse(p, datTree, varargin{:});

tTree = char(string(p.Results.datTree));
files = cellstr(string(p.Results.fileIDs(:)));
if isempty(files)
    files = unique(ctx.BWA_AQ2.fts1.fileID, 'stable');
end

BWA_use = bwaReduceByFiles(ctx.BWA_AQ2, files, ctx.filesTbl, ctx.eventsTbl);
memberlist = bwaComputeMemberlistForBWA(BWA_use, ctx.opts, ctx.DB);
dataType = bwaDatTreeToDataType(tTree);
if isempty(dataType)
    metadata = {};
else
    metadata = AddLimitingRuleset({}, 'percentiles', dataType, p.Results.percentiles);
end
[figs, means_alln, means_celln] = checkBySlice(BWA_use.fts1, ctx.opts, tTree, ctx.timeFrameRangeEpochs, memberlist, metadata);
out = struct('figures', {figs}, 'means_alln', {means_alln}, 'means_celln', {means_celln});
end

function ctx = get_derived_context(result, caps)
if isfield(result, 'derived') && isstruct(result.derived) && isfield(result.derived, 'context') ...
        && isstruct(result.derived.context) && isfield(result.derived.context, 'BWA_AQ2')
    ctx = result.derived.context;
    if ~isfield(ctx, 'filesTbl'), ctx.filesTbl = table(); end
    if ~isfield(ctx, 'eventsTbl'), ctx.eventsTbl = table(); end
    return
end
ctx = rebuild_derived_context(result, caps);
end

function ctx = rebuild_derived_context(result, caps)
if ~isfield(result, 'pathIn') || ~exist(char(string(result.pathIn)), 'dir')
    error('plot_local_analysis_file_metrics:MissingContext', ...
        'File-level legacy plotting requires a result with pathIn or saved derived context.');
end
pathIn = char(string(result.pathIn));
loaded = loadClusterOutputs(pathIn, 'attachLegacyMeta', false, 'attachEventMeta', false);
if ~isfield(loaded, 'hasLegacyPayload') || ~logical(loaded.hasLegacyPayload)
    error('plot_local_analysis_file_metrics:MissingLegacyPayload', ...
        ['File-level legacy plotting requires BWA_AQ2.mat because it rebuilds non-tabular legacy context. ' ...
         'Available payloads: %s. ' ...
         'Canonical summary plots remain available.'], char(string(caps.summary)));
end

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
