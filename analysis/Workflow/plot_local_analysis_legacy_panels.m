function out = plot_local_analysis_legacy_panels(source, varargin)
%plot_local_analysis_legacy_panels Plot workflow-level legacy BAviolin figures.
%
% This is plot-only and does not recompute metrics. It consumes the
% precomputed payload from run_local_analysis_compute(...).
%
% Usage
%   out = plot_local_analysis_legacy_panels(result)
%   out = plot_local_analysis_legacy_panels(result, 'metricFamilies', {'shape'})
%   out = plot_local_analysis_legacy_panels(result, 'metricFamilies', {'shape','recruitment'})
%   out = plot_local_analysis_legacy_panels(result, 'metricList', {'basic_area'})
%   out = plot_local_analysis_legacy_panels(result, 'metricList', {'cell.mean_perimeter_um'})
%   out = plot_local_analysis_legacy_panels(result, 'metricList', {'cell.frac_propagative_events'})
%   out = plot_local_analysis_legacy_panels(result, 'metricList', {'cell.mean_frac_outside_origin_cell'})
%   out = plot_local_analysis_legacy_panels(result, 'metricList', {'cell.mean_soma_centroid_dist_um'})
%   out = plot_local_analysis_legacy_panels(result, 'metricList', {'cell.mean_neighbor_entry_delay_sec'})
%   list_local_legacy_metric_keys(result)  % one-command list of available metric names
%
% Input
%   source : result struct from run_local_analysis_compute(...) OR a path/struct
%            accepted by load_local_analysis_result().
%
% Name-value options
%   'metricFamilies' (char|string|cell, default {})
%       Biological-family selector. Supported values are:
%         - 'activity'
%         - 'shape'
%         - 'kinetics'
%         - 'intensity'
%         - 'propagation'
%         - 'containment'
%         - 'recruitment'
%         - 'soma'
%       Slice-level synchrony metrics are not rendered by this wrapper.
%
%   'metricList' (cell|string, default {})
%       Unified selector across standard, extras, and cell-metric panels.
%       Empty => all enabled panels.
%       Non-empty => only matching panels are plotted (single selector for all types).
%       If both 'metricFamilies' and 'metricList' are supplied, their union is used.
%
%   'includeExtras' (logical, default true)
%       Include extra derived BAviolin-style panels (footprint, propagation,
%       centroid) when available.
%
%   'includeCellMetrics' (logical, default true)
%       Include cell-level derived metric BAviolin panels.
%
%   'style' (struct, default struct())
%       Plot style override struct merged with default style.
%
%   'statsConfig' (struct, default struct())
%       Statistics configuration passed into BAviolin subplot stats helpers.
%
%   'doByCellOnly' (logical, default false)
%       If true, emphasize/plot cell-normalized branch only where supported.
%
% Output
%   out struct with fields:
%     - out.baviolin      : figures for standard BAviolin metrics
%     - out.extras        : figures for footprint/propagation/centroid panels
%     - out.cell_metrics  : figures for cell-level derived metrics
%     - out.skipped       : per-metric errors if a panel failed
%
% Notes
%   - Epoch definitions are resolved from result.derived.context or epochMeta.
%   - The inset slice spider reuses result.spiderCombineEpochNames /
%     result.spiderCombinedEpochLabel when available, so legacy BAviolin
%     figures can follow the same combined post-manipulation framework as the main
%     spider workflow.
%   - This wrapper calls plotBAviolinPanels() for each available metric payload.
%   - Common cell metric keys include:
%       cell.mean_perimeter_um
%       cell.mean_circularity
%       cell.frac_propagative_events
%       cell.mean_propagation_distance_um
%       cell.mean_frac_outside_origin_cell
%       cell.frac_escape_events_gt_threshold
%       cell.mean_soma_centroid_dist_um
%       cell.mean_soma_start_dist_um
%       cell.mean_soma_end_dist_um
%       cell.frac_events_with_neighbor_recruitment
%       cell.mean_neighbor_entry_delay_sec
%       cell.mean_neighbor_peak_overlap_frac_of_cell
%       cell.mean_neighbor_partitions_per_event
userSetIncludeExtras = hasNameValueArg(varargin, 'includeExtras');
userSetIncludeCellMetrics = hasNameValueArg(varargin, 'includeCellMetrics');

result = load_local_analysis_result(source);
caps = bwaDescribeResultCapabilities(result);
if ~isfield(result, 'derived') || ~isstruct(result.derived) || ...
        ~isfield(result.derived, 'standard') || ~isstruct(result.derived.standard)
    error('plot_local_analysis_legacy_panels:MissingStandard', ...
        ['Legacy-panel plotting requires result.derived.standard. ' ...
         'This payload is only available from a legacy-enriched compute run. ' ...
         'Available payloads: %s. ' ...
         'Canonical summary plots remain available.'], char(string(caps.summary)));
end

p = inputParser;
p.StructExpand = false;
addParameter(p, 'metricList', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'metricFamilies', {}, @(x) iscell(x) || isstring(x) || ischar(x));
addParameter(p, 'extraMetricList', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'includeExtras', true, @islogical); % footprint/propagation/centroid
addParameter(p, 'includeCellMetrics', true, @islogical);
addParameter(p, 'style', struct(), @isstruct);
addParameter(p, 'statsConfig', struct(), @isstruct);
addParameter(p, 'doByCellOnly', false, @islogical);
parse(p, varargin{:});

includeExtras = p.Results.includeExtras;
includeCellMetrics = p.Results.includeCellMetrics;
extraMetricList = cellstr(string(p.Results.extraMetricList(:)))';
reqAll = cellstr(string(p.Results.metricList(:)))';

spiderCombineEpochNames = {};
if isfield(result, 'spiderCombineEpochNames') && ~isempty(result.spiderCombineEpochNames)
    spiderCombineEpochNames = cellstr(string(result.spiderCombineEpochNames(:)))';
end
spiderCombinedEpochLabel = 'CNO_post';
if isfield(result, 'spiderCombinedEpochLabel') && ~isempty(result.spiderCombinedEpochLabel)
    spiderCombinedEpochLabel = char(string(result.spiderCombinedEpochLabel));
end

stdOut = result.derived.standard;
if iscell(p.Results.metricFamilies)
    familyReq = cellstr(string(p.Results.metricFamilies(:)))';
elseif isstring(p.Results.metricFamilies)
    familyReq = cellstr(p.Results.metricFamilies(:))';
elseif isempty(p.Results.metricFamilies)
    familyReq = {};
else
    familyReq = cellstr(string(p.Results.metricFamilies));
end
if ~isempty(familyReq)
    reqAll = unique([reqAll, expandLegacyMetricFamilies(result, stdOut, familyReq)], 'stable');
    if isempty(p.Results.metricList)
        reqAll = suppressLegacyFamilyDuplicates(reqAll);
    end
end
tfre = getEpochSpec(result);

out = struct();
out.baviolin = struct();
out.extras = struct();
out.cell_metrics = struct();
out.skipped = struct();

allNames = {};
if isfield(stdOut, 'baviolin') && isstruct(stdOut.baviolin)
    allNames = fieldnames(stdOut.baviolin);
end
names = allNames;
keepExtras = {};
keepCellMetrics = {};

if ~isempty(reqAll)
    names = resolveMetricKeys(allNames, reqAll);
    keepExtras = resolveExtraKeys(availableExtraKeys(stdOut), reqAll);
    keepCellMetrics = resolveCellMetricKeys(result, reqAll);

    if ~userSetIncludeExtras
        includeExtras = ~isempty(keepExtras);
    end
    if ~userSetIncludeCellMetrics
        includeCellMetrics = ~isempty(keepCellMetrics);
    end

    if isempty(names) && isempty(keepExtras) && isempty(keepCellMetrics)
        warning('plot_local_analysis_legacy_panels:NoMetricMatch', ...
            'No metricList entries matched standard/extras/cell metrics. Requested: %s', strjoin(reqAll, ', '));
    end

    if isempty(extraMetricList)
        extraMetricList = keepExtras;
    end
end

if isfield(stdOut, 'baviolin') && isstruct(stdOut.baviolin)

    for i = 1:numel(names)
        key = names{i};
        comp = stdOut.baviolin.(key);
        try
            if (~isfield(comp, 'filesTbl') || ~istable(comp.filesTbl) || isempty(comp.filesTbl)) && ...
                    isfield(result, 'canonical') && isstruct(result.canonical) && ...
                    isfield(result.canonical, 'loaded') && isstruct(result.canonical.loaded) && ...
                    isfield(result.canonical.loaded, 'files') && istable(result.canonical.loaded.files)
                comp.filesTbl = result.canonical.loaded.files;
            end
            memberlist = inferMemberlist(comp);
            [h_fig, means_alln, means_celln] = plotBAviolinPanels(comp, tfre, memberlist, ...
                'style', p.Results.style, ...
                'statsConfig', p.Results.statsConfig, ...
                'doByCellOnly', p.Results.doByCellOnly, ...
                'spiderCombineEpochNames', spiderCombineEpochNames, ...
                'spiderCombinedEpochLabel', spiderCombinedEpochLabel);
            out.baviolin.(key) = struct('figure', h_fig, 'means_alln', means_alln, 'means_celln', means_celln);
        catch ME
            out.skipped.(key) = ME.message;
            warning('plot_local_analysis_legacy_panels:PlotFailed', ...
                'Failed plotting %s: %s', key, ME.message);
        end
    end
else
    warning('plot_local_analysis_legacy_panels:NoBaviolin', ...
        'No standard.baviolin payload found.');
end

if includeExtras
    keepE = {};
    if ~isempty(extraMetricList)
        keepE = resolveExtraKeys(availableExtraKeys(stdOut), extraMetricList);
        if isempty(keepE)
            warning('plot_local_analysis_legacy_panels:NoExtraMetricMatch', ...
                'No extra metrics matched request: %s', strjoin(cellstr(string(extraMetricList(:)))', ', '));
        end
    end
    out.extras = plotExtraDerivedPanels(stdOut, tfre, p.Results.style, p.Results.statsConfig, p.Results.doByCellOnly, keepE, spiderCombineEpochNames, spiderCombinedEpochLabel);
end
if includeCellMetrics
    out.cell_metrics = plotCellMetricPanels(result, tfre, p.Results.style, p.Results.statsConfig, keepCellMetrics, spiderCombineEpochNames, spiderCombinedEpochLabel);
end
end

function tfre = getEpochSpec(result)
if isfield(result, 'derived') && isfield(result.derived, 'context') && ...
        isstruct(result.derived.context) && isfield(result.derived.context, 'timeFrameRangeEpochs')
    tfre = result.derived.context.timeFrameRangeEpochs;
    return
end
if isfield(result, 'derived') && isfield(result.derived, 'epochMeta') && ...
        isstruct(result.derived.epochMeta) && isfield(result.derived.epochMeta, 'sourceFile')
    rootDir = '';
    if isfield(result.derived.epochMeta, 'rootDir')
        rootDir = char(string(result.derived.epochMeta.rootDir));
    end
    sourceFile = char(string(result.derived.epochMeta.sourceFile));
    if ~isempty(rootDir) && isfolder(rootDir)
        if isfield(result, 'presetName') && ~isempty(result.presetName)
            tfre = loadEpochSpecFromConfig(rootDir, 'presetName', result.presetName);
        else
            tfre = loadEpochSpecFromConfig(rootDir);
        end
        return
    end
    if ~isempty(sourceFile) && isfolder(sourceFile)
        if isfield(result, 'presetName') && ~isempty(result.presetName)
            tfre = loadEpochSpecFromConfig(sourceFile, 'presetName', result.presetName);
        else
            tfre = loadEpochSpecFromConfig(sourceFile);
        end
        return
    end
end
error('plot_local_analysis_legacy_panels:MissingEpochSpec', ...
    'Could not resolve epoch spec from result.derived.context or epochMeta.');
end

function memberlist = inferMemberlist(comp)
if isfield(comp, 'memberlist') && ~isempty(comp.memberlist)
    memberlist = double(comp.memberlist(:)');
    return
end
n = 0;
if isfield(comp, 'dffAll')
    n = numel(comp.dffAll);
end
memberlist = 1:max(1,n);
end

function out = plotExtraDerivedPanels(stdOut, tfre, style, statsConfig, doByCellOnly, keepExtras, spiderCombineEpochNames, spiderCombinedEpochLabel)
out = struct();
if nargin < 6 || isempty(keepExtras), keepExtras = {}; end
if nargin < 7, spiderCombineEpochNames = {}; end
if nargin < 8, spiderCombinedEpochLabel = 'CNO_post'; end

if (isempty(keepExtras) || any(strcmp('footprint', keepExtras))) && isfield(stdOut, 'footprint') && isstruct(stdOut.footprint) && isfield(stdOut.footprint, 'dffAll')
    try
        memberlist = inferMemberlist(stdOut.footprint);
        [h_fig, means_alln, means_celln] = plotBAviolinPanels(stdOut.footprint, tfre, memberlist, ...
            'style', style, 'statsConfig', statsConfig, 'doByCellOnly', doByCellOnly, ...
            'spiderCombineEpochNames', spiderCombineEpochNames, 'spiderCombinedEpochLabel', spiderCombinedEpochLabel);
        out.footprint = struct('figure', h_fig, 'means_alln', means_alln, 'means_celln', means_celln);
    catch ME
        out.footprint_error = ME.message;
    end
end

if isfield(stdOut, 'propagation') && isstruct(stdOut.propagation)
    if (isempty(keepExtras) || any(strcmp('propagation_distance', keepExtras))) && isfield(stdOut.propagation, 'distance') && isstruct(stdOut.propagation.distance)
        try
            memberlist = inferMemberlist(stdOut.propagation.distance);
            [h_fig, means_alln, means_celln] = plotBAviolinPanels(stdOut.propagation.distance, tfre, memberlist, ...
                'style', style, 'statsConfig', statsConfig, 'doByCellOnly', doByCellOnly, ...
            'spiderCombineEpochNames', spiderCombineEpochNames, 'spiderCombinedEpochLabel', spiderCombinedEpochLabel);
            out.propagation_distance = struct('figure', h_fig, 'means_alln', means_alln, 'means_celln', means_celln);
        catch ME
            out.propagation_distance_error = ME.message;
        end
    end
    if (isempty(keepExtras) || any(strcmp('propagation_angle', keepExtras))) && isfield(stdOut.propagation, 'angle') && isstruct(stdOut.propagation.angle)
        try
            memberlist = inferMemberlist(stdOut.propagation.angle);
            [h_fig, means_alln, means_celln] = plotBAviolinPanels(stdOut.propagation.angle, tfre, memberlist, ...
                'style', style, 'statsConfig', statsConfig, 'doByCellOnly', doByCellOnly, ...
            'spiderCombineEpochNames', spiderCombineEpochNames, 'spiderCombinedEpochLabel', spiderCombinedEpochLabel);
            out.propagation_angle = struct('figure', h_fig, 'means_alln', means_alln, 'means_celln', means_celln);
        catch ME
            out.propagation_angle_error = ME.message;
        end
    end
end

if (isempty(keepExtras) || any(strcmp('centroid', keepExtras))) && isfield(stdOut, 'centroid') && isstruct(stdOut.centroid) && ...
        isfield(stdOut.centroid, 'distComp') && isstruct(stdOut.centroid.distComp)
    try
        memberlist = inferMemberlist(stdOut.centroid.distComp);
        [h_fig, means_alln, means_celln] = plotBAviolinPanels(stdOut.centroid.distComp, tfre, memberlist, ...
            'style', style, 'statsConfig', statsConfig, 'doByCellOnly', doByCellOnly, ...
            'spiderCombineEpochNames', spiderCombineEpochNames, 'spiderCombinedEpochLabel', spiderCombinedEpochLabel);
        out.centroid = struct('figure', h_fig, 'means_alln', means_alln, 'means_celln', means_celln);
    catch ME
        out.centroid_error = ME.message;
    end
end
end

function out = plotCellMetricPanels(result, tfre, style, statsConfig, keepCellMetrics, spiderCombineEpochNames, spiderCombinedEpochLabel)
out = struct();
if nargin < 5, keepCellMetrics = {}; end
if nargin < 6, spiderCombineEpochNames = {}; end
if nargin < 7, spiderCombinedEpochLabel = 'CNO_post'; end
if ~isfield(result, 'metrics') || ~isstruct(result.metrics) || ...
        ~isfield(result.metrics, 'registry') || ~istable(result.metrics.registry) || ...
        ~isfield(result.metrics, 'tables') || ~isstruct(result.metrics.tables)
    return
end

reg = result.metrics.registry;
filesMeta = table();
if isfield(result, 'canonical') && isstruct(result.canonical)
    if isfield(result.canonical, 'loaded') && isstruct(result.canonical.loaded) && ...
            isfield(result.canonical.loaded, 'files') && istable(result.canonical.loaded.files)
        filesMeta = result.canonical.loaded.files;
    elseif isfield(result.canonical, 'subset') && isstruct(result.canonical.subset) && ...
            isfield(result.canonical.subset, 'files') && istable(result.canonical.subset.files)
        filesMeta = result.canonical.subset.files;
    end
end
for i = 1:height(reg)
    if ~strcmpi(string(reg.unit_level(i)), "cell")
        continue
    end
    metricID = char(string(reg.metric_id(i)));
    if ~isempty(keepCellMetrics)
        mNorm = normalizeKey(metricID);
        kk = cellfun(@normalizeKey, keepCellMetrics, 'UniformOutput', false);
        hit = false;
        for hk = 1:numel(kk)
            if contains(mNorm, kk{hk}) || contains(kk{hk}, mNorm)
                hit = true;
                break
            end
        end
        if ~hit
            continue
        end
    end
    tname = char(string(reg.source_table(i)));
    rawCol = char(string(reg.value_column_raw(i)));
    normCol = char(string(reg.value_column_norm(i)));
    if ~isfield(result.metrics.tables, tname)
        continue
    end
    T = result.metrics.tables.(tname);
    if ~istable(T) || isempty(T)
        continue
    end

    try
        comp = buildCellMetricComp(T, tfre, metricID, rawCol, normCol, filesMeta);
        if isempty(comp)
            continue
        end
        [h_fig, means_alln, means_celln] = plotBAviolinPanels(comp, tfre, 1, ...
            'style', style, 'statsConfig', statsConfig, ...
            'doByCellOnly', true, 'showSliceSpider', true, ...
            'spiderCombineEpochNames', spiderCombineEpochNames, ...
            'spiderCombinedEpochLabel', spiderCombinedEpochLabel);
        key = matlab.lang.makeValidName(strrep(metricID, '.', '_'));
        out.(key) = struct('figure', h_fig, 'means_alln', means_alln, 'means_celln', means_celln);
    catch ME
        key = matlab.lang.makeValidName(strrep(metricID, '.', '_'));
        out.([key '_error']) = ME.message;
    end
end
end

function comp = buildCellMetricComp(T, tfre, metricID, rawCol, normCol, filesMeta)
comp = struct();
need = {'epoch_idx','epoch_name'};
if ~all(ismember(need, T.Properties.VariableNames))
    return
end
if ~ismember(rawCol, T.Properties.VariableNames)
    return
end
if ~ismember('primary_cell_uid', T.Properties.VariableNames)
    if ismember('cell_uid', T.Properties.VariableNames)
        T.primary_cell_uid = string(T.cell_uid);
    else
        return
    end
end

raw = double(T.(rawCol));
ep = double(T.epoch_idx);
cells = string(T.primary_cell_uid);
cellNum = double(categorical(cells));
keepR = isfinite(raw) & isfinite(ep) & ep > 0 & strlength(cells) > 0;
if ~any(keepR)
    comp = [];
    return
end

if ~isempty(normCol) && ismember(normCol, T.Properties.VariableNames)
    normv = double(T.(normCol));
else
    normv = computeBaselineLog2PerCell(raw, ep, cells);
end
keepN = isfinite(normv) & isfinite(ep) & ep > 0 & strlength(cells) > 0;

comp.datTxt = metricID;
comp.avgs_byCell = raw(keepR);
comp.tCat_byCell = categorical(string(T.epoch_name(keepR)), string(tfre.epochNames), 'Ordinal', true);
comp.cells_byCell = cellNum(keepR);

comp.avgsn_byCell = normv(keepN);
comp.tCatn_byCell = categorical(string(T.epoch_name(keepN)), string(tfre.epochNames), 'Ordinal', true);
comp.cellsn_byCell = cellNum(keepN);

% Keep fold-scatter vectors aligned with normalized rows.
comp.outDffs_c = raw(keepN);
comp.outMember_c = cellNum(keepN);
comp.outEpochIndex_c = ep(keepN);

% Spider-compatible fields (table-native cell metric rows treated as per-slice observations).
comp.normByEvent = normv(keepN);
comp.epMembIdx = ep(keepN);
if ismember('file_id', T.Properties.VariableNames)
    comp.fileID = string(T.file_id(keepN));
else
    comp.fileID = repmat("", nnz(keepN), 1);
end
comp.filesTbl = buildFilesTblForLabels(T, keepN, filesMeta);
end

function ft = buildFilesTblForLabels(T, keepN, filesMeta)
want = {'file_id','exp_date','slice_num'};
ft = table();
if ismember('file_id', T.Properties.VariableNames)
    fkeep = unique(string(T.file_id(keepN)));
else
    return
end

if all(ismember(want, T.Properties.VariableNames))
    ft = unique(T(:, want), 'rows');
else
    ft = unique(T(:, intersect({'file_id'}, T.Properties.VariableNames)), 'rows');
end

if ~isempty(filesMeta) && istable(filesMeta) && all(ismember(want, filesMeta.Properties.VariableNames))
    fm = filesMeta(:, want);
    fm.file_id = string(fm.file_id);
    fm = unique(fm, 'rows');
    if ~ismember('file_id', ft.Properties.VariableNames)
        ft = table(fkeep, 'VariableNames', {'file_id'});
    end
    ft.file_id = string(ft.file_id);
    ft = outerjoin(ft, fm, 'Keys', 'file_id', 'Type', 'left', 'MergeKeys', true);
end

if ~ismember('exp_date', ft.Properties.VariableNames)
    ft.exp_date = repmat("", height(ft), 1);
end
if ~ismember('slice_num', ft.Properties.VariableNames)
    ft.slice_num = repmat("", height(ft), 1);
end
ft.file_id = string(ft.file_id);
ft = ft(ismember(ft.file_id, fkeep), :);
ft = unique(ft(:, want), 'rows');
end

function y = computeBaselineLog2PerCell(raw, ep, cells)
y = NaN(size(raw));
u = unique(cells, 'stable');
for i = 1:numel(u)
    idx = cells == u(i);
    b = raw(idx & ep == 1);
    b = b(isfinite(b) & b > 0);
    if isempty(b)
        continue
    end
    b0 = mean(b, 'omitnan');
    x = raw(idx);
    ok = isfinite(x) & x > 0;
    z = NaN(size(x));
    z(ok) = log2(x(ok) ./ b0);
    y(idx) = z;
end
end


function namesOut = resolveMetricKeys(available, requested)
available = cellstr(string(available(:)))';
req = cellstr(string(requested(:)))';
if isempty(req)
    namesOut = available;
    return
end

map = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(available)
    k = available{i};
    addAlias(map, k, k);
    addAlias(map, normalizeKey(k), k);
    addAlias(map, normalizeKey(strrep(k, '_', '.')), k);
    addAlias(map, normalizeKey(strrep(k, '.', '')), k);
end

% Common legacy/canonical aliases.
aliasPairs = {
    'area', 'basic_area';
    'basic.area', 'basic_area';
    'duration', 'curve_duration';
    'curve.duration', 'curve_duration';
    'rise19', 'curve_rise19';
    'curve.rise19', 'curve_rise19';
    'fall91', 'curve_fall91';
    'curve.fall91', 'curve_fall91';
    'decaytau', 'curve_decayTau';
    'curve.decaytau', 'curve_decayTau';
    'dffauc', 'curve_dffAUC';
    'curve.dffauc', 'curve_dffAUC';
    'dffmax2', 'curve_dffMax2';
    'curve.dffmax2', 'curve_dffMax2';
    'maxpropspeed', 'propagation_maxPropSpeed';
    'propagation.maxpropspeed', 'propagation_maxPropSpeed';
};
for i = 1:size(aliasPairs,1)
    src = normalizeKey(aliasPairs{i,1});
    dst = aliasPairs{i,2};
    if any(strcmpi(dst, available))
        addAlias(map, src, dst);
    end
end

picked = {};
for i = 1:numel(req)
    r = req{i};
    cand = normalizeKey(r);
    if isKey(map, cand)
        key = map(cand);
        if ~any(strcmp(key, picked))
            picked{end+1} = key; %#ok<AGROW>
        end
    end
end


if isempty(picked)
    reqN = cellfun(@normalizeKey, req, 'UniformOutput', false);
    for i = 1:numel(available)
        a = available{i};
        an = normalizeKey(a);
        for j = 1:numel(reqN)
            if contains(an, reqN{j}) || contains(reqN{j}, an)
                if ~any(strcmp(a, picked))
                    picked{end+1} = a; %#ok<AGROW>
                end
                break
            end
        end
    end
end

if isempty(picked)
    namesOut = {};
else
    namesOut = picked;
end
end

function out = normalizeKey(in)
out = lower(char(string(in)));
out = strrep(out, ' ', '');
out = strrep(out, '-', '');
out = strrep(out, '_', '');
out = strrep(out, '.', '');
out = strrep(out, '/', '');
end

function addAlias(map, alias, canonical)
alias = char(string(alias));
canonical = char(string(canonical));
if isempty(alias)
    return
end
if ~isKey(map, alias)
    map(alias) = canonical;
end
end


function tf = hasNameValueArg(args, name)
tf = false;
if isempty(args)
    return
end
if ~iscell(args)
    return
end
nm = lower(char(string(name)));
for i = 1:2:numel(args)
    if ischar(args{i}) || isstring(args{i})
        if strcmpi(char(string(args{i})), nm)
            tf = true;
            return
        end
    end
end
end


function out = filterExtraPanels(out, requested)
if ~isstruct(out) || isempty(fieldnames(out))
    return
end
allk = fieldnames(out);
keep = resolveExtraKeys(allk, requested);
for i = 1:numel(allk)
    if ~any(strcmp(allk{i}, keep))
        out = rmfield(out, allk{i});
    end
end
end

function keys = resolveExtraKeys(available, requested)
av = cellstr(string(available(:)))';
req = cellstr(string(requested(:)))';
if isempty(req)
    keys = av;
    return
end

map = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(av)
    addAlias(map, av{i}, av{i});
    addAlias(map, normalizeKey(av{i}), av{i});
    addAlias(map, normalizeKey(strrep(av{i}, '_', '.')), av{i});
end

% explicit aliases
explicit = {
    'centroid', 'centroid';
    'centroiddistcomp', 'centroid';
    'centroid.distcomp', 'centroid';
    'distcomp', 'centroid';
    'footprint', 'footprint';
    'propagation', 'propagation_distance';
    'propagationdistance', 'propagation_distance';
    'propagation.distance', 'propagation_distance';
    'propdistance', 'propagation_distance';
    'angle', 'propagation_angle';
    'propagationangle', 'propagation_angle';
    'propagation.angle', 'propagation_angle';
    'propangle', 'propagation_angle';
};
for i = 1:size(explicit,1)
    src = normalizeKey(explicit{i,1});
    dst = explicit{i,2};
    if any(strcmp(dst, av))
        addAlias(map, src, dst);
    end
end

picked = {};
for i = 1:numel(req)
    r = normalizeKey(req{i});
    if isKey(map, r)
        k = map(r);
        if ~any(strcmp(k, picked)), picked{end+1}=k; end %#ok<AGROW>
    end
end

if isempty(picked)
    % contains fallback
    reqN = cellfun(@normalizeKey, req, 'UniformOutput', false);
    for i = 1:numel(av)
        an = normalizeKey(av{i});
        for j = 1:numel(reqN)
            if contains(an, reqN{j}) || contains(reqN{j}, an)
                if ~any(strcmp(av{i}, picked)), picked{end+1}=av{i}; end %#ok<AGROW>
                break
            end
        end
    end
end

keys = picked;
end


function keys = availableExtraKeys(stdOut)
keys = {};
if isfield(stdOut, 'footprint') && isstruct(stdOut.footprint) && isfield(stdOut.footprint, 'dffAll')
    keys{end+1} = 'footprint'; %#ok<AGROW>
end
if isfield(stdOut, 'propagation') && isstruct(stdOut.propagation)
    if isfield(stdOut.propagation, 'distance') && isstruct(stdOut.propagation.distance)
        keys{end+1} = 'propagation_distance'; %#ok<AGROW>
    end
    if isfield(stdOut.propagation, 'angle') && isstruct(stdOut.propagation.angle)
        keys{end+1} = 'propagation_angle'; %#ok<AGROW>
    end
end
if isfield(stdOut, 'centroid') && isstruct(stdOut.centroid) && ...
        isfield(stdOut.centroid, 'distComp') && isstruct(stdOut.centroid.distComp)
    keys{end+1} = 'centroid'; %#ok<AGROW>
end
end


function keys = resolveCellMetricKeys(result, requested)
keys = {};
if ~isfield(result, 'metrics') || ~isstruct(result.metrics) || ...
        ~isfield(result.metrics, 'registry') || ~istable(result.metrics.registry)
    return
end
R = result.metrics.registry;
if ~all(ismember({'metric_id','unit_level'}, R.Properties.VariableNames))
    return
end
ids = cellstr(string(R.metric_id(strcmpi(string(R.unit_level), "cell"))))';
if isempty(ids)
    return
end
req = cellstr(string(requested(:)))';
if isempty(req)
    keys = ids;
    return
end
picked = {};
reqN = cellfun(@normalizeKey, req, 'UniformOutput', false);
for i = 1:numel(ids)
    in = normalizeKey(ids{i});
    for j = 1:numel(reqN)
        if contains(in, reqN{j}) || contains(reqN{j}, in)
            if ~any(strcmp(ids{i}, picked)), picked{end+1}=ids{i}; end %#ok<AGROW>
            break
        end
    end
end
keys = picked;
end

function req = expandLegacyMetricFamilies(result, stdOut, families)
if iscell(families)
    families = lower(string(families(:)))';
else
    families = lower(string(families));
end
families = families(strlength(families) > 0);
if isempty(families)
    req = {};
    return
end

allStd = {};
if isfield(stdOut, 'baviolin') && isstruct(stdOut.baviolin)
    allStd = fieldnames(stdOut.baviolin)';
end
allExtra = availableExtraKeys(stdOut);
allCell = {};
if isfield(result, 'metrics') && isstruct(result.metrics) && ...
        isfield(result.metrics, 'registry') && istable(result.metrics.registry) && ...
        all(ismember({'metric_id','unit_level'}, result.metrics.registry.Properties.VariableNames))
    R = result.metrics.registry;
    m = strcmpi(string(R.unit_level), "cell");
    allCell = cellstr(unique(string(R.metric_id(m)), 'stable'))';
    allCell = filterPreferredLegacyGroupCellKeys(allCell);
end

req = {};
allReq = [allStd, allExtra, allCell];
if any(families == "all")
    req = unique(allReq, 'stable');
    return
end
for i = 1:numel(allReq)
    family = bwaResolveMetricFamily(allReq{i}, result.metrics.registry);
    if any(family == families)
        req = unique([req, allReq(i)], 'stable');
    end
end
end

function req = suppressLegacyFamilyDuplicates(req)
req = cellstr(string(req(:)))';
if isempty(req)
    return
end
conceptsFromStandard = strings(0,1);
for i = 1:numel(req)
    key = string(req{i});
    if startsWith(key, "cell.")
        continue
    end
    conceptsFromStandard(end+1,1) = bwaMetricConceptKey(key); %#ok<AGROW>
end
conceptsFromStandard = unique(conceptsFromStandard, 'stable');
keep = true(1, numel(req));
for i = 1:numel(req)
    key = string(req{i});
    if ~startsWith(key, "cell.")
        continue
    end
    if any(bwaMetricConceptKey(key) == conceptsFromStandard)
        keep(i) = false;
    end
end
req = req(keep);
end

function keys = filterPreferredLegacyGroupCellKeys(keys)
if isempty(keys)
    return
end
keys = string(keys(:));
pairs = [ ...
    "cell.mean_soma_centroid_dist_px", "cell.mean_soma_centroid_dist_um"; ...
    "cell.mean_soma_start_dist_px", "cell.mean_soma_start_dist_um"; ...
    "cell.mean_soma_end_dist_px", "cell.mean_soma_end_dist_um"; ...
    "cell.mean_soma_net_distance_change_px", "cell.mean_soma_net_distance_change_um"; ...
    "cell.mean_perimeter_px", "cell.mean_perimeter_um"; ...
    "cell.mean_propagation_distance_px", "cell.mean_propagation_distance_um"; ...
    "cell.mean_neighbor_entry_delay_frames", "cell.mean_neighbor_entry_delay_sec"];
keep = true(size(keys));
for i = 1:size(pairs,1)
    raw = pairs(i,1); cal = pairs(i,2);
    if any(keys == cal)
        keep(keys == raw) = false;
    end
end
keys = cellstr(keys(keep))';
end
