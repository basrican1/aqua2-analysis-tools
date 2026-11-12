function out = build_event_cell_adjacency_graph(source, varargin)
%BUILD_EVENT_CELL_ADJACENCY_GRAPH Diagnostic multicell event graph.
%   OUT = BUILD_EVENT_CELL_ADJACENCY_GRAPH(RESULT, 'eventUID', '<file>_E###')
%   builds one event-level graph where:
%     - nodes are participating cells
%     - node geometry is event footprint intersected with each cell mask
%     - edges connect partition geometries that touch or nearly touch
%     - timing comes from thresholded entry timing when available

p = inputParser;
p.StructExpand = false;
addRequired(p, 'source');
addParameter(p, 'eventUID', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'eventID', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'maxGapPx', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'showFigure', true, @islogical);
addParameter(p, 'visible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
parse(p, source, varargin{:});

[result, BWA_AQ2, filesTbl, sourcePath] = localResolveAdjacencySource(source);
eventUID = string(p.Results.eventUID);
if strlength(eventUID) == 0
    if isfinite(p.Results.eventID)
        eventUID = localEventUIDFromResultIndex(result, round(double(p.Results.eventID)));
    else
        error('build_event_cell_adjacency_graph:MissingEvent', 'Provide eventUID or eventID.');
    end
end

eventID = localEventIndexFromUIDAdj(BWA_AQ2, eventUID);
fileID = localEventFileIDAdj(BWA_AQ2, eventID);
eventPartition = localResolveEventPartitionAdj(result, eventUID);
if isempty(eventPartition)
    error('build_event_cell_adjacency_graph:MissingPartition', 'No event_partition rows found for %s.', char(eventUID));
end

[imageSz, spatialRes, frameRateS] = localEventImagingMetaAdj(BWA_AQ2, eventID);
footMask = localEventFootprintMaskAdj(BWA_AQ2, eventID, imageSz);
[masks, filesTbl] = localResolveCellMasksAdj(BWA_AQ2, filesTbl, sourcePath, eventID, fileID);
if isempty(masks)
    error('build_event_cell_adjacency_graph:MissingMasks', 'Cell masks unavailable for %s.', char(eventUID));
end

nodes = localBuildNodeTable(eventPartition, masks, fileID, filesTbl, footMask, spatialRes);
edges = localBuildEdgeTable(nodes, double(p.Results.maxGapPx), spatialRes);

G = graph(double(edges.source_idx), double(edges.target_idx), edges.min_gap_um, string(nodes.label));
out = struct();
out.event_uid = eventUID;
out.event_id = eventID;
out.file_id = fileID;
out.primary_cell_uid = localPrimaryCellUIDAdj(nodes);
out.nodes = nodes;
out.edges = edges;
out.graph = G;
out.footprint_mask = footMask;
out.image_size = imageSz;
out.spatial_res = spatialRes;
out.frame_rate_s = frameRateS;

if p.Results.showFigure
    out.figure = localPlotAdjacencyGraph(out, char(string(p.Results.visible)));
end
end

function [result, BWA_AQ2, filesTbl, sourcePath] = localResolveAdjacencySource(source)
result = struct();
filesTbl = table();
sourcePath = "";
if isstruct(source) && isfield(source, 'metrics')
    result = source;
    if isfield(source, 'pathIn')
        sourcePath = string(source.pathIn);
    elseif isfield(source, 'source') && isstruct(source.source) && isfield(source.source, 'pathIn')
        sourcePath = string(source.source.pathIn);
    end
    if isfield(source, 'canonical') && isstruct(source.canonical) && isfield(source.canonical, 'subset') && ...
            isstruct(source.canonical.subset) && isfield(source.canonical.subset, 'files') && istable(source.canonical.subset.files)
        filesTbl = source.canonical.subset.files;
    elseif isfield(source, 'derived') && isstruct(source.derived) && isfield(source.derived, 'context') && ...
            isstruct(source.derived.context) && isfield(source.derived.context, 'filesTbl') && istable(source.derived.context.filesTbl)
        filesTbl = source.derived.context.filesTbl;
    end
else
    error('build_event_cell_adjacency_graph:BadSource', 'Expected a local-analysis result struct.');
end
if strlength(sourcePath) == 0
    error('build_event_cell_adjacency_graph:MissingPath', 'Result must carry pathIn so BWA_AQ2 can be loaded.');
end
loaded = loadClusterOutputs(char(sourcePath), 'attachLegacyMeta', true, 'attachEventMeta', true);
BWA_AQ2 = loaded.BWA_AQ2;
if isempty(filesTbl) && isfield(loaded, 'files') && istable(loaded.files)
    filesTbl = loaded.files;
end
end

function eventUID = localEventUIDFromResultIndex(result, eventID)
eventUID = "";
if ~(isfield(result, 'metrics') && isstruct(result.metrics) && isfield(result.metrics, 'tables') && ...
        isstruct(result.metrics.tables) && isfield(result.metrics.tables, 'event_canonical') && ...
        istable(result.metrics.tables.event_canonical))
    return
end
T = result.metrics.tables.event_canonical;
if eventID >= 1 && eventID <= height(T) && ismember('event_uid', T.Properties.VariableNames)
    eventUID = string(T.event_uid(eventID));
end
end

function T = localResolveEventPartitionAdj(result, eventUID)
T = table();
if ~(isfield(result, 'metrics') && isstruct(result.metrics) && isfield(result.metrics, 'tables') && ...
        isstruct(result.metrics.tables) && isfield(result.metrics.tables, 'event_partition') && ...
        istable(result.metrics.tables.event_partition))
    return
end
P = result.metrics.tables.event_partition;
if ~ismember('event_uid', P.Properties.VariableNames)
    return
end
T = P(string(P.event_uid) == string(eventUID), :);
end

function eventID = localEventIndexFromUIDAdj(BWA_AQ2, eventUID)
eventID = NaN;
fileIDs = localEventFileIDsAdj(BWA_AQ2);
tok = regexp(char(eventUID), '^(.*)_E(\d+)$', 'tokens', 'once');
if isempty(tok)
    error('build_event_cell_adjacency_graph:BadEventUID', 'Expected eventUID like <file_id>_E<number>.');
end
fileID = string(tok{1});
localID = str2double(tok{2});
idx = find(fileIDs == fileID);
if numel(idx) < localID
    error('build_event_cell_adjacency_graph:EventUIDNotFound', 'eventUID %s not found in BWA_AQ2.', char(eventUID));
end
eventID = idx(localID);
end

function fileIDs = localEventFileIDsAdj(BWA_AQ2)
fileIDs = strings(0,1);
fts = BWA_AQ2.fts1;
if isfield(fts, 'basic') && isfield(fts.basic, 'fileID') && ~isempty(fts.basic.fileID)
    fileIDs = string(fts.basic.fileID);
elseif isfield(fts, 'fileID') && ~isempty(fts.fileID)
    fileIDs = string(fts.fileID);
end
end

function fileID = localEventFileIDAdj(BWA_AQ2, eventID)
fileIDs = localEventFileIDsAdj(BWA_AQ2);
fileID = string(fileIDs(eventID));
end

function [imageSz, spatialRes, frameRateS] = localEventImagingMetaAdj(BWA_AQ2, eventID)
opts = BWA_AQ2.opts;
imageSz = double(opts.sz(1,:));
spatialRes = double(opts.spatialRes(1));
frameRateS = NaN;
if isfield(opts, 'frameRate') && ~isempty(opts.frameRate)
    frameRateS = double(opts.frameRate(1));
end
fileIDs = localCanonCellstrAdj(BWA_AQ2.fts1.fileID);
thisID = fileIDs{eventID};
rowIDs = unique(fileIDs, 'stable');
row = find(strcmp(rowIDs, thisID), 1, 'first');
if size(opts.sz,1) >= row
    imageSz = double(opts.sz(row,:));
end
if numel(opts.spatialRes) >= row
    spatialRes = double(opts.spatialRes(row));
end
if isfield(opts, 'frameRate') && numel(opts.frameRate) >= row
    frameRateS = double(opts.frameRate(row));
end
if numel(imageSz) < 3
    imageSz(3) = 1;
end
end

function c = localCanonCellstrAdj(v)
if iscell(v)
    c = cell(size(v));
    for i = 1:numel(v)
        c{i} = char(string(v{i}));
    end
elseif isstring(v)
    c = cellstr(v(:));
elseif ischar(v)
    c = {v};
else
    c = cellstr(string(v(:)));
end
end

function footMask = localEventFootprintMaskAdj(BWA_AQ2, eventID, imageSz)
sz = double(imageSz(:)');
sx = round(sz(1)); sy = round(sz(2));
sz3 = 1;
if numel(sz) >= 3 && isfinite(sz(3)) && sz(3) > 0
    sz3 = round(sz(3));
end
fts = BWA_AQ2.fts1;
footMask = false(sx, sy);
t0 = round(double(fts.loc.t0(eventID)));
t1 = round(double(fts.loc.t1(eventID)));
for frameNow = t0:t1
    footMask = footMask | localExplorerFrameMaskAdj(fts.loc.xSpaTemp{eventID}, frameNow, sx, sy, sz3);
end
end

function mask = localExplorerFrameMaskAdj(xSpaTemp_ID, frameNum, sx, sy, sz3)
limitStart = sx * sy * sz3 * (frameNum - 1) + 1;
limitEnd = limitStart + sx * sy * sz3;
linMovie = double(xSpaTemp_ID(:));
linFrame = linMovie(linMovie >= limitStart & linMovie <= limitEnd) - limitStart + 1;
linFrame = linFrame(isfinite(linFrame) & linFrame >= 1 & linFrame <= sx * sy);
mask = false(sx, sy);
mask(linFrame) = true;
end

function [masks, filesTbl] = localResolveCellMasksAdj(BWA_AQ2, filesTbl, sourcePath, eventID, fileID)
masks = {};
allMasks = BWA_AQ2.fts1.region.cell.mask;
if istable(filesTbl) && ~isempty(filesTbl) && all(ismember({'file_id','n_cells'}, filesTbl.Properties.VariableNames))
    totalCells = sum(double(filesTbl.n_cells), 'omitnan');
    if isfinite(totalCells) && totalCells == numel(allMasks)
        masks = allMasks;
        return
    end
end
fileIDs = localEventFileIDsAdj(BWA_AQ2);
if string(fileID) == string(fileIDs(1))
    masks = allMasks;
    return
end
if strlength(sourcePath) == 0 || ~(istable(filesTbl) && all(ismember({'file_id','filename'}, filesTbl.Properties.VariableNames)))
    return
end
idx = find(string(filesTbl.file_id) == string(fileID), 1, 'first');
if isempty(idx)
    return
end
fname = char(string(filesTbl.filename(idx)));
d = dir(fullfile(char(sourcePath), '**', fname));
if isempty(d)
    return
end
s = load(fullfile(d(1).folder, d(1).name));
masks = localExtractCellMasksAdj(s);
end

function masks = localExtractCellMasksAdj(s)
masks = {};
candidates = {
    {'BWA_AQ2','fts1','region','cell','mask'}
    {'res','fts1','region','cell','mask'}
    {'datOut','fts1','region','cell','mask'}
    {'fts1','region','cell','mask'}
    };
for i = 1:numel(candidates)
    try
        v = s;
        for j = 1:numel(candidates{i})
            v = v.(candidates{i}{j});
        end
        if iscell(v)
            masks = v;
            return
        end
    catch
    end
end
end

function nodes = localBuildNodeTable(P, masks, fileID, filesTbl, footMask, spatialRes)
cellCol = localFirstPresentAdj(P.Properties.VariableNames, {'partition_cell_uid','cell_uid'});
n = height(P);
nodes = table('Size', [n 13], ...
    'VariableTypes', {'double','string','logical','double','double','double','double','double','double','double','double','cell','cell'}, ...
    'VariableNames', {'node_idx','cell_uid','is_primary_cell','entry_frame','delay_frames','order_rank','cell_index','centroid_x_px','centroid_y_px','centroid_x_um','centroid_y_um','partition_mask','boundary_xy_px'});
delayCol = localFirstPresentAdj(P.Properties.VariableNames, {'cell_entry_delay_from_primary_cov_thresh_frames','cell_entry_delay_from_primary_frames'});
orderVals = NaN(n,1);
delayVals = NaN(n,1);
if ~isempty(delayCol)
    delayVals = double(P.(delayCol));
    nonPrimary = find(~logical(P.is_primary_cell) & isfinite(delayVals));
    [~, ord] = sort(delayVals(nonPrimary), 'ascend');
    orderVals(nonPrimary(ord)) = 1:numel(ord);
end
for i = 1:n
    uid = string(P.(cellCol)(i));
    cellIdx = localCellIndexFromUIDAdj(uid, fileID, filesTbl);
    mask = false(size(footMask));
    if isfinite(cellIdx) && cellIdx >= 1 && cellIdx <= numel(masks)
        cm = logical(masks{cellIdx});
        mask = localMaskLikeAdj(cm, size(footMask));
        mask = mask & footMask;
    end
    bxy = zeros(0,2);
    cx = NaN; cy = NaN;
    if any(mask(:))
        B = bwboundaries(mask, 8, 'noholes');
        [~, ib] = max(cellfun(@(x) size(x,1), B));
        b = B{ib};
        bxy = [b(:,2), b(:,1)];
        [r,c] = find(mask);
        cx = mean(c); cy = mean(r);
    end
    nodes.node_idx(i) = i;
    nodes.cell_uid(i) = uid;
    nodes.is_primary_cell(i) = logical(P.is_primary_cell(i));
    if ismember('cell_entry_frame_cov_thresh', P.Properties.VariableNames)
        nodes.entry_frame(i) = double(P.cell_entry_frame_cov_thresh(i));
    elseif ismember('cell_first_overlap_frame', P.Properties.VariableNames)
        nodes.entry_frame(i) = double(P.cell_first_overlap_frame(i));
    end
    nodes.delay_frames(i) = delayVals(i);
    nodes.order_rank(i) = orderVals(i);
    nodes.cell_index(i) = cellIdx;
    nodes.centroid_x_px(i) = cx;
    nodes.centroid_y_px(i) = cy;
    nodes.centroid_x_um(i) = cx * spatialRes;
    nodes.centroid_y_um(i) = cy * spatialRes;
    nodes.partition_mask{i} = mask;
    nodes.boundary_xy_px{i} = bxy;
end
nodes.label = strings(height(nodes),1);
for i = 1:height(nodes)
    nodes.label(i) = localCellLabelAdj(nodes.cell_uid(i));
end
end

function edges = localBuildEdgeTable(nodes, maxGapPx, spatialRes)
edges = table('Size', [0 8], ...
    'VariableTypes', {'double','double','double','double','double','double','logical','logical'}, ...
    'VariableNames', {'source_idx','target_idx','source_delay_frames','target_delay_frames','min_gap_px','min_gap_um','touches','time_forward'});
n = height(nodes);
for i = 1:n
    Mi = nodes.partition_mask{i};
    if isempty(Mi) || ~any(Mi(:))
        continue
    end
    for j = (i+1):n
        Mj = nodes.partition_mask{j};
        if isempty(Mj) || ~any(Mj(:))
            continue
        end
        [gapPx, touches] = localPartitionGapAdj(Mi, Mj);
        if ~(isfinite(gapPx) && gapPx <= maxGapPx)
            continue
        end
        di = nodes.delay_frames(i);
        dj = nodes.delay_frames(j);
        timeForward = isfinite(di) && isfinite(dj) && di <= dj;
        sourceIdx = i;
        targetIdx = j;
        sourceDelay = di;
        targetDelay = dj;
        if isfinite(di) && isfinite(dj) && dj < di
            sourceIdx = j;
            targetIdx = i;
            sourceDelay = dj;
            targetDelay = di;
            timeForward = true;
        end
        edges(end+1,:) = {sourceIdx, targetIdx, sourceDelay, targetDelay, gapPx, gapPx * spatialRes, touches, timeForward}; %#ok<AGROW>
    end
end
end

function [gapPx, touches] = localPartitionGapAdj(M1, M2)
touches = any(imdilate(M1, ones(3,3)) & M2, 'all') || any(imdilate(M2, ones(3,3)) & M1, 'all');
if touches
    gapPx = 0;
    return
end
[r1,c1] = find(bwperim(M1));
[r2,c2] = find(bwperim(M2));
if isempty(r1) || isempty(r2)
    gapPx = inf;
    return
end
gapPx = min(pdist2([c1 r1], [c2 r2]), [], 'all');
end

function maskOut = localMaskLikeAdj(maskIn, imageSize2D)
maskOut = false(imageSize2D(1), imageSize2D(2));
sx = min(size(maskIn,1), imageSize2D(1));
sy = min(size(maskIn,2), imageSize2D(2));
maskOut(1:sx,1:sy) = maskIn(1:sx,1:sy);
end

function name = localFirstPresentAdj(names, candidates)
name = '';
for i = 1:numel(candidates)
    if any(strcmp(names, candidates{i}))
        name = candidates{i};
        return
    end
end
end

function idx = localCellIndexFromUIDAdj(cellUID, fileID, filesTbl)
idx = NaN;
uid = char(string(cellUID));
fid = char(string(fileID));
prefix = [fid '_C'];
if ~startsWith(uid, prefix)
    return
end
localIdx = str2double(extractAfter(string(uid), strlength(string(prefix))));
if ~isfinite(localIdx)
    return
end
idx = localIdx;
if istable(filesTbl) && ~isempty(filesTbl) && all(ismember({'file_id','n_cells'}, filesTbl.Properties.VariableNames))
    fileIDs = string(filesTbl.file_id);
    row = find(fileIDs == string(fileID), 1, 'first');
    if isempty(row)
        idx = NaN;
        return
    end
    priorCounts = double(filesTbl.n_cells(1:(row-1)));
    priorCounts = priorCounts(isfinite(priorCounts));
    idx = sum(priorCounts) + localIdx;
end
end

function uid = localPrimaryCellUIDAdj(nodes)
uid = "";
hit = find(nodes.is_primary_cell, 1, 'first');
if ~isempty(hit)
    uid = nodes.cell_uid(hit);
end
end

function lbl = localCellLabelAdj(cellUID)
lbl = string(cellUID);
tok = regexp(char(string(cellUID)), '_(C\d+)$', 'tokens', 'once');
if ~isempty(tok)
    lbl = string(tok{1});
end
end

function hFig = localPlotAdjacencyGraph(out, visible)
hFig = figure('Name', 'EventCellAdjacencyGraph', 'Color', 'w', 'NumberTitle', 'off', 'Visible', visible);
ax = axes('Parent', hFig);
hold(ax, 'on');
spatialRes = double(out.spatial_res);
sx = size(out.footprint_mask, 1);
sy = size(out.footprint_mask, 2);
bg = zeros(sx, sy, 3);
for c = 1:3
    plane = bg(:,:,c);
    plane(out.footprint_mask) = 0.72;
    bg(:,:,c) = plane;
end
image(ax, [0 sy * spatialRes], [0 sx * spatialRes], bg);
axis(ax, 'image');
set(ax, 'YDir', 'reverse');
xlabel(ax, 'x (um)');
ylabel(ax, 'y (um)');
delayVals = double(out.nodes.delay_frames);
validDelay = isfinite(delayVals) & ~out.nodes.is_primary_cell;
timingMap = parula(256);
delayScale = 1;
delayLabel = 'Thresholded Entry Delay (frames)';
if isfield(out, 'frame_rate_s') && isfinite(out.frame_rate_s) && out.frame_rate_s > 0
    delayScale = double(out.frame_rate_s);
    delayLabel = 'Thresholded Entry Delay (s)';
end
delayValsDisplay = delayVals * delayScale;
if any(validDelay)
    delayMin = 0;
    delayMax = max(delayValsDisplay(validDelay));
else
    delayMin = 0;
    delayMax = 1;
end
for i = 1:height(out.nodes)
    b = out.nodes.boundary_xy_px{i};
    if isempty(b)
        continue
    end
    color = [0.20 0.85 0.30];
    if isfinite(out.nodes.delay_frames(i))
        color = localDelayColorAdj(out.nodes.delay_frames(i) * delayScale, delayMin, delayMax, timingMap);
    end
    plot(ax, b(:,1) * spatialRes, b(:,2) * spatialRes, '-', 'Color', color, 'LineWidth', 1.5);
    if isfinite(out.nodes.centroid_x_px(i)) && isfinite(out.nodes.centroid_y_px(i))
        text(ax, out.nodes.centroid_x_um(i), out.nodes.centroid_y_um(i), char(out.nodes.label(i)), ...
            'Color', color, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
end
for i = 1:height(out.edges)
    s = out.edges.source_idx(i);
    t = out.edges.target_idx(i);
    if ~all(isfinite([out.nodes.centroid_x_um(s), out.nodes.centroid_y_um(s), out.nodes.centroid_x_um(t), out.nodes.centroid_y_um(t)]))
        continue
    end
    x1 = out.nodes.centroid_x_um(s);
    y1 = out.nodes.centroid_y_um(s);
    x2 = out.nodes.centroid_x_um(t);
    y2 = out.nodes.centroid_y_um(t);
    plot(ax, [x1, x2], [y1, y2], '-', 'Color', [0.95 0.55 0.10], 'LineWidth', 1.5);
    if logical(out.edges.time_forward(i))
        dx = x2 - x1;
        dy = y2 - y1;
        if hypot(dx, dy) > 0
            qx = x1 + 0.75 * dx;
            qy = y1 + 0.75 * dy;
            quiver(ax, qx, qy, 0.18 * dx, 0.18 * dy, 0, ...
                'Color', [0.95 0.55 0.10], 'LineWidth', 1.2, ...
                'MaxHeadSize', 2.0, 'AutoScale', 'off');
        end
    end
end
title(ax, sprintf('%s | %s', char(out.file_id), char(out.event_uid)));
if any(validDelay)
    colormap(ax, timingMap);
    caxis(ax, [delayMin delayMax + (delayMax == delayMin)]);
    cb = colorbar(ax);
    ylabel(cb, delayLabel);
    localAddTimingStripAdj(hFig, cb, out.nodes, delayMin, delayMax + (delayMax == delayMin), delayScale, timingMap);
end
hold(ax, 'off');
end

function color = localDelayColorAdj(v, vmin, vmax, cmap)
if ~(isfinite(v) && isfinite(vmin) && isfinite(vmax))
    color = [0.20 0.85 0.30];
    return
end
if vmax <= vmin
    idx = round(size(cmap,1) * 0.5);
else
    t = (v - vmin) / max(vmax - vmin, eps);
    idx = 1 + round(t * (size(cmap,1) - 1));
end
idx = min(max(idx, 1), size(cmap,1));
color = cmap(idx, :);
end

function localAddTimingStripAdj(hFig, cb, nodes, delayMin, delayMax, delayScale, cmap)
if ~(ishandle(hFig) && isgraphics(cb))
    return
end
delayVals = double(nodes.delay_frames) * delayScale;
valid = isfinite(delayVals);
if ~any(valid)
    return
end
cbPos = cb.Position;
stripWidth = 0.22;
gap = 0.015;
x0 = min(cbPos(1) + cbPos(3) + gap, 0.98 - stripWidth);
ax2 = axes('Parent', hFig, ...
    'Units', 'normalized', ...
    'Position', [x0 cbPos(2) stripWidth cbPos(4)], ...
    'Color', 'none', ...
    'XColor', 'none', ...
    'YColor', 'none', ...
    'Box', 'off', ...
    'YDir', 'normal', ...
    'HitTest', 'off');
hold(ax2, 'on');
xlim(ax2, [0 1.8]);
ylim(ax2, [delayMin delayMax]);
labelXBase = 0.55;
lastY = -Inf;
staggerIdx = 0;
dyThresh = max((delayMax - delayMin) * 0.04, eps);
xOffsets = [0.00 0.35 0.70 1.05];
validIdx = find(isfinite(delayVals));
[~, sortOrd] = sort(delayVals(validIdx), 'ascend');
validIdx = validIdx(sortOrd);
for kk = 1:numel(validIdx)
    i = validIdx(kk);
    y = delayVals(i);
    color = localDelayColorAdj(y, delayMin, delayMax, cmap);
    plot(ax2, [0.05 0.22], [y y], '-', 'Color', color, 'LineWidth', 1.5, 'Clipping', 'off');
    if abs(y - lastY) <= dyThresh
        staggerIdx = mod(staggerIdx + 1, numel(xOffsets));
    else
        staggerIdx = 0;
    end
    lastY = y;
    text(ax2, labelXBase + xOffsets(staggerIdx + 1), y, char(nodes.label(i)), ...
        'Color', color, ...
        'FontSize', 8, ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'middle', ...
        'Interpreter', 'none', ...
        'Clipping', 'off');
end
hold(ax2, 'off');
end
