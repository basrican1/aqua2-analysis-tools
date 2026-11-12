function out = build_event_propagation_skeleton(source, varargin)
%BUILD_EVENT_PROPAGATION_SKELETON First-pass multicell propagation skeleton.
%   OUT = BUILD_EVENT_PROPAGATION_SKELETON(RESULT, 'eventUID', '<file>_E###')
%   decomposes the simplified event skeleton into node groups and segments,
%   then assigns participating cells to the nearest segment as a first-pass
%   trunk/branch diagnostic.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'source');
addParameter(p, 'eventUID', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'eventID', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'showFigure', true, @islogical);
addParameter(p, 'visible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
addParameter(p, 'nodeRadiusPx', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
parse(p, source, varargin{:});

S = build_event_structure_diagnostic(source, ...
    'eventUID', p.Results.eventUID, ...
    'eventID', p.Results.eventID, ...
    'showFigure', false, ...
    'delayUnits', 'sec');

skelMask = logical(S.skeleton_mask);
nodePts = [S.branch_points_px; S.end_points_px];
[nodeMask, nodeTbl] = localBuildNodeGroups(skelMask, nodePts, round(double(p.Results.nodeRadiusPx)));
[segMask, segTbl] = localBuildSegments(skelMask, nodeMask, S.pixel_map.delay_sec, S.pixel_map.spatial_res);
cellAssign = localAssignCellsToSegments(S.adjacency.nodes, segTbl, S.pixel_map.spatial_res);

out = struct();
out.event_uid = S.event_uid;
out.event_id = S.event_id;
out.file_id = S.file_id;
out.structure = S;
out.node_mask = nodeMask;
out.nodes = nodeTbl;
out.segment_mask = segMask;
out.segments = segTbl;
out.cell_segment_assignment = cellAssign;

if p.Results.showFigure
    out.figure = localPlotPropagationSkeleton(out, char(string(p.Results.visible)));
end
end

function [nodeMask, nodeTbl] = localBuildNodeGroups(skelMask, nodePts, nodeRadiusPx)
nodeMask = false(size(skelMask));
if isempty(nodePts)
    nodeTbl = table('Size', [0 5], ...
        'VariableTypes', {'double','double','double','double','double'}, ...
        'VariableNames', {'node_id','x_px','y_px','degree_hint','n_pixels'});
    return
end
tmp = false(size(skelMask));
xy = round(nodePts);
xy = xy(all(isfinite(xy),2), :);
xy(:,1) = min(max(xy(:,1),1), size(skelMask,2));
xy(:,2) = min(max(xy(:,2),1), size(skelMask,1));
lin = sub2ind(size(skelMask), xy(:,2), xy(:,1));
tmp(lin) = true;
if nodeRadiusPx > 0
    tmp = imdilate(tmp, strel('disk', nodeRadiusPx, 0));
end
tmp = tmp & skelMask;
cc = bwconncomp(tmp, 8);
nodeMask = false(size(skelMask));
nodeTbl = table('Size', [cc.NumObjects 5], ...
    'VariableTypes', {'double','double','double','double','double'}, ...
    'VariableNames', {'node_id','x_px','y_px','degree_hint','n_pixels'});
for i = 1:cc.NumObjects
    idx = cc.PixelIdxList{i};
    nodeMask(idx) = true;
    [r,c] = ind2sub(size(skelMask), idx);
    nodeTbl.node_id(i) = i;
    nodeTbl.x_px(i) = mean(c);
    nodeTbl.y_px(i) = mean(r);
    nodeTbl.degree_hint(i) = 0;
    nodeTbl.n_pixels(i) = numel(idx);
end
end

function [segMask, segTbl] = localBuildSegments(skelMask, nodeMask, delaySec, spatialRes)
segMask = skelMask & ~nodeMask;
cc = bwconncomp(segMask, 8);
segTbl = table('Size', [cc.NumObjects 9], ...
    'VariableTypes', {'double','double','double','double','double','double','double','double','cell'}, ...
    'VariableNames', {'segment_id','x_px','y_px','length_px','length_um','mean_delay_sec','min_delay_sec','max_delay_sec','pixel_idx'});
for i = 1:cc.NumObjects
    idx = cc.PixelIdxList{i};
    [r,c] = ind2sub(size(segMask), idx);
    vals = delaySec(idx);
    vals = vals(isfinite(vals));
    segTbl.segment_id(i) = i;
    segTbl.x_px(i) = mean(c);
    segTbl.y_px(i) = mean(r);
    segTbl.length_px(i) = numel(idx);
    segTbl.length_um(i) = numel(idx) * spatialRes;
    segTbl.mean_delay_sec(i) = mean(vals, 'omitnan');
    segTbl.min_delay_sec(i) = min(vals, [], 'omitnan');
    segTbl.max_delay_sec(i) = max(vals, [], 'omitnan');
    segTbl.pixel_idx{i} = idx;
end
end

function T = localAssignCellsToSegments(nodes, segTbl, spatialRes)
T = table('Size', [height(nodes) 6], ...
    'VariableTypes', {'string','double','double','double','double','double'}, ...
    'VariableNames', {'cell_label','segment_id','cell_x_um','cell_y_um','segment_x_um','segment_y_um'});
if isempty(segTbl)
    T.cell_label = string(nodes.label);
    T.segment_id(:) = NaN;
    T.cell_x_um = double(nodes.centroid_x_um);
    T.cell_y_um = double(nodes.centroid_y_um);
    T.segment_x_um(:) = NaN;
    T.segment_y_um(:) = NaN;
    return
end
segXY = [double(segTbl.x_px) * spatialRes, double(segTbl.y_px) * spatialRes];
for i = 1:height(nodes)
    T.cell_label(i) = string(nodes.label(i));
    T.cell_x_um(i) = double(nodes.centroid_x_um(i));
    T.cell_y_um(i) = double(nodes.centroid_y_um(i));
    if ~(isfinite(T.cell_x_um(i)) && isfinite(T.cell_y_um(i)))
        T.segment_id(i) = NaN;
        T.segment_x_um(i) = NaN;
        T.segment_y_um(i) = NaN;
        continue
    end
    d = hypot(segXY(:,1) - T.cell_x_um(i), segXY(:,2) - T.cell_y_um(i));
    [~, k] = min(d);
    T.segment_id(i) = double(segTbl.segment_id(k));
    T.segment_x_um(i) = segXY(k,1);
    T.segment_y_um(i) = segXY(k,2);
end
end

function hFig = localPlotPropagationSkeleton(out, visible)
hFig = figure('Name', 'EventPropagationSkeleton', 'Color', 'w', 'NumberTitle', 'off', 'Visible', visible);
ax = axes('Parent', hFig);
hold(ax, 'on');
spatialRes = double(out.structure.pixel_map.spatial_res);
footMask = out.structure.pixel_map.footprint_mask;
sx = size(footMask, 1);
sy = size(footMask, 2);
bg = zeros(sx, sy, 3);
for c = 1:3
    plane = bg(:,:,c);
    plane(footMask) = 0.72;
    bg(:,:,c) = plane;
end
image(ax, [0 sy * spatialRes], [0 sx * spatialRes], bg);
axis(ax, 'image');
set(ax, 'YDir', 'reverse');
xlabel(ax, 'x (um)');
ylabel(ax, 'y (um)');

delayMap = out.structure.pixel_map.delay_sec;
valid = isfinite(delayMap);
if any(valid(:))
    hIm = imagesc(ax, [0 sy * spatialRes], [0 sx * spatialRes], delayMap, 'AlphaData', 0.85 * valid);
    set(hIm, 'AlphaData', 0.85 * valid);
    colormap(ax, turbo(256));
    delayMax = max(delayMap(valid), [], 'omitnan');
    if ~(isfinite(delayMax) && delayMax > 0)
        delayMax = 1;
    end
    caxis(ax, [0 delayMax]);
    cb = colorbar(ax);
    ylabel(cb, 'Post-Origin First Activation Delay (s)');
end

colors = lines(max(height(out.segments), 1));
for i = 1:height(out.segments)
    idx = out.segments.pixel_idx{i};
    [r,c] = ind2sub(size(out.segment_mask), idx);
    plot(ax, c * spatialRes, r * spatialRes, '.', 'Color', colors(i,:), 'MarkerSize', 7);
    text(ax, out.segments.x_px(i) * spatialRes, out.segments.y_px(i) * spatialRes, sprintf('S%d', i), ...
        'Color', colors(i,:), 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

if ~isempty(out.nodes)
    plot(ax, out.nodes.x_px * spatialRes, out.nodes.y_px * spatialRes, 'wo', ...
        'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    for i = 1:height(out.nodes)
        text(ax, out.nodes.x_px(i) * spatialRes, out.nodes.y_px(i) * spatialRes, sprintf('N%d', i), ...
            'Color', 'w', 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    end
end

for i = 1:height(out.cell_segment_assignment)
    if ~(isfinite(out.cell_segment_assignment.cell_x_um(i)) && isfinite(out.cell_segment_assignment.segment_x_um(i)))
        continue
    end
    plot(ax, [out.cell_segment_assignment.cell_x_um(i), out.cell_segment_assignment.segment_x_um(i)], ...
        [out.cell_segment_assignment.cell_y_um(i), out.cell_segment_assignment.segment_y_um(i)], ...
        ':', 'Color', [1 1 1 0.6], 'LineWidth', 1.0);
    text(ax, out.cell_segment_assignment.cell_x_um(i), out.cell_segment_assignment.cell_y_um(i), ...
        char(out.cell_segment_assignment.cell_label(i)), 'Color', [1 1 1], ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

title(ax, sprintf('%s | %s', char(out.file_id), char(out.event_uid)));
hold(ax, 'off');
end
