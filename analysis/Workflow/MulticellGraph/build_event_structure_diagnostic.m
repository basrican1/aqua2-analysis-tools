function out = build_event_structure_diagnostic(source, varargin)
%BUILD_EVENT_STRUCTURE_DIAGNOSTIC Combined multicell structure diagnostic.
%   OUT = BUILD_EVENT_STRUCTURE_DIAGNOSTIC(RESULT, 'eventUID', '<file>_E###')
%   overlays the post-origin pixel activation map, event-within-cell
%   partitions, directed cell adjacency edges, and simple skeleton-derived
%   branch/bottleneck candidates for one event.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'source');
addParameter(p, 'eventUID', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'eventID', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'delayUnits', 'sec', @(x) any(strcmpi(string(x), ["sec","frames"])));
addParameter(p, 'showFigure', true, @islogical);
addParameter(p, 'visible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
addParameter(p, 'simplifyRadiusPx', 2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'spurLengthPx', 8, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
parse(p, source, varargin{:});

pix = build_event_pixel_activation_map(source, ...
    'eventUID', p.Results.eventUID, ...
    'eventID', p.Results.eventID, ...
    'delayUnits', p.Results.delayUnits, ...
    'showFigure', false);

adj = build_event_cell_adjacency_graph(source, ...
    'eventUID', pix.event_uid, ...
    'showFigure', false);

[simplifiedMask, skelMask, branchPts, endpointPts, bottleneckPts, widthPx] = ...
    localExtractStructureFeatures(pix.footprint_mask, round(double(p.Results.simplifyRadiusPx)), round(double(p.Results.spurLengthPx)));

out = struct();
out.event_uid = pix.event_uid;
out.event_id = pix.event_id;
out.file_id = pix.file_id;
out.pixel_map = pix;
out.adjacency = adj;
out.simplified_mask = simplifiedMask;
out.skeleton_mask = skelMask;
out.branch_points_px = branchPts;
out.end_points_px = endpointPts;
out.bottleneck_points_px = bottleneckPts;
out.skeleton_width_px = widthPx;

if p.Results.showFigure
    out.figure = localPlotStructureDiagnostic(out, char(string(p.Results.visible)), char(string(p.Results.delayUnits)));
end
end

function [simpleMask, skelMask, branchPts, endpointPts, bottleneckPts, widthPx] = localExtractStructureFeatures(footMask, simplifyRadiusPx, spurLengthPx)
simpleMask = localSimplifyMask(footMask, simplifyRadiusPx);
if exist('bwskel', 'file') == 2
    skelMask = bwskel(simpleMask);
else
    skelMask = bwmorph(simpleMask, 'skel', Inf);
end
if spurLengthPx > 0
    skelMask = bwmorph(skelMask, 'spur', spurLengthPx);
end
branchMask = bwmorph(skelMask, 'branchpoints');
endMask = bwmorph(skelMask, 'endpoints');
dt = bwdist(~simpleMask);
widthPx = nan(size(dt));
widthPx(skelMask) = 2 * dt(skelMask);
branchPts = localMaskToPoints(branchMask);
endpointPts = localMaskToPoints(endMask);

skelWidths = widthPx(skelMask);
skelWidths = skelWidths(isfinite(skelWidths) & skelWidths > 0);
if isempty(skelWidths)
    bottleneckPts = zeros(0,2);
    return
end

function simpleMask = localSimplifyMask(mask, radiusPx)
simpleMask = logical(mask);
if ~any(simpleMask(:))
    return
end
simpleMask = imfill(simpleMask, 'holes');
if radiusPx > 0
    se = strel('disk', radiusPx, 0);
    simpleMask = imopen(simpleMask, se);
    simpleMask = imclose(simpleMask, se);
    simpleMask = imfill(simpleMask, 'holes');
end
cc = bwconncomp(simpleMask, 8);
if cc.NumObjects > 1
    sizes = cellfun(@numel, cc.PixelIdxList);
    [~, imax] = max(sizes);
    keep = false(size(simpleMask));
    keep(cc.PixelIdxList{imax}) = true;
    simpleMask = keep;
end
end
cutoff = prctile(skelWidths, 25);
bottleMask = skelMask & isfinite(widthPx) & widthPx <= cutoff;
cc = bwconncomp(bottleMask, 8);
bottleneckPts = zeros(cc.NumObjects, 2);
for i = 1:cc.NumObjects
    idx = cc.PixelIdxList{i};
    [r,c] = ind2sub(size(bottleMask), idx);
    bottleneckPts(i,:) = [mean(c), mean(r)];
end
end

function pts = localMaskToPoints(mask)
[r,c] = find(mask);
pts = [c, r];
end

function hFig = localPlotStructureDiagnostic(out, visible, delayUnits)
hFig = figure('Name', 'EventStructureDiagnostic', 'Color', 'w', 'NumberTitle', 'off', 'Visible', visible);
ax = axes('Parent', hFig);
hold(ax, 'on');
spatialRes = double(out.pixel_map.spatial_res);
footMask = out.pixel_map.footprint_mask;
simpleMask = out.simplified_mask;
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

switch lower(delayUnits)
    case 'frames'
        delayMap = out.pixel_map.delay_frames;
        cbLabel = 'Post-Origin First Activation Delay (frames)';
    otherwise
        delayMap = out.pixel_map.delay_sec;
        cbLabel = 'Post-Origin First Activation Delay (s)';
end
valid = isfinite(delayMap);
if any(valid(:))
    hIm = imagesc(ax, [0 sy * spatialRes], [0 sx * spatialRes], delayMap, 'AlphaData', 0.90 * valid);
    set(hIm, 'AlphaData', 0.90 * valid);
    colormap(ax, turbo(256));
    delayMin = 0;
    delayMax = max(delayMap(valid), [], 'omitnan');
    if ~(isfinite(delayMax) && delayMax > delayMin)
        delayMax = delayMin + 1;
    end
    caxis(ax, [delayMin delayMax]);
    cb = colorbar(ax);
    ylabel(cb, cbLabel);
end

if any(simpleMask(:))
    Bsimple = bwboundaries(simpleMask, 8, 'noholes');
    for i = 1:numel(Bsimple)
        b = Bsimple{i};
        plot(ax, b(:,2) * spatialRes, b(:,1) * spatialRes, '-', 'Color', [0.2 1 1], 'LineWidth', 1.2);
    end
end

nodes = out.adjacency.nodes;
for i = 1:height(nodes)
    b = nodes.boundary_xy_px{i};
    if isempty(b)
        continue
    end
    plot(ax, b(:,1) * spatialRes, b(:,2) * spatialRes, '-', 'Color', [1 1 1], 'LineWidth', 1.0);
    if isfinite(nodes.centroid_x_px(i)) && isfinite(nodes.centroid_y_px(i))
        text(ax, nodes.centroid_x_um(i), nodes.centroid_y_um(i), char(nodes.label(i)), ...
            'Color', [1 1 1], 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
end

edges = out.adjacency.edges;
for i = 1:height(edges)
    s = edges.source_idx(i);
    t = edges.target_idx(i);
    if ~all(isfinite([nodes.centroid_x_um(s), nodes.centroid_y_um(s), nodes.centroid_x_um(t), nodes.centroid_y_um(t)]))
        continue
    end
    x1 = nodes.centroid_x_um(s);
    y1 = nodes.centroid_y_um(s);
    x2 = nodes.centroid_x_um(t);
    y2 = nodes.centroid_y_um(t);
    plot(ax, [x1 x2], [y1 y2], '-', 'Color', [0.95 0.55 0.10], 'LineWidth', 1.5);
    if logical(edges.time_forward(i))
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

if ~isempty(out.branch_points_px)
    plot(ax, out.branch_points_px(:,1) * spatialRes, out.branch_points_px(:,2) * spatialRes, ...
        'p', 'MarkerSize', 11, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0 1]);
end
if ~isempty(out.bottleneck_points_px)
    plot(ax, out.bottleneck_points_px(:,1) * spatialRes, out.bottleneck_points_px(:,2) * spatialRes, ...
        'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0.6 1], 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
end
if ~isempty(out.end_points_px)
    plot(ax, out.end_points_px(:,1) * spatialRes, out.end_points_px(:,2) * spatialRes, ...
        's', 'MarkerSize', 7, 'MarkerEdgeColor', [0.2 1 1], 'MarkerFaceColor', 'none', 'LineWidth', 1.0);
end

title(ax, sprintf('%s | %s', char(out.file_id), char(out.event_uid)));
hold(ax, 'off');
end
