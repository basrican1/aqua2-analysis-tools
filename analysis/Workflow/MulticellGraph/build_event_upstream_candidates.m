function out = build_event_upstream_candidates(source, varargin)
%BUILD_EVENT_UPSTREAM_CANDIDATES Rank plausible upstream cells for one target.
%   OUT = BUILD_EVENT_UPSTREAM_CANDIDATES(RESULT, 'eventUID', '<file>_E###', ...
%       'targetCell', 'C23')
%   uses the directed event-cell adjacency graph plus the post-origin pixel
%   activation map to rank plausible immediate upstream parents of one cell.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'source');
addParameter(p, 'eventUID', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'eventID', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'targetCell', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'showFigure', true, @islogical);
addParameter(p, 'visible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
parse(p, source, varargin{:});

G = build_event_cell_adjacency_graph(source, ...
    'eventUID', p.Results.eventUID, ...
    'eventID', p.Results.eventID, ...
    'showFigure', false);
M = build_event_pixel_activation_map(source, ...
    'eventUID', G.event_uid, ...
    'delayUnits', 'sec', ...
    'showFigure', false);

targetCell = string(p.Results.targetCell);
if strlength(targetCell) == 0
    error('build_event_upstream_candidates:MissingTarget', 'Provide targetCell like C23.');
end
targetIdx = localResolveTargetNodeIdx(G.nodes, targetCell);
if isempty(targetIdx)
    error('build_event_upstream_candidates:TargetNotFound', 'Target cell %s not found in event %s.', char(targetCell), char(G.event_uid));
end

inMask = G.edges.target_idx == targetIdx;
E = G.edges(inMask, :);
rows = table('Size', [0 11], ...
    'VariableTypes', {'string','string','double','double','double','logical','double','double','double','double','double'}, ...
    'VariableNames', {'source_label','target_label','source_delay_frames','target_delay_frames', ...
    'delay_diff_frames','touches','min_gap_um','path_valid_frac','path_monotone_frac','path_mean_delay_step_sec','score'});

for i = 1:height(E)
    s = E.source_idx(i);
    t = E.target_idx(i);
    [validFrac, monotoneFrac, meanStepSec] = localDelayPathStats(M.delay_sec, M.footprint_mask, ...
        G.nodes.centroid_x_px(s), G.nodes.centroid_y_px(s), ...
        G.nodes.centroid_x_px(t), G.nodes.centroid_y_px(t));
    delayDiff = double(E.target_delay_frames(i) - E.source_delay_frames(i));
    score = localCandidateScore(double(E.min_gap_um(i)), logical(E.touches(i)), validFrac, monotoneFrac, delayDiff);
    rows(end+1,:) = {string(G.nodes.label(s)), string(G.nodes.label(t)), ...
        double(E.source_delay_frames(i)), double(E.target_delay_frames(i)), delayDiff, ...
        logical(E.touches(i)), double(E.min_gap_um(i)), validFrac, monotoneFrac, meanStepSec, score}; %#ok<AGROW>
end

if ~isempty(rows)
    rows = sortrows(rows, {'score','min_gap_um','delay_diff_frames'}, {'descend','ascend','ascend'});
end

out = struct();
out.event_uid = G.event_uid;
out.target_cell = string(G.nodes.label(targetIdx));
out.target_cell_uid = string(G.nodes.cell_uid(targetIdx));
out.target_idx = targetIdx;
out.graph = G;
out.pixel_map = M;
out.candidates = rows;

if p.Results.showFigure
    out.figure = localPlotUpstreamCandidates(out, char(string(p.Results.visible)));
end
end

function idx = localResolveTargetNodeIdx(nodes, targetCell)
idx = [];
targetCell = string(targetCell);
hit = find(string(nodes.label) == targetCell, 1, 'first');
if ~isempty(hit)
    idx = hit;
    return
end
hit = find(string(nodes.cell_uid) == targetCell, 1, 'first');
if ~isempty(hit)
    idx = hit;
end
end

function [validFrac, monotoneFrac, meanStepSec] = localDelayPathStats(delaySec, footMask, x1, y1, x2, y2)
validFrac = NaN;
monotoneFrac = NaN;
meanStepSec = NaN;
if ~all(isfinite([x1 y1 x2 y2]))
    return
end
n = max(2, ceil(hypot(x2 - x1, y2 - y1)));
xs = round(linspace(x1, x2, n));
ys = round(linspace(y1, y2, n));
sx = size(delaySec, 2);
sy = size(delaySec, 1);
keep = xs >= 1 & xs <= sx & ys >= 1 & ys <= sy;
xs = xs(keep);
ys = ys(keep);
if isempty(xs)
    return
end
lin = sub2ind(size(delaySec), ys, xs);
inside = footMask(lin);
validFrac = mean(double(inside));
vals = delaySec(lin);
vals = vals(isfinite(vals));
if numel(vals) < 2
    return
end
d = diff(vals);
monotoneFrac = mean(double(d >= -1e-9));
meanStepSec = mean(d, 'omitnan');
end

function score = localCandidateScore(minGapUm, touches, validFrac, monotoneFrac, delayDiffFrames)
score = 0;
if touches
    score = score + 2.0;
end
if isfinite(minGapUm)
    score = score - 0.15 * minGapUm;
end
if isfinite(validFrac)
    score = score + 2.0 * validFrac;
end
if isfinite(monotoneFrac)
    score = score + 2.5 * monotoneFrac;
end
if isfinite(delayDiffFrames)
    if delayDiffFrames >= 0
        score = score + 0.5;
    else
        score = score - 2.0;
    end
end
end

function hFig = localPlotUpstreamCandidates(out, visible)
hFig = figure('Name', 'EventUpstreamCandidates', 'Color', 'w', 'NumberTitle', 'off', 'Visible', visible);
ax = axes('Parent', hFig);
hold(ax, 'on');
spatialRes = double(out.pixel_map.spatial_res);
footMask = out.pixel_map.footprint_mask;
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

delayMap = out.pixel_map.delay_sec;
valid = isfinite(delayMap);
if any(valid(:))
    hIm = imagesc(ax, [0 sy * spatialRes], [0 sx * spatialRes], delayMap, 'AlphaData', 0.90 * valid);
    set(hIm, 'AlphaData', 0.90 * valid);
    colormap(ax, turbo(256));
    delayMax = max(delayMap(valid), [], 'omitnan');
    if ~(isfinite(delayMax) && delayMax > 0)
        delayMax = 1;
    end
    caxis(ax, [0 delayMax]);
    cb = colorbar(ax);
    ylabel(cb, 'Post-Origin First Activation Delay (s)');
end

nodes = out.graph.nodes;
targetIdx = out.target_idx;
for i = 1:height(nodes)
    b = nodes.boundary_xy_px{i};
    if isempty(b)
        continue
    end
    color = [1 1 1];
    lw = 1.0;
    if i == targetIdx
        color = [1 0.2 0.2];
        lw = 2.0;
    end
    plot(ax, b(:,1) * spatialRes, b(:,2) * spatialRes, '-', 'Color', color, 'LineWidth', lw);
    if isfinite(nodes.centroid_x_px(i)) && isfinite(nodes.centroid_y_px(i))
        text(ax, nodes.centroid_x_um(i), nodes.centroid_y_um(i), char(nodes.label(i)), ...
            'Color', color, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
end

for i = 1:height(out.candidates)
    s = find(string(nodes.label) == string(out.candidates.source_label(i)), 1, 'first');
    if isempty(s)
        continue
    end
    x1 = nodes.centroid_x_um(s);
    y1 = nodes.centroid_y_um(s);
    x2 = nodes.centroid_x_um(targetIdx);
    y2 = nodes.centroid_y_um(targetIdx);
    plot(ax, [x1 x2], [y1 y2], '-', 'Color', [0.95 0.55 0.10], 'LineWidth', 1.8);
    dx = x2 - x1;
    dy = y2 - y1;
    if hypot(dx, dy) > 0
        qx = x1 + 0.70 * dx;
        qy = y1 + 0.70 * dy;
        quiver(ax, qx, qy, 0.20 * dx, 0.20 * dy, 0, ...
            'Color', [0.95 0.55 0.10], 'LineWidth', 1.2, ...
            'MaxHeadSize', 2.0, 'AutoScale', 'off');
    end
end

title(ax, sprintf('%s | target %s', char(out.event_uid), char(out.target_cell)));
hold(ax, 'off');
end
