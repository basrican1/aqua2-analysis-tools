function out = build_event_corridor_bundle(source, varargin)
%BUILD_EVENT_CORRIDOR_BUNDLE Bundle top-ranked upstream corridors for an event.
%   OUT = BUILD_EVENT_CORRIDOR_BUNDLE(RESULT, 'eventUID', '<file>_E###')

p = inputParser;
p.StructExpand = false;
addRequired(p, 'source');
addParameter(p, 'eventUID', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'eventID', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'showFigure', true, @islogical);
addParameter(p, 'visible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
addParameter(p, 'sourceAnchorInteriorFrac', 0.9, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
addParameter(p, 'targetAnchorInteriorFrac', 0.75, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
addParameter(p, 'edgeBandFrac', 0.20, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
addParameter(p, 'turnPenalty', 1.0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'corridorRadiusPx', 6, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'corridorCostSlack', 1.35, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'corridorAdaptive', true, @islogical);
addParameter(p, 'targetWidthWeight', 80, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'maxCandidatesPerTarget', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
parse(p, source, varargin{:});

G = build_event_cell_adjacency_graph(source, ...
    'eventUID', p.Results.eventUID, ...
    'eventID', p.Results.eventID, ...
    'showFigure', false);
M = build_event_pixel_activation_map(source, ...
    'eventUID', G.event_uid, ...
    'showFigure', false, ...
    'delayUnits', 'sec');

corridors = cell(0,1);
bundleMask = zeros(size(M.footprint_mask));
rows = table('Size', [0 4], ...
    'VariableTypes', {'string','string','double','double'}, ...
    'VariableNames', {'source_label','target_label','score','n_path_pixels'});

for i = 1:height(G.nodes)
    if G.nodes.is_primary_cell(i)
        continue
    end
    U = build_event_upstream_candidates(source, ...
        'eventUID', G.event_uid, ...
        'targetCell', G.nodes.label(i), ...
        'showFigure', false);
    if isempty(U.candidates)
        continue
    end
    nTry = min(height(U.candidates), max(1, round(double(p.Results.maxCandidatesPerTarget))));
    bestC = [];
    bestScore = inf;
    bestSrc = "";
    bestTgt = string(G.nodes.label(i));
    bestRankScore = NaN;
    for k = 1:nTry
        src = string(U.candidates.source_label(k));
        tgt = string(U.candidates.target_label(k));
        C = build_event_pair_corridor(source, ...
            'eventUID', G.event_uid, ...
            'sourceCell', src, ...
            'targetCell', tgt, ...
            'showFigure', false, ...
            'sourceAnchorInteriorFrac', double(p.Results.sourceAnchorInteriorFrac), ...
            'targetAnchorInteriorFrac', double(p.Results.targetAnchorInteriorFrac), ...
            'edgeBandFrac', double(p.Results.edgeBandFrac), ...
            'turnPenalty', double(p.Results.turnPenalty), ...
            'corridorRadiusPx', double(p.Results.corridorRadiusPx), ...
            'corridorCostSlack', double(p.Results.corridorCostSlack), ...
            'corridorAdaptive', logical(p.Results.corridorAdaptive), ...
            'targetWidthWeight', double(p.Results.targetWidthWeight));
        qual = localBundleCorridorScore(C, double(U.candidates.score(k)));
        if qual < bestScore
            bestScore = qual;
            bestC = C;
            bestSrc = src;
            bestTgt = tgt;
            bestRankScore = double(U.candidates.score(k));
        end
        if isfield(C,'refine_outcome') && isstruct(C.refine_outcome) && isfield(C.refine_outcome,'acceptable') && C.refine_outcome.acceptable
            break
        end
    end
    if isempty(bestC)
        continue
    end
    corridors{end+1,1} = bestC; %#ok<AGROW>
    bundleMask = bundleMask + double(bestC.support_mask);
    rows(end+1,:) = {bestSrc, bestTgt, bestRankScore, size(bestC.path_xy_px,1)}; %#ok<AGROW>
end

out = struct();
out.event_uid = G.event_uid;
out.file_id = G.file_id;
out.graph = G;
out.pixel_map = M;
out.corridors = corridors;
out.corridor_table = rows;
out.bundle_mask = bundleMask;

if p.Results.showFigure
    out.figure = localPlotCorridorBundle(out, char(string(p.Results.visible)));
end
end

function score = localBundleCorridorScore(C, rankScore)
score = inf;
if isempty(C) || ~isfield(C,'refine_outcome') || ~isstruct(C.refine_outcome)
    return
end
accPenalty = 0;
if ~(isfield(C.refine_outcome,'acceptable') && C.refine_outcome.acceptable)
    accPenalty = 500;
end
st = struct();
if isfield(C.refine_outcome,'stats')
    st = C.refine_outcome.stats;
end
meanWidth = 0; nPath = size(C.path_xy_px,1);
if isfield(st,'meanWidth'), meanWidth = double(st.meanWidth); end
score = accPenalty + 10*meanWidth - 0.5*nPath - rankScore;
end

function hFig = localPlotCorridorBundle(out, visible)
hFig = figure('Name', 'EventCorridorBundle', 'Color', 'w', 'NumberTitle', 'off', 'Visible', visible);
ax = axes('Parent', hFig);
hold(ax, 'on');
spatialRes = double(out.pixel_map.spatial_res);
footMask = out.pixel_map.footprint_mask;
sx = size(footMask,1);
sy = size(footMask,2);
bg = zeros(sx, sy, 3);
for c = 1:3
    plane = bg(:,:,c);
    plane(footMask) = 0.72;
    bg(:,:,c) = plane;
end
image(ax, [0 sy*spatialRes], [0 sx*spatialRes], bg);
axis(ax, 'image');
set(ax, 'YDir', 'reverse');
xlabel(ax, 'x (um)');
ylabel(ax, 'y (um)');

delayMap = out.pixel_map.delay_sec;
valid = isfinite(delayMap);
if any(valid(:))
    hIm = imagesc(ax, [0 sy*spatialRes], [0 sx*spatialRes], delayMap, 'AlphaData', 0.65 * valid);
    set(hIm, 'AlphaData', 0.65 * valid);
    colormap(ax, turbo(256));
    clim(ax, [0 max(delayMap(valid), [], 'omitnan') + eps]);
end

if any(out.bundle_mask(:) > 0)
    hB = imagesc(ax, [0 sy*spatialRes], [0 sx*spatialRes], out.bundle_mask, 'AlphaData', 0.55 * (out.bundle_mask > 0));
    uistack(hB, 'top');
    [CC, RR] = meshgrid((1:sy) * spatialRes, (1:sx) * spatialRes);
    contour(ax, CC, RR, double(out.bundle_mask > 0), [1 1], ...
        'LineColor', [1 1 1], 'LineWidth', 1.2);
end

bundleColors = lines(max(numel(out.corridors), 1));
for i = 1:numel(out.corridors)
    C = out.corridors{i};
    if isempty(C.path_xy_px), continue; end
    x = C.path_xy_px(:,1) * spatialRes;
    y = C.path_xy_px(:,2) * spatialRes;
    plot(ax, x, y, '-', 'Color', [0 0 0], 'LineWidth', 3.0);
    plot(ax, x, y, '-', 'Color', bundleColors(i,:), 'LineWidth', 1.8);
    plot(ax, x, y, '.', 'Color', bundleColors(i,:), 'MarkerSize', 10);
    plot(ax, x(1), y(1), 'o', 'Color', bundleColors(i,:), 'MarkerFaceColor', bundleColors(i,:), 'MarkerSize', 5);
    plot(ax, x(end), y(end), 's', 'Color', bundleColors(i,:), 'MarkerFaceColor', bundleColors(i,:), 'MarkerSize', 5);
    if ~isempty(x)
        text(ax, x(end), y(end), sprintf('%s->%s', char(C.source_cell), char(C.target_cell)), ...
            'Color', bundleColors(i,:), 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
end

nodes = out.graph.nodes;
for i = 1:height(nodes)
    if isfinite(nodes.centroid_x_um(i)) && isfinite(nodes.centroid_y_um(i))
        text(ax, nodes.centroid_x_um(i), nodes.centroid_y_um(i), char(nodes.label(i)), ...
            'Color', [1 1 1], 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
end

title(ax, sprintf('%s | corridor bundle', char(out.event_uid)));
hold(ax, 'off');
end
