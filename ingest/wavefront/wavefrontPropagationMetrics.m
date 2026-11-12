function [out, debug] = wavefrontPropagationMetrics(propGrow, propShrink, varargin)
%WAVEFRONTPROPAGATIONMETRICS Prototype leading-edge propagation metrics.
%   [OUT, DEBUG] = WAVEFRONTPROPAGATIONMETRICS(PROPGROW, PROPSHRINK)
%   reconstructs the 4-direction wavefront implied by AQuA2 propagation
%   arrays and derives leading-edge path metrics by backward-tracking from
%   the furthest attained frontier point.
%
%   Name-Value:
%     'spatialRes'        : calibration in um/px (default 1)
%     'minValidFrames'    : minimum propagated frames required (default 3)
%     'continuityWeight'  : weight on frame-to-frame continuity in backward
%                           tracking (default 1.0)
%     'radialWeight'      : weight on outward radial distance in backward
%                           tracking (default 0.35)
%     'maxJumpPx'         : optional hard cap on backward jump size; larger
%                           candidates are rejected (default inf)
%
%   OUT fields:
%     wavefront_max_reach_px
%     wavefront_max_reach_um
%     wavefront_frontier_length_px
%     wavefront_frontier_length_um
%     wavefront_net_displacement_px
%     wavefront_net_displacement_um
%     wavefront_tortuosity
%     wavefront_peak_frame
%     wavefront_peak_vertex_index
%     wavefront_valid_frames
%     wavefront_angle_deg
%
%   DEBUG contains the reconstructed centers, 4-point wavefront vertices,
%   chosen backward-tracked frontier path, and validity masks.

p = inputParser;
p.StructExpand = false;
addParameter(p, 'spatialRes', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'minValidFrames', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'continuityWeight', 1.0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'radialWeight', 0.35, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'maxJumpPx', inf, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
parse(p, varargin{:});

out = emptyWavefrontMetrics();
debug = struct( ...
    'centers_px', zeros(0, 2), ...
    'vertices_px', zeros(0, 2, 4), ...
    'frameValid', false(0, 1), ...
    'frameActivity', zeros(0, 1), ...
    'origin_center_px', [NaN NaN], ...
    'peak_frontier_px', [NaN NaN], ...
    'frontier_track_px', zeros(0, 2), ...
    'frontier_track_valid', false(0, 1), ...
    'reach_px', zeros(0, 4));

if isempty(propGrow) || isempty(propShrink) || ~isnumeric(propGrow) || ~isnumeric(propShrink)
    return
end
if ~isequal(size(propGrow), size(propShrink)) || size(propGrow, 2) < 4
    return
end

nFrames = size(propGrow, 1);
centers = NaN(nFrames, 2);
vertices = NaN(nFrames, 2, 4);
frameValid = false(nFrames, 1);
frameActivity = zeros(nFrames, 1);

prevNorth = [0 0];
prevWest = [0 0];
prevSouth = [0 0];
prevEast = [0 0];
prevCenter = [0 0];

for eFrame = 1:nFrames
    deltaUp = propGrow(eFrame, 1) + propShrink(eFrame, 1);
    deltaDown = -propGrow(eFrame, 2) - propShrink(eFrame, 2);
    deltaLeft = -propGrow(eFrame, 3) - propShrink(eFrame, 3);
    deltaRight = propGrow(eFrame, 4) + propShrink(eFrame, 4);
    frameActivity(eFrame) = abs(deltaUp) + abs(deltaDown) + abs(deltaLeft) + abs(deltaRight);

    dxN = signedDelta(prevNorth(1), prevCenter(1), deltaLeft, deltaRight);
    dxS = signedDelta(prevSouth(1), prevCenter(1), deltaLeft, deltaRight);
    dyW = signedDelta(prevWest(2), prevCenter(2), deltaDown, deltaUp);
    dyE = signedDelta(prevEast(2), prevCenter(2), deltaDown, deltaUp);

    newNorth = [prevNorth(1) + dxN, prevNorth(2) + deltaUp];
    newWest = [prevWest(1) + deltaLeft, prevWest(2) + dyW];
    newSouth = [prevSouth(1) + dxS, prevSouth(2) + deltaDown];
    newEast = [prevEast(1) + deltaRight, prevEast(2) + dyE];
    pts = [newNorth; newWest; newSouth; newEast];

    [cx, cy] = safeCentroid(pts(:,1), pts(:,2));
    if isfinite(cx) && isfinite(cy) && frameActivity(eFrame) > 0
        centers(eFrame, :) = [cx cy];
        vertices(eFrame, :, :) = permute(pts, [3 2 1]);
        frameValid(eFrame) = true;
        prevCenter = [cx cy];
    end

    prevNorth = newNorth;
    prevWest = newWest;
    prevSouth = newSouth;
    prevEast = newEast;
end

debug.centers_px = centers;
debug.vertices_px = vertices;
debug.frameValid = frameValid;
debug.frameActivity = frameActivity;

validIdx = find(frameValid);
if numel(validIdx) < double(p.Results.minValidFrames)
    return
end

originIdx = validIdx(1);
originCenter = centers(originIdx, :);
debug.origin_center_px = originCenter;

reach = NaN(nFrames, 4);
peakReachPx = -inf;
peakFrame = NaN;
peakVertex = NaN;
peakPoint = [NaN NaN];
for k = 1:numel(validIdx)
    t = validIdx(k);
    pts = squeeze(vertices(t, :, :)).';
    if isempty(pts)
        continue
    end
    for i = 1:4
        if all(isfinite(pts(i, :)))
            reach(t, i) = hypot(pts(i,1) - originCenter(1), pts(i,2) - originCenter(2));
            if reach(t, i) > peakReachPx
                peakReachPx = reach(t, i);
                peakFrame = t;
                peakVertex = i;
                peakPoint = pts(i, :);
            end
        end
    end
end
debug.reach_px = reach;
debug.peak_frontier_px = peakPoint;

if ~(isfinite(peakReachPx) && isfinite(peakFrame) && isfinite(peakVertex))
    return
end

track = NaN(nFrames, 2);
trackValid = false(nFrames, 1);
track(peakFrame, :) = peakPoint;
trackValid(peakFrame) = true;

contW = double(p.Results.continuityWeight);
radW = double(p.Results.radialWeight);
maxJumpPx = double(p.Results.maxJumpPx);
for k = find(validIdx < peakFrame, 1, 'last'):-1:1
    t = validIdx(k);
    pts = squeeze(vertices(t, :, :)).';
    nextPoint = track(validIdx(k + 1), :);
    bestScore = -inf;
    bestPoint = [NaN NaN];
    for i = 1:4
        cand = pts(i, :);
        if ~all(isfinite(cand))
            continue
        end
        jumpPx = hypot(cand(1) - nextPoint(1), cand(2) - nextPoint(2));
        if jumpPx > maxJumpPx
            continue
        end
        radialPx = hypot(cand(1) - originCenter(1), cand(2) - originCenter(2));
        score = -contW * jumpPx + radW * radialPx;
        if score > bestScore
            bestScore = score;
            bestPoint = cand;
        end
    end
    if all(isfinite(bestPoint))
        track(t, :) = bestPoint;
        trackValid(t) = true;
    end
end

debug.frontier_track_px = track;
debug.frontier_track_valid = trackValid;

trackedIdx = find(trackValid);
if numel(trackedIdx) < 2
    out.wavefront_max_reach_px = peakReachPx;
    out.wavefront_max_reach_um = peakReachPx * p.Results.spatialRes;
    out.wavefront_peak_frame = peakFrame;
    out.wavefront_peak_vertex_index = peakVertex;
    out.wavefront_valid_frames = numel(validIdx);
    out.wavefront_angle_deg = atan2d(peakPoint(2) - originCenter(2), peakPoint(1) - originCenter(1));
    return
end

% Include the initial origin-to-frontier segment so frontier length answers
% "how far along its path did the leading edge travel from the event origin?"
frontierLengthPx = hypot(track(trackedIdx(1),1) - originCenter(1), ...
    track(trackedIdx(1),2) - originCenter(2));
for i = 1:(numel(trackedIdx) - 1)
    a = track(trackedIdx(i), :);
    b = track(trackedIdx(i + 1), :);
    frontierLengthPx = frontierLengthPx + hypot(b(1) - a(1), b(2) - a(2));
end

firstTracked = track(trackedIdx(1), :);
lastTracked = track(trackedIdx(end), :);
netDispPx = hypot(lastTracked(1) - firstTracked(1), lastTracked(2) - firstTracked(2));

out.wavefront_max_reach_px = peakReachPx;
out.wavefront_max_reach_um = peakReachPx * p.Results.spatialRes;
out.wavefront_frontier_length_px = frontierLengthPx;
out.wavefront_frontier_length_um = frontierLengthPx * p.Results.spatialRes;
out.wavefront_net_displacement_px = netDispPx;
out.wavefront_net_displacement_um = netDispPx * p.Results.spatialRes;
if peakReachPx > 0
    out.wavefront_tortuosity = frontierLengthPx / peakReachPx;
end
out.wavefront_peak_frame = peakFrame;
out.wavefront_peak_vertex_index = peakVertex;
out.wavefront_valid_frames = numel(validIdx);
out.wavefront_angle_deg = atan2d(peakPoint(2) - originCenter(2), peakPoint(1) - originCenter(1));
end

function out = emptyWavefrontMetrics()
out = struct( ...
    'wavefront_max_reach_px', NaN, ...
    'wavefront_max_reach_um', NaN, ...
    'wavefront_frontier_length_px', NaN, ...
    'wavefront_frontier_length_um', NaN, ...
    'wavefront_net_displacement_px', NaN, ...
    'wavefront_net_displacement_um', NaN, ...
    'wavefront_tortuosity', NaN, ...
    'wavefront_peak_frame', NaN, ...
    'wavefront_peak_vertex_index', NaN, ...
    'wavefront_valid_frames', 0, ...
    'wavefront_angle_deg', NaN);
end

function val = signedDelta(edgeCoord, centerCoord, negDelta, posDelta)
if edgeCoord > centerCoord
    val = posDelta;
elseif edgeCoord < centerCoord
    val = negDelta;
else
    val = 0;
end
end

function [cx, cy] = safeCentroid(x, y)
cx = NaN;
cy = NaN;
if numel(x) < 3 || numel(y) < 3
    return
end
try
    pg = polyshape(x, y, 'Simplify', false);
    [cx, cy] = centroid(pg);
catch
    cx = NaN;
    cy = NaN;
end
end
