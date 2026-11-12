function out = canonicalEventMetrics(evtPix, spatialSize, cellMaskOrPix, propGrow, propShrink)
% Compute canonical event geometry/propagation metrics from raw spatial traces.

out = struct( ...
    'perimeter_px', NaN, ...
    'circularity', NaN, ...
    'projected_area_npix', NaN, ...
    'propagation_distance_px', NaN, ...
    'propagation_angle_deg', NaN, ...
    'cell_overlap_npix', NaN, ...
    'cell_overlap_frac_of_event', NaN, ...
    'cell_overlap_frac_of_cell', NaN, ...
    'frac_outside_cell', NaN, ...
    'outside_cell_area_npix', NaN);

evtPix = forcePix(evtPix);
if ~isempty(evtPix)
    shp = computeShape(evtPix, spatialSize);
    out.perimeter_px = shp.perimeter_px;
    out.circularity = shp.circularity;
    out.projected_area_npix = shp.projected_area_npix;

    ov = computeCellOverlap(evtPix, cellMaskOrPix);
    out.cell_overlap_npix = ov.cell_overlap_npix;
    out.cell_overlap_frac_of_event = ov.cell_overlap_frac_of_event;
    out.cell_overlap_frac_of_cell = ov.cell_overlap_frac_of_cell;
    out.frac_outside_cell = ov.frac_outside_cell;
    out.outside_cell_area_npix = ov.outside_cell_area_npix;
end

[distPx, angleDeg] = computePropagation(propGrow, propShrink);
out.propagation_distance_px = distPx;
out.propagation_angle_deg = angleDeg;
end

function shp = computeShape(evtPix, spatialSize)
shp = struct('perimeter_px', NaN, 'circularity', NaN, 'projected_area_npix', NaN);
bw = spatialProjectionMask(evtPix, spatialSize);
if isempty(bw)
    return
end
projArea = nnz(bw);
shp.projected_area_npix = projArea;
if projArea == 0
    return
end
perim = boundaryPixelCount(bw);
shp.perimeter_px = perim;
if isfinite(perim) && perim > 0
    shp.circularity = min(1, 4 * pi * projArea / (perim ^ 2));
end
end

function ov = computeCellOverlap(evtPix, cellMaskOrPix)
ov = struct( ...
    'cell_overlap_npix', NaN, ...
    'cell_overlap_frac_of_event', NaN, ...
    'cell_overlap_frac_of_cell', NaN, ...
    'frac_outside_cell', NaN, ...
    'outside_cell_area_npix', NaN);

evtPix = forcePix(evtPix);
cellPix = [];
if isnumeric(cellMaskOrPix) || islogical(cellMaskOrPix)
    if islogical(cellMaskOrPix) || ~isvector(cellMaskOrPix)
        cellPix = find(cellMaskOrPix > 0);
    else
        cellPix = double(cellMaskOrPix(:));
    end
end
cellPix = forcePix(cellPix);
if isempty(evtPix) || isempty(cellPix)
    return
end

overlapN = numel(intersect(evtPix, cellPix));
ov.cell_overlap_npix = overlapN;
ov.cell_overlap_frac_of_event = overlapN / numel(evtPix);
ov.cell_overlap_frac_of_cell = overlapN / numel(cellPix);
ov.frac_outside_cell = 1 - ov.cell_overlap_frac_of_event;
ov.outside_cell_area_npix = numel(evtPix) - overlapN;
end

function [distPx, angleDeg] = computePropagation(propGrow, propShrink)
distPx = NaN;
angleDeg = NaN;
if isempty(propGrow) || isempty(propShrink) || ~isnumeric(propGrow) || ~isnumeric(propShrink)
    return
end
if ~isequal(size(propGrow), size(propShrink)) || size(propGrow,2) < 4
    return
end

nEventFrames = size(propGrow, 1);
prevNorth = [0 0];
prevWest = [0 0];
prevSouth = [0 0];
prevEast = [0 0];
prevCenter = [0 0];
distPx = 0;
cx = 0;
cy = 0;

for eFrame = 1:nEventFrames
    deltaUp = propGrow(eFrame, 1) + propShrink(eFrame, 1);
    deltaDown = -propGrow(eFrame, 2) - propShrink(eFrame, 2);
    deltaLeft = -propGrow(eFrame, 3) - propShrink(eFrame, 3);
    deltaRight = propGrow(eFrame, 4) + propShrink(eFrame, 4);

    dxN = signedDelta(prevNorth(1), prevCenter(1), deltaLeft, deltaRight);
    dxS = signedDelta(prevSouth(1), prevCenter(1), deltaLeft, deltaRight);
    dyW = signedDelta(prevWest(2), prevCenter(2), deltaDown, deltaUp);
    dyE = signedDelta(prevEast(2), prevCenter(2), deltaDown, deltaUp);

    newNorth = [prevNorth(1) + dxN, prevNorth(2) + deltaUp];
    newWest = [prevWest(1) + deltaLeft, prevWest(2) + dyW];
    newSouth = [prevSouth(1) + dxS, prevSouth(2) + deltaDown];
    newEast = [prevEast(1) + deltaRight, prevEast(2) + dyE];

    [testx, testy] = safeCentroid([newNorth(1) newWest(1) newSouth(1) newEast(1)], ...
        [newNorth(2) newWest(2) newSouth(2) newEast(2)]);
    if ~(isfinite(testx) && isfinite(testy))
        testx = prevCenter(1);
        testy = prevCenter(2);
    else
        distPx = distPx + hypot(testx - prevCenter(1), testy - prevCenter(2));
    end

    prevNorth = newNorth;
    prevWest = newWest;
    prevSouth = newSouth;
    prevEast = newEast;
    prevCenter = [testx testy];
    cx = testx;
    cy = testy;
end

angleDeg = atan2d(cy, cx);
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

function bw = spatialProjectionMask(evtPix, spatialSize)
bw = [];
evtPix = forcePix(evtPix);
if isempty(evtPix) || isempty(spatialSize) || numel(spatialSize) < 2
    return
end
sx = round(spatialSize(1));
sy = round(spatialSize(2));
sz = 1;
if numel(spatialSize) >= 3 && isfinite(spatialSize(3)) && spatialSize(3) > 0
    sz = round(spatialSize(3));
end
if ~(isfinite(sx) && isfinite(sy)) || sx <= 0 || sy <= 0 || sz <= 0
    return
end
try
    [ih, iw, ~] = ind2sub([sx sy sz], evtPix(:));
    lin2 = sub2ind([sx sy], ih, iw);
    bw = false(sx, sy);
    bw(lin2) = true;
catch
    bw = [];
end
end

function n = boundaryPixelCount(bw)
n = NaN;
if isempty(bw) || ~(isnumeric(bw) || islogical(bw))
    return
end
if exist('bwperim', 'file') == 2
    try
        n = nnz(bwperim(bw, 8));
        return
    catch
    end
end
kernel = [1 1 1; 1 0 1; 1 1 1];
neighbors = conv2(double(bw), kernel, 'same');
n = nnz(bw & neighbors < 8);
end

function p = forcePix(p)
if isempty(p)
    p = [];
    return
end
p = unique(double(p(:)));
p = p(isfinite(p) & p > 0);
end
