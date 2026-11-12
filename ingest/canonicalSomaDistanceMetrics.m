function out = canonicalSomaDistanceMetrics(evtPix, xSpaTempEvent, t0Frame, t1Frame, spatialSize, somaMask, somaCenter)
%CANONICALSOMADISTANCEMETRICS Compute event-to-soma centroid distances.
% Distances are reported in matrix-space pixels:
%   - full event footprint centroid to soma centroid
%   - onset-frame centroid to soma centroid
%   - terminal-frame centroid to soma centroid

out = struct( ...
    'soma_centroid_dist_px', NaN, ...
    'soma_start_dist_px', NaN, ...
    'soma_end_dist_px', NaN);

if nargin < 6
    somaMask = [];
end
if nargin < 7
    somaCenter = [];
end

somaCtr = chooseSomaCenter(somaMask, somaCenter);
if any(~isfinite(somaCtr(1:2)))
    return
end

spatialSize = normalizeSpatialSize(spatialSize);
evtCtr = pixelCentroid(evtPix, spatialSize(1:3));
if all(isfinite(evtCtr(1:2)))
    out.soma_centroid_dist_px = centroidDistance(evtCtr, somaCtr);
end

startPix = eventPixelsAtAbsoluteFrame(xSpaTempEvent, t0Frame, spatialSize);
startCtr = pixelCentroid(startPix, spatialSize(1:3));
if all(isfinite(startCtr(1:2)))
    out.soma_start_dist_px = centroidDistance(startCtr, somaCtr);
end

endPix = eventPixelsAtAbsoluteFrame(xSpaTempEvent, t1Frame, spatialSize);
endCtr = pixelCentroid(endPix, spatialSize(1:3));
if all(isfinite(endCtr(1:2)))
    out.soma_end_dist_px = centroidDistance(endCtr, somaCtr);
end
end

function d = centroidDistance(a, b)
d = hypot(a(1) - b(1), a(2) - b(2));
end

function ctr = chooseSomaCenter(somaMask, somaCenter)
ctr = [NaN NaN NaN];
if isnumeric(somaCenter) && ~isempty(somaCenter)
    somaCenter = double(somaCenter);
    if isvector(somaCenter)
        ctr(1:min(numel(somaCenter), 3)) = somaCenter(1:min(numel(somaCenter), 3));
    elseif size(somaCenter, 1) >= 1
        ctr(1:min(size(somaCenter, 2), 3)) = somaCenter(1, 1:min(size(somaCenter, 2), 3));
    end
end
if any(~isfinite(ctr(1:2))) && (isnumeric(somaMask) || islogical(somaMask))
    ctr = maskCenterLocal(somaMask);
end
end

function spatialSize = normalizeSpatialSize(spatialSize)
spatialSize = double(spatialSize(:)');
if isempty(spatialSize)
    spatialSize = [NaN NaN NaN NaN];
end
if numel(spatialSize) < 4
    spatialSize(end+1:4) = NaN;
end
end

function ctr = pixelCentroid(pix, spatialSize3)
ctr = [NaN NaN NaN];
if isempty(pix) || numel(spatialSize3) < 3 || any(~isfinite(spatialSize3(1:3)))
    return
end
sx = round(spatialSize3(1));
sy = round(spatialSize3(2));
sz = round(spatialSize3(3));
if any([sx sy sz] <= 0)
    return
end
pix = unique(double(pix(:)));
pix = pix(isfinite(pix) & pix > 0);
if isempty(pix)
    return
end
try
    [r, c, z] = ind2sub([sx sy sz], pix);
    ctr = [mean(r), mean(c), mean(z)];
catch
    ctr = [NaN NaN NaN];
end
end

function pix = eventPixelsAtAbsoluteFrame(xSpaTempEvent, absFrame, spatialSize)
pix = [];
if isempty(xSpaTempEvent) || ~isnumeric(xSpaTempEvent) || ~isfinite(absFrame)
    return
end
spatialSize = normalizeSpatialSize(spatialSize);
sx = round(spatialSize(1));
sy = round(spatialSize(2));
sz = round(spatialSize(3));
st = max(1, round(spatialSize(4)));
if any(~isfinite([sx sy sz st])) || any([sx sy sz st] <= 0)
    return
end
try
    [ih, iw, il, it] = ind2sub([sx, sy, sz, st], double(xSpaTempEvent(:)));
    keep = (it == round(absFrame));
    if any(keep)
        pix = unique(sub2ind([sx sy sz], ih(keep), iw(keep), il(keep)));
    end
catch
    pix = [];
end
end

function ctr = maskCenterLocal(mask)
ctr = [NaN NaN NaN];
if ~(isnumeric(mask) || islogical(mask))
    return
end
[r, c, z] = ind2sub(size(mask), find(mask > 0));
if isempty(r)
    return
end
if isempty(z)
    z = ones(size(r));
end
ctr = [mean(r), mean(c), mean(z)];
end
