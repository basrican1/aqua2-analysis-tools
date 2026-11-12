function xSpaTempRefined = reconstruct_xspatemp_refined(delayMapRefinedSec, domainMask, t0Frame, t1Frame, spatialSize, varargin)
%RECONSTRUCT_XSPATEMP_REFINED Rebuild a cumulative xSpaTemp-like voxel list.
%   The reconstruction uses a smoothed onset-frame field on the working grid
%   and emits cumulative per-frame masks constrained to the event domain.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'delayMapRefinedSec');
addRequired(p, 'domainMask');
addRequired(p, 't0Frame');
addRequired(p, 't1Frame');
addRequired(p, 'spatialSize');
addParameter(p, 'frameDtSec', NaN, @(x) isnumeric(x) && isscalar(x) && (isnan(x) || (isfinite(x) && x > 0)));
addParameter(p, 'smoothSigmaPx', 1.2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
parse(p, delayMapRefinedSec, domainMask, t0Frame, t1Frame, spatialSize, varargin{:});

delayMapRefinedSec = double(delayMapRefinedSec);
domainMask = logical(domainMask);
t0Frame = round(double(t0Frame));
t1Frame = round(double(t1Frame));
spatialSize = double(spatialSize);

if numel(spatialSize) < 4
    error('reconstruct_xspatemp_refined:BadSpatialSize', ...
        'spatialSize must contain [sx sy sz st].');
end
sx = round(spatialSize(1));
sy = round(spatialSize(2));
sz3 = max(1, round(spatialSize(3)));

frameDtSec = double(p.Results.frameDtSec);
if ~isfinite(frameDtSec) || frameDtSec <= 0
    vals = unique(sort(delayMapRefinedSec(isfinite(delayMapRefinedSec))));
    if numel(vals) >= 2
        d = diff(vals);
        d = d(d > 0);
        if ~isempty(d)
            frameDtSec = median(d);
        else
            frameDtSec = 1.0;
        end
    else
        frameDtSec = 1.0;
    end
end

onsetFrameMap = NaN(size(delayMapRefinedSec));
onsetFrameMap(domainMask & isfinite(delayMapRefinedSec)) = ...
    t0Frame + delayMapRefinedSec(domainMask & isfinite(delayMapRefinedSec)) ./ frameDtSec;
if double(p.Results.smoothSigmaPx) > 0
    onsetFrameMap = localSmoothOnsetFrameMap(onsetFrameMap, domainMask, double(p.Results.smoothSigmaPx));
end
onsetFrameMap(domainMask) = max(t0Frame, min(t1Frame, round(onsetFrameMap(domainMask))));

allIdx = zeros(0, 1);
prevMask = false(sx, sy);
for frameNow = t0Frame:t1Frame
    currMask = domainMask & isfinite(onsetFrameMap) & onsetFrameMap <= frameNow;
    currMask = imclose(currMask, strel('disk', 1, 0)) & domainMask;
    currMask = currMask | prevMask;
    prevMask = currMask;
    linXY = find(currMask);
    if isempty(linXY)
        continue
    end
    linearMovie = linXY + (frameNow - 1) .* (sx * sy * sz3);
    allIdx = [allIdx; linearMovie(:)]; %#ok<AGROW>
end
xSpaTempRefined = double(unique(allIdx(:)));
end

function onsetOut = localSmoothOnsetFrameMap(onsetIn, domainMask, sigmaPx)
onsetOut = onsetIn;
valid = domainMask & isfinite(onsetIn);
if ~any(valid(:))
    return
end
w = double(valid);
num = imgaussfilt(localNanToZero(onsetIn) .* w, sigmaPx, 'FilterDomain', 'spatial');
den = imgaussfilt(w, sigmaPx, 'FilterDomain', 'spatial');
keep = domainMask & den > 0;
onsetOut(keep) = num(keep) ./ den(keep);
med = medfilt2(localNanToZero(onsetOut), [3 3], 'symmetric');
blend = domainMask & isfinite(onsetOut);
onsetOut(blend) = 0.5 * onsetOut(blend) + 0.5 * med(blend);
end

function x = localNanToZero(x)
x(~isfinite(x)) = 0;
end
