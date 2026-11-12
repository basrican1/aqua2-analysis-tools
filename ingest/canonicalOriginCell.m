function out = canonicalOriginCell(xSpaTempEvent, t0Frame, t1Frame, spatialSize, candidateCells, cellMasks, varargin)
% Determine onset-based origin cell among candidate region.cell.memberIdx cells.
% The chooser uses a denoised union of the first few event frames and only
% accepts an onset winner when there is sufficient evidence; otherwise it
% falls back to the candidate with the strongest full-footprint overlap.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'xSpaTempEvent');
addRequired(p, 't0Frame', @isnumeric);
addRequired(p, 't1Frame', @isnumeric);
addRequired(p, 'spatialSize', @isnumeric);
addRequired(p, 'candidateCells', @isnumeric);
addRequired(p, 'cellMasks');
addParameter(p, 'maxEarlyFrames', 3, @isnumeric);
addParameter(p, 'frameWeights', [1 0.5 0.25], @isnumeric);
addParameter(p, 'minOnsetScore', 0.15, @isnumeric);
addParameter(p, 'minScoreMargin', 0.05, @isnumeric);
addParameter(p, 'minOnsetComponentSize', 4, @isnumeric);
parse(p, xSpaTempEvent, t0Frame, t1Frame, spatialSize, candidateCells, cellMasks, varargin{:});

candidateCells = double(candidateCells(:));
candidateCells = candidateCells(isfinite(candidateCells) & candidateCells > 0);
out = struct( ...
    'originCell', NaN, ...
    'originMode', "none", ...
    'scores', NaN(size(candidateCells)), ...
    'overlapScores', NaN(size(candidateCells)), ...
    'centroidDistPx', NaN(size(candidateCells)), ...
    'fullOverlapScores', NaN(size(candidateCells)), ...
    'rankedCells', candidateCells, ...
    'nEarlyFramesUsed', 0, ...
    'scoreMargin', NaN);

if isempty(candidateCells) || isempty(xSpaTempEvent) || ~isfinite(t0Frame) || ~isfinite(t1Frame)
    return
end

nEvtFrames = max(1, round(t1Frame - t0Frame + 1));
nEarly = min([p.Results.maxEarlyFrames, nEvtFrames, numel(p.Results.frameWeights)]);
if nEarly <= 0
    return
end
w = double(p.Results.frameWeights(:));
w = w(1:nEarly);
w = w / sum(w);
out.nEarlyFramesUsed = nEarly;

[framePixList, earlyCentroid, onsetPix] = onsetFramePixels(xSpaTempEvent, round(t0Frame), spatialSize, ...
    nEarly, w, p.Results.minOnsetComponentSize, cellMasks, candidateCells);
fullEvtPix = allEventPixels(xSpaTempEvent, spatialSize);
if isempty(onsetPix) && isempty(fullEvtPix)
    return
end

[score, ovScore, distScore, fullOvScore] = scoreCandidates(candidateCells, cellMasks, ...
    framePixList, onsetPix, fullEvtPix, earlyCentroid, w);
out.scores = score;
out.overlapScores = ovScore;
out.centroidDistPx = distScore;
out.fullOverlapScores = fullOvScore;

if all(~isfinite(score))
    return
end

[sortedScore, ord] = sort(score, 'descend', 'MissingPlacement', 'last');
sortedCells = candidateCells(ord);
out.rankedCells = sortedCells;
bestScore = sortedScore(1);
secondScore = NaN;
if numel(sortedScore) >= 2
    secondScore = sortedScore(2);
end
if isfinite(bestScore) && isfinite(secondScore)
    out.scoreMargin = bestScore - secondScore;
elseif isfinite(bestScore)
    out.scoreMargin = inf;
end

bestIdx = ord(1);
if isfinite(bestScore) && isfinite(ovScore(bestIdx)) && ...
        ovScore(bestIdx) >= p.Results.minOnsetScore && ...
        (isinf(out.scoreMargin) || out.scoreMargin >= p.Results.minScoreMargin)
    out.originCell = candidateCells(bestIdx);
    out.originMode = "onset";
    return
end

fallbackIdx = chooseFullOverlapFallback(fullOvScore, distScore, candidateCells);
if ~isempty(fallbackIdx) && isfinite(fullOvScore(fallbackIdx)) && fullOvScore(fallbackIdx) > 0
    out.originCell = candidateCells(fallbackIdx);
    out.originMode = "full_overlap";
end
end

function idx = chooseFullOverlapFallback(fullOvScore, distScore, candidateCells)
idx = [];
if isempty(fullOvScore) || all(~isfinite(fullOvScore))
    return
end
keys = [ -fullOvScore(:), distScore(:), candidateCells(:) ];
keys(~isfinite(keys(:,1)),1) = inf;
keys(~isfinite(keys(:,2)),2) = inf;
[~, ord] = sortrows(keys, [1 2 3]);
if isempty(ord)
    return
end
idx = ord(1);
end

function [score, ovScore, distScore, fullOvScore] = scoreCandidates(candidateCells, cellMasks, ...
    framePixList, onsetPix, fullEvtPix, earlyCentroid, w)
nCand = numel(candidateCells);
score = NaN(nCand,1);
ovScore = NaN(nCand,1);
distScore = NaN(nCand,1);
fullOvScore = NaN(nCand,1);

for ii = 1:nCand
    cc = candidateCells(ii);
    if ~(iscell(cellMasks) && numel(cellMasks) >= cc)
        continue
    end
    cellMask = cellMasks{cc};
    if ~(isnumeric(cellMask) || islogical(cellMask))
        continue
    end
    cellPix = find(cellMask > 0);
    if isempty(cellPix)
        continue
    end

    ov = 0;
    usedW = 0;
    for ff = 1:numel(framePixList)
        fp = framePixList{ff};
        if isempty(fp)
            continue
        end
        ov = ov + w(ff) * (numel(intersect(fp, cellPix)) / max(1, numel(fp)));
        usedW = usedW + w(ff);
    end
    if usedW > 0
        ov = ov / usedW;
    end
    ovScore(ii) = ov;

    if ~isempty(fullEvtPix)
        fullOvScore(ii) = numel(intersect(fullEvtPix, cellPix)) / max(1, numel(fullEvtPix));
    end

    cCenter = maskCenterLocal(cellMask);
    if all(isfinite(cCenter(1:2))) && all(isfinite(earlyCentroid(1:2)))
        distScore(ii) = hypot(cCenter(1) - earlyCentroid(1), cCenter(2) - earlyCentroid(2));
    end

    score(ii) = ov;
    if isfinite(distScore(ii))
        score(ii) = score(ii) - 1e-3 * distScore(ii);
    end
end
end

function [framePixList, ctr, onsetPix] = onsetFramePixels(xSpaTempEvent, t0Frame, spatialSize, nEarly, w, minComponentSize, cellMasks, candidateCells)
framePixList = cell(nEarly,1);
ctr = [NaN NaN NaN];
onsetPix = [];
if isempty(xSpaTempEvent) || ~isnumeric(xSpaTempEvent) || numel(spatialSize) < 3
    return
end
sx = round(spatialSize(1));
sy = round(spatialSize(2));
sz = round(spatialSize(3));
if ~(isfinite(sx) && isfinite(sy) && isfinite(sz)) || sx <= 0 || sy <= 0 || sz <= 0
    return
end

sumR = 0; sumC = 0; sumZ = 0; sumW = 0;
vox = double(xSpaTempEvent(:));
for ff = 1:nEarly
    absFrame = t0Frame + ff - 1;
    fp = eventPixelsAtFrame(vox, absFrame, [sx sy sz]);
    framePixList{ff} = fp;
    onsetPix = union(onsetPix, fp);
    if isempty(fp)
        continue
    end
    [r, c, z] = ind2sub([sx sy sz], fp);
    ww = w(ff);
    sumR = sumR + ww * mean(r);
    sumC = sumC + ww * mean(c);
    sumZ = sumZ + ww * mean(z);
    sumW = sumW + ww;
end
if ~isempty(onsetPix)
    onsetPix = denoiseOnsetPixels(onsetPix, [sx sy sz], minComponentSize, cellMasks, candidateCells);
    for ff = 1:nEarly
        framePixList{ff} = intersect(framePixList{ff}, onsetPix);
    end
end
if sumW > 0
    ctr = [sumR/sumW, sumC/sumW, sumZ/sumW];
elseif ~isempty(onsetPix)
    [r, c, z] = ind2sub([sx sy sz], onsetPix);
    ctr = [mean(r), mean(c), mean(z)];
end
end

function evtPix = allEventPixels(xSpaTempEvent, spatialSize)
evtPix = [];
if isempty(xSpaTempEvent) || ~isnumeric(xSpaTempEvent) || numel(spatialSize) < 3
    return
end
sx = round(spatialSize(1));
sy = round(spatialSize(2));
sz = round(spatialSize(3));
st = max(1, round(spatialSize(min(4, numel(spatialSize)))));
if ~(isfinite(sx) && isfinite(sy) && isfinite(sz) && isfinite(st)) || any([sx sy sz st] <= 0)
    return
end
try
    [ih, iw, iz, ~] = ind2sub([sx sy sz st], double(xSpaTempEvent(:)));
    evtPix = unique(sub2ind([sx sy sz], ih, iw, iz));
    evtPix = evtPix(evtPix > 0 & isfinite(evtPix));
catch
    evtPix = [];
end
end

function onsetPix = denoiseOnsetPixels(onsetPix, spatialDims, minComponentSize, cellMasks, candidateCells)
onsetPix = unique(double(onsetPix(:)));
onsetPix = onsetPix(isfinite(onsetPix) & onsetPix > 0);
if isempty(onsetPix) || ~isfinite(minComponentSize) || minComponentSize <= 1
    return
end

bw = false(spatialDims(1), spatialDims(2), spatialDims(3));
bw(onsetPix) = true;

cc = [];
if exist('bwconncomp', 'file') == 2
    try
        if spatialDims(3) > 1
            conn = 26;
        else
            conn = 8;
        end
        cc = bwconncomp(bw, conn);
    catch
        cc = [];
    end
end
if isempty(cc)
    return
end

cellPixAll = [];
for ii = 1:numel(candidateCells)
    ccid = candidateCells(ii);
    if iscell(cellMasks) && numel(cellMasks) >= ccid
        m = cellMasks{ccid};
        if isnumeric(m) || islogical(m)
            cellPixAll = union(cellPixAll, find(m > 0));
        end
    end
end

keepPix = [];
for ii = 1:cc.NumObjects
    pix = cc.PixelIdxList{ii};
    if numel(pix) >= minComponentSize
        keepPix = union(keepPix, pix);
        continue
    end
    if ~isempty(cellPixAll) && any(ismember(pix, cellPixAll))
        keepPix = union(keepPix, pix);
    end
end
if ~isempty(keepPix)
    onsetPix = keepPix(:);
end
end

function fp = eventPixelsAtFrame(xSpaTempEvent, absFrame, spatialDims)
fp = [];
sx = spatialDims(1);
sy = spatialDims(2);
sz = spatialDims(3);
frameStride = sx * sy * sz;
limitStart = frameStride * (absFrame - 1) + 1;
limitEnd = limitStart + frameStride - 1;
v = xSpaTempEvent(xSpaTempEvent >= limitStart & xSpaTempEvent <= limitEnd);
if isempty(v)
    return
end
frameLocal = v - limitStart + 1;
fp = unique(frameLocal(:));
fp = fp(fp > 0 & isfinite(fp));
end

function c = maskCenterLocal(m)
c = [NaN NaN NaN];
if isempty(m) || ~(isnumeric(m) || islogical(m))
    return
end
nd = ndims(m);
if nd <= 2
    [r, col] = find(m > 0);
    if ~isempty(r)
        c = [mean(r), mean(col), 1];
    end
else
    [r, col, z] = ind2sub(size(m), find(m > 0));
    if ~isempty(r)
        c = [mean(r), mean(col), mean(z)];
    end
end
end
