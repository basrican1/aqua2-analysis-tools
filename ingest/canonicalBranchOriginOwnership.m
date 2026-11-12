function out = canonicalBranchOriginOwnership(BWA_AQ2, eventID, varargin)
%CANONICALBRANCHORIGINOWNERSHIP Lightweight branch-origin ownership for ingest.
% Returns only the canonical branch origin anchor and owning cell, without
% running the full prototype track solver.

BWA_setupPaths(fileparts(mfilename('fullpath')));

p = inputParser;
p.StructExpand = false;
addRequired(p, 'BWA_AQ2', @isstruct);
addRequired(p, 'eventID', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'minBranchPixels', 5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'minBranchPersistence', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'dominantPremergeAreaRatio', 1.5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'candidateCells', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'cellMasks', {}, @(x) iscell(x) || isempty(x));
parse(p, BWA_AQ2, eventID, varargin{:});

out = struct( ...
    'valid', false, ...
    'originXY', [NaN NaN], ...
    'originFrame', NaN, ...
    'ownerCell', NaN, ...
    'ownerMode', "none", ...
    'baseOriginXY', [NaN NaN], ...
    'baseOriginFrame', NaN, ...
    'usedDominantOverride', false);

fts = BWA_AQ2.fts1;
[imageSz, spatialRes] = localEventImagingMetaCanonicalOwnership(BWA_AQ2, eventID);
geomAll = CalcEventFrameGeometry(fts.loc.xSpaTemp, fts.loc.t0, fts.loc.t1, imageSz, 'eventIDs', eventID);
G = geomAll{eventID};
W = CalcWavefrontMetrics(G, spatialRes, ...
    'minValidFrames', 1, ...
    'minBranchPixels', p.Results.minBranchPixels, ...
    'minBranchPersistence', p.Results.minBranchPersistence);
if ~isstruct(W) || any(~isfinite(W.branchOriginXY))
    return
end

out.baseOriginXY = double(W.branchOriginXY);
out.baseOriginFrame = double(W.branchOriginFrame);

W2 = localApplyDominantBranchOriginOverride(BWA_AQ2, eventID, G, W, double(p.Results.dominantPremergeAreaRatio));
W2 = localSnapBranchOriginToFrameMask(G, W2);
if isstruct(W2) && isfield(W2, 'branchOriginXY') && all(isfinite(W2.branchOriginXY))
    out.originXY = double(W2.branchOriginXY);
    out.originFrame = double(W2.branchOriginFrame);
    out.usedDominantOverride = any(abs(out.originXY - out.baseOriginXY) > 0) || out.originFrame ~= out.baseOriginFrame;
else
    out.originXY = out.baseOriginXY;
    out.originFrame = out.baseOriginFrame;
end

ownerCell = localResolveOwnerCellCanonical(out.originXY, p.Results.cellMasks, p.Results.candidateCells);
if ~(isfinite(ownerCell) && ownerCell > 0)
    ownerCell = localResolveOwnerCellCanonical(out.originXY, p.Results.cellMasks, []);
    out.ownerMode = "nearest_any";
else
    out.ownerMode = "candidate";
end
out.ownerCell = ownerCell;
out.valid = all(isfinite(out.originXY)) && isfinite(out.originFrame);
end

function [imageSz, spatialRes] = localEventImagingMetaCanonicalOwnership(BWA_AQ2, eventID)
opts = BWA_AQ2.opts;
imageSz = double(opts.sz(1,:));
spatialRes = double(opts.spatialRes(1));
fts = BWA_AQ2.fts1;
if isfield(fts, 'basic') && isfield(fts.basic, 'fileID') && ~isempty(fts.basic.fileID)
    fileIDs = string(fts.basic.fileID);
elseif isfield(fts, 'fileID') && ~isempty(fts.fileID)
    fileIDs = string(fts.fileID);
else
    fileIDs = strings(numel(fts.loc.t0), 1);
end
rowIDs = unique(fileIDs, 'stable');
if eventID <= numel(fileIDs)
    row = find(rowIDs == fileIDs(eventID), 1, 'first');
    if ~isempty(row)
        if size(opts.sz, 1) >= row
            imageSz = double(opts.sz(row,:));
        end
        if numel(opts.spatialRes) >= row && isfinite(opts.spatialRes(row))
            spatialRes = double(opts.spatialRes(row));
        end
    end
end
if numel(imageSz) < 3
    imageSz(3) = 1;
end
end

function W = localApplyDominantBranchOriginOverride(BWA_AQ2, eventID, G, W, minAreaRatio)
try
    D = validateDominantSourceBranchPrototype(BWA_AQ2, eventID, ...
        'Geometry', G, ...
        'Wavefront', W, ...
        'SourceLabel', "");
catch
    return
end
if ~isfinite(D.premergeFrame) || isempty(D.candidateTable) || ~isfinite(D.dominantCandidateIdx) || ...
        ~isstruct(D.dominantTrace) || ~D.dominantTrace.valid || any(~isfinite(D.dominantTrace.originXY))
    return
end
localRows = find(D.candidateTable.contains_baseline_origin);
if isempty(localRows)
    return
end
localIdx = localRows(1);
if D.dominantCandidateIdx == localIdx
    return
end
dominantArea = double(D.candidateTable.area_px(D.dominantCandidateIdx));
localArea = double(D.candidateTable.area_px(localIdx));
if ~(isfinite(dominantArea) && isfinite(localArea) && localArea > 0)
    return
end
if dominantArea / localArea < minAreaRatio
    return
end
W.branchOriginFrame = D.dominantTrace.originFrame;
W.branchOriginXY = D.dominantTrace.originXY;
if isfield(W, 'originXY')
    W.originXY = W.branchOriginXY;
end
end

function W = localSnapBranchOriginToFrameMask(G, W)
frames = G.frames(:);
if ~isfinite(W.branchOriginFrame) || any(~isfinite(W.branchOriginXY))
    return
end
k0 = find(frames == W.branchOriginFrame, 1, 'first');
if isempty(k0)
    return
end
W.branchOriginXY = localSnapPointToMaskCanonical(logical(G.frameMasks{k0}), W.branchOriginXY);
if isfield(W, 'originXY')
    W.originXY = W.branchOriginXY;
end
end

function pt = localSnapPointToMaskCanonical(mask, pt)
pt = round(double(pt));
if isempty(mask) || ~any(mask(:)) || numel(pt) < 2 || any(~isfinite(pt))
    pt = [NaN NaN];
    return
end
pt(1) = min(max(pt(1), 1), size(mask,2));
pt(2) = min(max(pt(2), 1), size(mask,1));
if mask(pt(2), pt(1))
    return
end
[r, c] = find(mask);
if isempty(r)
    pt = [NaN NaN];
    return
end
[~, ibest] = min(hypot(c - pt(1), r - pt(2)));
pt = [c(ibest), r(ibest)];
end

function ownerCell = localResolveOwnerCellCanonical(pt, cellMasks, candidateCells)
ownerCell = NaN;
if any(~isfinite(pt)) || ~iscell(cellMasks) || isempty(cellMasks)
    return
end
if isempty(candidateCells)
    candidateCells = 1:numel(cellMasks);
end
candidateCells = double(candidateCells(:)');
candidateCells = candidateCells(isfinite(candidateCells) & candidateCells >= 1 & candidateCells <= numel(cellMasks));
row = round(pt(2));
col = round(pt(1));
inside = [];
for ii = 1:numel(candidateCells)
    cc = candidateCells(ii);
    mask = logical(cellMasks{cc});
    if isempty(mask) || row < 1 || col < 1 || row > size(mask,1) || col > size(mask,2)
        continue
    end
    if mask(row, col)
        inside(end+1) = cc; %#ok<AGROW>
    end
end
if isempty(inside)
    ownerCell = localNearestCellByCentroidCanonical(pt, candidateCells, cellMasks);
else
    ownerCell = localNearestCellByCentroidCanonical(pt, inside, cellMasks);
end
end

function ownerCell = localNearestCellByCentroidCanonical(pt, candidateCells, cellMasks)
ownerCell = NaN;
bestD = inf;
for ii = 1:numel(candidateCells)
    cc = candidateCells(ii);
    mask = logical(cellMasks{cc});
    [r, c] = find(mask);
    if isempty(r)
        continue
    end
    ctr = [mean(c), mean(r)];
    d = hypot(ctr(1) - pt(1), ctr(2) - pt(2));
    if d < bestD
        bestD = d;
        ownerCell = cc;
    end
end
end
