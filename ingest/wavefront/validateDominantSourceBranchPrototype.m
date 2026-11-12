function out = validateDominantSourceBranchPrototype(source, eventID, varargin)
%VALIDATEDOMINANTSOURCEBRANCHPROTOTYPE Validate dominant source branch origins.
%   OUT = validateDominantSourceBranchPrototype(RESULT, EVENTID)
%   OUT = validateDominantSourceBranchPrototype(PATHIN, EVENTID)
%
%   This diagnostic is intentionally independent of the prototype graph path
%   solver. It answers four questions for one event:
%     1. Is the baseline branch origin inside the event footprint?
%     2. Is it inside a real active cluster at that frame?
%     3. Are there multiple persistent source-cluster candidates?
%     4. If so, what origin is implied by tracing the dominant cluster
%        backward within the event footprint?

BWA_setupPaths(fileparts(mfilename('fullpath')));

p = inputParser;
p.StructExpand = false;
addParameter(p, 'minComponentPixels', 10, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'maxComponentsPerFrame', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'minAreaFractionOfLargest', 0.2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
addParameter(p, 'mergePersistenceFrames', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'minTraceOverlapFraction', 0.05, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
addParameter(p, 'maxTraceRepJumpPx', 20, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'Geometry', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'Wavefront', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'SourceLabel', "", @(x) ischar(x) || isstring(x));
parse(p, varargin{:});

sourceLabel = string(p.Results.SourceLabel);
if ~isempty(p.Results.Geometry)
    G = p.Results.Geometry;
    if isempty(p.Results.Wavefront)
        [BWA_AQ2, fallbackLabel] = localResolveSource(source);
        [~, spatialRes] = localEventImagingMeta(BWA_AQ2, eventID);
        W = CalcWavefrontMetrics(G, spatialRes);
        if strlength(sourceLabel) == 0
            sourceLabel = fallbackLabel;
        end
    else
        W = p.Results.Wavefront;
    end
else
    [BWA_AQ2, fallbackLabel] = localResolveSource(source);
    fts = BWA_AQ2.fts1;
    [imageSz, spatialRes] = localEventImagingMeta(BWA_AQ2, eventID);
    geomAll = CalcEventFrameGeometry(fts.loc.xSpaTemp, fts.loc.t0, fts.loc.t1, imageSz, 'eventIDs', eventID);
    G = geomAll{eventID};
    if isempty(p.Results.Wavefront)
        W = CalcWavefrontMetrics(G, spatialRes);
    else
        W = p.Results.Wavefront;
    end
    if strlength(sourceLabel) == 0
        sourceLabel = fallbackLabel;
    end
end
frames = G.frames(:);

originIdx = find(frames == W.branchOriginFrame, 1, 'first');
peakIdx = find(frames == W.branchPeakFrame, 1, 'first');

frameInfo = repmat(localEmptyFrameInfo(), numel(frames), 1);
for k = 1:numel(frames)
    mask = logical(G.frameMasks{k});
    comps = localLimitComponents( ...
        localFrameComponents(mask, p.Results.minComponentPixels), ...
        p.Results.maxComponentsPerFrame, p.Results.minAreaFractionOfLargest);
    frameInfo(k).frame = frames(k);
    frameInfo(k).mask = mask;
    frameInfo(k).components = comps;
    frameInfo(k).nComponents = numel(comps);
end

nComp = reshape([frameInfo.nComponents], [], 1);
mergeIdx = NaN;
for k = max(originIdx + 1, 2):peakIdx
    if nComp(k - 1) > 1 && localHasPersistentMerge(nComp, k, p.Results.mergePersistenceFrames)
        mergeIdx = k;
        break
    end
end
premergeIdx = NaN;
if isfinite(mergeIdx) && mergeIdx > 1
    premergeIdx = mergeIdx - 1;
end

baseline = localValidateAnchor(frameInfo, originIdx, W.branchOriginXY);
candidateTable = table();
dominantCandidateIdx = NaN;
dominantTrace = struct( ...
    'valid', false, ...
    'candidateIdx', NaN, ...
    'traceFrameIdx', NaN, ...
    'traceFrame', NaN, ...
    'traceOriginXY', [NaN NaN], ...
    'traceLengthFrames', 0, ...
    'traceFrameIdxList', zeros(0,1), ...
    'traceRepXY', zeros(0,2));

if isfinite(premergeIdx) && ~isempty(frameInfo(premergeIdx).components)
    comps = frameInfo(premergeIdx).components;
    nCand = numel(comps);
    area = NaN(nCand, 1);
    repX = NaN(nCand, 1);
    repY = NaN(nCand, 1);
    persistence = NaN(nCand, 1);
    earliestFrame = NaN(nCand, 1);
    containsBaseline = false(nCand, 1);
    traceOriginX = NaN(nCand, 1);
    traceOriginY = NaN(nCand, 1);

    for i = 1:nCand
        area(i) = comps(i).area;
        repX(i) = comps(i).rep(1);
        repY(i) = comps(i).rep(2);
        containsBaseline(i) = localComponentContainsPoint(comps(i), W.branchOriginXY);
        trace = localTraceComponentBackward(frameInfo, premergeIdx, i, ...
            p.Results.minTraceOverlapFraction, p.Results.maxTraceRepJumpPx);
        persistence(i) = trace.nFrames;
        earliestFrame(i) = trace.originFrame;
        traceOriginX(i) = trace.originXY(1);
        traceOriginY(i) = trace.originXY(2);
        if ~dominantTrace.valid || area(i) > area(dominantCandidateIdx)
            dominantCandidateIdx = i;
            dominantTrace = trace;
            dominantTrace.valid = trace.valid;
            dominantTrace.candidateIdx = i;
        end
    end

    candidateTable = table((1:nCand)', area, repX, repY, persistence, earliestFrame, ...
        containsBaseline, traceOriginX, traceOriginY, ...
        'VariableNames', {'candidate_idx','area_px','rep_x','rep_y','trace_frames','trace_origin_frame', ...
        'contains_baseline_origin','trace_origin_x','trace_origin_y'});
end

out = struct();
out.sourceLabel = sourceLabel;
out.eventID = eventID;
out.geometry = G;
out.wavefront = W;
out.frames = frames;
out.originFrameIdx = originIdx;
out.peakFrameIdx = peakIdx;
out.baselineOrigin = baseline;
out.nComponentsByFrame = nComp;
out.premergeFrameIdx = premergeIdx;
out.premergeFrame = localFrameValue(frames, premergeIdx);
out.mergeFrameIdx = mergeIdx;
out.mergeFrame = localFrameValue(frames, mergeIdx);
out.candidateTable = candidateTable;
out.dominantCandidateIdx = dominantCandidateIdx;
out.dominantTrace = dominantTrace;
out.frameInfo = frameInfo;
end

function info = localEmptyFrameInfo()
info = struct('frame', NaN, 'mask', false(0,0), 'components', struct('mask', {}, 'pts', {}, 'centroid', {}, 'rep', {}, 'area', {}), 'nComponents', 0);
end

function out = localValidateAnchor(frameInfo, frameIdx, pt)
out = struct( ...
    'frameIdx', frameIdx, ...
    'frame', NaN, ...
    'xy', pt, ...
    'insideFootprint', false, ...
    'insideRawComponent', false, ...
    'rawComponentArea', NaN, ...
    'insideTrackedComponent', false, ...
    'trackedComponentIdx', NaN, ...
    'trackedComponentArea', NaN, ...
    'snappedXY', [NaN NaN]);
if isempty(frameIdx) || ~isfinite(frameIdx) || frameIdx < 1 || frameIdx > numel(frameInfo)
    return
end
mask = frameInfo(frameIdx).mask;
pt0 = localClampPoint(pt, size(mask, 2), size(mask, 1));
out.frame = frameInfo(frameIdx).frame;
out.xy = pt0;
out.snappedXY = localNearestMaskPixel(mask, pt0);
if isempty(pt0) || isempty(mask)
    return
end
out.insideFootprint = mask(pt0(2), pt0(1));
if out.insideFootprint
    cc = bwconncomp(mask, 8);
    pixIdx = sub2ind(size(mask), pt0(2), pt0(1));
    for i = 1:cc.NumObjects
        if any(cc.PixelIdxList{i} == pixIdx)
            out.insideRawComponent = true;
            out.rawComponentArea = numel(cc.PixelIdxList{i});
            break
        end
    end
end
comps = frameInfo(frameIdx).components;
for i = 1:numel(comps)
    if localComponentContainsPoint(comps(i), pt0)
        out.insideTrackedComponent = true;
        out.trackedComponentIdx = i;
        out.trackedComponentArea = comps(i).area;
        return
    end
end
end

function trace = localTraceComponentBackward(frameInfo, startIdx, compIdx, minOverlapFraction, maxRepJumpPx)
trace = struct( ...
    'valid', false, ...
    'originFrameIdx', NaN, ...
    'originFrame', NaN, ...
    'originXY', [NaN NaN], ...
    'nFrames', 0, ...
    'frameIdxList', zeros(0,1), ...
    'repXY', zeros(0,2));
if startIdx < 1 || startIdx > numel(frameInfo) || isempty(frameInfo(startIdx).components) || ...
        compIdx < 1 || compIdx > numel(frameInfo(startIdx).components)
    return
end

currComp = frameInfo(startIdx).components(compIdx);
frameIdxList = startIdx;
repXY = currComp.rep;
originComp = currComp;
for k = (startIdx - 1):-1:1
    prevComps = frameInfo(k).components;
    if isempty(prevComps)
        break
    end
    nextIdx = localChooseBackwardMatch(prevComps, originComp, minOverlapFraction, maxRepJumpPx);
    if ~isfinite(nextIdx)
        break
    end
    originComp = prevComps(nextIdx);
    frameIdxList(end+1,1) = k; %#ok<AGROW>
    repXY(end+1,:) = originComp.rep; %#ok<AGROW>
end

frameIdxList = flipud(frameIdxList(:));
repXY = flipud(repXY);
trace.valid = ~isempty(frameIdxList);
trace.originFrameIdx = frameIdxList(1);
trace.originFrame = frameInfo(frameIdxList(1)).frame;
trace.originXY = repXY(1,:);
trace.nFrames = numel(frameIdxList);
trace.frameIdxList = frameIdxList;
trace.repXY = repXY;
end

function idx = localChooseBackwardMatch(prevComps, currComp, minOverlapFraction, maxRepJumpPx)
idx = NaN;
if isempty(prevComps)
    return
end
bestScore = -inf;
for i = 1:numel(prevComps)
    overlap = nnz(prevComps(i).mask & currComp.mask);
    minArea = min(prevComps(i).area, currComp.area);
    overlapFrac = overlap / max(1, minArea);
    repJump = hypot(prevComps(i).rep(1) - currComp.rep(1), prevComps(i).rep(2) - currComp.rep(2));
    if overlapFrac <= 0 && repJump > maxRepJumpPx
        continue
    end
    if overlapFrac < minOverlapFraction && repJump > maxRepJumpPx
        continue
    end
    score = 100 * overlapFrac - repJump;
    if score > bestScore
        bestScore = score;
        idx = i;
    end
end
end

function tf = localComponentContainsPoint(comp, pt)
tf = false;
if isempty(comp) || any(~isfinite(pt))
    return
end
xy = round(pt);
if xy(2) < 1 || xy(2) > size(comp.mask, 1) || xy(1) < 1 || xy(1) > size(comp.mask, 2)
    return
end
tf = logical(comp.mask(xy(2), xy(1)));
end

function tf = localHasPersistentMerge(nComp, startIdx, persistenceFrames)
stopIdx = min(numel(nComp), startIdx + persistenceFrames - 1);
tf = startIdx <= numel(nComp) && all(nComp(startIdx:stopIdx) == 1);
end

function val = localFrameValue(frames, idx)
val = NaN;
if ~isempty(idx) && isfinite(idx)
    val = frames(idx);
end
end

function comps = localFrameComponents(mask, minPixels)
comps = struct('mask', {}, 'pts', {}, 'centroid', {}, 'rep', {}, 'area', {});
if isempty(mask) || ~any(mask(:))
    return
end
cc = bwconncomp(mask, 8);
for i = 1:cc.NumObjects
    if numel(cc.PixelIdxList{i}) < minPixels
        continue
    end
    compMask = false(size(mask));
    compMask(cc.PixelIdxList{i}) = true;
    [r, c] = find(compMask);
    ctr = [mean(c, 'omitnan'), mean(r, 'omitnan')];
    rep = localNearestPointToTarget([c, r], ctr);
    comps(end+1,1) = struct( ... %#ok<AGROW>
        'mask', compMask, ...
        'pts', [c, r], ...
        'centroid', ctr, ...
        'rep', rep, ...
        'area', numel(cc.PixelIdxList{i}));
end
end

function comps = localLimitComponents(comps, maxCount, minAreaFractionOfLargest)
if isempty(comps)
    return
end
[~, order] = sort([comps.area], 'descend');
comps = comps(order);
largestArea = comps(1).area;
minArea = max(1, minAreaFractionOfLargest * largestArea);
keep = [comps.area] >= minArea;
comps = comps(keep);
if isempty(comps)
    return
end
comps = comps(1:min(numel(comps), maxCount));
end

function pt = localNearestMaskPixel(mask, pt)
[r, c] = find(mask);
if isempty(r)
    pt = [NaN NaN];
    return
end
[~, ibest] = min(hypot(c - pt(1), r - pt(2)));
pt = [c(ibest), r(ibest)];
end

function pt = localNearestPointToTarget(pts, target)
pt = [NaN NaN];
if isempty(pts)
    return
end
[~, ibest] = min(hypot(pts(:,1) - target(1), pts(:,2) - target(2)));
pt = pts(ibest, :);
end

function pt = localClampPoint(pt, nCols, nRows)
pt = round(pt);
pt(1) = min(max(pt(1), 1), nCols);
pt(2) = min(max(pt(2), 1), nRows);
end

function [BWA_AQ2, label] = localResolveSource(source)
label = "";
if isstruct(source) && isfield(source, 'derived') && isstruct(source.derived) && ...
        isfield(source.derived, 'context') && isstruct(source.derived.context) && ...
        isfield(source.derived.context, 'BWA_AQ2')
    cand = source.derived.context.BWA_AQ2;
    if isstruct(cand) && isfield(cand, 'fts1') && isfield(cand, 'opts')
        BWA_AQ2 = cand;
        if isfield(source, 'pathIn')
            label = string(source.pathIn);
        end
        return
    end
end
if isstruct(source) && isfield(source, 'pathIn') && strlength(string(source.pathIn)) > 0
    loaded = loadClusterOutputs(char(string(source.pathIn)), 'attachLegacyMeta', true, 'attachEventMeta', true);
    BWA_AQ2 = loaded.BWA_AQ2;
    label = string(source.pathIn);
    return
end
if ischar(source) || isstring(source)
    loaded = loadClusterOutputs(char(string(source)), 'attachLegacyMeta', true, 'attachEventMeta', true);
    BWA_AQ2 = loaded.BWA_AQ2;
    label = string(source);
    return
end
error('validateDominantSourceBranchPrototype:BadSource', ...
    'Expected a result with pathIn, a cluster-output root, or a result carrying derived.context.BWA_AQ2.');
end

function [imageSz, spatialRes] = localEventImagingMeta(BWA_AQ2, eventID)
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
