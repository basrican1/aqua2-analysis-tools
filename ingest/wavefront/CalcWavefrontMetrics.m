function wavefront = CalcWavefrontMetrics(eventGeom, spatialRes, varargin)
%CALCWAVEFRONTMETRICS Derive branch-based wavefront metrics from event geometry.

p = inputParser;
p.StructExpand = false;
addParameter(p, 'minValidFrames', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'continuityWeight', 0.5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'enforceConnectedFrontier', true, @islogical);
addParameter(p, 'maxComponentJumpPx', 2.5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'minBranchPixels', 5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'minBranchPersistence', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'minBranchAreaFracOfPeak', 0.05, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
addParameter(p, 'dominantBranchFutureWindow', 8, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'branchOriginRewindFrames', 4, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'minRewindBranchPixels', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'edgeSlideIterations', 2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'extendAfterPeak', true, @islogical);
parse(p, varargin{:});

wavefront = emptyWavefront();
if isempty(eventGeom) || ~isstruct(eventGeom) || isempty(eventGeom.frames)
    return
end

frames = eventGeom.frames(:);
[onsetOrigin, onsetFrameIdx] = localInitialCentroid(eventGeom);
if any(~isfinite(onsetOrigin))
    return
end

candidateComps = localConnectedFrontierComponents(eventGeom, logical(p.Results.enforceConnectedFrontier));
rawValidIdx = find(~cellfun(@isempty, candidateComps));
[onsetPeakReachPx, onsetPeakFrameIdx, onsetPeakPoint] = localFindPeakFromOrigin(candidateComps, onsetOrigin, rawValidIdx, onsetFrameIdx);
if ~(isfinite(onsetPeakReachPx) && isfinite(onsetPeakFrameIdx))
    return
end

branchLineage = localTracePeakComponentLineage(eventGeom, onsetPeakFrameIdx, onsetPeakPoint, onsetFrameIdx);
branchLineageValid = all(isfinite(branchLineage), 2);

[branchOrigin, branchOriginFrameIdx] = localSelectBranchOrigin(eventGeom, branchLineage, branchLineageValid, ...
    onsetFrameIdx, onsetPeakFrameIdx, double(p.Results.minBranchPixels), double(p.Results.minBranchPersistence), ...
    double(p.Results.minBranchAreaFracOfPeak), double(p.Results.dominantBranchFutureWindow), onsetOrigin);
[branchOrigin, branchOriginFrameIdx] = localRewindBranchOrigin(eventGeom, branchOrigin, branchOriginFrameIdx, ...
    onsetFrameIdx, double(p.Results.branchOriginRewindFrames), double(p.Results.minRewindBranchPixels));

branchValidIdx = rawValidIdx(rawValidIdx >= branchOriginFrameIdx);
[peakReachPx, peakFrameIdx, peakPoint] = localFindPeakFromOrigin(candidateComps, branchOrigin, branchValidIdx, branchOriginFrameIdx);
if ~(isfinite(peakReachPx) && isfinite(peakFrameIdx))
    branchOrigin = onsetOrigin;
    branchOriginFrameIdx = onsetFrameIdx;
    peakReachPx = onsetPeakReachPx;
    peakFrameIdx = onsetPeakFrameIdx;
    peakPoint = onsetPeakPoint;
end

validated = localValidateForwardBranchToPeak(eventGeom, branchOrigin, branchOriginFrameIdx, peakFrameIdx, peakPoint, double(p.Results.maxComponentJumpPx));
peakPoint = localRefinePeakPointOnValidatedBoundary(validated, peakFrameIdx, peakPoint, branchOrigin);
peakReachPx = hypot(peakPoint(1) - branchOrigin(1), peakPoint(2) - branchOrigin(2));
backwardMasks = localTracePeakComponentMasks(eventGeom, peakFrameIdx, peakPoint, branchOriginFrameIdx);
corridor = localIntersectBranchCorridor(validated, backwardMasks, eventGeom);
trackInfo = localBuildWavefrontTracks(eventGeom, corridor, branchOriginFrameIdx, branchOrigin, peakFrameIdx, peakPoint);
track = trackInfo.primaryTrackXY;
track = localRefineTrackByEdgeSliding(eventGeom, corridor, track, branchOriginFrameIdx, peakFrameIdx, double(p.Results.edgeSlideIterations));
if p.Results.extendAfterPeak
    track = localExtendTrackAfterPeak(eventGeom, corridor, track, peakFrameIdx);
end
trackValid = all(isfinite(track), 2);
validIdx = find(trackValid & frames >= frames(branchOriginFrameIdx) & frames <= frames(peakFrameIdx));
if isempty(validIdx)
    validIdx = unique([branchOriginFrameIdx; peakFrameIdx], 'stable');
end

frontierLengthPx = 0;
if numel(validIdx) >= 2
    for i = 1:(numel(validIdx) - 1)
        a = track(validIdx(i), :);
        b = track(validIdx(i + 1), :);
        frontierLengthPx = frontierLengthPx + hypot(b(1) - a(1), b(2) - a(2));
    end
end
netDispPx = hypot(peakPoint(1) - branchOrigin(1), peakPoint(2) - branchOrigin(2));

wavefront.onsetOriginXY = onsetOrigin;
wavefront.onsetOriginFrame = frames(onsetFrameIdx);
wavefront.onsetPeakPointXY = onsetPeakPoint;
wavefront.onsetPeakFrame = frames(onsetPeakFrameIdx);
wavefront.branchOriginXY = branchOrigin;
wavefront.branchOriginFrame = frames(branchOriginFrameIdx);
wavefront.branchPeakPointXY = peakPoint;
wavefront.branchPeakFrame = frames(peakFrameIdx);
wavefront.originXY = branchOrigin;
wavefront.peakPointXY = peakPoint;
wavefront.trackXY = track;
wavefront.trackValid = trackValid;
wavefront.primaryTrackXY = track;
wavefront.primaryTrackValid = trackValid;
wavefront.secondaryTracks = {trackInfo.secondaryTracks{:}};
wavefront.frames = frames;
wavefront.peakFrame = frames(peakFrameIdx);
wavefront.validFrames = numel(find(trackValid));
wavefront.validatedBranch = corridor;
wavefront.wavefront_max_reach_px = peakReachPx;
wavefront.wavefront_max_reach_um = peakReachPx * spatialRes;
wavefront.wavefront_frontier_length_px = frontierLengthPx;
wavefront.wavefront_frontier_length_um = frontierLengthPx * spatialRes;
wavefront.wavefront_net_displacement_px = netDispPx;
wavefront.wavefront_net_displacement_um = netDispPx * spatialRes;
if peakReachPx > 0
    wavefront.wavefront_tortuosity = frontierLengthPx / peakReachPx;
end
end

function out = emptyWavefront()
out = struct( ...
    'originXY', [NaN NaN], ...
    'onsetOriginXY', [NaN NaN], ...
    'onsetOriginFrame', NaN, ...
    'onsetPeakPointXY', [NaN NaN], ...
    'onsetPeakFrame', NaN, ...
    'branchOriginXY', [NaN NaN], ...
    'branchOriginFrame', NaN, ...
    'branchPeakPointXY', [NaN NaN], ...
    'branchPeakFrame', NaN, ...
    'peakPointXY', [NaN NaN], ...
    'trackXY', zeros(0,2), ...
    'trackValid', false(0,1), ...
    'primaryTrackXY', zeros(0,2), ...
    'primaryTrackValid', false(0,1), ...
    'secondaryTracks', {cell(0,1)}, ...
    'frames', zeros(0,1), ...
    'peakFrame', NaN, ...
    'validFrames', 0, ...
    'validatedBranch', struct('comps', {cell(0,1)}, 'repXY', zeros(0,2), 'valid', false(0,1), 'activeMasks', {cell(0,1)}), ...
    'wavefront_max_reach_px', NaN, ...
    'wavefront_max_reach_um', NaN, ...
    'wavefront_frontier_length_px', NaN, ...
    'wavefront_frontier_length_um', NaN, ...
    'wavefront_net_displacement_px', NaN, ...
    'wavefront_net_displacement_um', NaN, ...
    'wavefront_tortuosity', NaN);
end

function [origin, idx] = localInitialCentroid(eventGeom)
origin = [NaN NaN];
idx = [];
if ~isstruct(eventGeom) || ~isfield(eventGeom, 'centroidXY') || isempty(eventGeom.centroidXY)
    return
end
centroids = double(eventGeom.centroidXY);
idx = find(all(isfinite(centroids), 2), 1, 'first');
if ~isempty(idx)
    origin = centroids(idx, :);
end
end

function [peakReachPx, peakFrameIdx, peakPoint] = localFindPeakFromOrigin(candidateComps, origin, validIdx, minFrameIdx)
peakReachPx = -inf;
peakFrameIdx = NaN;
peakPoint = [NaN NaN];
for idx = 1:numel(validIdx)
    k = validIdx(idx);
    if k < minFrameIdx
        continue
    end
    pts = localFlattenComponents(candidateComps{k});
    if isempty(pts)
        continue
    end
    radial = hypot(pts(:,1) - origin(1), pts(:,2) - origin(2));
    [rmax, imax] = max(radial);
    if isfinite(rmax) && rmax > peakReachPx
        peakReachPx = rmax;
        peakFrameIdx = k;
        peakPoint = pts(imax, :);
    end
end
end

function track = localBacktrace(eventGeom, validated, peakFrameIdx, peakPoint, startFrameIdx)
nFrames = numel(eventGeom.frames);
track = NaN(nFrames, 2);
track(peakFrameIdx, :) = peakPoint;
for k = (peakFrameIdx - 1):-1:startFrameIdx
    nextPt = track(k + 1, :);
    if any(~isfinite(nextPt))
        break
    end
    [pt, ok] = localSelectBackwardBranchPoint(eventGeom, validated, k, nextPt);
    if ~ok
        prevMask = logical(eventGeom.frameMasks{k});
        [pt, ok] = localSelectBackwardFootprintPoint(prevMask, nextPt);
    end
    if ~ok
        break
    end
    track(k, :) = pt;
end
end

function track = localForwardTrackToPeak(eventGeom, validated, startFrameIdx, startPoint, peakFrameIdx, peakPoint)
nFrames = numel(eventGeom.frames);
track = NaN(nFrames, 2);
track(startFrameIdx, :) = startPoint;
track(peakFrameIdx, :) = peakPoint;
if peakFrameIdx <= startFrameIdx + 1
    return
end

cands = cell(nFrames, 1);
costs = cell(nFrames, 1);
prevIdx = cell(nFrames, 1);
cands{startFrameIdx} = struct('pts', {startPoint}, 'rep', startPoint);
costs{startFrameIdx} = 0;
prevIdx{startFrameIdx} = NaN;

for k = (startFrameIdx + 1):peakFrameIdx
    if k == peakFrameIdx
        cand = struct('pts', {peakPoint}, 'rep', peakPoint);
    else
        cand = localForwardCandidateComponents(eventGeom, validated, k);
    end
    if isempty(cand)
        break
    end
    cands{k} = cand;
    prevCand = cands{k - 1};
    prevCost = costs{k - 1};
    if isempty(prevCand) || isempty(prevCost)
        break
    end
    prevRep = vertcat(prevCand.rep);
    currRep = vertcat(cand.rep);
    thisCost = inf(size(currRep,1), 1);
    thisPrev = NaN(size(currRep,1), 1);
    for j = 1:size(currRep, 1)
        jump = hypot(prevRep(:,1) - currRep(j,1), prevRep(:,2) - currRep(j,2));
        total = prevCost + jump;
        [best, ibest] = min(total);
        thisCost(j) = best;
        thisPrev(j) = ibest;
    end
    costs{k} = thisCost;
    prevIdx{k} = thisPrev;
end

if isempty(costs{peakFrameIdx})
    return
end
[~, idx] = min(costs{peakFrameIdx});
for k = peakFrameIdx:-1:(startFrameIdx + 1)
    if isempty(cands{k}) || ~isfinite(idx) || idx < 1 || idx > numel(cands{k})
        break
    end
    if k == peakFrameIdx
        track(k, :) = peakPoint;
    else
        prevPt = track(max(startFrameIdx, k - 1), :);
        if any(~isfinite(prevPt)) && k < peakFrameIdx && ~isempty(prevIdx{k + 1})
            prevPt = cands{k}(idx).rep;
        end
        nextPt = track(min(peakFrameIdx, k + 1), :);
        if any(~isfinite(nextPt))
            nextPt = peakPoint;
        end
        track(k, :) = localRepresentativePointForComponent(cands{k}(idx).pts, prevPt, nextPt);
    end
    idx = prevIdx{k}(idx);
end
track(startFrameIdx, :) = startPoint;
end

function out = localBuildWavefrontTracks(eventGeom, corridor, startFrameIdx, startPoint, peakFrameIdx, peakPoint)
nFrames = numel(eventGeom.frames);
nodesByFrame = cell(nFrames, 1);
for k = startFrameIdx:peakFrameIdx
    mask = localPreferredBackwardMask(eventGeom, corridor, k);
    nodesByFrame{k} = localForwardCandidateComponents(eventGeom, corridor, k);
    if isempty(nodesByFrame{k}) && any(mask(:))
        nodesByFrame{k} = localComponentsFromMask(mask);
    end
end

out = struct('primaryTrackXY', NaN(nFrames,2), 'secondaryTracks', {cell(0,1)});
if isempty(nodesByFrame{startFrameIdx}) || isempty(nodesByFrame{peakFrameIdx})
    out.primaryTrackXY = localForwardTrackToPeak(eventGeom, corridor, startFrameIdx, startPoint, peakFrameIdx, peakPoint);
    return
end

startNode = localNearestComponentIndex(nodesByFrame{startFrameIdx}, startPoint);
peakNode = localContainingOrNearestComponentIndex(nodesByFrame{peakFrameIdx}, peakPoint);

costs = cell(nFrames, 1);
prevIdx = cell(nFrames, 1);
for k = startFrameIdx:peakFrameIdx
    n = numel(nodesByFrame{k});
    costs{k} = inf(n,1);
    prevIdx{k} = NaN(n,1);
end
costs{startFrameIdx}(startNode) = 0;

for k = startFrameIdx:(peakFrameIdx - 1)
    if isempty(nodesByFrame{k}) || isempty(nodesByFrame{k+1})
        continue
    end
    for i = 1:numel(nodesByFrame{k})
        if ~isfinite(costs{k}(i))
            continue
        end
        for j = 1:numel(nodesByFrame{k+1})
            if ~localComponentsLinked(nodesByFrame{k}(i), nodesByFrame{k+1}(j))
                continue
            end
            edgeCost = hypot(nodesByFrame{k}(i).rep(1) - nodesByFrame{k+1}(j).rep(1), ...
                nodesByFrame{k}(i).rep(2) - nodesByFrame{k+1}(j).rep(2));
            edgeCost = edgeCost + nodesByFrame{k+1}(j).nodePenalty;
            total = costs{k}(i) + edgeCost;
            if total < costs{k+1}(j)
                costs{k+1}(j) = total;
                prevIdx{k+1}(j) = i;
            end
        end
    end
end

primaryNodes = localBacktrackNodePath(prevIdx, costs, startFrameIdx, peakFrameIdx, peakNode);
if isempty(primaryNodes)
    out.primaryTrackXY = localForwardTrackToPeak(eventGeom, corridor, startFrameIdx, startPoint, peakFrameIdx, peakPoint);
    return
end
out.primaryTrackXY = localNodePathToTrack(nodesByFrame, primaryNodes, startFrameIdx, peakFrameIdx, startPoint, peakPoint, nFrames);

sinks = localReachableSinkNodes(nodesByFrame, prevIdx, costs, peakFrameIdx, peakNode);
secondary = {};
for s = 1:numel(sinks)
    sink = sinks{s};
    nodePath = localBacktrackNodePath(prevIdx, costs, startFrameIdx, sink.frameIdx, sink.nodeIdx);
    if isempty(nodePath)
        continue
    end
    T = localNodePathToTrack(nodesByFrame, nodePath, startFrameIdx, sink.frameIdx, startPoint, nodesByFrame{sink.frameIdx}(sink.nodeIdx).rep, nFrames);
    secondary{end+1,1} = struct('trackXY', T, 'trackValid', all(isfinite(T),2), 'endFrame', eventGeom.frames(sink.frameIdx)); %#ok<AGROW>
end
out.secondaryTracks = secondary;
end

function idx = localNearestComponentIndex(nodes, pt)
idx = 1;
best = inf;
for i = 1:numel(nodes)
    d = hypot(nodes(i).rep(1) - pt(1), nodes(i).rep(2) - pt(2));
    if d < best
        best = d;
        idx = i;
    end
end
end

function idx = localContainingOrNearestComponentIndex(nodes, pt)
idx = localNearestComponentIndex(nodes, pt);
for i = 1:numel(nodes)
    pts = nodes(i).pts;
    if any(pts(:,1) == round(pt(1)) & pts(:,2) == round(pt(2)))
        idx = i;
        return
    end
end
end

function tf = localComponentsLinked(aNode, bNode)
d = localMinComponentDistance(aNode.pts, bNode.pts);
if aNode.kind == "new" && bNode.kind == "new"
    tf = d <= 10;
else
    tf = d <= 15;
end
end

function nodePath = localBacktrackNodePath(prevIdx, costs, startFrameIdx, endFrameIdx, endNodeIdx)
nodePath = [];
if endFrameIdx < startFrameIdx || isempty(costs{endFrameIdx}) || endNodeIdx < 1 || endNodeIdx > numel(costs{endFrameIdx}) || ~isfinite(costs{endFrameIdx}(endNodeIdx))
    return
end
nodePath = NaN(endFrameIdx - startFrameIdx + 1, 1);
nodePath(end) = endNodeIdx;
for k = endFrameIdx:-1:(startFrameIdx + 1)
    prev = prevIdx{k}(nodePath(k - startFrameIdx + 1));
    if ~isfinite(prev)
        nodePath = [];
        return
    end
    nodePath(k - startFrameIdx) = prev;
end
end

function track = localNodePathToTrack(nodesByFrame, nodePath, startFrameIdx, endFrameIdx, startPoint, endPoint, nFrames)
track = NaN(nFrames, 2);
track(startFrameIdx, :) = startPoint;
for k = startFrameIdx:endFrameIdx
    idx = nodePath(k - startFrameIdx + 1);
    if ~isfinite(idx)
        continue
    end
    if k == startFrameIdx
        track(k,:) = startPoint;
    elseif k == endFrameIdx
        track(k,:) = endPoint;
    else
        prevPt = track(k - 1, :);
        nextPt = nodesByFrame{k + 1}(nodePath(k - startFrameIdx + 2)).rep;
        track(k,:) = localRepresentativePointForComponent(nodesByFrame{k}(idx).pts, prevPt, nextPt);
    end
end
end

function sinks = localReachableSinkNodes(nodesByFrame, prevIdx, costs, peakFrameIdx, peakNode)
sinks = {};
for k = 1:numel(nodesByFrame)
    if isempty(nodesByFrame{k}) || isempty(costs{k})
        continue
    end
    for i = 1:numel(nodesByFrame{k})
        if ~isfinite(costs{k}(i))
            continue
        end
        if k == peakFrameIdx && i == peakNode
            continue
        end
        if k == numel(nodesByFrame) || isempty(nodesByFrame{k+1})
            sinks{end+1,1} = struct('frameIdx', k, 'nodeIdx', i, 'cost', costs{k}(i)); %#ok<AGROW>
        end
    end
end
if isempty(sinks)
    return
end
[~, order] = sort(cellfun(@(s) s.cost, sinks), 'ascend');
sinks = sinks(order(1:min(3, numel(order))));
end

function track = localRefineTrackByEdgeSliding(eventGeom, corridor, track, startFrameIdx, peakFrameIdx, nIter)
if nIter <= 0
    return
end
for iter = 1:nIter
    for k = (startFrameIdx + 1):(peakFrameIdx - 1)
        prevPt = track(k - 1, :);
        currPt = track(k, :);
        nextPt = track(k + 1, :);
        if any(~isfinite(prevPt)) || any(~isfinite(currPt)) || any(~isfinite(nextPt))
            continue
        end
        mask = localPreferredBackwardMask(eventGeom, corridor, k);
        if isempty(mask) || ~any(mask(:))
            continue
        end
        [bestPt, bestCost, ok] = localBestEdgeSlidingPoint(mask, prevPt, nextPt);
        if ~ok
            continue
        end
        currCost = hypot(currPt(1) - prevPt(1), currPt(2) - prevPt(2)) + ...
            hypot(nextPt(1) - currPt(1), nextPt(2) - currPt(2));
        if bestCost + 1e-6 < currCost
            track(k, :) = bestPt;
        end
    end
end
end

function [pt, ok] = localSelectBackwardBranchPoint(eventGeom, validated, k, nextPt)
pt = [NaN NaN];
ok = false;
prefMask = localPreferredBackwardMask(eventGeom, validated, k);
if isempty(prefMask) || ~any(prefMask(:))
    return
end

frontierMask = false(size(prefMask));
if isfield(eventGeom, 'newMasks') && numel(eventGeom.newMasks) >= k && ~isempty(eventGeom.newMasks{k})
    frontierMask = logical(eventGeom.newMasks{k}) & prefMask;
end
if any(frontierMask(:))
    [pt, ok] = localSelectBackwardFootprintPoint(frontierMask, nextPt);
    if ok
        return
    end
end

edgeMask = bwperim(prefMask, 8);
if any(edgeMask(:))
    [pt, ok] = localSelectBackwardFootprintPoint(edgeMask, nextPt);
    if ok
        return
    end
end

[pt, ok] = localSelectBackwardFootprintPoint(prefMask, nextPt);
end

function [pt, cost, ok] = localBestEdgeSlidingPoint(mask, prevPt, nextPt)
pt = [NaN NaN];
cost = inf;
ok = false;
if isempty(mask) || ~any(mask(:))
    return
end
edgeMask = bwperim(mask, 8);
[r, c] = find(edgeMask);
if isempty(r)
    [r, c] = find(mask);
end
if isempty(r)
    return
end
cand = [c, r];
costs = hypot(cand(:,1) - prevPt(1), cand(:,2) - prevPt(2)) + ...
    hypot(nextPt(1) - cand(:,1), nextPt(2) - cand(:,2));
[cost, ibest] = min(costs);
pt = cand(ibest, :);
ok = all(isfinite(pt));
end

function track = localExtendTrackAfterPeak(eventGeom, corridor, track, peakFrameIdx)
nFrames = numel(eventGeom.frames);
if peakFrameIdx >= nFrames || any(~isfinite(track(peakFrameIdx, :)))
    return
end
currPt = track(peakFrameIdx, :);
for k = (peakFrameIdx + 1):nFrames
    mask = localPreferredBackwardMask(eventGeom, corridor, k);
    if isempty(mask) || ~any(mask(:))
        mask = logical(eventGeom.frameMasks{k});
    end
    if isempty(mask) || ~any(mask(:))
        break
    end
    [pt, ok] = localSelectBackwardFootprintPoint(mask, currPt);
    if ~ok
        break
    end
    track(k, :) = pt;
    currPt = pt;
end
end

function cand = localForwardCandidateComponents(eventGeom, validated, k)
cand = struct('pts', {}, 'rep', {}, 'kind', {}, 'nodePenalty', {});
prefMask = localPreferredBackwardMask(eventGeom, validated, k);
if isempty(prefMask) || ~any(prefMask(:))
    return
end

frontierMask = false(size(prefMask));
if isfield(eventGeom, 'newMasks') && numel(eventGeom.newMasks) >= k && ~isempty(eventGeom.newMasks{k})
    frontierMask = logical(eventGeom.newMasks{k}) & prefMask;
end
edgeMask = bwperim(prefMask, 8);
cand = localAppendComponentNodes(cand, frontierMask, "new", 0.0);
cand = localAppendComponentNodes(cand, edgeMask & ~frontierMask, "edge", 8.0);
if isempty(cand)
    cand = localAppendComponentNodes(cand, prefMask, "bulk", 12.0);
end
end

function comps = localComponentsFromMask(mask)
comps = struct('pts', {}, 'rep', {}, 'kind', {}, 'nodePenalty', {});
comps = localAppendComponentNodes(comps, mask, "bulk", 0.0);
end

function comps = localAppendComponentNodes(comps, mask, kind, nodePenalty)
if isempty(mask) || ~any(mask(:))
    return
end
cc = bwconncomp(logical(mask), 8);
for i = 1:cc.NumObjects
    compMask = false(size(mask));
    compMask(cc.PixelIdxList{i}) = true;
    [r, c] = find(compMask);
    pts = [c, r];
    rep = [mean(c, 'omitnan'), mean(r, 'omitnan')];
    comps(end+1,1) = struct('pts', pts, 'rep', rep, 'kind', kind, 'nodePenalty', nodePenalty); %#ok<AGROW>
end
end

function lineage = localTracePeakComponentLineage(eventGeom, peakFrameIdx, peakPoint, firstFrameIdx)
nFrames = numel(eventGeom.frames);
lineage = NaN(nFrames, 2);
activeComps = cell(nFrames, 1);
for k = firstFrameIdx:peakFrameIdx
    activeComps{k} = localConnectedActiveComponents(logical(eventGeom.frameMasks{k}));
end
[pt, comp, ok] = localSelectComponentContainingPoint(activeComps{peakFrameIdx}, peakPoint);
if ~ok
    return
end
lineage(peakFrameIdx, :) = pt;
currComp = comp;
currPt = pt;
for k = (peakFrameIdx - 1):-1:firstFrameIdx
    comps = activeComps{k};
    if isempty(comps)
        continue
    end
    [pt, comp, ok] = localSelectBackwardLineageComponent(comps, currComp, currPt);
    if ~ok
        continue
    end
    lineage(k, :) = pt;
    currComp = comp;
    currPt = pt;
end
end

function activeMasks = localTracePeakComponentMasks(eventGeom, peakFrameIdx, peakPoint, firstFrameIdx)
nFrames = numel(eventGeom.frames);
activeMasks = cell(nFrames, 1);
activeComps = cell(nFrames, 1);
for k = firstFrameIdx:peakFrameIdx
    activeComps{k} = localConnectedActiveComponents(logical(eventGeom.frameMasks{k}));
end
[~, comp, ok] = localSelectComponentContainingPoint(activeComps{peakFrameIdx}, peakPoint);
if ~ok
    return
end
activeMasks{peakFrameIdx} = localMaskFromPoints(size(eventGeom.frameMasks{peakFrameIdx}), comp);
currComp = comp;
currPt = peakPoint;
for k = (peakFrameIdx - 1):-1:firstFrameIdx
    comps = activeComps{k};
    if isempty(comps)
        continue
    end
    [pt, comp, ok] = localSelectBackwardLineageComponent(comps, currComp, currPt);
    if ~ok
        continue
    end
    activeMasks{k} = localMaskFromPoints(size(eventGeom.frameMasks{k}), comp);
    currComp = comp;
    currPt = pt;
end
end

function corridor = localIntersectBranchCorridor(validated, backwardMasks, eventGeom)
nFrames = numel(eventGeom.frames);
activeMasks = cell(nFrames, 1);
repXY = NaN(nFrames, 2);
valid = false(nFrames, 1);
comps = cell(nFrames, 1);
for k = 1:nFrames
    fwd = false(size(eventGeom.frameMasks{k}));
    if isfield(validated, 'activeMasks') && numel(validated.activeMasks) >= k && ~isempty(validated.activeMasks{k})
        fwd = logical(validated.activeMasks{k});
    end
    bwd = false(size(eventGeom.frameMasks{k}));
    if numel(backwardMasks) >= k && ~isempty(backwardMasks{k})
        bwd = logical(backwardMasks{k});
    end
    mask = fwd & bwd;
    if ~any(mask(:))
        if any(fwd(:))
            mask = fwd;
        elseif any(bwd(:))
            mask = bwd;
        end
    end
    activeMasks{k} = mask;
    if any(mask(:))
        pts = localConnectedActiveComponents(mask);
        comps{k} = pts;
        valid(k) = true;
        [r, c] = find(mask);
        repXY(k, :) = [mean(c, 'omitnan'), mean(r, 'omitnan')];
    else
        comps{k} = {};
    end
end
corridor = struct('comps', {comps}, 'activeMasks', {activeMasks}, 'repXY', repXY, 'valid', valid);
end

function [branchOrigin, branchOriginFrameIdx] = localSelectBranchOrigin(eventGeom, track, trackValid, firstFrameIdx, peakFrameIdx, minBranchPixels, minBranchPersistence, minBranchAreaFracOfPeak, dominantBranchFutureWindow, onsetOrigin)
branchOrigin = onsetOrigin;
branchOriginFrameIdx = firstFrameIdx;
trackAreas = zeros(size(trackValid));
for k = firstFrameIdx:peakFrameIdx
    if ~trackValid(k) || any(~isfinite(track(k, :)))
        continue
    end
    trackAreas(k) = localClusterAreaAtPoint(eventGeom.frameMasks{k}, track(k, :));
end
peakArea = max(trackAreas(firstFrameIdx:peakFrameIdx), [], 'omitnan');
if ~isfinite(peakArea) || peakArea <= 0
    peakArea = minBranchPixels;
end
areaThreshold = max(minBranchPixels, ceil(minBranchAreaFracOfPeak * peakArea));
for k = firstFrameIdx:peakFrameIdx
    if ~trackValid(k) || any(~isfinite(track(k, :)))
        continue
    end
    if ~localDominantBranchSupport(trackAreas, trackValid, k, peakFrameIdx, areaThreshold, dominantBranchFutureWindow)
        continue
    end
    if ~localFramePersists(eventGeom, track, trackValid, trackAreas, k, peakFrameIdx, areaThreshold, minBranchPersistence)
        continue
    end
    [ctr, ok] = localClusterCentroidAtPoint(eventGeom.frameMasks{k}, track(k, :));
    if ok
        branchOrigin = ctr;
        branchOriginFrameIdx = k;
        return
    end
end
end

function [branchOrigin, branchOriginFrameIdx] = localRewindBranchOrigin(eventGeom, branchOrigin, branchOriginFrameIdx, onsetFrameIdx, rewindFrames, minRewindBranchPixels)
if rewindFrames <= 0 || branchOriginFrameIdx <= onsetFrameIdx
    return
end

[currComp, ok] = localComponentContainingPoint(eventGeom.frameMasks{branchOriginFrameIdx}, branchOrigin);
if ~ok || size(currComp, 1) < minRewindBranchPixels
    return
end

currRefPt = branchOrigin;
startIdx = max(onsetFrameIdx, branchOriginFrameIdx - rewindFrames);
for k = (branchOriginFrameIdx - 1):-1:startIdx
    comps = localConnectedActiveComponents(logical(eventGeom.frameMasks{k}));
    if isempty(comps)
        break
    end
    [pt, comp, ok] = localSelectBackwardLineageComponent(comps, currComp, currRefPt);
    if ~ok || size(comp, 1) < minRewindBranchPixels
        break
    end
    [ctr, okCtr] = localClusterCentroidAtPoint(eventGeom.frameMasks{k}, pt);
    if ~okCtr
        break
    end
    branchOrigin = ctr;
    branchOriginFrameIdx = k;
    currComp = comp;
    currRefPt = pt;
end
end

function tf = localDominantBranchSupport(trackAreas, trackValid, startIdx, peakFrameIdx, areaThreshold, futureWindow)
tf = false;
endIdx = min(peakFrameIdx, startIdx + futureWindow - 1);
if endIdx < startIdx
    return
end
winValid = trackValid(startIdx:endIdx);
winAreas = trackAreas(startIdx:endIdx);
if ~any(winValid)
    return
end
futureMax = max(winAreas(winValid), [], 'omitnan');
if ~isfinite(futureMax)
    futureMax = 0;
end
tf = futureMax >= areaThreshold;
end

function tf = localFramePersists(eventGeom, track, trackValid, trackAreas, startIdx, peakFrameIdx, areaThreshold, minBranchPersistence)
endIdx = min(peakFrameIdx, startIdx + minBranchPersistence - 1);
tf = endIdx >= startIdx;
for k = startIdx:endIdx
    if ~trackValid(k) || any(~isfinite(track(k, :)))
        tf = false;
        return
    end
    if trackAreas(k) < areaThreshold
        tf = false;
        return
    end
end
end

function area = localClusterAreaAtPoint(frameMask, pt)
area = 0;
[comp, ok] = localComponentContainingPoint(frameMask, pt);
if ok
    area = size(comp, 1);
end
end

function [ctr, ok] = localClusterCentroidAtPoint(frameMask, pt)
ctr = [NaN NaN];
ok = false;
[comp, ok] = localComponentContainingPoint(frameMask, pt);
if ~ok
    return
end
ctr = [mean(comp(:,1), 'omitnan'), mean(comp(:,2), 'omitnan')];
ok = all(isfinite(ctr));
end

function [comp, ok] = localComponentContainingPoint(frameMask, pt)
comp = zeros(0, 2);
ok = false;
if isempty(frameMask) || ~any(frameMask(:)) || any(~isfinite(pt))
    return
end
cc = bwconncomp(logical(frameMask), 8);
L = labelmatrix(cc);
col = round(pt(1));
row = round(pt(2));
if row < 1 || row > size(L, 1) || col < 1 || col > size(L, 2)
    return
end
lab = L(row, col);
if lab == 0
    return
end
pix = cc.PixelIdxList{lab};
[r, c] = ind2sub(size(frameMask), pix);
comp = [c, r];
ok = ~isempty(comp);
end

function candidateComps = localConnectedFrontierComponents(eventGeom, enforceConnectedFrontier)
nFrames = numel(eventGeom.frames);
candidateComps = cell(nFrames, 1);
for k = 1:nFrames
    pts = double(eventGeom.frontierPts{k});
    if isempty(pts)
        candidateComps{k} = {};
        continue
    end
    if ~enforceConnectedFrontier || k == 1
        candidateComps{k} = {pts};
        continue
    end
    prevMask = logical(eventGeom.frameMasks{k - 1});
    newMask = logical(eventGeom.newMasks{k});
    if ~any(prevMask(:)) || ~any(newMask(:))
        candidateComps{k} = {};
        continue
    end
    candidateComps{k} = localConnectedComponentsFromFrontier(newMask, prevMask);
end
end

function validated = localValidateForwardBranchToPeak(eventGeom, origin, startFrameIdx, peakFrameIdx, peakPoint, maxJumpPx)
nFrames = numel(eventGeom.frames);
validatedComps = cell(nFrames, 1);
activeMasks = cell(nFrames, 1);
repXY = NaN(nFrames, 2);
valid = false(nFrames, 1);

activeComps = cell(nFrames, 1);
for k = 1:nFrames
    activeComps{k} = localConnectedActiveComponents(logical(eventGeom.frameMasks{k}));
end

firstIdx = find(~cellfun(@isempty, activeComps) & ((1:nFrames)' >= startFrameIdx), 1, 'first');
if isempty(firstIdx)
    validated = struct('comps', {validatedComps}, 'activeMasks', {activeMasks}, 'repXY', repXY, 'valid', valid);
    return
end

[pt, comp, ok] = localSelectSeedActiveComponent(activeComps{firstIdx}, origin);
if ~ok
    validated = struct('comps', {validatedComps}, 'activeMasks', {activeMasks}, 'repXY', repXY, 'valid', valid);
    return
end
validatedComps{firstIdx} = {comp};
activeMasks{firstIdx} = localMaskFromPoints(size(eventGeom.frameMasks{firstIdx}), comp);
repXY(firstIdx, :) = pt;
valid(firstIdx) = true;
prevComp = comp;
prevPt = pt;

for k = (firstIdx + 1):min(nFrames, peakFrameIdx)
    comps = activeComps{k};
    if isempty(comps)
        continue
    end
    if k == peakFrameIdx
        [pt, comp, ok] = localSelectComponentContainingPoint(comps, peakPoint);
    else
        [pt, comp, ok] = localSelectForwardActiveComponent(comps, prevComp, prevPt, maxJumpPx);
    end
    if ~ok
        continue
    end
    validatedComps{k} = {comp};
    activeMasks{k} = localMaskFromPoints(size(eventGeom.frameMasks{k}), comp);
    repXY(k, :) = pt;
    valid(k) = true;
    prevComp = comp;
    prevPt = pt;
end

validated = struct('comps', {validatedComps}, 'activeMasks', {activeMasks}, 'repXY', repXY, 'valid', valid);
end

function peakPoint = localRefinePeakPointOnValidatedBoundary(validated, peakFrameIdx, peakPoint, origin)
if ~(isstruct(validated) && isfield(validated, 'activeMasks') && numel(validated.activeMasks) >= peakFrameIdx)
    return
end
mask = validated.activeMasks{peakFrameIdx};
if isempty(mask) || ~any(mask(:))
    return
end
edgeMask = bwperim(logical(mask), 8);
[r, c] = find(edgeMask);
if isempty(r)
    [r, c] = find(mask);
end
if isempty(r)
    return
end
pts = [c, r];
radial = hypot(pts(:,1) - origin(1), pts(:,2) - origin(2));
[~, imax] = max(radial);
if ~isempty(imax) && isfinite(imax)
    peakPoint = pts(imax, :);
end
end

function comps = localConnectedComponentsFromFrontier(newMask, prevMask)
comps = {};
cc = bwconncomp(newMask, 8);
if cc.NumObjects == 0
    return
end
touchMask = conv2(double(prevMask), ones(3), 'same') > 0;
for i = 1:cc.NumObjects
    pix = cc.PixelIdxList{i};
    if ~any(touchMask(pix))
        continue
    end
    compMask = false(size(newMask));
    compMask(pix) = true;
    [r, c] = find(compMask);
    comps{end+1,1} = [c, r]; %#ok<AGROW>
end
end

function comps = localConnectedActiveComponents(frameMask)
comps = {};
if ~any(frameMask(:))
    return
end
cc = bwconncomp(frameMask, 8);
for i = 1:cc.NumObjects
    compMask = false(size(frameMask));
    compMask(cc.PixelIdxList{i}) = true;
    [r, c] = find(compMask);
    comps{end+1,1} = [c, r]; %#ok<AGROW>
end
end

function [pt, comp, ok] = localSelectSeedActiveComponent(comps, origin)
pt = [NaN NaN];
comp = zeros(0, 2);
ok = false;
if isempty(comps)
    return
end

bestComp = NaN;
bestJump = inf;
for i = 1:numel(comps)
    pts = double(comps{i});
    jump = hypot(pts(:,1) - origin(1), pts(:,2) - origin(2));
    thisMin = min(jump);
    if thisMin < bestJump
        bestJump = thisMin;
        bestComp = i;
    end
end
if ~isfinite(bestComp)
    return
end
pts = double(comps{bestComp});
jump = hypot(pts(:,1) - origin(1), pts(:,2) - origin(2));
[~, ibest] = min(jump);
pt = pts(ibest, :);
comp = pts;
ok = all(isfinite(pt)) && ~isempty(comp);
end

function [pt, comp, ok] = localSelectForwardActiveComponent(comps, prevComp, prevPt, maxJumpPx)
pt = [NaN NaN];
comp = zeros(0, 2);
ok = false;
if isempty(comps)
    return
end
bestComp = NaN;
bestJump = inf;
prevComp = double(prevComp);
for i = 1:numel(comps)
    pts = double(comps{i});
    thisMin = localMinComponentDistance(pts, prevComp);
    if thisMin < bestJump
        bestJump = thisMin;
        bestComp = i;
    end
end
if ~isfinite(bestComp) || bestJump > maxJumpPx
    return
end
pts = double(comps{bestComp});
jump = hypot(pts(:,1) - prevPt(1), pts(:,2) - prevPt(2));
[~, ibest] = min(jump);
pt = pts(ibest, :);
comp = pts;
ok = all(isfinite(pt)) && ~isempty(comp);
end

function [pt, comp, ok] = localSelectBackwardLineageComponent(comps, nextComp, nextPt)
pt = [NaN NaN];
comp = zeros(0, 2);
ok = false;
if isempty(comps)
    return
end
bestComp = NaN;
bestOverlap = -inf;
bestJump = inf;
for i = 1:numel(comps)
    pts = double(comps{i});
    overlap = localComponentOverlapSize(pts, nextComp);
    jump = hypot(pts(:,1) - nextPt(1), pts(:,2) - nextPt(2));
    thisJump = min(jump);
    if overlap > bestOverlap || (overlap == bestOverlap && thisJump < bestJump)
        bestOverlap = overlap;
        bestComp = i;
        bestJump = thisJump;
    end
end
if ~isfinite(bestComp)
    return
end
pts = double(comps{bestComp});
jump = hypot(pts(:,1) - nextPt(1), pts(:,2) - nextPt(2));
[~, ibest] = min(jump);
pt = pts(ibest, :);
comp = pts;
ok = all(isfinite(pt)) && ~isempty(comp);
end

function mask = localPreferredBackwardMask(eventGeom, validated, k)
mask = false(size(eventGeom.frameMasks{k}));
if isfield(validated, 'activeMasks') && numel(validated.activeMasks) >= k && ~isempty(validated.activeMasks{k})
    mask = logical(validated.activeMasks{k});
end
end

function [pt, ok] = localSelectBackwardFootprintPoint(prevMask, nextPt)
pt = [NaN NaN];
ok = false;
if isempty(prevMask) || ~any(prevMask(:))
    return
end
edgeMask = bwperim(prevMask, 8);
[r, c] = find(edgeMask);
if isempty(r)
    [r, c] = find(prevMask);
end
if isempty(r)
    return
end
cand = [c, r];
jump = hypot(cand(:,1) - nextPt(1), cand(:,2) - nextPt(2));
[~, ibest] = min(jump);
pt = cand(ibest, :);
ok = all(isfinite(pt));
end

function [pt, ok] = localSelectPointTowardPeak(mask, prevPt, peakPt)
pt = [NaN NaN];
ok = false;
if isempty(mask) || ~any(mask(:))
    return
end
[r, c] = find(mask);
if isempty(r)
    return
end
cand = [c, r];
jump = hypot(cand(:,1) - prevPt(1), cand(:,2) - prevPt(2));
toPeak = hypot(cand(:,1) - peakPt(1), cand(:,2) - peakPt(2));
score = 0.6 * jump + 1.0 * toPeak;
[~, ibest] = min(score);
pt = cand(ibest, :);
ok = all(isfinite(pt));
end

function pt = localRepresentativePointForComponent(pts, prevPt, nextPt)
if isempty(pts)
    pt = [NaN NaN];
    return
end
if any(~isfinite(prevPt)) && any(~isfinite(nextPt))
    pt = [mean(pts(:,1), 'omitnan'), mean(pts(:,2), 'omitnan')];
    return
end
if any(~isfinite(prevPt))
    prevPt = nextPt;
end
if any(~isfinite(nextPt))
    nextPt = prevPt;
end
cost = hypot(pts(:,1) - prevPt(1), pts(:,2) - prevPt(2)) + ...
    hypot(nextPt(1) - pts(:,1), nextPt(2) - pts(:,2));
[~, ibest] = min(cost);
pt = pts(ibest, :);
end

function [pt, comp, ok] = localSelectComponentContainingPoint(comps, targetPt)
pt = [NaN NaN];
comp = zeros(0, 2);
ok = false;
if isempty(comps)
    return
end
for i = 1:numel(comps)
    pts = double(comps{i});
    tf = pts(:,1) == targetPt(1) & pts(:,2) == targetPt(2);
    if any(tf)
        pt = targetPt;
        comp = pts;
        ok = true;
        return
    end
end
end

function pts = localFlattenComponents(comps)
if isempty(comps)
    pts = zeros(0, 2);
    return
end
pts = vertcat(comps{:});
end

function d = localMinComponentDistance(aPts, bPts)
if isempty(aPts) || isempty(bPts)
    d = inf;
    return
end
d = inf;
for i = 1:size(aPts, 1)
    dx = bPts(:,1) - aPts(i,1);
    dy = bPts(:,2) - aPts(i,2);
    d = min(d, min(hypot(dx, dy)));
end
end

function n = localComponentOverlapSize(aPts, bPts)
n = 0;
if isempty(aPts) || isempty(bPts)
    return
end
[~, ia] = unique(aPts, 'rows');
aPts = aPts(sort(ia), :);
[~, ib] = unique(bPts, 'rows');
bPts = bPts(sort(ib), :);
n = size(intersect(aPts, bPts, 'rows'), 1);
end

function mask = localMaskFromPoints(sz, pts)
mask = false(sz);
if isempty(pts)
    return
end
pts = round(double(pts));
keep = pts(:,1) >= 1 & pts(:,1) <= sz(2) & pts(:,2) >= 1 & pts(:,2) <= sz(1);
pts = pts(keep, :);
if isempty(pts)
    return
end
idx = sub2ind(sz, pts(:,2), pts(:,1));
mask(idx) = true;
end
