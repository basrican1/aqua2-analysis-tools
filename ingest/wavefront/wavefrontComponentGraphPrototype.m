function out = wavefrontComponentGraphPrototype(source, eventID, varargin)
%WAVEFRONTCOMPONENTGRAPHPROTOTYPE Prototype component-graph wavefront tracker.
%   OUT = wavefrontComponentGraphPrototype(RESULT, EVENTID)
%   OUT = wavefrontComponentGraphPrototype(PATHIN, EVENTID)
%   OUT = wavefrontComponentGraphPrototype(BWA_AQ2, EVENTID)

BWA_setupPaths(fileparts(mfilename('fullpath')));

p = inputParser;
p.StructExpand = false;
addParameter(p, 'maxLinkDistPx', 16, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'showFigure', true, @islogical);
addParameter(p, 'figureName', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'showEdges', false, @islogical);
addParameter(p, 'showNodeLabels', false, @islogical);
addParameter(p, 'showReachability', true, @islogical);
addParameter(p, 'showAlternatives', false, @islogical);
addParameter(p, 'progressRewardWeight', 0.75, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'peakDirectionPenaltyWeight', 4.0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'newNodeReward', 6.0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'maxRepresentativeJumpPx', 16, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'requireSupportConnectivity', true, @islogical);
addParameter(p, 'minNewComponentPixels', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'minEdgeComponentPixels', 6, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'tailBacktrackFrames', 5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'refinePixelWalk', 'auto', @(x) (islogical(x) && isscalar(x)) || any(strcmpi(string(x), ["auto","on","off"])));
addParameter(p, 'useDominantPremergeOverride', true, @islogical);
addParameter(p, 'dominantPremergeAreaRatio', 1.5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'maxPeakDelayFrames', 6, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
parse(p, varargin{:});

[BWA_AQ2, sourceLabel, sourceResult] = localResolveSource(source);
fts = BWA_AQ2.fts1;
[imageSz, spatialRes] = localEventImagingMeta(BWA_AQ2, eventID);
geomAll = CalcEventFrameGeometry(fts.loc.xSpaTemp, fts.loc.t0, fts.loc.t1, imageSz, 'eventIDs', eventID);
G = geomAll{eventID};
Wbase = CalcWavefrontMetrics(G, spatialRes);
W = Wbase;
if p.Results.useDominantPremergeOverride
    W = localApplyDominantPremergeOverride(sourceResult, BWA_AQ2, eventID, G, W, ...
        double(p.Results.dominantPremergeAreaRatio));
end
W = localSnapWavefrontAnchorsToFrames(G, W);
W = localRefreshPeakForOriginOnPeakFrame(G, W);

frames = G.frames(:);
[nodesByFrame, startIdx, endIdx, nodePath, costs, reachable, usedStrictConnectivity, track, trackValid] = ...
    localBuildPrototypeTrackState(G, W, ...
        round(double(p.Results.minNewComponentPixels)), ...
        round(double(p.Results.minEdgeComponentPixels)), ...
        double(p.Results.maxLinkDistPx), ...
        double(p.Results.progressRewardWeight), ...
        double(p.Results.peakDirectionPenaltyWeight), ...
        double(p.Results.newNodeReward), ...
        logical(p.Results.requireSupportConnectivity), ...
        double(p.Results.maxRepresentativeJumpPx), ...
        round(double(p.Results.tailBacktrackFrames)));
if ~localTrackReachesPeak(trackValid, endIdx) && double(p.Results.maxPeakDelayFrames) > 0
    Wdelayed = localDelayPeakFrameToConnectedBoundary(G, W, round(double(p.Results.maxPeakDelayFrames)));
    if Wdelayed.branchPeakFrame ~= W.branchPeakFrame || any(Wdelayed.branchPeakPointXY ~= W.branchPeakPointXY)
        W = Wdelayed;
        [nodesByFrame, startIdx, endIdx, nodePath, costs, reachable, usedStrictConnectivity, track, trackValid] = ...
            localBuildPrototypeTrackState(G, W, ...
                round(double(p.Results.minNewComponentPixels)), ...
                round(double(p.Results.minEdgeComponentPixels)), ...
                double(p.Results.maxLinkDistPx), ...
                double(p.Results.progressRewardWeight), ...
                double(p.Results.peakDirectionPenaltyWeight), ...
                double(p.Results.newNodeReward), ...
                logical(p.Results.requireSupportConnectivity), ...
                double(p.Results.maxRepresentativeJumpPx), ...
                round(double(p.Results.tailBacktrackFrames)));
    end
end
if ~isempty(nodePath)
    W = localSyncPeakFromTrack(frames, endIdx, track, trackValid, W);
end
if p.Results.useDominantPremergeOverride && ~localHasUsableTrack(trackValid, startIdx, endIdx)
    W = localSnapWavefrontAnchorsToFrames(G, Wbase);
    W = localRefreshPeakForOriginOnPeakFrame(G, W);
    [nodesByFrame, startIdx, endIdx, nodePath, costs, reachable, usedStrictConnectivity, track, trackValid] = ...
        localBuildPrototypeTrackState(G, W, ...
            round(double(p.Results.minNewComponentPixels)), ...
            round(double(p.Results.minEdgeComponentPixels)), ...
            double(p.Results.maxLinkDistPx), ...
            double(p.Results.progressRewardWeight), ...
            double(p.Results.peakDirectionPenaltyWeight), ...
            double(p.Results.newNodeReward), ...
            logical(p.Results.requireSupportConnectivity), ...
            double(p.Results.maxRepresentativeJumpPx), ...
            round(double(p.Results.tailBacktrackFrames)));
    if ~isempty(nodePath)
        W = localSyncPeakFromTrack(frames, endIdx, track, trackValid, W);
    end
end
[W, track, trackValid, endIdx] = localLockPeakFrameAtFirstArrival(frames, track, trackValid, W, endIdx);
isMulticellEvent = localIsMulticellEvent(sourceResult, BWA_AQ2, eventID);
useRefinedPixelWalk = localResolveRefinePixelWalk(p.Results.refinePixelWalk, isMulticellEvent);
displayPath = [];
if useRefinedPixelWalk
    displayPath = localTrackDisplayPath(track, trackValid, G, startIdx, endIdx);
end
endNode = [];
if ~isempty(endIdx) && endIdx >= 1 && endIdx <= numel(nodesByFrame) && ~isempty(nodesByFrame{endIdx})
    endNode = localContainingOrNearestNode(nodesByFrame{endIdx}, W.branchPeakPointXY);
end
backReachable = localBackwardReachability(nodesByFrame, endIdx, endNode, double(p.Results.maxLinkDistPx), ...
    logical(p.Results.requireSupportConnectivity));

[protoFrontierLengthPx, protoReachPx] = localTrackMetrics(track, trackValid, W.branchOriginXY, W.branchPeakPointXY, frames, W.branchOriginFrame, W.branchPeakFrame);
W.wavefront_frontier_length_px = protoFrontierLengthPx;
W.wavefront_frontier_length_um = protoFrontierLengthPx * spatialRes;
W.wavefront_max_reach_px = protoReachPx;
W.wavefront_max_reach_um = protoReachPx * spatialRes;
W.wavefront_net_displacement_px = protoReachPx;
W.wavefront_net_displacement_um = protoReachPx * spatialRes;
if protoReachPx > 0
    W.wavefront_tortuosity = protoFrontierLengthPx / protoReachPx;
else
    W.wavefront_tortuosity = NaN;
end

secondary = localSecondarySinks(nodesByFrame, costs, reachable, endIdx, endNode);
alternatives = localFrameAlternatives(nodesByFrame, costs, reachable, backReachable, nodePath, startIdx, endIdx);

hFig = [];
hAx = [];
if p.Results.showFigure
    figName = char(string(p.Results.figureName));
    if isempty(strtrim(figName))
        figName = sprintf('WavefrontComponentGraph_%d', eventID);
    end
    hFig = figure('Name', figName, 'Color', 'w');
    hAx = axes('Parent', hFig);
    hold(hAx, 'on');
    axis(hAx, 'equal');
    set(hAx, 'YDir', 'reverse');
    xlabel(hAx, 'x (um)');
    ylabel(hAx, 'y (um)');
    title(hAx, sprintf('Prototype graph path | Event %d | branch origin=%d | branch peak=%d', ...
        eventID, round(W.branchOriginFrame), round(W.branchPeakFrame)));

    bg = false(imageSz(1), imageSz(2));
    for k = 1:numel(G.frameMasks)
        bg = bg | logical(G.frameMasks{k});
    end
    [r, c] = find(bg);
    if ~isempty(r)
        plot(hAx, c * spatialRes, r * spatialRes, '.', 'Color', [0.80 0.80 0.80], 'MarkerSize', 6);
    end

    if p.Results.showEdges
        localPlotEdges(hAx, nodesByFrame, startIdx, endIdx, double(p.Results.maxLinkDistPx), spatialRes);
    end

    for k = startIdx:endIdx
        if isempty(nodesByFrame{k})
            continue
        end
        for i = 1:numel(nodesByFrame{k})
            rep = nodesByFrame{k}(i).rep * spatialRes;
            nodeColor = localNodeColor(nodesByFrame{k}(i), p.Results.showReachability, reachable{k}(i), backReachable{k}(i));
            plot(hAx, rep(1), rep(2), 'o', 'Color', nodeColor, 'MarkerFaceColor', nodeColor, 'MarkerSize', 5);
            if p.Results.showNodeLabels
                text(hAx, rep(1), rep(2), sprintf('%d:%d', frames(k), i), ...
                    'Color', nodeColor, 'FontSize', 7, 'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'bottom');
            end
        end
    end

    if ~isempty(secondary)
        for s = 1:numel(secondary)
            T = localNodePathToTrack(nodesByFrame, secondary{s}.nodePath, startIdx, secondary{s}.endFrameIdx, ...
                W.branchOriginXY, nodesByFrame{secondary{s}.endFrameIdx}(secondary{s}.endNodeIdx).rep, numel(frames), ...
                double(p.Results.maxRepresentativeJumpPx));
            keep = all(isfinite(T), 2);
            if any(keep)
                pts = T(keep, :) * spatialRes;
                plot(hAx, pts(:,1), pts(:,2), '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.25);
            end
        end
    end

    if any(trackValid)
        if ~isempty(displayPath)
            dpts = displayPath * spatialRes;
            plot(hAx, dpts(:,1), dpts(:,2), '-', 'Color', [0.85 0.20 0.20], 'LineWidth', 2);
        end
        pts = track(trackValid, :) * spatialRes;
        plot(hAx, pts(:,1), pts(:,2), 'o', 'Color', [0.85 0.20 0.20], ...
            'MarkerFaceColor', [0.85 0.20 0.20], 'LineWidth', 1.5, 'MarkerSize', 4);
    end
    if p.Results.showAlternatives
        localPlotAlternatives(hAx, nodesByFrame, alternatives, frames, spatialRes);
    end
    plot(hAx, W.branchOriginXY(1) * spatialRes, W.branchOriginXY(2) * spatialRes, 'ks', ...
        'MarkerFaceColor', 'y', 'MarkerSize', 8, 'LineWidth', 1.25);
    plot(hAx, W.branchPeakPointXY(1) * spatialRes, W.branchPeakPointXY(2) * spatialRes, 'rp', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 11, 'LineWidth', 1.25);
end

out = struct();
out.figure = hFig;
out.axes = hAx;
out.sourceLabel = sourceLabel;
out.geometry = G;
out.wavefront = W;
out.frames = frames;
out.nodesByFrame = nodesByFrame;
out.trackXY = track;
out.trackValid = trackValid;
out.supportPathXY = displayPath;
out.displayPathXY = displayPath;
out.secondary = secondary;
out.spatialRes = spatialRes;
out.forwardReachable = reachable;
out.backwardReachable = backReachable;
out.nodePath = nodePath;
out.nodeCosts = costs;
out.alternatives = alternatives;
out.usedStrictConnectivity = usedStrictConnectivity;
out.isMulticellEvent = isMulticellEvent;
out.usedRefinedPixelWalk = useRefinedPixelWalk;
end

function [nodesByFrame, startIdx, endIdx, nodePath, costs, reachable, usedStrictConnectivity, track, trackValid] = ...
    localBuildPrototypeTrackState(G, W, minNewPixels, minEdgePixels, maxLinkDistPx, ...
    progressRewardWeight, peakDirectionPenaltyWeight, newNodeReward, requireSupportConnectivity, ...
    maxRepresentativeJumpPx, tailBacktrackFrames)
frames = G.frames(:);
startIdx = find(frames == W.branchOriginFrame, 1, 'first');
endIdx = find(frames == W.branchPeakFrame, 1, 'first');
nodesByFrame = cell(numel(frames), 1);
nodePath = cell(numel(frames), 1);
costs = cell(numel(frames), 1);
reachable = cell(numel(frames), 1);
usedStrictConnectivity = false;
track = NaN(numel(frames), 2);
trackValid = false(numel(frames), 1);
if isempty(startIdx) || isempty(endIdx) || endIdx < startIdx
    return
end
for k = startIdx:endIdx
    nodesByFrame{k} = localFrontierComponents(G, k, minNewPixels, minEdgePixels);
end
startNode = localNearestNode(nodesByFrame{startIdx}, W.branchOriginXY);
endNode = localContainingOrNearestNode(nodesByFrame{endIdx}, W.branchPeakPointXY);
[nodePath, costs, reachable, usedStrictConnectivity] = localSolveNodePathWithFallback( ...
    nodesByFrame, startIdx, endIdx, startNode, endNode, maxLinkDistPx, ...
    W.branchOriginXY, W.branchPeakPointXY, progressRewardWeight, peakDirectionPenaltyWeight, newNodeReward, ...
    requireSupportConnectivity);
if isempty(nodePath)
    track = localConnectedPixelWalkTrack(G, W, startIdx, endIdx);
    trackValid = all(isfinite(track), 2);
    return
end
nodeTrack = localNodePathToTrack(nodesByFrame, nodePath, startIdx, endIdx, W.branchOriginXY, W.branchPeakPointXY, numel(frames), ...
    maxRepresentativeJumpPx);
nodeTrack = localRefineTrackTail(nodesByFrame, nodePath, startIdx, endIdx, nodeTrack, W.branchPeakPointXY, ...
    maxRepresentativeJumpPx, tailBacktrackFrames);
nodeTrack = localSnapTrackToActiveMasks(G, nodeTrack, startIdx, endIdx);
nodeTrack = localRefineTrackPairwiseSupport(G, nodeTrack, startIdx, endIdx);
nodeTrack = localRefineTrackPairwiseSupport(G, nodeTrack, startIdx, endIdx);
nodeTrack = localEnforceRawTransitionLegality(G, nodeTrack, startIdx, endIdx);
nodeTrackValid = all(isfinite(nodeTrack), 2);

walkTrack = localConnectedPixelWalkTrack(G, W, startIdx, endIdx);
walkTrackValid = all(isfinite(walkTrack), 2);
if localHasUsableTrack(walkTrackValid, startIdx, endIdx) && ...
        localPreferWalkTrack(G, nodeTrack, nodeTrackValid, walkTrack, walkTrackValid, startIdx, endIdx)
    track = walkTrack;
    trackValid = walkTrackValid;
else
    track = nodeTrack;
    trackValid = nodeTrackValid;
end
[track, trackValid] = localDelayTerminalApproach(G, frames, track, trackValid, W);
track = localEnforceRawTransitionLegality(G, track, startIdx, endIdx);
trackValid = all(isfinite(track), 2);
end

function tf = localHasUsableTrack(trackValid, startIdx, endIdx)
tf = false;
if isempty(startIdx) || isempty(endIdx) || isempty(trackValid)
    return
end
idx = find(trackValid(:));
if isempty(idx)
    return
end
tf = any(idx >= startIdx & idx <= endIdx);
end

function tf = localPreferWalkTrack(G, nodeTrack, nodeTrackValid, walkTrack, walkTrackValid, startIdx, endIdx)
tf = false;
walkScore = localTrackPenaltyScore(G, walkTrack, walkTrackValid, startIdx, endIdx);
nodeScore = localTrackPenaltyScore(G, nodeTrack, nodeTrackValid, startIdx, endIdx);
if walkScore < nodeScore - 1e-6
    tf = true;
elseif abs(walkScore - nodeScore) <= 1e-6
    walkCoverage = nnz(walkTrackValid(startIdx:endIdx));
    nodeCoverage = nnz(nodeTrackValid(startIdx:endIdx));
    if walkCoverage > nodeCoverage
        tf = true;
    end
end
end

function score = localTrackPenaltyScore(G, track, trackValid, startIdx, endIdx)
score = inf;
if isempty(track) || isempty(trackValid) || isempty(startIdx) || isempty(endIdx) || endIdx < startIdx
    return
end
idx = startIdx:endIdx;
valid = trackValid(idx);
if ~any(valid)
    return
end
coveragePenalty = 5.0 * (numel(idx) - nnz(valid));
illegalPenalty = 0;
stepPenalty = 0;
reversalPenalty = 0;
prevVec = [NaN NaN];
for k = startIdx:(endIdx - 1)
    if ~(trackValid(k) && trackValid(k + 1))
        continue
    end
    p1 = round(track(k, :));
    p2 = round(track(k + 1, :));
    rawUnion = logical(G.frameMasks{k}) | logical(G.frameMasks{k + 1});
    if ~localSegmentInsideMask(rawUnion, p1, p2)
        illegalPenalty = illegalPenalty + 1000;
    end
    dv = double(p2) - double(p1);
    stepPenalty = stepPenalty + 0.25 * hypot(dv(1), dv(2));
    if all(isfinite(prevVec)) && all(isfinite(dv))
        if dot(prevVec, dv) < 0
            reversalPenalty = reversalPenalty + 25;
        end
    end
    prevVec = dv;
end
score = coveragePenalty + illegalPenalty + stepPenalty + reversalPenalty;
end

function [track, trackValid] = localDelayTerminalApproach(G, frames, track, trackValid, W)
if isempty(track) || isempty(trackValid) || isempty(frames) || ...
        ~isfinite(W.branchPeakFrame) || any(~isfinite(W.branchPeakPointXY))
    return
end
endIdx = find(frames == W.branchPeakFrame, 1, 'first');
if isempty(endIdx) || endIdx < 2 || endIdx >= numel(frames)
    return
end
peakPt = round(W.branchPeakPointXY);
prePt = round(track(endIdx - 1, :));
if any(~isfinite(prePt))
    return
end
maskPrevRaw = logical(G.frameMasks{endIdx - 1}) | logical(G.frameMasks{endIdx});
if localSegmentInsideMask(maskPrevRaw, prePt, peakPt)
    return
end
maskPrev = localRepairInteriorCracks(maskPrevRaw);
nextMask = localRepairInteriorCracks(logical(G.frameMasks{endIdx + 1}));
if ~any(nextMask(:))
    return
end
elbowMask = localRepairInteriorCracks(logical(G.frameMasks{endIdx}));
[elbowPt, latePt, ok] = localTerminalElbowAndDelay(maskPrev, elbowMask, nextMask, prePt, peakPt);
if ~ok
    return
end
track(endIdx, :) = elbowPt;
trackValid(endIdx) = true;
track(endIdx + 1, :) = latePt;
trackValid(endIdx + 1) = true;
end

function [elbowPt, latePt, ok] = localTerminalElbowAndDelay(maskPrev, elbowMask, nextMask, prePt, peakPt)
elbowPt = [NaN NaN];
latePt = [NaN NaN];
ok = false;
if isempty(maskPrev) || isempty(elbowMask) || isempty(nextMask) || ~any(elbowMask(:)) || ~any(nextMask(:))
    return
end
[er, ec] = find(bwperim(elbowMask, 8));
[nr, nc] = find(bwperim(nextMask, 8));
if isempty(er) || isempty(nr)
    [er, ec] = find(elbowMask);
    [nr, nc] = find(nextMask);
end
if isempty(er) || isempty(nr)
    return
end
elbowPts = [ec, er];
latePts = [nc, nr];
lateDist = hypot(latePts(:, 1) - peakPt(1), latePts(:, 2) - peakPt(2));
[~, lateOrd] = sort(lateDist, 'ascend');
latePts = latePts(lateOrd, :);
lateDist = lateDist(lateOrd);
bestCost = inf;
for j = 1:size(latePts, 1)
    lpt = latePts(j, :);
    if lateDist(j) >= hypot(prePt(1) - peakPt(1), prePt(2) - peakPt(2))
        continue
    end
    elbowDist = hypot(elbowPts(:, 1) - lpt(1), elbowPts(:, 2) - lpt(2));
    [~, elbowOrd] = sort(elbowDist, 'ascend');
    for ii = 1:numel(elbowOrd)
        ept = elbowPts(elbowOrd(ii), :);
        if ~localSegmentInsideMask(maskPrev, prePt, ept)
            continue
        end
        if ~localSegmentInsideMask(elbowMask | nextMask, ept, lpt)
            continue
        end
        cost = 3.0 * hypot(lpt(1) - peakPt(1), lpt(2) - peakPt(2));
        cost = cost + 0.4 * hypot(ept(1) - peakPt(1), ept(2) - peakPt(2));
        cost = cost + 0.05 * hypot(ept(1) - prePt(1), ept(2) - prePt(2));
        cost = cost + 0.02 * hypot(lpt(1) - ept(1), lpt(2) - ept(2));
        if cost < bestCost
            bestCost = cost;
            elbowPt = ept;
            latePt = lpt;
            ok = true;
        end
        if ok
            return
        end
    end
end
end

function tf = localTrackReachesPeak(trackValid, endIdx)
tf = false;
if isempty(endIdx) || isempty(trackValid) || endIdx < 1 || endIdx > numel(trackValid)
    return
end
tf = logical(trackValid(endIdx));
end

function [W, track, trackValid, endIdx] = localLockPeakFrameAtFirstArrival(frames, track, trackValid, W, endIdx)
if isempty(endIdx) || endIdx < 1 || endIdx > numel(frames) || isempty(track) || isempty(trackValid) || ...
        ~isfinite(W.branchPeakFrame) || any(~isfinite(W.branchPeakPointXY)) || any(~isfinite(W.branchOriginXY))
    return
end
if endIdx > size(track, 1) || ~trackValid(endIdx) || any(~isfinite(track(endIdx, :)))
    return
end
arrivalRadiusPx = 10;
reachSlackPx = 8;
minTrailingFrames = 3;
peakPt = double(W.branchPeakPointXY);
originPt = double(W.branchOriginXY);
finalReach = hypot(double(track(endIdx, 1)) - originPt(1), double(track(endIdx, 2)) - originPt(2));
bestIdx = NaN;
for k = 1:max(0, endIdx - minTrailingFrames)
    if ~trackValid(k) || any(~isfinite(track(k, :)))
        continue
    end
    peakDist = hypot(double(track(k, 1)) - peakPt(1), double(track(k, 2)) - peakPt(2));
    if peakDist > arrivalRadiusPx
        continue
    end
    reachK = hypot(double(track(k, 1)) - originPt(1), double(track(k, 2)) - originPt(2));
    if reachK < finalReach - reachSlackPx
        continue
    end
    if ~any(trackValid((k + 1):endIdx))
        continue
    end
    bestIdx = k;
    break
end
if ~isfinite(bestIdx) || bestIdx >= endIdx
    return
end
track((bestIdx + 1):end, :) = NaN;
trackValid((bestIdx + 1):end) = false;
W.branchPeakFrame = frames(bestIdx);
W.branchPeakPointXY = track(bestIdx, :);
W.peakPointXY = W.branchPeakPointXY;
endIdx = bestIdx;
end

function W = localApplyDominantPremergeOverride(sourceResult, BWA_AQ2, eventID, G, W, minAreaRatio)
try
    if isempty(sourceResult)
        src = BWA_AQ2;
    else
        src = sourceResult;
        if ~(isfield(src, 'derived') && isfield(src.derived, 'context') && ...
                isfield(src.derived.context, 'BWA_AQ2'))
            if ~isfield(src, 'derived') || ~isstruct(src.derived)
                src.derived = struct();
            end
            if ~isfield(src.derived, 'context') || ~isstruct(src.derived.context)
                src.derived.context = struct();
            end
            src.derived.context.BWA_AQ2 = BWA_AQ2;
        end
    end
    D = validateDominantSourceBranchPrototype(src, eventID, ...
        'Geometry', G, ...
        'Wavefront', W, ...
        'SourceLabel', "");
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
    W.originXY = W.branchOriginXY;
end
end


function W = localSnapWavefrontAnchorsToFrames(G, W)
frames = G.frames(:);
if isfinite(W.branchOriginFrame)
    k0 = find(frames == W.branchOriginFrame, 1, 'first');
    if ~isempty(k0)
        W.branchOriginXY = localSnapPointToMask(logical(G.frameMasks{k0}), W.branchOriginXY);
        W.originXY = W.branchOriginXY;
    end
end
if isfinite(W.branchPeakFrame)
    kp = find(frames == W.branchPeakFrame, 1, 'first');
    if ~isempty(kp)
        W.branchPeakPointXY = localSnapPointToBoundary(logical(G.frameMasks{kp}), W.branchPeakPointXY);
        W.peakPointXY = W.branchPeakPointXY;
    end
end
end

function W = localRefreshPeakForOriginOnPeakFrame(G, W)
frames = G.frames(:);
if ~isfinite(W.branchPeakFrame) || any(~isfinite(W.branchOriginXY))
    return
end
kp = find(frames == W.branchPeakFrame, 1, 'first');
if isempty(kp)
    return
end
mask = logical(G.frameMasks{kp});
if ~any(mask(:))
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
pts = [c, r];
radial = hypot(pts(:,1) - W.branchOriginXY(1), pts(:,2) - W.branchOriginXY(2));
[~, imax] = max(radial);
if isempty(imax) || ~isfinite(imax)
    return
end
W.branchPeakPointXY = pts(imax, :);
W.peakPointXY = W.branchPeakPointXY;
end

function Wout = localDelayPeakFrameToConnectedBoundary(G, W, maxDelayFrames)
Wout = W;
frames = G.frames(:);
startIdx = find(frames == W.branchOriginFrame, 1, 'first');
peakIdx0 = find(frames == W.branchPeakFrame, 1, 'first');
if isempty(startIdx) || isempty(peakIdx0) || peakIdx0 < startIdx || maxDelayFrames < 1 || any(~isfinite(W.branchOriginXY))
    return
end
originPoint = localSnapPointToMask(logical(G.frameMasks{startIdx}), W.branchOriginXY);
peakHint = round(W.branchPeakPointXY);
lastIdx = min(numel(frames), peakIdx0 + round(maxDelayFrames));
for kp = (peakIdx0 + 1):lastIdx
    unionMask = false(size(G.frameMasks{kp}));
    for kk = startIdx:kp
        unionMask = unionMask | logical(G.frameMasks{kk});
    end
    connectedMask = localConnectedRegionFromSeed(unionMask, originPoint);
    currMask = logical(G.frameMasks{kp});
    currConnected = connectedMask & currMask;
    if ~any(currConnected(:))
        continue
    end
    candidatePt = localSnapPointToBoundary(currConnected, peakHint);
    if isempty(candidatePt) || any(~isfinite(candidatePt))
        continue
    end
    if ~localMaskConnected(unionMask, originPoint, candidatePt)
        continue
    end
    Wout.branchPeakFrame = frames(kp);
    Wout.branchPeakPointXY = candidatePt;
    Wout.peakPointXY = candidatePt;
    return
end
end

function W = localSyncPeakFromTrack(frames, endIdx, track, trackValid, W)
if isempty(endIdx) || endIdx < 1 || endIdx > numel(frames) || isempty(track) || isempty(trackValid)
    return
end
if endIdx <= size(track, 1) && endIdx <= numel(trackValid) && trackValid(endIdx) && all(isfinite(track(endIdx,:)))
    if ~isfinite(W.branchPeakFrame)
        W.branchPeakFrame = frames(endIdx);
    end
    if any(~isfinite(W.branchPeakPointXY))
        W.branchPeakPointXY = track(endIdx,:);
    end
    W.peakPointXY = W.branchPeakPointXY;
end
end

function track = localSnapTrackToActiveMasks(G, track, startIdx, endIdx)
for k = startIdx:endIdx
    if k < 1 || k > numel(G.frameMasks) || any(~isfinite(track(k,:)))
        continue
    end
    track(k,:) = localSnapPointToMask(logical(G.frameMasks{k}), track(k,:));
end
end

function track = localRefineTrackPairwiseSupport(G, track, startIdx, endIdx)
for pass = 1:2
    changed = false;
    for k = startIdx:(endIdx - 1)
        if any(~isfinite(track(k,:))) || any(~isfinite(track(k + 1,:)))
            continue
        end
        mask1 = logical(G.frameMasks{k});
        mask2 = logical(G.frameMasks{k + 1});
        unionMask = mask1 | mask2;
        supportMask = localRepairInteriorCracks(unionMask);
        p1 = round(track(k,:));
        p2 = round(track(k + 1,:));
        if localSegmentInsideMask(supportMask, p1, p2)
            continue
        end
        seg = localMaskedShortestPath(supportMask, p1, p2);
        if isempty(seg)
            continue
        end
        [q1, q2, ok] = localFrameTransitionPair(seg, mask1, mask2, p1, p2, supportMask);
        if ~ok
            continue
        end
        if any(q1 ~= p1)
            track(k,:) = q1;
            changed = true;
        end
        if any(q2 ~= p2)
            track(k + 1,:) = q2;
            changed = true;
        end
    end
    if ~changed
        break
    end
end
end

function track = localEnforceRawTransitionLegality(G, track, startIdx, endIdx)
for k = startIdx:(endIdx - 1)
    if any(~isfinite(track(k,:))) || any(~isfinite(track(k + 1,:)))
        continue
    end
    mask1 = logical(G.frameMasks{k});
    mask2 = logical(G.frameMasks{k + 1});
    rawUnion = mask1 | mask2;
    p1 = round(track(k,:));
    p2 = round(track(k + 1,:));
    if localSegmentInsideMask(rawUnion, p1, p2)
        continue
    end
    seg = localMaskedShortestPath(rawUnion, p1, p2);
    if ~isempty(seg)
        [q1, q2, ok] = localFrameTransitionPair(seg, mask1, mask2, p1, p2, rawUnion);
        if ok
            track(k,:) = q1;
            track(k + 1,:) = q2;
            continue
        end
    end
    guidePt = localFutureGuidePoint(track, k + 1, endIdx, p2);
    q2 = localNearestLegalNextPoint(mask1, mask2, rawUnion, p1, guidePt);
    if all(isfinite(q2))
        track(k + 1,:) = q2;
        continue
    end
    track((k + 1):end, :) = NaN;
    return
end
end

function guidePt = localFutureGuidePoint(track, startIdx, endIdx, defaultPt)
guidePt = defaultPt;
if isempty(track) || startIdx > size(track, 1)
    return
end
for k = endIdx:-1:startIdx
    if all(isfinite(track(k,:)))
        guidePt = round(track(k,:));
        return
    end
end
end

function pt = localNearestLegalNextPoint(mask1, mask2, unionMask, p1, targetPt)
pt = [NaN NaN];
if isempty(mask2) || ~any(mask2(:)) || isempty(unionMask) || ~any(unionMask(:))
    return
end
newMask = mask2 & ~mask1;
edgeMask = bwperim(imfill(mask2, 'holes'), 8) & mask2;
tiers = {newMask, edgeMask, mask2};
for t = 1:numel(tiers)
    [r, c] = find(tiers{t});
    if isempty(r)
        continue
    end
    pts = [c, r];
    distTarget = hypot(pts(:,1) - targetPt(1), pts(:,2) - targetPt(2));
    distStep = hypot(pts(:,1) - p1(1), pts(:,2) - p1(2));
    progress = hypot(pts(:,1) - p1(1), pts(:,2) - p1(2));
    [~, ord] = sortrows([distTarget, distStep, -progress], [1 2 3]);
    pts = pts(ord, :);
    for i = 1:size(pts, 1)
        cand = pts(i, :);
        if localSegmentInsideMask(unionMask, p1, cand)
            pt = cand;
            return
        end
    end
end
end

function pt = localSnapPointToMask(mask, pt)
pt = localClampPoint(pt, size(mask,2), size(mask,1));
if any(~isfinite(pt))
    return
end
if mask(round(pt(2)), round(pt(1)))
    pt = round(pt);
    return
end
pt = localNearestMaskPixel(mask, pt);
end

function tf = localSegmentInsideMask(mask, p1, p2)
tf = false;
if isempty(mask) || ~any(mask(:)) || any(~isfinite(p1)) || any(~isfinite(p2))
    return
end
nRows = size(mask, 1);
nCols = size(mask, 2);
p1 = localClampPoint(p1, nCols, nRows);
p2 = localClampPoint(p2, nCols, nRows);
n = max(abs(p2 - p1));
if n < 1
    tf = localPointInMask(mask, p1);
    return
end
t = linspace(0, 1, max(2 * n + 1, 3));
xs = round(p1(1) + (p2(1) - p1(1)) .* t);
ys = round(p1(2) + (p2(2) - p1(2)) .* t);
xs = min(max(xs, 1), nCols);
ys = min(max(ys, 1), nRows);
idx = sub2ind([nRows, nCols], ys(:), xs(:));
tf = all(mask(idx));
end

function [q1, q2, ok] = localFrameTransitionPair(seg, mask1, mask2, p1, p2, unionMask)
q1 = p1;
q2 = p2;
ok = false;
if size(seg, 1) < 2
    return
end
in1 = false(size(seg, 1), 1);
in2 = false(size(seg, 1), 1);
for i = 1:size(seg, 1)
    in1(i) = localPointInMask(mask1, seg(i,:));
    in2(i) = localPointInMask(mask2, seg(i,:));
end
transition = find(in2, 1, 'first');
if isempty(transition)
    return
end
pre = find(in1(1:transition), 1, 'last');
if isempty(pre)
    return
end
i1cand = max(1, pre - 2):pre;
i2cand = transition:min(size(seg, 1), transition + 2);
bestCost = inf;
for a = i1cand
    for b = i2cand
        c1 = seg(a,:);
        c2 = seg(b,:);
        if ~localPointInMask(mask1, c1) || ~localPointInMask(mask2, c2)
            continue
        end
        if ~localSegmentInsideMask(unionMask, c1, c2)
            continue
        end
        cost = hypot(c1(1) - p1(1), c1(2) - p1(2)) + hypot(c2(1) - p2(1), c2(2) - p2(2));
        cost = cost + 0.25 * hypot(c2(1) - c1(1), c2(2) - c1(2));
        if cost < bestCost
            bestCost = cost;
            q1 = c1;
            q2 = c2;
            ok = true;
        end
    end
end
end

function pt = localSnapPointToBoundary(mask, pt)
filledMask = imfill(mask, 'holes');
edgeMask = bwperim(filledMask, 8);
[r, c] = find(edgeMask);
if isempty(r)
    pt = localSnapPointToMask(mask, pt);
    return
end
pt = localNearestPointToTarget([c, r], pt);
end

function [frontierLengthPx, reachPx] = localTrackMetrics(track, trackValid, originXY, peakXY, frames, originFrame, peakFrame)
frontierLengthPx = 0;
reachPx = hypot(peakXY(1) - originXY(1), peakXY(2) - originXY(2));
if isempty(track) || isempty(trackValid)
    return
end
keep = trackValid & frames >= originFrame & frames <= peakFrame;
idx = find(keep);
for i = 1:(numel(idx) - 1)
    a = track(idx(i), :);
    b = track(idx(i + 1), :);
    frontierLengthPx = frontierLengthPx + hypot(b(1) - a(1), b(2) - a(2));
end
end

function comps = localFrontierComponents(G, k, minNewPixels, minEdgePixels)
comps = struct('pts', {}, 'boundaryPts', {}, 'rep', {}, 'kind', {}, 'nodePenalty', {}, 'supportMask', {}, 'area', {});

newMask = false(size(G.frameMasks{k}));
if k <= numel(G.newMasks) && ~isempty(G.newMasks{k})
    newMask = logical(G.newMasks{k});
end
filledMask = imfill(logical(G.frameMasks{k}), 'holes');
edgeMask = bwperim(filledMask, 8);

comps = localAppendComponents(comps, newMask, edgeMask, filledMask, "new", 0.0, minNewPixels);
comps = localAppendComponents(comps, edgeMask & ~newMask, edgeMask, filledMask, "edge", 18.0, minEdgePixels);
end

function comps = localAppendComponents(comps, mask, edgeMask, supportMask, kind, nodePenalty, minPixels)
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
    pts = [c, r];
    boundaryMask = compMask & edgeMask;
    if ~any(boundaryMask(:))
        boundaryMask = bwperim(compMask, 8);
    end
    [rb, cb] = find(boundaryMask);
    boundaryPts = [cb, rb];
    if isempty(boundaryPts)
        boundaryPts = pts;
    end
    rep = [mean(boundaryPts(:,1), 'omitnan'), mean(boundaryPts(:,2), 'omitnan')];
    comps(end+1,1) = struct( ...
        'pts', pts, ...
        'boundaryPts', boundaryPts, ...
        'rep', rep, ...
        'kind', kind, ...
        'nodePenalty', nodePenalty, ...
        'supportMask', supportMask, ...
        'area', numel(cc.PixelIdxList{i})); %#ok<AGROW>
end
end

function [nodePath, costs, reachable, usedStrictConnectivity] = localSolveNodePathWithFallback(nodesByFrame, startIdx, endIdx, startNode, endNode, maxLinkDistPx, originXY, peakPointXY, progressRewardWeight, peakDirectionPenaltyWeight, newNodeReward, requireSupportConnectivity)
[nodePath, costs, reachable] = localSolveNodePath( ...
    nodesByFrame, startIdx, endIdx, startNode, endNode, maxLinkDistPx, peakPointXY, ...
    progressRewardWeight, peakDirectionPenaltyWeight, newNodeReward, true);
usedStrictConnectivity = true;
if isempty(nodePath) && ~requireSupportConnectivity
    [nodePath, costs, reachable] = localSolveNodePath( ...
        nodesByFrame, startIdx, endIdx, startNode, endNode, maxLinkDistPx, peakPointXY, ...
        progressRewardWeight, peakDirectionPenaltyWeight, newNodeReward, false);
    usedStrictConnectivity = false;
end
end

function idx = localNearestNode(nodes, pt)
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

function idx = localContainingOrNearestNode(nodes, pt)
idx = localNearestNode(nodes, pt);
for i = 1:numel(nodes)
    pts = round(nodes(i).pts);
    if any(pts(:,1) == round(pt(1)) & pts(:,2) == round(pt(2)))
        idx = i;
        return
    end
end
end

function [nodePath, costs, reachable] = localSolveNodePath(nodesByFrame, startIdx, endIdx, startNode, endNode, maxLinkDistPx, peakPointXY, progressRewardWeight, peakDirectionPenaltyWeight, newNodeReward, requireSupportConnectivity)
costs = cell(numel(nodesByFrame), 1);
prev = cell(numel(nodesByFrame), 1);
reachable = cell(numel(nodesByFrame), 1);
for k = startIdx:endIdx
    n = numel(nodesByFrame{k});
    costs{k} = inf(n, 1);
    prev{k} = NaN(n, 1);
    reachable{k} = false(n, 1);
end
costs{startIdx}(startNode) = 0;
reachable{startIdx}(startNode) = true;

for k = startIdx:(endIdx - 1)
    for i = 1:numel(nodesByFrame{k})
        if ~reachable{k}(i)
            continue
        end
        for j = 1:numel(nodesByFrame{k + 1})
            if ~localNodesLinked(nodesByFrame{k}(i), nodesByFrame{k + 1}(j), maxLinkDistPx, requireSupportConnectivity)
                continue
            end
            edgeCost = localTransitionCost( ...
                nodesByFrame{k}(i), nodesByFrame{k + 1}(j), peakPointXY, ...
                progressRewardWeight, peakDirectionPenaltyWeight, newNodeReward);
            total = costs{k}(i) + edgeCost;
            if total < costs{k + 1}(j)
                costs{k + 1}(j) = total;
                prev{k + 1}(j) = i;
                reachable{k + 1}(j) = true;
            end
        end
    end
end

function edgeCost = localTransitionCost(a, b, peakPointXY, progressRewardWeight, peakDirectionPenaltyWeight, newNodeReward)
stepDx = b.rep(1) - a.rep(1);
stepDy = b.rep(2) - a.rep(2);
stepDist = hypot(stepDx, stepDy);
edgeCost = stepDist + b.nodePenalty;
if b.kind == "new"
    edgeCost = edgeCost - newNodeReward;
end

distA = hypot(peakPointXY(1) - a.rep(1), peakPointXY(2) - a.rep(2));
distB = hypot(peakPointXY(1) - b.rep(1), peakPointXY(2) - b.rep(2));
progressTowardPeak = distA - distB;
edgeCost = edgeCost - progressRewardWeight * progressTowardPeak;

toPeakDx = peakPointXY(1) - a.rep(1);
toPeakDy = peakPointXY(2) - a.rep(2);
toPeakDist = hypot(toPeakDx, toPeakDy);
if stepDist > 0 && toPeakDist > 0
    cosTheta = (stepDx * toPeakDx + stepDy * toPeakDy) / (stepDist * toPeakDist);
    cosTheta = max(-1, min(1, cosTheta));
    edgeCost = edgeCost + peakDirectionPenaltyWeight * (1 - cosTheta);
end
end

nodePath = [];
if ~reachable{endIdx}(endNode)
    return
end
nodePath = NaN(endIdx - startIdx + 1, 1);
nodePath(end) = endNode;
for k = endIdx:-1:(startIdx + 1)
    nodePath(k - startIdx) = prev{k}(nodePath(k - startIdx + 1));
end
end

function reachable = localBackwardReachability(nodesByFrame, endIdx, endNode, maxLinkDistPx, requireSupportConnectivity)
reachable = cell(numel(nodesByFrame), 1);
for k = 1:numel(nodesByFrame)
    reachable{k} = false(numel(nodesByFrame{k}), 1);
end
if isempty(endIdx) || ~isfinite(endIdx) || endIdx < 1 || endIdx > numel(nodesByFrame) || ...
        isempty(endNode) || ~isfinite(endNode) || endNode < 1 || endNode > numel(reachable{endIdx})
    return
end
reachable{endIdx}(endNode) = true;
for k = (endIdx - 1):-1:1
    if isempty(nodesByFrame{k}) || isempty(nodesByFrame{k + 1})
        continue
    end
    for i = 1:numel(nodesByFrame{k})
        for j = 1:numel(nodesByFrame{k + 1})
            if ~reachable{k + 1}(j)
                continue
            end
            if localNodesLinked(nodesByFrame{k}(i), nodesByFrame{k + 1}(j), maxLinkDistPx, requireSupportConnectivity)
                reachable{k}(i) = true;
                break
            end
        end
    end
end
end

function tf = localNodesLinked(a, b, maxLinkDistPx, requireSupportConnectivity)
[d, aPt, bPt] = localClosestPair(a.boundaryPts, b.boundaryPts);
if a.kind == "new" && b.kind == "new"
    tf = d <= maxLinkDistPx;
else
    tf = d <= maxLinkDistPx * 1.5;
end
if ~tf
    return
end
if ~requireSupportConnectivity
    return
end
unionMask = logical(a.supportMask) | logical(b.supportMask);
tf = localMaskConnected(unionMask, aPt, bPt);
end

function [d, aPt, bPt] = localClosestPair(aPts, bPts)
d = inf;
aPt = [NaN NaN];
bPt = [NaN NaN];
if isempty(aPts) || isempty(bPts)
    return
end
for i = 1:size(aPts, 1)
    dx = bPts(:,1) - aPts(i,1);
    dy = bPts(:,2) - aPts(i,2);
    [dHere, idx] = min(hypot(dx, dy));
    if dHere < d
        d = dHere;
        aPt = aPts(i, :);
        bPt = bPts(idx, :);
    end
end
end

function localPlotEdges(ax, nodesByFrame, startIdx, endIdx, maxLinkDistPx, spatialRes)
for k = startIdx:(endIdx - 1)
    if isempty(nodesByFrame{k}) || isempty(nodesByFrame{k + 1})
        continue
    end
    for i = 1:numel(nodesByFrame{k})
        a = nodesByFrame{k}(i);
        for j = 1:numel(nodesByFrame{k + 1})
            b = nodesByFrame{k + 1}(j);
            if ~localNodesLinked(a, b, maxLinkDistPx, false)
                continue
            end
            plot(ax, [a.rep(1), b.rep(1)] * spatialRes, [a.rep(2), b.rep(2)] * spatialRes, ...
                '-', 'Color', [0.82 0.82 0.82], 'LineWidth', 0.5);
        end
    end
end
end

function color = localNodeColor(node, showReachability, fwd, back)
if node.kind == "new"
    base = [0.16 0.45 0.95];
else
    base = [0.10 0.65 0.30];
end
if ~showReachability
    color = base;
    return
end
if fwd && back
    color = [0.88 0.35 0.15];
elseif fwd
    color = [0.35 0.60 0.95];
elseif back
    color = [0.25 0.70 0.45];
else
    color = [0.70 0.70 0.70];
end
end

function track = localNodePathToTrack(nodesByFrame, nodePath, startIdx, endIdx, startPoint, endPoint, nFrames, maxRepresentativeJumpPx)
track = NaN(nFrames, 2);
if isempty(nodePath)
    return
end
track(startIdx, :) = startPoint;
for k = (startIdx + 1):(endIdx - 1)
    idx = nodePath(k - startIdx + 1);
    prevPt = track(k - 1, :);
    nextRep = nodesByFrame{k + 1}(nodePath(k - startIdx + 2)).rep;
    track(k, :) = localRepresentativePoint(nodesByFrame{k}(idx), prevPt, nextRep, maxRepresentativeJumpPx);
end
track(endIdx, :) = endPoint;
end

function pt = localRepresentativePoint(node, prevPt, nextPt, maxRepresentativeJumpPx)
pts = node.boundaryPts;
dPrev = hypot(pts(:,1) - prevPt(1), pts(:,2) - prevPt(2));
dNext = hypot(nextPt(1) - pts(:,1), nextPt(2) - pts(:,2));
cost = dPrev + 0.75 * dNext + 0.05 * (dPrev .^ 2);
cost = cost + 5.0 * max(dPrev - maxRepresentativeJumpPx, 0);
[~, ibest] = min(cost);
pt = pts(ibest, :);
end

function track = localRefineTrackTail(nodesByFrame, nodePath, startIdx, endIdx, track, peakPoint, maxRepresentativeJumpPx, tailBacktrackFrames)
if isempty(nodePath) || tailBacktrackFrames <= 0
    return
end
track(endIdx, :) = peakPoint;
tailStart = max(startIdx + 1, endIdx - tailBacktrackFrames);
for k = (endIdx - 1):-1:tailStart
    idx = nodePath(k - startIdx + 1);
    node = nodesByFrame{k}(idx);
    nextPt = track(k + 1, :);
    basePt = track(k, :);
    pts = node.boundaryPts;
    dNext = hypot(pts(:,1) - nextPt(1), pts(:,2) - nextPt(2));
    dBase = hypot(pts(:,1) - basePt(1), pts(:,2) - basePt(2));
    cost = 1.25 * dNext + 0.40 * dBase + 0.06 * (dNext .^ 2);
    cost = cost + 8.0 * max(dNext - maxRepresentativeJumpPx, 0);
    [~, ibest] = min(cost);
    track(k, :) = pts(ibest, :);
end
end

function displayPath = localTrackDisplayPath(track, trackValid, G, startIdx, endIdx)
displayPath = zeros(0, 2);
validIdx = find(trackValid(:)');
if isempty(validIdx)
    return
end
displayPath = track(validIdx(1), :);
for ii = 1:(numel(validIdx) - 1)
    k1 = validIdx(ii);
    k2 = validIdx(ii + 1);
    p1 = round(track(k1, :));
    p2 = round(track(k2, :));
    [seg, pathMask] = localTransitionSupportSegment(G, startIdx, endIdx, k1, k2, p1, p2, 2);
    if isempty(seg)
        seg = [p1; p2];
        pathMask = false(size(G.frameMasks{1}));
    end
    seg = localOrthogonalizeDisplaySegment(pathMask, seg);
    if ~isempty(displayPath) && isequal(displayPath(end, :), seg(1, :))
        seg = seg(2:end, :);
    end
    displayPath = [displayPath; seg]; %#ok<AGROW>
end
end

function [seg, pathMask] = localTransitionSupportSegment(G, startIdx, endIdx, k1, k2, p1, p2, maxExtraFrames)
seg = zeros(0, 2);
pathMask = false(size(G.frameMasks{1}));
lo = max(startIdx, k1);
for extra = 0:max(0, round(maxExtraFrames))
    hi = min(endIdx, k2 + extra);
    mask = false(size(G.frameMasks{1}));
    for k = lo:hi
        mask = mask | logical(G.frameMasks{k});
    end
    pathMask = localRepairInteriorCracks(mask);
    seg = localMaskedShortestPath(pathMask, p1, p2);
    if ~isempty(seg)
        return
    end
end
end

function maskOut = localRepairInteriorCracks(maskIn)
maskOut = logical(maskIn);
if isempty(maskOut) || ~any(maskOut(:))
    return
end
filledMask = imfill(maskOut, 'holes');
bridgeMask = imclose(maskOut, strel('square', 3));
bridgeMask = bridgeMask | imclose(maskOut, strel('line', 5, 0));
bridgeMask = bridgeMask | imclose(maskOut, strel('line', 5, 90));
bridgeMask = bridgeMask | imclose(maskOut, strel('line', 5, 45));
bridgeMask = bridgeMask | imclose(maskOut, strel('line', 5, 135));
maskOut = maskOut | (bridgeMask & filledMask);
end

function seg = localOrthogonalizeDisplaySegment(mask, seg)
if size(seg, 1) <= 1 || isempty(mask) || ~any(mask(:))
    return
end
out = seg(1, :);
for ii = 1:(size(seg, 1) - 1)
    a = seg(ii, :);
    b = seg(ii + 1, :);
    step = b - a;
    if all(abs(step) == [1 1])
        c1 = [b(1), a(2)];
        c2 = [a(1), b(2)];
        useC1 = localPointInMask(mask, c1);
        useC2 = localPointInMask(mask, c2);
        if useC1
            out(end + 1, :) = c1; %#ok<AGROW>
        elseif useC2
            out(end + 1, :) = c2; %#ok<AGROW>
        end
    end
    out(end + 1, :) = b; %#ok<AGROW>
end
seg = localUniqueConsecutiveRows(out);
end

function track = localConnectedPixelWalkTrack(G, W, startIdx, endIdx)
nFrames = numel(G.frames);
track = NaN(nFrames, 2);
if isempty(startIdx) || isempty(endIdx) || endIdx < startIdx
    return
end
originPoint = localSnapPointToMask(logical(G.frameMasks{startIdx}), W.branchOriginXY);
peakPoint = localSnapPointToMask(logical(G.frameMasks{endIdx}), W.branchPeakPointXY);
track(startIdx, :) = originPoint;
walkFailed = false;
futureMasks = cell(nFrames, 1);
futureMasks{endIdx} = logical(G.frameMasks{endIdx});
for k = (endIdx - 1):-1:startIdx
    futureMasks{k} = logical(G.frameMasks{k}) | futureMasks{k + 1};
end
for k = (startIdx + 1):endIdx
    prevPt = track(k - 1, :);
    if any(~isfinite(prevPt))
        walkFailed = true;
        break
    end
    prevMask = localRepairInteriorCracks(logical(G.frameMasks{k - 1}));
    currMask = localRepairInteriorCracks(logical(G.frameMasks{k}));
    unionMask = prevMask | currMask;
    reachableMask = localConnectedRegionFromSeed(unionMask, prevPt);
    if ~any(reachableMask(:))
        walkFailed = true;
        break
    end
    currReachable = reachableMask & currMask;
    if ~any(currReachable(:))
        walkFailed = true;
        break
    end
    peakConnected = localConnectedRegionFromSeed(localRepairInteriorCracks(futureMasks{k}), peakPoint);
    if any(peakConnected(:))
        currTowardPeak = currReachable & peakConnected;
        if any(currTowardPeak(:))
            currReachable = currTowardPeak;
        end
    end
    if k == endIdx
        if localPointInMask(currReachable, peakPoint)
            track(k, :) = peakPoint;
            continue
        else
            walkFailed = true;
            break
        end
    end
    prevPrevPt = [NaN NaN];
    if k - 2 >= 1
        prevPrevPt = track(k - 2, :);
    end
    nextPt = localSelectNextSupportPoint(G, k, endIdx, currReachable, prevPt, prevPrevPt, originPoint, peakPoint, futureMasks{k}, currMask);
    if any(~isfinite(nextPt))
        walkFailed = true;
        break
    end
    track(k, :) = nextPt;
end
if walkFailed || ~(endIdx <= size(track, 1) && all(isfinite(track(endIdx, :))) && ...
        all(round(track(endIdx, :)) == round(peakPoint)))
    track = localCumulativeMaskTrack(G, startIdx, endIdx, originPoint, peakPoint, false);
end
end

function pt = localSelectNextSupportPoint(G, k, endIdx, currReachable, prevPt, prevPrevPt, originPoint, peakPoint, futureMask, currMask)
pt = [NaN NaN];
if isempty(currReachable) || ~any(currReachable(:))
    return
end
futureMask = localRepairInteriorCracks(logical(futureMask));
distMap = localMaskedDistanceMap(futureMask, peakPoint);
corridorPath = localMaskedShortestPath(futureMask, prevPt, peakPoint);
peakLinkMask = [];
if k == endIdx - 1
    peakLinkMask = localRepairInteriorCracks(currMask | logical(futureMask));
end
newMask = false(size(currReachable));
if k <= numel(G.newMasks) && ~isempty(G.newMasks{k})
    newMask = logical(G.newMasks{k});
end
newMask = currReachable & newMask;
pt = localBestSupportCandidate(G, k, newMask, prevPt, prevPrevPt, originPoint, peakPoint, distMap, corridorPath, true, peakLinkMask);
if all(isfinite(pt))
    return
end
filledReachable = imfill(currReachable, 'holes');
edgeMask = bwperim(filledReachable, 8) & currReachable;
pt = localBestSupportCandidate(G, k, edgeMask, prevPt, prevPrevPt, originPoint, peakPoint, distMap, corridorPath, false, peakLinkMask);
if all(isfinite(pt))
    return
end
pt = localBestSupportCandidate(G, k, currReachable, prevPt, prevPrevPt, originPoint, peakPoint, distMap, corridorPath, false, peakLinkMask);
end

function pt = localBestSupportCandidate(G, k, mask, prevPt, prevPrevPt, originPoint, peakPoint, distMap, corridorPath, isPinkTier, peakLinkMask)
pt = [NaN NaN];
if isempty(mask) || ~any(mask(:))
    return
end
[pts, bonus] = localSupportCandidates(G, k, mask, originPoint, isPinkTier);
if isempty(pts)
    return
end
keep = localGeodesicBandMask(pts, corridorPath, 8);
if any(keep)
    pts = pts(keep, :);
    bonus = bonus(keep);
end
if ~isempty(peakLinkMask)
    keep = false(size(pts, 1), 1);
    for i = 1:size(pts, 1)
        keep(i) = localSegmentInsideMask(peakLinkMask, pts(i,:), peakPoint);
    end
    if any(keep)
        pts = pts(keep, :);
        bonus = bonus(keep);
    end
end
dPrev = hypot(pts(:,1) - prevPt(1), pts(:,2) - prevPt(2));
gPeak = localCandidateDistance(distMap, pts, peakPoint);
reversalPenalty = localImmediateReversalPenalty(pts, prevPt, prevPrevPt);
score = gPeak + 0.18 * dPrev + 0.012 * (dPrev .^ 2) + reversalPenalty - bonus;
[~, ibest] = min(score);
pt = pts(ibest, :);
end

function keep = localGeodesicBandMask(pts, pathPts, bandPx)
keep = true(size(pts, 1), 1);
if isempty(pts) || isempty(pathPts)
    return
end
for i = 1:size(pts, 1)
    d = hypot(pathPts(:,1) - pts(i,1), pathPts(:,2) - pts(i,2));
    if min(d) > bandPx
        keep(i) = false;
    end
end
end

function penalty = localImmediateReversalPenalty(pts, prevPt, prevPrevPt)
penalty = zeros(size(pts, 1), 1);
if any(~isfinite(prevPt)) || any(~isfinite(prevPrevPt))
    return
end
vPrev = prevPt - prevPrevPt;
prevNorm = hypot(vPrev(1), vPrev(2));
if prevNorm < 1
    return
end
for i = 1:size(pts, 1)
    vNow = pts(i,:) - prevPt;
    nowNorm = hypot(vNow(1), vNow(2));
    if nowNorm < 1
        continue
    end
    cosTheta = (vPrev(1) * vNow(1) + vPrev(2) * vNow(2)) / (prevNorm * nowNorm);
    if cosTheta < -0.2
        penalty(i) = 6.0 * (-cosTheta) + 0.08 * nowNorm;
    end
end
end

function [pts, bonus] = localSupportCandidates(G, k, mask, originPoint, isPinkTier)
pts = zeros(0, 2);
bonus = zeros(0, 1);
if isempty(mask) || ~any(mask(:))
    return
end
cc = bwconncomp(mask, 8);
for i = 1:cc.NumObjects
    compMask = false(size(mask));
    compMask(cc.PixelIdxList{i}) = true;
    [r, c] = find(compMask);
    compPts = [c, r];
    pts = [pts; compPts]; %#ok<AGROW>
    compBonus = zeros(size(compPts, 1), 1);
    if isPinkTier
        areaBonus = 0.06 * numel(cc.PixelIdxList{i});
        reachBonus = 0;
        if all(isfinite(originPoint))
            reachBonus = 0.08 * max(hypot(compPts(:,1) - originPoint(1), compPts(:,2) - originPoint(2)));
        end
        persistBonus = localPinkPersistenceBonus(G, k, compMask);
        compBonus(:) = areaBonus + reachBonus + persistBonus;
    end
    bonus = [bonus; compBonus]; %#ok<AGROW>
end
end

function bonus = localPinkPersistenceBonus(G, k, compMask)
bonus = 0;
if ~(isstruct(G) && isfield(G, 'newMasks') && isfield(G, 'frameMasks'))
    return
end
nFrames = numel(G.frameMasks);
for ahead = 1:2
    kk = k + ahead;
    if kk > nFrames
        break
    end
    nextSupport = false(size(compMask));
    if kk <= numel(G.newMasks) && ~isempty(G.newMasks{kk})
        nextSupport = logical(G.newMasks{kk});
    end
    nextSupport = localRepairInteriorCracks(nextSupport | logical(G.frameMasks{kk}));
    overlap = nnz(compMask & nextSupport);
    if overlap > 0
        bonus = bonus + (0.18 / ahead) * overlap;
    end
end
end

function d = localCandidateDistance(distMap, pts, peakPoint)
d = inf(size(pts, 1), 1);
if isempty(pts)
    return
end
if isempty(distMap)
    d = hypot(pts(:,1) - peakPoint(1), pts(:,2) - peakPoint(2));
    return
end
nRows = size(distMap, 1);
nCols = size(distMap, 2);
xs = min(max(round(pts(:,1)), 1), nCols);
ys = min(max(round(pts(:,2)), 1), nRows);
idx = sub2ind([nRows, nCols], ys, xs);
d = distMap(idx);
fallback = ~isfinite(d);
if any(fallback)
    d(fallback) = hypot(xs(fallback) - peakPoint(1), ys(fallback) - peakPoint(2));
end
end

function track = localCumulativeMaskTrack(G, startIdx, endIdx, originPoint, peakPoint, useFullEventUnion)
nFrames = numel(G.frames);
track = NaN(nFrames, 2);
if isempty(startIdx) || isempty(endIdx) || endIdx < startIdx || ...
        any(~isfinite(originPoint)) || any(~isfinite(peakPoint))
    return
end
mask = false(size(G.frameMasks{startIdx}));
if useFullEventUnion
    lo = 1;
    hi = numel(G.frameMasks);
else
    lo = startIdx;
    hi = endIdx;
end
for k = lo:hi
    mask = mask | logical(G.frameMasks{k});
end
mask = localRepairInteriorCracks(mask);
pathPts = localMaskedShortestPath(mask, originPoint, peakPoint);
if isempty(pathPts)
    return
end
pathPts = localUniqueConsecutiveRows(pathPts);
if isempty(pathPts)
    return
end
pathPos = 1;
track(startIdx, :) = originPoint;
for k = startIdx:endIdx
    frameMask = logical(G.frameMasks{k});
    if ~any(frameMask(:))
        continue
    end
    [bestPt, bestPos] = localLatestPathPointInMask(pathPts, pathPos, frameMask);
    if isempty(bestPt)
        [bestPt, bestPos] = localLatestPathPointInMask(pathPts, 1, frameMask);
    end
    if isempty(bestPt)
        continue
    end
    track(k, :) = bestPt;
    pathPos = max(pathPos, bestPos);
end
if localPointInMask(logical(G.frameMasks{endIdx}), peakPoint)
    track(endIdx, :) = peakPoint;
end
end

function rows = localUniqueConsecutiveRows(rows)
if size(rows, 1) <= 1
    return
end
keep = [true; any(diff(rows, 1, 1) ~= 0, 2)];
rows = rows(keep, :);
end

function [pt, idx] = localLatestPathPointInMask(pathPts, startPos, mask)
pt = [];
idx = [];
if isempty(pathPts) || isempty(mask) || ~any(mask(:))
    return
end
startPos = max(1, min(size(pathPts, 1), round(startPos)));
for ii = size(pathPts, 1):-1:startPos
    cand = pathPts(ii, :);
    if localPointInMask(mask, cand)
        pt = cand;
        idx = ii;
        return
    end
end
end

function pt = localLastPathPointInMask(pathPts, mask)
pt = [NaN NaN];
if isempty(pathPts) || isempty(mask)
    return
end
for ii = size(pathPts, 1):-1:1
    if localPointInMask(mask, pathPts(ii,:))
        pt = pathPts(ii,:);
        return
    end
end
end

function pt = localFirstPathPointInMask(pathPts, mask)
pt = [NaN NaN];
if isempty(pathPts) || isempty(mask)
    return
end
for ii = 1:size(pathPts, 1)
    if localPointInMask(mask, pathPts(ii,:))
        pt = pathPts(ii,:);
        return
    end
end
end

function distMap = localMaskedDistanceMap(mask, targetPt)
distMap = [];
if isempty(mask) || ~any(mask(:)) || any(~isfinite(targetPt))
    return
end
nRows = size(mask, 1);
nCols = size(mask, 2);
targetPt = localClampPoint(targetPt, nCols, nRows);
if ~mask(targetPt(2), targetPt(1))
    targetPt = localNearestMaskPixel(mask, targetPt);
    if isempty(targetPt)
        return
    end
end
distMap = inf(nRows, nCols);
queue = zeros(nnz(mask), 1);
head = 1;
tail = 1;
seedIdx = sub2ind([nRows, nCols], targetPt(2), targetPt(1));
queue(tail) = seedIdx;
distMap(seedIdx) = 0;
offsets = int32([-1 -1; 0 -1; 1 -1; -1 0; 1 0; -1 1; 0 1; 1 1]);
stepCost = [sqrt(2); 1; sqrt(2); 1; 1; sqrt(2); 1; sqrt(2)];
while head <= tail
    idx = queue(head);
    head = head + 1;
    [r, c] = ind2sub([nRows, nCols], idx);
    base = distMap(idx);
    for m = 1:8
        rr = r + offsets(m, 2);
        cc = c + offsets(m, 1);
        if rr < 1 || rr > nRows || cc < 1 || cc > nCols || ~mask(rr, cc)
            continue
        end
        nextIdx = sub2ind([nRows, nCols], rr, cc);
        cand = base + stepCost(m);
        if cand + 1e-9 < distMap(nextIdx)
            distMap(nextIdx) = cand;
            tail = tail + 1;
            queue(tail) = nextIdx;
        end
    end
end
end

function pathPts = localMaskedShortestPath(mask, startPt, endPt)
pathPts = zeros(0, 2);
if isempty(mask) || ~any(mask(:))
    return
end
nRows = size(mask, 1);
nCols = size(mask, 2);
startPt = localClampPoint(startPt, nCols, nRows);
endPt = localClampPoint(endPt, nCols, nRows);
if ~mask(startPt(2), startPt(1))
    startPt = localNearestMaskPixel(mask, startPt);
end
if ~mask(endPt(2), endPt(1))
    endPt = localNearestMaskPixel(mask, endPt);
end
if isempty(startPt) || isempty(endPt)
    return
end
startIdx = sub2ind([nRows, nCols], startPt(2), startPt(1));
endIdx = sub2ind([nRows, nCols], endPt(2), endPt(1));
if startIdx == endIdx
    pathPts = startPt;
    return
end
queue = zeros(nnz(mask), 1);
head = 1;
tail = 1;
queue(tail) = startIdx;
visited = false(nRows, nCols);
visited(startIdx) = true;
prev = zeros(nRows * nCols, 1, 'uint32');
offsets = int32([-1 -1; 0 -1; 1 -1; -1 0; 1 0; -1 1; 0 1; 1 1]);
found = false;
while head <= tail
    idx = queue(head);
    head = head + 1;
    [r, c] = ind2sub([nRows, nCols], idx);
    for m = 1:size(offsets, 1)
        rr = r + offsets(m, 2);
        cc = c + offsets(m, 1);
        if rr < 1 || rr > nRows || cc < 1 || cc > nCols || ~mask(rr, cc) || visited(rr, cc)
            continue
        end
        nextIdx = sub2ind([nRows, nCols], rr, cc);
        visited(nextIdx) = true;
        prev(nextIdx) = uint32(idx);
        tail = tail + 1;
        queue(tail) = nextIdx;
        if nextIdx == endIdx
            found = true;
            break
        end
    end
    if found
        break
    end
end
if ~found
    return
end
rev = zeros(0, 1, 'uint32');
idx = uint32(endIdx);
while idx ~= 0
    rev(end+1, 1) = idx; %#ok<AGROW>
    if idx == startIdx
        break
    end
    idx = prev(idx);
end
rev = flipud(rev);
[r, c] = ind2sub([nRows, nCols], double(rev));
pathPts = [c, r];
end

function regionMask = localConnectedRegionFromSeed(mask, seedPt)
regionMask = false(size(mask));
if isempty(mask) || ~any(mask(:)) || any(~isfinite(seedPt))
    return
end
nRows = size(mask, 1);
nCols = size(mask, 2);
seedPt = localClampPoint(seedPt, nCols, nRows);
if ~mask(seedPt(2), seedPt(1))
    seedPt = localNearestMaskPixel(mask, seedPt);
    if isempty(seedPt)
        return
    end
end
seedIdx = sub2ind([nRows, nCols], seedPt(2), seedPt(1));
queue = zeros(nnz(mask), 1);
head = 1;
tail = 1;
queue(tail) = seedIdx;
visited = false(nRows, nCols);
visited(seedIdx) = true;
offsets = int32([-1 -1; 0 -1; 1 -1; -1 0; 1 0; -1 1; 0 1; 1 1]);
while head <= tail
    idx = queue(head);
    head = head + 1;
    [r, c] = ind2sub([nRows, nCols], idx);
    for m = 1:size(offsets, 1)
        rr = r + offsets(m, 2);
        cc = c + offsets(m, 1);
        if rr < 1 || rr > nRows || cc < 1 || cc > nCols || ~mask(rr, cc) || visited(rr, cc)
            continue
        end
        nextIdx = sub2ind([nRows, nCols], rr, cc);
        visited(nextIdx) = true;
        tail = tail + 1;
        queue(tail) = nextIdx;
    end
end
regionMask = visited;
end

function pt = localNearestMaskPixel(mask, pt)
[r, c] = find(mask);
if isempty(r)
    pt = [];
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

function tf = localPointInMask(mask, pt)
tf = false;
if isempty(mask) || any(~isfinite(pt))
    return
end
nRows = size(mask, 1);
nCols = size(mask, 2);
pt = localClampPoint(pt, nCols, nRows);
tf = mask(pt(2), pt(1));
end

function tf = localMaskConnected(mask, startPt, endPt)
tf = false;
if isempty(mask) || ~any(mask(:)) || any(~isfinite(startPt)) || any(~isfinite(endPt))
    return
end
nRows = size(mask, 1);
nCols = size(mask, 2);
startPt = localClampPoint(startPt, nCols, nRows);
endPt = localClampPoint(endPt, nCols, nRows);
if ~mask(startPt(2), startPt(1)) || ~mask(endPt(2), endPt(1))
    return
end
if isequal(startPt, endPt)
    tf = true;
    return
end
queue = zeros(nnz(mask), 1);
head = 1;
tail = 1;
startIdx = sub2ind([nRows, nCols], startPt(2), startPt(1));
endIdx = sub2ind([nRows, nCols], endPt(2), endPt(1));
queue(tail) = startIdx;
visited = false(nRows, nCols);
visited(startIdx) = true;
offsets = int32([-1 -1; 0 -1; 1 -1; -1 0; 1 0; -1 1; 0 1; 1 1]);
while head <= tail
    idx = queue(head);
    head = head + 1;
    [r, c] = ind2sub([nRows, nCols], idx);
    for m = 1:size(offsets, 1)
        rr = r + offsets(m, 2);
        cc = c + offsets(m, 1);
        if rr < 1 || rr > nRows || cc < 1 || cc > nCols || ~mask(rr, cc) || visited(rr, cc)
            continue
        end
        nextIdx = sub2ind([nRows, nCols], rr, cc);
        if nextIdx == endIdx
            tf = true;
            return
        end
        visited(nextIdx) = true;
        tail = tail + 1;
        queue(tail) = nextIdx;
    end
end
end

function pt = localClampPoint(pt, nCols, nRows)
pt = round(pt);
pt(1) = min(max(pt(1), 1), nCols);
pt(2) = min(max(pt(2), 1), nRows);
end

function secondary = localSecondarySinks(nodesByFrame, costs, reachable, peakFrameIdx, peakNode)
secondary = {};
for k = 1:numel(nodesByFrame)
    if isempty(nodesByFrame{k})
        continue
    end
    if k == peakFrameIdx
        continue
    end
    if k < numel(nodesByFrame) && ~isempty(nodesByFrame{k + 1})
        continue
    end
    idx = find(reachable{k});
    for i = idx(:)'
        secondary{end+1,1} = struct('endFrameIdx', k, 'endNodeIdx', i, 'cost', costs{k}(i)); %#ok<AGROW>
    end
end
if isempty(secondary)
    return
end
[~, order] = sort(cellfun(@(s) s.cost, secondary), 'ascend');
secondary = secondary(order(1:min(3, numel(order))));
end

function alternatives = localFrameAlternatives(nodesByFrame, costs, forwardReachable, backwardReachable, nodePath, startIdx, endIdx)
alternatives = cell(numel(nodesByFrame), 1);
if isempty(nodePath)
    return
end
for k = startIdx:endIdx
    nodes = nodesByFrame{k};
    if isempty(nodes)
        continue
    end
    chosenIdx = nodePath(k - startIdx + 1);
    if ~isfinite(chosenIdx) || chosenIdx < 1 || chosenIdx > numel(nodes)
        continue
    end
    candidates = find(forwardReachable{k} & backwardReachable{k});
    candidates = candidates(candidates ~= chosenIdx);
    if isempty(candidates)
        continue
    end
    [~, order] = sort(costs{k}(candidates), 'ascend', 'MissingPlacement', 'last');
    altIdx = candidates(order(1));
    if ~isfinite(costs{k}(altIdx))
        continue
    end
    alternatives{k} = struct( ...
        'chosenIdx', chosenIdx, ...
        'chosenRep', nodes(chosenIdx).rep, ...
        'altIdx', altIdx, ...
        'altRep', nodes(altIdx).rep, ...
        'altCost', costs{k}(altIdx));
end
end

function localPlotAlternatives(ax, nodesByFrame, alternatives, frames, spatialRes)
for k = 1:numel(alternatives)
    alt = alternatives{k};
    if isempty(alt)
        continue
    end
    chosen = alt.chosenRep * spatialRes;
    other = alt.altRep * spatialRes;
    plot(ax, [chosen(1), other(1)], [chosen(2), other(2)], '--', ...
        'Color', [0.95 0.60 0.10], 'LineWidth', 1.0);
    plot(ax, other(1), other(2), 'o', ...
        'Color', [0.95 0.60 0.10], 'MarkerFaceColor', 'none', 'MarkerSize', 7, 'LineWidth', 1.25);
end
end

function [BWA_AQ2, label, sourceResult] = localResolveSource(source)
label = "";
sourceResult = [];
if isstruct(source) && isfield(source, 'derived') && isstruct(source.derived) && ...
        isfield(source.derived, 'context') && isstruct(source.derived.context) && ...
        isfield(source.derived.context, 'BWA_AQ2')
    sourceResult = source;
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
    sourceResult = source;
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
if isstruct(source) && isfield(source, 'fts1') && isfield(source, 'opts')
    BWA_AQ2 = source;
    return
end
error('wavefrontComponentGraphPrototype:BadSource', 'Unsupported source type.');
end

function tf = localResolveRefinePixelWalk(modeValue, isMulticellEvent)
if islogical(modeValue)
    tf = modeValue;
    return
end
modeValue = lower(char(string(modeValue)));
switch modeValue
    case 'on'
        tf = true;
    case 'off'
        tf = false;
    otherwise
        tf = isMulticellEvent;
end
end

function tf = localIsMulticellEvent(result, BWA_AQ2, eventID)
tf = false;
if isempty(result) || ~isstruct(result) || isempty(BWA_AQ2) || ~isstruct(BWA_AQ2)
    return
end
try
    if ~isfield(result, 'metrics') || ~isstruct(result.metrics) || ...
            ~isfield(result.metrics, 'tables') || ~isstruct(result.metrics.tables) || ...
            ~isfield(result.metrics.tables, 'event_partition')
        return
    end
    eventUID = localEventUIDFromGlobalIndex(BWA_AQ2, eventID);
    if strlength(eventUID) == 0
        return
    end
    T = result.metrics.tables.event_partition;
    if ~istable(T) || ~ismember('event_uid', T.Properties.VariableNames)
        return
    end
    keep = string(T.event_uid) == eventUID;
    if ~any(keep)
        return
    end
    if ismember('partition_cell_uid', T.Properties.VariableNames)
        cellIDs = unique(string(T.partition_cell_uid(keep)));
    elseif ismember('cell_uid', T.Properties.VariableNames)
        cellIDs = unique(string(T.cell_uid(keep)));
    else
        return
    end
    cellIDs = unique(cellIDs(cellIDs ~= ""));
    tf = numel(cellIDs) >= 2;
catch
    tf = false;
end
end

function eventUID = localEventUIDFromGlobalIndex(BWA_AQ2, eventID)
eventUID = "";
if eventID < 1
    return
end
fts = BWA_AQ2.fts1;
fileIDs = [];
if isfield(fts, 'basic') && isfield(fts.basic, 'fileID') && ~isempty(fts.basic.fileID)
    fileIDs = string(fts.basic.fileID);
elseif isfield(fts, 'fileID') && ~isempty(fts.fileID)
    fileIDs = string(fts.fileID);
end
if isempty(fileIDs) || eventID > numel(fileIDs)
    return
end
thisFile = fileIDs(eventID);
localEventID = sum(fileIDs(1:eventID) == thisFile);
eventUID = sprintf('%s_E%d', char(thisFile), localEventID);
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
