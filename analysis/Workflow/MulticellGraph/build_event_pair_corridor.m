function out = build_event_pair_corridor(source, varargin)
%BUILD_EVENT_PAIR_CORRIDOR Extract one delay-consistent source-target corridor.
%   OUT = BUILD_EVENT_PAIR_CORRIDOR(RESULT, 'eventUID', '<file>_E###', ...
%       'sourceCell', 'C21', 'targetCell', 'C23')

p = inputParser;
p.StructExpand = false;
addRequired(p, 'source');
addParameter(p, 'eventUID', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'eventID', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'sourceCell', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'targetCell', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'graph', struct(), @(x) isstruct(x) || isempty(x));
addParameter(p, 'pixelMap', struct(), @(x) isstruct(x) || isempty(x));
addParameter(p, 'showFigure', true, @islogical);
addParameter(p, 'visible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
addParameter(p, 'smoothSigma', 1.25, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'interpScale', 2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'turnPenalty', 0.8, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'anchorInteriorFrac', 0.75, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
addParameter(p, 'sourceAnchorInteriorFrac', NaN, @(x) isnumeric(x) && isscalar(x) && (isnan(x) || (isfinite(x) && x >= 0 && x <= 1)));
addParameter(p, 'targetAnchorInteriorFrac', NaN, @(x) isnumeric(x) && isscalar(x) && (isnan(x) || (isfinite(x) && x >= 0 && x <= 1)));
addParameter(p, 'edgeBandFrac', 0.15, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
addParameter(p, 'corridorRadiusPx', 6, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'corridorCostSlack', 1.35, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'corridorAdaptive', true, @islogical);
addParameter(p, 'targetWidthWeight', 80, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'maxRefineIter', 5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'verbose', false, @islogical);
addParameter(p, 'useFitMemory', true, @islogical);
addParameter(p, 'beamWidth', 4, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
parse(p, source, varargin{:});

if isstruct(p.Results.graph) && isfield(p.Results.graph, 'nodes') && isfield(p.Results.graph, 'edges')
    G = p.Results.graph;
else
    G = build_event_cell_adjacency_graph(source, ...
        'eventUID', p.Results.eventUID, ...
        'eventID', p.Results.eventID, ...
        'showFigure', false);
end
if isstruct(p.Results.pixelMap) && isfield(p.Results.pixelMap, 'delay_sec') && isfield(p.Results.pixelMap, 'footprint_mask')
    M = p.Results.pixelMap;
else
    M = build_event_pixel_activation_map(source, ...
        'eventUID', G.event_uid, ...
        'showFigure', false, ...
        'delayUnits', 'sec');
end

sIdx = localResolveNodeIdxPair(G.nodes, string(p.Results.sourceCell));
tIdx = localResolveNodeIdxPair(G.nodes, string(p.Results.targetCell));
if isempty(sIdx) || isempty(tIdx)
    error('build_event_pair_corridor:MissingCell', 'Source or target cell not found.');
end

srcMask = G.nodes.partition_mask{sIdx};
tgtMask = G.nodes.partition_mask{tIdx};
footMask = localFillSmallHoles(logical(M.footprint_mask), 64);
delaySec = M.delay_sec;
params = struct();
params.srcInteriorFrac = double(p.Results.sourceAnchorInteriorFrac);
if ~isfinite(params.srcInteriorFrac)
    params.srcInteriorFrac = double(p.Results.anchorInteriorFrac);
end
params.tgtInteriorFrac = double(p.Results.targetAnchorInteriorFrac);
if ~isfinite(params.tgtInteriorFrac)
    params.tgtInteriorFrac = double(p.Results.anchorInteriorFrac);
end
params.edgeBandFrac = double(p.Results.edgeBandFrac);
params.corridorCostSlack = double(p.Results.corridorCostSlack);
params.corridorRadiusPx = double(p.Results.corridorRadiusPx);
params.turnPenalty = double(p.Results.turnPenalty);
params.targetWidthWeight = double(p.Results.targetWidthWeight);
params.smoothSigma = double(p.Results.smoothSigma);
params.interpScale = round(double(p.Results.interpScale));
params.corridorAdaptive = logical(p.Results.corridorAdaptive);
maxIter = max(1, round(double(p.Results.maxRefineIter)));
verbose = logical(p.Results.verbose);
trialBase = localPairTrialBase(footMask, srcMask, tgtMask);
targetWidthLocal = localEstimateTargetWidth(trialBase, srcMask, tgtMask);
pairFeatures = localPairFeatureStruct(srcMask, tgtMask, delaySec, targetWidthLocal);
fitMemoryMeta = struct('used', false, 'nMatches', 0, 'distance', NaN, 'usedAccepted', false);
if p.Results.useFitMemory
    [params, fitMemoryMeta] = pair_corridor_fit_memory('warm_start', pairFeatures, params);
    if verbose && fitMemoryMeta.used
        fprintf('pair warm start: matches=%d distance=%.3f usedAccepted=%d\n', fitMemoryMeta.nMatches, fitMemoryMeta.distance, fitMemoryMeta.usedAccepted);
    end
end
outcome = struct('acceptable', false, 'reason', "init");
pathXY = zeros(0,2);
pathDelaySec = zeros(0,1);
supportMask = false(size(footMask));
srcPlotAnchorLin = NaN;
tgtPlotAnchorLin = NaN;
borderAnchorLin = NaN;
beamWidth = max(1, round(double(p.Results.beamWidth)));
if fitMemoryMeta.used
    if fitMemoryMeta.usedAccepted && isfinite(fitMemoryMeta.distance) && fitMemoryMeta.distance <= 1.5
        maxIter = min(maxIter, 3);
        beamWidth = min(beamWidth, 2);
    elseif isfinite(fitMemoryMeta.distance) && fitMemoryMeta.distance <= 2.5
        maxIter = min(maxIter, 4);
        beamWidth = min(beamWidth, 3);
    end
end
state0 = struct('pathXY', pathXY, 'pathDelaySec', pathDelaySec, 'supportMask', supportMask, ...
    'srcPlotAnchorLin', srcPlotAnchorLin, 'tgtPlotAnchorLin', tgtPlotAnchorLin, 'borderAnchorLin', borderAnchorLin, ...
    'outcome', outcome, 'params', params, 'score', inf);
bestScore = inf;
bestState = state0;
beams = state0;
for iter = 1:maxIter
    pool = repmat(state0, 0, 1);
    poolCount = 0;
    for ib = 1:numel(beams)
        candidates = localGenerateParamCandidates(beams(ib).params, iter);
        for ic = 1:numel(candidates)
            candParams = candidates(ic);
            [srcSeedMask, tgtSeedMask] = localInteriorSeeds(srcMask, tgtMask, footMask, candParams.srcInteriorFrac, candParams.tgtInteriorFrac);
            [srcAnchorLin, borderAnchorLin, tgtAnchorLin] = localPickCorridorAnchors(footMask, delaySec, srcMask, tgtMask, srcSeedMask, tgtSeedMask);
            [candPathXY, candPathDelaySec, candSupportMask, candSrcPlotAnchorLin, candTgtPlotAnchorLin] = localExtractCorridorPath(footMask, delaySec, srcSeedMask, tgtSeedMask, srcAnchorLin, borderAnchorLin, tgtAnchorLin, candParams.smoothSigma, candParams.interpScale, candParams.turnPenalty, candParams.edgeBandFrac, candParams.corridorRadiusPx, candParams.corridorCostSlack, candParams.corridorAdaptive, candParams.targetWidthWeight);
            candOutcome = localJudgeCorridorOutcome(candSupportMask, candPathXY, srcSeedMask, tgtSeedMask, srcMask, tgtMask, footMask, targetWidthLocal);
            candScore = localOutcomeObjective(candOutcome, candPathXY, targetWidthLocal);
            if verbose
                fprintf('pair refine iter %d beam %d cand %d score=%.2f reason=%s nPath=%d width=%.2f spill=%.2f srcOcc=%.2f tgtOcc=%.2f branches=%.0f\n', ...
                    iter, ib, ic, candScore, string(candOutcome.reason), size(candPathXY,1), ...
                    localGetOutcomeField(candOutcome,'meanWidth'), localGetOutcomeField(candOutcome,'targetSpillFrac'), ...
                    localGetOutcomeField(candOutcome,'sourceOccupancyFrac'), localGetOutcomeField(candOutcome,'targetOccupancyFrac'), ...
                    localGetOutcomeField(candOutcome,'branchPoints'));
            end
            poolCount = poolCount + 1;
            pool(poolCount) = struct('pathXY', candPathXY, 'pathDelaySec', candPathDelaySec, 'supportMask', candSupportMask, ...
                'srcPlotAnchorLin', candSrcPlotAnchorLin, 'tgtPlotAnchorLin', candTgtPlotAnchorLin, 'borderAnchorLin', borderAnchorLin, ...
                'outcome', candOutcome, 'params', candParams, 'score', candScore);
        end
    end
    if isempty(pool)
        break
    end
    [~, ord] = sort([pool.score], 'ascend');
    pool = pool(ord);
    pool = localUniqueParamStates(pool);
    beams = pool(1:min(beamWidth, numel(pool)));
    iterBestState = beams(1);
    iterBestScore = iterBestState.score;
    if verbose
        fprintf('pair refine iter %d best score=%.2f reason=%s kept=%d\n', iter, iterBestScore, string(iterBestState.outcome.reason), numel(beams));
    end
    if iterBestScore < bestScore
        bestScore = iterBestScore;
        bestState = iterBestState;
    else
        break
    end
    accIdx = find(arrayfun(@(b) b.outcome.acceptable, beams), 1, 'first');
    if ~isempty(accIdx)
        bestState = beams(accIdx);
        bestScore = bestState.score;
        break
    end
    if iter >= 3 && localBeamConverged(beams)
        if verbose
            fprintf('pair refine iter %d converged: beam family collapsed\n', iter);
        end
        break
    end
end
pathXY = bestState.pathXY;
pathDelaySec = bestState.pathDelaySec;
supportMask = bestState.supportMask;
srcPlotAnchorLin = bestState.srcPlotAnchorLin;
tgtPlotAnchorLin = bestState.tgtPlotAnchorLin;
borderAnchorLin = bestState.borderAnchorLin;
outcome = bestState.outcome;
params = bestState.params;
if string(outcome.reason) == "init"
    outcome = localJudgeCorridorOutcome(supportMask, pathXY, false(size(srcMask)), false(size(tgtMask)), srcMask, tgtMask, footMask, targetWidthLocal);
end

out = struct();
out.event_uid = G.event_uid;
out.file_id = G.file_id;
out.source_cell = string(G.nodes.label(sIdx));
out.target_cell = string(G.nodes.label(tIdx));
out.source_idx = sIdx;
out.target_idx = tIdx;
out.pixel_map = M;
out.graph = G;
out.path_xy_px = pathXY;
out.path_delay_sec = pathDelaySec;
out.support_mask = supportMask;
out.source_anchor_xy_px = localLinToXY(size(footMask), srcPlotAnchorLin);
out.border_anchor_xy_px = localLinToXY(size(footMask), borderAnchorLin);
out.target_anchor_xy_px = localLinToXY(size(footMask), tgtPlotAnchorLin);
out.refine_outcome = outcome;
out.fit_memory = fitMemoryMeta;

if p.Results.useFitMemory
    pair_corridor_fit_memory('record', pairFeatures, params, outcome, bestScore);
end

if p.Results.showFigure
    out.figure = localPlotPairCorridor(out, char(string(p.Results.visible)));
end
end

function beamsOut = localUniqueParamStates(beamsIn)
if isempty(beamsIn)
    beamsOut = beamsIn;
    return
end
keys = strings(1, numel(beamsIn));
for i = 1:numel(beamsIn)
    st = struct();
    if isfield(beamsIn(i), 'outcome') && isfield(beamsIn(i).outcome, 'stats')
        st = beamsIn(i).outcome.stats;
    end
    if ~isempty(fieldnames(st))
        keys(i) = sprintf('%s|%.1f|%.2f|%.2f|%.2f|%.0f', ...
            string(beamsIn(i).outcome.reason), ...
            localRoundOrZero(localGetOutcomeField(beamsIn(i).outcome, 'meanWidth'), 1), ...
            localRoundOrZero(localGetOutcomeField(beamsIn(i).outcome, 'targetSpillFrac'), 2), ...
            localRoundOrZero(localGetOutcomeField(beamsIn(i).outcome, 'sourceOccupancyFrac'), 2), ...
            localRoundOrZero(localGetOutcomeField(beamsIn(i).outcome, 'targetOccupancyFrac'), 2), ...
            localRoundOrZero(localGetOutcomeField(beamsIn(i).outcome, 'branchPoints'), 0));
    else
        p = beamsIn(i).params;
        keys(i) = sprintf('%.3f|%.3f|%.3f|%.3f|%.3f|%.3f|%.3f|%d', ...
            p.srcInteriorFrac, p.tgtInteriorFrac, p.edgeBandFrac, p.corridorCostSlack, ...
            p.corridorRadiusPx, p.turnPenalty, p.targetWidthWeight, p.interpScale);
    end
end
[~, ia] = unique(keys, 'stable');
beamsOut = beamsIn(sort(ia));
end

function v = localRoundOrZero(v, digits)
if ~(isfinite(v) && isscalar(v))
    v = 0;
    return
end
scale = 10^digits;
v = round(v * scale) / scale;
end

function tf = localBeamConverged(beams)
tf = false;
if numel(beams) < 2
    return
end
reasons = strings(1, numel(beams));
widths = nan(1, numel(beams));
spills = nan(1, numel(beams));
srcOcc = nan(1, numel(beams));
tgtOcc = nan(1, numel(beams));
for i = 1:numel(beams)
    reasons(i) = string(beams(i).outcome.reason);
    widths(i) = localGetOutcomeField(beams(i).outcome, 'meanWidth');
    spills(i) = localGetOutcomeField(beams(i).outcome, 'targetSpillFrac');
    srcOcc(i) = localGetOutcomeField(beams(i).outcome, 'sourceOccupancyFrac');
    tgtOcc(i) = localGetOutcomeField(beams(i).outcome, 'targetOccupancyFrac');
end
sameReason = numel(unique(reasons)) == 1;
narrowFamily = range(widths) <= 0.9 && range(spills) <= 0.04 && range(srcOcc) <= 0.04 && range(tgtOcc) <= 0.04;
tf = sameReason && narrowFamily;
end

function feat = localPairFeatureStruct(srcMask, tgtMask, delaySec, targetWidthLocal)
feat = struct();
feat.srcArea = nnz(srcMask);
feat.tgtArea = nnz(tgtMask);
feat.areaRatio = feat.srcArea / max(feat.tgtArea, 1);
feat.contactFrac = nnz(imdilate(srcMask, strel('disk', 1, 0)) & tgtMask) / max(min(feat.srcArea, feat.tgtArea), 1);
[srcR, srcC] = find(srcMask);
[tgtR, tgtC] = find(tgtMask);
if isempty(srcR) || isempty(tgtR)
    feat.centroidDist = NaN;
else
    feat.centroidDist = hypot(mean(srcC) - mean(tgtC), mean(srcR) - mean(tgtR));
end
srcVals = delaySec(srcMask & isfinite(delaySec));
tgtVals = delaySec(tgtMask & isfinite(delaySec));
feat.srcDelay = mean(srcVals, 'omitnan');
feat.tgtDelay = mean(tgtVals, 'omitnan');
feat.delayDiff = abs(feat.tgtDelay - feat.srcDelay);
feat.targetWidth = targetWidthLocal;
end

function mask = localFillSmallHoles(mask, maxHoleArea)
mask = logical(mask);
filled = imfill(mask, 'holes');
holes = filled & ~mask;
if ~any(holes(:))
    return
end
cc = bwconncomp(holes, 8);
keep = false(size(mask));
for i = 1:cc.NumObjects
    pix = cc.PixelIdxList{i};
    if numel(pix) <= maxHoleArea
        keep(pix) = true;
    end
end
mask = mask | keep;
end

function [srcSeedMask, tgtSeedMask] = localInteriorSeeds(srcMask, tgtMask, footMask, srcInteriorFrac, tgtInteriorFrac)
srcMask = logical(srcMask) & footMask;
tgtMask = logical(tgtMask) & footMask;
overlap = srcMask & tgtMask;
srcSeedMask = localInteriorCandidateMask(srcMask & ~overlap, srcInteriorFrac);
tgtSeedMask = localInteriorCandidateMask(tgtMask & ~overlap, tgtInteriorFrac);
if ~any(srcSeedMask(:)), srcSeedMask = srcMask & ~overlap; end
if ~any(tgtSeedMask(:)), tgtSeedMask = tgtMask & ~overlap; end
if ~any(srcSeedMask(:)), srcSeedMask = srcMask; end
if ~any(tgtSeedMask(:)), tgtSeedMask = tgtMask; end
end

function idx = localResolveNodeIdxPair(nodes, cellLabel)
idx = find(string(nodes.label) == string(cellLabel) | string(nodes.cell_uid) == string(cellLabel), 1, 'first');
end

function baseMask = localPairTrialBase(footMask, srcMask, tgtMask)
baseMask = footMask & (imdilate(srcMask, strel('disk', 10, 0)) | ...
    imdilate(tgtMask, strel('disk', 10, 0)) | ...
    ((bwdist(srcMask) + bwdist(tgtMask)) <= 24));
baseMask(srcMask | tgtMask) = true;
end

function outcome = localJudgeCorridorOutcome(supportMask, pathXY, srcSeedMask, tgtSeedMask, srcMask, tgtMask, footMask, targetWidth)
outcome = struct('acceptable', false, 'reason', "must_exist", 'stats', struct());
stats = struct('exists', any(supportMask(:)), ...
    'nPath', size(pathXY,1), ...
    'linksSourceInterior', any(supportMask(srcSeedMask)), ...
    'linksTargetInterior', any(supportMask(tgtSeedMask)), ...
    'touchesSource', any(supportMask(srcMask)), ...
    'touchesTarget', any(supportMask(tgtMask)), ...
    'crossesBoundary', any(supportMask(footMask & ~(srcMask | tgtMask))), ...
    'area', 0, ...
    'branchPoints', inf, ...
    'meanWidth', 0, ...
    'sourceSpillFrac', 1, ...
    'targetSpillFrac', 1, ...
    'sourceReachFrac', 0, ...
    'targetReachFrac', 0, ...
    'sourceOccupancyFrac', 1, ...
    'targetOccupancyFrac', 1);
if stats.exists
    simp = localCorridorSimplicity(supportMask, srcMask, tgtMask, ones(size(supportMask)));
    stats.area = simp.area;
    stats.branchPoints = simp.branchPoints;
    stats.meanWidth = simp.meanWidth;
    stats.sourceReachFrac = localInteriorReachFrac(supportMask, srcMask);
    stats.targetReachFrac = localInteriorReachFrac(supportMask, tgtMask);
    stats.sourceSpillFrac = localTargetSpillFrac(supportMask, tgtMask, srcMask);
    stats.targetSpillFrac = localTargetSpillFrac(supportMask, srcMask, tgtMask);
    stats.sourceOccupancyFrac = localMaskOccupancyFrac(supportMask, srcMask);
    stats.targetOccupancyFrac = localMaskOccupancyFrac(supportMask, tgtMask);
end
outcome.stats = stats;
if ~stats.exists || stats.nPath < 6
    outcome.reason = "must_exist";
    return
end
if ~(stats.linksSourceInterior && stats.linksTargetInterior && stats.touchesSource && stats.touchesTarget && stats.crossesBoundary && ...
        stats.sourceReachFrac >= 0.35 && stats.targetReachFrac >= 0.35)
    outcome.reason = "must_link_interiors";
    return
end
targetOccMax = max(0.45, min(0.80, 0.15 + 1.35 * stats.targetReachFrac));
sourceOccMax = max(0.50, min(0.90, 0.20 + 1.50 * stats.sourceReachFrac));
maxBranchPoints = 2 + double(stats.meanWidth >= max(6.0, 1.1 * targetWidth));
if stats.meanWidth > max(4.0, 1.6 * targetWidth) || stats.targetSpillFrac > 0.28 || stats.sourceSpillFrac > 0.28 || ...
        stats.branchPoints > maxBranchPoints || ...
        stats.targetOccupancyFrac > targetOccMax || stats.sourceOccupancyFrac > sourceOccMax
    outcome.reason = "too_broad";
    return
end
outcome.acceptable = true;
outcome.reason = "ok";
end

function params = localAdjustCorridorParams(params, outcome)
switch string(outcome.reason)
    case {"must_exist", "must_link_interiors"}
        params.srcInteriorFrac = max(0.20, params.srcInteriorFrac - 0.15);
        params.tgtInteriorFrac = max(0.20, params.tgtInteriorFrac - 0.10);
        params.edgeBandFrac = max(0.05, params.edgeBandFrac - 0.03);
        params.corridorCostSlack = min(1.8, params.corridorCostSlack + 0.08);
        params.corridorRadiusPx = min(12, params.corridorRadiusPx + 1);
        params.targetWidthWeight = max(20, params.targetWidthWeight - 10);
    case "too_broad"
        params.srcInteriorFrac = min(0.95, params.srcInteriorFrac + 0.05);
        params.tgtInteriorFrac = min(0.95, params.tgtInteriorFrac + 0.05);
        params.edgeBandFrac = min(0.25, params.edgeBandFrac + 0.02);
        params.corridorCostSlack = max(1.08, params.corridorCostSlack - 0.06);
        params.corridorRadiusPx = max(3, params.corridorRadiusPx - 1);
        params.targetWidthWeight = min(160, params.targetWidthWeight + 15);
    otherwise
end
end

function score = localOutcomeObjective(outcome, pathXY, targetWidth)
st = struct();
if isfield(outcome,'stats'), st = outcome.stats; end
meanWidth = 0; nPath = size(pathXY,1);
if isfield(st,'meanWidth'), meanWidth = double(st.meanWidth); end
pen = 0;
reason = string(outcome.reason);
switch reason
    case "ok"
        pen = 0;
    case "too_broad"
        pen = 200;
    otherwise
        pen = 400;
end
srcSpill = localGetOutcomeField(outcome,'sourceSpillFrac');
spill = localGetOutcomeField(outcome,'targetSpillFrac');
srcOcc = localGetOutcomeField(outcome,'sourceOccupancyFrac');
tgtOcc = localGetOutcomeField(outcome,'targetOccupancyFrac');
srcReach = localGetOutcomeField(outcome,'sourceReachFrac');
tgtReach = localGetOutcomeField(outcome,'targetReachFrac');
widthErr = abs(meanWidth - targetWidth);
sourceOccMax = max(0.50, min(0.90, 0.20 + 1.50 * srcReach));
targetOccMax = max(0.45, min(0.80, 0.15 + 1.35 * tgtReach));
srcOverflow = max(0, srcOcc - sourceOccMax);
tgtOverflow = max(0, tgtOcc - targetOccMax);
widthPenalty = 18 * widthErr;
if srcOverflow > 0 || tgtOverflow > 0
    widthPenalty = 10 * widthErr;
end
branchPoints = localGetOutcomeField(outcome,'branchPoints');
branchPenalty = 45 * branchPoints + 110 * max(0, branchPoints - 2).^2;
overfillPenalty = 1100 * tgtOverflow + 850 * srcOverflow + 280 * tgtOcc + 140 * srcOcc;
reachPenalty = 160 * max(0, 0.35 - srcReach) + 180 * max(0, 0.35 - tgtReach);
spillPenalty = 240 * srcSpill + 320 * spill;
score = pen + widthPenalty + spillPenalty + overfillPenalty + branchPenalty + reachPenalty - 0.35*nPath;
end

function v = localGetOutcomeField(outcome, fieldName)
v = 0;
if isfield(outcome,'stats') && isfield(outcome.stats, fieldName)
    v = double(outcome.stats.(fieldName));
end
end

function candidates = localGenerateParamCandidates(params, iter)
base = params;
candidates = base;
steps = [1, -1];
for s = steps
    p = base; p.corridorCostSlack = min(1.8, max(1.02, base.corridorCostSlack + s*0.08)); candidates(end+1) = p; %#ok<AGROW>
    p = base; p.edgeBandFrac = min(0.25, max(0.05, base.edgeBandFrac + s*0.03)); candidates(end+1) = p; %#ok<AGROW>
    p = base; p.srcInteriorFrac = min(0.95, max(0.20, base.srcInteriorFrac + s*0.10)); candidates(end+1) = p; %#ok<AGROW>
    p = base; p.tgtInteriorFrac = min(0.95, max(0.20, base.tgtInteriorFrac + s*0.10)); candidates(end+1) = p; %#ok<AGROW>
    p = base; p.corridorRadiusPx = min(12, max(3, base.corridorRadiusPx + s*1)); candidates(end+1) = p; %#ok<AGROW>
    p = base; p.targetWidthWeight = min(160, max(20, base.targetWidthWeight + s*15)); candidates(end+1) = p; %#ok<AGROW>
end
if iter == 1
    p = base; p.srcInteriorFrac = min(0.95, max(0.20, base.srcInteriorFrac - 0.20)); p.tgtInteriorFrac = min(0.95, max(0.20, base.tgtInteriorFrac - 0.15)); candidates(end+1)=p; %#ok<AGROW>
    p = base; p.corridorCostSlack = min(1.8, base.corridorCostSlack + 0.15); p.edgeBandFrac = max(0.05, base.edgeBandFrac - 0.05); candidates(end+1)=p; %#ok<AGROW>
end
end

function seedMask = localInteriorCandidateMask(mask, interiorFrac)
seedMask = false(size(mask));
if ~any(mask(:))
    return
end
d = bwdist(~mask);
maxd = max(d(mask), [], 'omitnan');
if ~(isfinite(maxd) && maxd > 0)
    seedMask = mask;
    return
end
seedMask = mask & (d >= interiorFrac * maxd);
if ~any(seedMask(:))
    levels = [0.60 0.50 0.40 0.30 0.20 0.10 0.05 0.0];
    for k = 1:numel(levels)
        seedMask = mask & (d >= levels(k) * maxd);
        if any(seedMask(:))
            break
        end
    end
end
if any(seedMask(:))
    tmp = imerode(seedMask, strel('disk', 1, 0));
    if any(tmp(:))
        seedMask = tmp;
    end
end
end

function [srcAnchorLin, borderAnchorLin, tgtAnchorLin] = localPickCorridorAnchors(footMask, delaySec, srcMask, tgtMask, srcSeedMask, tgtSeedMask)
srcAnchorLin = NaN;
borderAnchorLin = NaN;
tgtAnchorLin = NaN;
DS = bwdistgeodesic(footMask, srcSeedMask, 'quasi-euclidean');
DT = bwdistgeodesic(footMask, tgtSeedMask, 'quasi-euclidean');
if ~any(isfinite(DS(:))) || ~any(isfinite(DT(:)))
    return
end
srcVals = delaySec(srcSeedMask & isfinite(delaySec));
tgtVals = delaySec(tgtSeedMask & isfinite(delaySec));
srcDelay = mean(srcVals, 'omitnan');
tgtDelay = mean(tgtVals, 'omitnan');
if ~(isfinite(srcDelay) && isfinite(tgtDelay))
    srcDelay = 0;
    tgtDelay = max(delaySec(isfinite(delaySec)), [], 'omitnan');
end
midDelay = 0.5 * (srcDelay + tgtDelay);
delaySpan = max(abs(tgtDelay - srcDelay), eps);

srcCost = DT + abs(delaySec - srcDelay) ./ delaySpan;
srcCost(~srcSeedMask | ~isfinite(srcCost)) = inf;
[~, srcAnchorLin] = min(srcCost(:));

tgtCost = DS + abs(delaySec - tgtDelay) ./ delaySpan;
tgtCost(~tgtSeedMask | ~isfinite(tgtCost)) = inf;
[~, tgtAnchorLin] = min(tgtCost(:));

overlap = srcMask & tgtMask & footMask;
if any(overlap(:))
    midCost = abs(delaySec - midDelay);
    midCost(~overlap | ~isfinite(midCost)) = inf;
    [~, borderAnchorLin] = min(midCost(:));
    return
end

handoffMask = footMask & ~srcSeedMask & ~tgtSeedMask;
handoffCost = DS + DT + 2.0 * abs(delaySec - midDelay) ./ delaySpan;
handoffCost(~handoffMask | ~isfinite(handoffCost)) = inf;
if any(isfinite(handoffCost(:)))
    [~, borderAnchorLin] = min(handoffCost(:));
else
    bridgeCost = DS + DT;
    bridgeCost(~footMask | ~isfinite(bridgeCost)) = inf;
    [~, borderAnchorLin] = min(bridgeCost(:));
end
end

function [pathXY, pathDelaySec, supportMask, srcPlotAnchorLin, tgtPlotAnchorLin] = localExtractCorridorPath(footMask, delaySec, srcMask, tgtMask, srcAnchorLin, borderAnchorLin, tgtAnchorLin, smoothSigma, interpScale, turnPenalty, edgeBandFrac, corridorRadiusPx, corridorCostSlack, corridorAdaptive, targetWidthWeight)
pathXY = zeros(0,2);
pathDelaySec = zeros(0,1);
supportMask = false(size(footMask));
srcPlotAnchorLin = srcAnchorLin;
tgtPlotAnchorLin = tgtAnchorLin;
srcMask = logical(srcMask) & footMask;
tgtMask = logical(tgtMask) & footMask;
if ~any(srcMask(:)) || ~any(tgtMask(:))
    return
end
if ~(isfinite(srcAnchorLin) && isfinite(tgtAnchorLin))
    return
end

interpScale = max(1, round(interpScale));
if interpScale > 1
    workFoot = imresize(footMask, interpScale, 'nearest');
    workDelay = imresize(delaySec, interpScale, 'bilinear');
    workSrc = imresize(srcMask, interpScale, 'nearest');
    workTgt = imresize(tgtMask, interpScale, 'nearest');
    workSrcAnchor = localScaleLin(size(footMask), srcAnchorLin, interpScale, size(workFoot));
    workTgtAnchor = localScaleLin(size(footMask), tgtAnchorLin, interpScale, size(workFoot));
    workBorderAnchor = localScaleLin(size(footMask), borderAnchorLin, interpScale, size(workFoot));
else
    workFoot = footMask;
    workDelay = delaySec;
    workSrc = srcMask;
    workTgt = tgtMask;
    workSrcAnchor = srcAnchorLin;
    workTgtAnchor = tgtAnchorLin;
    workBorderAnchor = borderAnchorLin;
end

DS = bwdistgeodesic(workFoot, workSrc, 'quasi-euclidean');
DT = bwdistgeodesic(workFoot, workTgt, 'quasi-euclidean');
if ~any(isfinite(DS(:))) || ~any(isfinite(DT(:)))
    return
end

srcVals = workDelay(workSrc & isfinite(workDelay));
tgtVals = workDelay(workTgt & isfinite(workDelay));
srcDelay = mean(srcVals, 'omitnan');
tgtDelay = mean(tgtVals, 'omitnan');
if ~(isfinite(srcDelay) && isfinite(tgtDelay))
    srcDelay = 0;
    tgtDelay = max(workDelay(isfinite(workDelay)), [], 'omitnan');
end

baseCost = 1 + 0.20 * (DS + DT);
delaySpan = max(abs(tgtDelay - srcDelay), eps);
frac = DS ./ max(DS + DT, eps);
expected = srcDelay + frac * (tgtDelay - srcDelay);
delayPenalty = abs(workDelay - expected) ./ delaySpan;
delayPenalty(~isfinite(delayPenalty)) = 1;
dtInterior = bwdist(~workFoot);
maxInterior = max(dtInterior(workFoot), [], 'omitnan');
if ~(isfinite(maxInterior) && maxInterior > 0)
    maxInterior = 1;
end
interiorPenalty = 1 ./ max(dtInterior, 1.5);
edgeFrac = dtInterior ./ maxInterior;
edgeBandMask = workFoot & (edgeFrac <= edgeBandFrac);
weightImg = baseCost + 2.0 * delayPenalty + 2.5 * interiorPenalty + 6.0 * double(edgeBandMask);
weightImg(~workFoot | ~isfinite(weightImg)) = inf;
if smoothSigma > 0
    smoothWeight = localSmoothCostField(weightImg, workFoot, smoothSigma * interpScale / max(interpScale, 1));
else
    smoothWeight = weightImg;
end

Dsrc = bwdistgeodesic(workFoot, workSrc, 'quasi-euclidean');
Dtgt = bwdistgeodesic(workFoot, workTgt, 'quasi-euclidean');
if ~any(isfinite(Dsrc(:))) || ~any(isfinite(Dtgt(:)))
    return
end

joint = Dsrc + Dtgt;
bridgeVals = joint(workTgt & isfinite(joint));
if isempty(bridgeVals)
    bridgeVals = joint(isfinite(joint));
end
bestTotal = min(bridgeVals, [], 'omitnan');
if ~(isfinite(bestTotal) && bestTotal > 0)
    return
end
supportBand = workFoot & isfinite(Dsrc) & isfinite(Dtgt);
supportBand = supportBand & (joint <= bestTotal * 1.35);
supportBand = supportBand & (delayPenalty <= 0.95);
supportBand(workSrc | workTgt) = true;
if isfinite(workBorderAnchor) && workBorderAnchor >= 1 && workBorderAnchor <= numel(supportBand)
    supportBand(workBorderAnchor) = true;
end
supportBand = imclose(supportBand, strel('disk', max(1, round(3 * interpScale / 2)), 0));
supportBand = imfill(supportBand, 'holes');
supportBand = localKeepComponentTouchingMasks(supportBand, workSrc, workTgt);
routeMask = supportBand & ~edgeBandMask;
routeMask(workSrc | workTgt) = true;
routeMask = localKeepComponentTouchingMasks(routeMask, workSrc, workTgt);
[pathLinTest, srcPathLin, tgtPathLin] = localShortestPathBetweenMasks(routeMask, smoothWeight, workSrc, workTgt, turnPenalty);
if isempty(pathLinTest)
    routeMask = supportBand;
    [pathLinTest, srcPathLin, tgtPathLin] = localShortestPathBetweenMasks(routeMask, smoothWeight, workSrc, workTgt, turnPenalty);
end
if isempty(pathLinTest)
    return
end
srcPlotAnchorLin = localScaleBackLin(srcPathLin, size(workFoot), interpScale, size(footMask));
tgtPlotAnchorLin = localScaleBackLin(tgtPathLin, size(workFoot), interpScale, size(footMask));
[pathXY, pathDelaySec, supportMask] = localExtractInterpolatedBandPath(routeMask, smoothWeight, workDelay, workSrc, workTgt, pathLinTest, interpScale, size(footMask), corridorRadiusPx, corridorCostSlack, turnPenalty, corridorAdaptive, targetWidthWeight);
end

function smoothWeight = localSmoothCostField(weightImg, footMask, smoothSigma)
fillVal = max(weightImg(isfinite(weightImg)), [], 'omitnan');
if ~(isfinite(fillVal) && fillVal > 0)
    fillVal = 1;
end
work = weightImg;
work(~isfinite(work)) = fillVal * 2;
work = imgaussfilt(work, smoothSigma, 'FilterDomain', 'spatial');
smoothWeight = work;
smoothWeight(~footMask) = inf;
end

function [pathXY, pathDelaySec, supportMaskOut] = localExtractInterpolatedBandPath(supportBand, smoothWeight, delaySec, srcMask, tgtMask, pathLinFallback, interpScale, baseSize, corridorRadiusPx, corridorCostSlack, turnPenalty, corridorAdaptive, targetWidthWeight)
supportMaskOut = localExpandCorridorRegion(supportBand, smoothWeight, srcMask, tgtMask, corridorRadiusPx, corridorCostSlack, corridorAdaptive, targetWidthWeight);
pathXY = zeros(0,2);
pathDelaySec = zeros(0,1);
pathLin = localCenterlineFromSupport(supportMaskOut, smoothWeight, srcMask, tgtMask, turnPenalty);
if isempty(pathLin)
    pathLin = pathLinFallback;
    if isempty(pathLin)
        return
    end
end
[r,c] = ind2sub(size(supportMaskOut), pathLin);
pathXY = [c, r];
pathDelaySec = delaySec(pathLin);
if interpScale > 1
    pathXY = [(pathXY(:,1) - 0.5) ./ interpScale + 0.5, (pathXY(:,2) - 0.5) ./ interpScale + 0.5];
    pathXY = localUniqueRowsRounded(pathXY);
    supportMaskOut = imresize(supportMaskOut, baseSize, 'nearest');
else
    supportMaskOut = supportMaskOut;
end
end

function pathLin = localShortestPathOnMask(mask, weightImg, srcLin, tgtLin, turnPenalty)
pathLin = zeros(0,1);
if ~(isfinite(srcLin) && isfinite(tgtLin))
    return
end
seedLin = localNearestTruePixel(mask, srcLin);
targetLin = localNearestTruePixel(mask, tgtLin);
if ~(isfinite(seedLin) && isfinite(targetLin))
    return
end
if nargin < 5 || ~isfinite(turnPenalty)
    turnPenalty = 0;
end
pathLin = localDirectionalShortestPath(mask, weightImg, seedLin, targetLin, turnPenalty);
end

function linOut = localScaleLin(szIn, linIn, scale, szOut)
linOut = NaN;
if ~(isfinite(linIn) && linIn >= 1 && linIn <= prod(szIn))
    return
end
[r,c] = ind2sub(szIn, linIn);
r2 = min(max(round((r - 0.5) * scale + 0.5), 1), szOut(1));
c2 = min(max(round((c - 0.5) * scale + 0.5), 1), szOut(2));
linOut = sub2ind(szOut, r2, c2);
end

function linOut = localScaleBackLin(linIn, szIn, scale, szOut)
linOut = NaN;
if ~(isfinite(linIn) && linIn >= 1 && linIn <= prod(szIn))
    return
end
if scale <= 1
    linOut = linIn;
    return
end
[r,c] = ind2sub(szIn, linIn);
r2 = min(max(round((r - 0.5) ./ scale + 0.5), 1), szOut(1));
c2 = min(max(round((c - 0.5) ./ scale + 0.5), 1), szOut(2));
linOut = sub2ind(szOut, r2, c2);
end

function mask = localKeepComponentTouchingMasks(mask, srcMask, tgtMask)
if ~any(mask(:))
    return
end
cc = bwconncomp(mask, 8);
keep = false(size(mask));
for i = 1:cc.NumObjects
    pix = cc.PixelIdxList{i};
    if any(srcMask(pix)) && any(tgtMask(pix))
        keep(pix) = true;
    end
end
if any(keep(:))
    mask = keep;
end
end

function [pathLin, srcLin, tgtLin] = localShortestPathBetweenMasks(mask, weightImg, srcMask, tgtMask, turnPenalty)
pathLin = zeros(0,1);
srcLin = NaN;
tgtLin = NaN;
srcCandidates = find(mask & srcMask);
tgtCandidates = find(mask & tgtMask);
if isempty(srcCandidates) || isempty(tgtCandidates)
    return
end
D = graydist(weightImg, mask & srcMask, 'quasi-euclidean');
tgtValid = tgtCandidates(isfinite(D(tgtCandidates)));
if isempty(tgtValid)
    return
end
[~,k] = min(D(tgtValid));
tgtLin = tgtValid(k);
pathLin = localDirectionalShortestPathToMask(mask, weightImg, srcCandidates, tgtLin, turnPenalty);
if isempty(pathLin)
    return
end
srcLin = pathLin(1);
end

function pathLin = localDirectionalShortestPathToMask(mask, weightImg, srcCandidates, targetLin, turnPenalty)
pathLin = zeros(0,1);
[rows, cols] = find(mask);
if isempty(rows)
    return
end
linList = sub2ind(size(mask), rows, cols);
linToNode = zeros(numel(mask),1);
linToNode(linList) = 1:numel(linList);
targetNode = linToNode(targetLin);
if targetNode == 0
    return
end
seedNodes = linToNode(srcCandidates);
seedNodes = seedNodes(seedNodes > 0);
if isempty(seedNodes)
    return
end
ndir = 8;
offRC = [-1 0; -1 1; 0 1; 1 1; 1 0; 1 -1; 0 -1; -1 -1];
ndirState = ndir + 1;
nState = numel(linList) * ndirState;
dist = inf(nState,1);
prev = zeros(nState,1);
visited = false(nState,1);
for i = 1:numel(seedNodes)
    startState = seedNodes(i) * ndirState;
    dist(startState) = 0;
end
for iter = 1:nState
    cand = find(~visited & isfinite(dist));
    if isempty(cand)
        break
    end
    [~,k] = min(dist(cand));
    s = cand(k);
    visited(s) = true;
    node = ceil(s / ndirState);
    dirPrev = s - (node-1)*ndirState - 1;
    if node == targetNode
        break
    end
    r = rows(node); c = cols(node);
    for d = 1:ndir
        rr = r + offRC(d,1); cc = c + offRC(d,2);
        if rr < 1 || rr > size(mask,1) || cc < 1 || cc > size(mask,2)
            continue
        end
        linN = sub2ind(size(mask), rr, cc);
        nextNode = linToNode(linN);
        if nextNode == 0
            continue
        end
        stepLen = hypot(offRC(d,1), offRC(d,2));
        base = stepLen * weightImg(linN);
        if ~isfinite(base)
            continue
        end
        turn = 0;
        if dirPrev > 0
            delta = abs(d - dirPrev);
            delta = min(delta, ndir - delta);
            turn = turnPenalty * (delta / 2);
        end
        s2 = (nextNode-1)*ndirState + (d+1);
        alt = dist(s) + base + turn;
        if alt < dist(s2)
            dist(s2) = alt;
            prev(s2) = s;
        end
    end
end
endStates = (targetNode-1)*ndirState + (2:ndirState);
[bestDist, ibest] = min(dist(endStates));
if ~isfinite(bestDist)
    return
end
s = endStates(ibest);
pathNodes = zeros(0,1);
while s > 0
    node = ceil(s / ndirState);
    pathNodes(end+1,1) = node; %#ok<AGROW>
    s = prev(s);
end
pathNodes = flipud(pathNodes);
pathLin = linList(pathNodes);
pathLin = localUniqueStable(pathLin);
end

function xy = localUniqueRowsRounded(xy)
if isempty(xy)
    return
end
xy = round(xy, 3);
keep = [true; any(diff(xy,1,1) ~= 0, 2)];
xy = xy(keep,:);
end

function D = localGrayDistFromLin(weightImg, seedLin)
seedMask = false(size(weightImg));
seedMask(seedLin) = true;
D = graydist(weightImg, seedMask, 'quasi-euclidean');
end

function supportMask = localExpandCorridorRegion(baseMask, costImg, srcMask, tgtMask, radiusPx, costSlack, adaptive, targetWidthWeight)
supportMask = baseMask;
if ~any(baseMask(:))
    return
end
regionCosts = costImg(baseMask & isfinite(costImg));
if isempty(regionCosts)
    return
end
refCost = median(regionCosts, 'omitnan');
if ~(isfinite(refCost) && refCost > 0)
    refCost = max(regionCosts, [], 'omitnan');
end
targetWidth = localEstimateTargetWidth(baseMask, srcMask, tgtMask);
if ~adaptive
    supportMask = localCorridorCoreAtSlack(baseMask, costImg, srcMask, tgtMask, radiusPx, refCost, costSlack);
    if ~any(supportMask(:))
        supportMask = baseMask;
    end
    return
end
slacks = unique([linspace(costSlack, 1.02, 10), costSlack]);
slacks = sort(slacks, 'descend');
connectedMasks = {};
connectedStats = repmat(localCorridorSimplicity(false(size(baseMask)), srcMask, tgtMask, costImg), 0, 1);
for i = 1:numel(slacks)
    cand = localCorridorCoreAtSlack(baseMask, costImg, srcMask, tgtMask, radiusPx, refCost, slacks(i));
    if ~any(cand(:))
        continue
    end
    stats = localCorridorSimplicity(cand, srcMask, tgtMask, costImg);
    if ~stats.connected
        continue
    end
    connectedMasks{end+1} = cand; %#ok<AGROW>
    connectedStats(end+1,1) = stats; %#ok<AGROW>
end
if isempty(connectedMasks)
    return
end
breakIdx = numel(connectedMasks);
refStats = connectedStats(breakIdx);
lo = max(1, breakIdx - 3);
hi = breakIdx;
useIdx = [];
for i = lo:hi
    st = connectedStats(i);
    tooMuch = st.meanWidth > max(4.0, 1.75 * refStats.meanWidth) || st.area > max(refStats.area + 150, 1.8 * refStats.area);
    tooLittle = st.pathLength < max(8, 0.7 * refStats.pathLength) || st.meanWidth < 0.6 * refStats.meanWidth;
    widthErr = abs(st.meanWidth - targetWidth);
    widthOk = widthErr <= max(1.5, 0.6 * targetWidth);
    if ~(tooMuch || tooLittle) && widthOk
        useIdx(end+1) = i; %#ok<AGROW>
    end
end
if isempty(useIdx)
    scores = inf(numel(connectedMasks),1);
    for i = 1:numel(connectedMasks)
        st = connectedStats(i);
        scores(i) = st.branchPoints * 50 + st.area + targetWidthWeight * abs(st.meanWidth - targetWidth);
    end
    [~, bestIdx] = min(scores);
    useIdx = max(1,bestIdx-1):min(numel(connectedMasks),bestIdx+1);
end
accum = zeros(size(baseMask));
for i = useIdx
    accum = accum + double(connectedMasks{i});
end
voteThresh = max(1, ceil(0.5 * numel(useIdx)));
consensus = accum >= voteThresh;
consensus(srcMask | tgtMask) = true;
consensus = imopen(consensus, strel('disk', 1, 0));
consensus = imclose(consensus, strel('disk', max(1, round(radiusPx/4)), 0));
consensus = imfill(consensus, 'holes');
consensus = localKeepComponentTouchingMasks(consensus, srcMask, tgtMask);
if any(consensus(:))
    supportMask = consensus;
else
    supportMask = connectedMasks{breakIdx};
end
end

function targetWidth = localEstimateTargetWidth(baseMask, srcMask, tgtMask)
d = bwdist(~baseMask);
vals = [];
if any(baseMask & srcMask, 'all')
    vals = [vals; 2*d(baseMask & srcMask)]; %#ok<AGROW>
end
if any(baseMask & tgtMask, 'all')
    vals = [vals; 2*d(baseMask & tgtMask)]; %#ok<AGROW>
end
bridge = baseMask & ~srcMask & ~tgtMask;
if any(bridge(:))
    vals = [vals; 2*d(bridge)]; %#ok<AGROW>
end
vals = vals(isfinite(vals) & vals > 0);
if isempty(vals)
    targetWidth = 4;
    return
end
q25 = prctile(vals,25);
q50 = prctile(vals,50);
targetWidth = max(3, min(12, 0.6*q25 + 0.4*q50));
end

function mask = localCorridorCoreAtSlack(baseMask, costImg, srcMask, tgtMask, radiusPx, refCost, slack)
mask = false(size(baseMask));
costCore = baseMask & isfinite(costImg) & (costImg <= refCost * slack);
if ~any(costCore(:))
    return
end
mask = costCore;
mask = imclose(mask, strel('disk', max(1, round(radiusPx/3)), 0));
mask = imfill(mask, 'holes');
mask(srcMask | tgtMask) = true;
mask = localKeepComponentTouchingMasks(mask, srcMask, tgtMask);
end

function stats = localCorridorSimplicity(mask, srcMask, tgtMask, costImg)
stats = struct('connected', false, 'area', 0, 'branchPoints', inf, 'pathLength', 0, 'meanWidth', 0, 'meanCost', inf);
if ~any(mask(:))
    return
end
stats.area = nnz(mask);
pathLin = localCenterlineFromSupport(mask, costImg, srcMask, tgtMask, 0);
stats.connected = ~isempty(pathLin);
if ~stats.connected
    return
end
stats.pathLength = numel(pathLin);
dIn = bwdist(~mask);
skel = bwmorph(mask, 'skel', Inf);
skel = bwmorph(skel, 'spur', 3);
widthVals = 2 * dIn(skel);
widthVals = widthVals(isfinite(widthVals) & widthVals > 0);
if isempty(widthVals)
    widthVals = 2 * dIn(mask);
    widthVals = widthVals(isfinite(widthVals) & widthVals > 0);
end
if isempty(widthVals)
    stats.meanWidth = 0;
else
    stats.meanWidth = median(widthVals, 'omitnan');
end
if ~isempty(pathLin)
    stats.meanCost = mean(costImg(pathLin), 'omitnan');
end
branchMask = bwmorph(skel, 'branchpoints');
stats.branchPoints = nnz(branchMask);
end

function m = pathLinMask(mask, pathLin)
m = false(size(mask));
if ~isempty(pathLin)
    m(pathLin) = true;
end
end

function frac = localInteriorReachFrac(supportMask, cellMask)
frac = 0;
cellMask = logical(cellMask);
if ~any(cellMask(:))
    return
end
d = bwdist(~cellMask);
vals = d(cellMask);
vals = vals(isfinite(vals) & vals > 0);
if isempty(vals)
    return
end
maxd = max(vals);
if ~(isfinite(maxd) && maxd > 0)
    return
end
innerMask = cellMask & (d >= 0.5 * maxd);
if ~any(innerMask(:))
    innerMask = cellMask & (d >= 0.35 * maxd);
end
if any(innerMask(:))
    frac = nnz(supportMask & innerMask) / nnz(innerMask);
end
end

function frac = localTargetSpillFrac(supportMask, srcMask, tgtMask)
frac = 0;
tgtOnly = supportMask & tgtMask;
if ~any(tgtOnly(:))
    return
end
bridgeSeeds = imdilate(srcMask, strel('disk', 3, 0)) & tgtMask;
if ~any(bridgeSeeds(:))
    bridgeSeeds = tgtMask;
end
Dt = bwdistgeodesic(tgtMask, bridgeSeeds, 'quasi-euclidean');
allVals = Dt(tgtMask & isfinite(Dt));
if isempty(allVals)
    return
end
coreTh = min(prctile(allVals, 20), 4.0);
if ~(isfinite(coreTh) && coreTh > 0)
    coreTh = min(max(allVals), 4.0);
end
entryCore = tgtMask & isfinite(Dt) & (Dt <= coreTh);
if ~any(entryCore(:))
    frac = 1;
    return
end
farBand = tgtMask & isfinite(Dt) & (Dt > max(coreTh, 2.0));
if ~any(farBand(:))
    frac = nnz(tgtOnly & ~entryCore) / nnz(tgtOnly);
    return
end
frac = 0.5 * (nnz(tgtOnly & ~entryCore) / nnz(tgtOnly)) + 0.5 * (nnz(tgtOnly & farBand) / nnz(farBand));
end

function frac = localMaskOccupancyFrac(supportMask, cellMask)
frac = 0;
cellMask = logical(cellMask);
if ~any(cellMask(:))
    return
end
frac = nnz(supportMask & cellMask) / nnz(cellMask);
end

function pathLin = localCenterlineFromSupport(supportMask, weightImg, srcMask, tgtMask, turnPenalty)
pathLin = zeros(0,1);
if ~any(supportMask(:))
    return
end
skel = bwmorph(supportMask, 'skel', Inf);
skel = bwmorph(skel, 'spur', 3);
skel(srcMask | tgtMask) = true;
skel = localKeepComponentTouchingMasks(skel, srcMask, tgtMask);
[pathLin, ~, ~] = localShortestPathBetweenMasks(skel, weightImg, srcMask, tgtMask, turnPenalty);
end

function pathLin = localDirectionalShortestPath(mask, weightImg, seedLin, targetLin, turnPenalty)
pathLin = zeros(0,1);
[rows, cols] = find(mask);
if isempty(rows)
    return
end
linList = sub2ind(size(mask), rows, cols);
linToNode = zeros(numel(mask),1);
linToNode(linList) = 1:numel(linList);
seedNode = linToNode(seedLin);
targetNode = linToNode(targetLin);
if seedNode == 0 || targetNode == 0
    return
end
ndir = 8;
offRC = [-1 0; -1 1; 0 1; 1 1; 1 0; 1 -1; 0 -1; -1 -1];
ndirState = ndir + 1;
nState = numel(linList) * ndirState;
dist = inf(nState,1);
prev = zeros(nState,1);
visited = false(nState,1);
startState = seedNode * ndirState;
dist(startState) = 0;
for iter = 1:nState
    active = find(~visited & isfinite(dist), 1);
    if isempty(active)
        break
    end
    cand = find(~visited & isfinite(dist));
    [~,k] = min(dist(cand));
    s = cand(k);
    visited(s) = true;
    node = ceil(s / ndirState);
    dirPrev = s - (node-1)*ndirState;
    dirPrev = dirPrev - 1;
    if node == targetNode
        break
    end
    r = rows(node); c = cols(node);
    for d = 1:ndir
        rr = r + offRC(d,1);
        cc = c + offRC(d,2);
        if rr < 1 || rr > size(mask,1) || cc < 1 || cc > size(mask,2)
            continue
        end
        linN = sub2ind(size(mask), rr, cc);
        nextNode = linToNode(linN);
        if nextNode == 0
            continue
        end
        stepLen = hypot(offRC(d,1), offRC(d,2));
        base = stepLen * weightImg(linN);
        if ~isfinite(base)
            continue
        end
        turn = 0;
        if dirPrev > 0
            delta = abs(d - dirPrev);
            delta = min(delta, ndir - delta);
            turn = turnPenalty * (delta / 2);
        end
        s2 = (nextNode-1)*ndirState + (d+1);
        alt = dist(s) + base + turn;
        if alt < dist(s2)
            dist(s2) = alt;
            prev(s2) = s;
        end
    end
end
endStates = (targetNode-1)*ndirState + (2:ndirState);
[bestDist, ibest] = min(dist(endStates));
if ~isfinite(bestDist)
    return
end
s = endStates(ibest);
pathNodes = zeros(0,1);
while s > 0
    node = ceil(s / ndirState);
    pathNodes(end+1,1) = node; %#ok<AGROW>
    s = prev(s);
end
pathNodes = flipud(pathNodes);
pathLin = linList(pathNodes);
pathLin = localUniqueStable(pathLin);
end

function mask = localKeepConnectedComponent(mask, anchorLin)
if ~(isfinite(anchorLin) && anchorLin >= 1 && anchorLin <= numel(mask) && mask(anchorLin))
    return
end
cc = bwconncomp(mask, 8);
keep = false(size(mask));
for i = 1:cc.NumObjects
    if any(cc.PixelIdxList{i} == anchorLin)
        keep(cc.PixelIdxList{i}) = true;
        mask = keep;
        return
    end
end
end

function lin = localNearestTruePixel(mask, anchorLin)
lin = NaN;
if ~(isfinite(anchorLin) && anchorLin >= 1 && anchorLin <= numel(mask))
    return
end
if mask(anchorLin)
    lin = anchorLin;
    return
end
[ra, ca] = ind2sub(size(mask), anchorLin);
[r, c] = find(mask);
if isempty(r)
    return
end
[~, k] = min((r - ra).^2 + (c - ca).^2);
lin = sub2ind(size(mask), r(k), c(k));
end

function xy = localLinToXY(sz, lin)
if ~(isfinite(lin) && lin >= 1 && lin <= prod(sz))
    xy = [NaN NaN];
    return
end
[r,c] = ind2sub(sz, lin);
xy = [c r];
end

function vals = localUniqueStable(vals)
if isempty(vals)
    return
end
keep = [true; diff(vals(:)) ~= 0];
vals = vals(keep);
end

function nbr = localNeighborLin(sz, r, c)
rr = max(1, r-1):min(sz(1), r+1);
cc = max(1, c-1):min(sz(2), c+1);
[CC, RR] = meshgrid(cc, rr);
nbr = sub2ind(sz, RR(:), CC(:));
nbr = nbr(~(RR(:) == r & CC(:) == c));
end

function hFig = localPlotPairCorridor(out, visible)
hFig = figure('Name', 'EventPairCorridor', 'Color', 'w', 'NumberTitle', 'off', 'Visible', visible);
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
    hIm = imagesc(ax, [0 sy*spatialRes], [0 sx*spatialRes], delayMap, 'AlphaData', 0.85 * valid);
    set(hIm, 'AlphaData', 0.85 * valid);
    colormap(ax, turbo(256));
    clim(ax, [0 max(delayMap(valid), [], 'omitnan') + eps]);
    cb = colorbar(ax);
    ylabel(cb, 'Post-Origin First Activation Delay (s)');
end

nodes = out.graph.nodes;
for i = [out.source_idx out.target_idx]
    b = nodes.boundary_xy_px{i};
    if isempty(b), continue; end
    col = [1 1 1];
    if i == out.source_idx, col = [0.2 1 0.2]; end
    if i == out.target_idx, col = [1 0.2 0.2]; end
    plot(ax, b(:,1)*spatialRes, b(:,2)*spatialRes, '-', 'Color', col, 'LineWidth', 2.0);
    text(ax, nodes.centroid_x_um(i), nodes.centroid_y_um(i), char(nodes.label(i)), ...
        'Color', col, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

if any(out.support_mask(:))
    xvec = (0:sy-1) * spatialRes;
    yvec = (0:sx-1) * spatialRes;
    [XX, YY] = meshgrid(xvec, yvec);
    contour(ax, XX, YY, double(out.support_mask), [1 1], ...
        'Color', [1 1 1], 'LineWidth', 1.2);
end
if ~isempty(out.path_xy_px)
    plot(ax, out.path_xy_px(:,1)*spatialRes, out.path_xy_px(:,2)*spatialRes, '-', ...
        'Color', [1 1 1], 'LineWidth', 2.2);
end
if all(isfinite(out.source_anchor_xy_px))
    plot(ax, out.source_anchor_xy_px(1)*spatialRes, out.source_anchor_xy_px(2)*spatialRes, 'o', ...
        'Color', [0.2 1 0.2], 'MarkerFaceColor', [0.2 1 0.2], 'MarkerSize', 6);
end
if all(isfinite(out.target_anchor_xy_px))
    plot(ax, out.target_anchor_xy_px(1)*spatialRes, out.target_anchor_xy_px(2)*spatialRes, 's', ...
        'Color', [1 0.2 0.2], 'MarkerFaceColor', [1 0.2 0.2], 'MarkerSize', 6);
end
title(ax, sprintf('%s | %s -> %s', char(out.event_uid), char(out.source_cell), char(out.target_cell)));
hold(ax, 'off');
end
