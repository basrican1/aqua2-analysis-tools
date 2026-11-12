function out = plot_local_event_motion(source, eventID, varargin)
%PLOT_LOCAL_EVENT_MOTION Visualize centroid and wavefront tracks for one event.
%   OUT = PLOT_LOCAL_EVENT_MOTION(RESULT, EVENTID)
%   OUT = PLOT_LOCAL_EVENT_MOTION(PATHIN, EVENTID)
%   OUT = PLOT_LOCAL_EVENT_MOTION(BWA_AQ2, EVENTID)

p = inputParser;
p.StructExpand = false;
addParameter(p, 'backgroundMode', 'event_progression', @(x) any(strcmpi(string(x), ["none","event_footprint","event_progression"])));
addParameter(p, 'progressFrame', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'showFrameLabels', false, @islogical);
addParameter(p, 'figureName', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'wavefrontMode', 'local', @(x) any(strcmpi(string(x), ["local","prototype"])));
addParameter(p, 'eventUID', "", @(x) ischar(x) || isstring(x));
parse(p, varargin{:});

[BWA_AQ2, pathLabel] = resolveLocalEventMotionSource(source);
if ~isstruct(BWA_AQ2) || ~isfield(BWA_AQ2, 'fts1') || ~isfield(BWA_AQ2, 'opts')
    error('plot_local_event_motion:MissingLegacyPayload', ...
        'This visualizer requires BWA_AQ2 with fts1 and opts.');
end

fts = BWA_AQ2.fts1;
opts = BWA_AQ2.opts;
if strlength(string(p.Results.eventUID)) > 0
    eventID = localEventIndexFromUIDStatic(BWA_AQ2, string(p.Results.eventUID));
end
[imageSz, spatialRes] = localEventImagingMeta(BWA_AQ2, eventID);
geomAll = CalcEventFrameGeometry(fts.loc.xSpaTemp, fts.loc.t0, fts.loc.t1, imageSz, 'eventIDs', eventID);
G = geomAll{eventID};
W = CalcWavefrontMetrics(G, spatialRes);
if strcmpi(string(p.Results.wavefrontMode), "prototype") && exist('wavefrontComponentGraphPrototype', 'file') == 2
    proto = wavefrontComponentGraphPrototype(BWA_AQ2, eventID, ...
        'showFigure', false, ...
        'refinePixelWalk', 'on');
    W = applyPrototypeWavefrontStatic(W, proto);
end
cent = localAdaptiveCentroidForDisplayStatic(G, W);

figName = char(string(p.Results.figureName));
if isempty(strtrim(figName))
    figName = sprintf('LocalEventMotion_%d', eventID);
end
hFig = figure('Name', figName);
ax = axes('Parent', hFig); hold(ax, 'on'); axis(ax, 'equal'); grid(ax, 'on');
xlabel(ax, 'x (um)'); ylabel(ax, 'y (um)');
title(ax, sprintf('Event %d | branch origin=%d | branch peak=%d | reach=%.3f um | frontier=%.3f um | tortuosity=%.3f', ...
    eventID, round(W.branchOriginFrame), round(W.branchPeakFrame), ...
    W.wavefront_max_reach_um, W.wavefront_frontier_length_um, W.wavefront_tortuosity));

[bgImg, bgRef] = buildLocalMotionBackground(fts, G, imageSz, spatialRes, eventID, p.Results.backgroundMode, p.Results.progressFrame, W);
bgHandle = gobjects(0);
if ~isempty(bgImg)
    bgHandle = imshow(bgImg, bgRef, 'Parent', ax);
    uistack(bgHandle, 'bottom');
end

polyHandles = gobjects(0);
frameColors = parula(max(nnz(G.validFrontier), 2));
validIdx = find(G.validFrontier);
for ii = 1:numel(validIdx)
    k = validIdx(ii);
    pts = G.frontierPts{k};
    if isempty(pts)
        continue
    end
    pts = pts * spatialRes;
    c = frameColors(min(ii, size(frameColors,1)), :);
    polyHandles(end+1,1) = plot(ax, pts(:,1), pts(:,2), '.', 'Color', c, 'MarkerSize', 8); %#ok<AGROW>
    if p.Results.showFrameLabels && ~isempty(pts)
        ctr = G.newCentroidXY(k, :) * spatialRes;
        text(ax, ctr(1), ctr(2), sprintf(' %d', G.frames(k)), 'Color', c, 'FontSize', 8);
    end
end

centroidHandle = gobjects(0);
if isstruct(cent) && isfield(cent, 'cx') && ~isempty(cent.cx)
    centKeep = cent.ct >= W.branchOriginFrame;
    if any(centKeep)
        centPts = [W.originXY; [cent.cx(centKeep), cent.cy(centKeep)]];
        centFrames = [W.branchOriginFrame; cent.ct(centKeep)];
        keepUnique = [true; diff(centFrames) > 0 | any(diff(centPts,1,1) ~= 0, 2)];
        centPts = centPts(keepUnique, :);
        centroidHandle = plot(ax, centPts(:,1) * spatialRes, centPts(:,2) * spatialRes, '-o', ...
            'Color', [0.1 0.8 0.8], 'LineWidth', 1.5, 'MarkerSize', 3, ...
            'MarkerFaceColor', [0.1 0.8 0.8]);
    end
end

originHandle = plot(ax, W.originXY(1) * spatialRes, W.originXY(2) * spatialRes, 'ks', ...
    'MarkerFaceColor', 'y', 'MarkerSize', 8, 'LineWidth', 1.25);
peakHandle = plot(ax, W.peakPointXY(1) * spatialRes, W.peakPointXY(2) * spatialRes, 'rp', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 12, 'LineWidth', 1.5);
secondaryHandles = gobjects(0);
if isfield(W, 'secondaryTracks') && ~isempty(W.secondaryTracks)
    for i = 1:numel(W.secondaryTracks)
        T = W.secondaryTracks{i};
        if ~isstruct(T) || ~isfield(T, 'trackValid') || ~isfield(T, 'trackXY')
            continue
        end
        secIdx = find(T.trackValid & W.frames >= W.branchOriginFrame);
        if isempty(secIdx)
            continue
        end
        secPts = T.trackXY(secIdx, :) * spatialRes;
        hSec = plotTrackRuns(ax, W.frames(secIdx), secPts, [0.55 0.55 0.55]);
        if ~isempty(hSec)
            secondaryHandles(end+1,1) = hSec(1); %#ok<AGROW>
        end
    end
end
trackHandle = gobjects(0);
trackIdx = find(W.trackValid & W.frames >= W.branchOriginFrame);
if ~isempty(trackIdx)
    trackPts = W.trackXY(trackIdx, :) * spatialRes;
    trackHandle = plotTrackRuns(ax, W.frames(trackIdx), trackPts, [0.85 0.2 0.2]);
end

legendItems = gobjects(0);
legendLabels = {};
if ~isempty(centroidHandle)
    legendItems(end+1,1) = centroidHandle(1); %#ok<AGROW>
    legendLabels{end+1} = 'centroid track'; %#ok<AGROW>
end
if ~isempty(secondaryHandles)
    legendItems(end+1,1) = secondaryHandles(1); %#ok<AGROW>
    legendLabels{end+1} = 'secondary branch'; %#ok<AGROW>
end
legendItems(end+1:end+2,1) = [originHandle; peakHandle];
legendLabels(end+1:end+2) = {'branch origin','branch peak'};
if ~isempty(trackHandle)
    legendItems(end+1,1) = trackHandle(1); %#ok<AGROW>
    legendLabels{end+1} = 'wavefront track'; %#ok<AGROW>
end
legend(ax, legendItems, legendLabels, 'Location', 'bestoutside');

out = struct('figure', hFig, 'axes', ax, 'background', bgHandle, 'geometry', G, ...
    'centroid', cent, 'wavefront', W, 'sourceLabel', pathLabel);
end

function W = applyPrototypeWavefrontStatic(W, proto)
if ~isstruct(proto)
    return
end
copyFields = {'branchOriginXY','branchOriginFrame','originXY','branchPeakPointXY','branchPeakFrame','peakPointXY','peakFrame', ...
    'wavefront_max_reach_px','wavefront_max_reach_um','wavefront_frontier_length_px','wavefront_frontier_length_um', ...
    'wavefront_net_displacement_px','wavefront_net_displacement_um','wavefront_tortuosity'};
for i = 1:numel(copyFields)
    if isfield(proto.wavefront, copyFields{i})
        W.(copyFields{i}) = proto.wavefront.(copyFields{i});
    end
end
if isfield(proto, 'trackXY') && ~isempty(proto.trackXY)
    W.trackXY = proto.trackXY;
end
if isfield(proto, 'trackValid') && ~isempty(proto.trackValid)
    W.trackValid = proto.trackValid;
end
if isfield(proto, 'secondary') && ~isempty(proto.secondary)
    W.prototypeSecondary = proto.secondary;
end
end

function cent = localAdaptiveCentroidForDisplayStatic(G, W)
cent = struct('cx', G.centroidXY(:,1), 'cy', G.centroidXY(:,2), 'ct', G.frames(:));
if isempty(G) || ~isstruct(G) || isempty(G.frames)
    return
end

nFrames = numel(G.frames);
branchCx = NaN(nFrames, 1);
branchCy = NaN(nFrames, 1);
globalArea = zeros(nFrames, 1);
branchArea = zeros(nFrames, 1);
for k = 1:nFrames
    mask = logical(G.frameMasks{k});
    globalArea(k) = nnz(mask);
    refPt = W.originXY;
    if k > 1
        refPt = [cent.cx(k - 1), cent.cy(k - 1)];
    end
    [compCtr, compArea, ok] = localStableCentroidComponentStatic(mask, refPt);
    if ok
        branchArea(k) = compArea;
        branchCx(k) = compCtr(1);
        branchCy(k) = compCtr(2);
    end
end

lowArea = 40;
highArea = 250;
maxSparseJumpPx = 18;
warmupJumpPx = 6;
warmupFrames = 4;
idx0 = find(cent.ct == W.branchOriginFrame, 1, 'first');
if isempty(idx0)
    idx0 = 1;
end
cent.cx(idx0) = W.originXY(1);
cent.cy(idx0) = W.originXY(2);
for k = 1:nFrames
    if k < idx0
        continue
    end
    if k == idx0
        continue
    end
    if ~isfinite(branchCx(k)) || ~isfinite(branchCy(k))
        continue
    end
    areaTot = globalArea(k);
    frac = branchArea(k) / max(areaTot, 1);
    areaScore = min(max((areaTot - lowArea) / max(highArea - lowArea, 1), 0), 1);
    fracScore = min(max((frac - 0.35) / 0.45, 0), 1);
    wGlobal = areaScore * fracScore;
    prevUsed = [cent.cx(k - 1), cent.cy(k - 1)];
    cand = wGlobal * [cent.cx(k), cent.cy(k)] + (1 - wGlobal) * [branchCx(k), branchCy(k)];
    jumpCand = hypot(cand(1) - prevUsed(1), cand(2) - prevUsed(2));
    framesFromOrigin = k - idx0;
    maxJumpAllowed = maxSparseJumpPx;
    if framesFromOrigin >= 1 && framesFromOrigin <= warmupFrames
        alpha = (framesFromOrigin - 1) / max(warmupFrames - 1, 1);
        maxJumpAllowed = warmupJumpPx + alpha * (maxSparseJumpPx - warmupJumpPx);
    end
    if areaTot < highArea && jumpCand > maxSparseJumpPx
        cand = [branchCx(k), branchCy(k)];
        jumpComp = hypot(cand(1) - prevUsed(1), cand(2) - prevUsed(2));
        jumpCand = jumpComp;
    end
    if jumpCand > maxJumpAllowed
        direction = cand - prevUsed;
        denom = hypot(direction(1), direction(2));
        if denom > 0
            cand = prevUsed + (maxJumpAllowed / denom) * direction;
        else
            cand = prevUsed;
        end
    end
    cent.cx(k) = cand(1);
    cent.cy(k) = cand(2);
end
end

function [ctr, area, ok] = localStableCentroidComponentStatic(mask, refPt)
ctr = [NaN NaN];
area = 0;
ok = false;
if isempty(mask) || ~any(mask(:))
    return
end
cc = bwconncomp(mask, 8);
bestIdx = NaN;
bestCost = inf;
for i = 1:cc.NumObjects
    compMask = false(size(mask));
    compMask(cc.PixelIdxList{i}) = true;
    [r, c] = find(compMask);
    pts = [c, r];
    compCtr = [mean(c, 'omitnan'), mean(r, 'omitnan')];
    if any(~isfinite(refPt))
        cost = 0;
    else
        dCtr = hypot(compCtr(1) - refPt(1), compCtr(2) - refPt(2));
        dPix = min(hypot(pts(:,1) - refPt(1), pts(:,2) - refPt(2)));
        cost = 0.7 * dPix + 0.3 * dCtr;
    end
    if cost < bestCost
        bestCost = cost;
        bestIdx = i;
        ctr = compCtr;
        area = numel(cc.PixelIdxList{i});
    end
end
ok = isfinite(bestIdx);
end

function h = plotTrackRuns(ax, frames, pts, color)
h = gobjects(0);
if isempty(frames) || isempty(pts)
    return
end
runStart = 1;
for i = 2:(numel(frames) + 1)
    isBreak = (i > numel(frames)) || (frames(i) ~= frames(i - 1) + 1);
    if ~isBreak
        continue
    end
    idx = runStart:(i - 1);
    h(end + 1, 1) = plot(ax, pts(idx,1), pts(idx,2), '-o', ... %#ok<AGROW>
        'Color', color, 'MarkerFaceColor', color, ...
        'LineWidth', 2, 'MarkerSize', 4);
    runStart = i;
end
end

function eventID = localEventIndexFromUIDStatic(BWA_AQ2, eventUID)
eventID = NaN;
fts = BWA_AQ2.fts1;
if isfield(fts, 'basic') && isfield(fts.basic, 'fileID') && ~isempty(fts.basic.fileID)
    fileIDs = string(fts.basic.fileID);
elseif isfield(fts, 'fileID') && ~isempty(fts.fileID)
    fileIDs = string(fts.fileID);
else
    error('plot_local_event_motion:MissingFileIDs', 'Could not resolve eventUID without file IDs.');
end
tok = regexp(char(eventUID), '^(.*)_E(\d+)$', 'tokens', 'once');
if isempty(tok)
    error('plot_local_event_motion:BadEventUID', 'Expected eventUID like <file_id>_E<number>.');
end
fileID = string(tok{1});
localID = str2double(tok{2});
if ~isfinite(localID) || localID < 1
    error('plot_local_event_motion:BadEventUID', 'Invalid eventUID suffix in %s.', char(eventUID));
end
idx = find(fileIDs == fileID);
if numel(idx) < localID
    error('plot_local_event_motion:EventUIDNotFound', 'eventUID %s was not found in loaded BWA_AQ2.', char(eventUID));
end
eventID = idx(localID);
end

function [BWA_AQ2, label] = resolveLocalEventMotionSource(source)
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
if isstruct(source) && isfield(source, 'fts1') && isfield(source, 'opts')
    BWA_AQ2 = source;
    return
end
error('plot_local_event_motion:BadSource', 'Expected result, pathIn, or BWA_AQ2.');
end

function [img, R] = buildLocalMotionBackground(fts, G, imageSz, spatialRes, eventID, modeIn, progressFrame, W)
mode = lower(char(string(modeIn)));
sz = double(imageSz(:)');
sx = round(sz(1)); sy = round(sz(2)); sz3 = 1;
if numel(sz) >= 3 && isfinite(sz(3)) && sz(3) > 0
    sz3 = round(sz(3));
end
img = [];
R = [];
if strcmp(mode, 'none')
    return
end
img = zeros(sx, sy, 3);
footMask = localGeometryFootprintMaskStatic(G, [sx sy]);
if isempty(footMask)
    footMask = false(sx, sy);
    footMask(double(fts.loc.xSpa{eventID})) = true;
end
img = applyMaskColorLocal(img, footMask, [0.72 0.72 0.72]);
if strcmp(mode, 'event_progression')
    frameNow = double(progressFrame);
    if ~isfinite(frameNow)
        frameNow = double(W.peakFrame);
    end
    prevMask = false(sx, sy);
    currMask = false(sx, sy);
    if isstruct(G) && isfield(G, 'frames') && isfield(G, 'frameMasks') && ~isempty(G.frames)
        idxNow = find(double(G.frames(:)) == double(frameNow), 1, 'first');
        if ~isempty(idxNow)
            currMask = localMaskLikeImageSizeStatic(G.frameMasks{idxNow}, [sx sy]);
            idxPrev = max(idxNow - 1, 1);
            prevMask = localMaskLikeImageSizeStatic(G.frameMasks{idxPrev}, [sx sy]);
        end
    else
        t0 = double(fts.loc.t0(eventID));
        framePrev = max(t0, frameNow - 1);
        prevMask = localFrameMask(fts.loc.xSpaTemp{eventID}, framePrev, sx, sy, sz3);
        currMask = localFrameMask(fts.loc.xSpaTemp{eventID}, frameNow, sx, sy, sz3);
    end
    img = applyMaskColorLocal(img, prevMask & ~currMask, [0.35 0.00 0.60]);
    img = applyMaskColorLocal(img, prevMask & currMask, [0.00 0.35 0.90]);
    img = applyMaskColorLocal(img, currMask & ~prevMask, [0.85 0.20 0.80]);
end
R = imref2d([sx sy], spatialRes, spatialRes);
end

function footMask = localGeometryFootprintMaskStatic(G, imageSize2D)
footMask = [];
if ~(isstruct(G) && isfield(G, 'frameMasks') && ~isempty(G.frameMasks))
    return
end
footMask = false(imageSize2D(1), imageSize2D(2));
for k = 1:numel(G.frameMasks)
    Mk = G.frameMasks{k};
    if isempty(Mk)
        continue
    end
    footMask = footMask | logical(Mk);
end
end

function maskOut = localMaskLikeImageSizeStatic(maskIn, imageSize2D)
maskOut = false(imageSize2D(1), imageSize2D(2));
if isempty(maskIn)
    return
end
maskIn = logical(maskIn);
sx = min(size(maskIn, 1), imageSize2D(1));
sy = min(size(maskIn, 2), imageSize2D(2));
maskOut(1:sx, 1:sy) = maskIn(1:sx, 1:sy);
end

function [imageSz, spatialRes] = localEventImagingMeta(BWA_AQ2, eventID)
opts = BWA_AQ2.opts;
imageSz = double(opts.sz(1,:));
spatialRes = double(opts.spatialRes(1));
if ~(isstruct(BWA_AQ2) && isfield(BWA_AQ2, 'fts1') && isfield(BWA_AQ2.fts1, 'fileID'))
    return
end
fileIDs = canonCellstrLocal(BWA_AQ2.fts1.fileID);
if eventID < 1 || eventID > numel(fileIDs)
    return
end
thisID = fileIDs{eventID};
rowIDs = unique(fileIDs, 'stable');
row = find(strcmp(rowIDs, thisID), 1, 'first');
if isempty(row)
    return
end
if size(opts.sz, 1) >= row
    imageSz = double(opts.sz(row, :));
end
if numel(opts.spatialRes) >= row && isfinite(opts.spatialRes(row))
    spatialRes = double(opts.spatialRes(row));
end
if numel(imageSz) < 3
    imageSz(3) = 1;
end
end

function c = canonCellstrLocal(v)
if iscell(v)
    c = cellfun(@char, cellfun(@string, v, 'UniformOutput', false), 'UniformOutput', false);
elseif isstring(v)
    c = cellstr(v(:));
elseif ischar(v)
    c = {v};
else
    c = cellstr(string(v(:)));
end
end

function mask = localFrameMask(xSpaTemp_ID, frameNum, sx, sy, sz3)
limitStart = sx * sy * sz3 * (frameNum - 1) + 1;
limitEnd = limitStart + sx * sy * sz3;
linMovie = double(xSpaTemp_ID(:));
linFrame = linMovie(linMovie >= limitStart & linMovie <= limitEnd) - limitStart + 1;
linFrame = linFrame(isfinite(linFrame) & linFrame >= 1 & linFrame <= sx * sy);
mask = false(sx, sy);
mask(linFrame) = true;
end

function img = applyMaskColorLocal(img, mask, color)
for c = 1:3
    plane = img(:,:,c);
    plane(mask) = color(c);
    img(:,:,c) = plane;
end
end
