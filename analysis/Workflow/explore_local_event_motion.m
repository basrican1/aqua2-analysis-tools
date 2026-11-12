function out = explore_local_event_motion(source, varargin)
%EXPLORE_LOCAL_EVENT_MOTION Interactive centroid + wavefront event viewer.
%   OUT = EXPLORE_LOCAL_EVENT_MOTION(RESULT)
%   OUT = EXPLORE_LOCAL_EVENT_MOTION(RESULT, 'eventID', 40)
%   OUT = EXPLORE_LOCAL_EVENT_MOTION(PATHIN, ...)
%   OUT = EXPLORE_LOCAL_EVENT_MOTION(BWA_AQ2, ...)
%
%   This viewer is modeled on the older centroid explorer, but overlays both
%   centroid and wavefront tracks on top of the evolving event footprint.
%
%   Name-Value options:
%     'eventID'    : initial event ID to show
%     'eventUID'   : exact event UID (overrides eventID)
%     'figureName' : figure title/name
%     'visible'    : 'on' or 'off' (useful for tests and batch smoke runs)
%     'showTracking' : true/false to compute and display centroid/wavefront tracks
%     'wavefrontMode' : 'local' or 'prototype'

p = inputParser;
p.StructExpand = false;
addParameter(p, 'eventID', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'eventUID', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'figureName', 'LocalEventMotionExplorer', @(x) ischar(x) || isstring(x));
addParameter(p, 'visible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
addParameter(p, 'showTracking', false, @islogical);
addParameter(p, 'wavefrontMode', 'local', @(x) any(strcmpi(string(x), ["local","prototype"])));
parse(p, varargin{:});

[BWA_AQ2, sourceLabel, eventList, filesTbl] = resolveLocalEventMotionSourceLocal(source);
if ~isstruct(BWA_AQ2) || ~isfield(BWA_AQ2, 'fts1') || ~isfield(BWA_AQ2, 'opts')
    error('explore_local_event_motion:MissingLegacyPayload', ...
        'This visualizer requires BWA_AQ2 with fts1 and opts.');
end

fts = BWA_AQ2.fts1;
opts = BWA_AQ2.opts;
nEvents = numel(fts.loc.t0);
if isempty(eventList)
    eventList = reshape(1:nEvents, [], 1);
else
    eventList = reshape(unique(round(double(eventList(:))), 'stable'), [], 1);
    eventList = eventList(isfinite(eventList) & eventList >= 1 & eventList <= nEvents);
end
if isempty(eventList)
    error('explore_local_event_motion:EmptyEventSubset', ...
        'No valid events from the supplied result could be matched to BWA_AQ2.');
end
eventID = localResolveInitialEventID(eventList, p.Results.eventID);
if strlength(string(p.Results.eventUID)) > 0
    eventID = localEventIndexFromUID(BWA_AQ2, string(p.Results.eventUID));
    if ~any(eventList == eventID)
        error('explore_local_event_motion:EventUIDOutsideSubset', ...
            'eventUID %s is not present in the supplied result subset.', char(string(p.Results.eventUID)));
    end
end
eventPos = find(eventList == eventID, 1, 'first');

hFig = figure('Name', char(string(p.Results.figureName)), 'Color', 'w', ...
    'NumberTitle', 'off', 'Position', [80 80 900 760], ...
    'Visible', char(string(p.Results.visible)));
hAx = axes('Parent', hFig, 'Units', 'pixels', 'Position', [80 120 620 620]);

uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Event', ...
    'Position', [80 84 40 18], 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
hEvent = uicontrol('Parent', hFig, 'Style', 'slider', ...
    'Min', 1, 'Max', max(numel(eventList), 1), 'Value', eventPos, ...
    'SliderStep', [1/max(max(numel(eventList)-1,1),1) min(10/max(max(numel(eventList)-1,1),1), 1)], ...
    'Position', [120 88 420 18], 'Callback', @onEventChange);
hEventText = uicontrol('Parent', hFig, 'Style', 'text', ...
    'String', localEventLabel(BWA_AQ2, filesTbl, eventPos, numel(eventList), eventID), ...
    'Position', [550 84 330 18], 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');

uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Frame', ...
    'Position', [80 54 40 18], 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
hFrame = uicontrol('Parent', hFig, 'Style', 'slider', ...
    'Min', 1, 'Max', 2, 'Value', 1, 'SliderStep', [1 1], ...
    'Position', [120 58 420 18], 'Callback', @onFrameChange);
hFrameText = uicontrol('Parent', hFig, 'Style', 'text', 'String', '', ...
    'Position', [550 54 120 18], 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');

hPlay = uicontrol('Parent', hFig, 'Style', 'togglebutton', 'String', 'Play', ...
    'Position', [690 54 70 22], 'Callback', @onPlayToggle);
hInfo = uicontrol('Parent', hFig, 'Style', 'text', 'String', char(string(sourceLabel)), ...
    'Position', [80 20 780 18], 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');

state = struct();
state.figure = hFig;
state.axes = hAx;
state.controls = struct('event', hEvent, 'eventText', hEventText, ...
    'frame', hFrame, 'frameText', hFrameText, 'play', hPlay, 'info', hInfo);
state.BWA_AQ2 = BWA_AQ2;
state.fts = fts;
state.opts = opts;
state.filesTbl = filesTbl;
state.imageSz = double(opts.sz(1,:));
state.spatialRes = double(opts.spatialRes(1));
state.nEvents = nEvents;
state.eventList = eventList;
state.nSelectableEvents = numel(eventList);
state.eventID = eventID;
state.showTracking = logical(p.Results.showTracking);
state.wavefrontMode = char(string(p.Results.wavefrontMode));
state.sourceLabel = sourceLabel;
state.sourcePath = localResolveSourcePathExplorer(source, sourceLabel);
state.eventPartition = localResolveEventPartitionTableExplorer(source);
state.cellMaskCache = struct('fileIDs', strings(0,1), 'masks', {{}});
state.geometry = [];
state.wavefront = [];
state.centroid = [];
state.cellOverlay = struct('available', false, 'reason', "", 'cellUID', strings(0,1), ...
    'isPrimary', false(0,1), 'boundaries', {{}}, 'partitionBoundaries', {{}}, 'delayFrames', [], ...
    'orderRank', [], 'colors', [], 'labelXY', []);
guidata(hFig, state);

loadEventState(hFig, eventID);

out = struct('figure', hFig, 'axes', hAx, 'controls', state.controls, 'sourceLabel', sourceLabel);

    function onEventChange(src, ~)
        if ~ishandle(src)
            return
        end
        S = guidata(hFig);
        pos = min(max(round(get(src, 'Value')), 1), numel(S.eventList));
        eid = S.eventList(pos);
        loadEventState(hFig, eid);
    end

    function onFrameChange(~, ~)
        if ~ishandle(hFig)
            return
        end
        S = guidata(hFig);
        renderFrame(hFig, round(get(S.controls.frame, 'Value')));
    end

    function onPlayToggle(src, ~)
        if ~ishandle(src) || ~get(src, 'Value')
            if ishandle(src)
                set(src, 'String', 'Play');
            end
            return
        end
        set(src, 'String', 'Stop');
        while ishandle(src) && get(src, 'Value')
            S = guidata(hFig);
            v = round(get(S.controls.frame, 'Value'));
            vmax = round(get(S.controls.frame, 'Max'));
            if v >= vmax
                set(S.controls.frame, 'Value', get(S.controls.frame, 'Min'));
            else
                set(S.controls.frame, 'Value', v + 1);
            end
            renderFrame(hFig, round(get(S.controls.frame, 'Value')));
            drawnow;
            pause(0.12);
        end
        if ishandle(src)
            set(src, 'Value', 0, 'String', 'Play');
        end
    end
end

function loadEventState(hFig, eventID)
S = guidata(hFig);
eventID = round(double(eventID));
if isempty(S.eventList)
    eventID = min(max(eventID, 1), S.nEvents);
else
    if ~any(S.eventList == eventID)
        pos = min(max(eventID, 1), numel(S.eventList));
        eventID = S.eventList(pos);
    end
end
S.eventID = eventID;
[imageSz, spatialRes] = localEventImagingMetaExplorer(S.BWA_AQ2, eventID);
eventPos = find(S.eventList == eventID, 1, 'first');
if isempty(eventPos)
    eventPos = 1;
end
set(S.controls.event, 'Value', eventPos);
set(S.controls.eventText, 'String', localEventLabel(S.BWA_AQ2, S.filesTbl, eventPos, numel(S.eventList), eventID));
set(S.controls.info, 'String', sprintf('%s | loading event %d...', char(string(S.sourceLabel)), eventID));
set(S.controls.frame, 'Min', 1, 'Max', 1, 'Value', 1, 'SliderStep', [1 1]);
set(S.controls.frameText, 'String', 'loading...');
renderLoadingPreview(hFig, eventID, imageSz, spatialRes);
drawnow limitrate;

geomAll = CalcEventFrameGeometry(S.fts.loc.xSpaTemp, S.fts.loc.t0, S.fts.loc.t1, imageSz, 'eventIDs', eventID);
S.geometry = geomAll{eventID};
S.imageSz = imageSz;
S.spatialRes = spatialRes;
S.wavefront = [];
S.centroid = [];
if S.showTracking
    S.wavefront = CalcWavefrontMetrics(S.geometry, S.spatialRes);
    if strcmpi(string(S.wavefrontMode), "prototype") && exist('wavefrontComponentGraphPrototype', 'file') == 2
        proto = wavefrontComponentGraphPrototype(S.BWA_AQ2, eventID, ...
            'showFigure', false, ...
            'refinePixelWalk', 'on');
        S.wavefront = applyPrototypeWavefrontExplorer(S.wavefront, proto);
    end
    S.centroid = localAdaptiveCentroidForDisplay(S.geometry, S.wavefront);
end
[cellOverlay, S] = localBuildEventCellOverlayExplorer(S, eventID);
S.cellOverlay = cellOverlay;

frames = S.geometry.frames(:);
if isempty(frames)
    frames = 1;
end

function W = applyPrototypeWavefrontExplorer(W, proto)
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
set(S.controls.event, 'Value', eventPos);
set(S.controls.eventText, 'String', localEventLabel(S.BWA_AQ2, S.filesTbl, eventPos, numel(S.eventList), eventID));
set(S.controls.info, 'String', sprintf('%s | event %d ready | %s', ...
    char(string(S.sourceLabel)), eventID, localCellOverlayStatusExplorer(S.cellOverlay)));
set(S.controls.frame, 'Min', double(frames(1)), 'Max', double(frames(end)), 'Value', double(frames(1)));
if numel(frames) > 1
    set(S.controls.frame, 'SliderStep', [1/max(numel(frames)-1,1) min(10/max(numel(frames)-1,1),1)]);
else
    set(S.controls.frame, 'SliderStep', [1 1]);
end
guidata(hFig, S);
renderFrame(hFig, double(frames(1)));
end

function renderLoadingPreview(hFig, eventID, imageSz, spatialRes)
if ~ishandle(hFig)
    return
end
S = guidata(hFig);
ax = S.axes;
cla(ax);
[bgImg, bgRef] = buildExplorerFootprintPreview(S.fts, imageSz, spatialRes, eventID);
imshow(bgImg, bgRef, 'Parent', ax);
hold(ax, 'on');
axis(ax, 'equal');
xlabel(ax, 'x (um)');
ylabel(ax, 'y (um)');
title(ax, sprintf('Event %d | loading geometry...', eventID));
end

function renderFrame(hFig, frameNow)
S = guidata(hFig);
fts = S.fts;
opts = S.opts;
eventID = S.eventID;
frameNow = max(min(round(double(frameNow)), round(S.geometry.frames(end))), round(S.geometry.frames(1)));
set(S.controls.frame, 'Value', frameNow);
set(S.controls.frameText, 'String', sprintf('%d / %d', frameNow, round(S.geometry.frames(end))));

[bgImg, bgRef] = buildExplorerBackground(fts, S.geometry, S.imageSz, S.spatialRes, eventID, frameNow);
ax = S.axes;
cla(ax);
imshow(bgImg, bgRef, 'Parent', ax);
hold(ax, 'on');
axis(ax, 'equal');
xlabel(ax, 'x (um)');
ylabel(ax, 'y (um)');
if S.showTracking && ~isempty(S.wavefront)
    title(ax, sprintf(['Event %d | frame %d | branch origin=%d | branch peak=%d | reach=%.3f um | frontier=%.3f um | ' ...
        'tortuosity=%.3f'], eventID, frameNow, round(S.wavefront.branchOriginFrame), round(S.wavefront.branchPeakFrame), ...
        S.wavefront.wavefront_max_reach_um, S.wavefront.wavefront_frontier_length_um, ...
        S.wavefront.wavefront_tortuosity));
else
    title(ax, sprintf('Event %d | frame %d', eventID, frameNow));
end

spatialRes = S.spatialRes;
legendItems = gobjects(0);
legendLabels = {};

if isfield(S, 'cellOverlay') && isstruct(S.cellOverlay) && logical(S.cellOverlay.available)
    for i = 1:numel(S.cellOverlay.boundaries)
        b = [];
        if isfield(S.cellOverlay, 'partitionBoundaries') && numel(S.cellOverlay.partitionBoundaries) >= i
            b = S.cellOverlay.partitionBoundaries{i};
        end
        if isempty(b)
            b = S.cellOverlay.boundaries{i};
        end
        if isempty(b)
            continue
        end
        color = [0.20 0.85 0.30];
        lw = 1.25;
        if isfield(S.cellOverlay, 'colors') && size(S.cellOverlay.colors, 1) >= i
            color = S.cellOverlay.colors(i, :);
        end
        if i <= numel(S.cellOverlay.isPrimary) && S.cellOverlay.isPrimary(i)
            color = [1.00 0.85 0.10];
            lw = 1.8;
        end
        plot(ax, b(:,1) * spatialRes, b(:,2) * spatialRes, '-', ...
            'Color', color, 'LineWidth', lw, 'HandleVisibility', 'off');
        if isfield(S.cellOverlay, 'labelXY') && size(S.cellOverlay.labelXY, 1) >= i && ...
                all(isfinite(S.cellOverlay.labelXY(i,:)))
            if i <= numel(S.cellOverlay.isPrimary) && S.cellOverlay.isPrimary(i)
                lbl = 'P';
            else
                lbl = '';
                if isfield(S.cellOverlay, 'orderRank') && numel(S.cellOverlay.orderRank) >= i && ...
                        isfinite(S.cellOverlay.orderRank(i))
                    lbl = char(string(S.cellOverlay.orderRank(i)));
                elseif isfield(S.cellOverlay, 'delayFrames') && numel(S.cellOverlay.delayFrames) >= i && ...
                        isfinite(S.cellOverlay.delayFrames(i))
                    lbl = char(string(S.cellOverlay.delayFrames(i)));
                end
            end
            if ~isempty(lbl)
                text(ax, S.cellOverlay.labelXY(i,1) * spatialRes, S.cellOverlay.labelXY(i,2) * spatialRes, lbl, ...
                    'Color', color, 'FontSize', 8, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                    'Clipping', 'on', 'HandleVisibility', 'off');
            end
        end
    end
end

hLegCent = plot(ax, NaN, NaN, '-o', ...
    'Color', [0.10 0.80 0.80], 'LineWidth', 1.5, 'MarkerSize', 3, ...
    'MarkerFaceColor', [0.10 0.80 0.80], 'HandleVisibility', 'on');
hLegTrack = plot(ax, NaN, NaN, '-o', ...
    'Color', [0.85 0.20 0.20], 'LineWidth', 2, 'MarkerSize', 4, ...
    'MarkerFaceColor', [0.85 0.20 0.20], 'HandleVisibility', 'on');
hLegCurrent = plot(ax, NaN, NaN, 'rp', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 11, 'LineWidth', 1.25, 'HandleVisibility', 'on');
hLegSecondary = plot(ax, NaN, NaN, '-o', ...
    'Color', [0.55 0.55 0.55], 'LineWidth', 2, 'MarkerSize', 4, ...
    'MarkerFaceColor', [0.55 0.55 0.55], 'HandleVisibility', 'on');
hLegCell = plot(ax, NaN, NaN, '-', 'Color', [0.20 0.85 0.30], 'LineWidth', 1.25, 'HandleVisibility', 'on');
hLegPrimaryCell = plot(ax, NaN, NaN, '-', 'Color', [1.00 0.85 0.10], 'LineWidth', 1.8, 'HandleVisibility', 'on');
if S.showTracking && ~isempty(S.wavefront)
    legendItems = [hLegCent; hLegTrack; hLegCurrent; hLegSecondary; hLegCell; hLegPrimaryCell]; %#ok<AGROW>
    legendLabels = {'centroid track','branch track','current branch point','secondary branch', ...
        'participating cell','primary cell'};
else
    delete([hLegCent; hLegTrack; hLegCurrent; hLegSecondary]);
    legendItems = [hLegCell; hLegPrimaryCell]; %#ok<AGROW>
    legendLabels = {'participating cell','primary cell'};
end

if S.showTracking && ~isempty(S.wavefront)
    centKeep = S.centroid.ct >= S.wavefront.branchOriginFrame & S.centroid.ct <= frameNow;
    if any(centKeep)
        centPts = [S.wavefront.originXY; [S.centroid.cx(centKeep), S.centroid.cy(centKeep)]];
        centFrames = [S.wavefront.branchOriginFrame; S.centroid.ct(centKeep)];
        keepUnique = [true; diff(centFrames) > 0 | any(diff(centPts,1,1) ~= 0, 2)];
        centPts = centPts(keepUnique, :);
        plot(ax, centPts(:,1) * spatialRes, centPts(:,2) * spatialRes, '-o', ...
            'Color', [0.10 0.80 0.80], 'LineWidth', 1.5, 'MarkerSize', 3, ...
            'MarkerFaceColor', [0.10 0.80 0.80]);
    end
    if isfield(S.wavefront, 'secondaryTracks') && ~isempty(S.wavefront.secondaryTracks)
        for i = 1:numel(S.wavefront.secondaryTracks)
            T = S.wavefront.secondaryTracks{i};
            if ~isstruct(T) || ~isfield(T, 'trackValid') || ~isfield(T, 'trackXY')
                continue
            end
            secKeep = T.trackValid & (S.wavefront.frames >= S.wavefront.branchOriginFrame) & (S.wavefront.frames <= frameNow);
            if ~any(secKeep)
                continue
            end
            secPts = T.trackXY(secKeep, :) * spatialRes;
            plotTrackRunsExplorer(ax, S.wavefront.frames(secKeep), secPts, [0.55 0.55 0.55]);
        end
    end

    wfKeep = S.wavefront.trackValid & (S.wavefront.frames >= S.wavefront.branchOriginFrame) & (S.wavefront.frames <= frameNow);
    if any(wfKeep)
        pts = S.wavefront.trackXY(wfKeep, :) * spatialRes;
        plotTrackRunsExplorer(ax, S.wavefront.frames(wfKeep), pts, [0.85 0.20 0.20]);
        plot(ax, pts(end,1), pts(end,2), 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 11, 'LineWidth', 1.25);
    end

    hBranchOrigin = plot(ax, S.wavefront.originXY(1) * spatialRes, S.wavefront.originXY(2) * spatialRes, ...
        'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 8, 'LineWidth', 1.25);
    hBranchPeak = plot(ax, S.wavefront.peakPointXY(1) * spatialRes, S.wavefront.peakPointXY(2) * spatialRes, ...
        'r*', 'MarkerSize', 9, 'LineWidth', 1.25);
    legendItems(end+1:end+2,1) = [hBranchOrigin; hBranchPeak];
    legendLabels(end+1:end+2) = {'branch origin','branch peak'};
end

legend(ax, legendItems, legendLabels, 'Location', 'bestoutside');
end

function cent = localAdaptiveCentroidForDisplay(G, W)
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
    [compCtr, compArea, ok] = localStableCentroidComponent(mask, refPt);
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

function [ctr, area, ok] = localStableCentroidComponent(mask, refPt)
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

function [BWA_AQ2, label, eventList, filesTbl] = resolveLocalEventMotionSourceLocal(source)
label = "";
eventList = [];
filesTbl = table();
if isstruct(source) && isfield(source, 'derived') && isstruct(source.derived) && ...
        isfield(source.derived, 'context') && isstruct(source.derived.context) && ...
        isfield(source.derived.context, 'BWA_AQ2')
    cand = source.derived.context.BWA_AQ2;
    if isstruct(cand) && isfield(cand, 'fts1') && isfield(cand, 'opts')
        BWA_AQ2 = cand;
        if isfield(source, 'pathIn')
            label = string(source.pathIn);
        end
        eventList = localResultEventList(source, BWA_AQ2);
        filesTbl = localResolveFilesTable(source, table());
        return
    end
end
if isstruct(source) && isfield(source, 'pathIn') && strlength(string(source.pathIn)) > 0
    loaded = loadClusterOutputs(char(string(source.pathIn)), 'attachLegacyMeta', true, 'attachEventMeta', true);
    BWA_AQ2 = loaded.BWA_AQ2;
    label = string(source.pathIn);
    eventList = localResultEventList(source, BWA_AQ2);
    if isfield(loaded, 'files') && istable(loaded.files)
        filesTbl = localResolveFilesTable(source, loaded.files);
    else
        filesTbl = localResolveFilesTable(source, table());
    end
    return
end
if ischar(source) || isstring(source)
    loaded = loadClusterOutputs(char(string(source)), 'attachLegacyMeta', true, 'attachEventMeta', true);
    BWA_AQ2 = loaded.BWA_AQ2;
    label = string(source);
    if isfield(loaded, 'files') && istable(loaded.files)
        filesTbl = loaded.files;
    end
    return
end
if isstruct(source) && isfield(source, 'fts1') && isfield(source, 'opts')
    BWA_AQ2 = source;
    return
end
error('explore_local_event_motion:BadSource', 'Expected result, pathIn, or BWA_AQ2.');
end

function eventID = localResolveInitialEventID(eventList, requested)
if isempty(eventList)
    eventID = max(round(double(requested)), 1);
    return
end
requested = round(double(requested));
if any(eventList == requested)
    eventID = requested;
    return
end
pos = min(max(requested, 1), numel(eventList));
eventID = eventList(pos);
end

function txt = localEventLabel(BWA_AQ2, filesTbl, eventPos, nSelectable, eventID)
fileLabel = localShortFileLabel(BWA_AQ2, filesTbl, eventID);
localLabel = localFileEventLabel(BWA_AQ2, eventID);
txt = sprintf('Subset %d/%d | %s | %s', eventPos, nSelectable, fileLabel, localLabel);
end

function localLabel = localFileEventLabel(BWA_AQ2, eventID)
localLabel = 'E?';
if ~(isstruct(BWA_AQ2) && isfield(BWA_AQ2, 'fts1'))
    return
end
fts = BWA_AQ2.fts1;
if isfield(fts, 'basic') && isfield(fts.basic, 'fileID') && ~isempty(fts.basic.fileID)
    fileIDs = string(fts.basic.fileID);
elseif isfield(fts, 'fileID') && ~isempty(fts.fileID)
    fileIDs = string(fts.fileID);
else
    return
end
if eventID < 1 || eventID > numel(fileIDs)
    return
end
thisID = fileIDs(eventID);
localID = nnz(fileIDs(1:eventID) == thisID);
if isfinite(localID) && localID >= 1
    localLabel = sprintf('E%d', localID);
end
end

function fileLabel = localShortFileLabel(BWA_AQ2, filesTbl, eventID)
fileLabel = 'unknown_file';
if ~(isstruct(BWA_AQ2) && isfield(BWA_AQ2, 'fts1'))
    return
end
fts = BWA_AQ2.fts1;
if isfield(fts, 'basic') && isfield(fts.basic, 'fileID') && ~isempty(fts.basic.fileID)
    fileIDs = string(fts.basic.fileID);
elseif isfield(fts, 'fileID') && ~isempty(fts.fileID)
    fileIDs = string(fts.fileID);
else
    return
end
if eventID < 1 || eventID > numel(fileIDs)
    return
end
fileID = char(fileIDs(eventID));
if istable(filesTbl) && ~isempty(filesTbl) && ...
        all(ismember({'file_id','exp_date','slice_num'}, filesTbl.Properties.VariableNames))
    idx = find(string(filesTbl.file_id) == string(fileID), 1, 'first');
    if ~isempty(idx)
        labels = bwaMakeSliceLabels(filesTbl.file_id(idx), filesTbl(idx,:));
        cand = strtrim(char(string(labels{1})));
        if ~isempty(cand)
            fileLabel = cand;
            return
        end
    end
end
if exist('bwaFormatFileIDLabel', 'file') == 2
    cand = char(string(bwaFormatFileIDLabel(fileID)));
    cand = strtrim(cand);
    if ~isempty(cand)
        fileLabel = cand;
        return
    end
end
fileLabel = strtrim(fileID);
if isempty(fileLabel)
    fileLabel = 'unknown_file';
end
end

function filesTbl = localResolveFilesTable(result, fallbackTbl)
filesTbl = fallbackTbl;
if isstruct(result) && isfield(result, 'canonical') && isstruct(result.canonical) && ...
        isfield(result.canonical, 'subset') && isstruct(result.canonical.subset) && ...
        isfield(result.canonical.subset, 'files') && istable(result.canonical.subset.files) && ...
        ~isempty(result.canonical.subset.files)
    filesTbl = result.canonical.subset.files;
    return
end
if isstruct(result) && isfield(result, 'derived') && isstruct(result.derived) && ...
        isfield(result.derived, 'context') && isstruct(result.derived.context) && ...
        isfield(result.derived.context, 'filesTbl') && istable(result.derived.context.filesTbl) && ...
        ~isempty(result.derived.context.filesTbl)
    filesTbl = result.derived.context.filesTbl;
end
end

function eventList = localResultEventList(result, BWA_AQ2)
eventList = [];
if ~(isstruct(result) && isfield(result, 'metrics') && isstruct(result.metrics) && ...
        isfield(result.metrics, 'tables') && isstruct(result.metrics.tables))
    return
end
T = result.metrics.tables;
uids = strings(0,1);
tblNames = fieldnames(T);
for i = 1:numel(tblNames)
    tbl = T.(tblNames{i});
    if istable(tbl) && ismember('event_uid', tbl.Properties.VariableNames)
        uids = [uids; string(tbl.event_uid(:))]; %#ok<AGROW>
    end
end
uids = unique(uids(strlength(uids) > 0), 'stable');
if isempty(uids)
    return
end
eventList = NaN(numel(uids), 1);
keep = false(numel(uids), 1);
for i = 1:numel(uids)
    try
        eventList(i) = localEventIndexFromUID(BWA_AQ2, uids(i));
        keep(i) = isfinite(eventList(i));
    catch
        keep(i) = false;
    end
end
eventList = eventList(keep);
end

function sourcePath = localResolveSourcePathExplorer(source, fallbackLabel)
sourcePath = "";
if isstruct(source) && isfield(source, 'pathIn') && strlength(string(source.pathIn)) > 0
    sourcePath = string(source.pathIn);
elseif ischar(source) || isstring(source)
    sourcePath = string(source);
elseif strlength(string(fallbackLabel)) > 0
    sourcePath = string(fallbackLabel);
end
end

function T = localResolveEventPartitionTableExplorer(source)
T = table();
if isstruct(source) && isfield(source, 'metrics') && isstruct(source.metrics) && ...
        isfield(source.metrics, 'tables') && isstruct(source.metrics.tables) && ...
        isfield(source.metrics.tables, 'event_partition') && istable(source.metrics.tables.event_partition)
    T = source.metrics.tables.event_partition;
end
end

function [img, R] = buildExplorerBackground(fts, G, imageSz, spatialRes, eventID, frameNow)
sz = double(imageSz(:)');
sx = round(sz(1));
sy = round(sz(2));
sz3 = 1;
if numel(sz) >= 3 && isfinite(sz(3)) && sz(3) > 0
    sz3 = round(sz(3));
end
img = zeros(sx, sy, 3);
footMask = localGeometryFootprintMask(G, [sx sy]);
if isempty(footMask)
    footMask = false(sx, sy);
    footMask(double(fts.loc.xSpa{eventID})) = true;
end
img = applyMaskColorExplorer(img, footMask, [0.72 0.72 0.72]);

prevMask = false(sx, sy);
currMask = false(sx, sy);
if isstruct(G) && isfield(G, 'frames') && isfield(G, 'frameMasks') && ~isempty(G.frames)
    idxNow = find(double(G.frames(:)) == double(frameNow), 1, 'first');
    if ~isempty(idxNow)
        currMask = localMaskLikeImageSize(G.frameMasks{idxNow}, [sx sy]);
        idxPrev = max(idxNow - 1, 1);
        prevMask = localMaskLikeImageSize(G.frameMasks{idxPrev}, [sx sy]);
    end
else
    t0 = double(fts.loc.t0(eventID));
    framePrev = max(t0, frameNow - 1);
    prevMask = explorerFrameMask(fts.loc.xSpaTemp{eventID}, framePrev, sx, sy, sz3);
    currMask = explorerFrameMask(fts.loc.xSpaTemp{eventID}, frameNow, sx, sy, sz3);
end
img = applyMaskColorExplorer(img, prevMask & ~currMask, [0.35 0.00 0.60]);
img = applyMaskColorExplorer(img, prevMask & currMask, [0.00 0.35 0.90]);
img = applyMaskColorExplorer(img, currMask & ~prevMask, [0.85 0.20 0.80]);
R = imref2d([sx sy], spatialRes, spatialRes);
end

function [overlay, S] = localBuildEventCellOverlayExplorer(S, eventID)
overlay = struct('available', false, 'reason', "", 'cellUID', strings(0,1), ...
    'isPrimary', false(0,1), 'boundaries', {{}}, 'partitionBoundaries', {{}}, 'delayFrames', [], ...
    'orderRank', [], 'colors', [], 'labelXY', []);
T = S.eventPartition;
if ~(istable(T) && ~isempty(T) && ismember('event_uid', T.Properties.VariableNames))
    overlay.reason = "no_event_partition";
    return
end
eventUID = localEventUIDFromIndexExplorer(S.BWA_AQ2, eventID);
if strlength(eventUID) == 0
    overlay.reason = "missing_event_uid";
    return
end
keep = string(T.event_uid) == eventUID;
if ~any(keep)
    overlay.reason = "event_not_in_partition";
    return
end
cellCol = localFirstPresentExplorer(T.Properties.VariableNames, {'partition_cell_uid','cell_uid'});
if strlength(cellCol) == 0
    overlay.reason = "missing_cell_uid_column";
    return
end
cellUID = unique(string(T.(char(cellCol))(keep)), 'stable');
cellUID = cellUID(strlength(cellUID) > 0);
if isempty(cellUID)
    overlay.reason = "no_cells_for_event";
    return
end
isPrimary = false(size(cellUID));
if ismember('is_primary_cell', T.Properties.VariableNames)
    rows = T(keep, :);
    for i = 1:numel(cellUID)
        isPrimary(i) = any(string(rows.(char(cellCol))) == cellUID(i) & logical(rows.is_primary_cell));
    end
elseif ismember('primary_cell_uid', T.Properties.VariableNames)
    primaryUID = unique(string(T.primary_cell_uid(keep)));
    isPrimary = ismember(cellUID, primaryUID(primaryUID ~= ""));
end
fileID = localEventFileIDExplorer(S.BWA_AQ2, eventID);
[masks, S] = localResolveCellMasksForEventExplorer(S, eventID, fileID);
overlay.cellUID = cellUID;
overlay.isPrimary = isPrimary;
if isempty(masks)
    overlay.reason = "cell_masks_unavailable";
    return
end
boundaries = cell(numel(cellUID), 1);
partitionBoundaries = cell(numel(cellUID), 1);
labelXY = NaN(numel(cellUID), 2);
partMaskEvent = localEventFootprintMaskExplorer(S.BWA_AQ2, eventID);
for i = 1:numel(cellUID)
    cellIdx = localCellIndexFromUIDExplorer(cellUID(i), fileID, S.filesTbl);
    if ~isfinite(cellIdx) || cellIdx < 1 || cellIdx > numel(masks)
        continue
    end
    mask = logical(masks{cellIdx});
    if isempty(mask) || ~any(mask(:))
        continue
    end
    B = bwboundaries(mask, 8, 'noholes');
    if isempty(B)
        continue
    end
    [~, ib] = max(cellfun(@(x) size(x,1), B));
    b = B{ib};
    boundaries{i} = [b(:,2), b(:,1)];
    drawB = b;
    partMask = [];
    if ~isempty(partMaskEvent)
        partMask = localMaskLikeImageSize(mask, size(partMaskEvent));
        partMask = partMask & partMaskEvent;
    end
    if ~isempty(partMask) && any(partMask(:))
        Bp = bwboundaries(partMask, 8, 'noholes');
        if ~isempty(Bp)
            [~, ip] = max(cellfun(@(x) size(x,1), Bp));
            bp = Bp{ip};
            partitionBoundaries{i} = [bp(:,2), bp(:,1)];
            drawB = bp;
        end
    end
    labelXY(i,:) = [mean(drawB(:,2)), mean(drawB(:,1))];
end
rows = T(keep, :);
delayFrames = localResolveCellOverlayDelayFramesExplorer(rows, cellCol, cellUID);
orderRank = localResolveCellOverlayOrderExplorer(delayFrames, isPrimary);
colors = localResolveCellOverlayColorsExplorer(delayFrames, isPrimary);
overlay.boundaries = boundaries;
overlay.partitionBoundaries = partitionBoundaries;
overlay.labelXY = labelXY;
overlay.delayFrames = delayFrames;
overlay.orderRank = orderRank;
overlay.colors = colors;
overlay.available = any(~cellfun(@isempty, boundaries));
if ~overlay.available
    overlay.reason = "cell_masks_missing_for_event_cells";
end
end

function delayFrames = localResolveCellOverlayDelayFramesExplorer(rows, cellCol, cellUID)
delayFrames = NaN(numel(cellUID), 1);
if isempty(rows) || ~istable(rows)
    return
end
prefer = "";
if ismember('cell_entry_delay_from_primary_cov_thresh_frames', rows.Properties.VariableNames)
    prefer = "cell_entry_delay_from_primary_cov_thresh_frames";
elseif ismember('cell_entry_delay_from_primary_frames', rows.Properties.VariableNames)
    prefer = "cell_entry_delay_from_primary_frames";
end
if strlength(prefer) == 0
    return
end
rowCellUID = string(rows.(char(cellCol)));
vals = double(rows.(char(prefer)));
for i = 1:numel(cellUID)
    hit = find(rowCellUID == cellUID(i), 1, 'first');
    if ~isempty(hit)
        delayFrames(i) = vals(hit);
    end
end
end

function orderRank = localResolveCellOverlayOrderExplorer(delayFrames, isPrimary)
orderRank = NaN(numel(delayFrames), 1);
nonPrimary = find(~isPrimary(:));
if isempty(nonPrimary)
    return
end
vals = delayFrames(nonPrimary);
good = isfinite(vals);
if ~any(good)
    return
end
[~, ord] = sort(vals(good), 'ascend');
tmpIdx = nonPrimary(good);
tmpIdx = tmpIdx(ord);
orderRank(tmpIdx) = 1:numel(tmpIdx);
end

function colors = localResolveCellOverlayColorsExplorer(delayFrames, isPrimary)
colors = repmat([0.20 0.85 0.30], numel(delayFrames), 1);
if isempty(delayFrames)
    return
end
nonPrimary = find(~isPrimary(:));
vals = delayFrames(nonPrimary);
good = isfinite(vals);
if ~any(good)
    return
end
u = unique(vals(good));
if numel(u) == 1
    t = ones(nnz(good), 1) * 0.5;
else
    t = (vals(good) - min(u)) ./ max(max(u) - min(u), eps);
end
early = [0.10 0.70 0.95];
late = [0.95 0.35 0.20];
colors(nonPrimary(good), :) = early + (late - early) .* t;
colors(isPrimary(:), :) = repmat([1.00 0.85 0.10], nnz(isPrimary), 1);
end

function [masks, S] = localResolveCellMasksForEventExplorer(S, eventID, fileID)
masks = {};
if strlength(fileID) == 0
    return
end
cache = S.cellMaskCache;
if isfield(cache, 'fileIDs') && isfield(cache, 'masks') && ~isempty(cache.fileIDs)
    idx = find(string(cache.fileIDs) == string(fileID), 1, 'first');
    if ~isempty(idx)
        masks = cache.masks{idx};
        return
    end
end
masks = localCurrentLegacyMasksIfMatchingFile(S.BWA_AQ2, S.filesTbl, eventID, fileID);
if isempty(masks)
    masks = localLoadCellMasksFromLocalAQuAFile(S.sourcePath, S.filesTbl, fileID);
end
if ~isempty(masks)
    cache.fileIDs(end+1,1) = string(fileID);
    cache.masks{end+1,1} = masks;
    S.cellMaskCache = cache;
end
end

function masks = localCurrentLegacyMasksIfMatchingFile(BWA_AQ2, filesTbl, eventID, fileID)
masks = {};
if ~(isstruct(BWA_AQ2) && isfield(BWA_AQ2, 'fts1') && isfield(BWA_AQ2.fts1, 'region') && ...
        isfield(BWA_AQ2.fts1.region, 'cell') && isfield(BWA_AQ2.fts1.region.cell, 'mask'))
    return
end
allMasks = BWA_AQ2.fts1.region.cell.mask;
if ~iscell(allMasks) || isempty(allMasks)
    return
end
fileIDs = localEventFileIDsExplorer(BWA_AQ2);
if isempty(fileIDs) || eventID < 1 || eventID > numel(fileIDs)
    return
end
if istable(filesTbl) && ~isempty(filesTbl) && all(ismember({'file_id','n_cells'}, filesTbl.Properties.VariableNames))
    totalCells = sum(double(filesTbl.n_cells), 'omitnan');
    if isfinite(totalCells) && totalCells == numel(allMasks)
        masks = allMasks;
        return
    end
end
if string(fileID) == string(fileIDs(1))
    masks = allMasks;
end
end

function masks = localLoadCellMasksFromLocalAQuAFile(sourcePath, filesTbl, fileID)
masks = {};
if strlength(string(sourcePath)) == 0 || ~(istable(filesTbl) && ~isempty(filesTbl)) || ...
        ~all(ismember({'file_id','filename'}, filesTbl.Properties.VariableNames))
    return
end
idx = find(string(filesTbl.file_id) == string(fileID), 1, 'first');
if isempty(idx)
    return
end
fname = char(string(filesTbl.filename(idx)));
if isempty(fname)
    return
end
d = dir(fullfile(char(string(sourcePath)), '**', fname));
if isempty(d)
    return
end
try
    s = load(fullfile(d(1).folder, d(1).name));
    masks = localExtractCellMasksFromStructExplorer(s);
catch
    masks = {};
end
end

function masks = localExtractCellMasksFromStructExplorer(s)
masks = {};
candidates = {
    {'BWA_AQ2','fts1','region','cell','mask'}
    {'res','fts1','region','cell','mask'}
    {'datOut','fts1','region','cell','mask'}
    {'fts1','region','cell','mask'}
    };
for i = 1:numel(candidates)
    try
        v = localGetFieldPathExplorer(s, candidates{i});
        if iscell(v)
            masks = v;
            return
        end
    catch
    end
end
end

function v = localGetFieldPathExplorer(s, path)
v = s;
for i = 1:numel(path)
    v = v.(path{i});
end
end

function fileIDs = localEventFileIDsExplorer(BWA_AQ2)
fileIDs = strings(0,1);
if ~(isstruct(BWA_AQ2) && isfield(BWA_AQ2, 'fts1'))
    return
end
fts = BWA_AQ2.fts1;
if isfield(fts, 'basic') && isfield(fts.basic, 'fileID') && ~isempty(fts.basic.fileID)
    fileIDs = string(fts.basic.fileID);
elseif isfield(fts, 'fileID') && ~isempty(fts.fileID)
    fileIDs = string(fts.fileID);
end
end

function fileID = localEventFileIDExplorer(BWA_AQ2, eventID)
fileID = "";
fileIDs = localEventFileIDsExplorer(BWA_AQ2);
if isempty(fileIDs) || eventID < 1 || eventID > numel(fileIDs)
    return
end
fileID = string(fileIDs(eventID));
end

function eventUID = localEventUIDFromIndexExplorer(BWA_AQ2, eventID)
eventUID = "";
fileIDs = localEventFileIDsExplorer(BWA_AQ2);
if isempty(fileIDs) || eventID < 1 || eventID > numel(fileIDs)
    return
end
thisID = fileIDs(eventID);
localID = nnz(fileIDs(1:eventID) == thisID);
eventUID = string(sprintf('%s_E%d', char(thisID), localID));
end

function name = localFirstPresentExplorer(names, candidates)
name = "";
for i = 1:numel(candidates)
    if any(strcmp(names, candidates{i}))
        name = string(candidates{i});
        return
    end
end
end

function idx = localCellIndexFromUIDExplorer(cellUID, fileID, filesTbl)
idx = NaN;
uid = char(string(cellUID));
fid = char(string(fileID));
prefix = [fid '_C'];
if ~startsWith(uid, prefix)
    return
end
suffix = extractAfter(string(uid), strlength(string(prefix)));
localIdx = double(str2double(suffix));
if ~isfinite(localIdx)
    idx = NaN;
    return
end
idx = localIdx;
if istable(filesTbl) && ~isempty(filesTbl) && all(ismember({'file_id','n_cells'}, filesTbl.Properties.VariableNames))
    fileIDs = string(filesTbl.file_id);
    row = find(fileIDs == string(fileID), 1, 'first');
    if isempty(row)
        idx = NaN;
        return
    end
    priorCounts = double(filesTbl.n_cells(1:(row - 1)));
    priorCounts = priorCounts(isfinite(priorCounts));
    idx = sum(priorCounts) + localIdx;
end
end

function txt = localCellOverlayStatusExplorer(overlay)
txt = 'cell borders: unavailable';
if ~(isstruct(overlay) && isfield(overlay, 'available'))
    return
end
if logical(overlay.available)
    n = 0;
    if isfield(overlay, 'boundaries') && ~isempty(overlay.boundaries)
        n = nnz(~cellfun(@isempty, overlay.boundaries));
    elseif isfield(overlay, 'cellUID')
        n = numel(overlay.cellUID);
    end
    txt = sprintf('cell borders: %d loaded', n);
    return
end
if isfield(overlay, 'reason') && strlength(string(overlay.reason)) > 0
    txt = sprintf('cell borders: %s', char(string(overlay.reason)));
end
end

function footMask = localGeometryFootprintMask(G, imageSize2D)
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

function footMask = localEventFootprintMaskExplorer(BWA_AQ2, eventID)
footMask = [];
try
    [imageSz, ~] = localEventImagingMetaExplorer(BWA_AQ2, eventID);
catch
    return
end
sz = double(imageSz(:)');
if numel(sz) < 2 || any(~isfinite(sz(1:2)))
    return
end
sx = round(sz(1));
sy = round(sz(2));
sz3 = 1;
if numel(sz) >= 3 && isfinite(sz(3)) && sz(3) > 0
    sz3 = round(sz(3));
end
fts = BWA_AQ2.fts1;
footMask = false(sx, sy);
if isfield(fts, 'loc') && isfield(fts.loc, 'xSpaTemp') && numel(fts.loc.xSpaTemp) >= eventID && ...
        ~isempty(fts.loc.xSpaTemp{eventID}) && isfield(fts.loc, 't0') && isfield(fts.loc, 't1')
    t0 = round(double(fts.loc.t0(eventID)));
    t1 = round(double(fts.loc.t1(eventID)));
    for frameNow = t0:t1
        footMask = footMask | explorerFrameMask(fts.loc.xSpaTemp{eventID}, frameNow, sx, sy, sz3);
    end
elseif isfield(fts.loc, 'xSpa') && numel(fts.loc.xSpa) >= eventID && ~isempty(fts.loc.xSpa{eventID})
    footMask(double(fts.loc.xSpa{eventID})) = true;
end
end

function maskOut = localMaskLikeImageSize(maskIn, imageSize2D)
maskOut = false(imageSize2D(1), imageSize2D(2));
if isempty(maskIn)
    return
end
maskIn = logical(maskIn);
sx = min(size(maskIn, 1), imageSize2D(1));
sy = min(size(maskIn, 2), imageSize2D(2));
maskOut(1:sx, 1:sy) = maskIn(1:sx, 1:sy);
end

function [img, R] = buildExplorerFootprintPreview(fts, imageSz, spatialRes, eventID)
sz = double(imageSz(:)');
sx = round(sz(1));
sy = round(sz(2));
sz3 = 1;
if numel(sz) >= 3 && isfinite(sz(3)) && sz(3) > 0
    sz3 = round(sz(3));
end
img = zeros(sx, sy, 3);
footMask = false(sx, sy);
if isfield(fts, 'loc') && isfield(fts.loc, 'xSpaTemp') && numel(fts.loc.xSpaTemp) >= eventID && ...
        ~isempty(fts.loc.xSpaTemp{eventID}) && isfield(fts.loc, 't0') && isfield(fts.loc, 't1')
    t0 = round(double(fts.loc.t0(eventID)));
    t1 = round(double(fts.loc.t1(eventID)));
    for frameNow = t0:t1
        footMask = footMask | explorerFrameMask(fts.loc.xSpaTemp{eventID}, frameNow, sx, sy, sz3);
    end
else
    footMask(double(fts.loc.xSpa{eventID})) = true;
end
img = applyMaskColorExplorer(img, footMask, [0.72 0.72 0.72]);
R = imref2d([sx sy], spatialRes, spatialRes);
end

function mask = explorerFrameMask(xSpaTemp_ID, frameNum, sx, sy, sz3)
limitStart = sx * sy * sz3 * (frameNum - 1) + 1;
limitEnd = limitStart + sx * sy * sz3;
linMovie = double(xSpaTemp_ID(:));
linFrame = linMovie(linMovie >= limitStart & linMovie <= limitEnd) - limitStart + 1;
linFrame = linFrame(isfinite(linFrame) & linFrame >= 1 & linFrame <= sx * sy);
mask = false(sx, sy);
mask(linFrame) = true;
end

function img = applyMaskColorExplorer(img, mask, color)
for c = 1:3
    plane = img(:,:,c);
    plane(mask) = color(c);
    img(:,:,c) = plane;
end
end

function [imageSz, spatialRes] = localEventImagingMetaExplorer(BWA_AQ2, eventID)
opts = BWA_AQ2.opts;
imageSz = double(opts.sz(1,:));
spatialRes = double(opts.spatialRes(1));
if ~(isstruct(BWA_AQ2) && isfield(BWA_AQ2, 'fts1') && isfield(BWA_AQ2.fts1, 'fileID'))
    return
end
fileIDs = canonCellstrLocalExplorer(BWA_AQ2.fts1.fileID);
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

function c = canonCellstrLocalExplorer(v)
if iscell(v)
    c = cell(size(v));
    for i = 1:numel(v)
        c{i} = char(string(v{i}));
    end
elseif isstring(v)
    c = cellstr(v(:));
elseif ischar(v)
    c = {v};
else
    c = cellstr(string(v(:)));
end
end

function h = plotTrackRunsExplorer(ax, frames, pts, color)
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
        'Color', color, 'MarkerFaceColor', color, 'LineWidth', 2, 'MarkerSize', 4);
    runStart = i;
end
end

function eventID = localEventIndexFromUID(BWA_AQ2, eventUID)
eventID = NaN;
fts = BWA_AQ2.fts1;
if isfield(fts, 'basic') && isfield(fts.basic, 'fileID') && ~isempty(fts.basic.fileID)
    fileIDs = string(fts.basic.fileID);
elseif isfield(fts, 'fileID') && ~isempty(fts.fileID)
    fileIDs = string(fts.fileID);
else
    error('explore_local_event_motion:MissingFileIDs', 'Could not resolve eventUID without file IDs.');
end
tok = regexp(char(eventUID), '^(.*)_E(\d+)$', 'tokens', 'once');
if isempty(tok)
    error('explore_local_event_motion:BadEventUID', 'Expected eventUID like <file_id>_E<number>.');
end
fileID = string(tok{1});
localID = str2double(tok{2});
if ~isfinite(localID) || localID < 1
    error('explore_local_event_motion:BadEventUID', 'Invalid eventUID suffix in %s.', char(eventUID));
end
idx = find(fileIDs == fileID);
if numel(idx) < localID
    error('explore_local_event_motion:EventUIDNotFound', 'eventUID %s was not found in loaded BWA_AQ2.', char(eventUID));
end
eventID = idx(localID);
end
