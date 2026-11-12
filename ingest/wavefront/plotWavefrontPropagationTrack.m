function out = plotWavefrontPropagationTrack(source, varargin)
BWA_setupPaths(fileparts(mfilename('fullpath')));

%PLOTWAVEFRONTPROPAGATIONTRACK Visualize backward-tracked wavefront frontier.
%   OUT = PLOTWAVEFRONTPROPAGATIONTRACK(PROPAGATION, OPTS, EVENTID)
%   visualizes one event using explicit propagation/opts inputs.
%
%   OUT = PLOTWAVEFRONTPROPAGATIONTRACK(RES, EVENTID)
%   accepts an AQuA2-style res struct directly and extracts:
%     - propagation = res.fts1.propagation
%     - fts         = res.fts1
%     - opts        = res.opts
%
%   Name-Value:
%     'spatialRes'      : um/px override (default opts.spatialRes(1) or 1)
%     'fts'             : optional AQuA2 fts/fts1 struct for event-footprint
%                         and per-frame background overlay
%     'backgroundMode'  : 'none', 'event_footprint', or 'event_progression'
%                         (default auto)
%     'progressFrame'   : frame to highlight for event_progression; default
%                         is the peak wavefront frame
%     'showFrameLabels' : label valid frames (default true)
%     'showCenters'     : draw reconstructed centers (default true)
%     'figureName'      : figure title/name
%     plus any wavefrontPropagationMetrics name-value args

[propagation, opts, eventID, defaultFts, extraArgs] = resolveWavefrontInputs(source, varargin{:});

p = inputParser;
p.StructExpand = false;
addParameter(p, 'spatialRes', inferSpatialRes(opts), @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'fts', defaultFts, @isstruct);
addParameter(p, 'backgroundMode', 'auto', @(x) any(strcmpi(string(x), ["auto","none","event_footprint","event_progression"])));
addParameter(p, 'progressFrame', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'showFrameLabels', true, @islogical);
addParameter(p, 'showCenters', true, @islogical);
addParameter(p, 'figureName', sprintf('Wavefront_Event%d', eventID), @(x) ischar(x) || isstring(x));
parse(p, extraArgs{:});

if ~isstruct(propagation) || ~isfield(propagation, 'propGrow') || ~iscell(propagation.propGrow) || ...
        ~isfield(propagation, 'propShrink') || ~iscell(propagation.propShrink)
    error('plotWavefrontPropagationTrack:BadPropagation', ...
        'propagation must contain cell arrays propGrow and propShrink.');
end
if eventID < 1 || eventID > numel(propagation.propGrow)
    error('plotWavefrontPropagationTrack:BadEventID', 'eventID is out of range.');
end

metricArgs = forwardMetricArgs(extraArgs);
[metrics, debug] = wavefrontPropagationMetrics( ...
    propagation.propGrow{eventID}, propagation.propShrink{eventID}, ...
    'spatialRes', p.Results.spatialRes, metricArgs{:});

hFig = figure('Name', char(string(p.Results.figureName)));
ax = axes('Parent', hFig); hold(ax, 'on');
axis(ax, 'equal');
grid(ax, 'on');
xlabel(ax, 'x (um)');
ylabel(ax, 'y (um)');
title(ax, sprintf(['Event %d | reach=%.3f um | frontier=%.3f um | ' ...
    'tortuosity=%.3f | peak frame=%d'], ...
    eventID, metrics.wavefront_max_reach_um, metrics.wavefront_frontier_length_um, ...
    metrics.wavefront_tortuosity, metrics.wavefront_peak_frame));

backgroundMode = resolveBackgroundMode(p.Results.backgroundMode, p.Results.fts);
bgHandle = gobjects(0);
if strcmp(backgroundMode, 'event_footprint')
    [bgImg, bgRef] = buildEventFootprintBackground(p.Results.fts, opts, eventID, p.Results.spatialRes);
    if ~isempty(bgImg)
        bgHandle = imshow(bgImg, bgRef, 'Parent', ax);
        uistack(bgHandle, 'bottom');
    end
elseif strcmp(backgroundMode, 'event_progression')
    [bgImg, bgRef] = buildEventProgressionBackground( ...
        p.Results.fts, opts, eventID, p.Results.spatialRes, metrics, p.Results.progressFrame);
    if ~isempty(bgImg)
        bgHandle = imshow(bgImg, bgRef, 'Parent', ax);
        uistack(bgHandle, 'bottom');
    end
end

offsetPx = [0 0];
if ~strcmp(backgroundMode, 'none')
    offsetPx = inferEventSpatialOffsetPx(p.Results.fts, opts, eventID);
end

validIdx = find(debug.frameValid);
nValid = numel(validIdx);
frameColors = parula(max(nValid, 2));
polyHandles = gobjects(0);
centerHandles = gobjects(0);
textHandles = gobjects(0);

for k = 1:nValid
    t = validIdx(k);
    pts = squeeze(debug.vertices_px(t, :, :)).';
    if isempty(pts) || any(~isfinite(pts(:)))
        continue
    end
    pts = (pts + offsetPx) * p.Results.spatialRes;
    c = frameColors(min(k, size(frameColors,1)), :);
    polyHandles(end+1,1) = plot(ax, [pts(:,1); pts(1,1)], [pts(:,2); pts(1,2)], '-', ... %#ok<AGROW>
        'Color', c, 'LineWidth', 1.25);
    ctr = (debug.centers_px(t,:) + offsetPx) * p.Results.spatialRes;
    if p.Results.showCenters && all(isfinite(ctr))
        centerHandles(end+1,1) = plot(ax, ctr(1), ctr(2), 'o', ... %#ok<AGROW>
            'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'MarkerSize', 4);
    end
    if p.Results.showFrameLabels && all(isfinite(ctr))
        textHandles(end+1,1) = text(ax, ctr(1), ctr(2), ... %#ok<AGROW>
            sprintf(' %d', t), 'Color', c, 'FontSize', 8, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
end

originUm = (debug.origin_center_px + offsetPx) * p.Results.spatialRes;
peakUm = (debug.peak_frontier_px + offsetPx) * p.Results.spatialRes;
originHandle = plot(ax, originUm(1), originUm(2), 'ks', ...
    'MarkerFaceColor', 'y', 'MarkerSize', 8, 'LineWidth', 1.25);
peakHandle = plot(ax, peakUm(1), peakUm(2), 'rp', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 12, 'LineWidth', 1.5);

trackIdx = find(debug.frontier_track_valid);
trackHandle = gobjects(0);
if ~isempty(trackIdx)
    trackPts = (debug.frontier_track_px(trackIdx, :) + offsetPx) * p.Results.spatialRes;
    trackHandle = plot(ax, trackPts(:,1), trackPts(:,2), '-o', ...
        'Color', [0.85 0.2 0.2], 'MarkerFaceColor', [0.85 0.2 0.2], ...
        'LineWidth', 2, 'MarkerSize', 4);
end

legendItems = [originHandle; peakHandle];
legendLabels = {'origin center', 'peak frontier'};
if ~isempty(trackHandle)
    legendItems(end+1,1) = trackHandle(1);
    legendLabels{end+1} = 'backward frontier track';
end
legend(ax, legendItems, legendLabels, 'Location', 'bestoutside');

out = struct();
out.figure = hFig;
out.axes = ax;
out.metrics = metrics;
out.debug = debug;
out.handles = struct( ...
    'background', bgHandle, ...
    'polygons', polyHandles, ...
    'centers', centerHandles, ...
    'frameLabels', textHandles, ...
    'origin', originHandle, ...
    'peak', peakHandle, ...
    'track', trackHandle);
end

function [propagation, opts, eventID, defaultFts, extraArgs] = resolveWavefrontInputs(source, varargin)
defaultFts = struct();
if isAquaResStruct(source)
    res = source;
    if isempty(varargin)
        error('plotWavefrontPropagationTrack:MissingEventID', 'Event ID is required when passing res.');
    end
    eventID = double(varargin{1});
    propagation = res.fts1.propagation;
    defaultFts = res.fts1;
    opts = res.opts;
    extraArgs = varargin(2:end);
    return
end
if isstruct(source) && isfield(source, 'res') && isAquaResStruct(source.res)
    res = source.res;
    if isempty(varargin)
        error('plotWavefrontPropagationTrack:MissingEventID', 'Event ID is required when passing a loaded AQuA2 struct.');
    end
    eventID = double(varargin{1});
    propagation = res.fts1.propagation;
    defaultFts = res.fts1;
    opts = res.opts;
    extraArgs = varargin(2:end);
    return
end
if numel(varargin) < 2
    error('plotWavefrontPropagationTrack:MissingInputs', ...
        'Expected either (res, eventID) or (propagation, opts, eventID).');
end
arg2 = varargin{1};
arg3 = varargin{2};
propagation = source;
opts = arg2;
eventID = double(arg3);
extraArgs = varargin(3:end);
end

function tf = isAquaResStruct(S)
tf = isstruct(S) && isfield(S, 'fts1') && isstruct(S.fts1) && ...
    isfield(S.fts1, 'propagation') && isfield(S, 'opts');
end

function sr = inferSpatialRes(opts)
sr = 1;
try
    if isstruct(opts) && isfield(opts, 'spatialRes') && ~isempty(opts.spatialRes)
        v = double(opts.spatialRes(1));
        if isfinite(v) && v > 0
            sr = v;
        end
    end
catch
end
end

function args = forwardMetricArgs(inArgs)
args = {};
skip = ["spatialres","fts","backgroundmode","progressframe","showframelabels","showcenters","figurename"];
i = 1;
while i <= numel(inArgs)
    name = string(inArgs{i});
    if i == numel(inArgs)
        break
    end
    if any(lower(name) == skip)
        i = i + 2;
        continue
    end
    args(end+1:end+2) = inArgs(i:i+1); %#ok<AGROW>
    i = i + 2;
end
end

function mode = resolveBackgroundMode(modeIn, fts)
mode = lower(char(string(modeIn)));
if strcmp(mode, 'auto')
    if hasFootprintData(fts)
        mode = 'event_progression';
    else
        mode = 'none';
    end
end
end

function tf = hasFootprintData(fts)
tf = isstruct(fts) && isfield(fts, 'loc') && isstruct(fts.loc) && ...
    isfield(fts.loc, 'xSpa') && iscell(fts.loc.xSpa) && ...
    isfield(fts.loc, 'xSpaTemp') && iscell(fts.loc.xSpaTemp);
end

function [img, R] = buildEventFootprintBackground(fts, opts, eventID, spatialRes)
img = [];
R = [];
[sx, sy] = backgroundSize(opts);
if sx <= 0 || sy <= 0 || ~hasFootprintData(fts)
    return
end
linPix = double(fts.loc.xSpa{eventID}(:));
linPix = linPix(isfinite(linPix) & linPix >= 1 & linPix <= sx*sy);
img = zeros(sx, sy, 3);
if ~isempty(linPix)
    mask = false(sx, sy);
    mask(linPix) = true;
    for c = 1:3
        plane = img(:,:,c);
        plane(mask) = 0.75;
        img(:,:,c) = plane;
    end
end
R = imref2d([sx sy], spatialRes, spatialRes);
end

function [img, R] = buildEventProgressionBackground(fts, opts, eventID, spatialRes, metrics, progressFrame)
img = [];
R = [];
[sx, sy] = backgroundSize(opts);
if sx <= 0 || sy <= 0 || ~hasFootprintData(fts)
    return
end
img = zeros(sx, sy, 3);

fullMask = maskForEventFootprint(fts.loc.xSpa{eventID}, sx, sy);
if any(fullMask(:))
    img = applyMaskColor(img, fullMask, [0.72 0.72 0.72]);
end

frameNow = double(progressFrame);
if ~isfinite(frameNow)
    frameNow = metrics.wavefront_peak_frame;
end
if ~(isfinite(frameNow) && frameNow >= 1)
    frameNow = NaN;
end

if isfinite(frameNow)
    t0 = inferEventStartFrame(fts, eventID);
    framePrev = max(t0, frameNow - 1);
    sz3 = backgroundDepth(opts);
    prevMask = temporalFrameMask(fts.loc.xSpaTemp{eventID}, framePrev, [sx sy sz3]);
    currMask = temporalFrameMask(fts.loc.xSpaTemp{eventID}, frameNow, [sx sy sz3]);
    sharedMask = prevMask & currMask;
    prevOnlyMask = prevMask & ~currMask;
    currOnlyMask = currMask & ~prevMask;
    img = applyMaskColor(img, prevOnlyMask, [0.35 0.00 0.60]);
    img = applyMaskColor(img, sharedMask, [0.00 0.35 0.90]);
    img = applyMaskColor(img, currOnlyMask, [0.85 0.20 0.80]);
end

R = imref2d([sx sy], spatialRes, spatialRes);
end

function img = applyMaskColor(img, mask, color)
for c = 1:3
    plane = img(:,:,c);
    plane(mask) = color(c);
    img(:,:,c) = plane;
end
end

function mask = maskForEventFootprint(linPix, sx, sy)
mask = false(sx, sy);
linPix = double(linPix(:));
linPix = linPix(isfinite(linPix) & linPix >= 1 & linPix <= sx*sy);
mask(linPix) = true;
end

function frameMask = temporalFrameMask(xSpaTemp_ID, frameNum, imageSz)
sx = imageSz(1);
sy = imageSz(2);
if numel(imageSz) >= 3 && isfinite(imageSz(3)) && imageSz(3) > 0
    sz3 = imageSz(3);
else
    sz3 = 1;
end
limitStart = sx * sy * sz3 * (frameNum - 1) + 1;
limitEnd = limitStart + sx * sy * sz3;
linMovie = double(xSpaTemp_ID(:));
linFrame = linMovie(linMovie >= limitStart & linMovie <= limitEnd) - limitStart + 1;
linFrame = linFrame(isfinite(linFrame) & linFrame >= 1 & linFrame <= sx*sy);
frameMask = false(sx, sy);
frameMask(linFrame) = true;
end

function [sx, sy] = backgroundSize(opts)
sx = 0;
sy = 0;
try
    sz = double(opts.sz(1,:));
    if numel(sz) >= 2 && all(isfinite(sz(1:2))) && all(sz(1:2) > 0)
        sx = round(sz(1));
        sy = round(sz(2));
    end
catch
end
end

function sz3 = backgroundDepth(opts)
sz3 = 1;
try
    sz = double(opts.sz(1,:));
    if numel(sz) >= 3 && isfinite(sz(3)) && sz(3) > 0
        sz3 = round(sz(3));
    end
catch
end
end

function t0 = inferEventStartFrame(fts, eventID)
t0 = 1;
try
    if isfield(fts.loc, 't0') && numel(fts.loc.t0) >= eventID
        v = double(fts.loc.t0(eventID));
        if isfinite(v) && v >= 1
            t0 = round(v);
        end
    end
catch
end
end

function offsetPx = inferEventSpatialOffsetPx(fts, opts, eventID)
offsetPx = [0 0];
try
    if ~hasFootprintData(fts) || ~isstruct(opts) || ~isfield(opts, 'sz') || isempty(opts.sz)
        return
    end
    sz = double(opts.sz(1,:));
    if numel(sz) < 3 || any(~isfinite(sz(1:3))) || any(sz(1:3) <= 0)
        return
    end
    frame0 = inferEventStartFrame(fts, eventID);
    linFrame = temporalFrameLinearIdx(fts.loc.xSpaTemp{eventID}, frame0, round(sz(1)), round(sz(2)), round(sz(3)));
    if isempty(linFrame)
        linFrame = double(fts.loc.xSpa{eventID}(:));
        linFrame = linFrame(isfinite(linFrame) & linFrame >= 1 & linFrame <= round(sz(1))*round(sz(2)));
    end
    if isempty(linFrame)
        return
    end
    [row, col] = ind2sub([round(sz(1)) round(sz(2))], linFrame);
    offsetPx = [mean(col, 'omitnan') mean(row, 'omitnan')];
    if any(~isfinite(offsetPx))
        offsetPx = [0 0];
    end
catch
    offsetPx = [0 0];
end
end

function linFrame = temporalFrameLinearIdx(xSpaTemp_ID, frameNum, sx, sy, sz3)
limitStart = sx * sy * sz3 * (frameNum - 1) + 1;
limitEnd = limitStart + sx * sy * sz3;
linMovie = double(xSpaTemp_ID(:));
linFrame = linMovie(linMovie >= limitStart & linMovie <= limitEnd) - limitStart + 1;
linFrame = linFrame(isfinite(linFrame) & linFrame >= 1 & linFrame <= sx * sy);
end
