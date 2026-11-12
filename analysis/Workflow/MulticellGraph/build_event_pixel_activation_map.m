function out = build_event_pixel_activation_map(source, varargin)
%BUILD_EVENT_PIXEL_ACTIVATION_MAP Diagnostic first-activation map for one event.
%   OUT = BUILD_EVENT_PIXEL_ACTIVATION_MAP(RESULT, 'eventUID', '<file>_E###')
%   computes the first frame each pixel becomes active within the event
%   footprint and optionally plots the delay map in seconds or frames.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'source');
addParameter(p, 'eventUID', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'eventID', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'showFigure', true, @islogical);
addParameter(p, 'visible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
addParameter(p, 'delayUnits', 'sec', @(x) any(strcmpi(string(x), ["sec","frames"])));
parse(p, source, varargin{:});

[result, BWA_AQ2, sourcePath] = localResolvePixelSource(source);
eventUID = string(p.Results.eventUID);
if strlength(eventUID) == 0
    if isfinite(p.Results.eventID)
        eventUID = localEventUIDFromResultIndexPixel(result, round(double(p.Results.eventID)));
    else
        error('build_event_pixel_activation_map:MissingEvent', 'Provide eventUID or eventID.');
    end
end

eventID = localEventIndexFromUIDPixel(BWA_AQ2, eventUID);
fileID = localEventFileIDPixel(BWA_AQ2, eventID);
[imageSz, spatialRes, frameRateS] = localEventImagingMetaPixel(BWA_AQ2, eventID);
[firstFrameMap, t0, t1] = localFirstActivationFrameMapPixel(BWA_AQ2, eventID, imageSz);
zeroFrame = localBranchOriginFramePixel(result, BWA_AQ2, eventUID, eventID);
postOriginFrameMap = localPostOriginActivationFrameMapPixel(BWA_AQ2, eventID, imageSz, zeroFrame);
footMask = isfinite(postOriginFrameMap);
delayFrames = postOriginFrameMap - zeroFrame;
delaySec = delayFrames * frameRateS;

out = struct();
out.event_uid = eventUID;
out.event_id = eventID;
out.file_id = fileID;
out.path_in = sourcePath;
out.image_size = imageSz;
out.spatial_res = spatialRes;
out.frame_rate_s = frameRateS;
out.t0_frame = t0;
out.t1_frame = t1;
out.zero_frame = zeroFrame;
out.first_activation_frame = firstFrameMap;
out.post_origin_activation_frame = postOriginFrameMap;
out.delay_frames = delayFrames;
out.delay_sec = delaySec;
out.footprint_mask = footMask;

if p.Results.showFigure
    out.figure = localPlotPixelActivationMap(out, char(string(p.Results.visible)), char(string(p.Results.delayUnits)));
end
end

function [result, BWA_AQ2, sourcePath] = localResolvePixelSource(source)
result = struct();
sourcePath = "";
if isstruct(source) && isfield(source, 'metrics')
    result = source;
    if isfield(source, 'pathIn')
        sourcePath = string(source.pathIn);
    elseif isfield(source, 'source') && isstruct(source.source) && isfield(source.source, 'pathIn')
        sourcePath = string(source.source.pathIn);
    end
else
    error('build_event_pixel_activation_map:BadSource', 'Expected a local-analysis result struct.');
end
if strlength(sourcePath) == 0
    error('build_event_pixel_activation_map:MissingPath', 'Result must carry pathIn so BWA_AQ2 can be loaded.');
end
loaded = loadClusterOutputs(char(sourcePath), 'attachLegacyMeta', true, 'attachEventMeta', true);
BWA_AQ2 = loaded.BWA_AQ2;
end

function eventUID = localEventUIDFromResultIndexPixel(result, eventID)
eventUID = "";
if ~(isfield(result, 'metrics') && isstruct(result.metrics) && isfield(result.metrics, 'tables') && ...
        isstruct(result.metrics.tables) && isfield(result.metrics.tables, 'event_canonical') && ...
        istable(result.metrics.tables.event_canonical))
    return
end
T = result.metrics.tables.event_canonical;
if eventID >= 1 && eventID <= height(T) && ismember('event_uid', T.Properties.VariableNames)
    eventUID = string(T.event_uid(eventID));
end
end

function eventID = localEventIndexFromUIDPixel(BWA_AQ2, eventUID)
eventID = NaN;
fileIDs = localEventFileIDsPixel(BWA_AQ2);
tok = regexp(char(eventUID), '^(.*)_E(\d+)$', 'tokens', 'once');
if isempty(tok)
    error('build_event_pixel_activation_map:BadEventUID', 'Expected eventUID like <file_id>_E<number>.');
end
fileID = string(tok{1});
localID = str2double(tok{2});
idx = find(fileIDs == fileID);
if numel(idx) < localID
    error('build_event_pixel_activation_map:EventUIDNotFound', 'eventUID %s not found in BWA_AQ2.', char(eventUID));
end
eventID = idx(localID);
end

function fileIDs = localEventFileIDsPixel(BWA_AQ2)
fileIDs = strings(0,1);
fts = BWA_AQ2.fts1;
if isfield(fts, 'basic') && isfield(fts.basic, 'fileID') && ~isempty(fts.basic.fileID)
    fileIDs = string(fts.basic.fileID);
elseif isfield(fts, 'fileID') && ~isempty(fts.fileID)
    fileIDs = string(fts.fileID);
end
end

function fileID = localEventFileIDPixel(BWA_AQ2, eventID)
fileIDs = localEventFileIDsPixel(BWA_AQ2);
fileID = string(fileIDs(eventID));
end

function [imageSz, spatialRes, frameRateS] = localEventImagingMetaPixel(BWA_AQ2, eventID)
opts = BWA_AQ2.opts;
imageSz = double(opts.sz(1,:));
spatialRes = double(opts.spatialRes(1));
frameRateS = NaN;
if isfield(opts, 'frameRate') && ~isempty(opts.frameRate)
    frameRateS = double(opts.frameRate(1));
end
fileIDs = localCanonCellstrPixel(BWA_AQ2.fts1.fileID);
thisID = fileIDs{eventID};
rowIDs = unique(fileIDs, 'stable');
row = find(strcmp(rowIDs, thisID), 1, 'first');
if size(opts.sz,1) >= row
    imageSz = double(opts.sz(row,:));
end
if numel(opts.spatialRes) >= row
    spatialRes = double(opts.spatialRes(row));
end
if isfield(opts, 'frameRate') && numel(opts.frameRate) >= row
    frameRateS = double(opts.frameRate(row));
end
if numel(imageSz) < 3
    imageSz(3) = 1;
end
end

function c = localCanonCellstrPixel(v)
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

function zeroFrame = localBranchOriginFramePixel(result, BWA_AQ2, eventUID, eventID)
zeroFrame = NaN;
if isfield(result, 'metrics') && isstruct(result.metrics) && isfield(result.metrics, 'tables') && ...
        isstruct(result.metrics.tables) && isfield(result.metrics.tables, 'event_canonical') && ...
        istable(result.metrics.tables.event_canonical)
    T = result.metrics.tables.event_canonical;
    if ismember('event_uid', T.Properties.VariableNames) && ismember('branch_origin_frame', T.Properties.VariableNames)
        hit = find(string(T.event_uid) == string(eventUID), 1, 'first');
        if ~isempty(hit)
            zeroFrame = double(T.branch_origin_frame(hit));
        end
    end
end
if ~(isfinite(zeroFrame) && zeroFrame > 0) && isfield(BWA_AQ2.fts1, 'branch_origin') && ...
        isstruct(BWA_AQ2.fts1.branch_origin) && isfield(BWA_AQ2.fts1.branch_origin, 'frame')
    zeroFrame = double(BWA_AQ2.fts1.branch_origin.frame(eventID));
end
if ~(isfinite(zeroFrame) && zeroFrame > 0)
    zeroFrame = round(double(BWA_AQ2.fts1.loc.t0(eventID)));
end
end

function postOriginFrameMap = localPostOriginActivationFrameMapPixel(BWA_AQ2, eventID, imageSz, zeroFrame)
fts = BWA_AQ2.fts1;
sz = double(imageSz(:)');
sx = round(sz(1));
sy = round(sz(2));
sz3 = 1;
if numel(sz) >= 3 && isfinite(sz(3)) && sz(3) > 0
    sz3 = round(sz(3));
end
t1 = round(double(fts.loc.t1(eventID)));
postOriginFrameMap = nan(sx, sy);
for frameNow = max(round(double(zeroFrame)), 1):t1
    mask = localExplorerFrameMaskPixel(fts.loc.xSpaTemp{eventID}, frameNow, sx, sy, sz3);
    addMask = mask & ~isfinite(postOriginFrameMap);
    postOriginFrameMap(addMask) = frameNow;
end
end

function [firstFrameMap, t0, t1] = localFirstActivationFrameMapPixel(BWA_AQ2, eventID, imageSz)
fts = BWA_AQ2.fts1;
sz = double(imageSz(:)');
sx = round(sz(1));
sy = round(sz(2));
sz3 = 1;
if numel(sz) >= 3 && isfinite(sz(3)) && sz(3) > 0
    sz3 = round(sz(3));
end
t0 = round(double(fts.loc.t0(eventID)));
t1 = round(double(fts.loc.t1(eventID)));
firstFrameMap = nan(sx, sy);
for frameNow = t0:t1
    mask = localExplorerFrameMaskPixel(fts.loc.xSpaTemp{eventID}, frameNow, sx, sy, sz3);
    addMask = mask & ~isfinite(firstFrameMap);
    firstFrameMap(addMask) = frameNow;
end
end

function mask = localExplorerFrameMaskPixel(xSpaTemp_ID, frameNum, sx, sy, sz3)
limitStart = sx * sy * sz3 * (frameNum - 1) + 1;
limitEnd = limitStart + sx * sy * sz3;
linMovie = double(xSpaTemp_ID(:));
linFrame = linMovie(linMovie >= limitStart & linMovie <= limitEnd) - limitStart + 1;
linFrame = linFrame(isfinite(linFrame) & linFrame >= 1 & linFrame <= sx * sy);
mask = false(sx, sy);
mask(linFrame) = true;
end

function hFig = localPlotPixelActivationMap(out, visible, delayUnits)
hFig = figure('Name', 'EventPixelActivationMap', 'Color', 'w', 'NumberTitle', 'off', 'Visible', visible);
ax = axes('Parent', hFig);
hold(ax, 'on');
spatialRes = double(out.spatial_res);
sx = size(out.footprint_mask, 1);
sy = size(out.footprint_mask, 2);
bg = zeros(sx, sy, 3);
for c = 1:3
    plane = bg(:,:,c);
    plane(out.footprint_mask) = 0.72;
    bg(:,:,c) = plane;
end
image(ax, [0 sy * spatialRes], [0 sx * spatialRes], bg);
axis(ax, 'image');
set(ax, 'YDir', 'reverse');
xlabel(ax, 'x (um)');
ylabel(ax, 'y (um)');

switch lower(delayUnits)
    case 'frames'
        delayMap = out.delay_frames;
        cbLabel = 'First Activation Delay From Branch Origin (frames)';
    otherwise
        delayMap = out.delay_sec;
        cbLabel = 'First Activation Delay From Branch Origin (s)';
end

valid = isfinite(delayMap);
if any(valid(:))
    xData = [0 sy * spatialRes];
    yData = [0 sx * spatialRes];
    hIm = imagesc(ax, xData, yData, delayMap, 'AlphaData', 0.90 * valid);
    set(hIm, 'AlphaData', 0.90 * valid);
    colormap(ax, turbo(256));
    delayMin = 0;
    delayMax = max(delayMap(valid), [], 'omitnan');
    if ~(isfinite(delayMax) && delayMax > delayMin)
        delayMax = delayMin + 1;
    end
    caxis(ax, [delayMin delayMax]);
    cb = colorbar(ax);
    ylabel(cb, cbLabel);
end

title(ax, sprintf('%s | %s', char(out.file_id), char(out.event_uid)));
hold(ax, 'off');
end
