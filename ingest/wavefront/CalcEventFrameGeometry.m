function geom = CalcEventFrameGeometry(xSpaTemp, allt0, allt1, imageSz, varargin)
%CALCEVENTFRAMEGEOMETRY Build per-frame event geometry from xSpaTemp once.
%   GEOM = CALCEVENTFRAMEGEOMETRY(XSPATEMP, ALLT0, ALLT1, IMAGESZ)
%   returns one geometry struct per event with frame masks, new-pixel masks,
%   per-frame centroids, and frontier candidate pixels derived from xSpaTemp.
%
%   Optional Name-Value:
%     'eventIDs' : subset of event IDs to compute (default all)

p = inputParser;
p.StructExpand = false;
addParameter(p, 'eventIDs', [], @(x) isempty(x) || isnumeric(x));
parse(p, varargin{:});

nEvents = size(xSpaTemp, 2);
geom = cell(1, nEvents);
eventIDs = p.Results.eventIDs;
if isempty(eventIDs)
    eventIDs = 1:nEvents;
end

for eventID = eventIDs(:)'
    if eventID < 1 || eventID > nEvents
        continue
    end
    geom{eventID} = oneEventGeometry(xSpaTemp{eventID}, allt0(eventID), allt1(eventID), imageSz);
end
end

function G = oneEventGeometry(xSpaTemp_ID, t0, t1, imageSz)
frames = (double(t0):double(t1))';
nFrames = numel(frames);
frameMasks = cell(nFrames, 1);
newMasks = cell(nFrames, 1);
frameIdx = cell(nFrames, 1);
newIdx = cell(nFrames, 1);
centroidXY = NaN(nFrames, 2);
newCentroidXY = NaN(nFrames, 2);
frontierPts = cell(nFrames, 1);
validFrontier = false(nFrames, 1);

prevMask = [];
for k = 1:nFrames
    frameNum = frames(k);
    idx = localEventPixInFrame(frameNum, xSpaTemp_ID, imageSz);
    frameIdx{k} = idx;
    M = false(imageSz(1), imageSz(2));
    M(idx) = true;
    frameMasks{k} = M;

    if k == 1
        newM = M;
    else
        newM = M & ~prevMask;
    end
    prevMask = M;
    newMasks{k} = newM;
    newIdx{k} = find(newM);

    if any(M(:))
        [r, c] = find(M);
        centroidXY(k, :) = [mean(c), mean(r)];
    end
    if any(newM(:))
        [r, c] = find(newM);
        newCentroidXY(k, :) = [mean(c), mean(r)];
        frontierPts{k} = [c, r];
        validFrontier(k) = true;
    else
        frontierPts{k} = zeros(0, 2);
    end
end

firstFrontier = find(validFrontier, 1, 'first');
originXY = [NaN NaN];
if ~isempty(firstFrontier)
    originXY = newCentroidXY(firstFrontier, :);
end

G = struct();
G.frames = frames;
G.t0 = double(t0);
G.t1 = double(t1);
G.frameMasks = frameMasks;
G.newMasks = newMasks;
G.frameIdx = frameIdx;
G.newIdx = newIdx;
G.centroidXY = centroidXY;
G.newCentroidXY = newCentroidXY;
G.frontierPts = frontierPts;
G.validFrontier = validFrontier;
G.originXY = originXY;
G.firstFrontierIdx = firstFrontier;
end

function LinIndexIntoFrame = localEventPixInFrame(eventFrame, xSpaTemp_ID, imageSz)
xPix = imageSz(1);
yPix = imageSz(2);
zPix = imageSz(3);

limitStart = xPix * yPix * zPix * (eventFrame - 1) + 1;
limitEnd = limitStart + xPix * yPix * zPix;

LinIndexIntoMovie = xSpaTemp_ID(xSpaTemp_ID >= limitStart & xSpaTemp_ID <= limitEnd);
LinIndexIntoFrame = LinIndexIntoMovie - limitStart + 1;
LinIndexIntoFrame = double(LinIndexIntoFrame(:));
LinIndexIntoFrame = LinIndexIntoFrame(isfinite(LinIndexIntoFrame) & LinIndexIntoFrame >= 1 & LinIndexIntoFrame <= xPix * yPix);
end
