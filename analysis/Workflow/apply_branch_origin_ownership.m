function [canonicalOut, info] = apply_branch_origin_ownership(canonicalOut, BWA_AQ2, varargin)
%APPLY_BRANCH_ORIGIN_OWNERSHIP Deprecated local ownership rewrite.
% Canonical ingest in the cluster repo is now responsible for authoritative
% branch-origin ownership and primary-cell assignment. This helper remains
% only as a legacy diagnostic entry point and is not used by the normal
% local-analysis workflow.

ip = inputParser;
ip.addRequired('canonicalOut', @isstruct);
ip.addRequired('BWA_AQ2', @isstruct);
ip.addParameter('minBranchPixels', 5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
ip.addParameter('minBranchPersistence', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
ip.addParameter('wavefrontMode', 'prototype', @(x) any(strcmpi(string(x), ["prototype","local"])));
ip.parse(canonicalOut, BWA_AQ2, varargin{:});

info = struct('applied', false, 'nEvents', 0, 'nBranchOwnersResolved', 0, 'nOwnersReassigned', 0);
if ~(isfield(canonicalOut, 'subset') && isstruct(canonicalOut.subset) && ...
        isfield(canonicalOut.subset, 'events') && istable(canonicalOut.subset.events) && ...
        isfield(canonicalOut.subset, 'cells') && istable(canonicalOut.subset.cells) && ...
        isfield(canonicalOut.subset, 'event_cell_map') && istable(canonicalOut.subset.event_cell_map))
    return
end
if ~(isstruct(BWA_AQ2) && isfield(BWA_AQ2, 'fts1') && isfield(BWA_AQ2, 'opts') && ...
        isfield(BWA_AQ2.fts1, 'loc') && isfield(BWA_AQ2.fts1.loc, 'xSpaTemp') && ...
        isfield(BWA_AQ2.fts1.loc, 't0') && isfield(BWA_AQ2.fts1.loc, 't1') && ...
        isfield(BWA_AQ2.fts1, 'region') && isfield(BWA_AQ2.fts1.region, 'cell') && ...
        isfield(BWA_AQ2.fts1.region.cell, 'mask'))
    return
end

events = canonicalOut.subset.events;
cells = canonicalOut.subset.cells;
map = canonicalOut.subset.event_cell_map;
fts = BWA_AQ2.fts1;
opts = BWA_AQ2.opts;
fileIDsByEvent = localCanonCellstr(fts.fileID);
uniqueFiles = unique(fileIDsByEvent, 'stable');

heightE = height(events);
onsetOwner = repmat("", heightE, 1);
if ismember('primary_cell_uid', events.Properties.VariableNames)
    onsetOwner = string(events.primary_cell_uid);
else
    events.primary_cell_uid = repmat("", heightE, 1);
end
branchOwner = repmat("", heightE, 1);
onsetOriginX = NaN(heightE, 1); onsetOriginY = NaN(heightE, 1); onsetOriginFrame = NaN(heightE, 1);
branchOriginX = NaN(heightE, 1); branchOriginY = NaN(heightE, 1); branchOriginFrame = NaN(heightE, 1);
branchPeakX = NaN(heightE, 1); branchPeakY = NaN(heightE, 1); branchPeakFrame = NaN(heightE, 1);

for fi = 1:numel(uniqueFiles)
    thisFile = uniqueFiles{fi};
    eventIdx = find(strcmp(fileIDsByEvent, thisFile));
    if isempty(eventIdx)
        continue
    end
    fprintf('apply_branch_origin_ownership: file %d/%d %s (%d events)\n', ...
        fi, numel(uniqueFiles), thisFile, numel(eventIdx));
    drawnow;
    for ii = 1:numel(eventIdx)
        eventID = eventIdx(ii);
        [imageSzEvent, spatialResEvent] = localEventMetaForOwnership(BWA_AQ2, eventID, fileIDsByEvent);
        geomEvent = CalcEventFrameGeometry(fts.loc.xSpaTemp(eventID), fts.loc.t0(eventID), fts.loc.t1(eventID), imageSzEvent, 'eventIDs', 1);
        G = geomEvent{1};
        W = localOwnershipWavefront(BWA_AQ2, eventID, G, spatialResEvent, ip.Results);
        if ~isstruct(W) || any(~isfinite(W.branchOriginXY))
            continue
        end
        onsetOriginX(eventID) = localFieldXY(W, 'onsetOriginXY', 1);
        onsetOriginY(eventID) = localFieldXY(W, 'onsetOriginXY', 2);
        onsetOriginFrame(eventID) = localScalarField(W, 'onsetOriginFrame');
        branchOriginX(eventID) = localFieldXY(W, 'branchOriginXY', 1);
        branchOriginY(eventID) = localFieldXY(W, 'branchOriginXY', 2);
        branchOriginFrame(eventID) = localScalarField(W, 'branchOriginFrame');
        branchPeakX(eventID) = localFieldXY(W, 'branchPeakPointXY', 1);
        branchPeakY(eventID) = localFieldXY(W, 'branchPeakPointXY', 2);
        branchPeakFrame(eventID) = localScalarField(W, 'branchPeakFrame');
        branchOwner(eventID) = localResolveBranchOwnerCellUID(BWA_AQ2, cells, eventID, [branchOriginX(eventID) branchOriginY(eventID)]);
    end
end

events.onset_cell_uid = onsetOwner;
events.branch_origin_cell_uid = branchOwner;
events.onset_origin_x_px = onsetOriginX;
events.onset_origin_y_px = onsetOriginY;
events.onset_origin_frame = onsetOriginFrame;
events.branch_origin_x_px = branchOriginX;
events.branch_origin_y_px = branchOriginY;
events.branch_origin_frame = branchOriginFrame;
events.branch_peak_x_px = branchPeakX;
events.branch_peak_y_px = branchPeakY;
events.branch_peak_frame = branchPeakFrame;

if all(ismember({'event_uid','cell_uid','is_primary_cell'}, map.Properties.VariableNames))
    mapEvent = string(map.event_uid);
    mapCell = string(map.cell_uid);
    if ~ismember('onset_is_primary_cell', map.Properties.VariableNames)
        map.onset_is_primary_cell = logical(map.is_primary_cell);
    end
    eventIDs = string(events.event_uid);
    for k = 1:heightE
        newOwner = string(branchOwner(k));
        if strlength(newOwner) == 0
            continue
        end
        rows = mapEvent == eventIDs(k);
        if ~any(rows & mapCell == newOwner)
            continue
        end
        map.is_primary_cell(rows) = false;
        map.is_primary_cell(rows & mapCell == newOwner) = true;
        events.primary_cell_uid(k) = newOwner;
    end
end

canonicalOut.subset.events = events;
canonicalOut.subset.event_cell_map = map;
info.applied = true;
info.nEvents = height(events);
info.nBranchOwnersResolved = nnz(strlength(branchOwner) > 0);
info.nOwnersReassigned = nnz(string(events.primary_cell_uid) ~= onsetOwner);
end

function uid = localResolveBranchOwnerCellUID(BWA_AQ2, cellsTbl, eventID, pt)
uid = "";
if any(~isfinite(pt))
    return
end
fts = BWA_AQ2.fts1;
allMasks = fts.region.cell.mask;
if isempty(allMasks)
    return
end

candidateIdx = [];
if isfield(fts.region.cell, 'memberIdx') && size(fts.region.cell.memberIdx, 1) >= eventID
    candidateIdx = unique(double(fts.region.cell.memberIdx(eventID, :)));
    candidateIdx = candidateIdx(isfinite(candidateIdx) & candidateIdx >= 1 & candidateIdx <= numel(allMasks));
end
if isempty(candidateIdx)
    candidateIdx = (1:min(numel(allMasks), height(cellsTbl)))';
end

if istable(cellsTbl) && ismember('file_id', cellsTbl.Properties.VariableNames) && isfield(fts, 'fileID') && eventID <= numel(fts.fileID)
    fileID = string(fts.fileID{eventID});
    sameFile = find(string(cellsTbl.file_id) == fileID);
    if ~isempty(sameFile)
        candidateIdx = intersect(candidateIdx(:), sameFile(:), 'stable');
        if isempty(candidateIdx)
            candidateIdx = sameFile(:);
        end
    end
end

row = round(pt(2));
col = round(pt(1));
inside = [];
for ii = 1:numel(candidateIdx)
    idx = candidateIdx(ii);
    if idx < 1 || idx > numel(allMasks)
        continue
    end
    mask = logical(allMasks{idx});
    if isempty(mask) || row < 1 || row > size(mask, 1) || col < 1 || col > size(mask, 2)
        continue
    end
    if mask(row, col)
        inside(end+1,1) = idx; %#ok<AGROW>
    end
end

if isempty(inside)
    bestIdx = localNearestCellByCentroid(allMasks, candidateIdx, pt);
else
    bestIdx = localNearestCellByCentroid(allMasks, inside, pt);
end
if isempty(bestIdx) || bestIdx < 1 || bestIdx > height(cellsTbl) || ~ismember('cell_uid', cellsTbl.Properties.VariableNames)
    return
end
uid = string(cellsTbl.cell_uid(bestIdx));
end

function W = localOwnershipWavefront(BWA_AQ2, eventID, G, spatialResEvent, opts)
W = struct();
mode = string(opts.wavefrontMode);
if strcmpi(mode, "prototype")
    try
        proto = wavefrontComponentGraphPrototype(BWA_AQ2, eventID, ...
            'showFigure', false, ...
            'refinePixelWalk', 'off');
        if isstruct(proto) && isfield(proto, 'wavefront') && isstruct(proto.wavefront) && ...
                isfield(proto.wavefront, 'branchOriginXY') && any(isfinite(proto.wavefront.branchOriginXY))
            W = proto.wavefront;
            return
        end
    catch
    end
end
W = CalcWavefrontMetrics(G, spatialResEvent, ...
    'minValidFrames', 1, ...
    'minBranchPixels', opts.minBranchPixels, ...
    'minBranchPersistence', opts.minBranchPersistence);
end

function idx = localNearestCellByCentroid(allMasks, candidateIdx, pt)
idx = [];
bestD = inf;
for ii = 1:numel(candidateIdx)
    cidx = candidateIdx(ii);
    if cidx < 1 || cidx > numel(allMasks)
        continue
    end
    mask = logical(allMasks{cidx});
    [r, c] = find(mask);
    if isempty(r)
        continue
    end
    ctr = [mean(c), mean(r)];
    d = hypot(ctr(1) - pt(1), ctr(2) - pt(2));
    if d < bestD
        bestD = d;
        idx = cidx;
    end
end
end

function v = localScalarField(S, field)
v = NaN;
if isstruct(S) && isfield(S, field)
    x = S.(field);
    if isnumeric(x) && isscalar(x) && isfinite(x)
        v = double(x);
    end
end
end

function v = localFieldXY(S, field, col)
v = NaN;
if isstruct(S) && isfield(S, field)
    x = S.(field);
    if isnumeric(x) && numel(x) >= col && isfinite(x(col))
        v = double(x(col));
    end
end
end

function [imageSz, spatialRes] = localEventMetaForOwnership(BWA_AQ2, eventID, fileIDs)
opts = BWA_AQ2.opts;
imageSz = double(opts.sz(1,:));
spatialRes = double(opts.spatialRes(1));
if ~(isfield(BWA_AQ2, 'fts1') && isfield(BWA_AQ2.fts1, 'fileID') && eventID <= numel(BWA_AQ2.fts1.fileID))
    if numel(imageSz) < 3
        imageSz(3) = 1;
    end
    return
end
if nargin < 3 || isempty(fileIDs)
    fileIDs = localCanonCellstr(BWA_AQ2.fts1.fileID);
end
rowIDs = unique(fileIDs, 'stable');
row = find(strcmp(rowIDs, fileIDs{eventID}), 1, 'first');
if ~isempty(row)
    if size(opts.sz, 1) >= row
        imageSz = double(opts.sz(row,:));
    end
    if numel(opts.spatialRes) >= row && isfinite(opts.spatialRes(row))
        spatialRes = double(opts.spatialRes(row));
    end
end
if numel(imageSz) < 3
    imageSz(3) = 1;
end
end

function c = localCanonCellstr(v)
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
