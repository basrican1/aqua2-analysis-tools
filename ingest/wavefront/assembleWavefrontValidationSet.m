function out = assembleWavefrontValidationSet(source, varargin)
%ASSEMBLEWAVEFRONTVALIDATIONSET Build a compact event set for wavefront review.
%   OUT = assembleWavefrontValidationSet(RESULT)
%   OUT = assembleWavefrontValidationSet(RESULT, 'TopMulticell', 4, 'NSingleCell', 3)
%
%   Returns:
%     OUT.summaryTable   - one row per selected event with basic spread and motion metrics
%     OUT.multicellTable - candidate multicell events ranked by recruited cells then area
%     OUT.singleTable    - representative single-cell events across the area range
%     OUT.selectedEvents - selected local event IDs

BWA_setupPaths(fileparts(mfilename('fullpath')));

p = inputParser;
p.StructExpand = false;
addParameter(p, 'TopMulticell', 4, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'NSingleCell', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'IncludeEvents', [], @(x) isempty(x) || isnumeric(x) || iscell(x) || isstring(x) || ischar(x));
addParameter(p, 'RefinePixelWalk', 'auto', @(x) (islogical(x) && isscalar(x)) || any(strcmpi(string(x), ["auto","on","off"])));
parse(p, varargin{:});

[result, BWA_AQ2, sourceLabel] = localResolveValidationSource(source);
eventTbl = localBuildEventSummaryTable(result);

multicell = eventTbl(eventTbl.n_participating_cells >= 2, :);
if ~isempty(multicell)
    sortVars = intersect({'n_recruited_cells','n_participating_cells','area_um2','area_px2','local_event_id'}, multicell.Properties.VariableNames, 'stable');
    sortDir = repmat({'descend'}, 1, numel(sortVars));
    multicell = sortrows(multicell, sortVars, sortDir);
end

single = eventTbl(eventTbl.n_participating_cells <= 1, :);
singlePicked = localPickRepresentativeSingles(single, round(double(p.Results.NSingleCell)));

selectedRows = table();
if ~isempty(p.Results.IncludeEvents)
    includeRows = localResolveIncludedRows(eventTbl, p.Results.IncludeEvents);
    if ~isempty(includeRows)
        selectedRows = [selectedRows; includeRows]; %#ok<AGROW>
    end
end
if ~isempty(multicell)
    selectedRows = [selectedRows; multicell(1:min(height(multicell), round(double(p.Results.TopMulticell))), :)]; %#ok<AGROW>
end
if ~isempty(singlePicked)
    selectedRows = [selectedRows; singlePicked]; %#ok<AGROW>
end
if ~isempty(selectedRows)
    selectedRows = unique(selectedRows, 'rows', 'stable');
end

summary = localComputeSelectedMetrics(result, BWA_AQ2, selectedRows, p.Results.RefinePixelWalk);

out = struct();
out.sourceLabel = sourceLabel;
out.eventTable = eventTbl;
out.multicellTable = multicell;
out.singleTable = singlePicked;
out.selectedEvents = selectedRows;
out.summaryTable = summary;
end

function [result, BWA_AQ2, sourceLabel] = localResolveValidationSource(source)
sourceLabel = "";
if isstruct(source) && isfield(source, 'metrics') && isfield(source, 'pathIn')
    result = source;
    sourceLabel = string(result.pathIn);
elseif ischar(source) || isstring(source)
    S = load(char(string(source)), 'result');
    if ~isfield(S, 'result')
        error('assembleWavefrontValidationSet:MissingResult', 'Expected MAT file containing variable ''result''.');
    end
    result = S.result;
    if isfield(result, 'pathIn')
        sourceLabel = string(result.pathIn);
    else
        sourceLabel = string(source);
    end
else
    error('assembleWavefrontValidationSet:BadSource', 'Expected a local-analysis result or result MAT path.');
end

if isfield(result, 'derived') && isfield(result.derived, 'context') && isfield(result.derived.context, 'BWA_AQ2') && ...
        isstruct(result.derived.context.BWA_AQ2) && isfield(result.derived.context.BWA_AQ2, 'fts1') && isfield(result.derived.context.BWA_AQ2, 'opts')
    BWA_AQ2 = result.derived.context.BWA_AQ2;
elseif isfield(result, 'pathIn') && strlength(string(result.pathIn)) > 0
    loaded = loadClusterOutputs(char(string(result.pathIn)), 'attachLegacyMeta', true, 'attachEventMeta', true);
    BWA_AQ2 = loaded.BWA_AQ2;
else
    error('assembleWavefrontValidationSet:MissingAqua2', 'Could not resolve BWA_AQ2 from the supplied source.');
end
end

function T = localBuildEventSummaryTable(result)
P = result.metrics.tables.event_partition;
S = result.metrics.tables.event_shape;

eventIDFromUID = arrayfun(@localEventIdFromUID, string(P.event_uid));
cellCol = localFirstPresent(P.Properties.VariableNames, {'partition_cell_uid','cell_uid'});
if isempty(cellCol)
    error('assembleWavefrontValidationSet:MissingCellColumn', 'event_partition must include partition_cell_uid or cell_uid.');
end
cellIDs = string(P.(cellCol));
eventUID = string(P.event_uid);

[G, grp] = findgroups(eventUID);
nCells = splitapply(@(x) numel(unique(x(x ~= ""))), cellIDs, G);
nRecruited = zeros(max(G), 1);
if ismember('is_recruited_partition', P.Properties.VariableNames)
    isRec = logical(double(P.is_recruited_partition));
    nRecruited = splitapply(@(cid, rec) numel(unique(cid(rec & cid ~= ""))), cellIDs, isRec, G);
end

base = unique(P(:, intersect({'event_uid','file_id'}, P.Properties.VariableNames, 'stable')), 'rows', 'stable');
base.local_event_id = arrayfun(@localEventIdFromUID, string(base.event_uid));
base.n_participating_cells = nCells;
base.n_recruited_cells = nRecruited;

shapeKeep = intersect({'event_uid','area_um2','area_px2'}, S.Properties.VariableNames, 'stable');
if ~isempty(shapeKeep)
    base = outerjoin(base, unique(S(:, shapeKeep), 'rows', 'stable'), 'Keys', 'event_uid', 'MergeKeys', true, 'Type', 'left');
end

T = sortrows(base, 'local_event_id');
end

function Tpick = localPickRepresentativeSingles(singleTbl, nPick)
Tpick = singleTbl([],:);
if isempty(singleTbl) || nPick <= 0
    return
end
scoreCol = localFirstPresent(singleTbl.Properties.VariableNames, {'area_um2','area_px2','local_event_id'});
vals = double(singleTbl.(scoreCol));
[vals, order] = sort(vals, 'ascend');
singleTbl = singleTbl(order, :);
if numel(vals) <= nPick
    Tpick = singleTbl;
    return
end
q = linspace(0.2, 0.8, nPick);
idx = max(1, min(numel(vals), round(1 + q * (numel(vals) - 1))));
idx = unique(idx, 'stable');
while numel(idx) < nPick
    cand = setdiff(round(linspace(1, numel(vals), nPick)), idx, 'stable');
    idx = unique([idx, cand], 'stable');
end
idx = idx(1:min(nPick, numel(idx)));
Tpick = singleTbl(idx, :);
end

function summary = localComputeSelectedMetrics(result, BWA_AQ2, selectedRows, refinePixelWalk)
summary = table();
if isempty(selectedRows)
    return
end

protoSource = result;
if ~isfield(protoSource, 'derived') || ~isstruct(protoSource.derived)
    protoSource.derived = struct();
end
if ~isfield(protoSource.derived, 'context') || ~isstruct(protoSource.derived.context)
    protoSource.derived.context = struct();
end
protoSource.derived.context.BWA_AQ2 = BWA_AQ2;

fts = BWA_AQ2.fts1;
rows = cell(height(selectedRows), 1);
globalIndex = localGlobalEventIndex(BWA_AQ2);
for i = 1:height(selectedRows)
    hit = selectedRows(i, :);
    eventUID = string(hit.event_uid(1));
    if ~isKey(globalIndex, char(eventUID))
        continue
    end
    eventID = globalIndex(char(eventUID));
    [imageSz, spatialRes] = localEventImagingMeta(BWA_AQ2, eventID);
    Gtmp = CalcEventFrameGeometry(fts.loc.xSpaTemp(eventID), fts.loc.t0(eventID), fts.loc.t1(eventID), imageSz, 'eventIDs', 1);
    G = Gtmp{1};
    motion = CalcCentroidAndWavefrontMetrics(fts.loc.xSpaTemp(eventID), fts.loc.t0(eventID), fts.loc.t1(eventID), imageSz, spatialRes, 'eventIDs', 1);
    proto = wavefrontComponentGraphPrototype(protoSource, eventID, 'showFigure', false, 'refinePixelWalk', refinePixelWalk);

    cohort = "single";
    if hit.n_participating_cells(1) >= 2
        cohort = "multicell";
    end

    rows{i} = table( ...
        double(hit.local_event_id(1)), ...
        eventUID, ...
        string(hit.file_id(1)), ...
        string(cohort), ...
        double(hit.n_participating_cells(1)), ...
        double(hit.n_recruited_cells(1)), ...
        localMaybeScalar(hit, 'area_um2'), ...
        localMaybeScalar(hit, 'area_px2'), ...
        double(motion.centroid.distTot(1)), ...
        double(motion.centroid.netDisplacementUm(1)), ...
        double(motion.centroid.tortuosity(1)), ...
        double(proto.wavefront.branchOriginFrame), ...
        double(proto.wavefront.branchPeakFrame), ...
        double(proto.wavefront.wavefront_max_reach_um), ...
        double(proto.wavefront.wavefront_frontier_length_um), ...
        double(proto.wavefront.wavefront_tortuosity), ...
        logical(proto.isMulticellEvent), ...
        logical(proto.usedRefinedPixelWalk), ...
        'VariableNames', {'event_id','event_uid','file_id','cohort','n_participating_cells','n_recruited_cells', ...
            'area_um2','area_px2','centroid_dist_um','centroid_net_disp_um','centroid_tortuosity', ...
            'branch_origin_frame','branch_peak_frame','wavefront_reach_um','wavefront_frontier_um','wavefront_tortuosity', ...
            'is_multicell_gate','used_refined_pixel_walk'});
    rows{i}.notes = "";
    rows{i}.geometry_frames = numel(G.frames);
end
rows = rows(~cellfun(@isempty, rows));
if isempty(rows)
    return
end
summary = vertcat(rows{:});
end

function outRows = localResolveIncludedRows(eventTbl, includeSpec)
outRows = eventTbl([],:);
if isempty(includeSpec)
    return
end
if isnumeric(includeSpec)
    vals = double(includeSpec(:));
    keep = ismember(double(eventTbl.local_event_id), vals);
    outRows = eventTbl(keep, :);
    return
end
uids = string(includeSpec(:));
keep = ismember(string(eventTbl.event_uid), uids);
outRows = eventTbl(keep, :);
end

function v = localMaybeScalar(T, name)
v = NaN;
if ismember(name, T.Properties.VariableNames)
    x = T.(name)(1);
    if isnumeric(x)
        v = double(x);
    end
end
end

function name = localFirstPresent(names, candidates)
name = '';
for i = 1:numel(candidates)
    if any(strcmp(names, candidates{i}))
        name = candidates{i};
        return
    end
end
end

function eid = localEventIdFromUID(uid)
tok = regexp(char(uid), '_E(\d+)$', 'tokens', 'once');
if isempty(tok)
    eid = NaN;
else
    eid = str2double(tok{1});
end
end

function idxMap = localGlobalEventIndex(BWA_AQ2)
idxMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
fts = BWA_AQ2.fts1;
fileIDs = [];
if isfield(fts, 'basic') && isfield(fts.basic, 'fileID') && ~isempty(fts.basic.fileID)
    fileIDs = string(fts.basic.fileID);
elseif isfield(fts, 'fileID') && ~isempty(fts.fileID)
    fileIDs = string(fts.fileID);
end
if isempty(fileIDs)
    return
end
counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:numel(fileIDs)
    fid = char(fileIDs(i));
    if ~isKey(counts, fid)
        counts(fid) = 0;
    end
    counts(fid) = counts(fid) + 1;
    uid = sprintf('%s_E%d', fid, counts(fid));
    idxMap(uid) = i;
end
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
