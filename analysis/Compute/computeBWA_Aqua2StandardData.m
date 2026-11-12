function out = computeBWA_Aqua2StandardData(ctx, varargin)
% Compute standard derived metrics without plotting.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'ctx', @isstruct);
addParameter(p, 'dataLimitRules', {}, @(x) iscell(x) || isempty(x));
addParameter(p, 'includeEventMotion', false, @islogical);
parse(p, ctx, varargin{:});
ctx = p.Results.ctx;
dataLimitRules = p.Results.dataLimitRules;
includeEventMotion = p.Results.includeEventMotion;

if ~isstruct(ctx) || ~isfield(ctx, 'BWA_AQ2') || ~isfield(ctx, 'opts')
    error('computeBWA_Aqua2StandardData:BadContext', 'ctx must be output from prepareLegacyAqua2Context().');
end

BWA_AQ2 = ctx.BWA_AQ2;
opts = ctx.opts;
tfre = ctx.timeFrameRangeEpochs;
memberlist = ctx.memberlistRow3;
fts = BWA_AQ2.fts1;

filesTbl = table();
eventsTbl = table();
mapTbl = table();
if isfield(ctx, 'filesTbl'), filesTbl = ctx.filesTbl; end
if isfield(ctx, 'eventsTbl'), eventsTbl = ctx.eventsTbl; end
if isfield(ctx, 'eventCellMapTbl'), mapTbl = ctx.eventCellMapTbl; end

out = struct();
out.status = 'ok';
out.errors = struct();
out.baviolin = struct();

trees = { ...
    'basic.area', ...
    'curve.duration', ...
    'curve.rise19', ...
    'curve.fall91', ...
    'curve.decayTau', ...
    'curve.dffAUC', ...
    'curve.dffMax2', ...
    'propagation.maxPropSpeed'};
for i = 1:numel(trees)
    t = trees{i};
        key = matlab.lang.makeValidName(strrep(t, '.', '_'));
    try
        out.baviolin.(key) = computeBAviolinData(fts, opts, t, tfre, memberlist, ...
            'eventsTbl', eventsTbl, 'filesTbl', filesTbl, 'limitRules', dataLimitRules);
    catch ME
        out.errors.(key) = ME.message;
    end
end

if ~isempty(filesTbl) && ~isempty(eventsTbl) && ~isempty(mapTbl)
    try
        [eventsF, mapF] = filterCanonicalEventsAndMap(eventsTbl, mapTbl, dataLimitRules);
        out.frequency = computeNEventsPerMinCanonical(eventsTbl, mapTbl, filesTbl, tfre, ...
            'groupBy', {'file_id'}, 'areaPercentiles', [0 100]);
        if ~isempty(dataLimitRules)
            out.frequency = computeNEventsPerMinCanonical(eventsF, mapF, filesTbl, tfre, ...
                'groupBy', {'file_id'}, 'areaPercentiles', [0 100]);
        end
    catch ME
        out.errors.frequency = ME.message;
    end
end

try
    covThreshold = 0.1;
    [coverage, t0_subset, epMembIdx, cellIdx, incs] = CellCoverage_calc( ...
        fts, opts, tfre, memberlist, covThreshold, ctx.associatedCoverageTraces);
    out.coverage = struct();
    out.coverage.threshold = covThreshold;
    out.coverage.coverage = coverage(:);
    out.coverage.t0_subset = t0_subset(:);
    out.coverage.epoch_idx = epMembIdx(:);
    out.coverage.cell_idx = cellIdx(:);
    out.coverage.included = logical(incs(:));
catch ME
    out.errors.coverage = ME.message;
end

try
    TempFts = struct();
    TempFts.loc.t0 = fts.loc.t0;
    if isfield(fts, 'fileID')
        TempFts.fileID = fts.fileID;
    end
    TempFts.propagation.pCell = ctx.memberCoverage;
    out.footprint = computeBAviolinData(TempFts, opts, 'propagation.pCell', tfre, memberlist, ...
        'eventsTbl', eventsTbl, 'filesTbl', filesTbl, 'limitRules', dataLimitRules);
catch ME
    out.errors.footprint = ME.message;
end

try
    spatialRes = localFirstSpatialRes(opts);
    [propdistVal, angleVal] = calcpropdistAngleLocal(fts, spatialRes);
    TempFts = struct();
    TempFts.loc.t0 = fts.loc.t0;
    if isfield(fts, 'fileID')
        TempFts.fileID = fts.fileID;
    end
    TempFts.propagation.propdist = propdistVal;
    TempFts.propagation.propangle = angleVal;
    out.propagation = struct();
    out.propagation.distance_values = propdistVal(:);
    out.propagation.angle_values = angleVal(:);
    out.propagation.distance = computeBAviolinData(TempFts, opts, 'propagation.propdist', tfre, memberlist, ...
        'eventsTbl', eventsTbl, 'filesTbl', filesTbl, 'limitRules', dataLimitRules);
    out.propagation.angle = computeBAviolinData(TempFts, opts, 'propagation.propangle', tfre, memberlist, ...
        'eventsTbl', eventsTbl, 'filesTbl', filesTbl, 'limitRules', dataLimitRules);
catch ME
    out.errors.propagation = ME.message;
end

if includeEventMotion
    try
        motion = CalcCentroidAndWavefrontMetrics( ...
            fts.loc.xSpaTemp, fts.loc.t0, fts.loc.t1, opts.sz(1,:), opts.spatialRes(1), ...
            'fileIDByEvent', fts.fileID, ...
            'imageSzByFile', opts.sz, ...
            'spatialResByFile', opts.spatialRes);
        TempFts = struct();
        TempFts.loc.t0 = fts.loc.t0;
        if isfield(fts, 'fileID')
            TempFts.fileID = fts.fileID;
        end

        out.centroid = struct();
        out.centroid.allCentroids = motion.centroid.allCentroids;
        out.centroid.dist = motion.centroid.dist;
        out.centroid.jumps = motion.centroid.jumps;
        out.centroid.distTot = motion.centroid.distTot(:);
        out.centroid.netDisplacementUm = motion.centroid.netDisplacementUm(:);
        out.centroid.tortuosity = motion.centroid.tortuosity(:);
        TempFts.centroid.distTot = motion.centroid.distTot;
        TempFts.centroid.netDisplacementUm = motion.centroid.netDisplacementUm;
        TempFts.centroid.tortuosity = motion.centroid.tortuosity;
        out.centroid.distComp = computeBAviolinData(TempFts, opts, 'centroid.distTot', tfre, memberlist, ...
            'eventsTbl', eventsTbl, 'filesTbl', filesTbl, 'limitRules', dataLimitRules);
        out.centroid.netDisplacementComp = computeBAviolinData(TempFts, opts, 'centroid.netDisplacementUm', tfre, memberlist, ...
            'eventsTbl', eventsTbl, 'filesTbl', filesTbl, 'limitRules', dataLimitRules);
        out.centroid.tortuosityComp = computeBAviolinData(TempFts, opts, 'centroid.tortuosity', tfre, memberlist, ...
            'eventsTbl', eventsTbl, 'filesTbl', filesTbl, 'limitRules', dataLimitRules);

        wf = motion.wavefront;
        out.wavefront = struct();
        out.wavefront.event = wf;
        out.wavefront.max_reach_um = cellfun(@(s) localFieldOrNaN(s, 'wavefront_max_reach_um'), wf(:));
        out.wavefront.frontier_length_um = cellfun(@(s) localFieldOrNaN(s, 'wavefront_frontier_length_um'), wf(:));
        out.wavefront.net_displacement_um = cellfun(@(s) localFieldOrNaN(s, 'wavefront_net_displacement_um'), wf(:));
        out.wavefront.tortuosity = cellfun(@(s) localFieldOrNaN(s, 'wavefront_tortuosity'), wf(:));
        out.wavefront.branch_origin_frame = cellfun(@(s) localFieldOrNaN(s, 'branchOriginFrame'), wf(:));
        out.wavefront.branch_peak_frame = cellfun(@(s) localFieldOrNaN(s, 'branchPeakFrame'), wf(:));
        TempFts.wavefront.maxReach = out.wavefront.max_reach_um(:)';
        TempFts.wavefront.frontierLength = out.wavefront.frontier_length_um(:)';
        TempFts.wavefront.netDisplacement = out.wavefront.net_displacement_um(:)';
        TempFts.wavefront.tortuosity = out.wavefront.tortuosity(:)';
        out.wavefront.maxReachComp = computeBAviolinData(TempFts, opts, 'wavefront.maxReach', tfre, memberlist, ...
            'eventsTbl', eventsTbl, 'filesTbl', filesTbl, 'limitRules', dataLimitRules);
        out.wavefront.frontierLengthComp = computeBAviolinData(TempFts, opts, 'wavefront.frontierLength', tfre, memberlist, ...
            'eventsTbl', eventsTbl, 'filesTbl', filesTbl, 'limitRules', dataLimitRules);
        out.wavefront.netDisplacementComp = computeBAviolinData(TempFts, opts, 'wavefront.netDisplacement', tfre, memberlist, ...
            'eventsTbl', eventsTbl, 'filesTbl', filesTbl, 'limitRules', dataLimitRules);
        out.wavefront.tortuosityComp = computeBAviolinData(TempFts, opts, 'wavefront.tortuosity', tfre, memberlist, ...
            'eventsTbl', eventsTbl, 'filesTbl', filesTbl, 'limitRules', dataLimitRules);
    catch ME
        out.errors.centroid = ME.message;
        out.errors.wavefront = ME.message;
    end
end
end

function v = localFieldOrNaN(s, field)
v = NaN;
if isstruct(s) && isfield(s, field)
    x = s.(field);
    if isnumeric(x) && isscalar(x) && isfinite(x)
        v = double(x);
    end
end
end

function [propdistVal, angleVal] = calcpropdistAngleLocal(fts, spatialRes)
t0 = fts.loc.t0;
nEvents = size(t0, 2);
propdistVal = zeros(1, nEvents);
angleVal = zeros(1, nEvents);
for eventID = 1:nEvents
    [~, propdistVal(eventID), angleVal(eventID)] = calcPropagationLocal(fts, eventID, spatialRes);
end
end

function [prop, valStr_distance, angle] = calcPropagationLocal(fts, eventnum, spatialRes)
propGrowAll = fts.propagation.propGrow;
propGrow_i = propGrowAll{1, eventnum};
propShrinkAll = fts.propagation.propShrink;
propShrink_i = propShrinkAll{1, eventnum};

nEventFrames = size(propGrow_i, 1);
if nEventFrames == 0
    prop = [];
    valStr_distance = 0;
    angle = NaN;
    return
end
coordX = zeros(4, nEventFrames);
coordY = zeros(4, nEventFrames);
centerTrack = zeros(2, nEventFrames);
valStr_distance = 0;
for eFrame = 1:nEventFrames
    if eFrame == 1
        prevNorth = [0 0];
        prevWest = [0 0];
        prevSouth = [0 0];
        prevEast = [0 0];
        previousCenter = [0 0];
    else
        prevNorth = [coordX(1, eFrame-1) coordY(1, eFrame-1)];
        prevWest = [coordX(2, eFrame-1) coordY(2, eFrame-1)];
        prevSouth = [coordX(3, eFrame-1) coordY(3, eFrame-1)];
        prevEast = [coordX(4, eFrame-1) coordY(4, eFrame-1)];
        previousCenter = [centerTrack(1,eFrame-1) centerTrack(2,eFrame-1)];
    end

    delta_up = propGrow_i(eFrame, 1) + propShrink_i(eFrame, 1);
    delta_down = -propGrow_i(eFrame, 2) - propShrink_i(eFrame, 2);
    delta_left = -propGrow_i(eFrame, 3) - propShrink_i(eFrame, 3);
    delta_right = propGrow_i(eFrame, 4) + propShrink_i(eFrame, 4);

    if prevNorth(1) > previousCenter(1), dxN = delta_right;
    elseif prevNorth(1) < previousCenter(1), dxN = delta_left;
    else, dxN = 0; end

    if prevSouth(1) > previousCenter(1), dxS = delta_right;
    elseif prevSouth(1) < previousCenter(1), dxS = delta_left;
    else, dxS = 0; end

    if prevWest(2) > previousCenter(2), dyW = delta_up;
    elseif prevWest(2) < previousCenter(2), dyW = delta_down;
    else, dyW = 0; end

    if prevEast(2) > previousCenter(2), dyE = delta_up;
    elseif prevEast(2) < previousCenter(2), dyE = delta_down;
    else, dyE = 0; end

    newNorth = [prevNorth(1)+dxN prevNorth(2)+delta_up];
    newWest = [prevWest(1)+delta_left prevWest(2)+dyW];
    newSouth = [prevSouth(1)+dxS prevSouth(2)+delta_down];
    newEast = [prevEast(1)+delta_right prevEast(2)+dyE];

    polyX = [newNorth(1) newWest(1) newSouth(1) newEast(1)];
    polyY = [newNorth(2) newWest(2) newSouth(2) newEast(2)];
    polyin = polyshape(polyX, polyY, 'Simplify', false);
    [testx, testy] = centroid(polyin);
    if isnan(testx)
        cx = previousCenter(1);
        cy = previousCenter(2);
    else
        valStr_distance = valStr_distance + spatialRes * pdist([previousCenter(1) previousCenter(2); testx testy], 'euclidean');
        cx = testx;
        cy = testy;
    end
    centerTrack(1,eFrame) = cx;
    centerTrack(2,eFrame) = cy;
    coordX(:, eFrame) = [newNorth(1) newWest(1) newSouth(1) newEast(1)];
    coordY(:, eFrame) = [newNorth(2) newWest(2) newSouth(2) newEast(2)];
end
angle = atan2d(cy, cx);
prop = [];
end

function [eventsOut, mapOut] = filterCanonicalEventsAndMap(eventsTbl, mapTbl, rules)
eventsOut = eventsTbl;
mapOut = mapTbl;
if isempty(rules) || isempty(eventsTbl)
    return
end
pass = true(height(eventsTbl),1);
for i = 1:size(rules,1)
    if size(rules,2) < 3 || isempty(rules{i,3})
        continue
    end
    limitType = char(string(rules{i,1}));
    dataType = char(string(rules{i,2}));
    value = rules{i,3};
    col = mapDataTypeToEventCol(dataType);
    if isempty(col) || ~ismember(col, eventsTbl.Properties.VariableNames)
        continue
    end
    x = double(eventsTbl.(col));
    pass = pass & applyRuleToVector(x, limitType, value);
end
eventsOut = eventsTbl(pass,:);
if ismember('event_uid', eventsOut.Properties.VariableNames) && ismember('event_uid', mapTbl.Properties.VariableNames)
    mapOut = mapTbl(ismember(mapTbl.event_uid, eventsOut.event_uid), :);
end
end

function pass = applyRuleToVector(x, limitType, value)
pass = false(size(x));
finiteMask = isfinite(x);
if ~any(finiteMask), return; end
switch lower(limitType)
    case {'percentile','percentiles'}
        lo = prctile(x(finiteMask), value(1));
        hi = prctile(x(finiteMask), value(2));
    case {'threshold','thresholds'}
        lo = value(1); hi = value(2);
    otherwise
        lo = -inf; hi = inf;
end
pass = finiteMask & x >= lo & x <= hi;
end

function col = mapDataTypeToEventCol(dataType)
switch char(string(dataType))
    case 'area', col = 'area';
    case 'amplitude', col = 'dff_max2';
    case 'AUC', col = 'dff_auc';
    case 'duration', col = 'duration';
    case 'rise', col = 'rise19';
    case 'fall', col = 'fall91';
    case 'decayTau', col = 'decay_tau';
    case 'maxspeed', col = 'max_prop_speed';
    otherwise, col = '';
end
end

function spatialRes = localFirstSpatialRes(opts)
spatialRes = 1;
if isstruct(opts) && isfield(opts, 'spatialRes') && ~isempty(opts.spatialRes)
    sr = opts.spatialRes;
    if isnumeric(sr)
        spatialRes = double(sr(1));
    elseif iscell(sr) && ~isempty(sr{1})
        spatialRes = double(sr{1}(1));
    end
end
if ~isfinite(spatialRes) || spatialRes <= 0
    spatialRes = 1;
end
end
