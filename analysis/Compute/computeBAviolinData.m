function comp = computeBAviolinData(fts, opts, datTree, timeFrameRangeEpochs, memberlist, varargin)
% Compute-only stage for BAviolins.
% Returns normalized/event-by-cell data used by plotBAviolinPanels().

    p = inputParser;
    p.StructExpand = false;
    addRequired(p, 'fts', @isstruct);
    addRequired(p, 'opts', @isstruct);
    addRequired(p, 'datTree', @(x) ischar(x) || isstring(x));
    addRequired(p, 'timeFrameRangeEpochs', @isstruct);
    addRequired(p, 'memberlist', @isnumeric);
    addParameter(p, 'eventsTbl', table(), @(x) istable(x) || isempty(x));
    addParameter(p, 'filesTbl', table(), @(x) istable(x) || isempty(x));
    addParameter(p, 'limitRules', {}, @(x) iscell(x) || isempty(x));
    parse(p, fts, opts, datTree, timeFrameRangeEpochs, memberlist, varargin{:});

    fts = p.Results.fts;
    opts = p.Results.opts;
    datTree = char(string(p.Results.datTree));
    timeFrameRangeEpochs = p.Results.timeFrameRangeEpochs;
    memberlist = p.Results.memberlist;
    eventsTbl = p.Results.eventsTbl;
    filesTbl = p.Results.filesTbl;
    limitRules = p.Results.limitRules;

    dffAll = canonicalMetricFromEvents(datTree, eventsTbl, filesTbl);
    if isempty(dffAll)
        dffAll = FindFieldValue(fts, datTree);
    end
    dffRaw = dffAll;

    t0 = fts.loc.t0;
    fileID = [];
    if isfield(fts, 'fileID')
        fileID = fts.fileID;
    end
    memberlist = memberlist(:);
    t0 = t0(:);
    dffAll = dffAll(:);
    if ~isempty(limitRules)
        pass = limitMembershipByDataRuleSet(fts, limitRules);
        pass = logical(pass(:));
        n0 = min([numel(pass), numel(memberlist), numel(t0), numel(dffAll)]);
        if n0 > 0
            pass = pass(1:n0);
            memberlist = memberlist(1:n0);
            t0 = t0(1:n0);
            dffAll = dffAll(1:n0);
            if isfield(fts, 'fileID')
                fileID = string(fts.fileID(:));
                fileID = fileID(1:n0);
                fileID = fileID(pass);
            end
            memberlist = memberlist(pass);
            t0 = t0(pass);
            dffAll = dffAll(pass);
        end
    end
    lens = [numel(memberlist), numel(t0), numel(dffAll)];
    if ~isempty(fileID)
        lens(end+1) = numel(fileID);
    end
    n = min(lens);
    if isempty(n) || n == 0
        dffAll = zeros(0,1);
        t0 = zeros(0,1);
        memberlist = zeros(0,1);
        fileID = strings(0,1);
    else
        if numel(memberlist) ~= numel(t0) || numel(memberlist) ~= numel(dffAll) || (~isempty(fileID) && numel(fileID) ~= n)
            warning('computeBAviolinData:LengthMismatch', ...
                'Truncating vectors to common length (%d): memberlist=%d, t0=%d, metric=%d, fileID=%d', ...
                n, numel(memberlist), numel(t0), numel(dffAll), numel(fileID));
        end
        memberlist = memberlist(1:n);
        t0 = t0(1:n);
        dffAll = dffAll(1:n);
        if isempty(fileID)
            fileID = repmat("", n, 1);
        else
            fileID = string(fileID(:));
            fileID = fileID(1:n);
        end
    end
    eventIndex = (1:n)';

    n_nonfinite_input = sum(~isfinite(double(dffRaw(:))));
    negIndx = dffAll < 0;
    n_neg = sum(negIndx);
    if n_neg && shouldDropNegativeValues(datTree)
        warning('computeBAviolinData:NegativeValuesDropped', ...
            '%s: dropped %d negative values.', datTree, n_neg);
        dffAll(negIndx) = NaN;
    elseif n_neg && ~strcmpi(datTree, 'propagation.propangle')
        warning('computeBAviolinData:NegativeValuesKept', ...
            '%s: retained %d negative values (expected for this metric).', datTree, n_neg);
    end

    tOccuring = t0 * opts.frameRate(1) / 60;

    datTxt = extractAfter(datTree, '.');
    if isempty(datTxt)
        datTxt = datTree;
    end

    epMembIdx = epochMembershipFromFrames(t0, timeFrameRangeEpochs);
    nonEpoch = isnan(epMembIdx);
    tEpochs = timeFrameRangeEpochs.epochNames(epMembIdx(~nonEpoch));

    iEpochs = double(categorical(tEpochs, timeFrameRangeEpochs.epochNames, 'Ordinal', true));
    fileInEpoch = fileID(~nonEpoch);
    if isempty(iEpochs)
        T_summary = table([], [], strings(0,1), [], 'VariableNames', {'epoch', 'cell', 'file_id', 'mean_dff'});
    else
        T = table(iEpochs(:), memberlist(~nonEpoch), fileInEpoch(:), dffAll(~nonEpoch), ...
            'VariableNames', {'epoch', 'cell', 'file_id', 'dff'});
        T_summary = groupsummary(T, {'epoch','cell','file_id'}, 'mean', 'dff', ...
            'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
    end

    avgs = T_summary.mean_dff;
    if isempty(avgs)
        tCat_byCell = categorical([]);
        avgs_byCell = [];
        cells_byCell = [];
    else
        tCat = timeFrameRangeEpochs.epochNames(T_summary.epoch);
        tCat_byCell = tCat(~isnan(avgs));
        avgs_byCell = avgs(~isnan(avgs));
        cells_byCell = T_summary.cell(~isnan(avgs));
    end

    baselineCellVal = T_summary.mean_dff(T_summary.epoch == 1);
    baselineCellID = T_summary.cell(T_summary.epoch == 1);
    baselineFileID = string(T_summary.file_id(T_summary.epoch == 1));
    [outDFFs, normDFFs, outMember, outT0s, outEpochIndex, normByEvent] = normalizedByCellBaseLocal( ...
        dffAll, memberlist, fileID, t0, baselineCellVal, baselineCellID, baselineFileID, epMembIdx);
    tOccuringN = outT0s * opts.frameRate(1) / 60;

    [outDffs_c, avgsn_byCell, outMember_c, ~, outEpochIndex_c] = normalizedByCellBaseLocal( ...
        T_summary.mean_dff, T_summary.cell, string(T_summary.file_id), zeros(size(T_summary.mean_dff)), ...
        baselineCellVal, baselineCellID, baselineFileID, T_summary.epoch);
    if isempty(outEpochIndex_c)
        tCatn_byCell = categorical([]);
    else
        tCatn_byCell = timeFrameRangeEpochs.epochNames(outEpochIndex_c);
    end

    comp = struct();
    comp.datTxt = char(datTxt);
    comp.dffAll = dffAll;
    comp.normByEvent = normByEvent(:);
    comp.t0 = t0;
    comp.fileID = fileID(:);
    comp.filesTbl = filesTbl;
    comp.memberlist = memberlist(:);
    comp.eventIndex = eventIndex;
    comp.tOccuring = tOccuring;
    comp.epMembIdx = epMembIdx;
    comp.nonEpoch = nonEpoch;
    comp.tEpochs = tEpochs;
    comp.T_summary = T_summary;
    comp.avgs_byCell = avgs_byCell;
    comp.tCat_byCell = tCat_byCell;
    comp.cells_byCell = cells_byCell;

    comp.outDFFs = outDFFs;
    comp.normDFFs = normDFFs;
    comp.outMember = outMember;
    comp.outT0s = outT0s;
    comp.outEpochIndex = outEpochIndex;
    comp.tOccuringN = tOccuringN;

    comp.outDffs_c = outDffs_c;
    comp.avgsn_byCell = avgsn_byCell;
    comp.outMember_c = outMember_c;
    comp.outEpochIndex_c = outEpochIndex_c;
    comp.tCatn_byCell = tCatn_byCell;
    comp.cellsn_byCell = outMember_c;
    comp.n_negative_values = n_neg;
    comp.negative_values_dropped = (n_neg > 0) && shouldDropNegativeValues(datTree);
    comp.n_nonfinite_input = n_nonfinite_input;
    comp.n_non_epoch = sum(nonEpoch);
    comp.n_events_input = numel(dffRaw);
end

function tf = shouldDropNegativeValues(datTree)
% Metrics like propagation angle can legitimately be negative.
switch char(string(datTree))
    case 'propagation.propangle'
        tf = false;
    otherwise
        tf = true;
end
end

function dffAll = canonicalMetricFromEvents(datTree, eventsTbl, filesTbl)
dffAll = [];
if ~istable(eventsTbl) || isempty(eventsTbl)
    return
end

col = mapDatTreeToEventCol(datTree);
if isempty(col) || ~ismember(col, eventsTbl.Properties.VariableNames)
    return
end
if ~ismember('file_id', eventsTbl.Properties.VariableNames) || ~ismember('local_event_id', eventsTbl.Properties.VariableNames)
    return
end

E = eventsTbl;
if istable(filesTbl) && ~isempty(filesTbl) && ismember('file_id', filesTbl.Properties.VariableNames)
    fileOrder = cellstr(string(filesTbl.file_id));
    eFile = cellstr(string(E.file_id));
    rank = nan(height(E),1);
    for ii = 1:numel(fileOrder)
        rank(strcmp(eFile, fileOrder{ii})) = ii;
    end
    known = rank(~isnan(rank));
    if isempty(known)
        rank(isnan(rank)) = 1;
    else
        rank(isnan(rank)) = max(known) + 1;
    end
    E.file_rank = rank;
    E = sortrows(E, {'file_rank', 'local_event_id'});
else
    E = sortrows(E, {'file_id', 'local_event_id'});
end

dffAll = double(E.(col));
dffAll = dffAll(:);
end

function col = mapDatTreeToEventCol(datTree)
switch char(string(datTree))
    case 'basic.area'
        col = 'area_um2';
    case 'curve.duration'
        col = 'duration_sec';
    case 'curve.rise19'
        col = 'rise19_sec';
    case 'curve.fall91'
        col = 'fall91_sec';
    case 'curve.decayTau'
        col = 'decay_tau_sec';
    case 'curve.dffAUC'
        col = 'dff_auc_sec';
    case 'curve.dffMax2'
        col = 'dff_max2';
    case 'propagation.maxPropSpeed'
        col = 'max_prop_speed';
    otherwise
        col = '';
end
end

function epMembIdx = epochMembershipFromFrames(t0, timeFrameRangeEpochs)
    epMembIdx = NaN(size(t0));
    epFrames = timeFrameRangeEpochs.frameRange;
    nEpochs = size(epFrames, 1);
    for epochIndex = 1:nEpochs
        frameRange = epFrames(epochIndex,1:2);
        isEpochMember = t0 >= frameRange(1) & t0 <= frameRange(2);
        epMembIdx(isEpochMember) = epochIndex;
    end
end

function [outDffs, normDFFs, outMember, outT0s, outEpochs, normFull] = normalizedByCellBaseLocal(inDFFs, inMember, inFileID, inT0s, baselineVals, baselineCell, baselineFileID, epMembIdx)
    inDFFs = inDFFs(:)';
    inMember = inMember(:)';
    inFileID = string(inFileID(:))';
    inT0s = inT0s(:)';
    epMembIdx = epMembIdx(:)';
    baselineVals = baselineVals(:)';
    baselineCell = baselineCell(:)';
    baselineFileID = string(baselineFileID(:))';

    normDFFs = NaN(size(inDFFs));
    normFull = NaN(size(inDFFs));
    outMember = NaN(size(inMember));

    nBaseline = numel(baselineVals);
    for baseIdx = 1:nBaseline
        baseCell = baselineCell(baseIdx);
        baseFile = baselineFileID(baseIdx);
        belong = (inMember == baseCell) & (inFileID == baseFile);
        idx = find(belong);
        if isempty(idx)
            continue;
        end

        % Prefer baseline computed directly from this cell's baseline-epoch events.
        baseVal = bwaComputeBaselineValue(inDFFs, belong & epMembIdx == 1, 'log2_fold');
        if ~isfinite(baseVal)
            % Backward-compatible fallback if baseline epoch events are absent.
            baseVal = baselineVals(baseIdx);
        end

        vals = inDFFs(idx);
        normVals = bwaApplyTransform(vals, repmat(baseVal, size(vals)), 'log2_fold');

        normDFFs(idx) = normVals;
        normFull(idx) = normVals;
        outMember(idx) = inMember(idx);
    end

    keep = ~isnan(normDFFs);
    outDffs = inDFFs(keep);
    normDFFs = normDFFs(keep);
    outMember = outMember(keep);
    outT0s = inT0s(keep);
    outEpochs = epMembIdx(keep);
end
