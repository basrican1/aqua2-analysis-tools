function out = computeBWA_Aqua2PolarsBySliceData(BWA_AQ2, varargin)
% Compute metric tables used by BWA_Aqua2PolarsBySlice plotting.
%
% This function contains only data preparation (no plotting).

    p = inputParser;
    p.StructExpand = false;
    addRequired(p, 'BWA_AQ2', @isstruct);
    addParameter(p,'timeFrameRangeEpochs', [], @(x) isempty(x) || isstruct(x));
    addParameter(p,'timeRange', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p,'frameRange', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p,'EpochList', {}, @(x) isempty(x) || iscellstr(x) || isstring(x));
    addParameter(p,'memberlist', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p,'associatedCoverageTraces', {}, @(x) isempty(x) || iscell(x));
    addParameter(p,'memberCoverage', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p,'filesTbl', table(), @(x) istable(x) || isempty(x));
    addParameter(p,'eventsTbl', table(), @(x) istable(x) || isempty(x));
    addParameter(p,'dataLimitRules', {}, @(x) iscell(x) || isempty(x));
    addParameter(p,'includeStandardFeatures', true, @islogical);
    addParameter(p,'filePlots', true, @islogical);
    addParameter(p,'combineEpochNames', {}, @(x) ischar(x) || isstring(x) || iscell(x));
    addParameter(p,'combinedEpochLabel', 'post', @(x) ischar(x) || isstring(x));
    parse(p, BWA_AQ2, varargin{:});

    ctx = prepareLegacyAqua2Context(BWA_AQ2, ...
        'timeFrameRangeEpochs', p.Results.timeFrameRangeEpochs, ...
        'timeRange', p.Results.timeRange, ...
        'frameRange', p.Results.frameRange, ...
        'EpochList', p.Results.EpochList, ...
        'memberlist', p.Results.memberlist, ...
        'associatedCoverageTraces', p.Results.associatedCoverageTraces, ...
        'memberCoverage', p.Results.memberCoverage, ...
        'filesTbl', p.Results.filesTbl, ...
        'eventsTbl', p.Results.eventsTbl);
    BWA_AQ2 = ctx.BWA_AQ2;
    opts = ctx.opts;
    timeFrameRangeEpochs = ctx.timeFrameRangeEpochs;
    memberlist = ctx.memberlistRow3;
    memberCoverage = ctx.memberCoverage;
    associatedCoverageTraces = ctx.associatedCoverageTraces;
    combineEpochNames = string(cellstr(string(p.Results.combineEpochNames(:))));
    combineEpochNames = combineEpochNames(strlength(combineEpochNames) > 0);
    combinedEpochLabel = char(string(p.Results.combinedEpochLabel));

    allFileIDs = BWA_AQ2.fts1.fileID;
    uniqueFiles = unique(allFileIDs, "stable");
    nFiles = numel(uniqueFiles);
    subs = numSubplotsForN(nFiles);

    areaPercentiles = [20 100];
    if ~isempty(p.Results.dataLimitRules)
        areaPercentiles = [0 100];
    end
    ftsDataList = {...
        'basic.area', ...
        'curve.duration', ...
        'curve.rise19', ...
        'curve.fall91', ...
        'curve.decayTau', ...
        'curve.dffAUC', ...
        'curve.dffMax2', ...
        'propagation.maxPropSpeed'};

    labels = {};
    saveValA = zeros(nFiles, 0);
    saveValC = zeros(nFiles, 0);
    saveSigA = false(nFiles, 0);
    saveSigC = false(nFiles, 0);
    entityValuesA = cell(0, 1);
    entityValuesC = cell(0, 1);
    entityValuesAByFile = cell(0, 1);
    entityValuesCByFile = cell(0, 1);
    feature = 0;
    computeCache = struct();

    if p.Results.includeStandardFeatures
        for ff = 1:numel(ftsDataList)
            datTree = ftsDataList{ff};
            [labels, saveValA, saveValC, saveSigA, saveSigC, entityValuesA, entityValuesC, entityValuesAByFile, entityValuesCByFile, feature, computeCache] = appendBySliceMetric( ...
                labels, saveValA, saveValC, saveSigA, saveSigC, entityValuesA, entityValuesC, entityValuesAByFile, entityValuesCByFile, feature, ...
                BWA_AQ2.fts1, opts, datTree, timeFrameRangeEpochs, memberlist, ...
                areaPercentiles, p.Results.filePlots, '', computeCache, p.Results.dataLimitRules, ...
                combineEpochNames, combinedEpochLabel, p.Results.filesTbl, p.Results.eventsTbl);
        end
    end

    % Derived metrics beyond canonical event fields.
    try
        TempFts = struct();
        TempFts.loc.t0 = BWA_AQ2.fts1.loc.t0;
        TempFts.fileID = BWA_AQ2.fts1.fileID;
        TempFts.propagation.pCell = memberCoverage;
        [labels, saveValA, saveValC, saveSigA, saveSigC, entityValuesA, entityValuesC, entityValuesAByFile, entityValuesCByFile, feature, computeCache] = appendBySliceMetric( ...
            labels, saveValA, saveValC, saveSigA, saveSigC, entityValuesA, entityValuesC, entityValuesAByFile, entityValuesCByFile, feature, ...
            TempFts, opts, 'propagation.pCell', timeFrameRangeEpochs, memberlist, ...
            {}, p.Results.filePlots, 'footprint', computeCache, p.Results.dataLimitRules, ...
            combineEpochNames, combinedEpochLabel, p.Results.filesTbl, p.Results.eventsTbl);
    catch ME
        warning('computeBWA_Aqua2PolarsBySliceData:DerivedMetricFailed', ...
            'Failed adding derived spider metric footprint: %s', ME.message);
    end

    try
        spatialRes = localFirstSpatialRes(opts);
        [propdistVal, angleVal] = calcpropdistAngleLocal(BWA_AQ2.fts1, spatialRes);
        TempFts = struct();
        TempFts.loc.t0 = BWA_AQ2.fts1.loc.t0;
        TempFts.fileID = BWA_AQ2.fts1.fileID;
        TempFts.propagation.propdist = propdistVal;
        TempFts.propagation.propangle = angleVal;
        [labels, saveValA, saveValC, saveSigA, saveSigC, entityValuesA, entityValuesC, entityValuesAByFile, entityValuesCByFile, feature, computeCache] = appendBySliceMetric( ...
            labels, saveValA, saveValC, saveSigA, saveSigC, entityValuesA, entityValuesC, entityValuesAByFile, entityValuesCByFile, feature, ...
            TempFts, opts, 'propagation.propdist', timeFrameRangeEpochs, memberlist, ...
            {}, p.Results.filePlots, 'propdist', computeCache, p.Results.dataLimitRules, ...
            combineEpochNames, combinedEpochLabel, p.Results.filesTbl, p.Results.eventsTbl);
        [labels, saveValA, saveValC, saveSigA, saveSigC, entityValuesA, entityValuesC, entityValuesAByFile, entityValuesCByFile, feature, computeCache] = appendBySliceMetric( ...
            labels, saveValA, saveValC, saveSigA, saveSigC, entityValuesA, entityValuesC, entityValuesAByFile, entityValuesCByFile, feature, ...
            TempFts, opts, 'propagation.propangle', timeFrameRangeEpochs, memberlist, ...
            {}, p.Results.filePlots, 'propangle', computeCache, p.Results.dataLimitRules, ...
            combineEpochNames, combinedEpochLabel, p.Results.filesTbl, p.Results.eventsTbl);
    catch ME
        warning('computeBWA_Aqua2PolarsBySliceData:DerivedMetricFailed', ...
            'Failed adding derived spider metric propagation (dist/angle): %s', ME.message);
    end

    try
        nCent = centroidCompatibleEventCount(BWA_AQ2.fts1);
        xSpaTemp = BWA_AQ2.fts1.loc.xSpaTemp(1:nCent);
        t0 = BWA_AQ2.fts1.loc.t0(1:nCent);
        t1 = BWA_AQ2.fts1.loc.t1(1:nCent);
        fileID = BWA_AQ2.fts1.fileID(1:nCent);
        propGrow = BWA_AQ2.fts1.propagation.propGrow(1:nCent);
        allCentroids = CalcCentroids(xSpaTemp, t0, t1, opts.sz(1,:), propGrow);
        [~, ~, distTot] = DistFromCentroids(allCentroids, opts.spatialRes(1));
        TempFts = struct();
        TempFts.loc.t0 = t0;
        TempFts.fileID = fileID;
        TempFts.centroid.distTot = distTot;
        [labels, saveValA, saveValC, saveSigA, saveSigC, entityValuesA, entityValuesC, entityValuesAByFile, entityValuesCByFile, feature, computeCache] = appendBySliceMetric( ...
            labels, saveValA, saveValC, saveSigA, saveSigC, entityValuesA, entityValuesC, entityValuesAByFile, entityValuesCByFile, feature, ...
            TempFts, opts, 'centroid.distTot', timeFrameRangeEpochs, memberlist, ...
            {}, p.Results.filePlots, 'centroid', computeCache, p.Results.dataLimitRules, ...
            combineEpochNames, combinedEpochLabel, p.Results.filesTbl, p.Results.eventsTbl);
    catch ME
        warning('computeBWA_Aqua2PolarsBySliceData:DerivedMetricFailed', ...
            'Failed adding derived spider metric centroid: %s', ME.message);
    end

    [figs, FreqMeansA, FreqMeansC] = nEventsPerMin(BWA_AQ2.fts1, opts, memberlist, timeFrameRangeEpochs, areaPercentiles, 1, ...
        'plotFigures', p.Results.filePlots, 'dataLimitRules', p.Results.dataLimitRules);
    feature = feature + 1;
    labels{feature,1} = 'freq';
    saveValA(:, feature) = cellfun(@(x) combineEpochValue(x, timeFrameRangeEpochs, combineEpochNames), FreqMeansA);
    saveValC(:, feature) = cellfun(@(x) combineEpochValue(x, timeFrameRangeEpochs, combineEpochNames), FreqMeansC);
    saveSigA(:, feature) = false(nFiles, 1);
    saveSigC(:, feature) = false(nFiles, 1);
    entityValuesA{feature,1} = [];
    entityValuesC{feature,1} = [];
    entityValuesAByFile{feature,1} = cell(nFiles, 1);
    entityValuesCByFile{feature,1} = cell(nFiles, 1);
    if ~isempty(figs)
        close([figs{:}]);
    end

    out = struct();
    out.labels = labels;
    out.valuesA = saveValA;
    out.valuesC = saveValC;
    out.sigMaskA = saveSigA;
    out.sigMaskC = saveSigC;
    out.entityValuesA = entityValuesA;
    out.entityValuesC = entityValuesC;
    out.entityValuesAByFile = entityValuesAByFile;
    out.entityValuesCByFile = entityValuesCByFile;
    out.uniqueFiles = uniqueFiles;
    out.fileLabels = bwaMakeSliceLabels(uniqueFiles, p.Results.filesTbl);
    out.subplotGrid = subs;
    out.areaPercentiles = areaPercentiles;
    out.timeFrameRangeEpochs = timeFrameRangeEpochs;
    out.memberlist = memberlist;
    out.memberCoverage = memberCoverage;
    out.associatedCoverageTraces = associatedCoverageTraces;
    out.combineEpochNames = cellstr(combineEpochNames(:))';
    out.combinedEpochLabel = combinedEpochLabel;
end

function [labels, saveValA, saveValC, saveSigA, saveSigC, entityValuesA, entityValuesC, entityValuesAByFile, entityValuesCByFile, feature, cacheOut] = appendBySliceMetric( ...
    labels, saveValA, saveValC, saveSigA, saveSigC, entityValuesA, entityValuesC, entityValuesAByFile, entityValuesCByFile, feature, fts, opts, datTree, timeFrameRangeEpochs, memberlist, areaPercentiles, plotFigures, labelOverride, cacheIn, dataLimitRules, combineEpochNames, combinedEpochLabel, filesTbl, eventsTbl)

    if nargin < 13 || ~isstruct(cacheIn)
        cacheIn = struct();
    end
    cacheOut = cacheIn;

    if plotFigures
        [figs, means_alln, means_celln] = checkBySlice(fts, opts, datTree, ...
            timeFrameRangeEpochs, memberlist, areaPercentiles, 'plotFigures', plotFigures, 'limitRules', dataLimitRules);
        sig_alln = repmat({false}, numel(means_alln), 1);
        sig_celln = repmat({false}, numel(means_celln), 1);
        pooledA = [];
        pooledC = [];
        pooledAByFile = cell(numel(means_alln), 1);
        pooledCByFile = cell(numel(means_celln), 1);
    else
        [means_alln, means_celln, sig_alln, sig_celln, pooledA, pooledC, pooledAByFile, pooledCByFile, cacheOut] = checkBySliceComputeOnlyFast(fts, opts, datTree, ...
            timeFrameRangeEpochs, memberlist, areaPercentiles, cacheIn, dataLimitRules, combineEpochNames, filesTbl, eventsTbl);
        figs = {};
    end
    if ~isempty(figs)
        close([figs{:}]);
    end
    if nargin >= 12 && ~isempty(labelOverride)
        datTxt = char(string(labelOverride));
    else
        datTxt = extractAfter(datTree, '.');
        if isempty(datTxt)
            datTxt = datTree;
        end
        datTxt = char(datTxt);
    end
    feature = feature + 1;
    labels{feature,1} = datTxt; %#ok<AGROW>
    saveValA(:, feature) = cellfun(@(x) combineEpochValue(x, timeFrameRangeEpochs, combineEpochNames), means_alln); %#ok<AGROW>
    saveValC(:, feature) = cellfun(@(x) combineEpochValue(x, timeFrameRangeEpochs, combineEpochNames), means_celln); %#ok<AGROW>
    saveSigA(:, feature) = cellfun(@logical, sig_alln); %#ok<AGROW>
    saveSigC(:, feature) = cellfun(@logical, sig_celln); %#ok<AGROW>
    entityValuesA{feature,1} = pooledA; %#ok<AGROW>
    entityValuesC{feature,1} = pooledC; %#ok<AGROW>
    entityValuesAByFile{feature,1} = pooledAByFile; %#ok<AGROW>
    entityValuesCByFile{feature,1} = pooledCByFile; %#ok<AGROW>
end

function [meansA, meansC, sigA, sigC, pooledA, pooledC, pooledAByFile, pooledCByFile, cacheOut] = checkBySliceComputeOnlyFast(fts, opts, datTree, timeFrameRangeEpochs, memberlist, anaRulz, cacheIn, dataLimitRules, combineEpochNames, filesTbl, eventsTbl)
% Compute-only version of checkBySlice with cache reuse across metrics.

cacheOut = cacheIn;
[cacheKey, ~] = getOrBuildCache(fts, memberlist, anaRulz, cacheIn, dataLimitRules);
cacheOut.(cacheKey.name) = cacheKey.data;
rawVals = FindFieldValue(fts, datTree);
rawVals = rawVals(cacheKey.data.passMask);

nFiles = numel(cacheKey.data.uniqueFiles);
meansA = cell(nFiles, 1);
meansC = cell(nFiles, 1);
sigA = cell(nFiles, 1);
sigC = cell(nFiles, 1);
pooledA = [];
pooledC = [];
pooledAByFile = cell(nFiles, 1);
pooledCByFile = cell(nFiles, 1);
for fileNum = 1:nFiles
    eventIndexes = cacheKey.data.fileEventIdx{fileNum};
    cellmemberSubset = cacheKey.data.memberFiltered(eventIndexes);
    t0subset = cacheKey.data.t0Filtered(eventIndexes);
    fileIDsubset = string(cacheKey.data.fileIDsFiltered(eventIndexes));
    dataSubset = rawVals(eventIndexes);

    TempAQ2 = struct();
    TempAQ2.fts.loc.t0 = t0subset;
    TempAQ2.fts.fileID = fileIDsubset;
    TempAQ2.fts = PutFieldValue(dataSubset, datTree, TempAQ2.fts, 0);
    thisFileID = string(cacheKey.data.uniqueFiles{fileNum});
    if istable(filesTbl) && ~isempty(filesTbl) && ismember('file_id', filesTbl.Properties.VariableNames)
        fileMask = ismember(string(filesTbl.file_id), thisFileID);
        thisFilesTbl = filesTbl(fileMask, :);
    else
        thisFilesTbl = table();
    end
    if istable(eventsTbl) && ~isempty(eventsTbl) && ismember('file_id', eventsTbl.Properties.VariableNames)
        eventMask = ismember(string(eventsTbl.file_id), thisFileID);
        thisEventsTbl = eventsTbl(eventMask, :);
    else
        thisEventsTbl = table();
    end
    comp = computeBAviolinData(TempAQ2.fts, opts, datTree, timeFrameRangeEpochs, cellmemberSubset, ...
        'filesTbl', thisFilesTbl, 'eventsTbl', thisEventsTbl);
    [means_alln, means_celln] = computeEpochMeansFromComp(comp, timeFrameRangeEpochs);
    [sig_alln, sig_celln] = computeCnoSignificanceFromComp(comp, timeFrameRangeEpochs, combineEpochNames);
    [valsA, valsC] = extractCombinedEpochValues(comp, timeFrameRangeEpochs, combineEpochNames);
    meansA{fileNum} = means_alln;
    meansC{fileNum} = means_celln;
    sigA{fileNum} = sig_alln;
    sigC{fileNum} = sig_celln;
    pooledA = [pooledA; valsA(:)]; %#ok<AGROW>
    pooledC = [pooledC; valsC(:)]; %#ok<AGROW>
    pooledAByFile{fileNum} = valsA(:);
    pooledCByFile{fileNum} = valsC(:);
end
end

function [out, hit] = getOrBuildCache(fts, memberlist, anaRulz, cacheIn, dataLimitRules)
allFileIDs = fts.fileID;
allT0 = fts.loc.t0;
if isnumeric(anaRulz) && numel(anaRulz) == 2
    anaRulz = AddLimitingRuleset({}, 'percentiles', 'area', anaRulz);
end
if ~isempty(dataLimitRules)
    if isempty(anaRulz)
        anaRulz = dataLimitRules;
    else
        anaRulz = [anaRulz; dataLimitRules];
    end
end
ruleSig = rulesSignature(anaRulz);
fieldName = matlab.lang.makeValidName(['rule_' ruleSig]);
hit = false;
if isstruct(cacheIn) && isfield(cacheIn, fieldName)
    c = cacheIn.(fieldName);
    if isfield(c, 'nEvents') && c.nEvents == numel(allFileIDs)
        out.name = fieldName;
        out.data = c;
        hit = true;
        return
    end
end

pass = limitMembershipByDataRuleSet(fts, anaRulz);
if size(memberlist,1) == 3
    memberV = memberlist(3, :);
else
    memberV = memberlist;
end
memberV = memberV(:);
pass = pass(:);
allT0 = allT0(:);

allFileIDsCell = cellstr(string(allFileIDs(:)));
uniqueFiles = unique(allFileIDsCell, "stable");
fileIDsFiltered = allFileIDsCell(pass);

fileEventIdx = cell(numel(uniqueFiles), 1);
for i = 1:numel(uniqueFiles)
    thisFileID = uniqueFiles{i};
    fileEventIdx{i} = find(ismember(fileIDsFiltered, thisFileID));
end

c = struct();
c.nEvents = numel(allFileIDsCell);
c.passMask = pass;
c.uniqueFiles = uniqueFiles;
c.fileEventIdx = fileEventIdx;
c.fileIDsFiltered = fileIDsFiltered;
c.memberFiltered = memberV(pass);
c.t0Filtered = allT0(pass);

out.name = fieldName;
out.data = c;
end

function sig = rulesSignature(anaRulz)
if isempty(anaRulz)
    sig = 'none';
else
    % Current compute path uses either {} or area-percentile rule sets.
    % Keep key coarse and stable to avoid type-shape issues in signatures.
    sig = 'limited';
end
end

function [means_alln, means_celln] = computeEpochMeansFromComp(comp, timeFrameRangeEpochs)
nEpoch = numel(timeFrameRangeEpochs.epochNames);
means_alln = NaN(1, nEpoch);
means_celln = NaN(1, nEpoch);

if isfield(comp, 'normDFFs') && isfield(comp, 'outEpochIndex')
    for e = 1:nEpoch
        x = double(comp.normDFFs(comp.outEpochIndex == e));
        x = x(isfinite(x));
        if ~isempty(x)
            means_alln(e) = mean(x);
        end
    end
end

if isfield(comp, 'avgsn_byCell') && isfield(comp, 'outEpochIndex_c')
    for e = 1:nEpoch
        x = double(comp.avgsn_byCell(comp.outEpochIndex_c == e));
        x = x(isfinite(x));
        if ~isempty(x)
            means_celln(e) = mean(x);
        end
    end
end
end

function [sig_alln, sig_celln] = computeCnoSignificanceFromComp(comp, timeFrameRangeEpochs, combineEpochNames)
nEpoch = numel(timeFrameRangeEpochs.epochNames);
sig_alln = false;
sig_celln = false;
if nEpoch < 2
    return
end
idx = resolveCombinedEpochIdx(timeFrameRangeEpochs, combineEpochNames);
if isempty(idx)
    return
end
if isfield(comp, 'normDFFs') && isfield(comp, 'outEpochIndex')
    x = double(comp.normDFFs(ismember(comp.outEpochIndex, idx)));
    sig_alln = isSigVsZero(x);
end
if isfield(comp, 'avgsn_byCell') && isfield(comp, 'outEpochIndex_c')
    x = double(comp.avgsn_byCell(ismember(comp.outEpochIndex_c, idx)));
    sig_celln = isSigVsZero(x);
end
end

function [valsA, valsC] = extractCombinedEpochValues(comp, timeFrameRangeEpochs, combineEpochNames)
valsA = [];
valsC = [];
idx = resolveCombinedEpochIdx(timeFrameRangeEpochs, combineEpochNames);
if isempty(idx)
    return
end
if isfield(comp, 'normDFFs') && isfield(comp, 'outEpochIndex')
    maskA = ismember(double(comp.outEpochIndex(:)), idx);
    valsA = double(comp.normDFFs(maskA));
    valsA = valsA(isfinite(valsA));
end
if isfield(comp, 'avgsn_byCell') && isfield(comp, 'outEpochIndex_c') && isfield(comp, 'outMember_c')
    epochC = double(comp.outEpochIndex_c(:));
    memberC = double(comp.outMember_c(:));
    valsCAll = double(comp.avgsn_byCell(:));
    maskC = ismember(epochC, idx) & isfinite(memberC) & isfinite(valsCAll);
    if any(maskC)
        Tc = table(memberC(maskC), valsCAll(maskC), 'VariableNames', {'cell_id','value'});
        G = findgroups(Tc.cell_id);
        valsC = splitapply(@(x) mean(x, 'omitnan'), Tc.value, G);
        valsC = double(valsC(:));
        valsC = valsC(isfinite(valsC));
    end
end
end

function tf = isSigVsZero(x)
x = x(:);
x = x(isfinite(x));
tf = false;
if numel(x) < 3
    return
end
try
    p = signrank(x, 0, 'method', 'approximate');
    tf = isfinite(p) && p < 0.05;
catch
    try
        [~, p] = ttest(x, 0);
        tf = isfinite(p) && p < 0.05;
    catch
        tf = false;
    end
end
end

function v = combineEpochValue(x, timeFrameRangeEpochs, combineEpochNames)
idx = resolveCombinedEpochIdx(timeFrameRangeEpochs, combineEpochNames);
if isempty(x) || isempty(idx)
    v = NaN;
    return
end
idx = idx(idx >= 1 & idx <= numel(x));
if isempty(idx)
    v = NaN;
    return
end
v = mean(double(x(idx)), 'omitnan');
end

function idx = resolveCombinedEpochIdx(timeFrameRangeEpochs, combineEpochNames)
idx = [];
if ~isstruct(timeFrameRangeEpochs) || ~isfield(timeFrameRangeEpochs, 'epochNames') || isempty(timeFrameRangeEpochs.epochNames)
    return
end
epochNames = string(timeFrameRangeEpochs.epochNames(:));
if isempty(combineEpochNames)
    idx = find((1:numel(epochNames))' >= 2);
    return
end
for i = 1:numel(combineEpochNames)
    hit = find(epochNames == string(combineEpochNames(i)), 1, 'first');
    if ~isempty(hit)
        idx(end+1,1) = hit; %#ok<AGROW>
    end
end
idx = unique(idx, 'stable');
end

function [propdistVal, angleVal] = calcpropdistAngleLocal(fts, spatialRes)
if nargin < 2 || ~(isnumeric(spatialRes) && isfinite(spatialRes) && spatialRes > 0)
    spatialRes = 1;
end
t0 = fts.loc.t0;
nEvents = size(t0, 2);
propdistVal = zeros(1, nEvents);
angleVal = zeros(1, nEvents);
for eventID = 1:nEvents
    [~, propdistVal(eventID), angleVal(eventID)] = calcPropagationLocal(fts, eventID, spatialRes);
end
end

function [prop, valStr_distance, angle] = calcPropagationLocal(fts, eventnum, spatialRes)
if nargin < 3 || ~(isnumeric(spatialRes) && isfinite(spatialRes) && spatialRes > 0)
    spatialRes = 1;
end
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
        valStr_distance = valStr_distance + pdist([previousCenter(1) previousCenter(2); testx testy], 'euclidean');
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

function sr = localFirstSpatialRes(opts)
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

function n = centroidCompatibleEventCount(fts)
n = 0;
try
    counts = [ ...
        numel(fts.loc.t0), ...
        numel(fts.loc.t1), ...
        numel(fts.fileID), ...
        numel(fts.loc.xSpaTemp), ...
        numel(fts.propagation.propGrow)];
    counts = counts(counts > 0);
    if isempty(counts)
        n = 0;
    else
        n = min(counts);
    end
catch
    n = 0;
end
if n < 1
    error('computeBWA_Aqua2PolarsBySliceData:BadCentroidInputs', ...
        'Centroid metric could not determine a compatible event count.');
end
end
