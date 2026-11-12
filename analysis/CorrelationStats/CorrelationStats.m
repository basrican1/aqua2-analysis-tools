function out = CorrelationStats(source, varargin)
%CORRELATIONSTATS Compute metric correlation tables from local analysis outputs.
%
% Usage:
%   out = CorrelationStats(result)
%   out = CorrelationStats(result, 'unitLevel', 'slice', 'metricFamilies', {'shape','synchrony'})
%   out = CorrelationStats(result, 'unitLevel', 'cell', 'metricFamilies', {'shape','recruitment'}, 'epochName', 'baseline')
%   out = CorrelationStats(result, 'unitLevel', 'event', 'metricFamilies', {'shape','intensity'}, 'epochName', 'baseline')
%   out = CorrelationStats(result, 'unitLevel', 'slice', 'metricFamilies', {'shape','synchrony'}, ...
%       'dataView', 'norm', 'combineEpochNames', {'ExperimentalEpoch1','ExperimentalEpoch2'}, 'combinedEpochLabel', 'CombinedPostEpochs')
%
% Name-value:
%   'unitLevel'    : 'slice' (default), 'cell', 'event', or 'animal'
%   'metricFamilies' : 'all' (default), or any of
%                    {'activity','shape','kinetics','intensity','propagation',
%                     'containment','recruitment','synchrony','soma'}
%                    when 'all' is used, families default by unit level
%   'dataView'     : 'raw' (default), 'norm', or 'recommended'
%                    recommended is the preferred mode for mixed metric sets
%                    and repeated-post analyses. It chooses normalization per metric:
%                    canonical positive magnitudes -> log2_fold
%                    event/cell/recruitment/network or signed metrics -> delta
%   'epochName'    : optional epoch filter; default '' means include all epochs as separate rows
%   'combineEpochNames' : optional cellstr/string list of epochs to combine into one repeated-measure
%                         post epoch before pivoting (example: {'ExperimentalEpoch1','ExperimentalEpoch2'}).
%                         This is the recommended pattern when ExperimentalEpoch1/ExperimentalEpoch2 are
%                         repeated post-drug epochs rather than separate conditions.
%   'combinedEpochLabel': label to assign the combined epoch (default 'combined_post')
%   'combineEpochMethod': aggregation across the selected epochs, currently 'mean' (default)
%   'metricSelection' : 'default' (default), 'broad', or 'strict'
%                    broad   = exploratory; keep more secondary metrics
%                    default = balanced everyday surface
%                    strict  = presentation-ready; prune more redundancy/sparse metrics
%   'method'       : 'Spearman' (default) or 'Pearson'
%   'minN'         : minimum pairwise observations for a reported correlation (default 3)
%   'includeOverlappingDerived' : internal legacy option for overlapping non-primary metrics (default false)
%
% Output fields:
%   out.longTable   : tidy table of unit/epoch/metric/value rows
%   out.wideTable   : one row per unit+epoch, one column per metric
%   out.metricNames : correlation metric column names used in wideTable
%   out.R           : correlation matrix
%   out.P           : p-value matrix
%   out.N           : pairwise observation-count matrix
%   out.pairTable   : long-form upper-triangle correlation table with FDR q-values

result = load_local_analysis_result(source);

p = inputParser;
p.StructExpand = false;
addParameter(p, 'unitLevel', 'slice', @(x) any(strcmpi(string(x), ["cell","slice","event","animal"])));
addParameter(p, 'metricFamilies', 'all', @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'dataView', 'raw', @(x) any(strcmpi(string(x), ["raw","norm","recommended"])));
addParameter(p, 'epochName', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'combineEpochNames', {}, @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'combinedEpochLabel', 'combined_post', @(x) ischar(x) || isstring(x));
addParameter(p, 'combineEpochMethod', 'mean', @(x) any(strcmpi(string(x), "mean")));
addParameter(p, 'metricSelection', 'default', @(x) any(strcmpi(string(x), ["default","broad","strict"])));
addParameter(p, 'method', 'Spearman', @(x) any(strcmpi(string(x), ["Spearman","Pearson"])));
addParameter(p, 'minN', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 2);
addParameter(p, 'includeOverlappingDerived', false, @islogical);
parse(p, varargin{:});

unitLevel = lower(char(string(p.Results.unitLevel)));
dataView = lower(char(string(p.Results.dataView)));
metricFamilies = normalizeMetricFamilies(p.Results.metricFamilies);
if isempty(metricFamilies)
    metricFamilies = {'all'};
end
selectionMode = lower(char(string(p.Results.metricSelection)));
method = lower(char(string(p.Results.method)));
epochName = string(p.Results.epochName);
combineEpochNames = string(cellstr(string(p.Results.combineEpochNames(:))));
combineEpochNames = combineEpochNames(strlength(combineEpochNames) > 0);
combinedEpochLabel = string(p.Results.combinedEpochLabel);
combineEpochMethod = lower(char(string(p.Results.combineEpochMethod)));

reg = table();
if isfield(result, 'metrics') && isstruct(result.metrics) && isfield(result.metrics, 'registry') && istable(result.metrics.registry)
    reg = result.metrics.registry;
end

L = cell(0,1);
metaTables = struct('canonical', struct(), 'derived', struct(), 'registry', struct());
metaTables.event = struct();

if isfield(result, 'canonical') && isstruct(result.canonical)
    [Lc, metaTables.canonical] = collectCanonicalLong(result, unitLevel, dataView);
    L = [L; Lc]; %#ok<AGROW>
    [Ld, metaTables.derived] = collectDerivedLong(result, unitLevel, dataView, p.Results.includeOverlappingDerived);
    L = [L; Ld]; %#ok<AGROW>
end

if strcmp(unitLevel, 'event')
    [Le, metaTables.event] = collectEventLocalLong(result, dataView);
    L = [L; Le]; %#ok<AGROW>
end

Lr = collectRegistryLong(result, reg, unitLevel, dataView);
L = [L; Lr]; %#ok<AGROW>

if isempty(L)
    error('CorrelationStats:NoRows', 'No metric rows were available for unitLevel=%s.', unitLevel);
end

Tlong = vertcat(L{:});
if ~isempty(combineEpochNames)
    if strlength(epochName) > 0 && epochName ~= combinedEpochLabel
        error('CorrelationStats:ConflictingEpochOptions', ...
            'epochName=%s conflicts with combineEpochNames. Use epochName=''%s'' or leave epochName empty.', ...
            char(epochName), char(combinedEpochLabel));
    end
    Tlong = combineLongAcrossEpochs(Tlong, combineEpochNames, combinedEpochLabel, combineEpochMethod, 'value');
    epochName = combinedEpochLabel;
elseif strlength(epochName) > 0 && ismember('epoch_name', Tlong.Properties.VariableNames)
    Tlong = Tlong(string(Tlong.epoch_name) == epochName, :);
end
Tlong = enforceRequestedUnitLevel(Tlong, unitLevel);
if isempty(Tlong)
    error('CorrelationStats:EmptyAfterFilter', 'No rows remain after filtering by epochName=%s.', char(epochName));
end

allowed = resolveAllowedMetricIDs(metricFamilies, unitLevel, reg, metaTables, selectionMode, p.Results.includeOverlappingDerived);
if ~isempty(allowed)
    Tlong = Tlong(ismember(string(Tlong.metric_id), string(allowed)), :);
end
if isempty(Tlong)
    error('CorrelationStats:NoMetricRows', 'No metrics matched the requested metricFamilies for unitLevel=%s.', unitLevel);
end
Tlong = filterPreferredRegistryMetricRows(Tlong);

[Twide, metricNames, metricIDs] = longToWide(Tlong);
denLong = collectDenominatorLong(result, reg, unitLevel, dataView, epochName);
if ~isempty(combineEpochNames)
    denLong = combineLongAcrossEpochs(denLong, combineEpochNames, combinedEpochLabel, combineEpochMethod, 'sum');
end
[Dwide, denomMetricNames, denomMetricIDs] = longToWideDenominator(denLong);
if isempty(metricNames)
    error('CorrelationStats:NoWideMetrics', 'No metric columns were available after pivoting.');
end

X = NaN(height(Twide), numel(metricNames));
for i = 1:numel(metricNames)
    xi = Twide.(metricNames{i});
    if ~isreal(xi)
        xr = real(double(xi));
        xr(abs(imag(double(xi))) > 0) = NaN;
        xi = xr;
    end
    X(:,i) = double(xi);
end

keepMetric = true(1, numel(metricNames));
for i = 1:numel(metricNames)
    xi = X(:,i);
    xi = xi(isfinite(xi));
    if isempty(xi)
        keepMetric(i) = false;
    elseif numel(unique(xi)) < 2
        keepMetric(i) = false;
    end
end
if ~all(keepMetric)
    dropped = string(metricIDs(~keepMetric));
    warning('CorrelationStats:DroppedDegenerateMetrics', ...
        'Dropping %d metrics with non-real, empty, or zero-variance values in the selected view: %s', ...
        nnz(~keepMetric), strjoin(cellstr(dropped), ', '));
    X = X(:, keepMetric);
    metricNames = metricNames(keepMetric);
    metricIDs = metricIDs(keepMetric);
end
if numel(metricNames) < 2
    error('CorrelationStats:TooFewMetrics', 'Fewer than 2 usable metrics remain after filtering non-real/empty values.');
end

[R, P] = corr(X, 'Rows', 'pairwise', 'Type', method);
N = computePairwiseCounts(X);
R(N < p.Results.minN) = NaN;
P(N < p.Results.minN) = NaN;

pairTable = buildPairTable(metricIDs, R, P, N, p.Results.minN);
if ~isempty(pairTable)
    pairTable.q_value = bhAdjust(pairTable.p_value);
    pairTable = sortrows(pairTable, {'q_value','p_value','abs_rho'}, {'ascend','ascend','descend'});
end

out = struct();
out.source = 'CorrelationStats';
out.unitLevel = unitLevel;
out.metricFamilies = metricFamilies;
out.dataView = dataView;
out.epochName = char(epochName);
out.combineEpochNames = cellstr(combineEpochNames(:))';
out.combinedEpochLabel = char(combinedEpochLabel);
out.combineEpochMethod = combineEpochMethod;
out.method = method;
out.metricSelection = selectionMode;
out.minN = p.Results.minN;
out.longTable = Tlong;
out.wideTable = Twide;
out.denomWideTable = Dwide;
out.denomMetricNames = denomMetricNames(:)';
out.denomMetricIDs = cellstr(denomMetricIDs(:))';
out.metricNames = metricNames(:)';
out.metricIDs = cellstr(metricIDs(:))';
out.R = R;
out.P = P;
out.N = N;
out.pairTable = pairTable;
out.allowedMetricIDs = allowed(:)';
unitLevelsByMetric = strings(size(metricIDs));
for i = 1:numel(metricIDs)
    hit = find(string(Tlong.metric_id) == string(metricIDs(i)), 1, 'first');
    if ~isempty(hit) && ismember('unit_level', Tlong.Properties.VariableNames)
        unitLevelsByMetric(i) = string(Tlong.unit_level(hit));
    end
end
out.normalizationProvenance = bwaResolveNormalizationProvenance(result, dataView, string(metricIDs(:)), unitLevelsByMetric(:));
end

function Tout = combineLongAcrossEpochs(Tin, epochNames, combinedLabel, methodName, valueVar)
Tout = Tin;
if isempty(Tin) || ~ismember('epoch_name', Tin.Properties.VariableNames)
    return
end
mask = ismember(string(Tin.epoch_name), string(epochNames(:)));
Tin = Tin(mask, :);
if isempty(Tin)
    Tout = Tin;
    return
end
groupCols = {'unit_level','entity_id','animal_id','slice_num','file_id','metric_id'};
groupCols = groupCols(ismember(groupCols, Tin.Properties.VariableNames));
if ~ismember(valueVar, Tin.Properties.VariableNames)
    error('CorrelationStats:MissingCombineValue', ...
        'Cannot combine epochs because value column %s is missing.', valueVar);
end
[G, Gtbl] = localGroupRows(Tin, groupCols);
switch lower(methodName)
    case 'mean'
        vals = accumarray(G, double(Tin.(valueVar)), [], @meanOmitNaN, NaN);
    otherwise
        error('CorrelationStats:UnknownCombineMethod', ...
            'Unsupported combineEpochMethod: %s', methodName);
end
if strcmp(valueVar, 'denominator_value')
    vals = accumarray(G, double(Tin.(valueVar)), [], @sum, NaN);
end
Tout = Gtbl;
Tout.epoch_idx = repmat(NaN, height(Tout), 1);
Tout.epoch_name = repmat(string(combinedLabel), height(Tout), 1);
Tout.(valueVar) = vals;
end

function families = normalizeMetricFamilies(families)
families = normalizeStringListLocal(families);
out = {};
for i = 1:numel(families)
    g = lower(char(string(families{i})));
    if strcmp(g, 'all')
        out = {'all'};
        families = out;
        return
    else
        out{end+1} = g; %#ok<AGROW>
    end
end
families = unique(out, 'stable');
end

function out = normalizeStringListLocal(in)
if iscell(in)
    out = cellstr(string(in(:)))';
elseif isstring(in)
    out = cellstr(in(:))';
else
    out = cellstr(string(in));
end
out = out(~cellfun(@isempty, out));
end

function allowed = resolveAllowedMetricIDs(families, unitLevel, reg, metaTables, selectionMode, includeOverlappingDerived)
P = correlationMetricPolicy(unitLevel, selectionMode);
if any(strcmp(families, 'all'))
    families = P.defaultFamilies;
end
allowed = strings(0,1);
for i = 1:numel(families)
    fam = string(families{i});
    ids = selectMetricIDsByFamily(fam, unitLevel, reg, metaTables);
    if fam == "recruitment" && any(strcmp(unitLevel, {'slice','animal'}))
        ids = [ids; dynamicSliceRecruitmentMetricIDs()]; %#ok<AGROW>
    end
    if isempty(ids)
        if ~ismember(fam, supportedMetricFamilies())
            warning('CorrelationStats:UnknownMetricFamily', 'Ignoring unknown metric family: %s', fam);
        end
    else
        allowed = [allowed; string(ids(:))]; %#ok<AGROW>
    end
end
allowed = unique(allowed, 'stable');
allowed = applyCorrelationMetricSelection(allowed, unitLevel, selectionMode);
if ~includeOverlappingDerived
    % Derived collector already filters overlapping metrics, so no extra action needed.
end
end

function families = supportedFamiliesForUnit(unitLevel)
P = correlationMetricPolicy(unitLevel, 'default');
families = P.defaultFamilies;
end

function allowed = applyCorrelationMetricSelection(allowed, unitLevel, selectionMode)
allowed = string(allowed(:));
P = correlationMetricPolicy(unitLevel, selectionMode);
if ~isempty(P.excludeIDs)
    allowed = allowed(~ismember(allowed, string(P.excludeIDs(:))));
end
allowed = unique(allowed, 'stable');
end

function ids = selectMetricIDsByFamily(familyName, unitLevel, reg, metaTables)
ids = strings(0,1);
if ~isempty(reg) && all(ismember({'metric_id','unit_level'}, reg.Properties.VariableNames))
    metricIDs = string(reg.metric_id);
    levels = lower(string(reg.unit_level));
    for i = 1:height(reg)
        id = metricIDs(i);
        lvl = levels(i);
        switch unitLevel
            case 'event'
                continue
            case 'cell'
                if lvl ~= "cell"
                    continue
                end
            case {'slice','animal'}
                if ~(lvl == "cell" || lvl == "slice")
                    continue
                end
            otherwise
                continue
        end
        if bwaResolveMetricFamily(id, reg) == familyName
            if familyName == "synchrony" && lvl ~= "slice"
                continue
            end
            ids(end+1,1) = id; %#ok<AGROW>
        end
    end
end

if ismember(unitLevel, {'event','cell','slice','animal'})
    ids = [ids; selectMetaTableMetricIDs(metaTables.canonical, familyName)]; %#ok<AGROW>
end
if strcmp(unitLevel, 'event')
    ids = [ids; selectMetaTableMetricIDs(metaTables.event, familyName)]; %#ok<AGROW>
end
ids = unique(ids, 'stable');
end

function ids = selectMetaTableMetricIDs(metaStruct, familyName)
ids = strings(0,1);
if ~isstruct(metaStruct) || isempty(fieldnames(metaStruct))
    return
end
fns = string(fieldnames(metaStruct));
for i = 1:numel(fns)
    if bwaResolveMetricFamily(fns(i)) == familyName
        ids(end+1,1) = fns(i); %#ok<AGROW>
    end
end
end

function families = supportedMetricFamilies()
families = ["activity","shape","kinetics","intensity","propagation","containment","recruitment","synchrony","soma"];
end

function [rows, meta] = collectCanonicalLong(result, unitLevel, dataView)
rows = cell(0,1);
meta = struct();
canonical = result.canonical;
if strcmp(dataView, 'norm') && isfield(result, 'canonical_norm') && isstruct(result.canonical_norm)
    canonical = result.canonical_norm;
end
if strcmp(unitLevel, 'cell')
    freqStruct = canonical.frequency;
    if strcmp(dataView, 'recommended')
        freqStruct = selectCanonicalFrequencyForCorrelation(result, 'frequency');
    end
    if isfield(freqStruct, 'ratesByCellEpoch') && ...
            istable(freqStruct.ratesByCellEpoch)
        Tf = freqStruct.ratesByCellEpoch;
        vcol = 'rate_events_per_min';
        if ismember(vcol, Tf.Properties.VariableNames)
            rows{end+1,1} = makeLongFromTable(Tf, 'frequency', vcol, 'cell', {'cell_uid'}); %#ok<AGROW>
            meta.frequency = true;
        end
    end
    if isfield(canonical, 'metrics') && isstruct(canonical.metrics)
        mNames = fieldnames(result.canonical.metrics);
        mNames = filterPreferredCanonicalMetricIDs(mNames);
        for i = 1:numel(mNames)
            m = mNames{i};
            mOut = canonical.metrics.(m);
            if strcmp(dataView, 'recommended')
                mOut = selectCanonicalMetricForCorrelation(result, m);
            end
            if isfield(mOut, 'eventValues') && ~isempty(mOut.eventValues)
                rows{end+1,1} = metricOutEventValuesToLong(mOut, m, unitLevel); %#ok<AGROW>
                meta.(m) = true;
            end
        end
    end
    rows = rows(~cellfun(@isempty, rows));
    return
end
if strcmp(unitLevel, 'event')
    if isfield(canonical, 'metrics') && isstruct(canonical.metrics)
        mNames = fieldnames(result.canonical.metrics);
        mNames = filterPreferredCanonicalMetricIDs(mNames);
        for i = 1:numel(mNames)
            m = mNames{i};
            mOut = canonical.metrics.(m);
            if strcmp(dataView, 'recommended')
                mOut = selectCanonicalMetricForCorrelation(result, m);
            end
            if isfield(mOut, 'eventValues') && ~isempty(mOut.eventValues)
                rows{end+1,1} = metricOutEventValuesToLong(mOut, m, unitLevel); %#ok<AGROW>
                meta.(m) = true;
            end
        end
    end
    rows = rows(~cellfun(@isempty, rows));
    return
end
    if isfield(canonical, 'frequency') && isstruct(canonical.frequency) && isfield(canonical.frequency, 'summaryByGroupEpoch')
        freqSummary = canonical.frequency.summaryByGroupEpoch;
        if strcmp(dataView, 'recommended')
            freqStruct = selectCanonicalFrequencyForCorrelation(result, 'frequency');
            if isfield(freqStruct, 'summaryByGroupEpoch')
                freqSummary = freqStruct.summaryByGroupEpoch;
            end
        end
        rows{end+1,1} = canonicalSummaryToLong(freqSummary, 'frequency', 'mean_rate_events_per_min', unitLevel); %#ok<AGROW>
        meta.frequency = true;
    end
    if isfield(canonical, 'metrics') && isstruct(canonical.metrics)
        mNames = fieldnames(result.canonical.metrics);
        mNames = filterPreferredCanonicalMetricIDs(mNames);
        for i = 1:numel(mNames)
            m = mNames{i};
            meta.(m) = true;
            mOut = canonical.metrics.(m);
            if strcmp(dataView, 'recommended')
                mOut = selectCanonicalMetricForCorrelation(result, m);
            end
            if strcmp(unitLevel, 'cell') && isfield(mOut, 'eventValues') && ~isempty(mOut.eventValues)
                rows{end+1,1} = metricOutEventValuesToLong(mOut, m, unitLevel); %#ok<AGROW>
            else
                S = mOut.summaryByGroupEpoch;
                vcol = inferSummaryValueColumn(S, m);
                rows{end+1,1} = canonicalSummaryToLong(S, m, vcol, unitLevel); %#ok<AGROW>
            end
        end
    end
    rows = rows(~cellfun(@isempty, rows));
    return
end

function ids = filterPreferredCanonicalMetricIDs(ids)
ids = cellstr(string(ids(:)));
if isempty(ids)
    return
end
skip = strings(0,1);
pairs = { ...
    'area_px2','area_um2'; ...
    'duration_frames','duration_sec'; ...
    'rise19_frames','rise19_sec'; ...
    'fall91_frames','fall91_sec'; ...
    'decay_tau_frames','decay_tau_sec'; ...
    'dff_auc','dff_auc_sec'; ...
    'perimeter_px','perimeter_um'; ...
    'propagation_distance_px','propagation_distance_um'; ...
    'max_prop_speed','max_prop_speed_um_per_sec'};
for i = 1:size(pairs,1)
    raw = pairs{i,1};
    cal = pairs{i,2};
    if any(strcmp(ids, raw)) && any(strcmp(ids, cal))
        skip(end+1,1) = string(raw); %#ok<AGROW>
    end
end
ids = ids(~ismember(string(ids), skip));
end

function [rows, meta] = collectDerivedLong(result, unitLevel, dataView, includeOverlappingDerived)
rows = cell(0,1);
meta = struct();
if ~isfield(result, 'derived') || ~isstruct(result.derived) || ~isfield(result.derived, 'standard') || ~isstruct(result.derived.standard)
    return
end
filesTbl = table();
canonical = result.canonical;
if strcmp(dataView, 'norm') && isfield(result, 'canonical_norm') && isstruct(result.canonical_norm)
    canonical = result.canonical_norm;
end
if isfield(canonical, 'subset') && isfield(canonical.subset, 'files')
    filesTbl = canonical.subset.files;
end
epochNames = {};
if isfield(result.derived, 'context') && isstruct(result.derived.context) && ...
        isfield(result.derived.context, 'timeFrameRangeEpochs') && isfield(result.derived.context.timeFrameRangeEpochs, 'epochNames')
    epochNames = result.derived.context.timeFrameRangeEpochs.epochNames;
end
normCfg = struct('mode','none','transform','log2_fold','baselineEpochIdx',1);
if strcmp(dataView, 'norm') && isfield(result, 'normalization') && isstruct(result.normalization)
    if isfield(result.normalization, 'mode'), normCfg.mode = char(string(result.normalization.mode)); end
    if isfield(result.normalization, 'transform'), normCfg.transform = char(string(result.normalization.transform)); end
    if isfield(result.normalization, 'baselineEpochIdx') && ~isempty(result.normalization.baselineEpochIdx)
        normCfg.baselineEpochIdx = double(result.normalization.baselineEpochIdx);
    end
end
entries = bwaCollectDerivedSummaryEntries(result.derived.standard, includeOverlappingDerived);
for i = 1:height(entries)
    key = char(entries.key(i));
    comp = entries.comp{i};
    try
        mOut = bwaCompToMetricOut(comp, filesTbl, epochNames, char(entries.title(i)), normCfg);
        if any(strcmp(unitLevel, {'cell','event'})) && isfield(mOut, 'eventValues') && ~isempty(mOut.eventValues)
            rows{end+1,1} = metricOutEventValuesToLong(mOut, key, unitLevel); %#ok<AGROW>
        else
            S = mOut.summaryByGroupEpoch;
            vcol = inferSummaryValueColumn(S, char(entries.title(i)));
            rows{end+1,1} = canonicalSummaryToLong(S, key, vcol, unitLevel); %#ok<AGROW>
        end
        meta.(key) = true;
    catch
    end
end
rows = rows(~cellfun(@isempty, rows));
end

function Tlong = metricOutEventValuesToLong(metricOut, metricID, unitLevel)
Tlong = table();
if ~any(strcmp(unitLevel, {'cell','event'})) || ~isfield(metricOut, 'eventValues') || isempty(metricOut.eventValues)
    return
end
T = metricOut.eventValues;
if ~istable(T) || isempty(T)
    return
end
valueCol = inferEventValueColumn(T, metricOut);
if ~ismember(valueCol, T.Properties.VariableNames)
    return
end
switch unitLevel
    case 'cell'
        if ~ismember('primary_cell_uid', T.Properties.VariableNames)
            return
        end
        Tlong = aggregateTableToLevel(T, metricID, valueCol, 'cell', cellGroupColsForTable(T));
    case 'event'
        if ~ismember('event_uid', T.Properties.VariableNames)
            return
        end
        Tlong = makeLongFromTable(T, metricID, valueCol, 'event', {'event_uid'});
end
end

function valueCol = inferEventValueColumn(T, metricOut)
valueCol = "";
if isfield(metricOut, 'metricCol') && strlength(string(metricOut.metricCol)) > 0
    cand = char(string(metricOut.metricCol));
    if ismember(cand, T.Properties.VariableNames)
        valueCol = cand;
        return
    end
end
if ismember('value', T.Properties.VariableNames)
    valueCol = 'value';
    return
end
cand = setdiff(T.Properties.VariableNames, {'event_uid','file_id','epoch_idx','epoch_name','primary_cell_uid','animal_id','slice_num'});
numCand = cand(varfun(@(x) isnumeric(x) || islogical(x), T(:, cand), 'OutputFormat', 'uniform'));
if ~isempty(numCand)
    valueCol = numCand{1};
end
end

function Tlong = canonicalSummaryToLong(S, metricID, valueCol, unitLevel)
Tlong = table();
if ~istable(S) || isempty(S) || ~ismember(valueCol, S.Properties.VariableNames)
    return
end
switch unitLevel
    case 'slice'
        req = {'animal_id','slice_num','file_id','epoch_idx','epoch_name',valueCol};
        req = req(ismember(req, S.Properties.VariableNames));
        T = S(:, req);
        T.entity_id = string(T.file_id);
        T.unit_level = repmat("slice", height(T), 1);
    case 'animal'
        if ~ismember('animal_id', S.Properties.VariableNames)
            return
        end
        keep = {'animal_id','epoch_idx','epoch_name',valueCol};
        keep = keep(ismember(keep, S.Properties.VariableNames));
        T0 = S(:, keep);
        [G, Gtbl] = localGroupRows(T0, {'animal_id','epoch_idx','epoch_name'});
        vals = accumarray(G, double(T0.(valueCol)), [], @meanOmitNaN, NaN);
        T = Gtbl;
        T.(valueCol) = vals;
        T.entity_id = string(T.animal_id);
        T.unit_level = repmat("animal", height(T), 1);
        T.slice_num = repmat("", height(T), 1);
        T.file_id = repmat("", height(T), 1);
    otherwise
        return
end
T.metric_id = repmat(string(metricID), height(T), 1);
T.value = double(T.(valueCol));
Tlong = T(:, {'unit_level','entity_id','animal_id','slice_num','file_id','epoch_idx','epoch_name','metric_id','value'});
end

function rows = collectRegistryLong(result, reg, unitLevel, dataView)
rows = cell(0,1);
if isempty(reg) || ~all(ismember({'metric_id','unit_level','source_table','value_column_raw','value_column_norm'}, reg.Properties.VariableNames)) || ...
        ~isfield(result, 'metrics') || ~isstruct(result.metrics) || ~isfield(result.metrics, 'tables')
    return
end
for i = 1:height(reg)
    metricID = char(string(reg.metric_id(i)));
    sourceTable = char(string(reg.source_table(i)));
    [regRow, tablesStruct] = selectRegistrySourceForCorrelation(result, reg, i, dataView);
    if ~isfield(tablesStruct, sourceTable)
        continue
    end
    T = tablesStruct.(sourceTable);
    if ~istable(T) || isempty(T)
        continue
    end
    rawCol = char(string(regRow.value_column_raw));
    normCol = char(string(regRow.value_column_norm));
    ycol = rawCol;
    if (strcmp(dataView, 'norm') || strcmp(dataView, 'recommended')) && strlength(string(normCol)) > 0 && ismember(normCol, T.Properties.VariableNames)
        ycol = normCol;
    end
    if ~ismember(ycol, T.Properties.VariableNames)
        continue
    end
    lvl = lower(char(string(reg.unit_level(i))));
    switch unitLevel
        case 'cell'
            if ~strcmp(lvl, 'cell') || ~ismember('primary_cell_uid', T.Properties.VariableNames)
                continue
            end
            rows{end+1,1} = makeLongFromTable(T, metricID, ycol, 'cell', {'primary_cell_uid'}); %#ok<AGROW>
        case 'slice'
            if strcmp(lvl, 'cell')
                rows{end+1,1} = aggregateTableToLevel(T, metricID, ycol, 'slice', {'animal_id','slice_num','file_id','epoch_idx','epoch_name'}); %#ok<AGROW>
            elseif strcmp(lvl, 'slice')
                rows{end+1,1} = makeLongFromTable(T, metricID, ycol, 'slice', {'file_id'}); %#ok<AGROW>
            end
        case 'animal'
            if strcmp(lvl, 'cell')
                rows{end+1,1} = aggregateTableToLevel(T, metricID, ycol, 'animal', {'animal_id','epoch_idx','epoch_name'}); %#ok<AGROW>
            elseif strcmp(lvl, 'slice')
                rows{end+1,1} = aggregateTableToLevel(T, metricID, ycol, 'animal', {'animal_id','epoch_idx','epoch_name'}); %#ok<AGROW>
            end
    end
end

function [regRow, tablesStruct] = selectRegistrySourceForCorrelation(result, reg, idx, dataView)
regRow = reg(idx, :);
tablesStruct = result.metrics.tables;
if ~strcmp(dataView, 'recommended') || ~isfield(result, 'metrics_by_transform') || ~isstruct(result.metrics_by_transform)
    return
end
metricID = char(string(reg.metric_id(idx)));
metricUnitLevel = char(string(reg.unit_level(idx)));
tr = bwaPreferredTransformForMetricID(metricID, metricUnitLevel);
field = matlab.lang.makeValidName(tr);
if ~isfield(result.metrics_by_transform, field)
    return
end
M = result.metrics_by_transform.(field);
if isfield(M, 'tables')
    tablesStruct = M.tables;
end
if isfield(M, 'registry') && istable(M.registry)
    hit = find(string(M.registry.metric_id) == string(metricID), 1, 'first');
    if ~isempty(hit)
        regRow = M.registry(hit, :);
    end
end
end
if any(strcmp(unitLevel, {'slice','animal'}))
    dynRows = collectDynamicSliceRecruitmentLong(result, unitLevel, dataView);
    rows = [rows; dynRows]; %#ok<AGROW>
end
rows = rows(~cellfun(@isempty, rows));
end

function T = enforceRequestedUnitLevel(T, unitLevel)
if isempty(T) || ~ismember('unit_level', T.Properties.VariableNames)
    return
end
want = string(unitLevel);
T = T(string(T.unit_level) == want, :);
end

function Tlong = makeLongFromTable(T, metricID, ycol, levelName, entityCols)
Tlong = table();
req = [entityCols {'animal_id','slice_num','file_id','epoch_idx','epoch_name', ycol}];
req = unique(req, 'stable');
req = req(ismember(req, T.Properties.VariableNames));
if ~all(ismember({'epoch_idx','epoch_name', ycol}, req))
    return
end
T0 = T(:, req);
entityCol = entityCols{1};
T0.entity_id = string(T0.(entityCol));
if ~ismember('animal_id', T0.Properties.VariableNames), T0.animal_id = repmat("", height(T0), 1); end
if ~ismember('slice_num', T0.Properties.VariableNames), T0.slice_num = repmat("", height(T0), 1); end
if ~ismember('file_id', T0.Properties.VariableNames), T0.file_id = repmat("", height(T0), 1); end
T0.unit_level = repmat(string(levelName), height(T0), 1);
T0.metric_id = repmat(string(metricID), height(T0), 1);
T0.value = double(T0.(ycol));
Tlong = T0(:, {'unit_level','entity_id','animal_id','slice_num','file_id','epoch_idx','epoch_name','metric_id','value'});
end

function Tlong = aggregateTableToLevel(T, metricID, ycol, levelName, groupCols)
Tlong = table();
if ~all(ismember([groupCols {ycol}], T.Properties.VariableNames))
    return
end
T0 = T(:, [groupCols {ycol}]);
[G, Gtbl] = localGroupRows(T0, groupCols);
vals = accumarray(G, double(T0.(ycol)), [], @meanOmitNaN, NaN);
T = Gtbl;
T.value = vals;
switch levelName
    case 'event'
        T.entity_id = string(T.event_uid);
        if ~ismember('slice_num', T.Properties.VariableNames), T.slice_num = repmat("", height(T), 1); end
        if ~ismember('file_id', T.Properties.VariableNames), T.file_id = repmat("", height(T), 1); end
    case 'cell'
        T.entity_id = string(T.primary_cell_uid);
        if ~ismember('slice_num', T.Properties.VariableNames), T.slice_num = repmat("", height(T), 1); end
        if ~ismember('file_id', T.Properties.VariableNames), T.file_id = repmat("", height(T), 1); end
    case 'slice'
        T.entity_id = string(T.file_id);
        if ~ismember('slice_num', T.Properties.VariableNames), T.slice_num = repmat("", height(T), 1); end
        if ~ismember('file_id', T.Properties.VariableNames), T.file_id = repmat("", height(T), 1); end
    case 'animal'
        T.entity_id = string(T.animal_id);
        T.slice_num = repmat("", height(T), 1);
        T.file_id = repmat("", height(T), 1);
end
T.unit_level = repmat(string(levelName), height(T), 1);
T.metric_id = repmat(string(metricID), height(T), 1);
if ~ismember('animal_id', T.Properties.VariableNames), T.animal_id = repmat("", height(T), 1); end
Tlong = T(:, {'unit_level','entity_id','animal_id','slice_num','file_id','epoch_idx','epoch_name','metric_id','value'});
end

function cols = cellGroupColsForTable(T)
cols = {'animal_id','slice_num','file_id','primary_cell_uid','epoch_idx','epoch_name'};
cols = cols(ismember(cols, T.Properties.VariableNames));
if ~all(ismember({'primary_cell_uid','epoch_idx','epoch_name'}, cols))
    cols = {};
end
end

function y = meanOmitNaN(x)
x = double(x);
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = mean(x);
end
end

function [rows, meta] = collectEventLocalLong(result, dataView)
rows = cell(0,1);
meta = struct('ids', strings(0,1));
metricContainer = result.metrics;
if strcmp(dataView, 'recommended') && isfield(result, 'metrics_by_transform') && isstruct(result.metrics_by_transform) && ...
        isfield(result.metrics_by_transform, 'delta')
    metricContainer = result.metrics_by_transform.delta;
end
if ~isfield(metricContainer, 'tables') || ~isstruct(metricContainer.tables)
    return
end
tables = metricContainer.tables;
tblNames = fieldnames(tables);
for i = 1:numel(tblNames)
    tname = tblNames{i};
    if ~startsWith(tname, 'event_')
        continue
    end
    T = tables.(tname);
    if ~istable(T) || isempty(T) || ~ismember('event_uid', T.Properties.VariableNames)
        continue
    end
    valueCols = pickEventLocalValueColumns(T, dataView);
    for j = 1:numel(valueCols)
        vcol = valueCols{j};
        metricID = ['event.' stripEventLocalColumnPrefix(vcol)];
        rows{end+1,1} = makeLongFromTable(T, metricID, vcol, 'event', {'event_uid'}); %#ok<AGROW>
        meta.ids(end+1,1) = string(metricID); %#ok<AGROW>
    end
end
rows = rows(~cellfun(@isempty, rows));
meta.ids = unique(meta.ids, 'stable');
end

function metricOut = selectCanonicalMetricForCorrelation(result, metricID)
metricOut = result.canonical.metrics.(metricID);
tr = bwaPreferredTransformForMetricID(metricID, 'canonical');
field = matlab.lang.makeValidName(tr);
if isfield(result, 'canonical_norms') && isstruct(result.canonical_norms) && ...
        isfield(result.canonical_norms, field) && isfield(result.canonical_norms.(field), 'metrics') && ...
        isfield(result.canonical_norms.(field).metrics, metricID)
    metricOut = result.canonical_norms.(field).metrics.(metricID);
elseif isfield(result, 'canonical_norm') && isfield(result.canonical_norm, 'metrics') && isfield(result.canonical_norm.metrics, metricID)
    metricOut = result.canonical_norm.metrics.(metricID);
end
end

function freqOut = selectCanonicalFrequencyForCorrelation(result, metricID)
freqOut = result.canonical.frequency;
tr = bwaPreferredTransformForMetricID(metricID, 'canonical');
field = matlab.lang.makeValidName(tr);
if isfield(result, 'canonical_norms') && isstruct(result.canonical_norms) && ...
        isfield(result.canonical_norms, field) && isfield(result.canonical_norms.(field), 'frequency')
    freqOut = result.canonical_norms.(field).frequency;
elseif isfield(result, 'canonical_norm') && isfield(result.canonical_norm, 'frequency')
    freqOut = result.canonical_norm.frequency;
end
end

function cols = pickEventLocalValueColumns(T, dataView)
cols = {};
vars = string(T.Properties.VariableNames);
skip = ["event_uid","animal_id","slice_num","file_id","primary_cell_uid","partition_cell_uid","epoch_idx","epoch_name","cell_uid", ...
    "frame_rate_s","spatial_res","n_frames","is_primary_cell"];
preferredPairs = [ ...
    "soma_centroid_dist_px","soma_centroid_dist_um"; ...
    "soma_start_dist_px","soma_start_dist_um"; ...
    "soma_end_dist_px","soma_end_dist_um"; ...
    "soma_net_distance_change_px","soma_net_distance_change_um"; ...
    "perimeter_px","perimeter_um"; ...
    "propagation_distance_px","propagation_distance_um"; ...
    "outside_cell_area_npix","outside_cell_area_um2"; ...
    "cell_entry_delay_from_primary_frames","cell_entry_delay_from_primary_sec"; ...
    "cell_overlap_duration_frames","cell_overlap_duration_sec"];
complementBase = ["cell_overlap_frac_of_event", "frac_outside_cell"];
redundantEventBase = ["soma_start_dist_um","soma_end_dist_um", ...
    "soma_start_dist_px","soma_end_dist_px"];
for i = 1:numel(vars)
    vn = vars(i);
    if any(vn == skip)
        continue
    end
    if contains(vn, "threshold")
        continue
    end
    base = regexprep(char(vn), '_norm$', '');
    base = string(base);
    if any(base == redundantEventBase)
        continue
    end
    if strcmp(dataView, 'norm')
        if endsWith(vn, "_norm")
            if ismember(base, complementBase) && any(ismember(vars, complementBase + "_norm"))
                if base == "cell_overlap_frac_of_event"
                    continue
                end
            end
            if any(preferredPairs(:,1) == base & ismember(preferredPairs(:,2) + "_norm", vars))
                continue
            end
            cols{end+1} = char(vn); %#ok<AGROW>
            continue
        end
        if any(vars == vn + "_norm")
            continue
        end
        continue
    end
    if ismember(vn, complementBase) && any(ismember(complementBase, vars))
        if vn == "cell_overlap_frac_of_event"
            continue
        end
    end
    if (endsWith(vn, "_px") || endsWith(vn, "_frames") || endsWith(vn, "_npix")) && any(preferredPairs(:,1) == vn & ismember(preferredPairs(:,2), vars))
        continue
    end
    x = T.(char(vn));
    if isnumeric(x) || islogical(x)
        cols{end+1} = char(vn); %#ok<AGROW>
    end
end
cols = unique(cols, 'stable');
end

function out = stripEventLocalColumnPrefix(vcol)
out = char(string(vcol));
out = regexprep(out, '_norm$', '');
end

function [Twide, metricNames, metricIDs] = longToWide(Tlong)
Tlong = Tlong(isfinite(Tlong.value), :);
if isempty(Tlong)
    Twide = table();
    metricNames = {};
    metricIDs = strings(0,1);
    return
end
keyCols = {'unit_level','entity_id','epoch_name'};
[Gu, Ukey] = localGroupRows(Tlong, keyCols);
firstIdx = accumarray(Gu, (1:height(Tlong))', [], @min);
U = Tlong(firstIdx, {'unit_level','entity_id','animal_id','slice_num','file_id','epoch_idx','epoch_name'});
metricIDs = unique(string(Tlong.metric_id), 'stable');
metricNames = matlab.lang.makeValidName(cellstr(metricIDs));
X = NaN(height(U), numel(metricIDs));
for j = 1:numel(metricIDs)
    idx = string(Tlong.metric_id) == metricIDs(j);
    Gj = Gu(idx);
    vals = double(Tlong.value(idx));
    % duplicate rows for same unit+metric are averaged
    xj = accumarray(Gj, vals, [height(U),1], @mean, NaN);
    X(:,j) = xj;
end
Twide = U;
for j = 1:numel(metricIDs)
    Twide.(metricNames{j}) = X(:,j);
end
end

function [G, grpTbl] = localGroupRows(T, cols)
n = height(T);
if n == 0
    G = zeros(0,1);
    grpTbl = T(:, cols);
    return
end
K = strings(n,1);
for i = 1:numel(cols)
    K = K + "|" + toKeyLocal(T.(cols{i}));
end
[~,~,G] = unique(K, 'stable');
G = full(double(G(:)));
firstIdx = accumarray(G, (1:n)', [], @min);
grpTbl = T(firstIdx, cols);
end

function s = toKeyLocal(v)
if isstring(v)
    s = v(:);
elseif iscell(v)
    s = strings(numel(v),1);
    for i = 1:numel(v)
        s(i) = string(v{i});
    end
elseif isnumeric(v) || islogical(v)
    s = string(v(:));
else
    s = string(v(:));
end
end

function valueCol = inferSummaryValueColumn(S, metricName)
valueCol = ['mean_' metricName];
if ismember(valueCol, S.Properties.VariableNames)
    return
end
cand = S.Properties.VariableNames(startsWith(S.Properties.VariableNames, 'mean_'));
if isempty(cand)
    error('CorrelationStats:NoSummaryValueColumn', 'No mean_* column found for metric %s.', metricName);
end
valueCol = cand{1};
end

function N = computePairwiseCounts(X)
p = size(X,2);
N = zeros(p,p);
for i = 1:p
    xi = isfinite(X(:,i));
    for j = i:p
        n = nnz(xi & isfinite(X(:,j)));
        N(i,j) = n;
        N(j,i) = n;
    end
end
end

function T = buildPairTable(metricIDs, R, P, N, minN)
metricIDs = string(metricIDs(:));
rows = [];
ma = strings(0,1); mb = strings(0,1); rho = []; pval = []; nobs = [];
for i = 1:numel(metricIDs)
    for j = i+1:numel(metricIDs)
        if ~isfinite(R(i,j)) || ~isfinite(P(i,j)) || N(i,j) < minN
            continue
        end
        ma(end+1,1) = metricIDs(i); %#ok<AGROW>
        mb(end+1,1) = metricIDs(j); %#ok<AGROW>
        rho(end+1,1) = R(i,j); %#ok<AGROW>
        pval(end+1,1) = P(i,j); %#ok<AGROW>
        nobs(end+1,1) = N(i,j); %#ok<AGROW>
    end
end
if isempty(ma)
    T = table();
    return
end
T = table(ma, mb, rho, abs(rho), pval, nobs, 'VariableNames', ...
    {'metric_a','metric_b','rho','abs_rho','p_value','n_obs'});
end

function q = bhAdjust(p)
p = double(p(:));
q = NaN(size(p));
keep = isfinite(p);
if ~any(keep)
    return
end
pv = p(keep);
[m, ord] = sort(pv);
n = numel(m);
qq = m .* n ./ (1:n)';
qq = flipud(cummin(flipud(qq)));
qq = min(qq, 1);
qv = NaN(size(pv));
qv(ord) = qq;
q(keep) = qv;
end

function T = filterPreferredRegistryMetricRows(T)
if isempty(T) || ~ismember('metric_id', T.Properties.VariableNames)
    return
end
ids = string(T.metric_id);
pairs = [ ...
    "cell.mean_soma_centroid_dist_px", "cell.mean_soma_centroid_dist_um"; ...
    "cell.mean_soma_start_dist_px", "cell.mean_soma_start_dist_um"; ...
    "cell.mean_soma_end_dist_px", "cell.mean_soma_end_dist_um"; ...
    "cell.mean_perimeter_px", "cell.mean_perimeter_um"; ...
    "cell.mean_propagation_distance_px", "cell.mean_propagation_distance_um"; ...
    "cell.mean_neighbor_entry_delay_frames", "cell.mean_neighbor_entry_delay_sec"];
for i = 1:size(pairs,1)
    rawID = pairs(i,1);
    calID = pairs(i,2);
    if any(ids == calID)
        T = T(ids ~= rawID, :);
        ids = string(T.metric_id);
    end
end
end

function rows = collectDenominatorLong(result, reg, unitLevel, dataView, epochName)
rows = table();
if ~isfield(result, 'metrics') || ~isstruct(result.metrics) || ~isfield(result.metrics, 'tables')
    return
end
parts = {};
if ~isempty(reg)
    for i = 1:height(reg)
        metricID = char(string(reg.metric_id(i)));
        sourceTable = char(string(reg.source_table(i)));
        if ~isfield(result.metrics.tables, sourceTable)
            continue
        end
        denomCol = denominatorColumnForMetric(metricID, sourceTable);
        if strlength(string(denomCol)) == 0
            continue
        end
        T = result.metrics.tables.(sourceTable);
        if ~istable(T) || isempty(T) || ~ismember(denomCol, T.Properties.VariableNames)
            continue
        end
        lvl = lower(char(string(reg.unit_level(i))));
        switch unitLevel
            case 'cell'
                if ~strcmp(lvl, 'cell') || ~ismember('primary_cell_uid', T.Properties.VariableNames)
                    continue
                end
                parts{end+1,1} = makeDenLongFromTable(T, metricID, denomCol, 'cell', {'primary_cell_uid'}); %#ok<AGROW>
            case 'slice'
                if strcmp(lvl, 'cell')
                    parts{end+1,1} = aggregateDenTableToLevel(T, metricID, denomCol, 'slice', {'animal_id','slice_num','file_id','epoch_idx','epoch_name'}); %#ok<AGROW>
                elseif strcmp(lvl, 'slice')
                    parts{end+1,1} = makeDenLongFromTable(T, metricID, denomCol, 'slice', {'file_id'}); %#ok<AGROW>
                end
            case 'animal'
                if any(strcmp(lvl, {'cell','slice'}))
                    parts{end+1,1} = aggregateDenTableToLevel(T, metricID, denomCol, 'animal', {'animal_id','epoch_idx','epoch_name'}); %#ok<AGROW>
                end
        end
    end
end
if any(strcmp(unitLevel, {'slice','animal'}))
    dyn = collectDynamicSliceRecruitmentDenominatorLong(result, unitLevel, dataView);
    for i = 1:numel(dyn)
        parts{end+1,1} = dyn{i}; %#ok<AGROW>
    end
end
if ~isempty(parts)
    rows = vertcat(parts{:});
    if strlength(epochName) > 0 && ismember('epoch_name', rows.Properties.VariableNames)
        rows = rows(string(rows.epoch_name) == epochName, :);
    end
    rows = filterPreferredRegistryMetricRows(rows);
end
end

function denomCol = denominatorColumnForMetric(metricID, sourceTable)
metricID = string(metricID);
sourceTable = string(sourceTable);
denomCol = "";
switch metricID
    case "cell.frac_multicell_events"
        denomCol = "n_total_events";
    case "cell.frac_propagative_events"
        denomCol = "n_total_events";
    case "cell.frac_events_with_neighbor_recruitment"
        denomCol = "n_primary_events";
    case "slice.frac_primary_events_with_neighbor_recruitment"
        denomCol = "n_primary_events";
end
if sourceTable == "cell_soma" && any(metricID == ["cell.frac_soma_events","cell.frac_toward_events","cell.frac_away_events","cell.frac_traverse_events","cell.frac_net_toward_soma_events","cell.frac_net_away_soma_events","cell.frac_net_stationary_soma_events"])
    denomCol = "n_total_events";
end
if sourceTable == "cell_boundary" && metricID == "cell.frac_escape_events_gt_threshold"
    denomCol = "n_total_events";
end
end

function Tlong = makeDenLongFromTable(T, metricID, dcol, levelName, entityCols)
Tlong = table();
req = [entityCols {'animal_id','slice_num','file_id','epoch_idx','epoch_name', char(dcol)}];
req = unique(req, 'stable');
req = req(ismember(req, T.Properties.VariableNames));
if ~all(ismember({'epoch_idx','epoch_name', char(dcol)}, req))
    return
end
T0 = T(:, req);
entityCol = entityCols{1};
T0.entity_id = string(T0.(entityCol));
if ~ismember('animal_id', T0.Properties.VariableNames), T0.animal_id = repmat("", height(T0), 1); end
if ~ismember('slice_num', T0.Properties.VariableNames), T0.slice_num = repmat("", height(T0), 1); end
if ~ismember('file_id', T0.Properties.VariableNames), T0.file_id = repmat("", height(T0), 1); end
T0.unit_level = repmat(string(levelName), height(T0), 1);
T0.metric_id = repmat(string(metricID), height(T0), 1);
T0.denominator_value = double(T0.(dcol));
Tlong = T0(:, {'unit_level','entity_id','animal_id','slice_num','file_id','epoch_idx','epoch_name','metric_id','denominator_value'});
end

function Tlong = aggregateDenTableToLevel(T, metricID, dcol, levelName, groupCols)
Tlong = table();
if ~all(ismember([groupCols {char(dcol)}], T.Properties.VariableNames))
    return
end
T0 = T(:, [groupCols {char(dcol)}]);
[G, Gtbl] = localGroupRows(T0, groupCols);
vals = accumarray(G, double(T0.(dcol)), [], @mean, NaN);
T = Gtbl;
T.denominator_value = vals;
switch levelName
    case 'slice'
        T.entity_id = string(T.file_id);
        if ~ismember('slice_num', T.Properties.VariableNames), T.slice_num = repmat("", height(T), 1); end
        if ~ismember('file_id', T.Properties.VariableNames), T.file_id = repmat("", height(T), 1); end
    case 'animal'
        T.entity_id = string(T.animal_id);
        T.slice_num = repmat("", height(T), 1);
        T.file_id = repmat("", height(T), 1);
end
T.unit_level = repmat(string(levelName), height(T), 1);
T.metric_id = repmat(string(metricID), height(T), 1);
if ~ismember('animal_id', T.Properties.VariableNames), T.animal_id = repmat("", height(T), 1); end
Tlong = T(:, {'unit_level','entity_id','animal_id','slice_num','file_id','epoch_idx','epoch_name','metric_id','denominator_value'});
end

function ids = dynamicSliceRecruitmentMetricIDs()
ids = string({ ...
    'slice.frac_primary_events_with_neighbor_recruitment'; ...
    'slice.mean_neighbor_entry_delay_sec'; ...
    'slice.mean_neighbor_peak_overlap_frac_of_cell'; ...
    'slice.mean_neighbor_partitions_per_event'});
end

function rows = collectDynamicSliceRecruitmentLong(result, unitLevel, dataView)
rows = cell(0,1);
transform = resolveDynamicTransform(result, dataView, 'slice.frac_primary_events_with_neighbor_recruitment', 'slice');
T = bwaBuildDynamicSliceRecruitmentTable(result, 'transform', transform);
if ~istable(T) || isempty(T)
    return
end

spec = { ...
    'slice.frac_primary_events_with_neighbor_recruitment','frac_primary_events_with_neighbor_recruitment','frac_primary_events_with_neighbor_recruitment_norm'; ...
    'slice.mean_neighbor_entry_delay_sec','mean_neighbor_entry_delay_sec','mean_neighbor_entry_delay_sec_norm'; ...
    'slice.mean_neighbor_peak_overlap_frac_of_cell','mean_neighbor_peak_overlap_frac_of_cell','mean_neighbor_peak_overlap_frac_of_cell_norm'; ...
    'slice.mean_neighbor_partitions_per_event','mean_neighbor_partitions_per_event','mean_neighbor_partitions_per_event_norm'};
for i = 1:size(spec,1)
    metricID = spec{i,1};
    rawCol = spec{i,2};
    normCol = spec{i,3};
    ycol = rawCol;
    if (strcmp(dataView, 'norm') || strcmp(dataView, 'recommended')) && ismember(normCol, T.Properties.VariableNames)
        ycol = normCol;
    end
    switch unitLevel
        case 'slice'
            rows{end+1,1} = makeLongFromTable(T, metricID, ycol, 'slice', {'file_id'}); %#ok<AGROW>
        case 'animal'
            rows{end+1,1} = aggregateTableToLevel(T, metricID, ycol, 'animal', {'animal_id','epoch_idx','epoch_name'}); %#ok<AGROW>
    end
end
rows = rows(~cellfun(@isempty, rows));
end

function rows = collectDynamicSliceRecruitmentDenominatorLong(result, unitLevel, dataView)
rows = cell(0,1);
transform = resolveDynamicTransform(result, dataView, 'slice.frac_primary_events_with_neighbor_recruitment', 'slice');
T = bwaBuildDynamicSliceRecruitmentTable(result, 'transform', transform);
if ~istable(T) || isempty(T) || ~ismember('n_primary_events', T.Properties.VariableNames)
    return
end
switch unitLevel
    case 'slice'
        rows{end+1,1} = makeDenLongFromTable(T, 'slice.frac_primary_events_with_neighbor_recruitment', 'n_primary_events', 'slice', {'file_id'}); %#ok<AGROW>
    case 'animal'
        rows{end+1,1} = aggregateDenTableToLevel(T, 'slice.frac_primary_events_with_neighbor_recruitment', 'n_primary_events', 'animal', {'animal_id','epoch_idx','epoch_name'}); %#ok<AGROW>
end
rows = rows(~cellfun(@isempty, rows));
end

function transform = resolveDynamicTransform(result, dataView, metricID, unitLevel)
switch lower(char(string(dataView)))
    case 'recommended'
        transform = bwaPreferredTransformForMetricID(metricID, unitLevel);
    case 'norm'
        transform = 'delta';
        if isstruct(result) && isfield(result, 'normalization') && isstruct(result.normalization) && ...
                isfield(result.normalization, 'transform')
            transform = lower(char(string(result.normalization.transform)));
        end
    otherwise
        transform = 'delta';
end
end

function [Twide, metricNames, metricIDs] = longToWideDenominator(Tlong)
Twide = table(); metricNames = {}; metricIDs = strings(0,1);
if isempty(Tlong)
    return
end
Tlong = Tlong(isfinite(Tlong.denominator_value), :);
if isempty(Tlong)
    return
end
keys = {'unit_level','entity_id','epoch_name'};
keys = keys(ismember(keys, Tlong.Properties.VariableNames));
metricIDs = unique(string(Tlong.metric_id), 'stable');
[G, ~] = localGroupRows(Tlong, keys);
firstIdx = accumarray(G, (1:height(Tlong))', [], @min);
Twide = Tlong(firstIdx, {'unit_level','entity_id','animal_id','slice_num','file_id','epoch_idx','epoch_name'});
metricNames = cell(1, numel(metricIDs));
for j = 1:numel(metricIDs)
    idx = string(Tlong.metric_id) == metricIDs(j);
    Tj = Tlong(idx, {'unit_level','entity_id','epoch_name','denominator_value'});
    if j == 1
        vname = matlab.lang.makeValidName(char(metricIDs(j)));
    else
        vname = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(char(metricIDs(j))), metricNames(1:j-1));
    end
    metricNames{j} = vname;
    Tj.Properties.VariableNames{end} = vname;
    Twide = outerjoin(Twide, Tj, 'Keys', {'unit_level','entity_id','epoch_name'}, 'MergeKeys', true, 'Type', 'left');
end
end
