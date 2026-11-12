function [normOut, meta] = normalizeCanonicalOutput(canonicalOut, varargin)
% Normalize canonical summary outputs using baseline-referenced transforms.
%
% Options:
%   mode: none | by_cell_baseline_mean | by_global_baseline_mean
%   transform: log2_fold | delta | percent_change
%   baselineEpoch: epoch index (numeric) or epoch name (char/string)

p = inputParser;
p.StructExpand = false;
addRequired(p, 'canonicalOut', @isstruct);
addParameter(p, 'mode', 'none', @(x) ischar(x) || isstring(x));
addParameter(p, 'transform', 'log2_fold', @(x) ischar(x) || isstring(x));
addParameter(p, 'baselineEpoch', 1);
parse(p, canonicalOut, varargin{:});

mode = lower(char(string(p.Results.mode)));
transform = lower(char(string(p.Results.transform)));

meta = struct();
meta.enabled = ~strcmp(mode, 'none');
meta.mode = mode;
meta.transform = transform;
meta.baselineEpoch = p.Results.baselineEpoch;
meta.baselineEpochIdx = NaN;
meta.notes = '';

normOut = canonicalOut;
if strcmp(mode, 'none')
    meta.notes = 'No normalization applied.';
    return
end

if ~isfield(canonicalOut, 'frequency') || ~isfield(canonicalOut.frequency, 'epochTable')
    meta.notes = 'Missing canonical frequency/epochTable.';
    return
end
baseIdx = resolveBaselineEpochIdx(canonicalOut.frequency.epochTable, p.Results.baselineEpoch);
meta.baselineEpochIdx = baseIdx;

if isfield(canonicalOut, 'frequency') && isfield(canonicalOut.frequency, 'ratesByCellEpoch')
    normOut.frequency = normalizeFrequency(canonicalOut.frequency, mode, transform, baseIdx);
end

if isfield(canonicalOut, 'metrics') && isstruct(canonicalOut.metrics)
    mNames = fieldnames(canonicalOut.metrics);
    for i = 1:numel(mNames)
        m = mNames{i};
        normOut.metrics.(m) = normalizeMetric(canonicalOut.metrics.(m), mode, transform, baseIdx);
    end
end

meta.notes = 'Normalization applied to canonical frequency and metric summaries.';
end

function out = normalizeFrequency(freqOut, mode, transform, baseIdx)
out = freqOut;
T = freqOut.ratesByCellEpoch;
if isempty(T) || ~ismember('rate_events_per_min', T.Properties.VariableNames)
    return
end
vals = double(T.rate_events_per_min);
base = baselineVector(T, vals, mode, transform, baseIdx, 'cell_uid');
normVals = applyTransform(vals, base, transform);
T.rate_events_per_min = normVals;
out.ratesByCellEpoch = T;

grpCols = {'epoch_idx','epoch_name'};
if isfield(freqOut, 'groupBy') && ~isempty(freqOut.groupBy)
    grpCols = unique([cellstr(string(freqOut.groupBy(:)')) grpCols], 'stable');
end
grpCols = grpCols(ismember(grpCols, T.Properties.VariableNames));
[G, grp] = groupRows(T, grpCols);
v = T.rate_events_per_min;
nGrp = max([G;0]);
meanVals = NaN(nGrp,1); medVals = NaN(nGrp,1); stdVals = NaN(nGrp,1); nVals = zeros(nGrp,1);
for k = 1:nGrp
    x = v(G == k);
    x = x(isfinite(x));
    nVals(k) = numel(x);
    if isempty(x), continue; end
    meanVals(k) = mean(x);
    medVals(k) = median(x);
    stdVals(k) = std(x);
end
S = grp;
S.mean_rate_events_per_min = meanVals;
S.median_rate_events_per_min = medVals;
S.std_rate_events_per_min = stdVals;
S.n_cells = nVals;
S.sem_rate_events_per_min = stdVals ./ sqrt(max(nVals,1));
out.summaryByGroupEpoch = S;
out.normalization = struct('mode', mode, 'transform', transform, 'baselineEpochIdx', baseIdx);
end

function out = normalizeMetric(metricOut, mode, transform, baseIdx)
out = metricOut;
metricCol = metricOut.metricCol;
didEventLevel = false;
if isfield(metricOut, 'eventValues') && ~isempty(metricOut.eventValues)
    T = metricOut.eventValues;
    if ismember(metricCol, T.Properties.VariableNames)
        vals = double(T.(metricCol));
        keyCol = selectCellKeyColumn(T);
        base = baselineVector(T, vals, mode, transform, baseIdx, keyCol);
        normVals = applyTransform(vals, base, transform);
        T.(metricCol) = normVals;
        out.eventValues = T;
        didEventLevel = true;

        grpCols = {'epoch_idx','epoch_name'};
        if isfield(metricOut, 'groupBy') && ~isempty(metricOut.groupBy)
            grpCols = unique([cellstr(string(metricOut.groupBy(:)')) grpCols], 'stable');
        end
        grpCols = grpCols(ismember(grpCols, T.Properties.VariableNames));
        [G, grp] = groupRows(T, grpCols);
        v = T.(metricCol);
        nGrp = max([G;0]);
        meanVals = NaN(nGrp,1); medVals = NaN(nGrp,1); stdVals = NaN(nGrp,1); nVals = zeros(nGrp,1);
        for k = 1:nGrp
            x = v(G == k);
            x = x(isfinite(x));
            nVals(k) = numel(x);
            if isempty(x), continue; end
            meanVals(k) = mean(x);
            medVals(k) = median(x);
            stdVals(k) = std(x);
        end
        S = grp;
        S.(['mean_' metricCol]) = meanVals;
        S.(['median_' metricCol]) = medVals;
        S.(['std_' metricCol]) = stdVals;
        S.n_events = nVals;
        S.sem = stdVals ./ sqrt(max(nVals,1));
        out.summaryByGroupEpoch = S;
    end
end

if ~didEventLevel && isfield(metricOut, 'summaryByGroupEpoch') && ~isempty(metricOut.summaryByGroupEpoch)
    out.summaryByGroupEpoch = normalizeMetricSummaryOnly(metricOut.summaryByGroupEpoch, metricCol, mode, transform, baseIdx);
end
out.normalization = struct('mode', mode, 'transform', transform, 'baselineEpochIdx', baseIdx);
end

function base = baselineVector(T, vals, mode, transform, baseIdx, keyCol)
base = NaN(size(vals));
if isempty(vals) || ~ismember('epoch_idx', T.Properties.VariableNames)
    return
end

if strcmp(mode, 'by_global_baseline_mean')
    b = bwaComputeBaselineValue(vals, T.epoch_idx == baseIdx, transform);
    base(:) = b;
    return
end

if strcmp(mode, 'by_cell_baseline_mean') && ismember(keyCol, T.Properties.VariableNames)
    keys = string(T.(keyCol));
    u = unique(keys, 'stable');
    for i = 1:numel(u)
        k = u(i);
        if strlength(k) == 0
            continue
        end
        idx = keys == k;
        b = bwaComputeBaselineValue(vals, idx & T.epoch_idx == baseIdx, transform);
        base(idx) = b;
    end
else
    b = bwaComputeBaselineValue(vals, T.epoch_idx == baseIdx, transform);
    base(:) = b;
end
end

function keyCol = selectCellKeyColumn(T)
% Prefer globally unique cell identifiers when available.
if ismember('cell_uid', T.Properties.VariableNames)
    keyCol = 'cell_uid';
elseif ismember('primary_cell_uid', T.Properties.VariableNames)
    keyCol = 'primary_cell_uid';
else
    keyCol = '';
end
end

function out = applyTransform(vals, base, transform)
out = bwaApplyTransform(vals, base, transform);
end

function S = normalizeMetricSummaryOnly(S, metricCol, mode, transform, baseIdx)
meanCol = ['mean_' metricCol];
medCol = ['median_' metricCol];
if ~ismember('epoch_idx', S.Properties.VariableNames) || ~ismember(meanCol, S.Properties.VariableNames)
    return
end

groupCols = {};
for c = {'animal_id','slice_num','file_id'}
    if ismember(c{1}, S.Properties.VariableNames)
        groupCols{end+1} = c{1}; %#ok<AGROW>
    end
end

if strcmp(mode, 'by_global_baseline_mean') || isempty(groupCols)
    b = bwaComputeBaselineValue(S.(meanCol), S.epoch_idx == baseIdx, transform);
    S.(meanCol) = applyTransform(S.(meanCol), b*ones(height(S),1), transform);
    if ismember(medCol, S.Properties.VariableNames)
        S.(medCol) = applyTransform(S.(medCol), b*ones(height(S),1), transform);
    end
    if ismember(['std_' metricCol], S.Properties.VariableNames)
        S.(['std_' metricCol])(:) = NaN;
    end
    if ismember('sem', S.Properties.VariableNames)
        S.sem(:) = NaN;
    end
    return
end

[G, grp] = groupRows(S, groupCols);
bv = NaN(height(S),1);
for i = 1:height(grp)
    idx = G == i;
    b = bwaComputeBaselineValue(S.(meanCol), idx & S.epoch_idx == baseIdx, transform);
    bv(idx) = b;
end
S.(meanCol) = applyTransform(S.(meanCol), bv, transform);
if ismember(medCol, S.Properties.VariableNames)
    S.(medCol) = applyTransform(S.(medCol), bv, transform);
end
if ismember(['std_' metricCol], S.Properties.VariableNames)
    S.(['std_' metricCol])(:) = NaN;
end
if ismember('sem', S.Properties.VariableNames)
    S.sem(:) = NaN;
end
end

function idx = resolveBaselineEpochIdx(epochTable, baselineEpoch)
idx = 1;
if isnumeric(baselineEpoch) && isscalar(baselineEpoch)
    idx = baselineEpoch;
    return
end
if isstring(baselineEpoch) || ischar(baselineEpoch)
    b = string(baselineEpoch);
    if ismember('epoch_name', epochTable.Properties.VariableNames)
        names = string(epochTable.epoch_name);
        hit = find(strcmpi(names, b), 1, 'first');
        if ~isempty(hit) && ismember('epoch_idx', epochTable.Properties.VariableNames)
            idx = epochTable.epoch_idx(hit);
        end
    end
end
end

function [G, grpTbl] = groupRows(T, cols)
n = height(T);
if n == 0
    G = zeros(0,1);
    grpTbl = T(:, cols);
    return
end
K = strings(n,1);
for i = 1:numel(cols)
    K = K + "|" + toKey(T.(cols{i}));
end
[~,~,G] = unique(K, 'stable');
G = full(double(G(:)));
firstIdx = accumarray(G, (1:n)', [], @min);
grpTbl = T(firstIdx, cols);
end

function s = toKey(v)
if isstring(v)
    s = v(:);
elseif iscell(v)
    s = strings(numel(v),1);
    for i = 1:numel(v)
        s(i) = string(v{i});
    end
elseif isnumeric(v) || islogical(v) || iscategorical(v)
    s = string(v(:));
else
    s = string(v);
    s = s(:);
end
end
