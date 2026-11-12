function out = adaptCanonicalFrequencyFromCellMetric(cellFreqTable, templateFreq, varargin)
% Build canonical.frequency-style struct from authoritative cell_frequency table.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'cellFreqTable', @istable);
addRequired(p, 'templateFreq', @isstruct);
addParameter(p, 'valueColumn', 'rate_events_per_min', @(x) ischar(x) || isstring(x));
parse(p, cellFreqTable, templateFreq, varargin{:});

valueCol = char(string(p.Results.valueColumn));
if ~ismember(valueCol, cellFreqTable.Properties.VariableNames)
    error('adaptCanonicalFrequencyFromCellMetric:MissingValueColumn', ...
        'Value column %s missing from cell frequency table.', valueCol);
end

T = cellFreqTable;
if ismember('primary_cell_uid', T.Properties.VariableNames) && ~ismember('cell_uid', T.Properties.VariableNames)
    T.cell_uid = string(T.primary_cell_uid);
end
if ~ismember('cell_uid', T.Properties.VariableNames)
    error('adaptCanonicalFrequencyFromCellMetric:MissingCellUID', ...
        'cell_frequency table must include primary_cell_uid or cell_uid.');
end
T.cell_uid = string(T.cell_uid);

rates = T;
if ismember('rate_events_per_min_norm', rates.Properties.VariableNames)
    rates.rate_events_per_min_norm = [];
end
rates.rate_events_per_min = double(T.(valueCol));

out = templateFreq;
out.ratesByCellEpoch = rates;
out.summaryByGroupEpoch = summarizeRates(rates, templateFreq);
out.epochTable = pickEpochTable(T, templateFreq);

if ~isfield(out, 'groupBy') || isempty(out.groupBy)
    out.groupBy = {'file_id'};
end
if ~iscell(out.groupBy)
    out.groupBy = cellstr(string(out.groupBy));
end
end

function S = summarizeRates(rates, templateFreq)
grpCols = {'file_id'};
if isfield(templateFreq, 'groupBy') && ~isempty(templateFreq.groupBy)
    grpCols = cellstr(string(templateFreq.groupBy(:)))';
end
for c = grpCols
    if ~ismember(c{1}, rates.Properties.VariableNames)
        rates.(c{1}) = repmat("", height(rates), 1);
    end
end
sumCols = unique([grpCols {'epoch_idx','epoch_name'}], 'stable');
[G, grp] = groupRows(rates, sumCols);
vals = double(rates.rate_events_per_min);
nGrp = max([G;0]);
meanRate = NaN(nGrp,1);
medianRate = NaN(nGrp,1);
stdRate = NaN(nGrp,1);
nCells = zeros(nGrp,1);
for gg = 1:nGrp
    x = vals(G == gg);
    x = x(isfinite(x));
    nCells(gg) = numel(x);
    if isempty(x)
        continue
    end
    meanRate(gg) = mean(x);
    medianRate(gg) = median(x);
    stdRate(gg) = std(x);
end
S = grp;
S.mean_rate_events_per_min = meanRate;
S.median_rate_events_per_min = medianRate;
S.std_rate_events_per_min = stdRate;
S.n_cells = nCells;
S.sem_rate_events_per_min = stdRate ./ sqrt(max(nCells,1));
end

function epochTable = pickEpochTable(T, templateFreq)
if isfield(templateFreq, 'epochTable') && istable(templateFreq.epochTable)
    epochTable = templateFreq.epochTable;
elseif ismember('epoch_idx', T.Properties.VariableNames) && ismember('epoch_name', T.Properties.VariableNames)
    epochTable = unique(T(:, {'epoch_idx','epoch_name'}), 'rows', 'stable');
else
    epochTable = table();
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
