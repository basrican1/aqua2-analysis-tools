function mOut = bwaBuildMetricOutFromTable(T, metricCol, epochTable, varargin)
% Convert a metric table into canonical-like summary/event structs.
%
% Usage:
%   mOut = bwaBuildMetricOutFromTable(T, metricCol, epochTable)
%   mOut = bwaBuildMetricOutFromTable(T, metricCol, epochTable, ...
%       ''entityKeys'', {''primary_cell_uid''})

p = inputParser;
p.StructExpand = false;
addRequired(p, 'T', @istable);
addRequired(p, 'metricCol', @(x) ischar(x) || isstring(x));
addRequired(p, 'epochTable', @(x) isempty(x) || istable(x));
addParameter(p, 'entityKeys', {}, @(x) iscell(x) || isstring(x));
parse(p, T, metricCol, epochTable, varargin{:});
metricCol = char(string(metricCol));
entityKeys = cellstr(string(p.Results.entityKeys));

req = {'epoch_idx','epoch_name'};
if ~all(ismember(req, T.Properties.VariableNames))
    error('Metric table missing epoch columns.');
end
for i = 1:numel(entityKeys)
    ek = entityKeys{i};
    if ~ismember(ek, T.Properties.VariableNames)
        error('Metric table missing entity column: %s', ek);
    end
end

if isempty(epochTable) || ~istable(epochTable) || ~ismember('epoch_idx', epochTable.Properties.VariableNames)
    epochTable = unique(T(:, {'epoch_idx','epoch_name'}), 'rows', 'stable');
end

grp = {'animal_id','slice_num','file_id'};
grp = grp(ismember(grp, T.Properties.VariableNames));
sumCols = unique([grp entityKeys {'epoch_idx','epoch_name'}], 'stable');
[G, Gtbl] = localGroupRows(T, sumCols);
nG = max([G;0]);
mu = NaN(nG,1); md = NaN(nG,1); sd = NaN(nG,1); nn = zeros(nG,1);
vals = double(T.(metricCol));
for k = 1:nG
    x = vals(G == k);
    x = x(isfinite(x));
    nn(k) = numel(x);
    if isempty(x)
        continue
    end
    mu(k) = mean(x);
    md(k) = median(x);
    sd(k) = std(x);
end

S = Gtbl;
S.(['mean_' metricCol]) = mu;
S.(['median_' metricCol]) = md;
S.(['std_' metricCol]) = sd;
S.n_events = nn;
S.sem = sd ./ sqrt(max(nn,1));

mOut = struct();
mOut.metricCol = metricCol;
mOut.summaryByGroupEpoch = S;
mOut.eventValues = T;
mOut.epochTable = epochTable;
mOut.groupBy = grp;
mOut.entityKeys = entityKeys;
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
elseif isnumeric(v) || islogical(v)
    s = string(v(:));
else
    s = string(v(:));
end
end
