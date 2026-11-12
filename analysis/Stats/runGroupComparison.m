function out = runGroupComparison(values, groups, varargin)
% Flexible group-comparison engine for plotting and summaries.
%
% Supports auto test selection across:
%   - 2 groups: ttest2/ranksum (unpaired), ttest/signrank (paired)
%   - >2 groups: anova1/kruskalwallis (unpaired), friedman (paired, complete blocks)
%
% Name-value:
%   'method'      : 'auto' (default), 'anova', 'kruskal', 'ttest2',
%                   'ttest_paired', 'ranksum', 'signrank', 'friedman'
%   'paired'      : logical (default false)
%   'parametric'  : logical (default true)
%   'pairID'      : vector matching values (required for paired tests)
%   'groupOrder'  : desired group order (cellstr/string/categorical)

p = inputParser;
p.StructExpand = false;
addRequired(p, 'values', @isnumeric);
addRequired(p, 'groups');
addParameter(p, 'method', 'auto', @(x) ischar(x) || isstring(x));
addParameter(p, 'paired', false, @islogical);
addParameter(p, 'parametric', true, @islogical);
addParameter(p, 'pairID', [], @(x) isnumeric(x) || isstring(x) || iscell(x) || isempty(x));
addParameter(p, 'groupOrder', [], @(x) iscell(x) || isstring(x) || iscategorical(x) || isempty(x));
parse(p, values, groups, varargin{:});

y = values(:);
g = normalizeGroup(grouprow(groups), p.Results.groupOrder);
pairID = normalizePairID(p.Results.pairID, numel(y));

keep = isfinite(y) & ~ismissing(g);
if ~isempty(pairID)
    keep = keep & ~ismissing(pairID);
end
y = y(keep);
g = g(keep);
if ~isempty(pairID)
    pairID = pairID(keep);
end

groupNames = categories(g);
nGroups = numel(groupNames);
if nGroups < 2
    error('runGroupComparison:NotEnoughGroups', 'Need at least 2 non-empty groups.');
end

method = lower(char(string(p.Results.method)));
if strcmp(method, 'auto')
    method = chooseAutoMethod(nGroups, p.Results.paired, p.Results.parametric);
end

out = struct();
out.method = method;
out.groupNames = groupNames;
out.nGroups = nGroups;
out.nObs = numel(y);
out.summaryTable = groupSummary(y, g, groupNames);

switch method
    case 'anova'
        [pval, tbl, stats] = anova1(y, g, 'off');
        [c, m, ~, gnames] = multcompare(stats, 'Display', 'off');
        out = fillResult(out, 'One-way ANOVA', pval, tbl, c, m, gnames);

    case 'kruskal'
        [pval, tbl, stats] = kruskalwallis(y, g, 'off');
        [c, m, ~, gnames] = multcompare(stats, 'Display', 'off');
        out = fillResult(out, 'Kruskal-Wallis', pval, tbl, c, m, gnames);

    case 'ttest2'
        [c, m, pval, tbl, gnames] = runTwoGroupUnpaired(y, g, true);
        out = fillResult(out, 'Two-sample t-test', pval, tbl, c, m, gnames);

    case 'ranksum'
        [c, m, pval, tbl, gnames] = runTwoGroupUnpaired(y, g, false);
        out = fillResult(out, 'Wilcoxon rank-sum', pval, tbl, c, m, gnames);

    case 'ttest_paired'
        [c, m, pval, tbl, gnames] = runTwoGroupPaired(y, g, pairID, true);
        out = fillResult(out, 'Paired t-test', pval, tbl, c, m, gnames);

    case 'signrank'
        [c, m, pval, tbl, gnames] = runTwoGroupPaired(y, g, pairID, false);
        out = fillResult(out, 'Wilcoxon signed-rank', pval, tbl, c, m, gnames);

    case 'friedman'
        [pval, tbl, stats, gnames] = runFriedman(y, g, pairID);
        [c, m] = multcompare(stats, 'Display', 'off');
        out = fillResult(out, 'Friedman', pval, tbl, c, m, gnames);

    otherwise
        error('runGroupComparison:UnknownMethod', 'Unknown method: %s', method);
end
end

function method = chooseAutoMethod(nGroups, paired, parametric)
if nGroups == 2
    if paired
        if parametric
            method = 'ttest_paired';
        else
            method = 'signrank';
        end
    else
        if parametric
            method = 'ttest2';
        else
            method = 'ranksum';
        end
    end
else
    if paired
        method = 'friedman';
    else
        if parametric
            method = 'anova';
        else
            method = 'kruskal';
        end
    end
end
end

function g = grouprow(groups)
if iscategorical(groups)
    g = groups(:);
elseif iscell(groups) || isstring(groups)
    g = categorical(groups(:));
else
    g = categorical(groups(:));
end
end

function g = normalizeGroup(g, groupOrder)
if isempty(groupOrder)
    return
end
ord = cellstr(string(groupOrder(:)));
g = categorical(cellstr(string(g)), ord, 'Ordinal', true);
end

function pairID = normalizePairID(pairIDIn, n)
if isempty(pairIDIn)
    pairID = [];
    return
end
if isnumeric(pairIDIn)
    pairID = string(pairIDIn(:));
elseif isstring(pairIDIn)
    pairID = pairIDIn(:);
elseif iscell(pairIDIn)
    pairID = string(pairIDIn(:));
else
    error('runGroupComparison:BadPairID', 'Unsupported pairID type.');
end
if numel(pairID) ~= n
    error('runGroupComparison:PairIDSize', 'pairID must match number of values.');
end
end

function out = fillResult(out, testName, pval, tbl, c, m, gnames)
out.testName = testName;
out.pValue = pval;
out.omnibusTable = tbl;
out.pairwiseC = c;
out.groupStatsM = m;
out.groupNames = cellstr(string(gnames(:)));
end

function T = groupSummary(y, g, gnames)
n = numel(gnames);
means = NaN(n,1);
sem = NaN(n,1);
ns = zeros(n,1);
for i = 1:n
    yi = y(g == gnames{i});
    yi = yi(isfinite(yi));
    ns(i) = numel(yi);
    if ns(i) > 0
        means(i) = mean(yi);
    end
    if ns(i) > 1
        sem(i) = std(yi) / sqrt(ns(i));
    else
        sem(i) = NaN;
    end
end
T = table(string(gnames(:)), ns, means, sem, 'VariableNames', {'group', 'n', 'mean', 'sem'});
end

function [c, m, pval, tbl, gnames] = runTwoGroupUnpaired(y, g, parametric)
gnames = categories(g);
if numel(gnames) ~= 2
    error('runGroupComparison:TwoGroupNeeded', 'Selected method requires exactly 2 groups.');
end
g1 = y(g == gnames{1});
g2 = y(g == gnames{2});
if parametric
    [~, pval, ci] = ttest2(g1, g2);
    methodName = 'ttest2';
else
    pval = ranksum(g1, g2);
    ci = [NaN NaN];
    methodName = 'ranksum';
end
diffMean = mean(g1, 'omitnan') - mean(g2, 'omitnan');
c = [1 2 ci(1) diffMean ci(2) pval];
m = mkGroupStatsForTwo(g1, g2);
tbl = {methodName, 'pValue'; 'overall', pval};
end

function [c, m, pval, tbl, gnames] = runTwoGroupPaired(y, g, pairID, parametric)
if isempty(pairID)
    error('runGroupComparison:PairIDRequired', 'Paired test requires pairID.');
end
gnames = categories(g);
if numel(gnames) ~= 2
    error('runGroupComparison:TwoGroupNeeded', 'Selected method requires exactly 2 groups.');
end
[x1, x2] = pairedVectors(y, g, pairID, gnames);
if isempty(x1)
    error('runGroupComparison:NoPairs', 'No complete pairs found for paired test.');
end
if parametric
    [~, pval, ci] = ttest(x1, x2);
    methodName = 'paired_ttest';
else
    pval = signrank(x1, x2);
    ci = [NaN NaN];
    methodName = 'signrank';
end
diffMean = mean(x1 - x2, 'omitnan');
c = [1 2 ci(1) diffMean ci(2) pval];
m = mkGroupStatsForTwo(x1, x2);
tbl = {methodName, 'pValue'; 'overall', pval};
end

function [pval, tbl, stats, gnames] = runFriedman(y, g, pairID)
if isempty(pairID)
    error('runGroupComparison:PairIDRequired', 'Friedman test requires pairID.');
end
gnames = categories(g);
ids = unique(pairID, 'stable');
nID = numel(ids);
nG = numel(gnames);
X = NaN(nID, nG);
for i = 1:nID
    for j = 1:nG
        mask = pairID == ids(i) & g == gnames{j};
        vals = y(mask);
        if ~isempty(vals)
            X(i,j) = mean(vals, 'omitnan');
        end
    end
end
complete = all(isfinite(X), 2);
Xc = X(complete, :);
if isempty(Xc)
    error('runGroupComparison:FriedmanNoComplete', 'No complete paired blocks found for Friedman test.');
end
[pval, tbl, stats] = friedman(Xc, 1, 'off');
end

function [x1, x2] = pairedVectors(y, g, pairID, gnames)
ids = unique(pairID, 'stable');
x1 = [];
x2 = [];
for i = 1:numel(ids)
    a = y(pairID == ids(i) & g == gnames{1});
    b = y(pairID == ids(i) & g == gnames{2});
    if ~isempty(a) && ~isempty(b)
        x1(end+1,1) = mean(a, 'omitnan'); %#ok<AGROW>
        x2(end+1,1) = mean(b, 'omitnan'); %#ok<AGROW>
    end
end
end

function m = mkGroupStatsForTwo(g1, g2)
n1 = numel(g1);
n2 = numel(g2);
se1 = std(g1, 'omitnan') / max(sqrt(n1), 1);
se2 = std(g2, 'omitnan') / max(sqrt(n2), 1);
m = [1 mean(g1, 'omitnan') se1; 2 mean(g2, 'omitnan') se2];
end
