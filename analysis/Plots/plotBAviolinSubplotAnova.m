function out = plotBAviolinSubplotAnova(values, groups, varargin)
% Plot BAviolin ANOVA chart and table.
% Can be used standalone or in a provided figure+axes pair.

p = inputParser;
p.StructExpand = false;
addParameter(p, 'figure', [], @(x) isempty(x) || ishghandle(x, 'figure'));
addParameter(p, 'axesPlot', [], @(x) isempty(x) || ishghandle(x, 'axes'));
addParameter(p, 'axesTable', [], @(x) isempty(x) || ishghandle(x, 'axes'));
addParameter(p, 'figureName', 'BA_ANOVA', @(x) ischar(x) || isstring(x));
addParameter(p, 'style', struct(), @isstruct);
addParameter(p, 'statsConfig', struct(), @isstruct);
addParameter(p, 'pairID', [], @(x) isnumeric(x) || isstring(x) || iscell(x) || isempty(x));
parse(p, varargin{:});
style = mergeBWAPlotStyle(getBWAPlotStyle(), p.Results.style);

[h_fig, axPlot] = resolvePlotAxes('figure', p.Results.figure, 'axes', p.Results.axesPlot, 'name', p.Results.figureName);
axTable = p.Results.axesTable;
if isempty(axTable)
    figure(h_fig);
    axPlot.Position = [0.08 0.12 0.40 0.78];
    axTable = axes('Parent', h_fig, 'Position', [0.52 0.12 0.40 0.78]);
end

% Compute stats by configurable design (paired/unpaired, 2-group/multi-group).
statsCfg = p.Results.statsConfig;
if ~isfield(statsCfg, 'pairID') && ~isempty(p.Results.pairID)
    statsCfg.pairID = p.Results.pairID;
end
if isempty(fieldnames(statsCfg))
    statsRes = runGroupComparison(values, groups);
else
    nv = structToNV(statsCfg);
    statsRes = runGroupComparison(values, groups, nv{:});
end
pval = statsRes.pValue;

% Draw directly on target axes.
figure(h_fig); axes(axPlot); cla(axPlot);
boxplot(axPlot, values, groups, 'Symbol', '.');
set(axPlot, 'FontSize', style.fontSizeAxis);
ylabel(axPlot, 'value');
title(axPlot, sprintf('%s (p=%.3g)', statsRes.testName, pval));
h_an = findall(axPlot);

figure(h_fig); axes(axTable);
posT = get(axTable, 'Position');
delete(axTable);

uitableHandles = findall(h_fig, 'Type', 'uitable');
for i = 1:numel(uitableHandles)
    handlePos = get(uitableHandles(i), 'Position');
    if max(abs(posT - handlePos)) < 1e-10
        delete(uitableHandles(i));
        break;
    end
end
sumTbl = statsRes.summaryTable;
sumTbl.pValue = repmat(pval, height(sumTbl), 1);
sumTbl = normalizeUitableTable(sumTbl);
ht = uitable('Parent', h_fig, 'Units', 'normalized', 'Position', posT, ...
    'Data', table2cell(sumTbl), 'ColumnName', sumTbl.Properties.VariableNames, ...
    'RowName', {}, 'FontSize', style.fontSizeAxis);

out = struct('figure', h_fig, 'axesPlot', axPlot, 'uitable', ht, ...
    'stats', statsRes, 'anovaGraphics', h_an, ...
    'tableCell', {table2cell(sumTbl)}, 'pvalue', pval);
end

function nv = structToNV(s)
if isempty(s)
    nv = {};
    return
end
f = fieldnames(s);
nv = cell(1, 2*numel(f));
for i = 1:numel(f)
    nv{2*i-1} = f{i};
    nv{2*i} = s.(f{i});
end
end

function T = normalizeUitableTable(T)
vars = T.Properties.VariableNames;
for i = 1:numel(vars)
    v = T.(vars{i});
    if isstring(v)
        T.(vars{i}) = cellstr(v);
    elseif iscategorical(v)
        T.(vars{i}) = cellstr(string(v));
    end
end
end
