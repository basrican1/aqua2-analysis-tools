function out = plotBAviolinSubplotMultcompare(statsIn, varargin)
% Plot BAviolin multcompare subplot.
% Can be used standalone or in a parent figure/axes.

p = inputParser;
p.StructExpand = false;
addParameter(p, 'figure', [], @(x) isempty(x) || ishghandle(x, 'figure'));
addParameter(p, 'axes', [], @(x) isempty(x) || ishghandle(x, 'axes'));
addParameter(p, 'figureName', 'BA_Multcompare', @(x) ischar(x) || isstring(x));
addParameter(p, 'style', struct(), @isstruct);
parse(p, varargin{:});
style = mergeBWAPlotStyle(getBWAPlotStyle(), p.Results.style);

[h_fig, h_ax] = resolvePlotAxes('figure', p.Results.figure, 'axes', p.Results.axes, 'name', p.Results.figureName);
figure(h_fig); axes(h_ax); cla(h_ax);
set(h_ax, 'FontSize', style.fontSizeAxis);

[c, m, gnames, testName] = getPairwiseFromStats(statsIn);

if isempty(c) || isempty(m)
    text(h_ax, 0.5, 0.5, 'No pairwise comparisons available', ...
        'HorizontalAlignment', 'center', 'FontSize', style.fontSizeAxis);
    axis(h_ax, 'off');
    out = struct('figure', h_fig, 'axes', h_ax, 'c', [], 'm', [], ...
        'graphics', findall(h_ax), 'groupNames', {gnames});
    return
end

% Draw directly on target axes.
nCmp = size(c, 1);
y = 1:nCmp;
lo = c(:,3);
mu = c(:,4);
hi = c(:,5);
hold(h_ax, 'on');
for i = 1:nCmp
    line(h_ax, [lo(i) hi(i)], [y(i) y(i)], ...
        'Color', style.colors.groupA, 'LineWidth', style.lineWidthMain);
end
plot(h_ax, mu, y, 'o', ...
    'MarkerSize', style.markerSize, ...
    'MarkerEdgeColor', style.colors.groupC, ...
    'MarkerFaceColor', style.colors.groupC, ...
    'LineStyle', 'none');
xline(h_ax, 0, '--', 'Color', style.colors.baseline, 'LineWidth', style.lineWidthBaseline);
hold(h_ax, 'off');

cmpLabels = strings(nCmp,1);
for i = 1:nCmp
    g1 = c(i,1);
    g2 = c(i,2);
    if g1 >= 1 && g1 <= numel(gnames) && g2 >= 1 && g2 <= numel(gnames)
        cmpLabels(i) = string(gnames{g1}) + " - " + string(gnames{g2});
    else
        cmpLabels(i) = "cmp " + i;
    end
end
set(h_ax, 'YTick', y, 'YTickLabel', cellstr(cmpLabels), 'YDir', 'reverse');
xlabel(h_ax, 'Difference in means (95% CI)');
ylabel(h_ax, 'Group comparison');
title(h_ax, ['Multiple comparisons: ' testName]);

out = struct('figure', h_fig, 'axes', h_ax, 'c', c, 'm', m, ...
    'graphics', findall(h_ax), 'groupNames', {gnames});
end

function [c, m, gnames, testName] = getPairwiseFromStats(statsIn)
c = [];
m = [];
gnames = {};
testName = '';

if isstruct(statsIn) && isfield(statsIn, 'pairwiseC') && isfield(statsIn, 'groupStatsM')
    c = statsIn.pairwiseC;
    m = statsIn.groupStatsM;
    if isfield(statsIn, 'groupNames')
        gnames = cellstr(string(statsIn.groupNames));
    else
        gnames = {};
    end
    if isfield(statsIn, 'testName')
        testName = char(string(statsIn.testName));
    end
    return
end

% Backward compatibility: raw stats struct from anova/kruskal/friedman.
if isstruct(statsIn)
    [c, m, ~, gnames] = multcompare(statsIn, 'Display', 'off');
    testName = 'legacy';
end
end
