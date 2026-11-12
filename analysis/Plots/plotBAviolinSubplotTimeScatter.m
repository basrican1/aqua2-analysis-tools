function out = plotBAviolinSubplotTimeScatter(dffAll, tOccuring, memberlist, datTxt, timeFrameRangeEpochs, varargin)
% Plot BAviolin time scatter subplot.
% Can be used standalone or in a parent figure/axes.

p = inputParser;
p.StructExpand = false;
addParameter(p, 'figure', [], @(x) isempty(x) || ishghandle(x, 'figure'));
addParameter(p, 'axes', [], @(x) isempty(x) || ishghandle(x, 'axes'));
addParameter(p, 'style', struct(), @isstruct);
parse(p, varargin{:});
style = mergeBWAPlotStyle(getBWAPlotStyle(), p.Results.style);

[h_fig, h_ax] = resolvePlotAxes('figure', p.Results.figure, 'axes', p.Results.axes, ...
    'name', ['BA_Time_' char(string(datTxt))]);
figure(h_fig); axes(h_ax); cla(h_ax);
set(h_ax, 'FontSize', style.fontSizeAxis);
colors = get(h_ax, 'ColorOrder');
scatterHandle = gscatter(h_ax, tOccuring, dffAll, memberlist, [], [], [], 'off');
for ii = 1:numel(scatterHandle)
    try
        scatterHandle(ii).MarkerSize = style.scatterSize;
    catch
    end
end
xlabel(h_ax, 'time (m)');
ylabel(h_ax, bwaMetricLabelLatex(datTxt, inferLabelMode(datTxt)), 'Interpreter', 'latex');

yl = ylim(h_ax);
nEpochs = size(timeFrameRangeEpochs.frameRange, 1);
labelHandles = gobjects(0, 1);
for eidx = 1:nEpochs
    tstart = timeFrameRangeEpochs.timeRange(eidx,1);
    tend = timeFrameRangeEpochs.timeRange(eidx,2);
    xregion(h_ax, tstart, tend, 'FaceColor', colors(eidx, :));
    txt = timeFrameRangeEpochs.epochNames{eidx};
    xCenter = (tstart + tend) / 2;
    labelHandles(end+1,1) = text(xCenter, 0.95 * yl(2), txt, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
        'Tag', 'xregionLabel', 'FontSize', style.fontSizeEpochLabel, 'Color', style.colors.text, 'Parent', h_ax); %#ok<AGROW>
end

out = struct('figure', h_fig, 'axes', h_ax, 'scatter', scatterHandle, 'epochLabels', labelHandles);
end

function mode = inferLabelMode(datTxt)
dt = lower(strtrim(char(string(datTxt))));
if startsWith(dt, 'norm ')
    mode = 'norm';
else
    mode = 'raw';
end
end
