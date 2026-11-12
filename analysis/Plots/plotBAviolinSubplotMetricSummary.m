function out = plotBAviolinSubplotMetricSummary(comp, timeFrameRangeEpochs, datTxt, varargin)
% Plot compact epoch summary for one BAviolin metric (normalized domain).
% Shows mean+-SEM by epoch for event-level and cell-level normalized values.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'comp', @isstruct);
addRequired(p, 'timeFrameRangeEpochs', @isstruct);
addRequired(p, 'datTxt', @(x) ischar(x) || isstring(x));
addParameter(p, 'figure', [], @(x) isempty(x) || ishghandle(x, 'figure'));
addParameter(p, 'axes', [], @(x) isempty(x) || ishghandle(x, 'axes'));
addParameter(p, 'style', struct(), @isstruct);
parse(p, comp, timeFrameRangeEpochs, datTxt, varargin{:});
style = mergeBWAPlotStyle(getBWAPlotStyle(), p.Results.style);

[h_fig, h_ax] = resolvePlotAxes('figure', p.Results.figure, 'axes', p.Results.axes, ...
    'name', ['BA_Summary_' char(string(datTxt))]);
figure(h_fig); axes(h_ax); cla(h_ax);
set(h_ax, 'FontSize', style.fontSizeAxis);

nEpoch = numel(timeFrameRangeEpochs.epochNames);
x = 1:nEpoch;

muE = NaN(1, nEpoch); seE = NaN(1, nEpoch);
if isfield(comp, 'normDFFs') && isfield(comp, 'outEpochIndex')
    vals = double(comp.normDFFs(:));
    ep = double(comp.outEpochIndex(:));
    k = isfinite(vals) & isfinite(ep) & ep >= 1 & ep <= nEpoch;
    vals = vals(k); ep = ep(k);
    for e = 1:nEpoch
        ve = vals(ep == e);
        ve = ve(isfinite(ve));
        if ~isempty(ve)
            muE(e) = mean(ve);
            seE(e) = std(ve, 0) / sqrt(numel(ve));
        end
    end
end

muC = NaN(1, nEpoch); seC = NaN(1, nEpoch);
if isfield(comp, 'avgsn_byCell') && isfield(comp, 'tCatn_byCell')
    vc = double(comp.avgsn_byCell(:));
    tc = string(comp.tCatn_byCell(:));
    names = string(timeFrameRangeEpochs.epochNames(:));
    for e = 1:nEpoch
        ve = vc(tc == names(e));
        ve = ve(isfinite(ve));
        if ~isempty(ve)
            muC(e) = mean(ve);
            seC(e) = std(ve, 0) / sqrt(numel(ve));
        end
    end
end

hold(h_ax, 'on');
yline(h_ax, 0, '-', 'Color', [0 0 0], 'LineWidth', style.lineWidthBaseline);
he = errorbar(h_ax, x, muE, seE, '-o', 'Color', [0.35 0.35 0.35], ...
    'LineWidth', max(1, style.lineWidthMain-0.2), 'MarkerSize', 4, 'DisplayName', 'event mean+-SEM');
hc = errorbar(h_ax, x, muC, seC, '-o', 'Color', style.colors.groupC, ...
    'LineWidth', style.lineWidthMain, 'MarkerSize', 5, 'DisplayName', 'cell mean+-SEM');
hold(h_ax, 'off');

xlim(h_ax, [0.5 nEpoch + 0.5]);
xticks(h_ax, x);
xticklabels(h_ax, timeFrameRangeEpochs.epochNames);
xtickangle(h_ax, 20);
ylabel(h_ax, '${\sf log}_2{\sf FC}$', 'Interpreter', 'latex');
title(h_ax, sprintf('%s: epoch summary', char(string(datTxt))), 'Interpreter', 'none');
grid(h_ax, 'on');
legend(h_ax, [he hc], 'Location', 'southoutside', 'Interpreter', 'none', 'Box', 'off');

out = struct('figure', h_fig, 'axes', h_ax, 'eventMean', muE, 'eventSEM', seE, ...
    'cellMean', muC, 'cellSEM', seC);
end
