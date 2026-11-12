function out = plotBAviolinSubplotFoldScatter(dffs, foldChanges, epMembIdx, datTxt, varargin)
% Plot BAviolin value-vs-fold scatter subplot.

p = inputParser;
p.StructExpand = false;
addParameter(p, 'figure', [], @(x) isempty(x) || ishghandle(x, 'figure'));
addParameter(p, 'axes', [], @(x) isempty(x) || ishghandle(x, 'axes'));
addParameter(p, 'figureName', ['BACmpFld_' char(string(datTxt))], @(x) ischar(x) || isstring(x));
addParameter(p, 'style', struct(), @isstruct);
parse(p, varargin{:});
style = mergeBWAPlotStyle(getBWAPlotStyle(), p.Results.style);

[h_fig, h_ax] = resolvePlotAxes('figure', p.Results.figure, 'axes', p.Results.axes, 'name', p.Results.figureName);
figure(h_fig); axes(h_ax); cla(h_ax);
set(h_ax, 'FontSize', style.fontSizeAxis);

epList = unique(epMembIdx(~isnan(epMembIdx))); epList = sort(epList(:));
hold(h_ax, 'on');
for k = 1:numel(epList)
    epID = epList(k);
    theseDFFs = dffs(epMembIdx == epID);
    theseFolds = foldChanges(epMembIdx == epID);
    scatter(h_ax, theseDFFs, theseFolds, 12, 'Marker', '.', 'CData', epochColor(epID));
end
ylims = ylim(h_ax);
yregion(ylims(1), 0);
hold(h_ax, 'off');

suffix = '';
if strlength(string(datTxt)) >= 2 && strcmp(datTxt(end-1:end), '_c')
    suffix = 'cell_norm';
else
    suffix = 'norm';
end
xlabel(h_ax, bwaMetricLabelLatex(datTxt, 'raw'), 'Interpreter', 'latex');
ylabel(h_ax, bwaMetricLabelLatex(datTxt, suffix), 'Interpreter', 'latex');

out = struct('figure', h_fig, 'axes', h_ax, 'epList', epList);
end

function c = epochColor(epID)
pal = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13];
idx = max(1, min(size(pal,1), round(double(epID))));
c = pal(idx,:);
end
