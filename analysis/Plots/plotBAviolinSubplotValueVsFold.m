function out = plotBAviolinSubplotValueVsFold(c1, m1, c2, m2, datTxt, varargin)
% Plot BAviolin mean value vs mean fold-change summary subplot.

p = inputParser;
p.StructExpand = false;
addParameter(p, 'figure', [], @(x) isempty(x) || ishghandle(x, 'figure'));
addParameter(p, 'axes', [], @(x) isempty(x) || ishghandle(x, 'axes'));
addParameter(p, 'figureName', ['BA_ValueVsFold_' char(string(datTxt))], @(x) ischar(x) || isstring(x));
addParameter(p, 'style', struct(), @isstruct);
addParameter(p, 'groupNames', {}, @(x) iscell(x) || isstring(x));
parse(p, varargin{:});
style = mergeBWAPlotStyle(getBWAPlotStyle(), p.Results.style);

[h_fig, h_ax] = resolvePlotAxes('figure', p.Results.figure, 'axes', p.Results.axes, 'name', p.Results.figureName);
figure(h_fig); axes(h_ax); cla(h_ax); hold(h_ax, 'on');
set(h_ax, 'FontSize', style.fontSizeAxis);

hw1 = halfwidthFromMultCmpLocal(c1, m1);
hw2 = halfwidthFromMultCmpLocal(c2, m2);
for ii = 1:length(m1)
    c = epochColorFromGroup(ii, p.Results.groupNames);
    errorbar(m1(ii,1), m2(ii,1), hw2(ii), hw2(ii), hw1(ii), hw1(ii), 'o', ...
        'LineWidth', style.lineWidthError, 'MarkerSize', 6, 'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
end
hold(h_ax, 'off');

if strlength(string(datTxt)) >= 2 && strcmp(datTxt(end-1:end), '_c')
    ymode = 'cell_norm';
else
    ymode = 'norm';
end
xlabel(h_ax, bwaMetricLabelLatex(datTxt, 'raw'), 'Interpreter', 'latex');
ylabel(h_ax, bwaMetricLabelLatex(datTxt, ymode), 'Interpreter', 'latex');

out = struct('figure', h_fig, 'axes', h_ax, 'hwX', hw1, 'hwY', hw2);
end

function c = epochColorFromGroup(idx, groupNames)
if iscell(groupNames) || isstring(groupNames)
    g = string(groupNames);
    if idx <= numel(g)
        name = lower(strtrim(char(g(idx))));
        if contains(name, 'base')
            c = [0.00 0.45 0.74]; return
        elseif contains(name, 'cno1') || contains(name, 'drug1') || contains(name, 'epoch2')
            c = [0.85 0.33 0.10]; return
        elseif contains(name, 'cno2') || contains(name, 'drug2') || contains(name, 'epoch3')
            c = [0.93 0.69 0.13]; return
        end
    end
end
pal = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13];
idx = max(1, min(size(pal,1), round(double(idx))));
c = pal(idx,:);
end

function hw = halfwidthFromMultCmpLocal(c, m)
lower = c(:,3); upper = c(:,5);
margin = (upper - lower) / 2;
se = m(:,2);
crit = margin ./ se;

ng = length(m);
i12 = sub2ind([ng ng], c(:,1), c(:,2));
d = zeros(ng, ng);
d(i12) = m(:,2);
sum1 = sum(sum(d));
d = d + d';
sum2 = sum(d);
if ng > 2
    w = ((ng-1) * sum2 - sum1) ./ ((ng-1) * (ng-2));
else
    w = sum1 * ones(2, 1) / 2;
end
hw = crit .* w(:);
end
