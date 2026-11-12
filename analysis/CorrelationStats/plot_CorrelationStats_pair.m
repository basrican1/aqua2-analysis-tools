function out = plot_CorrelationStats_pair(source, metricA, metricB, varargin)
%PLOT_CORRELATIONSTATS_PAIR Scatter drill-down for a selected metric pair.
%
% Usage:
%   out = plot_CorrelationStats_pair(corrOut, 'area_um2', 'frequency')
%   out = plot_CorrelationStats_pair(result, 'area_um2', 'slice.frac_active_cells', ...
%       'unitLevel', 'slice', 'epochName', 'baseline')
%
% Name-value:
%   'figureName'        : figure name override
%   'showFit'           : add least-squares fit line (default true)
%   'labelPoints'       : label points with entity_id (default false)
%   'colorByDenominator': 'auto' (default), false, or one of {'x','y','min','max'}
%   plus any CorrelationStats name-value args when source is not already a CorrelationStats output

if isstruct(source) && isfield(source, 'source') && strcmp(source.source, 'CorrelationStats')
    corrOut = source;
else
    corrOut = CorrelationStats(source, varargin{:});
    varargin = {};
end

p = inputParser;
p.StructExpand = false;
addParameter(p, 'figureName', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'showFit', true, @islogical);
addParameter(p, 'labelPoints', false, @islogical);
addParameter(p, 'colorByDenominator', 'auto', @(x) islogical(x) || ischar(x) || isstring(x));
parse(p, varargin{:});

metricA = char(string(metricA));
metricB = char(string(metricB));
varA = findMetricVar(corrOut, metricA);
varB = findMetricVar(corrOut, metricB);

T = corrOut.wideTable;
x = double(T.(varA));
y = double(T.(varB));
keep = isfinite(x) & isfinite(y);
if nnz(keep) < 2
    error('plot_CorrelationStats_pair:InsufficientData', ...
        'Need at least two finite paired observations for %s vs %s.', metricA, metricB);
end
x = x(keep);
y = y(keep);
T = T(keep, :);

hFig = figure('Color', 'w', 'Visible', 'on');
if strlength(string(p.Results.figureName)) > 0
    set(hFig, 'Name', char(string(p.Results.figureName)));
else
    set(hFig, 'Name', sprintf('Correlation Pair: %s vs %s', metricA, metricB));
end
ax = axes('Parent', hFig);
[cvals, cLabel] = resolveDenominatorColor(corrOut, T, metricA, metricB, p.Results.colorByDenominator);
if isempty(cvals)
    scatter(ax, x, y, 42, 'filled', 'MarkerFaceAlpha', 0.78);
else
    scatter(ax, x, y, 42, cvals, 'filled', 'MarkerFaceAlpha', 0.82);
    cb = colorbar(ax);
    ylabel(cb, cLabel, 'Interpreter', 'none');
    colormap(ax, parula);
end
hold(ax, 'on');
if p.Results.showFit && numel(x) >= 2
    pp = polyfit(x, y, 1);
    xx = linspace(min(x), max(x), 100);
    yy = polyval(pp, xx);
    plot(ax, xx, yy, '-', 'LineWidth', 1.6, 'Color', [0.15 0.15 0.15]);
end
if p.Results.labelPoints && ismember('entity_id', T.Properties.VariableNames)
    ids = string(T.entity_id);
    for i = 1:numel(x)
        text(ax, x(i), y(i), [' ' char(ids(i))], 'Interpreter', 'none', 'FontSize', 8, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    end
end
xlabel(ax, metricA, 'Interpreter', 'none');
ylabel(ax, metricB, 'Interpreter', 'none');
box(ax, 'off');
grid(ax, 'on');
figure(hFig);

[rho, pval] = corr(x, y, 'Rows', 'complete', 'Type', corrOut.method);
title(ax, sprintf('%s vs %s | %s rho=%.3f, p=%.3g, n=%d', ...
    metricA, metricB, capitalizeMethod(corrOut.method), rho, pval, numel(x)), 'Interpreter', 'none');

out = struct();
out.figure = hFig;
out.axes = ax;
out.metricA = metricA;
out.metricB = metricB;
out.x = x;
out.y = y;
keepCols = intersect({'unit_level','entity_id','animal_id','slice_num','file_id','epoch_idx','epoch_name'}, T.Properties.VariableNames);
if ~isempty(cvals)
    T.denominator_color = cvals;
    keepCols = [keepCols {'denominator_color'}];
end
out.table = T(:, keepCols);
out.rho = rho;
out.p_value = pval;
out.method = corrOut.method;
end

function varName = findMetricVar(corrOut, metricID)
metricID = string(metricID);
ids = string(corrOut.metricIDs(:));
idx = find(ids == metricID, 1);
if isempty(idx)
    error('plot_CorrelationStats_pair:MetricNotFound', ...
        'Metric %s not found in CorrelationStats output.', char(metricID));
end
varName = corrOut.metricNames{idx};
end

function out = capitalizeMethod(in)
out = char(string(in));
if ~isempty(out)
    out(1) = upper(out(1));
end
end

function [cvals, cLabel] = resolveDenominatorColor(corrOut, T, metricA, metricB, colorOpt)
cvals = [];
cLabel = '';
if isequal(colorOpt, false)
    return
end
mode = lower(char(string(colorOpt)));
if strcmp(mode, 'true')
    mode = 'auto';
end
if ~isfield(corrOut, 'denomWideTable') || ~istable(corrOut.denomWideTable) || isempty(corrOut.denomWideTable)
    return
end
Dx = fetchDenominator(corrOut, T, metricA);
Dy = fetchDenominator(corrOut, T, metricB);
if strcmp(mode, 'auto')
    if isFractionMetric(metricA) && isFractionMetric(metricB)
        if ~isempty(Dx) && any(isfinite(Dx)) && ~isempty(Dy) && any(isfinite(Dy))
            mode = 'min';
        elseif ~isempty(Dx) && any(isfinite(Dx))
            mode = 'x';
        elseif ~isempty(Dy) && any(isfinite(Dy))
            mode = 'y';
        else
            return
        end
    elseif isFractionMetric(metricA) && ~isempty(Dx) && any(isfinite(Dx))
        mode = 'x';
    elseif isFractionMetric(metricB) && ~isempty(Dy) && any(isfinite(Dy))
        mode = 'y';
    else
        return
    end
end
switch mode
    case 'x'
        if isempty(Dx) || ~any(isfinite(Dx)), return; end
        cvals = Dx; cLabel = sprintf('denominator: %s', metricA);
    case 'y'
        if isempty(Dy) || ~any(isfinite(Dy)), return; end
        cvals = Dy; cLabel = sprintf('denominator: %s', metricB);
    case 'min'
        if isempty(Dx) || isempty(Dy) || ~any(isfinite(Dx)) || ~any(isfinite(Dy)), return; end
        cvals = min([Dx Dy], [], 2); cLabel = 'denominator: min(x,y)';
    case 'max'
        if isempty(Dx) || isempty(Dy) || ~any(isfinite(Dx)) || ~any(isfinite(Dy)), return; end
        cvals = max([Dx Dy], [], 2); cLabel = 'denominator: max(x,y)';
    otherwise
        return
end
if isempty(cvals) || ~any(isfinite(cvals))
    cvals = [];
    cLabel = '';
end
end

function d = fetchDenominator(corrOut, Tplot, metricID)
d = [];
if ~isfield(corrOut, 'denomMetricIDs') || ~isfield(corrOut, 'denomMetricNames')
    return
end
ids = string(corrOut.denomMetricIDs(:));
idx = find(ids == string(metricID), 1);
if isempty(idx)
    return
end
vname = corrOut.denomMetricNames{idx};
D = corrOut.denomWideTable;
keys = intersect({'unit_level','entity_id','animal_id','slice_num','file_id','epoch_idx','epoch_name'}, D.Properties.VariableNames);
keys = keys(ismember(keys, Tplot.Properties.VariableNames));
if ~ismember(vname, D.Properties.VariableNames) || isempty(keys)
    return
end
rhs = D(:, [keys {vname}]);
rhs.Properties.VariableNames{end} = 'denominator_value';
J = outerjoin(Tplot(:, keys), rhs, 'Keys', keys, 'MergeKeys', true, 'Type', 'left');
d = double(J.denominator_value);
end

function tf = isFractionMetric(metricID)
id = lower(char(string(metricID)));
tf = startsWith(id, 'cell.frac_') || startsWith(id, 'slice.frac_') || contains(id, '.frac_');
end
