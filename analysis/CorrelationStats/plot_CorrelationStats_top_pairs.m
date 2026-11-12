function out = plot_CorrelationStats_top_pairs(corrOut, summaryOut, varargin)
%PLOT_CORRELATIONSTATS_TOP_PAIRS Plot top positive/negative pairs in tiled figures.
%
% Usage:
%   out = plot_CorrelationStats_top_pairs(corrOut, S)
%   out = plot_CorrelationStats_top_pairs(corrOut, S, 'topK', 6)
%
% Name-value:
%   'topK'               : max rows from S.positive / S.negative (default: all)
%   'which'              : 'both' (default), 'positive', or 'negative'
%   'showFit'            : add fit line (default true)
%   'labelPoints'        : label entity_id (default false)
%   'colorByDenominator' : 'auto' (default), false, or one of {'x','y','min','max'}
%   'figurePrefix'       : optional figure-name prefix
%
% Output:
%   out.positiveFigure
%   out.negativeFigure
%   out.positiveAxes
%   out.negativeAxes

if ~(isstruct(corrOut) && isfield(corrOut, 'source') && strcmp(corrOut.source, 'CorrelationStats'))
    error('plot_CorrelationStats_top_pairs:BadCorrelationInput', ...
        'First input must be a CorrelationStats output struct.');
end
if ~(isstruct(summaryOut) && isfield(summaryOut, 'positive') && isfield(summaryOut, 'negative'))
    error('plot_CorrelationStats_top_pairs:BadSummaryInput', ...
        'Second input must be a summarize_CorrelationStats_top output struct.');
end

p = inputParser;
p.StructExpand = false;
addParameter(p, 'topK', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1));
addParameter(p, 'which', 'both', @(x) any(strcmpi(string(x), ["both","positive","negative"])));
addParameter(p, 'showFit', true, @islogical);
addParameter(p, 'labelPoints', false, @islogical);
addParameter(p, 'colorByDenominator', 'auto', @(x) islogical(x) || ischar(x) || isstring(x));
addParameter(p, 'figurePrefix', '', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});

Tpos = summaryOut.positive;
Tneg = summaryOut.negative;
if ~isempty(p.Results.topK)
    k = double(p.Results.topK);
    if ~isempty(Tpos), Tpos = Tpos(1:min(height(Tpos), k), :); end
    if ~isempty(Tneg), Tneg = Tneg(1:min(height(Tneg), k), :); end
end

whichMode = lower(char(string(p.Results.which)));
out = struct('positiveFigure', [], 'negativeFigure', [], 'positiveAxes', gobjects(0,1), 'negativeAxes', gobjects(0,1));
if any(strcmp(whichMode, {'both','positive'}))
    [out.positiveFigure, out.positiveAxes] = plotPairSet(corrOut, Tpos, 'Positive', p.Results);
end
if any(strcmp(whichMode, {'both','negative'}))
    [out.negativeFigure, out.negativeAxes] = plotPairSet(corrOut, Tneg, 'Negative', p.Results);
end
end

function [hFig, axArr] = plotPairSet(corrOut, Tset, signLabel, opts)
hFig = [];
axArr = gobjects(0,1);
if isempty(Tset)
    return
end
n = height(Tset);
[nr, nc] = bestGrid(n);
prefix = char(string(opts.figurePrefix));
if strlength(string(prefix)) > 0
    fName = sprintf('%s %s Correlation Pairs', prefix, signLabel);
else
    fName = sprintf('%s Correlation Pairs', signLabel);
end
hFig = figure('Color', 'w', 'Visible', 'on', 'Name', fName);
tl = tiledlayout(hFig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s Correlation Pairs', signLabel), 'Interpreter', 'none');
axArr = gobjects(n,1);
for i = 1:n
    ax = nexttile(tl);
    axArr(i) = ax;
    metricA = char(string(Tset.metric_a(i)));
    metricB = char(string(Tset.metric_b(i)));
    drawPairIntoAxes(ax, corrOut, metricA, metricB, opts, Tset(i,:));
end
figure(hFig);
end

function drawPairIntoAxes(ax, corrOut, metricA, metricB, opts, rowT)
varA = findMetricVar(corrOut, metricA);
varB = findMetricVar(corrOut, metricB);
T = corrOut.wideTable;
x = double(T.(varA));
y = double(T.(varB));
keep = isfinite(x) & isfinite(y);
if nnz(keep) < 2
    text(ax, 0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center', 'Interpreter', 'none');
    axis(ax, 'off');
    return
end
x = x(keep);
y = y(keep);
T = T(keep, :);
[cvals, cLabel] = resolveDenominatorColor(corrOut, T, metricA, metricB, opts.colorByDenominator);
if isempty(cvals)
    scatter(ax, x, y, 24, 'filled', 'MarkerFaceAlpha', 0.78);
else
    scatter(ax, x, y, 24, cvals, 'filled', 'MarkerFaceAlpha', 0.82);
    cb = colorbar(ax);
    ylabel(cb, cLabel, 'Interpreter', 'none');
    colormap(ax, parula);
end
hold(ax, 'on');
if opts.showFit && numel(x) >= 2
    pp = polyfit(x, y, 1);
    xx = linspace(min(x), max(x), 100);
    yy = polyval(pp, xx);
    plot(ax, xx, yy, '-', 'LineWidth', 1.1, 'Color', [0.15 0.15 0.15]);
end
if opts.labelPoints && ismember('entity_id', T.Properties.VariableNames)
    ids = string(T.entity_id);
    for k = 1:numel(x)
        text(ax, x(k), y(k), [' ' char(ids(k))], 'Interpreter', 'none', 'FontSize', 7, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    end
end
xlabel(ax, abbreviateMetric(metricA), 'Interpreter', 'none');
ylabel(ax, abbreviateMetric(metricB), 'Interpreter', 'none');
box(ax, 'off');
grid(ax, 'on');
title(ax, sprintf('rho=%.3f, q=%.3g, n=%d', double(rowT.rho(1)), double(rowT.q_value(1)), getObsN(rowT)), ...
    'Interpreter', 'none', 'FontSize', 9);
end

function n = getObsN(T)
if ismember('n_obs', T.Properties.VariableNames)
    n = double(T.n_obs(1));
elseif ismember('n', T.Properties.VariableNames)
    n = double(T.n(1));
else
    n = NaN;
end
end

function out = abbreviateMetric(metricID)
id = char(string(metricID));
switch id
    case 'cell.frac_events_with_neighbor_recruitment'
        out = 'frac neighbor recruit';
    case 'cell.mean_neighbor_entry_delay_sec'
        out = 'neighbor delay (s)';
    case 'cell.mean_neighbor_partitions_per_event'
        out = 'neighbor parts/event';
    case 'cell.mean_neighbor_peak_overlap_frac_of_cell'
        out = 'neighbor peak frac';
    case 'cell.frac_propagative_events'
        out = 'frac propagative';
    case 'cell.mean_propagation_distance_um'
        out = 'prop dist (um)';
    case 'cell.mean_perimeter_um'
        out = 'perimeter (um)';
    otherwise
        out = strrep(id, 'cell.', '');
        out = strrep(out, 'slice.', '');
        out = strrep(out, '_', ' ');
end
end

function varName = findMetricVar(corrOut, metricID)
metricID = string(metricID);
ids = string(corrOut.metricIDs(:));
idx = find(ids == metricID, 1);
if isempty(idx)
    error('plot_CorrelationStats_top_pairs:MetricNotFound', ...
        'Metric %s not found in CorrelationStats output.', char(metricID));
end
varName = corrOut.metricNames{idx};
end

function [cvals, cLabel] = resolveDenominatorColor(corrOut, Tplot, metricA, metricB, colorOpt)
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
Dx = fetchDenominator(corrOut, Tplot, metricA);
Dy = fetchDenominator(corrOut, Tplot, metricB);
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

function [nr, nc] = bestGrid(n)
nc = ceil(sqrt(n));
nr = ceil(n / nc);
end
