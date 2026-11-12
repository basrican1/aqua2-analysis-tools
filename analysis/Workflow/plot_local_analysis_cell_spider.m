function out = plot_local_analysis_cell_spider(source, cellID, varargin)
%PLOT_LOCAL_ANALYSIS_CELL_SPIDER Event-level spider drill-down for one cell.
%
%   out = plot_local_analysis_cell_spider(result, cellID)
%   out = plot_local_analysis_cell_spider(result, cellID, ...
%       'metricFamilies', {'shape','kinetics','intensity','recruitment'}, ...
%       'dataView', 'recommended', ...
%       'combineEpochNames', {'ExperimentalEpoch1','ExperimentalEpoch2'}, ...
%       'combinedEpochLabel', 'CombinedPostEpochs')
%
% This is a cell-exploration figure, not a cross-metric absolute-comparison
% figure. Each spoke is scaled independently using the selected cell's event
% values across baseline and combined post-manipulation, so spoke radius is only
% interpretable within a metric, not across different metrics.
%
% The plot shows:
%   - baseline event distribution on each spoke (gray jittered points)
%   - combined post-manipulation event distribution on each spoke (red jittered points)
%   - baseline mean polygon
%   - post-manipulation mean polygon
%   - optional SEM ribbons per condition
%   - per-spoke baseline-vs-post significance stars
%
% Branch semantics are not used here; this figure operates directly on one
% cell's event-level distributions. With 'dataView'='recommended', canonical
% metrics use their preferred transform (often log2 fold for positive-valued
% magnitude metrics) while local event metrics currently use delta.
%
% Name-value:
%   'metricFamilies'     : event-level families, default all non-synchrony families
%   'metricSelection'    : 'default' | 'broad' | 'strict'
%   'dataView'           : 'raw' | 'norm' | 'recommended' (default)
%   'baselineEpochName'  : 'baseline' by default
%   'combineEpochNames'  : post epochs to combine, default from result if present
%   'combinedEpochLabel' : label for combined post epochs, default from result
%   'showBands'          : true (SEM ribbons)
%   'showPoints'         : true
%   'showScaleInLabels'  : true
%   'sigAlpha'           : 0.05
%   'jitterFrac'         : 0.18 angular jitter fraction of spoke spacing
%   'robustPercentiles'  : [5 95] robust per-spoke scaling range
%   'style'              : partial style override struct
%   'figure'             : optional target figure
%   'axes'               : optional target axes
%
% Output:
%   out.figure
%   out.axes
%   out.metricTable      one row per plotted metric with counts, scales, and p-values
%   out.cellID
%   out.fileID
%
% See also: plot_local_analysis_spider, plot_local_analysis_spider_mean,
%           CorrelationStats

result = load_local_analysis_result(source);

p = inputParser;
p.StructExpand = false;
addRequired(p, 'cellID', @(x) ischar(x) || isstring(x));
addParameter(p, 'metricFamilies', defaultCellSpiderFamilies(), @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'metricSelection', 'default', @(x) any(strcmpi(string(x), ["default","broad","strict"])));
addParameter(p, 'dataView', 'recommended', @(x) any(strcmpi(string(x), ["raw","norm","recommended"])));
addParameter(p, 'baselineEpochName', 'baseline', @(x) ischar(x) || isstring(x));
addParameter(p, 'combineEpochNames', inferDefaultCombineEpochs(result), @(x) iscell(x) || isstring(x));
addParameter(p, 'combinedEpochLabel', inferDefaultCombineLabel(result), @(x) ischar(x) || isstring(x));
addParameter(p, 'showBands', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'showPoints', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'showScaleInLabels', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'sigAlpha', 0.05, @isnumeric);
addParameter(p, 'jitterFrac', 0.18, @isnumeric);
addParameter(p, 'robustPercentiles', [5 95], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'style', struct(), @isstruct);
addParameter(p, 'figure', [], @(x) isempty(x) || ishghandle(x, 'figure'));
addParameter(p, 'axes', [], @(x) isempty(x) || ishghandle(x, 'axes'));
parse(p, cellID, varargin{:});

cellID = string(p.Results.cellID);
metricFamilies = normalizeFamilies(p.Results.metricFamilies);
selectionMode = lower(char(string(p.Results.metricSelection)));
dataView = lower(char(string(p.Results.dataView)));
baselineEpoch = string(p.Results.baselineEpochName);
combineEpochNames = string(p.Results.combineEpochNames(:));
combinedEpochLabel = string(p.Results.combinedEpochLabel);
style = mergeBWAPlotStyle(getBWAPlotStyle(), p.Results.style);

rows = collectCellSpiderRows(result, cellID, metricFamilies, selectionMode, dataView, baselineEpoch, combineEpochNames);
if isempty(rows)
    error('plot_local_analysis_cell_spider:NoMetrics', ...
        'No event-level metrics were available for cell %s with the requested settings.', cellID);
end
rows = enrichCellSpiderRows(rows, p.Results.robustPercentiles, p.Results.sigAlpha);
rows = rows(arrayfun(@(r) ~isempty(r.baseline) || ~isempty(r.post), rows));

metricTable = buildMetricSummaryTable(rows);
[baselineMean, postMean, baselineSem, postSem, labels, fileID] = buildScaledSpiderArrays(rows, p.Results.showScaleInLabels);
prov = resolveCellSpiderProvenance(result, rows, dataView);
headerText = buildHeaderText(prov, cellID, fileID, baselineEpoch, combinedEpochLabel);

[hFig, hAx] = resolvePlotAxesLocal(p.Results.figure, p.Results.axes, ['CellSpider_' char(cellID)]);
cla(hAx);
axes(hAx); %#ok<LAXES>
hold(hAx, 'on');
axis(hAx, 'equal');
axis(hAx, 'off');
set(hAx, 'Color', 'none');

nMetrics = numel(labels);
theta = pi/2 - 2*pi*(0:nMetrics-1)/nMetrics;
drawSpiderGrid(hAx, theta, style);

if p.Results.showBands
    drawBand(hAx, theta, baselineMean, baselineSem, [0.35 0.35 0.35], 0.14);
    drawBand(hAx, theta, postMean, postSem, style.colors.groupC, 0.12);
end

if p.Results.showPoints
    drawEventPoints(hAx, theta, rows, p.Results.jitterFrac, [0.35 0.35 0.35], style.colors.groupC, style);
end

hBase = plotSpiderPolygon(hAx, theta, baselineMean, [0.35 0.35 0.35], '--', 1.6);
hPost = plotSpiderPolygon(hAx, theta, postMean, style.colors.groupC, '-', 2.1);
drawMetricLabels(hAx, theta, labels, style);
drawSignificanceStars(hAx, theta, rows, style, p.Results.sigAlpha);

text(hAx, 0, 1.34, headerText, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontWeight', 'bold', 'FontSize', inferHeaderFontSize(style), ...
    'Color', style.colors.text, 'Interpreter', 'none');
text(hAx, 0, -1.31, 'Per-spoke scaling: robust range across baseline + combined post-CNO for this cell', ...
    'HorizontalAlignment', 'center', 'FontSize', max(style.fontSizeAxis-1,8), 'Color', [0.25 0.25 0.25], ...
    'Interpreter', 'none');
legend(hAx, [hBase hPost], {'baseline mean', char(combinedEpochLabel) + " mean"}, 'Location', 'southoutside', ...
    'Interpreter', 'none', 'FontSize', style.fontSizeLegend, 'Box', 'off');

out = struct();
out.figure = hFig;
out.axes = hAx;
out.metricTable = metricTable;
out.cellID = char(cellID);
out.fileID = char(fileID);
out.rows = rows;
out.normalizationProvenance = prov;
end

function families = normalizeFamilies(in)
if iscell(in)
    families = lower(string(in(:)))';
else
    families = lower(string(in));
end
families = families(strlength(families) > 0);
if isempty(families)
    families = defaultCellSpiderFamilies();
end
families(families == "all") = [];
if isempty(families)
    families = defaultCellSpiderFamilies();
end
end

function families = defaultCellSpiderFamilies()
families = ["shape","kinetics","intensity","propagation","containment","recruitment","soma"];
end

function combineEpochNames = inferDefaultCombineEpochs(result)
combineEpochNames = {"CNO1","CNO2"};
if isstruct(result) && isfield(result, 'spiderCombineEpochNames') && ~isempty(result.spiderCombineEpochNames)
    combineEpochNames = cellstr(string(result.spiderCombineEpochNames(:)));
end
end

function label = inferDefaultCombineLabel(result)
label = 'CNO_post';
if isstruct(result) && isfield(result, 'spiderCombinedEpochLabel') && strlength(string(result.spiderCombinedEpochLabel)) > 0
    label = char(string(result.spiderCombinedEpochLabel));
end
end

function rows = collectCellSpiderRows(result, cellID, metricFamilies, selectionMode, dataView, baselineEpoch, combineEpochNames)
rows = struct('metricID', {}, 'label', {}, 'fileID', {}, 'baseline', {}, 'post', {}, ...
    'conceptKey', {}, 'scaleLo', {}, 'scaleHi', {}, 'pValue', {}, 'isSig', {}, 'sigMethod', {});
policy = correlationMetricPolicy('event', selectionMode);
excludeIDs = string(policy.excludeIDs(:));
excludeIDs = setdiff(excludeIDs, ["event.cell_peak_overlap_frac_of_cell", "event.cell_peak_overlap_frac_of_event_frame"], 'stable');
canonicalConcepts = strings(0,1);

cRows = collectCanonicalRows(result, cellID, dataView, baselineEpoch, combineEpochNames, excludeIDs, metricFamilies);
rows = [rows, cRows]; %#ok<AGROW>
canonicalConcepts = string(arrayfun(@(r) r.conceptKey, cRows, 'UniformOutput', false));
eRows = collectEventRows(result, cellID, dataView, baselineEpoch, combineEpochNames, excludeIDs, metricFamilies, canonicalConcepts);
rows = [rows, eRows]; %#ok<AGROW>
rows = rows(arrayfun(@(r) ~isempty(r.baseline) || ~isempty(r.post), rows));
end

function rows = collectCanonicalRows(result, cellID, dataView, baselineEpoch, combineEpochNames, excludeIDs, metricFamilies)
rows = struct('metricID', {}, 'label', {}, 'fileID', {}, 'baseline', {}, 'post', {}, ...
    'conceptKey', {}, 'scaleLo', {}, 'scaleHi', {}, 'pValue', {}, 'isSig', {}, 'sigMethod', {});
if ~isfield(result, 'canonical') || ~isfield(result.canonical, 'metrics') || ~isstruct(result.canonical.metrics)
    return
end
mNames = fieldnames(result.canonical.metrics);
mNames = filterPreferredCanonicalMetricIDs(mNames);
for i = 1:numel(mNames)
    metricID = string(mNames{i});
    if any(metricID == excludeIDs)
        continue
    end
    if ~ismember(bwaResolveMetricFamily(metricID), metricFamilies)
        continue
    end
    mOut = selectCanonicalMetricForCellSpider(result, metricID, dataView);
    if ~isfield(mOut, 'eventValues') || ~istable(mOut.eventValues) || isempty(mOut.eventValues)
        continue
    end
    T = mOut.eventValues;
    [b, p, fileID] = extractCellEpochValues(T, cellID, baselineEpoch, combineEpochNames, metricID);
    if isempty(b) && isempty(p)
        continue
    end
    rows(end+1) = makeRow(metricID, fileID, b, p); %#ok<AGROW>
end
end

function rows = collectEventRows(result, cellID, dataView, baselineEpoch, combineEpochNames, excludeIDs, metricFamilies, canonicalConcepts)
rows = struct('metricID', {}, 'label', {}, 'fileID', {}, 'baseline', {}, 'post', {}, ...
    'conceptKey', {}, 'scaleLo', {}, 'scaleHi', {}, 'pValue', {}, 'isSig', {}, 'sigMethod', {});
metricContainer = selectMetricContainer(result, dataView, 'delta');
if ~isfield(metricContainer, 'tables') || ~isstruct(metricContainer.tables)
    return
end
for tname = fieldnames(metricContainer.tables)'
    tblName = tname{1};
    if ~startsWith(tblName, 'event_')
        continue
    end
    T = metricContainer.tables.(tblName);
    if ~istable(T) || isempty(T) || ~ismember('event_uid', T.Properties.VariableNames)
        continue
    end
    valueCols = pickEventLocalValueColumnsForCellSpider(T, dataView);
    for j = 1:numel(valueCols)
        vcol = valueCols{j};
        metricID = "event." + string(stripEventLocalColumnPrefix(vcol));
        if any(metricID == excludeIDs)
            continue
        end
        if ~ismember(bwaResolveMetricFamily(metricID), metricFamilies)
            continue
        end
        ck = metricConceptKey(metricID);
        if any(ck == canonicalConcepts)
            continue
        end
        [b, p, fileID] = extractCellEpochValues(T, cellID, baselineEpoch, combineEpochNames, metricID, vcol);
        if isempty(b) && isempty(p)
            continue
        end
        rows(end+1) = makeRow(metricID, fileID, b, p); %#ok<AGROW>
    end
end
end

function metricContainer = selectMetricContainer(result, dataView, preferredTransform)
metricContainer = result.metrics;
if strcmp(dataView, 'recommended') || strcmp(dataView, 'norm')
    field = matlab.lang.makeValidName(preferredTransform);
    if isfield(result, 'normalized') && isfield(result.normalized, 'metrics') && isfield(result.normalized.metrics, field)
        metricContainer = result.normalized.metrics.(field);
    elseif isfield(result, 'metrics_by_transform') && isfield(result.metrics_by_transform, field)
        metricContainer = result.metrics_by_transform.(field);
    end
end
end

function mOut = selectCanonicalMetricForCellSpider(result, metricID, dataView)
mOut = result.canonical.metrics.(char(metricID));
if strcmp(dataView, 'raw')
    return
end
if strcmp(dataView, 'norm')
    if isfield(result, 'canonical_norm') && isfield(result.canonical_norm, 'metrics') && isfield(result.canonical_norm.metrics, char(metricID))
        mOut = result.canonical_norm.metrics.(char(metricID));
    end
    return
end
tr = bwaPreferredTransformForMetricID(metricID, 'canonical');
field = matlab.lang.makeValidName(tr);
if isfield(result, 'normalized') && isfield(result.normalized, 'canonical') && isfield(result.normalized.canonical, field) && ...
        isfield(result.normalized.canonical.(field), 'metrics') && isfield(result.normalized.canonical.(field).metrics, char(metricID))
    mOut = result.normalized.canonical.(field).metrics.(char(metricID));
elseif isfield(result, 'canonical_norms') && isfield(result.canonical_norms, field) && ...
        isfield(result.canonical_norms.(field), 'metrics') && isfield(result.canonical_norms.(field).metrics, char(metricID))
    mOut = result.canonical_norms.(field).metrics.(char(metricID));
end
end

function [baselineVals, postVals, fileID] = extractCellEpochValues(T, cellID, baselineEpoch, combineEpochNames, metricID, varargin)
baselineVals = [];
postVals = [];
fileID = "";
if nargin > 5
    valCol = varargin{1};
else
    valCol = inferValueColumnFromTable(T, metricID);
end
if isempty(valCol) || ~ismember(valCol, T.Properties.VariableNames)
    return
end
cellCol = resolveCellColumn(T);
if isempty(cellCol)
    return
end
T = T(string(T.(cellCol)) == string(cellID), :);
if isempty(T)
    return
end
if ismember('file_id', T.Properties.VariableNames)
    fileID = string(T.file_id(1));
end
if ~ismember('epoch_name', T.Properties.VariableNames)
    return
end
Tb = T(string(T.epoch_name) == baselineEpoch, :);
Tp = T(ismember(string(T.epoch_name), combineEpochNames), :);
if isempty(Tb) && isempty(Tp)
    return
end
baselineVals = aggregateEventRows(Tb, valCol);
postVals = aggregateEventRows(Tp, valCol);
end

function vals = aggregateEventRows(T, valCol)
vals = [];
if isempty(T) || ~ismember(valCol, T.Properties.VariableNames)
    return
end
if ~ismember('event_uid', T.Properties.VariableNames)
    x = double(T.(valCol));
    vals = x(isfinite(x));
    return
end
[G,~] = findgroups(string(T.event_uid));
vals = splitapply(@meanOmitNaN, double(T.(valCol)), G);
vals = vals(isfinite(vals));
end

function row = makeRow(metricID, fileID, baselineVals, postVals)
row = struct();
row.metricID = char(metricID);
row.label = formatCellSpiderLabel(metricID);
row.fileID = char(fileID);
row.baseline = baselineVals(:)';
row.post = postVals(:)';
row.conceptKey = char(metricConceptKey(metricID));
row.scaleLo = NaN;
row.scaleHi = NaN;
row.pValue = NaN;
row.isSig = false;
row.sigMethod = '';
end

function rows = enrichCellSpiderRows(rows, prc, sigAlpha)
for i = 1:numel(rows)
    allVals = [rows(i).baseline(:); rows(i).post(:)];
    allVals = allVals(isfinite(allVals));
    if isempty(allVals)
        continue
    end
    lo = prctile(allVals, prc(1));
    hi = prctile(allVals, prc(2));
    if ~isfinite(lo), lo = min(allVals); end
    if ~isfinite(hi), hi = max(allVals); end
    if hi <= lo
        lo = min(allVals);
        hi = max(allVals);
    end
    rows(i).scaleLo = lo;
    rows(i).scaleHi = hi;
    [pVal, method] = computeRowPValue(rows(i).baseline, rows(i).post);
    rows(i).pValue = pVal;
    rows(i).sigMethod = method;
    rows(i).isSig = isfinite(pVal) && pVal < sigAlpha;
end
end

function [pVal, method] = computeRowPValue(baselineVals, postVals)
pVal = NaN;
method = '';
b = double(baselineVals(:));
p = double(postVals(:));
b = b(isfinite(b));
p = p(isfinite(p));
if numel(b) < 2 || numel(p) < 2
    return
end
try
    pVal = ranksum(b, p);
    method = 'ranksum';
catch
    try
        [~, pVal] = ttest2(b, p, 'Vartype', 'unequal');
        method = 'ttest2';
    catch
        pVal = NaN;
        method = '';
    end
end
end

function tbl = buildMetricSummaryTable(rows)
n = numel(rows);
tbl = table('Size', [n, 11], ...
    'VariableTypes', {'string','string','double','double','double','double','double','double','double','double','string'}, ...
    'VariableNames', {'metric_id','label','n_baseline','n_post','mean_baseline','mean_post','delta_mean','scale_lo','scale_hi','p_value','sig_method'});
for i = 1:n
    tbl.metric_id(i) = string(rows(i).metricID);
    tbl.label(i) = string(rows(i).label);
    tbl.n_baseline(i) = numel(rows(i).baseline);
    tbl.n_post(i) = numel(rows(i).post);
    tbl.mean_baseline(i) = meanOmitNaN(rows(i).baseline);
    tbl.mean_post(i) = meanOmitNaN(rows(i).post);
    tbl.delta_mean(i) = tbl.mean_post(i) - tbl.mean_baseline(i);
    tbl.scale_lo(i) = rows(i).scaleLo;
    tbl.scale_hi(i) = rows(i).scaleHi;
    tbl.p_value(i) = rows(i).pValue;
    tbl.sig_method(i) = string(rows(i).sigMethod);
end
end

function [mB, mP, sB, sP, labels, fileID] = buildScaledSpiderArrays(rows, showScaleInLabels)
n = numel(rows);
mB = nan(1,n); mP = nan(1,n); sB = nan(1,n); sP = nan(1,n);
labels = cell(1,n);
fileID = "";
for i = 1:n
    labels{i} = decorateLabel(rows(i), showScaleInLabels);
    if strlength(fileID) == 0 && strlength(string(rows(i).fileID)) > 0
        fileID = string(rows(i).fileID);
    end
    [bScaled, pScaled] = scaleRowValues(rows(i));
    mB(i) = meanOmitNaN(bScaled);
    mP(i) = meanOmitNaN(pScaled);
    sB(i) = semOmitNaN(bScaled);
    sP(i) = semOmitNaN(pScaled);
end
end

function label = decorateLabel(row, showScaleInLabels)
label = row.label;
if ~showScaleInLabels || ~isfinite(row.scaleLo) || ~isfinite(row.scaleHi)
    return
end
label = sprintf('%s\n[%s, %s]', label, formatScaleNumber(row.scaleLo), formatScaleNumber(row.scaleHi));
end

function txt = formatScaleNumber(x)
if ~isfinite(x)
    txt = 'NaN';
elseif abs(x) >= 100
    txt = sprintf('%.0f', x);
elseif abs(x) >= 10
    txt = sprintf('%.1f', x);
else
    txt = sprintf('%.2f', x);
end
end

function [hFig, hAx] = resolvePlotAxesLocal(hFig, hAx, figName)
if isempty(hFig) || ~ishghandle(hFig, 'figure')
    hFig = figure('Name', figName, 'Color', 'w', 'Units', 'pixels', 'Position', [80 80 1100 900]);
else
    figure(hFig);
end
if isempty(hAx) || ~ishghandle(hAx, 'axes')
    hAx = axes('Parent', hFig, 'Position', [0.08 0.14 0.84 0.70]);
end
end

function fs = inferHeaderFontSize(style)
fs = 11;
if isstruct(style)
    if isfield(style, 'fontSizeTitle') && ~isempty(style.fontSizeTitle)
        fs = max(style.fontSizeTitle - 1, 10);
        return
    end
    if isfield(style, 'fontSizeAxis') && ~isempty(style.fontSizeAxis)
        fs = max(style.fontSizeAxis + 1, 10);
        return
    end
end
end

function drawSpiderGrid(hAx, theta, style)
rings = [0.25 0.5 0.75 1.0];
for r = rings
    t = linspace(0, 2*pi, 300);
    plot(hAx, r*cos(t), r*sin(t), 'Color', [0.88 0.88 0.88], 'LineWidth', 0.8);
end
for j = 1:numel(theta)
    plot(hAx, [0 cos(theta(j))], [0 sin(theta(j))], 'Color', [0.84 0.84 0.84], 'LineWidth', 0.8);
end
text(hAx, 0, 0, '0', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', max(style.fontSizeAxis-2,8), 'Color', [0.4 0.4 0.4]);
end

function drawBand(hAx, theta, mu, se, color, alphaVal)
if all(~isfinite(mu))
    return
end
lo = max(mu - se, 0);
hi = min(mu + se, 1);
[xh,yh] = pol2cartWrap(theta, hi);
[xl,yl] = pol2cartWrap(theta, lo);
patch('Parent', hAx, 'XData', [xh fliplr(xl)], 'YData', [yh fliplr(yl)], ...
    'FaceColor', color, 'FaceAlpha', alphaVal, 'EdgeColor', 'none');
end

function drawEventPoints(hAx, theta, rows, jitterFrac, colorB, colorP, style)
n = numel(rows);
if n == 0
    return
end
step = 2*pi / n;
for j = 1:n
    [bScaled, pScaled] = scaleRowValues(rows(j));
    if ~isempty(bScaled)
        jt = theta(j) + (rand(size(bScaled)) - 0.5) * step * jitterFrac - step*0.12;
        [x,y] = pol2cart(jt, bScaled(:)');
        scatter(hAx, x, y, style.scatterSize*5, 'MarkerFaceColor', colorB, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.35);
    end
    if ~isempty(pScaled)
        jt = theta(j) + (rand(size(pScaled)) - 0.5) * step * jitterFrac + step*0.12;
        [x,y] = pol2cart(jt, pScaled(:)');
        scatter(hAx, x, y, style.scatterSize*5, 'MarkerFaceColor', colorP, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.30);
    end
end
end

function h = plotSpiderPolygon(hAx, theta, mu, color, lineStyle, lw)
[x,y] = pol2cartWrap(theta, mu);
h = patch('Parent', hAx, 'XData', x, 'YData', y, 'FaceColor', color, 'FaceAlpha', 0.10, ...
    'EdgeColor', color, 'LineStyle', lineStyle, 'LineWidth', lw);
end

function drawMetricLabels(hAx, theta, labels, style)
for j = 1:numel(labels)
    [x,y] = pol2cart(theta(j), 1.17);
    text(hAx, x, y, labels{j}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Interpreter', 'none', 'FontSize', max(style.fontSizeAxis-1,8), 'Color', style.colors.text);
end
end

function drawSignificanceStars(hAx, theta, rows, style, ~)
for j = 1:numel(rows)
    if ~rows(j).isSig
        continue
    end
    [x,y] = pol2cart(theta(j), 1.05);
    text(hAx, x, y, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Interpreter', 'none', 'FontSize', max(style.fontSizeAxis+8, 16), 'FontWeight', 'bold', ...
        'Color', style.colors.groupC);
end
end

function txt = buildHeaderText(prov, cellID, fileID, baselineEpoch, postLabel)
cellLabel = compactCellID(cellID);
fileLabel = bwaFormatFileIDLabel(char(fileID));
if strlength(string(fileLabel)) == 0
    fileLabel = char(string(fileID));
end
normLabel = string(prov.compactLabel);
txt = sprintf('%s | %s | %s vs %s | %s', cellLabel, fileLabel, char(string(baselineEpoch)), char(string(postLabel)), normLabel);
end

function out = compactCellID(cellID)
out = char(string(cellID));
tok = regexp(out, '(_C\d+)$', 'tokens', 'once');
if ~isempty(tok)
    out = tok{1}(2:end);
end
end

function prov = resolveCellSpiderProvenance(result, rows, dataView)
metricIDs = strings(0,1);
unitLevels = strings(0,1);
for i = 1:numel(rows)
    metricIDs(end+1,1) = string(rows(i).metricID); %#ok<AGROW>
    if startsWith(metricIDs(end), "event.")
        unitLevels(end+1,1) = "event"; %#ok<AGROW>
    else
        unitLevels(end+1,1) = "canonical"; %#ok<AGROW>
    end
end
prov = bwaResolveNormalizationProvenance(result, dataView, metricIDs, unitLevels);
end

function [x,y] = pol2cartWrap(theta, r)
theta = theta(:)';
r = r(:)';
[x,y] = pol2cart(theta, r);
x(end+1) = x(1);
y(end+1) = y(1);
end

function [bScaled, pScaled] = scaleRowValues(row)
allVals = [row.baseline(:); row.post(:)];
allVals = allVals(isfinite(allVals));
if isempty(allVals)
    bScaled = [];
    pScaled = [];
    return
end
lo = row.scaleLo;
hi = row.scaleHi;
if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
    lo = min(allVals);
    hi = max(allVals);
end
if hi <= lo
    bScaled = 0.5 * ones(size(row.baseline));
    pScaled = 0.5 * ones(size(row.post));
else
    bScaled = min(max((row.baseline - lo) ./ (hi - lo), 0), 1);
    pScaled = min(max((row.post - lo) ./ (hi - lo), 0), 1);
end
bScaled = bScaled(isfinite(bScaled));
pScaled = pScaled(isfinite(pScaled));
end

function v = meanOmitNaN(x)
x = double(x(:));
x = x(isfinite(x));
if isempty(x)
    v = NaN;
else
    v = mean(x);
end
end

function v = semOmitNaN(x)
x = double(x(:));
x = x(isfinite(x));
if numel(x) < 2
    v = 0;
else
    v = std(x, 0) / sqrt(numel(x));
end
end

function col = resolveCellColumn(T)
for cand = {'primary_cell_uid','cell_uid'}
    if ismember(cand{1}, T.Properties.VariableNames)
        col = cand{1};
        return
    end
end
col = '';
end

function valueCol = inferValueColumnFromTable(T, metricID)
valueCol = '';
mid = char(string(metricID));
mid = regexprep(mid, '^event\.', '');
cands = [mid + "_norm", mid, "value"];
for i = 1:numel(cands)
    if ismember(char(cands(i)), T.Properties.VariableNames)
        valueCol = char(cands(i));
        return
    end
end
numCols = T.Properties.VariableNames(varfun(@(x) isnumeric(x) || islogical(x), T, 'OutputFormat', 'uniform'));
numCols = setdiff(numCols, {'epoch_idx'});
if ~isempty(numCols)
    valueCol = numCols{1};
end
end

function valueCols = pickEventLocalValueColumnsForCellSpider(T, dataView)
valueCols = {};
vars = string(T.Properties.VariableNames);
skip = ["event_uid","animal_id","slice_num","file_id","primary_cell_uid","partition_cell_uid","epoch_idx","epoch_name","cell_uid", ...
    "frame_rate_s","spatial_res","n_frames","is_primary_cell"];
preferredPairs = [ ...
    "soma_centroid_dist_px","soma_centroid_dist_um"; ...
    "soma_start_dist_px","soma_start_dist_um"; ...
    "soma_end_dist_px","soma_end_dist_um"; ...
    "soma_net_distance_change_px","soma_net_distance_change_um"; ...
    "perimeter_px","perimeter_um"; ...
    "propagation_distance_px","propagation_distance_um"; ...
    "outside_cell_area_npix","outside_cell_area_um2"; ...
    "cell_entry_delay_from_primary_frames","cell_entry_delay_from_primary_sec"; ...
    "cell_overlap_duration_frames","cell_overlap_duration_sec"];
complementBase = ["cell_overlap_frac_of_event", "frac_outside_cell"];
redundantBase = ["soma_start_dist_um","soma_end_dist_um","soma_start_dist_px","soma_end_dist_px", ...
    "perimeter_um","perimeter_px","circularity","propagation_distance_um","propagation_distance_px"];
for i = 1:numel(vars)
    vn = vars(i);
    if any(vn == skip) || contains(vn, "threshold")
        continue
    end
    base = string(regexprep(char(vn), '_norm$', ''));
    if any(base == redundantBase)
        continue
    end
    if strcmp(dataView, 'raw')
        if endsWith(vn, "_norm")
            continue
        end
        if ismember(vn, complementBase) && any(ismember(complementBase, vars)) && vn == "cell_overlap_frac_of_event"
            continue
        end
        if (endsWith(vn, "_px") || endsWith(vn, "_frames") || endsWith(vn, "_npix")) && any(preferredPairs(:,1) == vn & ismember(preferredPairs(:,2), vars))
            continue
        end
    else
        if ~endsWith(vn, "_norm")
            if any(ismember(vars, vn + "_norm"))
                continue
            end
            continue
        end
        if ismember(base, complementBase) && any(ismember(vars, complementBase + "_norm")) && base == "cell_overlap_frac_of_event"
            continue
        end
        if any(preferredPairs(:,1) == base & ismember(preferredPairs(:,2) + "_norm", vars))
            continue
        end
    end
    x = T.(char(vn));
    if isnumeric(x) || islogical(x)
        valueCols{end+1} = char(vn); %#ok<AGROW>
    end
end
valueCols = unique(valueCols, 'stable');
end

function out = stripEventLocalColumnPrefix(vcol)
out = char(string(vcol));
out = regexprep(out, '_norm$', '');
end

function ids = filterPreferredCanonicalMetricIDs(ids)
ids = cellstr(string(ids(:)));
if isempty(ids), return; end
skip = strings(0,1);
pairs = { ...
    'area_px2','area_um2'; ...
    'duration_frames','duration_sec'; ...
    'rise19_frames','rise19_sec'; ...
    'fall91_frames','fall91_sec'; ...
    'decay_tau_frames','decay_tau_sec'; ...
    'dff_auc','dff_auc_sec'; ...
    'perimeter_px','perimeter_um'; ...
    'propagation_distance_px','propagation_distance_um'; ...
    'max_prop_speed','max_prop_speed_um_per_sec'};
for i = 1:size(pairs,1)
    raw = pairs{i,1}; cal = pairs{i,2};
    if any(strcmp(ids, raw)) && any(strcmp(ids, cal))
        skip(end+1,1) = string(raw); %#ok<AGROW>
    end
end
ids = ids(~ismember(string(ids), skip));
end

function key = metricConceptKey(metricID)
id = lower(char(string(metricID)));
id = regexprep(id, '^(cell|slice|event)\.', '');
switch id
    case {'area_um2','area_px2'}
        key = "area";
    case {'duration_sec','duration_frames'}
        key = "duration";
    case {'rise19_sec','rise19_frames'}
        key = "rise19";
    case {'fall91_sec','fall91_frames'}
        key = "fall91";
    case {'decay_tau_sec','decay_tau_frames'}
        key = "decay_tau";
    case {'dff_auc_sec','dff_auc'}
        key = "dff_auc";
    case {'dff_max2'}
        key = "dff_max2";
    case {'perimeter_um','perimeter_px'}
        key = "perimeter";
    case {'circularity'}
        key = "circularity";
    case {'propagation_distance_um','propagation_distance_px'}
        key = "propagation_distance";
    case {'max_prop_speed_um_per_sec','max_prop_speed'}
        key = "max_prop_speed";
    case {'soma_centroid_dist_um','soma_centroid_dist_px'}
        key = "soma_centroid_dist";
    case {'soma_net_distance_change_um','soma_net_distance_change_px'}
        key = "soma_net_distance_change";
    case {'spread_count'}
        key = "spread_count";
    case {'is_multicell'}
        key = "is_multicell";
    case {'frac_outside_cell','outside_cell_area_um2','outside_cell_area_npix'}
        key = "outside_cell";
    case {'cell_entry_delay_from_primary_sec','cell_entry_delay_from_primary_frames'}
        key = "neighbor_delay";
    case {'cell_peak_overlap_frac_of_cell'}
        key = "neighbor_peak";
    case {'cell_peak_overlap_frac_of_event_frame'}
        key = "neighbor_peak_event";
    case {'cell_overlap_duration_sec','cell_overlap_duration_frames'}
        key = "overlap_duration";
    otherwise
        key = string(id);
end
end

function label = formatCellSpiderLabel(metricID)
id = lower(char(string(metricID)));
id = regexprep(id, '^(cell|slice|event)[._]', '');
switch id
    case 'frequency'
        label = sprintf('freq\n(events/min)');
    case 'area_um2'
        label = sprintf('area\n(um^2)');
    case 'duration_sec'
        label = sprintf('duration\n(s)');
    case 'rise19_sec'
        label = sprintf('rise19\n(s)');
    case 'fall91_sec'
        label = sprintf('fall91\n(s)');
    case 'decay_tau_sec'
        label = sprintf('decay tau\n(s)');
    case 'dff_auc_sec'
        label = sprintf('dF/F AUC\n((dF/F)*s)');
    case 'dff_max2'
        label = 'dF/F max';
    case 'perimeter_um'
        label = sprintf('perimeter\n(um)');
    case 'circularity'
        label = 'circularity';
    case 'propagation_distance_um'
        label = sprintf('prop dist\n(um)');
    case 'max_prop_speed_um_per_sec'
        label = sprintf('max prop spd\n(um/s)');
    case 'spread_count'
        label = 'cells/event';
    case 'is_multicell'
        label = 'multicell';
    case 'soma_centroid_dist_um'
        label = sprintf('centroid-soma d\n(um)');
    case 'soma_net_distance_change_um'
        label = sprintf('net soma Δd\n(um)');
    case 'frac_outside_cell'
        label = 'outside frac';
    case 'outside_cell_area_um2'
        label = sprintf('outside area\n(um^2)');
    case 'cell_entry_delay_from_primary_sec'
        label = sprintf('nbr delay\n(s)');
    case 'cell_peak_overlap_frac_of_cell'
        label = 'nbr peak';
    case 'cell_peak_overlap_frac_of_event_frame'
        label = 'peak of event';
    case 'cell_overlap_duration_sec'
        label = sprintf('overlap dur\n(s)');
    otherwise
        label = strrep(id, '_', ' ');
end
end
