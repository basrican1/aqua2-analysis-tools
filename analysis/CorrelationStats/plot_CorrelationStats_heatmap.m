function out = plot_CorrelationStats_heatmap(source, varargin)
%PLOT_CORRELATIONSTATS_HEATMAP Plot a correlation heatmap from CorrelationStats output.
%
% Usage:
%   out = plot_CorrelationStats_heatmap(corrOut)
%   out = plot_CorrelationStats_heatmap(result, 'unitLevel', 'slice', 'metricFamilies', {'shape','synchrony'})
%
% Name-value:
%   'cluster'         : true/false, reorder metrics by average-linkage clustering (default true)
%   'showLabels'      : true/false, show metric labels on both axes (default true)
%   'showValues'      : true/false, overlay numeric rho values (default false)
%   'minAbsR'         : optional filter; drop metrics with all |rho| below this threshold (default 0)
%   'figureName'      : figure name (default 'Correlation Heatmap')
%   'colorLimits'     : 1x2 numeric, default [-1 1]
%   plus any CorrelationStats name-value arguments when source is not already a CorrelationStats output

if isstruct(source) && isfield(source, 'source') && strcmp(source.source, 'CorrelationStats')
    corrOut = source;
else
    corrOut = CorrelationStats(source, varargin{:});
    varargin = {};
end

p = inputParser;
p.StructExpand = false;
addParameter(p, 'cluster', true, @islogical);
addParameter(p, 'showLabels', true, @islogical);
addParameter(p, 'showValues', false, @islogical);
addParameter(p, 'showGroupBars', true, @islogical);
addParameter(p, 'showRightGroupBar', false, @islogical); % ignored; retained for compatibility
addParameter(p, 'minAbsR', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'figureName', 'Correlation Heatmap', @(x) ischar(x) || isstring(x));
addParameter(p, 'colorLimits', [-1 1], @(x) isnumeric(x) && numel(x) == 2);
parse(p, varargin{:});

R = double(corrOut.R);
metricIDs = string(corrOut.metricIDs(:));
if isempty(R) || isempty(metricIDs)
    error('plot_CorrelationStats_heatmap:EmptyCorrelation', 'CorrelationStats output is empty.');
end

keep = true(size(metricIDs));
if p.Results.minAbsR > 0
    for i = 1:numel(metricIDs)
        r = R(i,:);
        r(i) = NaN;
        keep(i) = any(abs(r(isfinite(r))) >= p.Results.minAbsR);
    end
end
R = R(keep, keep);
metricIDs = metricIDs(keep);
if isempty(metricIDs)
    error('plot_CorrelationStats_heatmap:NoMetricsAfterFilter', 'No metrics remain after minAbsR filtering.');
end

order = 1:numel(metricIDs);
if p.Results.cluster && numel(metricIDs) > 2
    order = computeMetricOrder(R);
end
R = R(order, order);
metricIDs = metricIDs(order);
labels = cellstr(string(metricIDs));
labels = abbreviateCorrelationMetricLabels(labels);
groups = classifyCorrelationMetricGroups(metricIDs);

hFig = figure('Name', char(string(p.Results.figureName)), 'Color', 'w', 'Visible', 'on');
ax = axes('Parent', hFig);
imagesc(ax, R, p.Results.colorLimits);
axis(ax, 'image');
colormap(ax, parula(256));
cb = colorbar(ax);
cb.Label.String = sprintf('%s correlation', capitalizeMethod(corrOut.method));
mainTitle = annotation(hFig, 'textbox', [0.08 0.95 0.82 0.04], ...
    'String', buildTitle(corrOut), 'Interpreter', 'none', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold', ...
    'Color', [0.15 0.15 0.15], 'FitBoxToText', 'off');

if p.Results.showLabels
    set(ax, 'XTick', 1:numel(metricIDs), 'YTick', 1:numel(metricIDs), ...
        'XTickLabel', labels, 'YTickLabel', labels, 'TickLabelInterpreter', 'none');
    xtickangle(ax, 45);
else
    set(ax, 'XTick', [], 'YTick', []);
end

if p.Results.showValues
    addValueText(ax, R);
end

groupHandles = struct();
if p.Results.showGroupBars
    groupHandles = addGroupBoxes(hFig, ax, groups);
    updateGroupBoxLayout(hFig, ax, groupHandles);
    hFig.SizeChangedFcn = @(src, evt) updateGroupBoxLayout(src, ax, groupHandles); %#ok<NASGU>
end
figure(hFig);

out = struct();
out.figure = hFig;
out.axes = ax;
out.colorbar = cb;
out.titleAnnotation = mainTitle;
out.groupHandles = groupHandles;
out.metricIDs = cellstr(metricIDs(:))';
out.labels = labels;
out.groups = cellstr(groups(:))';
out.R = R;
out.order = order(:)';
end

function order = computeMetricOrder(R)
D = 1 - abs(R);
D(~isfinite(D)) = 1;
D(1:size(D,1)+1:end) = 0;
% Numerical cleanup for clustering helpers that require a proper distance form.
D = 0.5 * (D + D.');
D(D < 0) = 0;
Y = squareformSymmetric(D);
if isempty(Y)
    order = 1:size(R,1);
    return
end
Z = linkage(Y, 'average');
try
    order = optimalleaforder(Z, Y);
catch
    try
        order = optimalleaforder(Z, D);
    catch
        order = 1:size(R,1);
    end
end
end

function y = squareformSymmetric(D)
n = size(D,1);
if n < 2
    y = [];
    return
end
y = zeros(n*(n-1)/2,1);
k = 0;
for i = 1:n-1
    for j = i+1:n
        k = k + 1;
        y(k) = D(i,j);
    end
end
end

function labels = abbreviateCorrelationMetricLabels(labels)
for i = 1:numel(labels)
    labels{i} = abbreviateOne(labels{i});
end
end

function out = abbreviateOne(in)
key = lower(char(string(in)));
key = strrep(key, 'cell.', '');
key = strrep(key, 'slice.', '');
switch key
    case 'area_um2'
        out = 'area';
    case 'duration_sec'
        out = 'duration';
    case 'rise19_sec'
        out = 'rise19';
    case 'fall91_sec'
        out = 'fall91';
    case 'decay_tau_sec'
        out = 'decay tau';
    case 'dff_auc_sec'
        out = 'dF/F AUC';
    case 'dff_max2'
        out = 'dF/F max';
    case 'perimeter_um'
        out = 'perimeter';
    case 'propagation_distance_um'
        out = 'prop dist';
    case 'max_prop_speed_um_per_sec'
        out = 'max prop spd';
    case 'frequency'
        out = 'freq';
    case 'frac_events_with_neighbor_recruitment'
        out = 'cell frac recruited';
    case 'mean_neighbor_entry_delay_sec'
        out = 'nbr delay';
    case 'mean_neighbor_peak_overlap_frac_of_cell'
        out = 'nbr peak frac';
    case 'mean_neighbor_partitions_per_event'
        out = 'cell nbr parts/event';
    case 'frac_primary_events_with_neighbor_recruitment'
        out = 'slice frac recruited';
    case 'frac_active_cells'
        out = 'frac active';
    case 'coincidence_rate'
        out = 'coincidence';
    case 'burst_fraction'
        out = 'burst frac';
    case 'mean_pairwise_cell_correlation'
        out = 'pair corr';
    otherwise
        out = strrep(key, '_', ' ');
    end
end

function addValueText(ax, R)
n = size(R,1);
for i = 1:n
    for j = 1:n
        if ~isfinite(R(i,j))
            continue
        end
        c = 'k';
        if abs(R(i,j)) > 0.55
            c = 'w';
        end
        text(ax, j, i, sprintf('%.2f', R(i,j)), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Color', c, 'FontSize', 8, 'Interpreter', 'none');
    end
end
end

function groups = classifyCorrelationMetricGroups(metricIDs)
groups = strings(numel(metricIDs),1);
for i = 1:numel(metricIDs)
    id = lower(char(string(metricIDs(i))));
    if startsWith(id, 'cell.')
        if contains(id, 'neighbor') || contains(id, 'recruitment')
            groups(i) = "recruitment";
        else
            groups(i) = "cell";
        end
    elseif startsWith(id, 'slice.')
        if contains(id, 'coincidence') || contains(id, 'burst') || contains(id, 'pairwise') || contains(id, 'active_cells')
            groups(i) = "network";
        elseif contains(id, 'neighbor') || contains(id, 'recruitment')
            groups(i) = "recruitment";
        else
            groups(i) = "slice";
        end
    elseif startsWith(id, 'event.')
        groups(i) = "event";
    elseif startsWith(id, 'event_local.')
        groups(i) = "event";
    elseif startsWith(id, 'derived_')
        groups(i) = "derived";
    elseif any(strcmp(id, {'footprint','propagation.distance','propagation.angle','centroid.distcomp'}))
        groups(i) = "derived";
    else
        groups(i) = "canonical";
    end
end
end

function H = addGroupBoxes(hFig, ax, groups)
groups = string(groups(:));
H = struct();
if isempty(groups)
    return
end
palette = groupPalette();
[runs, labels] = contiguousGroupRuns(groups);
H.boxes = gobjects(0,1);
H.text = gobjects(0,1);
H.runs = runs;
H.labels = labels;
for i = 1:size(runs,1)
    key = char(labels(i));
    color = [0.5 0.5 0.5];
    if isfield(palette, matlab.lang.makeValidName(key))
        color = palette.(matlab.lang.makeValidName(key));
    end
    H.boxes(end+1,1) = annotation(hFig, 'rectangle', [0 0 0.1 0.1], ...
        'FaceColor', 'none', 'Color', color, 'LineWidth', 1.6); %#ok<AGROW>
    H.text(end+1,1) = annotation(hFig, 'textbox', [0 0 0.1 0.04], ...
        'String', char(labels(i)), 'Interpreter', 'none', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'Color', color, 'FitBoxToText', 'off'); %#ok<AGROW>
end
end

function palette = groupPalette()
palette = struct();
palette.canonical = [0.16 0.38 0.85];
palette.derived = [0.70 0.45 0.12];
palette.cell = [0.05 0.55 0.48];
palette.recruitment = [0.85 0.48 0.08];
palette.network = [0.78 0.16 0.20];
palette.event = [0.47 0.30 0.78];
palette.eventLocal = palette.event;
palette.slice = [0.45 0.45 0.45];
end

function txt = buildTitle(corrOut)
parts = {sprintf('%s correlation', capitalizeMethod(corrOut.method)), ...
    sprintf('unit=%s', corrOut.unitLevel), ...
    sprintf('view=%s', corrOut.dataView)};
if isfield(corrOut, 'epochName') && ~isempty(corrOut.epochName)
    parts{end+1} = sprintf('epoch=%s', corrOut.epochName);
end
txt = strjoin(parts, ' | ');
end

function out = capitalizeMethod(in)
out = char(string(in));
if isempty(out)
    return
end
out(1) = upper(out(1));
end

function updateGroupBoxLayout(hFig, ax, H)
if ~isstruct(H) || ~isfield(H, 'boxes') || ~ishghandle(ax)
    return
end
drawnow limitrate nocallbacks;
axPos = getAxesBoxPosition(ax);
n = 0;
if isfield(H, 'runs') && ~isempty(H.runs)
    n = H.runs(end,2);
end
if n < 1
    return
end
leftGap = 0.012;
boxW = 0.024;
textW = 0.095;
for i = 1:size(H.runs,1)
    r0 = H.runs(i,1);
    r1 = H.runs(i,2);
    yTop = axPos(2) + axPos(4) * (1 - (r0-1)/n);
    yBot = axPos(2) + axPos(4) * (1 - r1/n);
    posBox = [max(0.005, axPos(1) - leftGap - boxW), yBot, boxW, max(0.002, yTop - yBot)];
    set(H.boxes(i), 'Position', posBox);
    posText = [max(0.001, posBox(1) - textW - 0.006), yBot + 0.5*(yTop-yBot) - 0.018, textW, 0.036];
    set(H.text(i), 'Position', posText);
end
end

function pos = getAxesBoxPosition(ax)
pos = get(ax, 'Position');
try
    ip = get(ax, 'InnerPosition');
    if isnumeric(ip) && numel(ip) == 4 && all(isfinite(ip)) && all(ip(3:4) > 0)
        pos = ip;
    end
catch
end
end

function [runs, labels] = contiguousGroupRuns(groups)
groups = string(groups(:));
if isempty(groups)
    runs = zeros(0,2);
    labels = strings(0,1);
    return
end
startIdx = 1;
runs = zeros(0,2);
labels = strings(0,1);
for i = 2:numel(groups)
    if groups(i) ~= groups(i-1)
        runs(end+1,:) = [startIdx, i-1]; %#ok<AGROW>
        labels(end+1,1) = groups(i-1); %#ok<AGROW>
        startIdx = i;
    end
end
runs(end+1,:) = [startIdx, numel(groups)];
labels(end+1,1) = groups(end);
end
