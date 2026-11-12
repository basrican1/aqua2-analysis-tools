function h = plotMetricBySliceCanonical(metricOut, varargin)
% Plot output of computeMetricBySliceCanonical().
%
% Options:
%   'comparisonLevel' one of: overall | by_animal | by_slice | by_cell
%   'groupLabel' color-group column (default auto)
%   'traceLabel' trace column, 'auto' (default), 'none', or explicit column
%   'stat' 'mean' (default) or 'median'

p = inputParser;
p.StructExpand = false;
addRequired(p, 'metricOut', @isstruct);
addParameter(p, 'comparisonLevel', 'by_slice', @(x) ischar(x) || isstring(x));
addParameter(p, 'cellInclusion', 'complete_cells_only', @(x) ischar(x) || isstring(x));
addParameter(p, 'groupLabel', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'traceLabel', 'auto', @(x) ischar(x) || isstring(x));
addParameter(p, 'stat', 'mean', @(x) ischar(x) || isstring(x));
addParameter(p, 'showError', true, @islogical);
addParameter(p, 'showTraces', true, @islogical);
addParameter(p, 'showGroupMean', true, @islogical);
addParameter(p, 'showGroupMeanError', true, @islogical);
addParameter(p, 'meanLineWidth', 2.8, @isnumeric);
addParameter(p, 'meanMode', 'overall', @(x) ischar(x) || isstring(x)); % overall|group
addParameter(p, 'meanColor', [0 0 0], @isnumeric);
addParameter(p, 'yLabelSuffix', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'showSignificance', true, @islogical);
addParameter(p, 'sigAlpha', 0.05, @isnumeric);
addParameter(p, 'sigMethod', 'auto', @(x) ischar(x) || isstring(x)); % auto|paired|unpaired
parse(p, metricOut, varargin{:});

comparisonLevel = lower(char(string(p.Results.comparisonLevel)));

metricCol = metricOut.metricCol;
stat = lower(char(string(p.Results.stat)));
if strcmp(stat, 'median')
    ycol = ['median_' metricCol];
else
    ycol = ['mean_' metricCol];
end

[S, gcol, tvals, tname, ycol] = buildPlotTable(metricOut, comparisonLevel, ...
    p.Results.groupLabel, p.Results.traceLabel, p.Results.cellInclusion, ycol);
if isempty(S)
    warning('No summary rows to plot.');
    h = [];
    return
end
if ~ismember(ycol, S.Properties.VariableNames)
    error('Expected summary column not found: %s', ycol);
end

hName = ['Canonical Metric: ' metricCol];
if isfield(metricOut, 'titleName') && ~isempty(metricOut.titleName)
    hName = char(string(metricOut.titleName));
end
h = figure('Name', hName, 'Color', 'w', 'Visible', 'on');
hold on;
groups = unique(string(S.(gcol)), 'stable');
if isempty(groups)
    groups = "all";
    S.(gcol) = repmat("all", height(S), 1);
end
groupColors = lines(max(1, numel(groups)));
TM = buildTraceEpochTable(S, tvals, gcol, ycol);
showTraces = p.Results.showTraces && ~strcmp(comparisonLevel, 'overall');
if showTraces
    traces = unique(tvals, 'stable');
    for i = 1:numel(traces)
        t = traces(i);
        R = S(tvals == t, :);
        if isempty(R)
            continue
        end
        [x, y] = aggregateByEpoch(R.epoch_idx, R.(ycol));
        g = string(R.(gcol)(1));
        gi = find(groups == g, 1, 'first');
        if isempty(gi), gi = 1; end
        c = groupColors(gi, :);
        plot(x, y, '-o', 'LineWidth', 1.2, 'Color', c, 'HandleVisibility', 'off');

        if p.Results.showError && ismember('sem', R.Properties.VariableNames) && strcmp(tname, 'none')
            [~, e] = aggregateByEpoch(R.epoch_idx, R.sem);
            errorbar(x, y, e, 'LineStyle', 'none', 'CapSize', 6, 'Color', c, 'HandleVisibility', 'off');
        end
    end
end

meanMode = lower(char(string(p.Results.meanMode)));
if p.Results.showGroupMean && strcmp(meanMode, 'group')
    for i = 1:numel(groups)
        g = groups(i);
        Rg = TM(string(TM.group) == g, :);
        [uEpoch, ~, gidx] = unique(Rg.epoch_idx, 'stable');
        m = accumarray(gidx, Rg.value, [], @mean, NaN);
        s = accumarray(gidx, Rg.value, [], @std, NaN);
        n = accumarray(gidx, ones(size(gidx)), [], @sum, NaN);
        sem = s ./ sqrt(max(n,1));
        c = groupColors(i,:);
        plot(uEpoch, m, '-o', 'Color', c, 'LineWidth', p.Results.meanLineWidth, ...
            'DisplayName', char(g));
        if p.Results.showGroupMeanError
            errorbar(uEpoch, m, sem, 'LineStyle', 'none', 'CapSize', 8, ...
                'Color', c, 'LineWidth', max(1, p.Results.meanLineWidth-1), 'HandleVisibility', 'off');
        end
    end
elseif p.Results.showGroupMean
    [uEpoch, ~, gidx] = unique(TM.epoch_idx, 'stable');
    m = accumarray(gidx, TM.value, [], @mean, NaN);
    s = accumarray(gidx, TM.value, [], @std, NaN);
    n = accumarray(gidx, ones(size(gidx)), [], @sum, NaN);
    sem = s ./ sqrt(max(n,1));
    c = p.Results.meanColor;
    plot(uEpoch, m, '-o', 'Color', c, 'LineWidth', p.Results.meanLineWidth, ...
        'DisplayName', 'Overall mean');
    if p.Results.showGroupMeanError
        errorbar(uEpoch, m, sem, 'LineStyle', 'none', 'CapSize', 8, ...
            'Color', c, 'LineWidth', max(1, p.Results.meanLineWidth-1), 'HandleVisibility', 'off');
    end
end

if p.Results.showSignificance
    try
        baselineIdx = min(double(metricOut.epochTable.epoch_idx));
        sig = bwaComputeEpochVsBaselineSig(TM, baselineIdx, char(string(p.Results.sigMethod)), p.Results.sigAlpha);
        addEpochSigStars(sig, p.Results.sigAlpha);
    catch
    end
end

% Always include color key for groups.
for i = 1:numel(groups)
    plot(NaN, NaN, '-', 'Color', groupColors(i,:), 'LineWidth', 2.0, ...
        'DisplayName', ['Animal: ' char(groups(i))]);
end

epochNames = metricOut.epochTable.epoch_name;
xticks(metricOut.epochTable.epoch_idx);
xticklabels(cellstr(epochNames));
xtickangle(30);
xlabel('Epoch');
ySuffix = char(string(p.Results.yLabelSuffix));
ylabel([prettyMetricLabel(metricCol) ySuffix], 'Interpreter', 'none');
title(sprintf('%s | traces: %s', hName, strrep(tname, '_', '\_')), 'Interpreter', 'none');
lg = legend('Location', 'best');
set(lg, 'Interpreter', 'none');
grid on;
box on;
hold off;
end

function [S, gcol, tvals, tname, ycol] = buildPlotTable(metricOut, comparisonLevel, groupLabel, traceLabel, cellInclusion, ycol)
S = metricOut.summaryByGroupEpoch;
if isempty(S)
    gcol = 'file_id';
    tvals = strings(0,1);
    tname = 'none';
    return
end

if isempty(groupLabel)
    gcol = 'animal_id';
else
    gcol = char(string(groupLabel));
end
if ~ismember(gcol, S.Properties.VariableNames)
    if ismember('animal_id', S.Properties.VariableNames)
        gcol = 'animal_id';
    else
        gcol = 'file_id';
    end
end

switch comparisonLevel
    case 'overall'
        tname = 'none';
        tvals = repmat("all", height(S), 1);
    case 'by_animal'
        tname = 'animal_id';
        if ismember('animal_id', S.Properties.VariableNames)
            tvals = string(S.animal_id);
            gcol = 'animal_id';
        else
            tvals = string(S.(gcol));
            tname = gcol;
        end
    case 'by_slice'
        if ismember('animal_id', S.Properties.VariableNames) && ismember('slice_num', S.Properties.VariableNames)
            tname = 'animal_slice';
            tvals = string(S.animal_id) + "|" + string(S.slice_num);
            if isempty(groupLabel)
                gcol = 'animal_id';
            end
        elseif ismember('file_id', S.Properties.VariableNames)
            tname = 'file_id';
            tvals = string(S.file_id);
        else
            tname = 'none';
            tvals = repmat("all", height(S), 1);
        end
    case 'by_cell'
        [S, gcol, ycol] = metricByCellSummary(metricOut, gcol, cellInclusion);
        if isempty(S)
            tname = 'none';
            tvals = strings(0,1);
            return
        end
        tname = 'primary_cell_uid';
        tvals = string(S.primary_cell_uid);
    otherwise
        req = char(string(traceLabel));
        if strcmpi(req, 'none') || isempty(req)
            tname = 'none';
            tvals = repmat("all", height(S), 1);
        elseif strcmpi(req, 'auto')
            if ismember('file_id', S.Properties.VariableNames)
                tname = 'file_id';
                tvals = string(S.file_id);
            else
                tname = gcol;
                tvals = string(S.(gcol));
            end
        elseif ismember(req, S.Properties.VariableNames)
            tname = req;
            tvals = string(S.(req));
        else
            tname = 'none';
            tvals = repmat("all", height(S), 1);
        end
end
end

function [S, gcol, ycol] = metricByCellSummary(metricOut, gcol, cellInclusion)
S = table();
ycol = ['mean_' metricOut.metricCol];
if ~isfield(metricOut, 'eventValues') || isempty(metricOut.eventValues)
    warning('plotMetricBySliceCanonical:NoEventValues', 'Cell-level plotting unavailable (eventValues missing).');
    return
end
T = metricOut.eventValues;
if ~ismember('primary_cell_uid', T.Properties.VariableNames)
    warning('plotMetricBySliceCanonical:NoPrimaryCell', 'Cell-level plotting unavailable (primary_cell_uid missing).');
    return
end
mcol = metricOut.metricCol;
keep = isfinite(double(T.(mcol))) & T.epoch_idx > 0 & string(T.primary_cell_uid) ~= "";
T = T(keep, :);
if isempty(T)
    return
end
if ~ismember(gcol, T.Properties.VariableNames)
    if ismember('animal_id', T.Properties.VariableNames)
        gcol = 'animal_id';
    else
        gcol = 'file_id';
    end
end
S = groupsummary(T, {gcol,'primary_cell_uid','epoch_idx','epoch_name'}, 'mean', mcol);
S.(ycol) = S.(['mean_' mcol]);
S = expandCellEpochGrid(S, gcol, ycol, metricOut.epochTable);
S = applyCellInclusionPolicy(S, gcol, 'primary_cell_uid', ycol, metricOut.epochTable, cellInclusion);
S.sem = zeros(height(S),1);
end

function S = expandCellEpochGrid(Sin, gcol, ycol, epochTable)
if isempty(Sin) || isempty(epochTable)
    S = Sin;
    return
end
eIdx = epochTable.epoch_idx(:);
if ismember('epoch_name', epochTable.Properties.VariableNames)
    eName = string(epochTable.epoch_name(:));
else
    eName = "E" + string(eIdx);
end
grp = unique(Sin(:, {gcol, 'primary_cell_uid'}), 'rows', 'stable');
S = table();
for i = 1:height(grp)
    G = repmat(grp(i,:), numel(eIdx), 1);
    G.epoch_idx = eIdx;
    G.epoch_name = eName;
    J = outerjoin(G, Sin(:, {gcol,'primary_cell_uid','epoch_idx',ycol}), ...
        'Keys', {gcol,'primary_cell_uid','epoch_idx'}, 'Type', 'left', 'MergeKeys', true);
    if ~ismember('epoch_name', J.Properties.VariableNames)
        J.epoch_name = eName;
    end
    S = [S; J(:, {gcol,'primary_cell_uid','epoch_idx','epoch_name',ycol})]; %#ok<AGROW>
end
end

function S = applyCellInclusionPolicy(S, gcol, cellCol, ycol, epochTable, cellInclusion)
mode = lower(char(string(cellInclusion)));
if isempty(S)
    return
end
switch mode
    case 'available_case'
        return
    otherwise % complete_cells_only
        nEpoch = numel(unique(epochTable.epoch_idx));
        keys = unique(S(:, {gcol, cellCol}), 'rows', 'stable');
        keep = false(height(S),1);
        for i = 1:height(keys)
            idx = string(S.(gcol)) == string(keys.(gcol)(i)) & string(S.(cellCol)) == string(keys.(cellCol)(i));
            nPresent = sum(isfinite(double(S.(ycol)(idx))));
            if nPresent == nEpoch
                keep(idx) = true;
            end
        end
        S = S(keep, :);
end
end

function [x, y] = aggregateByEpoch(epochIdx, vals)
epochIdx = double(epochIdx(:));
vals = double(vals(:));
keep = isfinite(epochIdx) & isfinite(vals);
epochIdx = epochIdx(keep);
vals = vals(keep);
if isempty(epochIdx)
    x = zeros(0,1);
    y = zeros(0,1);
    return
end
[x, ~, gidx] = unique(epochIdx, 'stable');
y = accumarray(double(gidx), vals, [], @mean, NaN);
[x, ord] = sort(x);
y = y(ord);
end

function TM = buildTraceEpochTable(S, tvals, gcol, ycol)
tr = unique(tvals, 'stable');
TM = table(strings(0,1), strings(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', {'trace','group','epoch_idx','value'});
for i = 1:numel(tr)
    R = S(tvals == tr(i), :);
    if isempty(R)
        continue
    end
    [x, y] = aggregateByEpoch(R.epoch_idx, R.(ycol));
    t = repmat(tr(i), numel(x), 1);
    g = repmat(string(R.(gcol)(1)), numel(x), 1);
    TM = [TM; table(t, g, x, y, 'VariableNames', TM.Properties.VariableNames)]; %#ok<AGROW>
end
end

function addEpochSigStars(sig, alpha)
if isempty(sig)
    return
end
yl = ylim;
dy = 0.04 * max(eps, yl(2) - yl(1));
yStar = yl(2) - 0.6*dy;
yTxt = yl(2) - 1.8*dy;
for i = 1:height(sig)
    if sig.is_sig(i)
        text(sig.epoch_idx(i), yStar, '*', 'Color', [0 0 0], 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Interpreter', 'none');
        if isfinite(sig.p_value(i))
            text(sig.epoch_idx(i), yTxt, sprintf('p=%.3g', sig.p_value(i)), ...
                'Color', [0.2 0.2 0.2], 'FontSize', 8, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Interpreter', 'none');
        end
    end
end
ylim([yl(1), yl(2) + 2.5*dy]);
txt = sprintf('Significance vs baseline (alpha=%.3g)', alpha);
titleObj = get(gca, 'Title');
if isgraphics(titleObj)
    set(titleObj, 'String', sprintf('%s | %s', char(titleObj.String), txt), 'Interpreter', 'none');
end
end

function lbl = prettyMetricLabel(metricCol)
mc = char(string(metricCol));
switch lower(mc)
    case 'area_px2', lbl = 'area (px^2)';
    case 'area_um2', lbl = 'area (um^2)';
    case 'duration_frames', lbl = 'duration (frames)';
    case 'duration_sec', lbl = 'duration (s)';
    case 'rise19_frames', lbl = 'rise19 (frames)';
    case 'rise19_sec', lbl = 'rise19 (s)';
    case 'fall91_frames', lbl = 'fall91 (frames)';
    case 'fall91_sec', lbl = 'fall91 (s)';
    case 'decay_tau_frames', lbl = 'decay tau (frames)';
    case 'decay_tau_sec', lbl = 'decay tau (s)';
    case 'perimeter_px', lbl = 'perimeter (px)';
    case 'perimeter_um', lbl = 'perimeter (um)';
    case 'propagation_distance_px', lbl = 'propagation distance (px)';
    case 'propagation_distance_um', lbl = 'propagation distance (um)';
    case 'mean_soma_centroid_dist_px', lbl = 'mean soma centroid dist (px)';
    case 'mean_soma_centroid_dist_um', lbl = 'mean soma centroid dist (um)';
    case 'mean_soma_start_dist_px', lbl = 'mean soma start dist (px)';
    case 'mean_soma_start_dist_um', lbl = 'mean soma start dist (um)';
    case 'mean_soma_end_dist_px', lbl = 'mean soma end dist (px)';
    case 'mean_soma_end_dist_um', lbl = 'mean soma end dist (um)';
    case 'mean_neighbor_entry_delay_frames', lbl = 'mean neighbor entry delay (frames)';
    case 'mean_neighbor_entry_delay_sec', lbl = 'mean neighbor entry delay (s)';
    case 'mean_perimeter_px', lbl = 'mean perimeter (px)';
    case 'mean_perimeter_um', lbl = 'mean perimeter (um)';
    case 'mean_propagation_distance_px', lbl = 'mean propagation distance (px)';
    case 'mean_propagation_distance_um', lbl = 'mean propagation distance (um)';
    otherwise, lbl = mc;
end
end
