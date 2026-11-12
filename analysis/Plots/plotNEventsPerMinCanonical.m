function h = plotNEventsPerMinCanonical(freqOut, varargin)
% Plot output of computeNEventsPerMinCanonical().
%
% Options:
%   'comparisonLevel' one of: overall | by_animal | by_slice | by_cell
%   'groupLabel' column name for color grouping
%   'traceLabel' column name for individual traces, 'auto' (default), 'none', or explicit
%   'stat' one of: mean (default) or median
%   'showError' logical (default true)

p = inputParser;
p.StructExpand = false;
addRequired(p, 'freqOut', @isstruct);
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
parse(p, freqOut, varargin{:});

comparisonLevel = lower(char(string(p.Results.comparisonLevel)));

stat = lower(char(string(p.Results.stat)));
if strcmp(stat, 'median')
    ycol = 'median_rate_events_per_min';
else
    ycol = 'mean_rate_events_per_min';
end
[S, gcol, tvals, tname] = buildFreqPlotTable(freqOut, comparisonLevel, p.Results.groupLabel, p.Results.traceLabel, p.Results.cellInclusion);
if isempty(S)
    warning('No summary rows to plot.');
    h = [];
    return
end

h = figure('Name', 'Canonical Frequency Summary', 'Color', 'w', 'Visible', 'on');
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

        if p.Results.showError && ismember('sem_rate_events_per_min', R.Properties.VariableNames) && strcmp(tname, 'none')
            [~, e] = aggregateByEpoch(R.epoch_idx, R.sem_rate_events_per_min);
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
        baselineIdx = min(double(freqOut.epochTable.epoch_idx));
        sig = bwaComputeEpochVsBaselineSig(TM, baselineIdx, char(string(p.Results.sigMethod)), p.Results.sigAlpha);
        addEpochSigStars(sig, p.Results.sigAlpha);
    catch
    end
end
for i = 1:numel(groups)
    plot(NaN, NaN, '-', 'Color', groupColors(i,:), 'LineWidth', 2.0, ...
        'DisplayName', ['Animal: ' char(groups(i))]);
end

epochNames = freqOut.epochTable.epoch_name;
xticks(freqOut.epochTable.epoch_idx);
xticklabels(cellstr(epochNames));
xtickangle(30);
xlabel('Epoch');
ySuffix = char(string(p.Results.yLabelSuffix));
ylabel(['Events / min' ySuffix], 'Interpreter', 'none');
title(sprintf('Event Frequency by Epoch | traces: %s', strrep(tname, '_', '\_')), 'Interpreter', 'none');
lg = legend('Location', 'best');
set(lg, 'Interpreter', 'none');
grid on;
box on;
hold off;
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

function [S, gcol, tvals, tname] = buildFreqPlotTable(freqOut, level, groupLabel, traceLabel, cellInclusion)
if strcmp(level, 'by_cell') && isfield(freqOut, 'ratesByCellEpoch') && ~isempty(freqOut.ratesByCellEpoch)
    S = buildFreqByCellRows(freqOut.ratesByCellEpoch, freqOut.epochTable, cellInclusion);
    S.mean_rate_events_per_min = S.rate_events_per_min;
    S.median_rate_events_per_min = S.rate_events_per_min;
    S.sem_rate_events_per_min = zeros(height(S),1);
else
    S = freqOut.summaryByGroupEpoch;
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

switch lower(level)
    case 'overall'
        tname = 'none';
        tvals = repmat("all", height(S), 1);
    case 'by_animal'
        if ismember('animal_id', S.Properties.VariableNames)
            gcol = 'animal_id';
            tname = 'animal_id';
            tvals = string(S.animal_id);
        else
            tname = gcol;
            tvals = string(S.(gcol));
        end
    case 'by_slice'
        if ismember('animal_id', S.Properties.VariableNames) && ismember('slice_num', S.Properties.VariableNames)
            gcol = 'animal_id';
            tname = 'animal_slice';
            tvals = string(S.animal_id) + "|" + string(S.slice_num);
        elseif ismember('file_id', S.Properties.VariableNames)
            tname = 'file_id';
            tvals = string(S.file_id);
        else
            tname = 'none';
            tvals = repmat("all", height(S), 1);
        end
    case 'by_cell'
        if ismember('cell_uid', S.Properties.VariableNames)
            tname = 'cell_uid';
            tvals = string(S.cell_uid);
        else
            tname = 'none';
            tvals = repmat("all", height(S), 1);
        end
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

function S = buildFreqByCellRows(Tin, epochTable, cellInclusion)
if isempty(Tin)
    S = Tin;
    return
end
grpCols = {};
for c = {'animal_id','slice_num','file_id'}
    if ismember(c{1}, Tin.Properties.VariableNames)
        grpCols{end+1} = c{1}; %#ok<AGROW>
    end
end
if ~ismember('cell_uid', Tin.Properties.VariableNames)
    S = Tin;
    return
end
if isempty(epochTable) || ~ismember('epoch_idx', epochTable.Properties.VariableNames)
    S = Tin;
    return
end
eIdx = epochTable.epoch_idx(:);
if ismember('epoch_name', epochTable.Properties.VariableNames)
    eName = string(epochTable.epoch_name(:));
else
    eName = "E" + string(eIdx);
end
keys = unique(Tin(:, [grpCols {'cell_uid'}]), 'rows', 'stable');
S = table();
for i = 1:height(keys)
    G = repmat(keys(i,:), numel(eIdx), 1);
    G.epoch_idx = eIdx;
    G.epoch_name = eName;
    J = outerjoin(G, Tin, 'Keys', [grpCols {'cell_uid','epoch_idx'}], 'Type', 'left', 'MergeKeys', true);
    if ~ismember('rate_events_per_min', J.Properties.VariableNames)
        J.rate_events_per_min = NaN(height(J),1);
    end
    % Frequency missing cell-epoch combos are true zeros.
    J.rate_events_per_min(~isfinite(J.rate_events_per_min)) = 0;
    if ~ismember('epoch_name', J.Properties.VariableNames)
        J.epoch_name = eName;
    end
    colsKeep = [grpCols {'cell_uid','epoch_idx','epoch_name','rate_events_per_min'}];
    colsKeep = colsKeep(ismember(colsKeep, J.Properties.VariableNames));
    S = [S; J(:, colsKeep)]; %#ok<AGROW>
end

mode = lower(char(string(cellInclusion)));
if strcmp(mode, 'complete_cells_only')
    nEpoch = numel(unique(eIdx));
    keep = false(height(S),1);
    keys2 = unique(S(:, [grpCols {'cell_uid'}]), 'rows', 'stable');
    for i = 1:height(keys2)
        idx = true(height(S),1);
        for c = [grpCols {'cell_uid'}]
            idx = idx & string(S.(c{1})) == string(keys2.(c{1})(i));
        end
        if sum(idx) == nEpoch
            keep(idx) = true;
        end
    end
    S = S(keep,:);
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
titleObj = get(gca, 'Title');
if isgraphics(titleObj)
    set(titleObj, 'String', sprintf('%s | Significance vs baseline (alpha=%.3g)', ...
        char(titleObj.String), alpha), 'Interpreter', 'none');
end
end
