function out = plotBWA_Aqua2PolarsBySlice(polData, varargin)
% Plot spider/polar panels from precomputed BWA_Aqua2PolarsBySlice data.
%
% Inputs:
%   polData: output struct from computeBWA_Aqua2PolarsBySliceData()
%
% Name-Value options:
%   'spiderInputScale'      : 'fold' (default) or 'log2'
%   'spiderSymmetricLimit'  : [] auto, scalar, or 1xN vector
%   'spiderScaleMode'       : 'per_metric' (default) or 'global'
%   'spiderMetricClipPct'   : percentile for auto limits (default 95)
%   'branchMode'            : 'C' (default), 'E', or 'both'
%   'comparisonLevel'       : '' (auto), 'by_slice', 'by_event',
%                             'by_cell', or 'by_branch'
%   'aggregateEBy'          : 'slice' (default) or 'events'
%   'aggregateCBy'          : 'slice' (default) or 'cells'
%   'showStatsFigure'       : true/false (default false)
%   'statsFigureName'       : figure name for stats table
%   'exportCompanionTable'  : true/false (default true)
%   'companionTableVarName' : base workspace variable name
%   'figureName'            : figure name (default 'BAPol')
%   'showLegend'           : true/false (default true)

    p = inputParser;
    p.StructExpand = false;
    addRequired(p, 'polData', @isstruct);
    addParameter(p,'spiderInputScale', 'fold', @(x) ischar(x) || isstring(x));
    addParameter(p,'spiderSymmetricLimit', [], @isnumeric);
    addParameter(p,'spiderScaleMode', 'per_metric', @(x) any(strcmpi(string(x), ["per_metric","global"])));
    addParameter(p,'spiderMetricClipPct', 95, @isnumeric);
    addParameter(p,'branchMode', 'C', @(x) any(strcmpi(string(x), ["E","C","both"])));
    addParameter(p,'comparisonLevel', '', @(x) any(strcmpi(string(x), ["","by_slice","by_event","by_cell","by_branch"])));
    addParameter(p,'aggregateEBy', 'slice', @(x) any(strcmpi(string(x), ["slice","events"])));
    addParameter(p,'aggregateCBy', 'slice', @(x) any(strcmpi(string(x), ["slice","cells"])));
    addParameter(p,'showStatsFigure', false, @islogical);
    addParameter(p,'statsFigureName', 'BAPol_stats', @(x) ischar(x) || isstring(x));
    addParameter(p,'exportCompanionTable', true, @islogical);
    addParameter(p,'companionTableVarName', 'BWA_spiderFoldTable', @(x) ischar(x) || isstring(x));
    addParameter(p,'figureName', 'BAPol', @(x) ischar(x) || isstring(x));
    addParameter(p,'showLegend', true, @islogical);
    addParameter(p,'style', struct(), @isstruct);
    parse(p, polData, varargin{:});
    style = mergeBWAPlotStyle(getBWAPlotStyle(), p.Results.style);
    branchMode = normalizeSpiderBranchMode(p.Results.branchMode);
    [aggregateEBy, aggregateCBy] = resolveSpiderSliceAggregation(branchMode, p.Results.comparisonLevel, p.Results.aggregateEBy, p.Results.aggregateCBy);

    needFields = {'labels','valuesA','valuesC','uniqueFiles','subplotGrid'};
    for ii = 1:numel(needFields)
        if ~isfield(polData, needFields{ii})
            error('plotBWA_Aqua2PolarsBySlice:MissingField', ...
                'polData is missing required field: %s', needFields{ii});
        end
    end

    nFiles = numel(polData.uniqueFiles);
    if nFiles == 0
        warning('plotBWA_Aqua2PolarsBySlice:NoFiles', 'No file IDs found; skipping plot.');
        out = struct('figure', [], 'spider', [], 'companionTable', table());
        return;
    end

    fName = char(string(p.Results.figureName));
    h_pol = findall(0, 'Type', 'figure', 'Name', fName);
    if isempty(h_pol)
        figPos = style.figurePositionPolars;
        if numel(polData.subplotGrid) >= 1
            nRows = max(1, polData.subplotGrid(1));
        else
            nRows = 1;
        end
        figPos(4) = max(figPos(4), 240 * nRows + 140);
        h_pol = figure('Position', figPos, 'Name', fName);
    else
        figure(h_pol);
    end

    [panelValuesA, panelValuesC, aggInfo] = resolvePanelBranchValues(polData, aggregateEBy, aggregateCBy);
    companionTbl = buildCompanionTable(panelValuesA, panelValuesC, polData.labels, ...
        polData.uniqueFiles, p.Results.spiderInputScale, aggInfo);

    limForPanels = resolvePanelLimitFromPolars(panelValuesA, panelValuesC, ...
        p.Results.spiderInputScale, p.Results.spiderSymmetricLimit, numel(polData.labels), ...
        p.Results.spiderScaleMode, p.Results.spiderMetricClipPct);

    if p.Results.exportCompanionTable
        assignin('base', char(string(p.Results.companionTableVarName)), companionTbl);
        fprintf('[plotBWA_Aqua2PolarsBySlice] Companion table exported to base workspace variable: %s\n', ...
            char(string(p.Results.companionTableVarName)));
    end

    tileSpacing = style.tileSpacing;
    tilePadding = style.tilePadding;
    if strcmpi(char(string(tileSpacing)), 'none')
        tileSpacing = 'compact';
    end
    if strcmpi(char(string(tilePadding)), 'loose')
        tilePadding = 'compact';
    end
    tl = tiledlayout(h_pol, polData.subplotGrid(1), polData.subplotGrid(2), ...
        'TileSpacing', tileSpacing, 'Padding', tilePadding);
    firstAx = [];
    for pp = 1:nFiles
        ax = nexttile(tl);
        if isempty(firstAx)
            firstAx = ax;
        end
        thisFileID = char(string(polData.uniqueFiles(pp)));
        thisFileLabel = thisFileID;
        if isfield(polData, 'fileLabels') && iscell(polData.fileLabels) && numel(polData.fileLabels) >= pp
            thisFileLabel = char(string(polData.fileLabels{pp}));
        end
        sigA = [];
        sigC = [];
        if isfield(polData, 'sigMaskA') && size(polData.sigMaskA, 1) >= pp
            sigA = logical(polData.sigMaskA(pp, :));
        end
        if isfield(polData, 'sigMaskC') && size(polData.sigMaskC, 1) >= pp
            sigC = logical(polData.sigMaskC(pp, :));
        end
        plotPolarsSpiderSubplot(panelValuesA(pp,:), panelValuesC(pp,:), polData.labels, ...
            thisFileLabel, ...
            'inputScale', p.Results.spiderInputScale, ...
            'symmetricLimit', limForPanels, ...
            'branchMode', branchMode, ...
            'sigMaskA', sigA, ...
            'sigMaskC', sigC, ...
            'style', style, ...
            'figure', h_pol, ...
            'axes', ax);
    end

    if strcmpi(char(string(p.Results.spiderInputScale)), 'log2')
        scaleTxt = upper(char(string(p.Results.spiderScaleMode)));
        if strcmp(branchMode, 'both')
            title(tl, sprintf('Spider Plots (log2FC, baseline=0; %s scaling; Event-weighted + Cell-weighted)', scaleTxt), ...
                'Interpreter', 'none');
        elseif strcmp(branchMode, 'E')
            title(tl, sprintf('Spider Plots (log2FC, baseline=0; %s scaling; Event-weighted only)', scaleTxt), ...
                'Interpreter', 'none');
        else
            title(tl, sprintf('Spider Plots (log2FC, baseline=0; %s scaling; Cell-weighted only)', scaleTxt), ...
                'Interpreter', 'none');
        end
    else
        scaleTxt = upper(char(string(p.Results.spiderScaleMode)));
        if strcmp(branchMode, 'both')
            title(tl, sprintf('Spider Plots (log2 view from fold, baseline=0; %s scaling; Event-weighted + Cell-weighted)', scaleTxt), ...
                'Interpreter', 'none');
        elseif strcmp(branchMode, 'E')
            title(tl, sprintf('Spider Plots (log2 view from fold, baseline=0; %s scaling; Event-weighted only)', scaleTxt), ...
                'Interpreter', 'none');
        else
            title(tl, sprintf('Spider Plots (log2 view from fold, baseline=0; %s scaling; Cell-weighted only)', scaleTxt), ...
                'Interpreter', 'none');
        end
    end

    if p.Results.showLegend && ~isempty(firstAx) && isgraphics(firstAx)
        hold(firstAx, 'on');
        d1 = gobjects(1);
        d2 = plot(firstAx, NaN, NaN, '-', 'Color', style.colors.groupC, 'LineWidth', style.lineWidthMain);
        d3 = plot(firstAx, NaN, NaN, '-', 'Color', style.colors.baseline, 'LineWidth', style.lineWidthBaseline);
        if strcmp(branchMode, 'both') || strcmp(branchMode, 'E')
            d1 = plot(firstAx, NaN, NaN, '-', 'Color', style.colors.groupA, 'LineWidth', style.lineWidthMain);
        end
        hold(firstAx, 'off');
        if strcmp(branchMode, 'both')
            legend(firstAx, [d1 d2 d3], {'Event-weighted (E)', 'Cell-weighted (C)', 'Baseline (0)'}, ...
                'FontSize', style.fontSizeLegend, 'Location', style.legendLocation, 'Interpreter', 'none');
        elseif strcmp(branchMode, 'E')
            legend(firstAx, [d1 d3], {'Event-weighted (E)', 'Baseline (0)'}, ...
                'FontSize', style.fontSizeLegend, 'Location', style.legendLocation, 'Interpreter', 'none');
        else
            legend(firstAx, [d2 d3], {'Cell-weighted (C)', 'Baseline (0)'}, ...
                'FontSize', style.fontSizeLegend, 'Location', style.legendLocation, 'Interpreter', 'none');
        end
    end

    out = struct();
    out.figure = h_pol;
    out.spider = [];
    out.companionTable = companionTbl;
    out.summaryTable = companionTbl;
    out.statsSummary = buildSpiderSliceStatsSummary( ...
        polData.labels, polData.uniqueFiles, polData, panelValuesA, panelValuesC, ...
        branchMode, aggregateEBy, aggregateCBy);
    out.statsInfo = struct( ...
        'baseline_reference', 'log2FC = 0', ...
        'comparison_level', string(resolveSpiderSliceComparisonLevelLabel(p.Results.comparisonLevel, aggregateEBy, aggregateCBy)), ...
        'aggregate_E_by', string(aggregateEBy), ...
        'aggregate_C_by', string(aggregateCBy), ...
        'sig_source', string(branchMode), ...
        'requested_method', "precomputed_per_slice", ...
        'alpha', NaN, ...
        'n_files', nFiles);
    out.aggregationInfo = aggInfo;
    out.panelLimitsLog2 = limForPanels;
    out.statsFigure = [];
    if p.Results.showStatsFigure
        out.statsFigure = displaySpiderSliceStatsFigure(out.statsSummary, char(string(p.Results.statsFigureName)), out.statsInfo);
    end
end

function mode = normalizeSpiderBranchMode(modeIn)
mode = upper(char(string(modeIn)));
switch lower(strtrim(mode))
        case 'e'
        mode = 'E';
    case 'c'
        mode = 'C';
    case 'both'
        mode = 'both';
    otherwise
        error('plotBWA_Aqua2PolarsBySlice:BadBranchMode', ...
            'branchMode must be ''E'', ''C'', or ''both''.');
end
end

function [aggregateEBy, aggregateCBy] = resolveSpiderSliceAggregation(branchMode, comparisonLevel, aggregateEIn, aggregateCIn)
aggregateEBy = lower(char(string(aggregateEIn)));
aggregateCBy = lower(char(string(aggregateCIn)));
level = lower(strtrim(char(string(comparisonLevel))));
if isempty(level)
    return
end
switch level
    case 'by_slice'
        aggregateEBy = 'slice';
        aggregateCBy = 'slice';
    case 'by_event'
        aggregateEBy = 'events';
        aggregateCBy = 'slice';
    case 'by_cell'
        aggregateEBy = 'slice';
        aggregateCBy = 'cells';
    case 'by_branch'
        aggregateEBy = 'events';
        aggregateCBy = 'cells';
    otherwise
        error('plotBWA_Aqua2PolarsBySlice:BadComparisonLevel', ...
            'comparisonLevel must be '''', ''by_slice'', ''by_event'', ''by_cell'', or ''by_branch''.');
end
if strcmpi(branchMode, 'E') && strcmp(aggregateCBy, 'cells')
    aggregateCBy = 'slice';
elseif strcmpi(branchMode, 'C') && strcmp(aggregateEBy, 'events')
    aggregateEBy = 'slice';
end
end


function label = resolveSpiderSliceComparisonLevelLabel(comparisonLevel, aggregateEBy, aggregateCBy)
level = lower(strtrim(char(string(comparisonLevel))));
if ~isempty(level)
    label = string(level);
    return
end
if strcmpi(aggregateEBy, 'events') && strcmpi(aggregateCBy, 'cells')
    label = "by_branch";
elseif strcmpi(aggregateEBy, 'events')
    label = "by_event";
elseif strcmpi(aggregateCBy, 'cells')
    label = "by_cell";
else
    label = "by_slice";
end
end

function [valuesA, valuesC, aggInfo] = resolvePanelBranchValues(polData, aggregateEBy, aggregateCBy)
valuesA = double(polData.valuesA);
valuesC = double(polData.valuesC);
usedEntityA = false(1, size(valuesA, 2));
usedEntityC = false(1, size(valuesC, 2));
if strcmpi(aggregateEBy, 'events')
    [valuesA, usedEntityA] = applyEntityAggregationByFile(valuesA, polData, 'entityValuesAByFile');
end
if strcmpi(aggregateCBy, 'cells')
    [valuesC, usedEntityC] = applyEntityAggregationByFile(valuesC, polData, 'entityValuesCByFile');
end
aggInfo = struct();
aggInfo.requestedAggregateEBy = string(aggregateEBy);
aggInfo.requestedAggregateCBy = string(aggregateCBy);
aggInfo.usedEntityA = usedEntityA;
aggInfo.usedEntityC = usedEntityC;
aggInfo.usedAggregateEBy = repmat(string("slice"), 1, size(valuesA, 2));
aggInfo.usedAggregateCBy = repmat(string("slice"), 1, size(valuesC, 2));
aggInfo.usedAggregateEBy(usedEntityA) = "events";
aggInfo.usedAggregateCBy(usedEntityC) = "cells";
end

function [values, usedEntity] = applyEntityAggregationByFile(values, polData, fieldName)
usedEntity = false(1, size(values, 2));
if ~isstruct(polData) || ~isfield(polData, fieldName) || ~iscell(polData.(fieldName))
    return
end
entityByFile = polData.(fieldName);
if numel(entityByFile) ~= size(values, 2)
    return
end
for j = 1:size(values, 2)
    perMetric = entityByFile{j};
    if ~iscell(perMetric)
        continue
    end
    for i = 1:min(numel(perMetric), size(values, 1))
        v = double(perMetric{i});
        v = v(isfinite(v));
        if ~isempty(v)
            values(i, j) = mean(v, 'omitnan');
            usedEntity(j) = true;
        end
    end
end
end

function T = buildCompanionTable(valuesA, valuesC, labels, fileIDs, inputScale, aggInfo)
valuesA = double(valuesA);
valuesC = double(valuesC);
nFiles = size(valuesA, 1);
nMetrics = size(valuesA, 2);
labels = cellstr(string(labels(:)))';
fileIDs = cellstr(string(fileIDs(:)));
if numel(labels) ~= nMetrics
    error('plotBWA_Aqua2PolarsBySlice:LabelMismatch', 'labels length must equal number of metric columns.');
end
if numel(fileIDs) ~= nFiles
    error('plotBWA_Aqua2PolarsBySlice:FileMismatch', 'fileIDs length must equal number of file rows.');
end

scale = lower(char(string(inputScale)));
switch scale
    case 'fold'
        foldA = valuesA;
        foldC = valuesC;
        log2A = NaN(size(foldA));
        log2C = NaN(size(foldC));
        mA = isfinite(foldA) & foldA > 0;
        mC = isfinite(foldC) & foldC > 0;
        log2A(mA) = log2(foldA(mA));
        log2C(mC) = log2(foldC(mC));
    case 'log2'
        log2A = valuesA;
        log2C = valuesC;
        foldA = 2.^log2A;
        foldC = 2.^log2C;
    otherwise
        error('plotBWA_Aqua2PolarsBySlice:BadScale', 'spiderInputScale must be ''fold'' or ''log2''.');
end

nRows = nFiles * nMetrics * 2;
fileCol = strings(nRows,1);
metricCol = strings(nRows,1);
groupCol = strings(nRows,1);
foldCol = NaN(nRows,1);
log2Col = NaN(nRows,1);
r = 0;
for i = 1:nFiles
    for j = 1:nMetrics
        r = r + 1;
        fileCol(r) = string(fileIDs{i});
        metricCol(r) = string(labels{j});
        groupCol(r) = "E";
        foldCol(r) = foldA(i,j);
        log2Col(r) = log2A(i,j);

        r = r + 1;
        fileCol(r) = string(fileIDs{i});
        metricCol(r) = string(labels{j});
        groupCol(r) = "C";
        foldCol(r) = foldC(i,j);
        log2Col(r) = log2C(i,j);
    end
end
T = table(fileCol, metricCol, groupCol, foldCol, log2Col, ...
    'VariableNames', {'file_id', 'metric', 'group', 'value_fold', 'value_log2'});
if nargin >= 6 && isstruct(aggInfo)
    aggCol = strings(nRows, 1);
    entityCol = false(nRows, 1);
    r = 0;
    for i = 1:nFiles
        for j = 1:nMetrics
            r = r + 1;
            aggCol(r) = aggInfo.usedAggregateEBy(j);
            entityCol(r) = aggInfo.usedEntityA(j);
            r = r + 1;
            aggCol(r) = aggInfo.usedAggregateCBy(j);
            entityCol(r) = aggInfo.usedEntityC(j);
        end
    end
    T.aggregate_by = aggCol;
    T.entity_backed = entityCol;
end
end

function T = buildSpiderSliceStatsSummary(labels, fileIDs, polData, valuesA, valuesC, branchMode, aggregateEBy, aggregateCBy)
labels = string(labels(:));
fileIDs = string(fileIDs(:));
if isfield(polData, 'fileLabels') && iscell(polData.fileLabels) && numel(polData.fileLabels) == numel(fileIDs)
    sliceLabels = string(polData.fileLabels(:));
else
    sliceLabels = fileIDs;
end
metricFamily = string(bwaResolveMetricFamily(labels));
metricFamily = metricFamily(:);
if numel(metricFamily) ~= numel(labels)
    metricFamily = repmat("", numel(labels), 1);
end
branchList = resolveStatsBranches(branchMode);

T = table();
for b = 1:numel(branchList)
    branch = branchList(b);
    if branch == "E"
        vals = double(valuesA);
        sig = false(size(vals));
        nVals = nan(size(vals));
        if isfield(polData, 'sigMaskA') && isnumeric(polData.sigMaskA) && all(size(polData.sigMaskA) == size(vals))
            sig = logical(polData.sigMaskA);
        end
        aggBy = repmat(string(aggregateEBy), numel(labels), 1);
        entityBacked = false(numel(labels), 1);
        if isfield(polData, 'entityValuesAByFile') && iscell(polData.entityValuesAByFile)
            entityBacked = metricEntityBacked(polData.entityValuesAByFile, size(vals, 1));
            aggBy(entityBacked) = "events";
            nVals = metricEntityCounts(polData.entityValuesAByFile, size(vals, 1));
        end
    else
        vals = double(valuesC);
        sig = false(size(vals));
        nVals = nan(size(vals));
        if isfield(polData, 'sigMaskC') && isnumeric(polData.sigMaskC) && all(size(polData.sigMaskC) == size(vals))
            sig = logical(polData.sigMaskC);
        end
        aggBy = repmat(string(aggregateCBy), numel(labels), 1);
        entityBacked = false(numel(labels), 1);
        if isfield(polData, 'entityValuesCByFile') && iscell(polData.entityValuesCByFile)
            entityBacked = metricEntityBacked(polData.entityValuesCByFile, size(vals, 1));
            aggBy(entityBacked) = "cells";
            nVals = metricEntityCounts(polData.entityValuesCByFile, size(vals, 1));
        end
    end

    nMetrics = min([ ...
        numel(labels), ...
        numel(metricFamily), ...
        numel(aggBy), ...
        numel(entityBacked), ...
        size(vals, 2), ...
        size(nVals, 2), ...
        size(sig, 2)]);
    if nMetrics <= 0
        continue
    end

    labelsUse = labels(1:nMetrics);
    metricFamilyUse = metricFamily(1:nMetrics);
    aggByUse = aggBy(1:nMetrics);
    entityBackedUse = entityBacked(1:nMetrics);
    vals = vals(:, 1:nMetrics);
    nVals = nVals(:, 1:nMetrics);
    sig = sig(:, 1:nMetrics);

    Tbranch = table(labelsUse, metricFamilyUse, repmat(branch, nMetrics, 1), aggByUse, entityBackedUse, ...
        'VariableNames', {'metric','metric_family','group','aggregate_by','entity_backed'});
    for i = 1:numel(sliceLabels)
        base = matlab.lang.makeValidName(char(sliceLabels(i)));
        Tbranch.(sprintf('%s_log2', base)) = vals(i, :)';
        Tbranch.(sprintf('%s_n', base)) = nVals(i, :)';
        Tbranch.(sprintf('%s_sig', base)) = sig(i, :)';
        Tbranch.(sprintf('%s_direction', base)) = significanceDirection(vals(i, :)', sig(i, :)');
    end
    T = [T; Tbranch]; %#ok<AGROW>
end
T = sortSpiderSliceStatsTable(T);
end

function T = sortSpiderSliceStatsTable(T)
if ~istable(T) || isempty(T) || ~ismember('metric', T.Properties.VariableNames)
    return
end
[order, familiesOrdered] = bwaOrderSpiderMetricsByFamily(string(T.metric));
T = T(order, :);
if ismember('metric_family', T.Properties.VariableNames)
    T.metric_family = familiesOrdered;
end
end

function branches = resolveStatsBranches(branchMode)
switch upper(char(string(branchMode)))
    case 'E'
        branches = "E";
    case 'C'
        branches = "C";
    otherwise
        branches = ["E"; "C"];
end
end

function tf = metricEntityBacked(entityByFile, nFiles)
nMetrics = numel(entityByFile);
tf = false(nMetrics, 1);
for j = 1:nMetrics
    perMetric = entityByFile{j};
    if ~iscell(perMetric)
        continue
    end
    for i = 1:min(numel(perMetric), nFiles)
        v = double(perMetric{i});
        v = v(isfinite(v));
        if ~isempty(v)
            tf(j) = true;
            break
        end
    end
end
end

function nVals = metricEntityCounts(entityByFile, nFiles)
nMetrics = numel(entityByFile);
nVals = nan(nFiles, nMetrics);
for j = 1:nMetrics
    perMetric = entityByFile{j};
    if ~iscell(perMetric)
        continue
    end
    for i = 1:min(numel(perMetric), nFiles)
        v = double(perMetric{i});
        v = v(isfinite(v));
        if ~isempty(v)
            nVals(i, j) = numel(v);
        end
    end
end
end

function d = significanceDirection(vals, sig)
vals = double(vals(:));
sig = logical(sig(:));
d = repmat("ns", numel(vals), 1);
d(sig & vals > 0) = "increase";
d(sig & vals < 0) = "decrease";
d(sig & vals == 0) = "significant_zero";
end

function hFig = displaySpiderSliceStatsFigure(T, figName, statsInfo)
if isempty(T) || ~istable(T)
    hFig = [];
    return
end
fName = char(string(figName));
old = findall(0, 'Type', 'figure', 'Name', fName);
if ~isempty(old)
    try
        delete(old);
    catch
    end
end
hFig = figure('Name', fName, 'NumberTitle', 'off', 'Color', 'w', 'Position', [180 180 1200 420]);
ax = axes('Parent', hFig, 'Position', [0 0 1 1], 'Visible', 'off'); %#ok<NASGU>

infoTxt = '';
if nargin >= 3 && isstruct(statsInfo)
    infoTxt = sprintf('Comparison: %s | E: %s | C: %s | baseline: %s | source: %s | method: %s', ...
        char(string(statsInfo.comparison_level)), ...
        char(string(statsInfo.aggregate_E_by)), ...
        char(string(statsInfo.aggregate_C_by)), ...
        char(string(statsInfo.baseline_reference)), ...
        char(string(statsInfo.sig_source)), ...
        char(string(statsInfo.requested_method)));
end
if ~isempty(infoTxt)
    uicontrol('Parent', hFig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.02 0.92 0.96 0.06], 'BackgroundColor', 'w', ...
        'HorizontalAlignment', 'left', 'String', infoTxt, 'FontWeight', 'bold');
end

data = table2cell(T);
for r = 1:size(data,1)
    for c = 1:size(data,2)
        if isstring(data{r,c})
            data{r,c} = char(data{r,c});
        elseif islogical(data{r,c})
            data{r,c} = double(data{r,c});
        end
    end
end

uitable('Parent', hFig, ...
    'Data', data, ...
    'ColumnName', T.Properties.VariableNames, ...
    'Units', 'normalized', ...
    'Position', [0.02 0.05 0.96 0.85]);
end

function lims = resolvePanelLimitFromPolars(valuesA, valuesC, inputScale, userLimit, nMetrics, scaleMode, clipPct)
if ~isempty(userLimit)
    if isscalar(userLimit)
        lims = repmat(abs(userLimit), 1, nMetrics);
    else
        if numel(userLimit) ~= nMetrics
            error('plotBWA_Aqua2PolarsBySlice:BadLimitSize', ...
                'spiderSymmetricLimit vector must match number of metrics.');
        end
        lims = abs(reshape(userLimit, 1, []));
    end
    lims(lims == 0) = 1;
    return
end

log2A = convertToLog2(double(valuesA), inputScale);
log2C = convertToLog2(double(valuesC), inputScale);

clipPct = double(clipPct);
if ~isfinite(clipPct)
    clipPct = 95;
end
clipPct = min(max(clipPct, 1), 100);

mode = lower(char(string(scaleMode)));
switch mode
    case 'global'
        vv = [log2A(:); log2C(:)];
        vv = abs(vv(isfinite(vv)));
        if isempty(vv)
            lims = ones(1, nMetrics);
        else
            mx = prctile(vv, clipPct);
            if mx <= 0 || ~isfinite(mx)
                mx = 1;
            end
            mx = ceil(mx * 2) / 2;
            lims = repmat(mx, 1, nMetrics);
        end
    otherwise % per_metric
        lims = NaN(1, nMetrics);
        for j = 1:nMetrics
            vv = [log2A(:,j); log2C(:,j)];
            vv = abs(vv(isfinite(vv)));
            if isempty(vv)
                lims(j) = 1;
            else
                mx = prctile(vv, clipPct);
                if mx <= 0 || ~isfinite(mx)
                    mx = 1;
                end
                lims(j) = ceil(mx * 2) / 2;
            end
        end
end
lims(lims == 0 | ~isfinite(lims)) = 1;
end

function log2X = convertToLog2(X, inputScale)
scale = lower(char(string(inputScale)));
switch scale
    case 'log2'
        log2X = X;
    case 'fold'
        log2X = NaN(size(X));
        m = isfinite(X) & X > 0;
        log2X(m) = log2(X(m));
    otherwise
        error('plotBWA_Aqua2PolarsBySlice:BadScale', 'spiderInputScale must be ''fold'' or ''log2''.');
end
end
