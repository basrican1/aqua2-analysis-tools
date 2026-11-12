function out = plotBAviolinSubplotSliceSpider(comp, timeFrameRangeEpochs, datTxt, varargin)
% Plot BA-violin companion spider with slice axes and epoch/post log2FC amplitudes.
% Supports standalone mode or embedding via provided figure/axes handles.
%
% By default the legacy wrapper passes the same combined post-drug epoch
% selection stored in result.spiderCombineEpochNames/result.spiderCombinedEpochLabel.
% If combineEpochNames is non-empty, the inset spider collapses the selected
% post epochs into one branch and compares that branch against the zero line.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'comp', @isstruct);
addRequired(p, 'timeFrameRangeEpochs', @isstruct);
addRequired(p, 'datTxt', @(x) ischar(x) || isstring(x));
addParameter(p, 'epochPair', [2 3], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'combineEpochNames', {}, @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'combinedEpochLabel', 'CNO_post', @(x) ischar(x) || isstring(x));
addParameter(p, 'branchMode', 'both', @(x) any(strcmpi(string(x), ["cell_only","c","both"])));
addParameter(p, 'figure', [], @(x) isempty(x) || ishghandle(x, 'figure'));
addParameter(p, 'axes', [], @(x) isempty(x) || ishghandle(x, 'axes'));
addParameter(p, 'style', struct(), @isstruct);
addParameter(p, 'sigAlpha', 0.05, @isnumeric);
parse(p, comp, timeFrameRangeEpochs, datTxt, varargin{:});

if ~isfield(comp, 'normByEvent') || ~isfield(comp, 'epMembIdx') || ~isfield(comp, 'fileID')
    warning('plotBAviolinSubplotSliceSpider:MissingFields', ...
        'Required comp fields missing (normByEvent/epMembIdx/fileID).');
    out = struct('figure', [], 'axes', [], 'spider', [], 'labels', {{}}, ...
        'valuesA', [], 'valuesC', []);
    return
end

normVals = double(comp.normByEvent(:));
epIdx = double(comp.epMembIdx(:));
fileID = string(comp.fileID(:));
n = min([numel(normVals), numel(epIdx), numel(fileID)]);
if isempty(n) || n == 0
    warning('plotBAviolinSubplotSliceSpider:NoData', 'No event-level data available for slice spider.');
    out = struct('figure', [], 'axes', [], 'spider', [], 'labels', {{}}, ...
        'valuesA', [], 'valuesC', []);
    return
end
normVals = normVals(1:n);
epIdx = epIdx(1:n);
fileID = fileID(1:n);

mask = isfinite(normVals) & ~isnan(epIdx) & (strlength(fileID) > 0);
normVals = normVals(mask);
epIdx = epIdx(mask);
fileID = fileID(mask);

if isempty(normVals)
    warning('plotBAviolinSubplotSliceSpider:NoFiniteData', ...
        'No finite normalized values found for slice spider.');
    out = struct('figure', [], 'axes', [], 'spider', [], 'labels', {{}}, ...
        'valuesA', [], 'valuesC', []);
    return
end

uFiles = unique(fileID, 'stable');
nFiles = numel(uFiles);
filesTbl = table();
if isfield(comp, 'filesTbl') && istable(comp.filesTbl)
    filesTbl = comp.filesTbl;
end
labels = makeSliceLabels(uFiles, filesTbl);

combineEpochNames = string(cellstr(string(p.Results.combineEpochNames(:))));
combineEpochNames = combineEpochNames(strlength(combineEpochNames) > 0);
combinedEpochLabel = char(string(p.Results.combinedEpochLabel));

epochPair = double(p.Results.epochPair(:))';
nEpoch = size(timeFrameRangeEpochs.frameRange, 1);
epochPair(epochPair < 1) = 1;
epochPair(epochPair > nEpoch) = nEpoch;
if nEpoch == 1
    epochPair = [1 1];
elseif numel(unique(epochPair)) == 1
    epochPair(2) = min(nEpoch, epochPair(1) + 1);
end
combineIdx = resolveCombinedEpochIdx(timeFrameRangeEpochs, combineEpochNames);
useCombined = ~isempty(combineIdx);
effectiveBranchMode = normalizeLegacySpiderBranchMode(p.Results.branchMode);

valuesA = nan(1, nFiles);
valuesC = nan(1, nFiles);
sigA = false(1, nFiles);
sigC = false(1, nFiles);
for i = 1:nFiles
    fmask = (fileID == uFiles(i));
    if useCombined
        xC = normVals(fmask & ismember(epIdx, combineIdx));
        valuesA(i) = 0;
        valuesC(i) = mean(xC, 'omitnan');
        sigA(i) = false;
        sigC(i) = isSignificantVsZero(xC, p.Results.sigAlpha);
    else
        xA = normVals(fmask & epIdx == epochPair(1));
        xC = normVals(fmask & epIdx == epochPair(2));
        valuesA(i) = mean(xA, 'omitnan');
        valuesC(i) = mean(xC, 'omitnan');
        sigA(i) = isSignificantVsZero(xA, p.Results.sigAlpha);
        sigC(i) = isSignificantVsZero(xC, p.Results.sigAlpha);
    end
end
if useCombined
    effectiveBranchMode = 'C';
end

styleIn = mergeBWAPlotStyle(getBWAPlotStyle(), p.Results.style);
styleSp = styleIn;
styleSp.colors.groupA = epochColor(epochPair(1));
styleSp.colors.groupC = epochColor(epochPair(min(2, numel(epochPair))));
styleSp.colors.baseline = [0 0 0];

titleID = sprintf('%s slice log2FC', char(string(datTxt)));
if nFiles < 3
    [h_fig, h_ax] = resolvePlotAxes('figure', p.Results.figure, 'axes', p.Results.axes, 'name', titleID);
    plotLineFallback(h_fig, h_ax, valuesA, valuesC, labels, timeFrameRangeEpochs, epochPair, effectiveBranchMode, styleSp, combineEpochNames, combinedEpochLabel);
    out = struct();
    out.figure = h_fig;
    out.axes = h_ax;
    out.spider = [];
    out.embedded = ~(isempty(p.Results.figure) && isempty(p.Results.axes));
    out.labels = labels;
    out.valuesA = valuesA;
    out.valuesC = valuesC;
    out.fileIDs = cellstr(uFiles);
    out.epochPair = epochPair;
    out.combineEpochNames = cellstr(combineEpochNames(:))';
    out.combinedEpochLabel = combinedEpochLabel;
    return
end

sp = plotPolarsSpiderSubplot(valuesA, valuesC, labels, titleID, ...
    'inputScale', 'log2', ...
    'branchMode', effectiveBranchMode, ...
    'sigMaskA', sigA, ...
    'sigMaskC', sigC, ...
    'sigColorA', styleSp.colors.groupA, ...
    'sigColorC', styleSp.colors.groupC, ...
    'sigAlpha', p.Results.sigAlpha, ...
    'figure', p.Results.figure, ...
    'axes', p.Results.axes, ...
    'style', styleSp);

if ~isempty(sp.axes)
    if isfield(timeFrameRangeEpochs, 'epochNames') && numel(timeFrameRangeEpochs.epochNames) >= max(epochPair)
        eNames = string(timeFrameRangeEpochs.epochNames);
        try
            if useCombined
            hl = legend(sp.axes, {combinedEpochLabel, 'zero'}, ...
                'Location', 'southoutside', 'Interpreter', 'none');
            nudgeLegendRight(hl, 0.06);
        elseif strcmpi(char(string(effectiveBranchMode)), 'both')
                hl = legend(sp.axes, {char(eNames(epochPair(1))), char(eNames(epochPair(2))), 'zero'}, ...
                    'Location', 'southoutside', 'Interpreter', 'none');
                nudgeLegendRight(hl, 0.06);
            else
                hl = legend(sp.axes, {char(eNames(epochPair(2))), 'zero'}, ...
                    'Location', 'southoutside', 'Interpreter', 'none');
                nudgeLegendRight(hl, 0.06);
            end
        catch
        end
    end
end

out = sp;
out.labels = labels;
out.valuesA = valuesA;
out.valuesC = valuesC;
out.sigA = sigA;
out.sigC = sigC;
out.fileIDs = cellstr(uFiles);
out.epochPair = epochPair;
out.combineEpochNames = cellstr(combineEpochNames(:))';
out.combinedEpochLabel = combinedEpochLabel;
end

function c = epochColor(idx)
pal = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13];
idx = max(1, min(size(pal,1), round(double(idx))));
c = pal(idx,:);
end

function mode = normalizeLegacySpiderBranchMode(modeIn)
mode = lower(char(string(modeIn)));
switch mode
    case 'cell_only'
        mode = 'C';
    case 'c'
        mode = 'C';
    case 'both'
        mode = 'both';
    otherwise
        error('plotBAviolinSubplotSliceSpider:BadBranchMode', ...
            'branchMode must be ''cell_only'', ''C'', or ''both''.');
end
end

function labels = makeSliceLabels(fileIDs, filesTbl)
labels = bwaMakeSliceLabels(fileIDs, filesTbl);
end

function plotLineFallback(hFig, hAx, valuesA, valuesC, labels, timeFrameRangeEpochs, epochPair, branchMode, styleSp, combineEpochNames, combinedEpochLabel)
figure(hFig); axes(hAx); cla(hAx);
x = 1:numel(labels);
if strcmpi(char(string(branchMode)), 'both')
    plot(hAx, x, valuesA, '-o', 'LineWidth', 1.8, 'MarkerSize', 5, 'Color', styleSp.colors.groupA); hold(hAx, 'on');
    plot(hAx, x, valuesC, '-o', 'LineWidth', 1.8, 'MarkerSize', 5, 'Color', styleSp.colors.groupC);
    legVals = epochLegend(timeFrameRangeEpochs, epochPair, true, combineEpochNames, combinedEpochLabel);
else
    plot(hAx, x, valuesC, '-o', 'LineWidth', 1.8, 'MarkerSize', 5, 'Color', styleSp.colors.groupC); hold(hAx, 'on');
    legVals = epochLegend(timeFrameRangeEpochs, epochPair, false, combineEpochNames, combinedEpochLabel);
end
yline(hAx, 0, '-', 'Color', styleSp.colors.baseline, 'LineWidth', 2.5);
xticks(hAx, x); xticklabels(hAx, labels); xtickangle(hAx, 30);
ylabel(hAx, 'log2FC', 'Interpreter', 'none');
title(hAx, 'Slice log2FC (fallback: <3 slices)', 'Interpreter', 'none');
grid(hAx, 'on');
hl = legend(hAx, legVals, 'Location', 'southoutside', 'Interpreter', 'none');
nudgeLegendRight(hl, 0.06);
end

function legVals = epochLegend(timeFrameRangeEpochs, epochPair, bothMode, combineEpochNames, combinedEpochLabel)
legVals = {'epoch', 'zero'};
if nargin >= 4 && ~isempty(combineEpochNames)
    legVals = {combinedEpochLabel, 'zero'};
    return
end
if isfield(timeFrameRangeEpochs, 'epochNames') && numel(timeFrameRangeEpochs.epochNames) >= max(epochPair)
    eNames = string(timeFrameRangeEpochs.epochNames);
    if bothMode
        legVals = {char(eNames(epochPair(1))), char(eNames(epochPair(2))), 'zero'};
    else
        legVals = {char(eNames(epochPair(2))), 'zero'};
    end
end
end

function nudgeLegendRight(hLeg, frac)
if isempty(hLeg) || ~isgraphics(hLeg) || ~isfinite(frac), return; end
try
    oldUnits = hLeg.Units;
    hLeg.Units = 'normalized';
    pos = hLeg.Position;
    pos(1) = min(0.98 - pos(3), pos(1) + frac * (1 - pos(3)));
    hLeg.Position = pos;
    hLeg.Units = oldUnits;
catch
end
end

function tf = isSignificantVsZero(x, alpha)
x = double(x(:));
x = x(isfinite(x));
tf = false;
if numel(x) < 3
    return
end
try
    p = signrank(x, 0, 'method', 'approximate');
    tf = isfinite(p) && p < alpha;
catch
    try
        [~, p] = ttest(x, 0);
        tf = isfinite(p) && p < alpha;
    catch
        tf = false;
    end
end
end

function idx = resolveCombinedEpochIdx(timeFrameRangeEpochs, combineEpochNames)
idx = [];
if ~isstruct(timeFrameRangeEpochs) || ~isfield(timeFrameRangeEpochs, 'epochNames') || isempty(timeFrameRangeEpochs.epochNames)
    return
end
epochNames = string(timeFrameRangeEpochs.epochNames(:));
if isempty(combineEpochNames)
    return
end
for i = 1:numel(combineEpochNames)
    hit = find(epochNames == string(combineEpochNames(i)), 1, 'first');
    if ~isempty(hit)
        idx(end+1,1) = hit; %#ok<AGROW>
    end
end
idx = unique(idx, 'stable');
end
