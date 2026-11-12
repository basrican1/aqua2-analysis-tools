function [sp, companionTbl, figHandles] = plotSpiderLog2FC(valuesA, valuesC, labels, fileIDs, varargin)
% Plot spider charts in log2-fold space with a zero-reference ring.
%
% Inputs:
%   valuesA, valuesC : [nFiles x nMetrics] values in fold or log2 scale
%   labels           : metric labels (cellstr/string)
%   fileIDs          : file IDs (cellstr/string)
%
% Options:
%   'inputScale'     : 'fold' (default) or 'log2'
%   'symmetricLimit' : [] auto, scalar, or 1xN vector for +/- axis limits
%
% Outputs:
%   sp           : array of spider_plot_class objects
%   companionTbl : long table with fold/log2 values per file+metric+group
%   figHandles   : temporary figure handles used for tile composition

p = inputParser;
p.StructExpand = false;
addRequired(p, 'valuesA', @isnumeric);
addRequired(p, 'valuesC', @isnumeric);
addRequired(p, 'labels', @(x) iscell(x) || isstring(x));
addRequired(p, 'fileIDs', @(x) iscell(x) || isstring(x));
addParameter(p, 'inputScale', 'fold', @(x) ischar(x) || isstring(x));
addParameter(p, 'symmetricLimit', [], @isnumeric);
addParameter(p, 'style', struct(), @isstruct);
parse(p, valuesA, valuesC, labels, fileIDs, varargin{:});
style = mergeBWAPlotStyle(getBWAPlotStyle(), p.Results.style);

if ~isequal(size(valuesA), size(valuesC))
    error('plotSpiderLog2FC:SizeMismatch', 'valuesA and valuesC must have same size.');
end
[nFiles, nMetrics] = size(valuesA);

labels = cellstr(string(labels(:)))';
if numel(labels) ~= nMetrics
    error('plotSpiderLog2FC:LabelCountMismatch', 'labels must match number of metric columns.');
end
fileIDs = cellstr(string(fileIDs(:)));
if numel(fileIDs) ~= nFiles
    error('plotSpiderLog2FC:FileCountMismatch', 'fileIDs must match number of rows.');
end

inputScale = lower(char(string(p.Results.inputScale)));
switch inputScale
    case 'fold'
        foldA = valuesA;
        foldC = valuesC;
        badA = isfinite(foldA) & foldA <= 0;
        badC = isfinite(foldC) & foldC <= 0;
        if any(badA(:)) || any(badC(:))
            warning('plotSpiderLog2FC:NonPositiveFold', ...
                'Detected %d non-positive fold values; corresponding log2 values set to NaN.', ...
                nnz(badA) + nnz(badC));
        end
        log2A = log2Safe(foldA);
        log2C = log2Safe(foldC);
    case 'log2'
        log2A = valuesA;
        log2C = valuesC;
        foldA = 2.^log2A;
        foldC = 2.^log2C;
    otherwise
        error('plotSpiderLog2FC:BadScale', 'inputScale must be ''fold'' or ''log2''.');
end

lims = resolveSymmetricLimits(log2A, log2C, p.Results.symmetricLimit, nMetrics);
mins = -lims;
maxs = lims;

sp = repmat(spider_plot_class([]), 1, nFiles);
figHandles = gobjects(0);
for ii = 1:nFiles
    thisA = log2A(ii, :);
    thisC = log2C(ii, :);
    baseline = zeros(1, nMetrics);
    P = [thisA; thisC; baseline];

    figHandles(ii) = figure; %#ok<AGROW>
    sp(ii) = spider_plot_class(P, ...
        'AxesLabels', labels, ...
        'AxesLimits', [mins; maxs], ...
        'AxesWebType', 'circular', ...
        'FillOption', {'on', 'on', 'off'}, ...
        'FillTransparency', [0.2, 0.2, 1], ...
        'LineWidth', [style.lineWidthMain, style.lineWidthMain, style.lineWidthBaseline], ...
        'Color', [style.colors.groupA; style.colors.groupC; style.colors.baseline], ...
        'Marker', {'o', 'o', 'none'}, ...
        'MarkerSize', [style.markerSize, style.markerSize, 1], ...
        'AxesTickFormat', {'%.2f'});
    title(bwaFormatFileIDLabel(fileIDs{ii}), 'Interpreter', 'none');
end

companionTbl = makeCompanionTable(fileIDs, labels, foldA, foldC, log2A, log2C);
end

function out = log2Safe(in)
out = NaN(size(in));
pos = in > 0 & isfinite(in);
out(pos) = log2(in(pos));
end

function lims = resolveSymmetricLimits(a, c, userLim, nMetrics)
if isempty(userLim)
    allVals = [a(:); c(:)];
    allVals = allVals(isfinite(allVals));
    if isempty(allVals)
        lims = ones(1, nMetrics);
        return
    end
    mx = max(abs(allVals));
    if mx <= 0
        mx = 1;
    end
    mx = ceil(mx * 2) / 2; % round up to nearest 0.5
    lims = repmat(mx, 1, nMetrics);
elseif isscalar(userLim)
    lims = repmat(abs(userLim), 1, nMetrics);
else
    if numel(userLim) ~= nMetrics
        error('plotSpiderLog2FC:BadLimitSize', 'symmetricLimit vector must match number of metrics.');
    end
    lims = abs(reshape(userLim, 1, []));
end

lims(lims == 0) = 1;
end

function T = makeCompanionTable(fileIDs, labels, foldA, foldC, log2A, log2C)
nFiles = numel(fileIDs);
nMetrics = numel(labels);
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
        groupCol(r) = "A";
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
end
