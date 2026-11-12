function out = plotPolarsSpiderSubplot(valuesA, valuesC, labels, fileID, varargin)
% Plot one spider/polar subplot.
% Can be used standalone or embedded in a provided figure/axes.
%
% Name-value:
%   'inputScale'    : 'fold' (default) or 'log2'
%   'symmetricLimit': [] auto, scalar, or 1xN vector
%   'branchMode'    : 'C' (default), 'E', or 'both'
%   'figure'        : target figure (optional)
%   'axes'          : target axes (optional)

p = inputParser;
p.StructExpand = false;
addRequired(p, 'valuesA', @isnumeric);
addRequired(p, 'valuesC', @isnumeric);
addRequired(p, 'labels', @(x) iscell(x) || isstring(x));
addRequired(p, 'fileID', @(x) ischar(x) || isstring(x));
addParameter(p, 'inputScale', 'fold', @(x) ischar(x) || isstring(x));
addParameter(p, 'symmetricLimit', [], @isnumeric);
addParameter(p, 'branchMode', 'C', @(x) any(strcmpi(string(x), ["E","C","both"])));
addParameter(p, 'figure', [], @(x) isempty(x) || ishghandle(x, 'figure'));
addParameter(p, 'axes', [], @(x) isempty(x) || ishghandle(x, 'axes'));
addParameter(p, 'style', struct(), @isstruct);
addParameter(p, 'sigMaskA', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
addParameter(p, 'sigMaskC', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
addParameter(p, 'sigColorA', [], @isnumeric);
addParameter(p, 'sigColorC', [], @isnumeric);
addParameter(p, 'sigAlpha', 0.05, @isnumeric);
parse(p, valuesA, valuesC, labels, fileID, varargin{:});
style = mergeBWAPlotStyle(getBWAPlotStyle(), p.Results.style);

valuesA = reshape(valuesA, 1, []);
valuesC = reshape(valuesC, 1, []);
labels = cellstr(string(labels(:)))';
fileID = char(string(fileID));
fileLabel = bwaFormatFileIDLabel(fileID);

[h_fig, h_ax] = resolvePlotAxes('figure', p.Results.figure, 'axes', p.Results.axes, ...
    'name', ['BAPol_' fileLabel]);

spObj = drawSpiderIntoAxes(h_fig, h_ax, valuesA, valuesC, labels, fileID, ...
    'inputScale', p.Results.inputScale, ...
    'symmetricLimit', p.Results.symmetricLimit, ...
    'branchMode', p.Results.branchMode, ...
    'sigMaskA', p.Results.sigMaskA, ...
    'sigMaskC', p.Results.sigMaskC, ...
    'sigColorA', p.Results.sigColorA, ...
    'sigColorC', p.Results.sigColorC, ...
    'sigAlpha', p.Results.sigAlpha, ...
    'style', style);
set(h_ax, 'FontSize', style.fontSizeAxis);

out = struct();
out.figure = h_fig;
out.axes = h_ax;
out.spider = spObj;
out.embedded = ~(isempty(p.Results.figure) && isempty(p.Results.axes));
end

function spObj = drawSpiderIntoAxes(hFig, hAx, valuesA, valuesC, labels, fileID, varargin)
p = inputParser;
p.StructExpand = false;
addParameter(p, 'inputScale', 'fold', @(x) ischar(x) || isstring(x));
addParameter(p, 'symmetricLimit', [], @isnumeric);
addParameter(p, 'branchMode', 'C', @(x) any(strcmpi(string(x), ["E","C","both"])));
addParameter(p, 'style', struct(), @isstruct);
addParameter(p, 'sigMaskA', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
addParameter(p, 'sigMaskC', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
addParameter(p, 'sigColorA', [], @isnumeric);
addParameter(p, 'sigColorC', [], @isnumeric);
addParameter(p, 'sigAlpha', 0.05, @isnumeric);
parse(p, varargin{:});
style = p.Results.style;
labelsAxes = abbreviateSpiderLabels(labels);
branchMode = normalizeSpiderBranchMode(p.Results.branchMode);

% Convert to log2-domain for plotting.
inputScale = lower(char(string(p.Results.inputScale)));
switch inputScale
    case 'fold'
        log2A = nan(size(valuesA));
        log2C = nan(size(valuesC));
        mA = isfinite(valuesA) & valuesA > 0;
        mC = isfinite(valuesC) & valuesC > 0;
        log2A(mA) = log2(valuesA(mA));
        log2C(mC) = log2(valuesC(mC));
    case 'log2'
        log2A = valuesA;
        log2C = valuesC;
    otherwise
        error('plotPolarsSpiderSubplot:BadScale', 'inputScale must be ''fold'' or ''log2''.');
end

nMetrics = numel(labels);
lims = resolveSymmetricLimit([log2A(:); log2C(:)], p.Results.symmetricLimit, nMetrics);
mins = -lims;
maxs = lims;

% spider_plot/polyshape requires finite numeric coordinates.
% Replace any non-finite values with baseline-equivalent zero in log2 space.
badC = false(size(log2C));
badA = false(size(log2A));
if strcmp(branchMode, 'both') || strcmp(branchMode, 'E')
    badA = ~isfinite(log2A);
end
if strcmp(branchMode, 'both') || strcmp(branchMode, 'C')
    badC = ~isfinite(log2C);
end
if any(badA) || any(badC)
    nBad = nnz(badA) + nnz(badC);
    warning('plotPolarsSpiderSubplot:NonFiniteReplaced', ...
        'File %s has %d non-finite spider values; replacing with 0 for plotting.', fileID, nBad);
    log2A(badA) = 0;
    log2C(badC) = 0;
end

if strcmp(branchMode, 'both')
    P = [log2A; log2C; zeros(1, nMetrics)];
    fillOpt = {'on', 'on', 'off'};
    fillAlpha = [0.2, 0.2, 1];
    lineW = [style.lineWidthMain, style.lineWidthMain, style.lineWidthBaseline];
    colors = [style.colors.groupA; style.colors.groupC; style.colors.baseline];
    markers = {'o', 'o', 'none'};
    markerSz = [style.markerSize, style.markerSize, 1];
elseif strcmp(branchMode, 'E')
    P = [log2A; zeros(1, nMetrics)];
    fillOpt = {'on', 'off'};
    fillAlpha = [0.2, 1];
    lineW = [style.lineWidthMain, style.lineWidthBaseline];
    colors = [style.colors.groupA; style.colors.baseline];
    markers = {'o', 'none'};
    markerSz = [style.markerSize, 1];
else
    P = [log2C; zeros(1, nMetrics)];
    fillOpt = {'on', 'off'};
    fillAlpha = [0.2, 1];
    lineW = [style.lineWidthMain, style.lineWidthBaseline];
    colors = [style.colors.groupC; style.colors.baseline];
    markers = {'o', 'none'};
    markerSz = [style.markerSize, 1];
end
if ~isreal(P)
    warning('plotPolarsSpiderSubplot:ComplexValues', ...
        'File %s has complex spider values; using real part only.', fileID);
    P = real(P);
end
if any(~isfinite(P(:)))
    warning('plotPolarsSpiderSubplot:ResidualNonFinite', ...
        'File %s still has non-finite values after cleanup; forcing to 0.', fileID);
    P(~isfinite(P)) = 0;
end
if any(~isfinite(mins)) || any(~isfinite(maxs))
    mins = -ones(1, nMetrics);
    maxs = ones(1, nMetrics);
end

% Clamp plotted values to spoke limits to avoid spider_plot shading errors.
for jj = 1:nMetrics
    P(:,jj) = min(max(P(:,jj), mins(jj)), maxs(jj));
end

figure(hFig); axes(hAx); cla(hAx);
fileLabel = bwaFormatFileIDLabel(fileID);
try
    spObj = spider_plot(P, ...
        'AxesLabels', labelsAxes, ...
        'AxesLimits', [mins; maxs], ...
        'AxesHandle', hAx, ...
        'AxesWebType', 'circular', ...
        'FillOption', fillOpt, ...
        'FillTransparency', fillAlpha, ...
        'LineWidth', lineW, ...
        'Color', colors, ...
        'Marker', markers, ...
        'MarkerSize', markerSz, ...
        'AxesTickFormat', {'%.2f'});
catch ME
    warning('plotPolarsSpiderSubplot:FillRetry', ...
        'File %s spider fill failed (%s). Retrying with FillOption off.', fileID, ME.message);
    spObj = spider_plot(P, ...
        'AxesLabels', labelsAxes, ...
        'AxesLimits', [mins; maxs], ...
        'AxesHandle', hAx, ...
        'AxesWebType', 'circular', ...
        'FillOption', repmat({'off'}, 1, size(P,1)), ...
        'LineWidth', lineW, ...
        'Color', colors, ...
        'Marker', markers, ...
        'MarkerSize', markerSz, ...
        'AxesTickFormat', {'%.2f'});
end
title(hAx, fileLabel, 'Interpreter', 'none');
try
    tt = get(hAx, 'Title');
    set(tt, 'Units', 'normalized');
    pos = get(tt, 'Position');
    pos(2) = max(pos(2), 1.06);
    set(tt, 'Position', pos);
catch
end
% Ensure underscores are rendered literally for all spider axis text.
forceTextInterpreterNone(hFig, hAx);
plotSignificanceStars(hAx, mins, maxs, p.Results.sigMaskA, p.Results.sigMaskC, ...
    p.Results.sigColorA, p.Results.sigColorC, style, branchMode);
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
        error('plotPolarsSpiderSubplot:BadBranchMode', ...
            'branchMode must be ''E'', ''C'', or ''both''.');
end
end

function lims = resolveSymmetricLimit(vals, userLim, nMetrics)
if isempty(userLim)
    vals = vals(isfinite(vals));
    if isempty(vals)
        lims = ones(1, nMetrics);
        return
    end
    mx = max(abs(vals));
    if mx <= 0
        mx = 1;
    end
    mx = ceil(mx * 2) / 2;
    lims = repmat(mx, 1, nMetrics);
elseif isscalar(userLim)
    lims = repmat(abs(userLim), 1, nMetrics);
else
    lims = abs(reshape(userLim, 1, []));
    if numel(lims) ~= nMetrics
        error('plotPolarsSpiderSubplot:BadLimitSize', ...
            'symmetricLimit vector must match number of metrics.');
    end
end
lims(lims == 0) = 1;
end


function forceTextInterpreterNone(hFig, hAx)
% spider_plot often creates TeX text labels; force literal underscores.
th = findall(hAx, 'Type', 'text');
for k = 1:numel(th)
    try
        set(th(k), 'Interpreter', 'none');
        set(th(k), 'HorizontalAlignment', 'center');
    catch
    end
end
end

function plotSignificanceStars(hAx, mins, maxs, sigA, sigC, colorA, colorC, style, branchMode)
nMetrics = numel(maxs);
if isempty(sigA), sigA = false(1, nMetrics); end
if isempty(sigC), sigC = false(1, nMetrics); end
sigA = logical(reshape(sigA, 1, []));
sigC = logical(reshape(sigC, 1, []));
if numel(sigA) < nMetrics, sigA(end+1:nMetrics) = false; end
if numel(sigC) < nMetrics, sigC(end+1:nMetrics) = false; end
sigA = sigA(1:nMetrics);
sigC = sigC(1:nMetrics);
if isempty(colorA), colorA = style.colors.groupA; end
if isempty(colorC), colorC = style.colors.groupC; end

% Place stars in normalized axes coordinates for stable distance to outer ring.
cx = 0.5;
cy = 0.5;
rBase = 0.405;  % very near the outer ring
offT = 0.010;   % slight tangential split if both branches are significant
offR = 0;
hold(hAx, 'on');
starSize = max(style.fontSizeAxis + 8, round(style.fontSizeAxis * 1.75));
for j = 1:nMetrics
    theta = pi/2 - 2*pi*(j-1)/nMetrics;
    ux = cos(theta); uy = sin(theta);
    tx = -sin(theta); ty = cos(theta);
    if (strcmp(branchMode, 'both') || strcmp(branchMode, 'E')) && sigA(j)
        xA = (rBase + offR)*ux + offT*tx;
        yA = (rBase + offR)*uy + offT*ty;
        text(hAx, cx + xA, cy + yA, '*', 'Color', colorA, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Interpreter', 'none', 'Clipping', 'off', 'Units', 'normalized', ...
            'FontSize', starSize);
    end
    if (strcmp(branchMode, 'both') || strcmp(branchMode, 'C')) && sigC(j)
        xC = (rBase + offR)*ux - offT*tx;
        yC = (rBase + offR)*uy - offT*ty;
        text(hAx, cx + xC, cy + yC, '*', 'Color', colorC, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Interpreter', 'none', 'Clipping', 'off', 'Units', 'normalized', ...
            'FontSize', starSize);
    end
end
hold(hAx, 'off');
end
% Keep interpreter override local to this spider axes only (avoid clobbering LaTeX labels in other subplots).

function labelsOut = abbreviateSpiderLabels(labelsIn)
labelsIn = cellstr(string(labelsIn(:)))';
labelsOut = labelsIn;
for ii = 1:numel(labelsIn)
    labelsOut{ii} = abbreviateSpiderLabel(labelsIn{ii});
end
end

function out = abbreviateSpiderLabel(in)
if contains(string(in), newline)
    out = char(string(in));
    return
end
key = lower(strtrim(char(string(in))));
key = regexprep(key, '^(cell|slice|event)[._]', '');
switch key
    case 'freq'
        out = 'freq';
    case 'area'
        out = 'area';
    case 'duration'
        out = 'duration';
    case 'rise19'
        out = 'rise19';
    case 'fall91'
        out = 'fall91';
    case 'decaytau'
        out = 'decay tau';
    case 'dffauc'
        out = 'dF/F AUC';
    case 'dffmax2'
        out = 'dF/F max';
    case 'maxpropspeed'
        out = 'max prop spd';
    case 'footprint'
        out = 'footprint';
    case 'propdist'
        out = 'legacy prop dist';
    case 'propangle'
        out = 'prop angle';
    case 'centroid'
        out = 'centroid travel';
    case 'frac_multicell_events'
        out = 'multicell frac';
    case 'mean_cells_per_multicell_event'
        out = 'cells/event';
    case 'frac_soma_events'
        out = 'soma frac';
    case 'frac_toward_events'
        out = 'toward frac';
    case 'frac_away_events'
        out = 'away frac';
    case 'frac_traverse_events'
        out = 'traverse frac';
    case 'frac_net_toward_soma_events'
        out = 'net toward';
    case 'frac_net_away_soma_events'
        out = 'net away';
    case 'mean_soma_net_distance_change_um'
        out = 'net soma Δd';
    case 'mean_soma_centroid_dist_um'
        out = 'centroid-soma d';
    case 'mean_soma_start_dist_um'
        out = 'start-soma d';
    case 'mean_soma_end_dist_um'
        out = 'end-soma d';
    case 'mean_perimeter_um'
        out = 'perimeter';
    case 'mean_circularity'
        out = 'circularity';
    case 'frac_propagative_events'
        out = 'frac propagative';
    case 'mean_propagation_distance_um'
        out = 'prop dist';
    case 'mean_frac_outside_origin_cell'
        out = 'outside frac';
    case 'frac_escape_events_gt_threshold'
        out = 'escape frac';
    case 'frac_events_with_neighbor_recruitment'
        out = 'frac recruited';
    case 'mean_neighbor_entry_delay_sec'
        out = 'nbr delay';
    case 'mean_neighbor_peak_overlap_frac_of_cell'
        out = 'nbr peak';
    case 'mean_neighbor_partitions_per_event'
        out = 'nbr parts/event';
    case 'frac_active_cells'
        out = 'frac active';
    case 'coincidence_rate'
        out = 'coincidence';
    case 'burst_fraction'
        out = 'burst frac';
    case 'mean_pairwise_cell_correlation'
        out = 'pair corr';
    case 'frac_primary_events_with_neighbor_recruitment'
        out = 'frac recruited';
    otherwise
        out = strrep(key, '_', ' ');
        out = regexprep(out, '\\bmean\\b', '');
        out = regexprep(out, '\\bfrac\\b', 'frac');
        out = strtrim(regexprep(out, '\\s+', ' '));
    end
end
