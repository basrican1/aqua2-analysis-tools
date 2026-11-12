function out = plot_local_analysis_spider_mean(source, varargin)
% Plot a single BAPol-style spider that is the mean across all slices/files.
% Supports per-metric significance vs baseline (0 in log2 space), with
% significance markers on spoke edges and optional spoke error bars.
%
% Usage:
%   out = plot_local_analysis_spider_mean(result)
%   out = plot_local_analysis_spider_mean(result, 'branchMode', 'both')
%   out = plot_local_analysis_spider_mean(result, 'addSigMarkers', true, 'sigAlpha', 0.05)
%
% Name-Value:
%   'spiderInputScale' : '' (auto), 'fold', or 'log2'
%   'branchMode'       : 'C' (default), 'E', or 'both'
%   'comparisonLevel'  : '' (auto), 'by_slice', 'by_event', 'by_cell',
%                        or 'by_branch'
%                        by_slice  = both branches aggregate across slices
%                        by_event  = E branch pools post events; C stays by slice
%                        by_cell   = C branch pools post cells; E stays by slice
%                        by_branch = E pools post events and C pools post cells
%   'metricFamilies'   : 'all' (default), or any of:
%                        'activity','shape','kinetics','intensity',
%                        'propagation','containment','recruitment',
%                        'synchrony','soma'
%   'addSigMarkers'    : true/false, draw '*' at significant spokes (default true)
%   'sigAlpha'         : significance threshold (default 0.05)
%   'sigMethod'        : 'auto' (default), 'signrank', or 'ttest'
%   'sigSource'        : 'auto' (default), 'C', 'E', or 'both'
%   'aggregateEBy'     : low-level override, 'slice' (default) or 'events'
%   'aggregateCBy'     : low-level override, 'slice' (default) or 'cells'
%   'addErrorBars'     : true/false, draw spoke error bars from slice values (default false)
%   'addErrorBand'     : true/false, draw mean+-error uncertainty band (default true)
%   'showEntityPoints' : true/false, overlay pooled event/cell points on spokes (default false)
%   'entityPointSize'  : scatter size for entity points (default 18)
%   'entityPointAlpha' : marker face alpha for entity points (default 0.35)
%   'errorBandAlpha'    : band opacity in [0,1] (default 0.28)
%   'errorBandColor'    : RGB triplet (default [0 0 0])
%   'errorSource'      : 'auto' (default), 'A', 'C', or 'both'
%   'errorType'        : 'sem' (default) or 'sd'
%   'spiderScaleMode'  : 'per_metric' (default) or 'global'
%   'spiderMetricClipPct': percentile for auto limits (default 95)
%   'spiderSymmetricLimit': [] auto, scalar, or 1xN vector
%   plus any plotBWA_Aqua2PolarsBySlice name-value args
%
% Spider branch semantics:
%   E = event-weighted branch. For canonical metrics, E pools events within the
%       slice/file before baseline-vs-post summarization. For cell/recruitment
%       metrics, E uses the closest event-level analog.
%   C = cell-weighted branch. C first computes one summary per cell, then
%       averages those cell summaries within the slice/file.
%   Network note: network metrics are slice-level by construction, so only C is
%       independently meaningful.
%
% Epoch combination:
%   The mean spider uses the same combined post-epoch payload as
%   plot_local_analysis_spider. By default this is the mean of the selected
%   post epochs after baseline. In the local compute workflow you can make
%   this explicit with:
%       'spiderCombineEpochNames', {'ExperimentalEpoch1','ExperimentalEpoch2'}
%       'spiderCombinedEpochLabel', 'CombinedPostEpochs'

result = load_local_analysis_result(source);
caps = bwaDescribeResultCapabilities(result);

p = inputParser;
p.StructExpand = false;
p.KeepUnmatched = true;
addParameter(p, 'spiderInputScale', '', @(x) ischar(x) || isstring(x)); % '' => auto
addParameter(p, 'branchMode', 'C', @(x) any(strcmpi(string(x), ["E","C","both"])));
addParameter(p, 'comparisonLevel', '', @(x) any(strcmpi(string(x), ["","by_slice","by_event","by_cell","by_branch"])));
addParameter(p, 'metricFamilies', 'all', @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'addSigMarkers', true, @islogical);
addParameter(p, 'sigAlpha', 0.05, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x < 1);
addParameter(p, 'sigMethod', 'auto', @(x) any(strcmpi(string(x), ["auto","signrank","ttest"])));
addParameter(p, 'sigSource', 'auto', @(x) any(strcmpi(string(x), ["auto","E","C","both"])));
addParameter(p, 'aggregateEBy', 'slice', @(x) any(strcmpi(string(x), ["slice","events"])));
addParameter(p, 'aggregateCBy', 'slice', @(x) any(strcmpi(string(x), ["slice","cells"])));
addParameter(p, 'addErrorBars', false, @islogical);
addParameter(p, 'addErrorBand', true, @islogical);
addParameter(p, 'showEntityPoints', false, @islogical);
addParameter(p, 'entityPointSize', 18, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'entityPointAlpha', 0.35, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
addParameter(p, 'errorBandAlpha', 0.28, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
addParameter(p, 'errorBandColor', [0 0 0], @(x) isnumeric(x) && numel(x)==3);
addParameter(p, 'errorSource', 'auto', @(x) any(strcmpi(string(x), ["auto","A","C","both"])));
addParameter(p, 'errorType', 'sem', @(x) any(strcmpi(string(x), ["sem","sd"])));
addParameter(p, 'spiderScaleMode', 'per_metric', @(x) any(strcmpi(string(x), ["per_metric","global"])));
addParameter(p, 'spiderMetricClipPct', 95, @isnumeric);
addParameter(p, 'spiderSymmetricLimit', [], @isnumeric);
addParameter(p, 'showStatsFigure', true, @islogical);
addParameter(p, 'statsFigureName', 'BAPol_Mean_Stats', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
spiderScale = char(string(p.Results.spiderInputScale));
extraArgs = struct2nv(p.Unmatched);
branchMode = normalizeSpiderBranchMode(p.Results.branchMode);
if isSynchronyOnlySelection(p.Results.metricFamilies)
    branchMode = 'C';
end
[aggregateEBy, aggregateCBy] = resolveSpiderMeanAggregation( ...
    branchMode, p.Results.comparisonLevel, p.Results.aggregateEBy, p.Results.aggregateCBy);

needEntityValues = strcmpi(aggregateEBy, 'events') || strcmpi(aggregateCBy, 'cells');
[polData, result, meta] = ensure_local_polars_payload(result, 'requireEntityValues', needEntityValues);
polData = augment_polars_with_cell_metrics(result, polData);
polData = filter_spider_metric_families(polData, 'metricFamilies', p.Results.metricFamilies);

need = {'labels','valuesA','valuesC'};
if ~hasRequiredFields(polData, need)
    throwSpiderCapabilityError('plot_local_analysis_spider_mean:MissingPolarsPayload', caps, meta.message);
end

valuesA = double(polData.valuesA);
valuesC = double(polData.valuesC);
if isempty(valuesA) || isempty(valuesC)
    error('plot_local_analysis_spider_mean:EmptyValues', 'Polar values are empty.');
end

if isempty(strtrim(spiderScale))
    spiderScale = inferSpiderScale(result, polData);
else
    spiderScale = bwaValidateSpiderScale(spiderScale, polData, 'plot_local_analysis_spider_mean');
end

valsALog2 = convertToLog2(valuesA, spiderScale);
valsCLog2 = convertToLog2(valuesC, spiderScale);
[branchData, aggInfo] = buildSpiderMeanBranchData(polData, valsALog2, valsCLog2, spiderScale, ...
    aggregateEBy, aggregateCBy);
branchData.labels = relabelSpiderMeanMetrics(branchData.labels, branchMode, resolveComparisonLevelLabel(p.Results.comparisonLevel, aggregateEBy, aggregateCBy));
[branchData, polData] = pruneSpiderMeanBranchData(branchData, polData, branchMode, ...
    resolveComparisonLevelLabel(p.Results.comparisonLevel, aggregateEBy, aggregateCBy));
meanALog2 = branchData.meanA;
meanCLog2 = branchData.meanC;

if ~any(isfinite(meanCLog2)) && ~any(isfinite(meanALog2))
    error('plot_local_analysis_spider_mean:NoFiniteValues', 'No finite values available to plot mean spider.');
end

nMetrics = numel(polData.labels);
lims = [];
if p.Results.addErrorBars || p.Results.addErrorBand || ~isempty(p.Results.spiderSymmetricLimit)
    lims = resolvePanelLimitFromPolars(valuesA, valuesC, spiderScale, p.Results.spiderSymmetricLimit, ...
        nMetrics, p.Results.spiderScaleMode, p.Results.spiderMetricClipPct);
end

sigSource = resolveSigSource(p.Results.sigSource, branchMode);
stats = computeMetricSignificance(branchData, sigSource, p.Results.sigMethod, p.Results.sigAlpha);

meanPol = struct();
meanPol.labels = decorateSpiderMeanLabelsWithUnits(branchData.labels);
meanPol.valuesA = meanALog2;
meanPol.valuesC = meanCLog2;
meanPol.uniqueFiles = "ALL_SLICES_MEAN";
meanPol.subplotGrid = [1 1];
meanPol.fileLabels = {'ALL_SLICES_MEAN'};

passArgs = extraArgs;
if ~bwaHasNameValueArg(passArgs, 'spiderInputScale')
    passArgs = [passArgs {'spiderInputScale', spiderScale}];
end
if ~bwaHasNameValueArg(passArgs, 'branchMode')
    passArgs = [passArgs {'branchMode', branchMode}];
end
if ~bwaHasNameValueArg(passArgs, 'figureName')
    passArgs = [passArgs {'figureName', 'BAPol_Mean'}];
end
if ~bwaHasNameValueArg(passArgs, 'exportCompanionTable')
    passArgs = [passArgs {'exportCompanionTable', false}];
end
if ~bwaHasNameValueArg(passArgs, 'spiderScaleMode')
    passArgs = [passArgs {'spiderScaleMode', p.Results.spiderScaleMode}];
end
if ~bwaHasNameValueArg(passArgs, 'spiderMetricClipPct')
    passArgs = [passArgs {'spiderMetricClipPct', p.Results.spiderMetricClipPct}];
end
if ~isempty(p.Results.spiderSymmetricLimit) && ~bwaHasNameValueArg(passArgs, 'spiderSymmetricLimit')
    passArgs = [passArgs {'spiderSymmetricLimit', p.Results.spiderSymmetricLimit}];
end
if ~bwaHasNameValueArg(passArgs, 'showLegend')
    passArgs = [passArgs {'showLegend', false}];
end

out = plotBWA_Aqua2PolarsBySlice(meanPol, passArgs{:});
out.meanValuesA = meanALog2;
out.meanValuesC = meanCLog2;
out.labels = meanPol.labels;
out.stats = stats;
out.statsSummary = buildStatsSummaryTable(stats.table);
out.statsInfo = struct('baseline_reference', 'log2FC = 0', ...
    'sig_source', upper(char(string(sigSource))), ...
    'comparison_level', string(resolveComparisonLevelLabel(p.Results.comparisonLevel, aggregateEBy, aggregateCBy)), ...
    'aggregate_E_by', string(aggInfo.aggregateEBy), ...
    'aggregate_C_by', string(aggInfo.aggregateCBy), ...
    'requested_method', lower(char(string(p.Results.sigMethod))), ...
    'alpha', p.Results.sigAlpha);
limsPlot = lims;
if isfield(out, 'panelLimitsLog2') && ~isempty(out.panelLimitsLog2)
    limsPlot = double(out.panelLimitsLog2);
end
out.panelLimitsLog2 = limsPlot;
out.statsFigure = [];
if p.Results.showStatsFigure
    out.statsFigure = displayStatsSummaryFigure(out.statsSummary, char(string(p.Results.statsFigureName)), out.statsInfo);
end

ax = findPrimaryAxes(out.figure);
if ~isempty(ax) && isgraphics(ax)
    style = getBWAPlotStyle();
    errSources = resolveSpiderErrorSources(p.Results.errorSource, branchMode);
    bandHandles = gobjects(0);
    errHandles = gobjects(0);
    for ii = 1:numel(errSources)
        src = errSources{ii};
        if strcmp(src, 'A')
            valsLog2 = branchData.distA;
            bandColor = style.colors.groupA;
        else
            valsLog2 = branchData.distC;
            bandColor = style.colors.groupC;
        end
        [muLog2, errLog2] = computeMetricError(valsLog2, p.Results.errorType);
        if p.Results.addErrorBand
            if strcmpi(char(string(p.Results.errorSource)), 'auto') || strcmpi(char(string(p.Results.errorSource)), 'both')
                useColor = bandColor;
            else
                useColor = reshape(double(p.Results.errorBandColor), 1, 3);
            end
            hBand = drawSpiderErrorBand(ax, muLog2, errLog2, limsPlot, ...
                'alpha', p.Results.errorBandAlpha, ...
                'color', useColor);
            bandHandles = [bandHandles; hBand(:)]; %#ok<AGROW>
        end
        if p.Results.addErrorBars
            hErr = drawSpiderErrorBars(ax, muLog2, errLog2, limsPlot);
            errHandles = [errHandles; hErr(:)]; %#ok<AGROW>
        end
    end
    out.errorBand = bandHandles;
    out.errorBars = errHandles;
    out.entityPoints = gobjects(0);

    if p.Results.showEntityPoints
        entitySources = resolveSpiderErrorSources('auto', branchMode);
        pointHandles = gobjects(0);
        for ii = 1:numel(entitySources)
            src = entitySources{ii};
            if strcmp(src, 'A')
                valsLog2 = branchData.distA;
                pointColor = style.colors.groupA;
            else
                valsLog2 = branchData.distC;
                pointColor = style.colors.groupC;
            end
            hPts = drawSpiderEntityPoints(ax, valsLog2, limsPlot, ...
                'color', pointColor, ...
                'markerSize', p.Results.entityPointSize, ...
                'markerAlpha', p.Results.entityPointAlpha);
            pointHandles = [pointHandles; hPts(:)]; %#ok<AGROW>
        end
        out.entityPoints = pointHandles;
    end

    if p.Results.addSigMarkers
        hSig = drawSpiderSignificanceMarkers(ax, stats.significant, nMetrics);
        out.sigMarkers = hSig;
    else
        out.sigMarkers = gobjects(0);
    end
else
    out.errorBand = gobjects(0);
    out.errorBars = gobjects(0);
    out.entityPoints = gobjects(0);
    out.sigMarkers = gobjects(0);
end
out.statsAnnotation = addStatsInfoAnnotation(out.figure, out.statsInfo, stats);
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
        error('plot_local_analysis_spider_mean:BadBranchMode', ...
            'branchMode must be ''E'', ''C'', or ''both''.');
end
end

function sigSource = resolveSigSource(sigSourceIn, branchMode)
sigSource = lower(char(string(sigSourceIn)));
if strcmp(sigSource, 'auto')
    sigSource = lower(char(string(branchMode)));
end
end

function [aggregateEBy, aggregateCBy] = resolveSpiderMeanAggregation(branchMode, comparisonLevel, aggregateEIn, aggregateCIn)
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
        error('plot_local_analysis_spider_mean:BadComparisonLevel', ...
            'comparisonLevel must be '''', ''by_slice'', ''by_event'', ''by_cell'', or ''by_branch''.');
end

if strcmpi(branchMode, 'E') && strcmp(aggregateCBy, 'cells')
    aggregateCBy = 'slice';
elseif strcmpi(branchMode, 'C') && strcmp(aggregateEBy, 'events')
    aggregateEBy = 'slice';
end
end

function label = resolveComparisonLevelLabel(comparisonLevel, aggregateEBy, aggregateCBy)
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

function sources = resolveSpiderErrorSources(errorSourceIn, branchMode)
src = lower(char(string(errorSourceIn)));
if strcmp(src, 'auto')
    switch upper(char(string(branchMode)))
        case 'E'
            src = 'a';
        case 'C'
            src = 'c';
        otherwise
            src = 'both';
    end
end
switch upper(src)
    case 'A'
        sources = {'A'};
    case 'C'
        sources = {'C'};
    otherwise
        sources = {'A','C'};
end
end

function labelsOut = relabelSpiderMeanMetrics(labelsIn, branchMode, comparisonLevel)
labelsOut = cellstr(string(labelsIn(:)))';
mode = upper(char(string(branchMode)));
level = lower(char(string(comparisonLevel)));
useEventLabels = (strcmp(mode, 'E') && strcmp(level, 'by_event')) || ...
    (strcmpi(mode, 'both') && strcmp(level, 'by_branch'));
if ~useEventLabels
    return
end
for i = 1:numel(labelsOut)
    labelsOut{i} = relabelEventModeMetric(labelsOut{i});
end
end

function [branchData, polData] = pruneSpiderMeanBranchData(branchData, polData, branchMode, comparisonLevel)
keep = true(1, numel(branchData.labels));
mode = upper(char(string(branchMode)));
level = lower(char(string(comparisonLevel)));

switch level
    case 'by_event'
        if strcmp(mode, 'E')
            keep = branchData.usedEntityA;
        end
    case 'by_cell'
        if strcmp(mode, 'C')
            keep = branchData.usedEntityC;
        end
    case 'by_branch'
        if strcmp(mode, 'both')
            keep = branchData.usedEntityA & branchData.usedEntityC;
        elseif strcmp(mode, 'E')
            keep = branchData.usedEntityA;
        elseif strcmp(mode, 'C')
            keep = branchData.usedEntityC;
        end
end

labels = string(branchData.labels(:));
if any(keep)
    concept = strings(size(labels));
    for i = 1:numel(labels)
        if keep(i)
            concept(i) = bwaMetricConceptKey(labels(i));
        end
    end
    seen = strings(0,1);
    for i = 1:numel(labels)
        if ~keep(i)
            continue
        end
        key = concept(i);
        if strlength(key) == 0
            key = labels(i);
        end
        if any(seen == key)
            keep(i) = false;
        else
            seen(end+1,1) = key; %#ok<AGROW>
        end
    end
end

branchData = subsetSpiderMeanBranchData(branchData, keep);
polData = subsetSpiderMeanPolData(polData, keep);
end

function branchData = subsetSpiderMeanBranchData(branchData, keep)
branchData.labels = branchData.labels(keep);
branchData.distA = branchData.distA(keep);
branchData.distC = branchData.distC(keep);
branchData.meanA = branchData.meanA(keep);
branchData.meanC = branchData.meanC(keep);
branchData.usedEntityA = branchData.usedEntityA(keep);
branchData.usedEntityC = branchData.usedEntityC(keep);
end

function polData = subsetSpiderMeanPolData(polData, keep)
polData.labels = cellstr(string(polData.labels(:)));
polData.labels = polData.labels(keep);
if isfield(polData, 'valuesA') && size(polData.valuesA, 2) == numel(keep)
    polData.valuesA = polData.valuesA(:, keep);
end
if isfield(polData, 'valuesC') && size(polData.valuesC, 2) == numel(keep)
    polData.valuesC = polData.valuesC(:, keep);
end
for f = ["entityValuesA","entityValuesC"]
    if isfield(polData, char(f)) && iscell(polData.(char(f))) && numel(polData.(char(f))) == numel(keep)
        polData.(char(f)) = polData.(char(f))(keep);
    end
end
for f = ["sigMaskA","sigMaskC","panelLimitsLog2"]
    if isfield(polData, char(f)) && isnumeric(polData.(char(f))) && size(polData.(char(f)), 2) == numel(keep)
        polData.(char(f)) = polData.(char(f))(:, keep);
    end
end
end

function labelOut = relabelEventModeMetric(labelIn)
label = char(string(labelIn));
switch lower(label)
    case 'cell_mean_perimeter_um'
        labelOut = 'perimeter_um';
    case 'cell_mean_circularity'
        labelOut = 'circularity';
    case 'cell_frac_propagative_events'
        labelOut = 'frac_propagative_events';
    case 'cell_mean_propagation_distance_um'
        labelOut = 'propagation_distance_um';
    case 'cell_mean_frac_outside_origin_cell'
        labelOut = 'frac_outside_origin_cell';
    case 'cell_frac_escape_events_gt_threshold'
        labelOut = 'frac_escape_events_gt_threshold';
    case 'cell_frac_events_with_neighbor_recruitment'
        labelOut = 'frac_events_with_neighbor_recruitment';
    case 'cell_mean_neighbor_entry_delay_sec'
        labelOut = 'neighbor_entry_delay_sec';
    case 'cell_mean_neighbor_peak_overlap_frac_of_cell'
        labelOut = 'neighbor_peak_overlap_frac_of_cell';
    case 'cell_mean_neighbor_partitions_per_event'
        labelOut = 'neighbor_partitions_per_event';
    case 'cell_frac_multicell_events'
        labelOut = 'frac_multicell_events';
    case 'cell_mean_cells_per_multicell_event'
        labelOut = 'cells_per_multicell_event';
    case 'cell_frac_soma_events'
        labelOut = 'frac_soma_events';
    case 'cell_frac_toward_events'
        labelOut = 'frac_toward_events';
    case 'cell_frac_away_events'
        labelOut = 'frac_away_events';
    case 'cell_frac_traverse_events'
        labelOut = 'frac_traverse_events';
    case 'cell_frac_net_toward_soma_events'
        labelOut = 'frac_net_toward_soma_events';
    case 'cell_frac_net_away_soma_events'
        labelOut = 'frac_net_away_soma_events';
    case 'cell_mean_soma_net_distance_change_um'
        labelOut = 'soma_net_distance_change_um';
    case 'cell_mean_soma_centroid_dist_um'
        labelOut = 'soma_centroid_dist_um';
    case 'cell_mean_soma_start_dist_um'
        labelOut = 'soma_start_dist_um';
    case 'cell_mean_soma_end_dist_um'
        labelOut = 'soma_end_dist_um';
    otherwise
        if startsWith(lower(label), 'cell_mean_')
            labelOut = extractAfter(label, strlength('cell_mean_'));
        elseif startsWith(lower(label), 'cell_frac_')
            labelOut = ['frac_' char(extractAfter(label, strlength('cell_frac_')))];
        elseif startsWith(lower(label), 'cell_')
            labelOut = char(extractAfter(label, strlength('cell_')));
        else
            labelOut = label;
        end
end
end

function throwSpiderCapabilityError(errID, caps, rebuildMessage)
if nargin < 3
    rebuildMessage = '';
end
if caps.hasPolarsPayload
    detail = "Spider payload exists but is incomplete.";
elseif caps.hasLegacyPayload
    detail = "Legacy payload is available, but spider payload rebuild failed.";
else
    detail = "This result is canonical-only and does not include the legacy payload needed to rebuild mean spider data.";
end
msg = sprintf(['Mean spider plotting requires either a saved result.derived.polars payload or access to BWA_AQ2.mat for rebuild. ' ...
    '%s Available payloads: %s. ' ...
    'Canonical summary plots remain available.'], detail, char(string(caps.summary)));
if strlength(string(rebuildMessage)) > 0
    msg = sprintf('%s Rebuild detail: %s', msg, char(string(rebuildMessage)));
end
error(errID, '%s', msg);
end

function ax = findPrimaryAxes(hFig)
ax = [];
if isempty(hFig) || ~isgraphics(hFig)
    return
end
allAx = findall(hFig, 'Type', 'axes');
if isempty(allAx)
    return
end
for k = 1:numel(allAx)
    a = allAx(k);
    try
        if strcmpi(get(a, 'Visible'), 'on')
            ax = a;
            return
        end
    catch
    end
end
ax = allAx(1);
end

function h = drawSpiderErrorBand(ax, muLog2, errLog2, lims, varargin)
p = inputParser;
p.StructExpand = false;
addParameter(p, 'alpha', 0.28, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
addParameter(p, 'color', [0 0 0], @(x) isnumeric(x) && numel(x)==3);
parse(p, varargin{:});
alphaVal = double(p.Results.alpha);
bandColor = reshape(double(p.Results.color), 1, 3);
muLog2 = reshape(muLog2, 1, []);
errLog2 = reshape(errLog2, 1, []);
lims = reshape(lims, 1, []);
n = numel(muLog2);
if numel(errLog2) ~= n || numel(lims) ~= n
    h = gobjects(0);
    return
end
lo = muLog2 - errLog2;
hi = muLog2 + errLog2;
for j = 1:n
    if isfinite(lims(j)) && lims(j) > 0
        lo(j) = max(lo(j), -lims(j));
        hi(j) = min(hi(j), lims(j));
    end
end
rhoLo = NaN(1,n);
rhoHi = NaN(1,n);
for j = 1:n
    if ~isfinite(lo(j)) || ~isfinite(hi(j)) || ~isfinite(lims(j)) || lims(j) <= 0
        continue
    end
    rhoLo(j) = mapLog2ToRadius(lo(j), lims(j));
    rhoHi(j) = mapLog2ToRadius(hi(j), lims(j));
end
th = arrayfun(@(k) spokeTheta(k, n), 1:n);
valid = isfinite(rhoLo) & isfinite(rhoHi);
if nnz(valid) < 3
    h = gobjects(0);
    return
end
th = th(valid);
rhoLo = rhoLo(valid);
rhoHi = rhoHi(valid);
thWrap = [th, th(1)];
rhoHiWrap = [rhoHi, rhoHi(1)];
rhoLoWrap = [rhoLo, rhoLo(1)];
[xHi, yHi] = pol2cart(thWrap, rhoHiWrap);
[xLo, yLo] = pol2cart(fliplr(thWrap), fliplr(rhoLoWrap));
x = [xHi, xLo];
y = [yHi, yLo];
hold(ax, 'on');
h = patch(ax, x, y, bandColor, 'FaceAlpha', alphaVal, 'EdgeColor', 'none', 'HandleVisibility', 'off');
uistack(h, 'bottom');
hold(ax, 'off');
end

function h = drawSpiderErrorBars(ax, muLog2, errLog2, lims)
muLog2 = reshape(muLog2, 1, []);
errLog2 = reshape(errLog2, 1, []);
lims = reshape(lims, 1, []);
n = numel(muLog2);
if numel(errLog2) ~= n || numel(lims) ~= n
    h = gobjects(0);
    return
end
h = gobjects(0);
hold(ax, 'on');
for j = 1:n
    if ~isfinite(muLog2(j)) || ~isfinite(errLog2(j)) || ~isfinite(lims(j)) || lims(j) <= 0
        continue
    end
    lo = max(muLog2(j) - errLog2(j), -lims(j));
    hi = min(muLog2(j) + errLog2(j), lims(j));
    r1 = mapLog2ToRadius(lo, lims(j));
    r2 = mapLog2ToRadius(hi, lims(j));
    th = spokeTheta(j, n);
    [x, y] = pol2cart([th th], [r1 r2]);
    hh = line(ax, x, y, 'Color', [0 0 0], 'LineWidth', 1.5, 'LineStyle', '-');
    h(end+1) = hh; %#ok<AGROW>
end
hold(ax, 'off');
end

function h = drawSpiderEntityPoints(ax, valsLog2, lims, varargin)
p = inputParser;
p.StructExpand = false;
addParameter(p, 'color', [0 0 0], @(x) isnumeric(x) && numel(x)==3);
addParameter(p, 'markerSize', 18, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'markerAlpha', 0.35, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
parse(p, varargin{:});

[distCell, n] = normalizeMetricDistributions(valsLog2);
lims = reshape(lims, 1, []);
if numel(lims) ~= n
    h = gobjects(0);
    return
end

h = gobjects(0);
hold(ax, 'on');
for j = 1:n
    v = double(distCell{j});
    v = v(isfinite(v));
    if isempty(v) || ~isfinite(lims(j)) || lims(j) <= 0
        continue
    end
    v = max(min(v, lims(j)), -lims(j));
    rho = arrayfun(@(x) mapLog2ToRadius(x, lims(j)), v(:));
    th0 = spokeTheta(j, n);
    th = th0 + localSpokeJitter(numel(rho));
    [x, y] = pol2cart(th, rho(:)');
    hh = scatter(ax, x, y, p.Results.markerSize, ...
        'Marker', 'o', ...
        'MarkerFaceColor', reshape(double(p.Results.color), 1, 3), ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', double(p.Results.markerAlpha), ...
        'HandleVisibility', 'off');
    h(end+1,1) = hh; %#ok<AGROW>
end
hold(ax, 'off');
end

function dth = localSpokeJitter(n)
if n <= 1
    dth = 0;
    return
end
span = deg2rad(4);
dth = linspace(-span, span, n);
end

function h = drawSpiderSignificanceMarkers(ax, sigMask, nMetrics)
sigMask = logical(sigMask(:));
if numel(sigMask) < nMetrics
    sigMask(end+1:nMetrics) = false;
end
h = gobjects(0);
hold(ax, 'on');
for j = 1:nMetrics
    if ~sigMask(j)
        continue
    end
    th = spokeTheta(j, nMetrics);
    [xx, yy] = pol2cart(th, 1.07);
    hh = text(ax, xx, yy, '*', 'Color', [0.85 0 0], 'FontWeight', 'bold', ...
        'FontSize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Interpreter', 'none', 'Clipping', 'off');
    h(end+1) = hh; %#ok<AGROW>
end
hold(ax, 'off');
end

function th = spokeTheta(idx, n)
% Match spider_plot orientation: first spoke at top, clockwise.
th = pi/2 - 2*pi*(idx-1)/n;
end

function r = mapLog2ToRadius(v, lim)
% Match spider_plot scaling: map [-lim, +lim] into [rho_offset, 1].
axesInterval = 3;
axesOffset = 1;
fullInterval = axesInterval + 1;
rhoOffset = axesOffset / fullInterval;
r = (v + lim) / (2*lim);
r = max(0, min(1, r));
r = r * (1 - rhoOffset) + rhoOffset;
end

function [mu, err] = computeMetricError(valsLog2, errType)
[distCell, nMetrics] = normalizeMetricDistributions(valsLog2);
mu = NaN(1, nMetrics);
err = NaN(1, nMetrics);
for j = 1:nMetrics
    v = distCell{j};
    v = v(isfinite(v));
    if isempty(v)
        continue
    end
    mu(j) = mean(v);
    s = std(v, 0);
    if strcmpi(char(string(errType)), 'sd')
        err(j) = s;
    else
        err(j) = s / sqrt(numel(v));
    end
end
end

function [branchData, aggInfo] = buildSpiderMeanBranchData(polData, valsALog2, valsCLog2, spiderScale, aggregateEBy, aggregateCBy)
branchData = buildSpiderMeanBranchDataStruct(polData, valsALog2, valsCLog2, spiderScale, aggregateEBy, aggregateCBy);
aggInfo = struct('aggregateEBy', string(aggregateEBy), 'aggregateCBy', string(aggregateCBy));
end

function branchData = buildSpiderMeanBranchDataStruct(polData, valsALog2, valsCLog2, spiderScale, aggregateEBy, aggregateCBy)
labels = cellstr(string(polData.labels(:)));
[distA, usedEntityA] = metricDistributionsForBranch(polData, valsALog2, spiderScale, 'A', aggregateEBy);
[distC, usedEntityC] = metricDistributionsForBranch(polData, valsCLog2, spiderScale, 'C', aggregateCBy);
branchData = struct();
branchData.labels = labels;
branchData.distA = distA;
branchData.distC = distC;
branchData.meanA = meanByMetricDistribution(distA);
branchData.meanC = meanByMetricDistribution(distC);
branchData.usedEntityA = usedEntityA;
branchData.usedEntityC = usedEntityC;
branchData.aggregateEBy = string(aggregateEBy);
branchData.aggregateCBy = string(aggregateCBy);
end

function [distCell, usedEntity] = metricDistributionsForBranch(polData, fallbackLog2, spiderScale, branchName, aggregateBy)
labels = cellstr(string(polData.labels(:)));
distCell = cell(1, numel(labels));
usedEntity = false(1, numel(labels));
switch upper(char(string(branchName)))
    case 'A'
        entityField = 'entityValuesA';
        useEntity = strcmpi(char(string(aggregateBy)), 'events');
    otherwise
        entityField = 'entityValuesC';
        useEntity = strcmpi(char(string(aggregateBy)), 'cells');
end
entityVals = {};
if isfield(polData, entityField) && iscell(polData.(entityField)) && numel(polData.(entityField)) == numel(labels)
    entityVals = polData.(entityField);
end
for j = 1:numel(labels)
    v = [];
    if useEntity && ~isempty(entityVals)
        raw = entityVals{j};
        if ~isempty(raw)
            v = convertToLog2(double(raw(:)), spiderScale);
            if any(isfinite(v))
                usedEntity(j) = true;
            end
        end
    end
    if isempty(v)
        v = double(fallbackLog2(:, j));
    end
    distCell{j} = v(isfinite(v));
end
end

function mu = meanByMetricDistribution(distCell)
mu = NaN(1, numel(distCell));
for j = 1:numel(distCell)
    v = double(distCell{j});
    v = v(isfinite(v));
    if ~isempty(v)
        mu(j) = mean(v, 'omitnan');
    end
end
end

function [distCell, nMetrics] = normalizeMetricDistributions(valsLog2)
if iscell(valsLog2)
    distCell = valsLog2;
    nMetrics = numel(valsLog2);
    return
end
valsLog2 = double(valsLog2);
if isempty(valsLog2)
    distCell = cell(1, 0);
    nMetrics = 0;
    return
end
nMetrics = size(valsLog2, 2);
distCell = cell(1, nMetrics);
for j = 1:nMetrics
    distCell{j} = valsLog2(:, j);
end
end

function S = computeMetricSignificance(branchData, sigSource, sigMethod, alpha)
labels = cellstr(string(branchData.labels(:)))';
sources = resolveStatSources(sigSource);
rows = table();
sigMaskUnion = false(1, numel(labels));
bySource = struct();

for s = 1:numel(sources)
    src = sources{s};
    if strcmpi(src, 'E')
        distVals = branchData.distA;
    else
        distVals = branchData.distC;
    end

    [distVals, nMetrics] = normalizeMetricDistributions(distVals);
    n = zeros(1, nMetrics);
    pvals = NaN(1, nMetrics);
    methodUsed = strings(1, nMetrics);
    sigMask = false(1, nMetrics);

    for j = 1:nMetrics
        v = double(distVals{j});
        v = v(isfinite(v));
        n(j) = numel(v);
        if n(j) < 2
            methodUsed(j) = "insufficient_n";
            continue
        end

        [pval, used] = oneSampleP(v, char(sigMethod));
        pvals(j) = pval;
        methodUsed(j) = string(used);
        sigMask(j) = isfinite(pval) && (pval < alpha);
    end

    sigMaskUnion = sigMaskUnion | sigMask;
    bySource.(upper(src)) = struct( ...
        'significant', sigMask, ...
        'p_values', pvals, ...
        'n', n, ...
        'method_used', methodUsed);
    Tsrc = table(repmat(string(upper(src)), numel(labels), 1), string(labels(:)), n(:), methodUsed(:), pvals(:), sigMask(:), ...
        'VariableNames', {'branch','metric','n','method_used','p_value','significant'});
    rows = [rows; Tsrc]; %#ok<AGROW>
end

S = struct();
S.alpha = alpha;
S.source = upper(char(string(sigSource)));
S.requested_method = lower(char(sigMethod));
S.significant = sigMaskUnion;
S.table = rows;
S.bySource = bySource;
end

function T = buildStatsSummaryTable(Tlong)
if isempty(Tlong) || ~istable(Tlong)
    T = Tlong;
    return
end
cols = intersect({'branch','metric','n','method_used','p_value','significant'}, Tlong.Properties.VariableNames, 'stable');
Tlong = Tlong(:, cols);
if ~ismember('branch', Tlong.Properties.VariableNames)
    T = Tlong;
    return
end

metricFamily = string(bwaResolveMetricFamily(string(Tlong.metric)));
metricFamily = metricFamily(:);
if numel(metricFamily) ~= height(Tlong)
    metricFamily = repmat("", height(Tlong), 1);
end
Tlong.metric_family = metricFamily;

branches = unique(string(Tlong.branch), 'stable');
if numel(branches) < 2
    T = movevars(Tlong, 'metric_family', 'After', 'metric');
    T = sortSpiderStatsTable(T);
    return
end

metrics = unique(string(Tlong.metric), 'stable');
metricFamilies = string(bwaResolveMetricFamily(metrics));
metricFamilies = metricFamilies(:);
if numel(metricFamilies) ~= numel(metrics)
    metricFamilies = repmat("", numel(metrics), 1);
end
T = table(metrics, metricFamilies, 'VariableNames', {'metric','metric_family'});
valueVars = setdiff(string(Tlong.Properties.VariableNames), ["branch","metric"], 'stable');
for b = 1:numel(branches)
    Tb = Tlong(string(Tlong.branch) == branches(b), :);
    Tb = removevars(Tb, 'branch');
    renameMap = cell(1, 2*numel(valueVars));
    for j = 1:numel(valueVars)
        renameMap{2*j-1} = char(valueVars(j));
        renameMap{2*j} = char(valueVars(j) + "_" + branches(b));
    end
    Tb = renamevars(Tb, renameMap(1:2:end), renameMap(2:2:end));
    T = outerjoin(T, Tb, 'Keys', 'metric', 'MergeKeys', true, 'Type', 'left');
end
if ismember('metric_family_Tb', T.Properties.VariableNames)
    T.metric_family = T.metric_family;
    T = removevars(T, 'metric_family_Tb');
end
T = sortSpiderStatsTable(T);
end

function T = sortSpiderStatsTable(T)
if ~istable(T) || isempty(T) || ~ismember('metric', T.Properties.VariableNames)
    return
end
[order, familiesOrdered] = bwaOrderSpiderMetricsByFamily(string(T.metric));
T = T(order, :);
if ismember('metric_family', T.Properties.VariableNames)
    T.metric_family = familiesOrdered;
end
end

function sources = resolveStatSources(sigSource)
switch upper(char(string(sigSource)))
    case 'BOTH'
        sources = {'E','C'};
    case 'E'
        sources = {'E'};
    otherwise
        sources = {'C'};
end
end

function [pval, used] = oneSampleP(v, method)
method = lower(char(string(method)));
pval = NaN;
used = method;

switch method
    case 'signrank'
        pval = trySignrank(v);
        used = 'signrank';
    case 'ttest'
        pval = tryTtest(v);
        used = 'ttest';
    otherwise
        pval = trySignrank(v);
        used = 'signrank';
        if ~isfinite(pval)
            pval = tryTtest(v);
            used = 'ttest';
        end
end
end

function labelsOut = decorateSpiderMeanLabelsWithUnits(labelsIn)
labelsIn = cellstr(string(labelsIn(:)))';
labelsOut = labelsIn;
for ii = 1:numel(labelsIn)
    labelsOut{ii} = decorateSpiderMeanLabel(labelsIn{ii});
end
end

function out = decorateSpiderMeanLabel(in)
base = abbreviateSpiderMeanBaseLabel(in);
unit = spiderMeanLabelUnit(in);
if strlength(string(unit)) == 0
    out = base;
else
    out = sprintf('%s\n(%s)', base, unit);
end
end

function out = abbreviateSpiderMeanBaseLabel(in)
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

function out = spiderMeanLabelUnit(in)
key = lower(strtrim(char(string(in))));
if startsWith(key, 'cell_')
    key = char(extractAfter(string(key), strlength("cell_")));
elseif startsWith(key, 'slice_')
    key = char(extractAfter(string(key), strlength("slice_")));
end
switch key
    case 'freq'
        out = 'events/min';
    case 'area'
        out = 'um^2';
    case {'duration','rise19','fall91','decaytau'}
        out = 's';
    case 'dffauc'
        out = 'dF/F*s';
    case 'maxpropspeed'
        out = 'um/s';
    case 'propdist'
        out = 'um';
    case 'propangle'
        out = 'deg';
    case 'centroid'
        out = 'um';
    case {'mean_soma_net_distance_change_um','mean_soma_centroid_dist_um', ...
          'mean_soma_start_dist_um','mean_soma_end_dist_um', ...
          'mean_perimeter_um','mean_propagation_distance_um'}
        out = 'um';
    case 'mean_neighbor_entry_delay_sec'
        out = 's';
    otherwise
        out = '';
end
end

function p = trySignrank(v)
p = NaN;
try
    p = signrank(v, 0);
catch
end
end

function p = tryTtest(v)
p = NaN;
try
    [~, p] = ttest(v, 0);
catch
end
end

function y = convertToLog2(x, inputScale)
x = double(x);
scale = lower(char(string(inputScale)));
switch scale
    case 'log2'
        y = x;
    case 'fold'
        y = NaN(size(x));
        m = isfinite(x) & (x > 0);
        y(m) = log2(x(m));
    otherwise
        error('plot_local_analysis_spider_mean:BadScale', 'spiderInputScale must be ''fold'' or ''log2''.');
end
end

function lims = resolvePanelLimitFromPolars(valuesA, valuesC, inputScale, userLimit, nMetrics, scaleMode, clipPct)
if ~isempty(userLimit)
    if isscalar(userLimit)
        lims = repmat(abs(userLimit), 1, nMetrics);
    else
        if numel(userLimit) ~= nMetrics
            error('plot_local_analysis_spider_mean:BadLimitSize', ...
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
    otherwise
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
lims(lims == 0) = 1;
end

function ok = hasRequiredFields(S, fields)
ok = isstruct(S);
if ~ok
    return
end
for i = 1:numel(fields)
    if ~isfield(S, fields{i}) || isempty(S.(fields{i}))
        ok = false;
        return
    end
end
end

function scale = inferSpiderScale(result, polData)
scale = 'log2';
if isfield(result, 'normalization') && isstruct(result.normalization) && isfield(result.normalization, 'transform')
    tr = lower(char(string(result.normalization.transform)));
    if ~strcmp(tr, 'log2_fold')
        warning('plot_local_analysis_spider_mean:NormTransformMismatch', ...
            ['Spider plots intentionally use fold/log2-style spider payload values, not the workflow''s default normalized view. ' ...
             'This result has normalizationTransform=%s, so spider outputs are not directly comparable to delta-based summary/correlation outputs.'], tr);
    end
end
scale = bwaValidateSpiderScale(scale, polData, 'plot_local_analysis_spider_mean');
end

function nv = struct2nv(S)
nv = {};
if ~isstruct(S)
    return
end
fn = fieldnames(S);
for i = 1:numel(fn)
    nv(end+1:end+2) = {fn{i}, S.(fn{i})}; %#ok<AGROW>
end
end

function hFig = displayStatsSummaryFigure(T, figName, statsInfo)
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
hFig = figure('Name', fName, 'NumberTitle', 'off', 'Color', 'w', 'Position', [150 150 900 360]);
ax = axes('Parent', hFig, 'Position', [0 0 1 1], 'Visible', 'off'); %#ok<NASGU>

infoTxt = '';
if nargin >= 3 && isstruct(statsInfo)
    infoTxt = sprintf('Comparison: %s | E: %s | C: %s | Test vs baseline: %s | source: %s | requested: %s | alpha: %.3g', ...
        char(string(statsInfo.comparison_level)), ...
        char(string(statsInfo.aggregate_E_by)), ...
        char(string(statsInfo.aggregate_C_by)), ...
        char(string(statsInfo.baseline_reference)), ...
        char(string(statsInfo.sig_source)), ...
        char(string(statsInfo.requested_method)), ...
        double(statsInfo.alpha));
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

function h = addStatsInfoAnnotation(hFig, statsInfo, stats)
h = [];
if isempty(hFig) || ~isgraphics(hFig) || ~isstruct(statsInfo)
    return
end
try
    nSig = 0;
    nTot = 0;
    if isstruct(stats) && isfield(stats, 'significant')
        m = logical(stats.significant(:));
        nSig = nnz(m);
        nTot = numel(m);
    end
    txt = sprintf('comparison=%s | E=%s | C=%s | one-sample test vs %s | source: %s | method: %s (auto may choose signrank/ttest) | alpha=%.3g | sig metrics: %d/%d', ...
        char(string(statsInfo.comparison_level)), ...
        char(string(statsInfo.aggregate_E_by)), ...
        char(string(statsInfo.aggregate_C_by)), ...
        char(string(statsInfo.baseline_reference)), ...
        char(string(statsInfo.sig_source)), ...
        char(string(statsInfo.requested_method)), ...
        double(statsInfo.alpha), nSig, nTot);
    h = annotation(hFig, 'textbox', [0.01 0.01 0.98 0.04], ...
        'String', txt, 'Interpreter', 'none', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 9, 'Color', [0.2 0.2 0.2]);
catch
end
end

function tf = isSynchronyOnlySelection(metricFamilies)
g = lower(string(metricFamilies));
g(g == "all") = [];
tf = ~isempty(g) && all(g == "synchrony");
end
