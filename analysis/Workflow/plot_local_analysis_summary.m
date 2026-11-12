function figs = plot_local_analysis_summary(source, varargin)
% Replot canonical summary figures using a saved result structure/path.
%
% dataView:
%   raw         -> raw values
%   norm        -> the compute run's selected normalizationTransform
%   recommended -> metric-aware normalization
%                  canonical positive magnitudes -> log2_fold
%                  cell/event/recruitment/network or signed metrics -> delta
%
% metricFamilies:
%   optional family filter across canonical, cell, and slice summary surfaces
%   supported values: activity, shape, kinetics, intensity, propagation,
%   containment, recruitment, synchrony, soma
result = load_local_analysis_result(source);
if ~isfield(result, 'canonical')
    error('plot_local_analysis_summary:MissingCanonical', 'Canonical summary output required.');
end
p = inputParser;
p.StructExpand = false;
addParameter(p, 'metricList', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'metricFamilies', 'all', @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'plotFrequencyArgs', {}, @iscell);
addParameter(p, 'plotMetricArgs', {}, @iscell);
addParameter(p, 'includeDerived', true, @islogical);
addParameter(p, 'includeOverlappingDerived', false, @islogical);
addParameter(p, 'comparisonLevel', 'by_slice', @(x) ischar(x) || isstring(x));
addParameter(p, 'cellInclusion', 'complete_cells_only', @(x) ischar(x) || isstring(x));
addParameter(p, 'dataView', 'raw', @(x) ischar(x) || isstring(x)); % raw|norm|recommended
addParameter(p, 'canonicalAggregation', 'cell_promoted', @(x) ischar(x) || isstring(x)); % cell_promoted|event_native
parse(p, varargin{:});
dataView = lower(char(string(p.Results.dataView)));
metricFamilies = normalizeMetricFamiliesLocal(p.Results.metricFamilies);
figs = struct();
figs.provenance = struct();
figs.canonicalAggregation = lower(char(string(p.Results.canonicalAggregation)));

canonical = result.canonical;
if strcmp(dataView, 'norm')
    if isfield(result, 'canonical_norm') && isstruct(result.canonical_norm)
        canonical = result.canonical_norm;
        if isfield(result, 'normalization') && isstruct(result.normalization) && ...
                isfield(result.normalization, 'mode') && strcmpi(string(result.normalization.mode), "none")
            warning('plot_local_analysis_summary:NormModeNone', ...
                'dataView=norm requested, but normalizationMode is ''none'' in this result.');
        end
    else
        warning('plot_local_analysis_summary:NoNormView', ...
            'Requested dataView=norm but canonical_norm is missing; using raw canonical output.');
    end
end
figs.provenance.overall = bwaResolveNormalizationProvenance(result, dataView);
metricArgs = bwaEnsureDefaultNameValueArgs(p.Results.plotMetricArgs, ...
    {'groupLabel','animal_id','comparisonLevel',char(string(p.Results.comparisonLevel)),'cellInclusion',char(string(p.Results.cellInclusion))});
if isfield(canonical, 'frequency')
    if familyMatchesMetricLocal("frequency", metricFamilies)
        freqArgs = bwaEnsureDefaultNameValueArgs(p.Results.plotFrequencyArgs, ...
            {'groupLabel','animal_id','comparisonLevel',char(string(p.Results.comparisonLevel)),'cellInclusion',char(string(p.Results.cellInclusion)),'meanMode','overall','meanColor',[0 0 0]});
        freqProv = bwaResolveNormalizationProvenance(result, dataView, "frequency", "canonical");
        figs.provenance.frequency = freqProv;
        freqArgs = setNameValueArg(freqArgs, 'yLabelSuffix', [' (' char(freqProv.compactLabel) ')']);
        freqOut = canonical.frequency;
        if strcmp(dataView, 'recommended')
            freqOut = selectCanonicalFrequencyForMetric(result, 'frequency');
        end
        figs.frequency = plotNEventsPerMinCanonical(freqOut, freqArgs{:});
    end
end

metrics = fieldnames(canonical.metrics);
metrics = filterPreferredMetricIDs(metrics);
if ~isempty(p.Results.metricList)
    metrics = intersect(metrics, cellstr(string(p.Results.metricList(:))), 'stable');
end
metricArgs = bwaEnsureDefaultNameValueArgs(metricArgs, {'meanMode','overall','meanColor',[0 0 0]});
figs.metrics = struct();
figs.provenance.metrics = struct();
figs.vocabulary = struct('canonical', struct(), 'cell_metrics', struct(), 'slice_metrics', struct());
suppressedPromotedCellIDs = strings(0,1);
for i = 1:numel(metrics)
    metricName = metrics{i};
    if ~familyMatchesMetricLocal(metricName, metricFamilies)
        continue
    end
    [metricOut, aggInfo] = bwaSelectCanonicalMetricSummary(result, metricName, ...
        'dataView', dataView, ...
        'canonicalAggregation', p.Results.canonicalAggregation);
    metricProv = bwaResolveNormalizationProvenance(result, dataView, string(metricName), "canonical");
    metricProv.requestedAggregation = string(p.Results.canonicalAggregation);
    metricProv.effectiveAggregation = string(aggInfo.effectiveAggregation);
    metricProv.promoted_metric_id = string(aggInfo.promoted_metric_id);
    figs.provenance.metrics.(matlab.lang.makeValidName(metricName)) = metricProv;
    figs.vocabulary.canonical.(matlab.lang.makeValidName(metricName)) = bwaLookupMetricVocabulary(metricName);
    metricArgsOne = setNameValueArg(metricArgs, 'yLabelSuffix', [' (' char(metricProv.compactLabel) ')']);
    figs.metrics.(metricName) = plotMetricBySliceCanonical(metricOut, metricArgsOne{:});
    promotedID = bwaCanonicalPromotedMetricID(metricName);
    if strlength(promotedID) > 0
        suppressedPromotedCellIDs(end+1,1) = promotedID; %#ok<AGROW>
    end
end
suppressedPromotedCellIDs = unique(suppressedPromotedCellIDs, 'stable');

figs.derived = struct();
figs.cell_metrics = struct();
figs.slice_metrics = struct();
figs.provenance.derived = bwaResolveNormalizationProvenance(result, dataView);
figs.provenance.cell_metrics = struct();
figs.provenance.slice_metrics = struct();
if p.Results.includeDerived
    if ~(isfield(result, 'derived') && isstruct(result.derived))
        warning('plot_local_analysis_summary:NoDerivedStruct', ...
            'includeDerived=true but result.derived is missing.');
    elseif ~(isfield(result.derived, 'standard') && isstruct(result.derived.standard))
        warning('plot_local_analysis_summary:NoDerivedStandard', ...
            'includeDerived=true but result.derived.standard is missing/non-struct.');
    else
        filesTbl = table();
        if isfield(canonical, 'subset') && isfield(canonical.subset, 'files')
            filesTbl = canonical.subset.files;
        end
        epochNames = {};
        if isfield(result.derived, 'context') && isstruct(result.derived.context) ...
                && isfield(result.derived.context, 'timeFrameRangeEpochs') ...
                && isfield(result.derived.context.timeFrameRangeEpochs, 'epochNames')
            epochNames = result.derived.context.timeFrameRangeEpochs.epochNames;
        end
        normCfg = struct('mode','none','transform','log2_fold','baselineEpochIdx',1);
        if strcmp(dataView, 'norm') && isfield(result, 'normalization') && isstruct(result.normalization)
            if isfield(result.normalization, 'mode'), normCfg.mode = char(string(result.normalization.mode)); end
            if isfield(result.normalization, 'transform'), normCfg.transform = char(string(result.normalization.transform)); end
            if isfield(result.normalization, 'baselineEpochIdx') && ~isempty(result.normalization.baselineEpochIdx)
                normCfg.baselineEpochIdx = double(result.normalization.baselineEpochIdx);
            end
        end
        [figs.derived, diag] = plotDerivedFromStandard(result.derived.standard, filesTbl, epochNames, metricArgs, ...
            p.Results.includeOverlappingDerived, normCfg, metricFamilies);
        if diag.nEntries == 0
            if isfield(result.derived.standard, 'errors') && isstruct(result.derived.standard.errors) && ~isempty(fieldnames(result.derived.standard.errors))
                warning('plot_local_analysis_summary:DerivedEmptyWithErrors', ...
                    'No derived summary entries discovered. standard.errors keys: %s', ...
                    strjoin(fieldnames(result.derived.standard.errors), ', '));
            else
                warning('plot_local_analysis_summary:DerivedEmpty', ...
                    'No plottable derived summary entries discovered in result.derived.standard.');
            end
        end
    end
end

if isfield(result, 'metrics') && isstruct(result.metrics) && ...
        isfield(result.metrics, 'registry') && istable(result.metrics.registry) && ...
        isfield(result.metrics, 'tables') && isstruct(result.metrics.tables)
    epochTable = table();
    if isfield(canonical, 'frequency') && isstruct(canonical.frequency) && isfield(canonical.frequency, 'epochTable')
        epochTable = canonical.frequency.epochTable;
    end
    reg = result.metrics.registry;
    for i = 1:height(reg)
        unitLevel = char(string(reg.unit_level(i)));
        if ~(strcmpi(unitLevel, 'cell') || strcmpi(unitLevel, 'slice'))
            continue
        end
        if strcmpi(string(reg.metric_id(i)), "cell.frequency")
            continue
        end
        if strcmpi(unitLevel, 'slice') && strcmpi(string(p.Results.comparisonLevel), 'by_cell')
            continue
        end
        tname = char(string(reg.source_table(i)));
        if ~isfield(result.metrics.tables, tname)
            continue
        end
        T = result.metrics.tables.(tname);
        if ~istable(T) || isempty(T)
            continue
        end
        metricID = char(string(reg.metric_id(i)));
        if ~familyMatchesMetricLocal(metricID, metricFamilies, reg)
            continue
        end
        if any(string(metricID) == suppressedPromotedCellIDs)
            continue
        end
        if shouldSkipRawMetricID(metricID, reg.metric_id)
            continue
        end
        rawCol = char(string(reg.value_column_raw(i)));
        normCol = char(string(reg.value_column_norm(i)));
        [metricTables, metricReg, rawCol, normCol] = selectMetricSchemaForView(result, reg, i, dataView);
        if isempty(metricTables) || ~isfield(metricTables, tname)
            continue
        end
        T = metricTables.(tname);
        if ~istable(T) || isempty(T)
            continue
        end
        ycol = rawCol;
        if (strcmp(dataView, 'norm') || strcmp(dataView, 'recommended')) && strlength(string(normCol)) > 0 && ismember(normCol, T.Properties.VariableNames)
            ycol = normCol;
        end
        if ~ismember(ycol, T.Properties.VariableNames)
            continue
        end
        try
            if strcmpi(unitLevel, 'cell')
                mOut = bwaBuildMetricOutFromTable(T, ycol, epochTable, 'entityKeys', {'primary_cell_uid'});
                vocabRow = bwaLookupMetricVocabulary(metricID, reg);
                mOut.titleName = resolveMetricDisplayName(metricID, vocabRow);
                metricProv = bwaResolveNormalizationProvenance(result, dataView, string(metricID), "cell");
                figs.provenance.cell_metrics.(matlab.lang.makeValidName(metricID)) = metricProv;
                figs.vocabulary.cell_metrics.(matlab.lang.makeValidName(metricID)) = vocabRow;
                metricArgsOne = setNameValueArg(metricArgs, 'yLabelSuffix', [' (' char(metricProv.compactLabel) ')']);
                figs.cell_metrics.(matlab.lang.makeValidName(metricID)) = plotMetricBySliceCanonical(mOut, metricArgsOne{:});
            else
                mOut = bwaBuildMetricOutFromTable(T, ycol, epochTable);
                vocabRow = bwaLookupMetricVocabulary(metricID, reg);
                mOut.titleName = resolveMetricDisplayName(metricID, vocabRow);
                metricProv = bwaResolveNormalizationProvenance(result, dataView, string(metricID), "slice");
                figs.provenance.slice_metrics.(matlab.lang.makeValidName(metricID)) = metricProv;
                figs.vocabulary.slice_metrics.(matlab.lang.makeValidName(metricID)) = vocabRow;
                metricArgsOne = setNameValueArg(metricArgs, 'yLabelSuffix', [' (' char(metricProv.compactLabel) ')']);
                figs.slice_metrics.(matlab.lang.makeValidName(metricID)) = plotMetricBySliceCanonical(mOut, metricArgsOne{:});
            end
        catch ME
            warning('plot_local_analysis_summary:CellMetricPlotFailed', ...
                'Failed plotting %s from metrics schema: %s', metricID, ME.message);
        end
    end
end
if ~strcmpi(string(p.Results.comparisonLevel), 'by_cell')
    dynSliceIDs = string({ ...
        'slice.frac_primary_events_with_neighbor_recruitment'; ...
        'slice.mean_neighbor_entry_delay_sec'; ...
        'slice.mean_neighbor_peak_overlap_frac_of_cell'; ...
        'slice.mean_neighbor_partitions_per_event'});
    for i = 1:numel(dynSliceIDs)
        metricID = char(dynSliceIDs(i));
        if ~familyMatchesMetricLocal(metricID, metricFamilies)
            continue
        end
        try
            [mOut, ok] = bwaBuildDynamicSliceMetricOut(result, metricID, dataView);
            if ~ok
                continue
            end
            vocabRow = bwaLookupMetricVocabulary(metricID);
            metricProv = bwaResolveNormalizationProvenance(result, dataView, string(metricID), "slice");
            figs.provenance.slice_metrics.(matlab.lang.makeValidName(metricID)) = metricProv;
            figs.vocabulary.slice_metrics.(matlab.lang.makeValidName(metricID)) = vocabRow;
            metricArgsOne = setNameValueArg(metricArgs, 'yLabelSuffix', [' (' char(metricProv.compactLabel) ')']);
            figs.slice_metrics.(matlab.lang.makeValidName(metricID)) = plotMetricBySliceCanonical(mOut, metricArgsOne{:});
        catch ME
            warning('plot_local_analysis_summary:DynamicSliceMetricPlotFailed', ...
                'Failed plotting dynamic slice metric %s: %s', metricID, ME.message);
        end
    end
end

function out = resolveMetricDisplayName(metricID, vocabRow)
out = char(string(metricID));
if istable(vocabRow) && ~isempty(vocabRow) && ismember('public_metric_id', vocabRow.Properties.VariableNames)
    val = string(vocabRow.public_metric_id(1));
    if strlength(val) > 0
        out = char(val);
    end
end
end

if (~isfield(figs, 'frequency') || isempty(figs.frequency)) && ...
        isempty(fieldnames(figs.metrics)) && isempty(fieldnames(figs.derived)) && ...
        isempty(fieldnames(figs.cell_metrics)) && isempty(fieldnames(figs.slice_metrics))
    warning('plot_local_analysis_summary:NoPlots', ...
        'No summary figures were created. Check canonical/derived outputs for non-empty rows.');
end
end

function ids = filterPreferredMetricIDs(ids)
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
    raw = pairs{i,1};
    cal = pairs{i,2};
    if any(strcmp(ids, raw)) && any(strcmp(ids, cal))
        skip(end+1,1) = string(raw); %#ok<AGROW>
    end
end
ids = ids(~ismember(string(ids), skip));
end

function families = normalizeMetricFamiliesLocal(in)
if iscell(in)
    families = lower(string(in(:)))';
else
    families = lower(string(in));
end
families = families(strlength(families) > 0);
if isempty(families)
    families = "all";
end
families = unique(families, 'stable');
end

function tf = familyMatchesMetricLocal(metricID, families, reg)
if nargin < 3
    reg = [];
end
if any(families == "all")
    tf = true;
    return
end
fam = bwaResolveMetricFamily(string(metricID), reg);
tf = any(fam == families);
end

function tf = shouldSkipRawMetricID(metricID, allMetricIDs)
allMetricIDs = string(allMetricIDs(:));
tf = false;
pairs = { ...
    'cell.mean_soma_centroid_dist_px','cell.mean_soma_centroid_dist_um'; ...
    'cell.mean_soma_start_dist_px','cell.mean_soma_start_dist_um'; ...
    'cell.mean_soma_end_dist_px','cell.mean_soma_end_dist_um'; ...
    'cell.mean_perimeter_px','cell.mean_perimeter_um'; ...
    'cell.mean_propagation_distance_px','cell.mean_propagation_distance_um'; ...
    'cell.mean_neighbor_entry_delay_frames','cell.mean_neighbor_entry_delay_sec'};
for i = 1:size(pairs,1)
    if strcmp(metricID, pairs{i,1}) && any(allMetricIDs == string(pairs{i,2}))
        tf = true;
        return
    end
end
end

function freqOut = selectCanonicalFrequencyForMetric(result, metricID)
freqOut = result.canonical.frequency;
tr = bwaPreferredTransformForMetricID(metricID, 'canonical');
if ~isfield(result, 'canonical_norms') || ~isstruct(result.canonical_norms)
    return
end
field = matlab.lang.makeValidName(tr);
if isfield(result.canonical_norms, field) && isfield(result.canonical_norms.(field), 'frequency')
    freqOut = result.canonical_norms.(field).frequency;
elseif isfield(result, 'canonical_norm') && isfield(result.canonical_norm, 'frequency')
    freqOut = result.canonical_norm.frequency;
end
end

function [tablesOut, regOut, rawCol, normCol] = selectMetricSchemaForView(result, reg, rowIdx, dataView)
tablesOut = result.metrics.tables;
regOut = reg;
rawCol = char(string(reg.value_column_raw(rowIdx)));
normCol = char(string(reg.value_column_norm(rowIdx)));
if ~strcmp(dataView, 'recommended')
    return
end
metricID = char(string(reg.metric_id(rowIdx)));
unitLevel = char(string(reg.unit_level(rowIdx)));
tr = bwaPreferredTransformForMetricID(metricID, unitLevel);
if isfield(result, 'metrics_by_transform') && isstruct(result.metrics_by_transform)
    field = matlab.lang.makeValidName(tr);
    if isfield(result.metrics_by_transform, field) && isfield(result.metrics_by_transform.(field), 'tables')
        tablesOut = result.metrics_by_transform.(field).tables;
        if isfield(result.metrics_by_transform.(field), 'registry') && istable(result.metrics_by_transform.(field).registry)
            reg2 = result.metrics_by_transform.(field).registry;
            hit = find(string(reg2.metric_id) == string(metricID), 1, 'first');
            if ~isempty(hit)
                rawCol = char(string(reg2.value_column_raw(hit)));
                normCol = char(string(reg2.value_column_norm(hit)));
                regOut = reg2;
            end
        end
        end
    end
end

function args = setNameValueArg(argsIn, key, value)
args = argsIn;
if ~iscell(args)
    args = {};
end
for i = 1:2:numel(args)
    if strcmpi(char(string(args{i})), char(string(key)))
        args{i+1} = value;
        return
    end
end
args = [args {key, value}];
end

function [out, diag] = plotDerivedFromStandard(S, filesTbl, epochNames, metricArgs, includeOverlappingDerived, normCfg, metricFamilies)
out = struct();
diag = struct('nEntries', 0, 'nPlotted', 0, 'nFailed', 0);
entries = bwaCollectDerivedSummaryEntries(S, includeOverlappingDerived);
if nargin >= 7 && ~any(metricFamilies == "all")
    keep = false(height(entries), 1);
    for i = 1:height(entries)
        keep(i) = familyMatchesMetricLocal(entries.title(i), metricFamilies);
    end
    entries = entries(keep, :);
end
diag.nEntries = height(entries);
for i = 1:height(entries)
    key = char(entries.key(i));
    comp = entries.comp{i};
    titleName = char(entries.title(i));
    try
        mOut = bwaCompToMetricOut(comp, filesTbl, epochNames, titleName, normCfg);
        out.(key) = plotMetricBySliceCanonical(mOut, metricArgs{:});
        diag.nPlotted = diag.nPlotted + 1;
    catch ME
        diag.nFailed = diag.nFailed + 1;
        warning('plot_local_analysis_summary:DerivedPlotFailed', ...
            'Failed derived summary plot for %s: %s', key, ME.message);
    end
end
end
