function out = compute_local_summary_significance(source, varargin)
% Compute epoch-vs-baseline significance tables without making plots.
%
% Preferred use after run_local_analysis_compute(...):
%   S = compute_local_summary_significance(result, ...)
%
% Key options:
%   'dataView'                  : 'raw' | 'norm'
%   'comparisonLevel'           : 'overall' | 'by_animal' | 'by_slice' | 'by_cell'
%   'cellInclusion'             : 'complete_cells_only' | 'available_case'
%   'metricList'                : subset of canonical metric names
%   'includeDerived'            : include derived summary comps
%   'includeOverlappingDerived' : include derived metrics that overlap canonical
%   'groupLabel'                : grouping variable override
%   'traceLabel'                : trace variable override
%   'stat'                      : 'mean' | 'median'
%   'sigAlpha'                  : significance alpha
%   'sigMethod'                 : 'auto' | 'paired' | 'unpaired'
%   'metricFamilies'            : optional family filter across canonical,
%                                 cell, and slice significance surfaces

result = load_local_analysis_result(source);
if ~isfield(result, 'canonical')
    error('compute_local_summary_significance:MissingCanonical', ...
        'Canonical summary output required.');
end

p = inputParser;
p.StructExpand = false;
addParameter(p, 'metricList', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'metricFamilies', 'all', @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'includeDerived', true, @islogical);
addParameter(p, 'includeOverlappingDerived', false, @islogical);
addParameter(p, 'comparisonLevel', 'by_slice', @(x) ischar(x) || isstring(x));
addParameter(p, 'cellInclusion', 'complete_cells_only', @(x) ischar(x) || isstring(x));
addParameter(p, 'dataView', 'raw', @(x) ischar(x) || isstring(x));
addParameter(p, 'groupLabel', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'traceLabel', 'auto', @(x) ischar(x) || isstring(x));
addParameter(p, 'stat', 'mean', @(x) ischar(x) || isstring(x));
addParameter(p, 'sigAlpha', 0.05, @isnumeric);
addParameter(p, 'sigMethod', 'auto', @(x) ischar(x) || isstring(x));
addParameter(p, 'canonicalAggregation', 'cell_promoted', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});

dataView = lower(char(string(p.Results.dataView)));
metricFamilies = normalizeMetricFamiliesLocal(p.Results.metricFamilies);
canonical = result.canonical;
if strcmp(dataView, 'norm')
    if isfield(result, 'canonical_norm') && isstruct(result.canonical_norm)
        canonical = result.canonical_norm;
    else
        warning('compute_local_summary_significance:NoNormView', ...
            'Requested dataView=norm but canonical_norm is missing; using raw canonical output.');
    end
end

commonArgs = {'comparisonLevel', char(string(p.Results.comparisonLevel)), ...
    'cellInclusion', char(string(p.Results.cellInclusion)), ...
    'groupLabel', char(string(p.Results.groupLabel)), ...
    'traceLabel', char(string(p.Results.traceLabel)), ...
    'stat', char(string(p.Results.stat)), ...
    'sigAlpha', p.Results.sigAlpha, ...
    'sigMethod', char(string(p.Results.sigMethod))};

out = struct();
out.source = 'compute_local_summary_significance';
out.dataView = dataView;
out.canonicalAggregation = lower(char(string(p.Results.canonicalAggregation)));
out.normalizationProvenance = bwaResolveNormalizationProvenance(result, dataView);
out.metricVocabulary = struct('canonical', struct(), 'cell_metrics', struct(), 'slice_metrics', struct());
out.comparisonLevel = char(string(p.Results.comparisonLevel));
out.cellInclusion = char(string(p.Results.cellInclusion));
out.stat = char(string(p.Results.stat));
out.sigAlpha = p.Results.sigAlpha;
out.sigMethod = char(string(p.Results.sigMethod));
out.frequency = struct();
out.metrics = struct();
out.derived = struct();
out.cell_metrics = struct();
out.slice_metrics = struct();

if isfield(canonical, 'frequency') && isstruct(canonical.frequency)
    if familyMatchesMetricLocal("frequency", metricFamilies)
        out.frequency = summarizeFrequencyEpochVsBaseline(canonical.frequency, commonArgs{:});
    end
end

metrics = fieldnames(canonical.metrics);
if ~isempty(p.Results.metricList)
    metrics = intersect(metrics, cellstr(string(p.Results.metricList(:))), 'stable');
end
suppressedPromotedCellIDs = strings(0,1);
for i = 1:numel(metrics)
    metricName = metrics{i};
    if ~familyMatchesMetricLocal(metricName, metricFamilies)
        continue
    end
    out.metricVocabulary.canonical.(matlab.lang.makeValidName(metricName)) = bwaLookupMetricVocabulary(metricName);
    [metricOut, aggInfo] = bwaSelectCanonicalMetricSummary(result, metricName, ...
        'dataView', dataView, ...
        'canonicalAggregation', p.Results.canonicalAggregation);
    out.metrics.(metricName) = summarizeMetricEpochVsBaseline(metricOut, commonArgs{:});
    out.metrics.(metricName).requestedAggregation = string(p.Results.canonicalAggregation);
    out.metrics.(metricName).effectiveAggregation = string(aggInfo.effectiveAggregation);
    out.metrics.(metricName).promoted_metric_id = string(aggInfo.promoted_metric_id);
    promotedID = bwaCanonicalPromotedMetricID(metricName);
    if strlength(promotedID) > 0
        suppressedPromotedCellIDs(end+1,1) = promotedID; %#ok<AGROW>
    end
end
suppressedPromotedCellIDs = unique(suppressedPromotedCellIDs, 'stable');

if p.Results.includeDerived && isfield(result, 'derived') && isstruct(result.derived) && ...
        isfield(result.derived, 'standard') && isstruct(result.derived.standard)
    filesTbl = table();
    if isfield(canonical, 'subset') && isfield(canonical.subset, 'files')
        filesTbl = canonical.subset.files;
    end
    epochNames = {};
    if isfield(result.derived, 'context') && isstruct(result.derived.context) && ...
            isfield(result.derived.context, 'timeFrameRangeEpochs') && ...
            isfield(result.derived.context.timeFrameRangeEpochs, 'epochNames')
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
    entries = bwaCollectDerivedSummaryEntries(result.derived.standard, p.Results.includeOverlappingDerived);
    for i = 1:height(entries)
        key = char(entries.key(i));
        titleName = char(entries.title(i));
        if ~familyMatchesMetricLocal(titleName, metricFamilies)
            continue
        end
        comp = entries.comp{i};
        try
            mOut = bwaCompToMetricOut(comp, filesTbl, epochNames, titleName, normCfg);
            out.derived.(key) = summarizeMetricEpochVsBaseline(mOut, commonArgs{:});
        catch ME
            out.derived.(key) = struct('status', 'error', 'error', ME.message);
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
        rawCol = char(string(reg.value_column_raw(i)));
        normCol = char(string(reg.value_column_norm(i)));
        ycol = rawCol;
        if strcmp(dataView, 'norm') && strlength(string(normCol)) > 0 && ismember(normCol, T.Properties.VariableNames)
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
                out.metricVocabulary.cell_metrics.(matlab.lang.makeValidName(metricID)) = vocabRow;
                out.cell_metrics.(matlab.lang.makeValidName(metricID)) = summarizeMetricEpochVsBaseline(mOut, commonArgs{:});
            else
                mOut = bwaBuildMetricOutFromTable(T, ycol, epochTable);
                vocabRow = bwaLookupMetricVocabulary(metricID, reg);
                mOut.titleName = resolveMetricDisplayName(metricID, vocabRow);
                out.metricVocabulary.slice_metrics.(matlab.lang.makeValidName(metricID)) = vocabRow;
                out.slice_metrics.(matlab.lang.makeValidName(metricID)) = summarizeMetricEpochVsBaseline(mOut, commonArgs{:});
            end
        catch ME
            if strcmpi(unitLevel, 'cell')
                out.cell_metrics.(matlab.lang.makeValidName(metricID)) = struct('status', 'error', 'error', ME.message);
            else
                out.slice_metrics.(matlab.lang.makeValidName(metricID)) = struct('status', 'error', 'error', ME.message);
            end
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
            out.metricVocabulary.slice_metrics.(matlab.lang.makeValidName(metricID)) = vocabRow;
            out.slice_metrics.(matlab.lang.makeValidName(metricID)) = summarizeMetricEpochVsBaseline(mOut, commonArgs{:});
        catch ME
            out.slice_metrics.(matlab.lang.makeValidName(metricID)) = struct('status', 'error', 'error', ME.message);
        end
    end
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

function out = resolveMetricDisplayName(metricID, vocabRow)
out = char(string(metricID));
if istable(vocabRow) && ~isempty(vocabRow) && ismember('public_metric_id', vocabRow.Properties.VariableNames)
    val = string(vocabRow.public_metric_id(1));
    if strlength(val) > 0
        out = char(val);
    end
end
end
end
