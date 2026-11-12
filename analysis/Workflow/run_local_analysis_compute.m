function result = run_local_analysis_compute(varargin)
% Primary local workflow entrypoint.
% Computes canonical + derived summary data without opening figures.
% Recommended workspace variable name:
%   result
%
% The returned result stores:
%   - one default normalized view, selected by 'normalizationTransform'
%     (default: 'delta'; recommended for new work)
%   - all cached normalized views, so downstream tools can choose per metric
%
% Recommended convention:
%   - leave the default normalized view as 'delta'
%   - use dataView='recommended' in downstream summary/correlation calls
%   - combine repeated post-drug epochs (for example ExperimentalEpoch1/ExperimentalEpoch2 -> CombinedPostEpochs)
%     for spider and correlation analyses
%
% Clear aliases:
%   result.normalized.default_transform
%   result.normalized.canonical.log2_fold / result.normalized.canonical.delta
%   result.normalized.metrics.log2_fold / result.normalized.metrics.delta
%
% Compatibility fields kept in parallel:
%   result.canonical_norm
%   result.metrics
%   result.canonical_norms
%   result.metrics_by_transform
%
% Preferred plotting entrypoints after compute:
%   plot_local_analysis_summary(result, ...)
%   plot_local_analysis_spider(result, ...)
%   plot_local_analysis_legacy_panels(result, ...)
%
% Spider post-epoch combine options:
%   'spiderCombineEpochNames'  : post epochs to combine before spider plotting
%                                (example {'ExperimentalEpoch1','ExperimentalEpoch2'}). Empty means all
%                                epochs after baseline.
%   'spiderCombinedEpochLabel' : label stored with the combined spider payload
%                                (default 'CombinedPostEpochs')

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

p = inputParser;
p.StructExpand = false;
addParameter(p, 'pathIn', pwd, @ischar);
addParameter(p, 'mode', 'preset', @(x) any(strcmpi(string(x), ["preset","epoch","genotype","treatment"])));
addParameter(p, 'presetName', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'epochPreset', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'recordingManifestPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'epochManifestPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'windowManifestPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'analysisSpecPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'figureSpecPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'comparisonSpecPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'statisticsPlanPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'projectManifestDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'validateProjectManifest', false, @islogical);
addParameter(p, 'epochIncludeRoles', {'baseline','response'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'respectInferenceStatus', true, @islogical);
addParameter(p, 'dropZeroRates', true, @islogical);
addParameter(p, 'groupValues', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'metrics', {'area_px2','duration_frames','rise19_frames','fall91_frames','decay_tau_frames','dff_auc','dff_max2','max_prop_speed','perimeter_px','circularity','propagation_distance_px'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'useQC', true, @islogical);
addParameter(p, 'qc_files', table(), @(x) istable(x) || isempty(x));
addParameter(p, 'qc_cells', table(), @(x) istable(x) || isempty(x));
addParameter(p, 'qc_events', table(), @(x) istable(x) || isempty(x));
addParameter(p, 'strictQC', false, @islogical);
addParameter(p, 'runStats', true, @islogical);
addParameter(p, 'pairedEpochStats', true, @islogical);
addParameter(p, 'epochForBetweenGroups', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'normalizationMode', 'by_cell_baseline_mean', @(x) ischar(x) || isstring(x));
addParameter(p, 'normalizationTransform', 'delta', @(x) ischar(x) || isstring(x));
addParameter(p, 'normalizationBaselineEpoch', 1);
addParameter(p, 'useNormalizedForStats', false, @islogical);
addParameter(p, 'spreadThresholdPct', 10, @isnumeric);
addParameter(p, 'propagationThresholdPx', 5, @isnumeric);
addParameter(p, 'escapeThresholdFrac', 0.10, @isnumeric);
addParameter(p, 'coincidenceWindowSec', 3, @isnumeric);
addParameter(p, 'burstBinSec', 5, @isnumeric);
addParameter(p, 'burstMinEvents', 3, @isnumeric);
addParameter(p, 'corrBinSec', 5, @isnumeric);
addParameter(p, 'useBranchOriginOwnership', false, @islogical);
addParameter(p, 'includeEventMotionDerived', false, @islogical);
addParameter(p, 'spiderCombineEpochNames', {}, @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'spiderCombinedEpochLabel', 'CNO_post', @(x) ischar(x) || isstring(x));
addParameter(p, 'dataLimitRules', {}, @(x) iscell(x) || isempty(x));
addParameter(p, 'dataLimitType', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'dataLimitPercentiles', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'saveQCIfCreated', false, @islogical);
addParameter(p, 'exportOutputs', true, @islogical);
addParameter(p, 'exportDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'exportSlim', true, @islogical);
parse(p, varargin{:});

opts = p.Results;
pathIn = opts.pathIn;
mode = lower(char(string(opts.mode)));
opts.dataLimitRules = resolveDataLimitRules(opts);
projectManifestPreflight = preflightSharedAnalysisProject(opts.projectManifestDir, ...
    'Required', logical(opts.validateProjectManifest));
sharedSpecs = loadSharedAnalysisSpecs(opts);

loaded = loadClusterOutputs(pathIn, 'attachLegacyMeta', false, 'attachEventMeta', false);
qc = struct('qc_files', table(), 'qc_cells', table(), 'qc_events', table());
qcMeta = struct('source', 'disabled');
if opts.useQC
    if ~isempty(opts.qc_files) || ~isempty(opts.qc_cells) || ~isempty(opts.qc_events)
        qc.qc_files = opts.qc_files;
        qc.qc_cells = opts.qc_cells;
        qc.qc_events = opts.qc_events;
        qcMeta.source = 'provided';
    else
        [qc, ~, qcMeta] = loadOrCreateQC(pathIn);
        if opts.saveQCIfCreated && isfield(qcMeta, 'source') && strcmpi(qcMeta.source, 'created')
            saveQC(pathIn, qc, 'versioned', true);
        end
    end
else
    qc = createQCTables(loaded.files, loaded.cells, loaded.events, 'defaultKeep', true, 'reason', 'qc_disabled_keep_all');
end

[canonicalOut, epochSpec, epochMeta, presetName] = run_local_canonical_by_mode(pathIn, mode, opts, qc, loaded, false);
if isempty(opts.dataLimitRules) && isfield(canonicalOut, 'dataLimitRules') && iscell(canonicalOut.dataLimitRules)
    opts.dataLimitRules = canonicalOut.dataLimitRules;
end

canonicalOut = addCalibratedCanonicalMetrics(canonicalOut, epochSpec);
branchOwnership = struct('applied', false, 'mode', 'disabled', 'nEvents', 0, ...
    'nBranchOwnersResolved', 0, 'nOwnersReassigned', 0);
if opts.useBranchOriginOwnership
    if localCanonicalHasAuthoritativeBranchOwnership(canonicalOut)
        branchOwnership.applied = true;
        branchOwnership.mode = 'canonical';
        if isfield(canonicalOut, 'subset') && isfield(canonicalOut.subset, 'events') && istable(canonicalOut.subset.events)
            branchOwnership.nEvents = height(canonicalOut.subset.events);
            if ismember('branch_origin_cell_uid', canonicalOut.subset.events.Properties.VariableNames)
                branchOwnership.nBranchOwnersResolved = nnz(string(canonicalOut.subset.events.branch_origin_cell_uid) ~= "");
            end
        end
    else
        branchOwnership.mode = 'missing_canonical';
        warning('run_local_analysis_compute:MissingCanonicalBranchOwnership', ...
            ['Canonical ingest does not contain authoritative branch-origin ownership. ' ...
             'Rerun cluster canonical ingest with the updated cluster repo rather than assigning ownership locally.']);
    end
    if isfield(canonicalOut, 'subset') && isfield(canonicalOut.subset, 'events') && ...
            isfield(canonicalOut.subset, 'event_cell_map') && isfield(canonicalOut.subset, 'files')
        try
            canonicalOut.frequency = computeNEventsPerMinCanonical( ...
                canonicalOut.subset.events, canonicalOut.subset.event_cell_map, canonicalOut.subset.files, epochSpec, ...
                'groupBy', {'file_id'}, 'areaPercentiles', [0 100]);
        catch
        end
    end
end
[ctx, standard, polars, derivedMeta] = tryDerivedCompute(loaded.BWA_AQ2, canonicalOut, epochSpec, opts.dataLimitRules, ...
    opts.includeEventMotionDerived, opts.spiderCombineEpochNames, opts.spiderCombinedEpochLabel);
[canonicalNorms, metricsByTransform, normalizationByTransform, canonicalOut, canonicalNorm, metricsSchema, normMeta] = ...
    buildAllNormalizedViews(canonicalOut, ctx, epochSpec, opts);

statsOut = struct();
if opts.runStats
    if logical(opts.respectInferenceStatus) && strcmpi(string(sharedSpecs.inference_status), "descriptive_only")
        statsOut = struct('status', 'suppressed_by_analysis_spec', ...
            'inference_status', 'descriptive_only', ...
            'message', sharedSpecs.inference_warning);
    else
        statsSource = canonicalOut;
        if logical(opts.useNormalizedForStats) && ~strcmpi(string(opts.normalizationMode), "none")
            statsSource = canonicalNorm;
        end
        statsOut = compute_local_mode_stats(statsSource, mode, opts);
    end
end

result = struct();
result.canonical = canonicalOut;
result.canonical_norm = canonicalNorm;
result.metrics = metricsSchema;
result.metrics_by_transform = metricsByTransform;
result.normalization = normMeta;
result.normalization_by_transform = normalizationByTransform;
result.canonical_norms = canonicalNorms;
result.normalized = struct();
result.normalized.default_transform = lower(char(string(opts.normalizationTransform)));
result.normalized.default_view = struct('canonical', canonicalNorm, 'metrics', metricsSchema);
result.normalized.canonical = canonicalNorms;
result.normalized.metrics = metricsByTransform;
result.normalized.meta = normalizationByTransform;
result.derived = struct();
result.derived.standard = standard;
result.derived.polars = polars;
result.derived.context = ctx;
result.derived.qc = qc;
result.derived.epochMeta = epochMeta;
result.derived.status = derivedMeta.status;
result.derived.error = derivedMeta.error;
if isfield(derivedMeta, 'errors')
    result.derived.errors = derivedMeta.errors;
end
result.metric_qc_report = build_metric_qc_report(canonicalOut, standard);
result.metricVocabulary = buildMetricVocabulary(result);
result.runMode = mode;
result.pathIn = pathIn;
result.ranCanonicalQC = qcMeta.source;
result.stats = statsOut;
result.sharedSpecs = sharedSpecs;
result.projectManifestPreflight = projectManifestPreflight;
result.dataLimitRules = opts.dataLimitRules;
result.spreadThresholdPct = opts.spreadThresholdPct;
result.propagationThresholdPx = opts.propagationThresholdPx;
result.escapeThresholdFrac = opts.escapeThresholdFrac;
result.coincidenceWindowSec = opts.coincidenceWindowSec;
result.burstBinSec = opts.burstBinSec;
result.burstMinEvents = opts.burstMinEvents;
result.corrBinSec = opts.corrBinSec;
result.spiderCombineEpochNames = cellstr(string(opts.spiderCombineEpochNames(:)))';
result.spiderCombinedEpochLabel = char(string(opts.spiderCombinedEpochLabel));
result.branchOwnership = branchOwnership;
if ~isempty(presetName)
    result.presetName = presetName;
end

function vocab = buildMetricVocabulary(result)
ids = strings(0,1);
if isfield(result, 'canonical') && isstruct(result.canonical) && ...
        isfield(result.canonical, 'metrics') && isstruct(result.canonical.metrics)
    ids = [ids; string(fieldnames(result.canonical.metrics))]; %#ok<AGROW>
end
if isfield(result, 'metrics') && isstruct(result.metrics) && ...
        isfield(result.metrics, 'registry') && istable(result.metrics.registry) && ...
        ismember('metric_id', result.metrics.registry.Properties.VariableNames)
    ids = [ids; string(result.metrics.registry.metric_id(:))]; %#ok<AGROW>
end
ids = unique(ids(strlength(ids) > 0), 'stable');
vocab = bwaMetricVocabularyCatalog(ids);
end

function [canonicalNorms, metricsByTransform, normalizationByTransform, canonicalOut, canonicalNorm, metricsSchema, normMeta] = ...
        buildAllNormalizedViews(canonicalOut, ctx, epochSpec, opts)
transforms = unique([{char(string(opts.normalizationTransform))} {'log2_fold','delta'}], 'stable');
canonicalNorms = struct();
metricsByTransform = struct();
normalizationByTransform = struct();

canonicalNorm = struct();
metricsSchema = struct();
normMeta = struct();
selectedTransform = lower(char(string(opts.normalizationTransform)));

for i = 1:numel(transforms)
    tr = lower(char(string(transforms{i})));
    [cn, nm] = normalizeCanonicalOutput(canonicalOut, ...
        'mode', opts.normalizationMode, ...
        'transform', tr, ...
        'baselineEpoch', opts.normalizationBaselineEpoch);
    ms = build_local_metrics_schema(canonicalOut, cn, ctx, epochSpec, ...
        'spreadThresholdPct', opts.spreadThresholdPct, ...
        'propagationThresholdPx', opts.propagationThresholdPx, ...
        'escapeThresholdFrac', opts.escapeThresholdFrac, ...
        'coincidenceWindowSec', opts.coincidenceWindowSec, ...
        'burstBinSec', opts.burstBinSec, ...
        'burstMinEvents', opts.burstMinEvents, ...
        'corrBinSec', opts.corrBinSec, ...
        'normalizationTransform', tr);
    [canonicalOutTmp, cn] = alignCanonicalFrequencyToCellMetric(canonicalOut, cn, ms);
    field = matlab.lang.makeValidName(tr);
    canonicalNorms.(field) = cn;
    metricsByTransform.(field) = ms;
    normalizationByTransform.(field) = nm;
    if strcmp(tr, selectedTransform)
        canonicalOut = canonicalOutTmp;
        canonicalNorm = cn;
        metricsSchema = ms;
        normMeta = nm;
    end
end

if isempty(fieldnames(canonicalNorms))
    [canonicalOut, canonicalNorm] = alignCanonicalFrequencyToCellMetric(canonicalOut, struct(), metricsSchema);
end
end

if opts.exportOutputs
    exportDir = getExportDir(opts.pathIn, mode, opts.exportDir);
    result.exportDir = exportDir;
    if ~exist(exportDir, 'dir')
        mkdir(exportDir);
    end
    exportNorm = bwaWriteNormalizationExportArtifacts(exportDir, result);
    exportCaps = bwaWriteCapabilityExportArtifacts(exportDir, result);
    exportVocab = bwaWriteMetricVocabularyArtifacts(exportDir, result);
    exportFrequency = bwaWritePost2FrequencyExports(exportDir, result);
    result.exportMeta = struct('normalization', exportNorm.bundle, 'capabilities', exportCaps.bundle, ...
        'metricVocabularyRows', height(exportVocab.table), 'frequency', exportFrequency);
    result_full = result;
    result_to_save = result;
    if opts.exportSlim
        result_to_save = slimResultForExport(result_to_save);
    end
    result = result_to_save; %#ok<NASGU>
    save(fullfile(exportDir, 'analysis_result.mat'), 'result', '-v7.3');
    result = result_full;
end
end

function exportDir = getExportDir(pathIn, mode, override)
if strlength(char(string(override))) > 0
    exportDir = char(string(override));
else
    ts = datestr(now, 'yyyymmdd_HHMMSS');
    exportDir = fullfile(char(string(pathIn)), 'derived', 'local_outputs', [mode '_' ts]);
end
end

function [ctx, standard, polars, meta] = tryDerivedCompute(BWA_AQ2, canonicalOut, epochSpec, dataLimitRules, includeEventMotionDerived, spiderCombineEpochNames, spiderCombinedEpochLabel)
ctx = struct();
standard = struct();
polars = struct();
meta = struct('status', 'skipped', 'error', '', 'errors', struct());

if isfield(epochSpec, 'assignmentTable') && istable(epochSpec.assignmentTable)
    meta.status = 'skipped_recording_relative_epochs';
    meta.error = ['Legacy-derived modules require one global frame range per epoch and cannot safely represent ' ...
        'recording-relative experimental activity windows. Canonical and cell-level metrics remain available.'];
    meta.errors.legacy = meta.error;
    return
end

if ~isstruct(BWA_AQ2) || isempty(fieldnames(BWA_AQ2)) || ~isfield(BWA_AQ2, 'fts1')
    meta.status = 'skipped_missing_legacy';
    meta.error = 'BWA_AQ2.mat not available; skipped legacy-derived metrics.';
    meta.errors.legacy = meta.error;
    return
end

try
    ctx = prepareLegacyAqua2Context(BWA_AQ2, ...
        'timeFrameRangeEpochs', epochSpec, ...
        'filesTbl', canonicalOut.subset.files, ...
        'eventsTbl', canonicalOut.subset.events, ...
        'eventCellMapTbl', canonicalOut.subset.event_cell_map);
catch ME
    meta.status = 'error';
    meta.error = ME.message;
    meta.errors.context = ME.message;
    return
end

okStandard = false;
okPolars = false;
try
    standard = computeBWA_Aqua2StandardData(ctx, ...
        'dataLimitRules', dataLimitRules, ...
        'includeEventMotion', includeEventMotionDerived);
    okStandard = true;
catch ME
    standard = struct('status', 'error', 'error', ME.message);
    meta.errors.standard = ME.message;
end

try
    polars = computeBWA_Aqua2PolarsBySliceData(ctx.BWA_AQ2, ...
        'timeFrameRangeEpochs', ctx.timeFrameRangeEpochs, ...
        'memberlist', ctx.memberlistRow3, ...
        'associatedCoverageTraces', ctx.associatedCoverageTraces, ...
        'memberCoverage', ctx.memberCoverage, ...
        'filesTbl', canonicalOut.subset.files, ...
        'eventsTbl', canonicalOut.subset.events, ...
        'dataLimitRules', dataLimitRules, ...
        'combineEpochNames', spiderCombineEpochNames, ...
        'combinedEpochLabel', spiderCombinedEpochLabel, ...
        'includeStandardFeatures', true, ...
        'filePlots', false);
    okPolars = true;
catch ME
    polars = struct('status', 'error', 'error', ME.message);
    meta.errors.polars = ME.message;
end

if okStandard && okPolars
    meta.status = 'ok';
elseif okStandard || okPolars
    meta.status = 'partial';
    meta.error = 'One or more derived modules failed; see derived.error or derived.status.';
else
    meta.status = 'error';
    meta.error = 'All derived modules failed; see derived.status/errors.';
end
end

function shared = loadSharedAnalysisSpecs(opts)
mapping = {
    'analysis_spec', 'analysisSpecPath';
    'figure_spec', 'figureSpecPath';
    'comparison_spec', 'comparisonSpecPath';
    'statistics_plan', 'statisticsPlanPath'};
shared = struct('paths', struct(), 'inference_status', '', 'inference_warning', '');
for ii = 1:size(mapping,1)
    outName = mapping{ii,1}; optionName = mapping{ii,2};
    pathValue = char(string(getOptionOr(opts, optionName, '')));
    shared.paths.(outName) = pathValue;
    shared.(outName) = struct();
    if isempty(strtrim(pathValue)), continue; end
    if ~isfile(pathValue)
        error('run_local_analysis_compute:MissingSharedSpec', '%s does not exist: %s', optionName, pathValue);
    end
    shared.(outName) = jsondecode(fileread(pathValue));
end
if isfield(shared.analysis_spec, 'statistics')
    stats = shared.analysis_spec.statistics;
    if isfield(stats, 'status'), shared.inference_status = char(string(stats.status)); end
    if isfield(stats, 'warning'), shared.inference_warning = char(string(stats.warning)); end
end
if isempty(shared.inference_status) && isfield(shared.statistics_plan, 'inference_status')
    shared.inference_status = char(string(shared.statistics_plan.inference_status));
end
end

function value = getOptionOr(opts, name, defaultValue)
if isstruct(opts) && isfield(opts, name) && ~isempty(opts.(name)), value = opts.(name); else, value = defaultValue; end
end

function rules = resolveDataLimitRules(opts)
rules = {};
if isfield(opts, 'dataLimitRules') && iscell(opts.dataLimitRules) && ~isempty(opts.dataLimitRules)
    rules = opts.dataLimitRules;
    return
end

dt = "";
if isfield(opts, 'dataLimitType') && ~isempty(opts.dataLimitType)
    dt = string(opts.dataLimitType);
end
pp = [];
if isfield(opts, 'dataLimitPercentiles') && ~isempty(opts.dataLimitPercentiles)
    pp = opts.dataLimitPercentiles;
end
if strlength(dt) > 0 && ~isempty(pp) && ~isNoOpPercentileLimit(pp)
    rules = AddLimitingRuleset({}, 'percentiles', char(dt), pp);
end
end

function tf = isNoOpPercentileLimit(pp)
pp = double(pp(:))';
tf = numel(pp) == 2 && pp(1) <= 0 && pp(2) >= 100;
end

function out = slimResultForExport(in)
% Remove large redundant fields while keeping summary plotting compatibility.
out = in;
stripped = {};

if isfield(out, 'canonical') && isstruct(out.canonical)
    if isfield(out.canonical, 'rates')
        out.canonical = rmfield(out.canonical, 'rates');
        stripped{end+1} = 'canonical.rates'; %#ok<AGROW>
    end
    if isfield(out.canonical, 'subset')
        out.canonical = rmfield(out.canonical, 'subset');
        stripped{end+1} = 'canonical.subset'; %#ok<AGROW>
    end
end

if isfield(out, 'canonical_norm') && isstruct(out.canonical_norm)
    if isfield(out.canonical_norm, 'rates')
        out.canonical_norm = rmfield(out.canonical_norm, 'rates');
        stripped{end+1} = 'canonical_norm.rates'; %#ok<AGROW>
    end
    if isfield(out.canonical_norm, 'subset')
        out.canonical_norm = rmfield(out.canonical_norm, 'subset');
        stripped{end+1} = 'canonical_norm.subset'; %#ok<AGROW>
    end
    if isfield(out, 'normalization') && isstruct(out.normalization) && ...
            isfield(out.normalization, 'mode') && strcmpi(string(out.normalization.mode), "none")
        out = rmfield(out, 'canonical_norm');
        stripped{end+1} = 'canonical_norm'; %#ok<AGROW>
    end
end

if isfield(out, 'derived') && isstruct(out.derived) && isfield(out.derived, 'context') && isstruct(out.derived.context)
    cmin = struct();
    if isfield(out.derived.context, 'timeFrameRangeEpochs')
        cmin.timeFrameRangeEpochs = out.derived.context.timeFrameRangeEpochs;
    end
    out.derived.context = cmin;
    stripped{end+1} = 'derived.context.(large fields)'; %#ok<AGROW>
end

if ~isfield(out, 'exportMeta') || ~isstruct(out.exportMeta)
    out.exportMeta = struct();
end
out.exportMeta.slim = true;
out.exportMeta.strippedFields = stripped;
end

function tf = localCanonicalHasAuthoritativeBranchOwnership(canonicalOut)
tf = false;
if ~(isstruct(canonicalOut) && isfield(canonicalOut, 'subset') && isstruct(canonicalOut.subset) && ...
        isfield(canonicalOut.subset, 'events') && istable(canonicalOut.subset.events))
    return
end
E = canonicalOut.subset.events;
need = {'branch_origin_x_px','branch_origin_y_px','branch_origin_frame','branch_origin_cell_uid','primary_cell_uid'};
if ~all(ismember(need, E.Properties.VariableNames))
    return
end
tf = any(isfinite(double(E.branch_origin_frame))) || any(string(E.branch_origin_cell_uid) ~= "");
end

function [canonicalOut, canonicalNorm] = alignCanonicalFrequencyToCellMetric(canonicalOut, canonicalNorm, metricsSchema)
if ~isstruct(metricsSchema) || ~isfield(metricsSchema, 'tables') || ~isstruct(metricsSchema.tables) || ...
        ~isfield(metricsSchema.tables, 'cell_frequency') || ~istable(metricsSchema.tables.cell_frequency) || ...
        isempty(metricsSchema.tables.cell_frequency) || ...
        ~isfield(canonicalOut, 'frequency') || ~isstruct(canonicalOut.frequency)
    return
end

cellFreq = metricsSchema.tables.cell_frequency;
canonicalOut.frequency = adaptCanonicalFrequencyFromCellMetric(cellFreq, canonicalOut.frequency, ...
    'valueColumn', 'rate_events_per_min');

if isstruct(canonicalNorm) && isfield(canonicalNorm, 'frequency') && isstruct(canonicalNorm.frequency) && ...
        ismember('rate_events_per_min_norm', cellFreq.Properties.VariableNames)
    canonicalNorm.frequency = adaptCanonicalFrequencyFromCellMetric(cellFreq, canonicalNorm.frequency, ...
        'valueColumn', 'rate_events_per_min_norm');
end
end
