function out = plot_local_analysis_spider(source, varargin)
% Render summary spider panels from saved or rebuilt polar compute data.
%
% Usage:
%   out = plot_local_analysis_spider(result)
%   out = plot_local_analysis_spider(result, 'metricFamilies', 'shape')
%   out = plot_local_analysis_spider(result, 'metricFamilies', 'recruitment')
%   out = plot_local_analysis_spider(result, 'metricFamilies', 'synchrony')
%   out = plot_local_analysis_spider(result, 'metricFamilies', {'shape','recruitment'})
%
% Name-Value:
%   'spiderInputScale' : '' (auto), 'fold', or 'log2'
%   'branchMode'       : 'C' (default), 'E', or 'both'
%   'comparisonLevel'  : '' (auto), 'by_slice', 'by_event', 'by_cell', or 'by_branch'
%   'aggregateEBy'     : low-level override, 'slice' (default) or 'events'
%   'aggregateCBy'     : low-level override, 'slice' (default) or 'cells'
%   'metricFamilies'   : 'all' (default), or any of:
%                        'activity','shape','kinetics','intensity',
%                        'propagation','containment','recruitment',
%                        'synchrony','soma'
%   plus any plotBWA_Aqua2PolarsBySlice name-value args
%
% Epoch combination:
%   The spider payload summarizes the selected post-drug epochs into one
%   post condition before plotting. This is the recommended mode for
%   repeated post-drug epochs such as ExperimentalEpoch1/ExperimentalEpoch2. By default this means the
%   mean of all epochs after baseline. In the local workflow you should
%   usually make this explicit at compute time with:
%       'spiderCombineEpochNames', {'ExperimentalEpoch1','ExperimentalEpoch2'}
%       'spiderCombinedEpochLabel', 'CombinedPostEpochs'
%
% Branch semantics:
%   E = event-weighted branch
%   C = cell-weighted branch
%   synchrony   : only C is meaningful; synchrony-only requests are coerced to C
%
% Full per-spoke E/C definitions are documented in LOCAL_METRIC_DICTIONARY.md
result = load_local_analysis_result(source);
caps = bwaDescribeResultCapabilities(result);

p = inputParser;
p.StructExpand = false;
p.KeepUnmatched = true;
addParameter(p, 'spiderInputScale', '', @(x) ischar(x) || isstring(x)); % '' => auto
addParameter(p, 'branchMode', 'C', @(x) any(strcmpi(string(x), ["E","C","both"])));
addParameter(p, 'comparisonLevel', '', @(x) any(strcmpi(string(x), ["","by_slice","by_event","by_cell","by_branch"])));
addParameter(p, 'aggregateEBy', 'slice', @(x) any(strcmpi(string(x), ["slice","events"])));
addParameter(p, 'aggregateCBy', 'slice', @(x) any(strcmpi(string(x), ["slice","cells"])));
addParameter(p, 'metricFamilies', 'all', @(x) ischar(x) || isstring(x) || iscell(x));
parse(p, varargin{:});
spiderScale = char(string(p.Results.spiderInputScale));
extraArgs = struct2nv(p.Unmatched);
needByFileEntities = strcmpi(string(p.Results.comparisonLevel), "by_event") || ...
    strcmpi(string(p.Results.comparisonLevel), "by_cell") || ...
    strcmpi(string(p.Results.comparisonLevel), "by_branch") || ...
    strcmpi(string(p.Results.aggregateEBy), "events") || ...
    strcmpi(string(p.Results.aggregateCBy), "cells");

[polData, result, meta] = ensure_local_polars_payload(result, ...
    'requireEntityValuesByFile', needByFileEntities);
polData = augment_polars_with_cell_metrics(result, polData);
polData = filter_spider_metric_families(polData, 'metricFamilies', p.Results.metricFamilies);
branchMode = normalizeSpiderBranchMode(p.Results.branchMode);
if isSynchronyOnlySelection(p.Results.metricFamilies)
    branchMode = 'C';
end
polData = attachSpiderFileLabels(result, polData);
polData = fillSpiderDefaults(polData);

need = {'labels','valuesA','valuesC','uniqueFiles','subplotGrid'};
if ~hasRequiredFields(polData, need)
    throwSpiderCapabilityError('plot_local_analysis_spider:MissingPolarsPayload', caps, meta.message);
end

if isempty(strtrim(spiderScale))
    spiderScale = inferSpiderScale(result, polData);
else
    spiderScale = bwaValidateSpiderScale(spiderScale, polData, 'plot_local_analysis_spider');
end

function polData = fillSpiderDefaults(polData)
if ~isstruct(polData)
    return
end
if (~isfield(polData, 'subplotGrid') || isempty(polData.subplotGrid)) && ...
        isfield(polData, 'uniqueFiles') && ~isempty(polData.uniqueFiles)
    n = numel(polData.uniqueFiles);
    polData.subplotGrid = numSubplotsForN(n);
end
end

passArgs = extraArgs;
if ~bwaHasNameValueArg(passArgs, 'spiderInputScale')
    passArgs = [passArgs {'spiderInputScale', spiderScale}];
end
if ~bwaHasNameValueArg(passArgs, 'branchMode')
    passArgs = [passArgs {'branchMode', branchMode}];
end
if ~bwaHasNameValueArg(passArgs, 'comparisonLevel')
    passArgs = [passArgs {'comparisonLevel', p.Results.comparisonLevel}];
end
if ~bwaHasNameValueArg(passArgs, 'aggregateEBy')
    passArgs = [passArgs {'aggregateEBy', p.Results.aggregateEBy}];
end
if ~bwaHasNameValueArg(passArgs, 'aggregateCBy')
    passArgs = [passArgs {'aggregateCBy', p.Results.aggregateCBy}];
end
if ~bwaHasNameValueArg(passArgs, 'figureName')
    passArgs = [passArgs {'figureName', spiderFigureNameForFamilies(p.Results.metricFamilies)}];
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
        error('plot_local_analysis_spider:BadBranchMode', ...
            'branchMode must be ''E'', ''C'', or ''both''.');
end
end
out = plotBWA_Aqua2PolarsBySlice(polData, passArgs{:});
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
% Spider summaries are baseline-referenced fold/log2 quantities.
% Default to log2 for robust signed visualization.
scale = 'log2';

if isfield(result, 'normalization') && isstruct(result.normalization) && isfield(result.normalization, 'transform')
    tr = lower(char(string(result.normalization.transform)));
    if ~strcmp(tr, 'log2_fold')
        warning('plot_local_analysis_spider:NormTransformMismatch', ...
            ['Spider plots intentionally use fold/log2-style spider payload values, not the workflow''s default normalized view. ' ...
             'This result has normalizationTransform=%s, so spider outputs are not directly comparable to delta-based summary/correlation outputs.'], tr);
    end
end

scale = bwaValidateSpiderScale(scale, polData, 'plot_local_analysis_spider');
end

function polData = attachSpiderFileLabels(result, polData)
if ~isstruct(polData) || ~isfield(polData, 'uniqueFiles') || isempty(polData.uniqueFiles)
    return
end
ids = string(polData.uniqueFiles(:));

F = table();
if isfield(result, 'canonical') && isstruct(result.canonical) && ...
        isfield(result.canonical, 'subset') && isstruct(result.canonical.subset) && ...
        isfield(result.canonical.subset, 'files') && istable(result.canonical.subset.files)
    F = result.canonical.subset.files;
elseif isfield(result, 'canonical') && isstruct(result.canonical) && ...
        isfield(result.canonical, 'files') && istable(result.canonical.files)
    F = result.canonical.files;
end
if isempty(F) && isfield(result, 'pathIn') && ~isempty(result.pathIn)
    try
        L = loadClusterOutputs(char(string(result.pathIn)), 'attachLegacyMeta', false, 'attachEventMeta', false);
        if isfield(L, 'files') && istable(L.files)
            F = L.files;
        end
    catch
    end
end

polData.fileLabels = bwaMakeSliceLabels(ids, F);
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

function name = spiderFigureNameForFamilies(metricFamilies)
name = 'BAPol';
try
    if iscell(metricFamilies)
        groups = string(metricFamilies(:));
    else
        groups = string(metricFamilies);
    end
    groups = groups(strlength(groups) > 0);
    if isempty(groups)
        return
    end
    groups = lower(groups);
    suffix = strjoin(cellstr(groups(:))', '_');
    suffix = regexprep(suffix, '[^a-z0-9_]+', '_');
    if strlength(suffix) > 0
        name = sprintf('BAPol_%s', char(suffix));
    end
catch
end
end

function tf = isSynchronyOnlySelection(metricFamilies)
g = lower(string(metricFamilies));
g(g == "all") = [];
tf = ~isempty(g) && all(g == "synchrony");
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
    detail = "This result is canonical-only and does not include the legacy payload needed to rebuild spider data.";
end
msg = sprintf(['Spider plotting requires either a saved result.derived.polars payload or access to BWA_AQ2.mat for rebuild. ' ...
    '%s Available payloads: %s. ' ...
    'Canonical summary plots remain available.'], detail, char(string(caps.summary)));
if strlength(string(rebuildMessage)) > 0
    msg = sprintf('%s Rebuild detail: %s', msg, char(string(rebuildMessage)));
end
error(errID, '%s', msg);
end
