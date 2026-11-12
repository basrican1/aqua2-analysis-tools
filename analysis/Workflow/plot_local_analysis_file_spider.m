function out = plot_local_analysis_file_spider(source, fileIDs, varargin)
% Plot spider panels for individual file(s) using saved or rebuilt polar data.
%
% Usage:
%   out = plot_local_analysis_file_spider(result, fileIDs)
%   out = plot_local_analysis_file_spider(result, fileIDs, 'metricFamilies', 'shape')
%   out = plot_local_analysis_file_spider(result, fileIDs, 'metricFamilies', 'recruitment')
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
%   Spider branch semantics:
%   E = event-weighted branch; C = cell-weighted branch. For network metrics,
%   only C is independently meaningful because the source metrics are
%   slice-level by construction.
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
if (~isfield(polData, 'subplotGrid') || isempty(polData.subplotGrid)) && ...
        isfield(polData, 'uniqueFiles') && ~isempty(polData.uniqueFiles)
    polData.subplotGrid = numSubplotsForN(numel(polData.uniqueFiles));
end
need = {'labels','valuesA','valuesC','uniqueFiles','subplotGrid'};
if ~hasRequiredFields(polData, need)
    throwSpiderCapabilityError('plot_local_analysis_file_spider:MissingPolarsPayload', caps, meta.message);
end
if nargin < 2 || isempty(fileIDs)
    fileIDs = polData.uniqueFiles(1);
end
subset = subsetPolarsDataLocal(polData, cellstr(string(fileIDs(:))));

if isempty(strtrim(spiderScale))
    spiderScale = inferSpiderScale(result, subset);
else
    spiderScale = bwaValidateSpiderScale(spiderScale, subset, 'plot_local_analysis_file_spider');
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
out = plotBWA_Aqua2PolarsBySlice(subset, passArgs{:});
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
        warning('plot_local_analysis_file_spider:NormTransformMismatch', ...
            ['Spider plots intentionally use fold/log2-style spider payload values, not the workflow''s default normalized view. ' ...
             'This result has normalizationTransform=%s, so spider outputs are not directly comparable to delta-based summary/correlation outputs.'], tr);
    end
end
scale = bwaValidateSpiderScale(scale, polData, 'plot_local_analysis_file_spider');
end

function polSub = subsetPolarsDataLocal(polData, fileIDs)
mask = ismember(cellstr(string(polData.uniqueFiles(:))), cellstr(string(fileIDs(:))));
if ~any(mask)
    error('plot_local_analysis_file_spider:NoMatch', 'No matching file IDs found in polar data.');
end
polSub = struct();
fields = fieldnames(polData);
for i = 1:numel(fields)
    f = fields{i};
    val = polData.(f);
    if ischar(val) || isstring(val)
        polSub.(f) = val;
        continue
    end
    if iscell(val) && numel(val) == numel(polData.uniqueFiles)
        polSub.(f) = val(mask);
        continue
    end
    if isnumeric(val) && size(val, 1) == numel(polData.uniqueFiles)
        polSub.(f) = val(mask, :);
        continue
    end
    if iscell(val) && numel(val) == numel(polData.labels)
        keepNested = false;
        subVal = val;
        for j = 1:numel(val)
            if iscell(val{j}) && numel(val{j}) == numel(polData.uniqueFiles)
                subVal{j} = val{j}(mask);
                keepNested = true;
            end
        end
        if keepNested
            polSub.(f) = subVal;
            continue
        end
    end
    polSub.(f) = val;
end
polSub.uniqueFiles = polData.uniqueFiles(mask);
polSub.subplotGrid = [1 numel(polSub.uniqueFiles)];
end

function polData = attachSpiderFileLabels(result, polData)
if ~isstruct(polData) || ~isfield(polData, 'uniqueFiles') || isempty(polData.uniqueFiles)
    return
end
ids = string(polData.uniqueFiles(:));
labels = repmat("", numel(ids), 1);

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

for i = 1:numel(ids)
    lbl = "";
    if ~isempty(F) && all(ismember({'file_id','exp_date','slice_num'}, F.Properties.VariableNames))
        idx = find(string(F.file_id) == ids(i), 1, 'first');
        if ~isempty(idx)
            d = regexprep(char(string(F.exp_date(idx))), '[^0-9]', '');
            if numel(d) >= 8
                d = d(end-5:end);
            end
            s = regexprep(char(string(F.slice_num(idx))), '[^0-9]', '');
            if isempty(s)
                s = char(string(F.slice_num(idx)));
            end
            if ~isempty(d) && ~isempty(s)
                lbl = string(sprintf('%s_S%s', d, s));
            end
        end
    end
    if strlength(lbl) == 0
        alt = string(bwaFormatFileIDLabel(char(ids(i))));
        if alt == string(ids(i))
            alt = "S" + string(i);
        end
        lbl = alt;
    end
    labels(i) = lbl;
end

for i = 1:numel(labels)
    if i == 1
        continue
    end
    if any(labels(1:i-1) == labels(i))
        k = 2;
        cand = labels(i) + "_" + string(k);
        while any(labels(1:i-1) == cand)
            k = k + 1;
            cand = labels(i) + "_" + string(k);
        end
        labels(i) = cand;
    end
end

polData.fileLabels = cellstr(labels);
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
        error('plot_local_analysis_file_spider:BadBranchMode', ...
            'branchMode must be ''E'', ''C'', or ''both''.');
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
    detail = "This result is canonical-only and does not include the legacy payload needed to rebuild file-level spider data.";
end
msg = sprintf(['File spider plotting requires either a saved result.derived.polars payload or access to BWA_AQ2.mat for rebuild. ' ...
    '%s Available payloads: %s. ' ...
    'Canonical summary plots remain available.'], detail, char(string(caps.summary)));
if strlength(string(rebuildMessage)) > 0
    msg = sprintf('%s Rebuild detail: %s', msg, char(string(rebuildMessage)));
end
error(errID, '%s', msg);
end
