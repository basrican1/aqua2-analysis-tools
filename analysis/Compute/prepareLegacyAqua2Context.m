function ctx = prepareLegacyAqua2Context(BWA_AQ2, varargin)
% Prepare shared legacy context used by plotting workflows.
% This function performs metadata/event attachment and epoch/member setup.

    p = inputParser;
    p.StructExpand = false;
    addRequired(p, 'BWA_AQ2', @isstruct);
    addParameter(p,'timeFrameRangeEpochs', [], @(x) isempty(x) || isstruct(x));
    addParameter(p,'timeRange', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p,'frameRange', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p,'EpochList', {}, @(x) isempty(x) || iscellstr(x) || isstring(x));
    addParameter(p,'memberlist', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p,'associatedCoverageTraces', {}, @(x) isempty(x) || iscell(x));
    addParameter(p,'memberCoverage', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p,'filesTbl', table(), @(x) istable(x) || isempty(x));
    addParameter(p,'eventsTbl', table(), @(x) istable(x) || isempty(x));
    addParameter(p,'eventCellMapTbl', table(), @(x) istable(x) || isempty(x));
    parse(p, BWA_AQ2, varargin{:});

    if isfield(BWA_AQ2, 'opts') && isfield(BWA_AQ2, 'DB')
        opts = BWA_AQ2.opts;
        DB = BWA_AQ2.DB;
    else
        if isempty(p.Results.filesTbl)
            error('BWA_AQ2 missing opts/DB. Provide canonical files table via ''filesTbl''.');
        end
        [opts, DB] = buildLegacyMetaFromCanonical(p.Results.filesTbl, BWA_AQ2.fts1);
    end

    if ~(isfield(BWA_AQ2,'fts1') && isfield(BWA_AQ2.fts1,'loc') && isfield(BWA_AQ2.fts1.loc,'t0') && isfield(BWA_AQ2.fts1,'fileID'))
        if isempty(p.Results.eventsTbl)
            error('BWA_AQ2 missing t0/t1/fileID. Provide canonical events table via ''eventsTbl''.');
        end
        BWA_AQ2 = attachEventFieldsFromCanonical(BWA_AQ2, p.Results.eventsTbl, p.Results.filesTbl);
    end

    % When a canonical file subset is provided, reduce the legacy object to
    % those files before any downstream metric computation so single-file QC
    % results and full-result per-file slices are derived from the same
    % underlying event set.
    if istable(p.Results.filesTbl) && ~isempty(p.Results.filesTbl) && ...
            ismember('file_id', p.Results.filesTbl.Properties.VariableNames) && ...
            isfield(BWA_AQ2, 'fts1') && isfield(BWA_AQ2.fts1, 'fileID')
        keepFiles = string(p.Results.filesTbl.file_id(:));
        haveFiles = unique(string(BWA_AQ2.fts1.fileID(:)), 'stable');
        if ~isempty(keepFiles) && (~isequal(haveFiles(:), keepFiles(:)) || numel(haveFiles) ~= numel(keepFiles))
            BWA_AQ2 = ReduceAQ2ByFile(BWA_AQ2, cellstr(keepFiles), ...
                'filesTbl', p.Results.filesTbl, ...
                'eventsTbl', p.Results.eventsTbl);
            if isfield(BWA_AQ2, 'opts') && isfield(BWA_AQ2, 'DB')
                opts = BWA_AQ2.opts;
                DB = BWA_AQ2.DB;
            end
        end
    end

    if isempty(p.Results.timeFrameRangeEpochs)
        framesIn = p.Results.frameRange;
        timesIn = p.Results.timeRange;
        epochNames = p.Results.EpochList;
        if isempty(epochNames)
            epochNames = '';
        end
        sz_t = opts.sz(1,4);
        frameRate = opts.frameRate(1);
        timeFrameRangeEpochs = frameTimeConversionsByEpoch(framesIn, timesIn, sz_t, frameRate, epochNames);
    else
        timeFrameRangeEpochs = p.Results.timeFrameRangeEpochs;
    end

    if isempty(p.Results.memberlist) || isempty(p.Results.associatedCoverageTraces) || isempty(p.Results.memberCoverage)
        [memberlist, memberCoverage, ~, associatedCoverageTraces] = eventBelongsCell(BWA_AQ2.fts1, opts, DB);
    else
        memberlist = p.Results.memberlist;
        associatedCoverageTraces = p.Results.associatedCoverageTraces;
        memberCoverage = p.Results.memberCoverage;
    end

    eventsTbl = p.Results.eventsTbl;
    if istable(eventsTbl) && ~isempty(eventsTbl)
        BWA_AQ2 = fillLegacyMetricsFromEvents(BWA_AQ2, eventsTbl, p.Results.filesTbl);
    end

    memberlistRow3 = memberlist;
    if size(memberlist, 1) == 3
        memberlistRow3 = memberlist(3,:);
    end

    ctx = struct();
    ctx.BWA_AQ2 = BWA_AQ2;
    ctx.opts = opts;
    ctx.DB = DB;
    ctx.timeFrameRangeEpochs = timeFrameRangeEpochs;
    ctx.memberlist = memberlist;
    ctx.memberlistRow3 = memberlistRow3;
    ctx.memberCoverage = memberCoverage;
    ctx.associatedCoverageTraces = associatedCoverageTraces;
    ctx.filesTbl = p.Results.filesTbl;
    ctx.eventsTbl = p.Results.eventsTbl;
    ctx.eventCellMapTbl = p.Results.eventCellMapTbl;
end

function BWA_AQ2 = fillLegacyMetricsFromEvents(BWA_AQ2, eventsTbl, filesTbl)
if ~isfield(BWA_AQ2, 'fts1')
    return
end
E = canonicalizeEvents(eventsTbl, BWA_AQ2.fts1);
E = addLegacyCalibratedEventColumns(E, filesTbl);
if isempty(E)
    return
end

mapping = { ...
    'basic.area', 'area_um2'; ...
    'curve.duration', 'duration_sec'; ...
    'curve.rise19', 'rise19_sec'; ...
    'curve.fall91', 'fall91_sec'; ...
    'curve.decayTau', 'decay_tau_sec'; ...
    'curve.dffAUC', 'dff_auc_sec'; ...
    'curve.dffMax2', 'dff_max2'; ...
    'propagation.maxPropSpeed', 'max_prop_speed'};

for ii = 1:size(mapping,1)
    path = mapping{ii,1};
    col = mapping{ii,2};
    if ~ismember(col, E.Properties.VariableNames)
        continue
    end
    val = rowNum(E.(col));
    if isempty(val)
        continue
    end
    % Prefer the calibrated canonical event values whenever they exist.
    % Legacy fts1 fields may already be populated in older/raw units, which
    % causes single-file rebuilds and full-result payloads to diverge.
    BWA_AQ2.fts1 = setNestedField(BWA_AQ2.fts1, path, val);
end
end

function E = canonicalizeEvents(eventsTbl, fts1)
if ~istable(eventsTbl) || isempty(eventsTbl)
    E = table();
    return
end
E = eventsTbl;
if nargin < 2 || ~isstruct(fts1) || ~isfield(fts1, 'fileID') || isempty(fts1.fileID)
    return
end
if ~all(ismember({'file_id','local_event_id'}, E.Properties.VariableNames))
    return
end

fileIDs = string(fts1.fileID(:));
if numel(fileIDs) ~= height(E)
    return
end

rowOrder = zeros(height(E), 1);
cursor = 0;
uFiles = unique(fileIDs, 'stable');
for i = 1:numel(uFiles)
    fid = uFiles(i);
    idxFts = find(fileIDs == fid);
    idxEvt = find(string(E.file_id) == fid);
    if numel(idxEvt) ~= numel(idxFts)
        return
    end
    [~, ord] = sort(double(E.local_event_id(idxEvt)));
    idxEvt = idxEvt(ord);
    rowOrder(cursor + (1:numel(idxEvt))) = idxEvt;
    cursor = cursor + numel(idxEvt);
end
if all(rowOrder > 0)
    E = E(rowOrder, :);
end
end

function v = rowNum(v)
if iscell(v)
    v = cellfun(@(x) double(x(1)), v);
else
    v = double(v(:));
end
v = v(:)';
end

function [exists, value] = nestedField(S, path)
parts = strsplit(path, '.');
value = S;
exists = true;
for i = 1:numel(parts)
    if ~isstruct(value) || ~isfield(value, parts{i})
        exists = false;
        value = [];
        return
    end
    value = value.(parts{i});
end
end

function S = setNestedField(S, path, val)
parts = strsplit(path, '.');
if numel(parts) == 1
    S.(parts{1}) = val;
    return
end
p = parts{1};
if ~isfield(S, p) || ~isstruct(S.(p))
    S.(p) = struct();
end
S.(p) = setNestedField(S.(p), strjoin(parts(2:end), '.'), val);
end

function E = addLegacyCalibratedEventColumns(E, filesTbl)
if ~istable(E) || isempty(E) || ~istable(filesTbl) || isempty(filesTbl)
    return
end
fileKeep = intersect({'file_id','frame_rate_s','spatial_res'}, filesTbl.Properties.VariableNames);
if ~ismember('file_id', fileKeep)
    return
end
E = outerjoin(E, unique(filesTbl(:, fileKeep), 'rows', 'stable'), 'Keys', {'file_id'}, 'Type', 'left', 'MergeKeys', true);
pairs = { ...
    'area_px2', 'area_um2', 'spatial_res', 'area'; ...
    'duration_frames', 'duration_sec', 'frame_rate_s', 'linear'; ...
    'rise19_frames', 'rise19_sec', 'frame_rate_s', 'linear'; ...
    'fall91_frames', 'fall91_sec', 'frame_rate_s', 'linear'; ...
    'decay_tau_frames', 'decay_tau_sec', 'frame_rate_s', 'linear'; ...
    'dff_auc', 'dff_auc_sec', 'frame_rate_s', 'linear'};
for i = 1:size(pairs,1)
    src = pairs{i,1}; dst = pairs{i,2}; cal = pairs{i,3}; mode = pairs{i,4};
    if ~ismember(src, E.Properties.VariableNames) || ~ismember(cal, E.Properties.VariableNames)
        continue
    end
    x = double(E.(src)); c = double(E.(cal)); y = NaN(size(x));
    ok = isfinite(x) & isfinite(c);
    switch mode
        case 'linear'
            y(ok) = x(ok) .* c(ok);
        case 'area'
            y(ok) = x(ok) .* (c(ok).^2);
    end
    E.(dst) = y;
end
end
