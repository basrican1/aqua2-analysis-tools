function varargout = pair_corridor_fit_memory(action, varargin)
%PAIR_CORRIDOR_FIT_MEMORY Persistent warm-start memory for pair corridor fitting.

switch lower(string(action))
    case "warm_start"
        [paramsOut, meta] = localWarmStart(varargin{:});
        varargout = {paramsOut, meta};
    case "record"
        localRecord(varargin{:});
        varargout = {};
    case "path"
        varargout = {localStorePath()};
    otherwise
        error('pair_corridor_fit_memory:UnknownAction', 'Unknown action: %s', string(action));
end
end

function [paramsOut, meta] = localWarmStart(featureStruct, paramsDefault)
paramsOut = paramsDefault;
meta = struct('used', false, 'nMatches', 0, 'distance', NaN, 'usedAccepted', false);
S = localLoadStore();
if isempty(S.entries)
    return
end
entries = S.entries;
accepted = arrayfun(@(e) isfield(e,'acceptable') && e.acceptable, entries);
if any(accepted)
    entries = entries(accepted);
    meta.usedAccepted = true;
else
    finiteScore = arrayfun(@(e) isfield(e,'score') && isfinite(e.score), entries);
    entries = entries(finiteScore);
    if isempty(entries)
        return
    end
    scores = arrayfun(@(e) double(e.score), entries);
    cutoff = min(prctile(scores, 35), 650);
    keepNear = scores <= cutoff;
    entries = entries(keepNear);
    if isempty(entries)
        return
    end
end
q = localFeatureVector(featureStruct);
M = nan(numel(entries), numel(q));
for i = 1:numel(entries)
    M(i,:) = localFeatureVector(entries(i).features);
end
scale = median(abs(M - median(M,1,'omitnan')), 1, 'omitnan');
scale(~isfinite(scale) | scale <= 0) = 1;
dist = sqrt(sum(((M - q) ./ scale).^2, 2, 'omitnan'));
scoreAdj = zeros(size(dist));
for i = 1:numel(entries)
    if isfield(entries(i),'acceptable') && entries(i).acceptable
        scoreAdj(i) = 0;
    elseif isfield(entries(i),'score') && isfinite(entries(i).score)
        scoreAdj(i) = min(2.0, 0.003 * double(entries(i).score));
    else
        scoreAdj(i) = 2.0;
    end
end
distEff = dist + scoreAdj;
[distEff, ord] = sort(distEff, 'ascend');
entries = entries(ord);
dist = dist(ord);
distEff = distEff(ord);
keep = isfinite(distEff) & distEff <= 5.5;
entries = entries(keep);
dist = dist(keep);
distEff = distEff(keep);
if isempty(entries)
    return
end
k = min(5, numel(entries));
entries = entries(1:k);
dist = dist(1:k);
distEff = distEff(1:k);
w = 1 ./ max(distEff, 0.35);
w = w ./ sum(w);
fields = fieldnames(paramsDefault);
for i = 1:numel(fields)
    f = fields{i};
    v0 = paramsDefault.(f);
    if islogical(v0)
        vals = arrayfun(@(e) double(e.params.(f)), entries);
        paramsOut.(f) = mean(vals .* w') >= 0.5;
    elseif isnumeric(v0) && isscalar(v0)
        vals = arrayfun(@(e) double(e.params.(f)), entries);
        paramsOut.(f) = sum(vals .* w');
    else
        paramsOut.(f) = v0;
    end
end
meta.nMatches = k;
meta.used = true;
end

function localRecord(featureStruct, params, outcome, score)
S = localLoadStore();
entry = struct();
entry.features = featureStruct;
entry.params = params;
entry.acceptable = isfield(outcome, 'acceptable') && logical(outcome.acceptable);
entry.reason = "";
if isfield(outcome, 'reason')
    entry.reason = string(outcome.reason);
end
entry.score = double(score);
entry.timestamp = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
S.entries = [S.entries; entry]; %#ok<AGROW>
maxKeep = 1000;
if numel(S.entries) > maxKeep
    S.entries = S.entries(max(1,end-maxKeep+1):end);
end
storePath = localStorePath();
entries = S.entries; %#ok<NASGU>
tmpPath = [storePath '.tmp'];
save(tmpPath, 'entries');
movefile(tmpPath, storePath, 'f');
end

function S = localLoadStore()
S = struct('entries', repmat(localEmptyEntry(), 0, 1));
storePath = localStorePath();
if ~isfile(storePath)
    return
end
try
    loaded = load(storePath, 'entries');
catch ME
    warning('pair_corridor_fit_memory:LoadFailed', 'Ignoring unreadable fit-memory store: %s', ME.message);
    return
end
if isfield(loaded, 'entries')
    S.entries = loaded.entries;
end
end

function p = localStorePath()
p = fullfile(fileparts(mfilename('fullpath')), '.pair_corridor_fit_memory.mat');
end

function e = localEmptyEntry()
e = struct('features', struct(), 'params', struct(), 'acceptable', false, 'reason', "", 'score', inf, 'timestamp', '');
end

function v = localFeatureVector(s)
fields = {'srcArea','tgtArea','areaRatio','centroidDist','contactFrac','delayDiff','targetWidth','srcDelay','tgtDelay'};
v = nan(1, numel(fields));
for i = 1:numel(fields)
    f = fields{i};
    if isfield(s, f)
        v(i) = double(s.(f));
    end
end
end
