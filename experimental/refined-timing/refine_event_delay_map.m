function out = refine_event_delay_map(source, preprocessed, varargin)
%REFINE_EVENT_DELAY_MAP Prototype BWA-only refined per-pixel timing.
%   OUT = REFINE_EVENT_DELAY_MAP(SOURCE, PREPROCESSED) consumes canonical
%   AQuA timing plus Project A domain masks and returns:
%     - delay_map_refined_sec
%     - delay_confidence
%     - delay_imputed_mask
%     - xSpaTemp_refined
%
% SOURCE may be:
%   - a local-analysis result carrying pathIn
%   - a cluster-output root path
%   - a struct carrying BWA_AQ2 directly
%
% PREPROCESSED may be:
%   - the loaded Project A bundle
%   - a path to a MAT file containing M_canonical/domain_masks/domain_hole_rows

p = inputParser;
p.StructExpand = false;
addRequired(p, 'source');
addRequired(p, 'preprocessed');
addParameter(p, 'eventUID', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'eventID', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'spaRatio', NaN, @(x) isnumeric(x) && isscalar(x) && (isnan(x) || (isfinite(x) && x >= 1)));
addParameter(p, 'frameDtSec', NaN, @(x) isnumeric(x) && isscalar(x) && (isnan(x) || (isfinite(x) && x > 0)));
addParameter(p, 'anchorWindowFrames', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'baselineFrames', 4, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'smoothSpanFrames', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'onsetThresholdFrac', 0.35, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x <= 1);
addParameter(p, 'onsetMinSnr', 1.5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'derivativeWeight', 0.50, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'regularizationSigmaPx', 1.0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'maxWithinBlockShiftFrames', 2.0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
parse(p, source, preprocessed, varargin{:});

ctx = localResolvePrototypeContext(p.Results.source, p.Results.preprocessed, ...
    string(p.Results.eventUID), double(p.Results.eventID));

MCanonical = ctx.bundle.M_canonical;
domainMasks = ctx.bundle.domain_masks;
coarseDelaySec = double(MCanonical.delay_sec);
baseMask = logical(domainMasks.mask_domain_base);
domainMask = logical(domainMasks.mask_domain_filled);
fillAddedMask = logical(domainMasks.mask_fill_added);

frameDtSec = double(p.Results.frameDtSec);
if ~isfinite(frameDtSec) || frameDtSec <= 0
    frameDtSec = localInferFrameDtSec(coarseDelaySec);
end

spaRatio = double(p.Results.spaRatio);
if ~isfinite(spaRatio) || spaRatio < 1
    spaRatio = localInferSpaRatio(coarseDelaySec, baseMask);
end
spaRatio = max(1, round(spaRatio));

t0Frame = round(double(ctx.event.t0));
t1Frame = round(double(ctx.event.t1));
nFrames = size(ctx.movie.stack, 3);
if t0Frame < 1 || t0Frame > nFrames || t1Frame < t0Frame
    error('refine_event_delay_map:BadEventFrames', ...
        'Bad event frames: t0=%d t1=%d nFrames=%d.', t0Frame, t1Frame, nFrames);
end

workFrames = t0Frame:t1Frame;
movieFrameRange = double(ctx.movie.frame_range(:)');
[inRange, workFrameIdx] = ismember(workFrames, movieFrameRange);
if ~all(inRange)
    error('refine_event_delay_map:EventOutsideMovieWindow', ...
        'Event frames %d:%d are not fully contained in movie frame range %d:%d.', ...
        t0Frame, t1Frame, min(movieFrameRange), max(movieFrameRange));
end
movieWin = double(ctx.movie.stack(:,:,workFrameIdx));
coarseFrameMap = t0Frame + coarseDelaySec ./ frameDtSec;
coarseFrameMap(~isfinite(coarseFrameMap)) = NaN;
delayFrameLocal = coarseFrameMap - t0Frame;
[blockLabelMap, blockStats] = localBuildApproxBlocks(baseMask, delayFrameLocal, spaRatio);

pixelDelayFrameLocal = NaN(size(coarseDelaySec));
pixelConfidence = zeros(size(coarseDelaySec));
rawOnsetFrameLocal = NaN(size(coarseDelaySec));
rawPeakSnr = zeros(size(coarseDelaySec));
rawBestScore = NaN(size(coarseDelaySec));
validRawMask = false(size(coarseDelaySec));

opts = struct();
opts.anchorWindowFrames = round(double(p.Results.anchorWindowFrames));
opts.baselineFrames = round(double(p.Results.baselineFrames));
opts.smoothSpanFrames = round(double(p.Results.smoothSpanFrames));
opts.onsetThresholdFrac = double(p.Results.onsetThresholdFrac);
opts.onsetMinSnr = double(p.Results.onsetMinSnr);
opts.derivativeWeight = double(p.Results.derivativeWeight);

for b = 1:numel(blockStats)
    pix = blockStats(b).pixel_idx;
    if isempty(pix)
        continue
    end
    coarseAnchorLocal = blockStats(b).coarse_frame_local;
    if ~isfinite(coarseAnchorLocal)
        continue
    end
    onsetLocal = NaN(numel(pix), 1);
    conf = zeros(numel(pix), 1);
    peakSnr = zeros(numel(pix), 1);
    bestScore = NaN(numel(pix), 1);
    bestFrame = NaN(numel(pix), 1);
    for i = 1:numel(pix)
        [r, c] = ind2sub(size(baseMask), pix(i));
        trace = reshape(movieWin(r, c, :), 1, []);
        est = estimate_pixel_onset_time(trace, coarseAnchorLocal, opts);
        onsetLocal(i) = est.onset_frame_local;
        conf(i) = est.confidence;
        peakSnr(i) = est.peak_snr;
        bestScore(i) = est.best_score;
        bestFrame(i) = est.best_frame_local;
    end
    candidate = isfinite(bestFrame) & isfinite(bestScore) & bestScore > 0;
    if any(candidate)
        rankResidual = localRankResidual(bestScore(candidate), double(p.Results.maxWithinBlockShiftFrames));
        candidateIdx = find(candidate);
        support = candidate & (isfinite(onsetLocal) | bestScore >= max(0.75, 0.5 * double(p.Results.onsetMinSnr)));
        supportIdx = find(support);
        supportRankResidual = rankResidual(ismember(candidateIdx, supportIdx));
        supportOnset = bestFrame(support);
        hasDirect = isfinite(onsetLocal(support));
        supportOnset(hasDirect) = onsetLocal(support);
        supportResidual = supportOnset - coarseAnchorLocal;
        supportResidual = max(-double(p.Results.maxWithinBlockShiftFrames), ...
            min(double(p.Results.maxWithinBlockShiftFrames), supportResidual));
        supportResidual = 0.5 * supportResidual + 0.5 * supportRankResidual;
        refinedLocal = coarseAnchorLocal + supportResidual;
        rawOnsetFrameLocal(pix(support)) = supportOnset;
        pixelDelayFrameLocal(pix(support)) = refinedLocal;
        pixelConfidence(pix(support)) = max(conf(support), min(0.6, bestScore(support) ./ 6));
        rawPeakSnr(pix(support)) = peakSnr(support);
        rawBestScore(pix(support)) = bestScore(support);
        validRawMask(pix(support)) = true;
    end
end

status = struct('success', true, 'code', "ok", 'message', "", ...
    'used_coarse_fallback', false, 'movie_alignment_method', string(ctx.movie.alignment.method), ...
    'movie_alignment_note', string(ctx.movie.alignment.note));

if ~any(validRawMask(:))
    status.success = false;
    status.code = "no_raw_support";
    status.message = "No usable raw per-pixel onset evidence; returned coarse-anchor fallback.";
    status.used_coarse_fallback = true;
    [pixelDelayFrameLocal, imputedFromMissing] = localImputeDelay(delayFrameLocal, ...
        baseMask & isfinite(delayFrameLocal), domainMask);
    delayImputedMask = domainMask & (~baseMask | ~isfinite(delayFrameLocal) | imputedFromMissing);
    pixelConfidence = zeros(size(coarseDelaySec));
    pixelConfidence(baseMask & domainMask) = 0.35;
    pixelConfidence(fillAddedMask) = 0.15;
    pixelConfidence(delayImputedMask & ~fillAddedMask) = min(pixelConfidence(delayImputedMask & ~fillAddedMask), 0.25);
else
    residualMap = pixelDelayFrameLocal - delayFrameLocal;
    residualMap = localMedianResidual(residualMap, validRawMask, domainMask);
    if double(p.Results.regularizationSigmaPx) > 0
        residualMap = localSmoothResidual(residualMap, pixelConfidence, validRawMask, ...
            domainMask, double(p.Results.regularizationSigmaPx));
    end
    residualMap = localSmoothResidual(residualMap, max(pixelConfidence, 0.25), domainMask, ...
        domainMask, max(1.5, 1.5 * double(p.Results.regularizationSigmaPx)));
    pixelDelayFrameLocal(domainMask) = delayFrameLocal(domainMask) + residualMap(domainMask);
    pixelDelayFrameLocal(domainMask) = 0.65 * pixelDelayFrameLocal(domainMask) + 0.35 * delayFrameLocal(domainMask);

    pixelDelayFrameLocal = localPreserveBlockAnchors(pixelDelayFrameLocal, delayFrameLocal, blockStats);
    imputeSeedMask = domainMask & validRawMask & isfinite(pixelDelayFrameLocal);
    [pixelDelayFrameLocal, imputedFromMissing] = localImputeDelay(pixelDelayFrameLocal, imputeSeedMask, domainMask);
    delayImputedMask = domainMask & (~validRawMask | fillAddedMask | imputedFromMissing);
    pixelConfidence = localImputeConfidence(pixelConfidence, validRawMask, domainMask);
    pixelConfidence(delayImputedMask) = min(pixelConfidence(delayImputedMask), 0.55);
    pixelConfidence(fillAddedMask) = min(pixelConfidence(fillAddedMask), 0.45);
end

pixelConfidence(~domainMask) = 0;
delayMapRefinedSec = NaN(size(coarseDelaySec));
delayMapRefinedSec(domainMask) = max(0, pixelDelayFrameLocal(domainMask) .* frameDtSec);
delayMapRefinedSec = localDenoiseDelayMap( ...
    delayMapRefinedSec, coarseDelaySec, pixelConfidence, validRawMask, domainMask, frameDtSec);
xSpaTempRefined = reconstruct_xspatemp_refined( ...
    delayMapRefinedSec, domainMask, t0Frame, t1Frame, ctx.event.spatial_size, ...
    'frameDtSec', frameDtSec);

stats = struct();
stats.event_uid = string(ctx.event.event_uid);
stats.event_id = double(ctx.event.event_id);
stats.file_id = string(ctx.event.file_id);
stats.file_row = double(ctx.event.file_row);
stats.spa_ratio = double(spaRatio);
stats.frame_dt_sec = double(frameDtSec);
stats.coarse_delay_pixels = nnz(baseMask & isfinite(coarseDelaySec));
stats.refined_delay_pixels = nnz(domainMask & isfinite(delayMapRefinedSec));
stats.raw_supported_pixels = nnz(validRawMask);
stats.imputed_pixels = nnz(delayImputedMask);
stats.fill_added_pixels = nnz(fillAddedMask);
stats.block_count = numel(blockStats);
stats.mean_confidence = mean(pixelConfidence(domainMask), 'omitnan');
stats.mean_peak_snr = mean(rawPeakSnr(validRawMask), 'omitnan');
stats.mean_best_score = mean(rawBestScore(validRawMask), 'omitnan');
stats.min_refined_delay_sec = min(delayMapRefinedSec(domainMask), [], 'omitnan');
stats.max_refined_delay_sec = max(delayMapRefinedSec(domainMask), [], 'omitnan');
stats.status_code = status.code;
stats.used_coarse_fallback = double(status.used_coarse_fallback);

out = struct();
out.event_uid = string(ctx.event.event_uid);
out.event_id = double(ctx.event.event_id);
out.source_label = string(ctx.source_label);
out.bundle_contract = localGetFieldOr(ctx.bundle, 'audit.bundle_contract', "");
out.status = status;
out.delay_map_refined_sec = delayMapRefinedSec;
out.delay_confidence = pixelConfidence;
out.delay_imputed_mask = delayImputedMask;
out.xSpaTemp_refined = xSpaTempRefined;
out.stats = stats;
out.raw = struct( ...
    'onset_frame_local', rawOnsetFrameLocal, ...
    'peak_snr', rawPeakSnr, ...
    'best_score', rawBestScore, ...
    'valid_mask', validRawMask);
out.inputs = struct( ...
    'M_canonical', MCanonical, ...
    'domain_masks', domainMasks, ...
    'domain_hole_rows', localGetFieldOr(ctx.bundle, 'domain_hole_rows', []));
out.anchor = struct( ...
    'delay_map_coarse_sec', coarseDelaySec, ...
    'delay_map_coarse_frame_local', delayFrameLocal, ...
    'block_label_map', blockLabelMap, ...
    'block_stats', blockStats);
out.event = ctx.event;
out.movie = struct( ...
    'frames', workFrames, ...
    'source_frame_range', movieFrameRange, ...
    'frame_dt_sec', frameDtSec, ...
    'stack', movieWin, ...
    'alignment', ctx.movie.alignment);
out.canonical = struct( ...
    'xSpaTemp', ctx.event.xSpaTemp, ...
    'delay_sec', coarseDelaySec, ...
    'footprint_mask', logical(MCanonical.footprint_mask));
end

function ctx = localResolvePrototypeContext(source, preprocessed, eventUID, eventID)
bundle = localResolveBundle(preprocessed);
[BWA_AQ2, sourceLabel, sourceResult, loadedStruct] = localResolveSource(source);
event = localResolveEvent(BWA_AQ2, sourceResult, bundle, eventUID, eventID);
movie = localResolveMovie(loadedStruct, BWA_AQ2, sourceResult, event.file_row, event.spatial_size);

ctx = struct();
ctx.bundle = bundle;
ctx.BWA_AQ2 = BWA_AQ2;
ctx.source_label = sourceLabel;
ctx.source_result = sourceResult;
ctx.loaded = loadedStruct;
ctx.event = event;
ctx.movie = movie;
end

function bundle = localResolveBundle(preprocessed)
if ischar(preprocessed) || isstring(preprocessed)
    bundle = load(char(string(preprocessed)));
else
    bundle = preprocessed;
end
if ~isstruct(bundle) || ~isfield(bundle, 'M_canonical') || ~isfield(bundle, 'domain_masks')
    error('refine_event_delay_map:BadPreprocessedBundle', ...
        'Expected preprocessed bundle with M_canonical and domain_masks.');
end
end

function [BWA_AQ2, label, sourceResult, loadedStruct] = localResolveSource(source)
label = "";
sourceResult = [];
loadedStruct = struct();
if isstruct(source) && isfield(source, 'derived') && isstruct(source.derived) && ...
        isfield(source.derived, 'context') && isstruct(source.derived.context) && ...
        isfield(source.derived.context, 'BWA_AQ2')
    cand = source.derived.context.BWA_AQ2;
    if isstruct(cand) && ~isfield(cand, 'opts') && isfield(source.derived.context, 'opts')
        cand.opts = source.derived.context.opts;
    end
    if isstruct(cand) && isfield(cand, 'fts1') && isfield(cand, 'opts')
        BWA_AQ2 = cand;
        sourceResult = source;
        if isfield(source, 'pathIn')
            label = string(source.pathIn);
        end
        return
    end
end
if isstruct(source) && isfield(source, 'pathIn') && strlength(string(source.pathIn)) > 0
    sourceResult = source;
    loadedStruct = loadClusterOutputs(char(string(source.pathIn)), 'attachLegacyMeta', true, 'attachEventMeta', true);
    BWA_AQ2 = loadedStruct.BWA_AQ2;
    label = string(source.pathIn);
    return
end
if ischar(source) || isstring(source)
    loadedStruct = loadClusterOutputs(char(string(source)), 'attachLegacyMeta', true, 'attachEventMeta', true);
    BWA_AQ2 = loadedStruct.BWA_AQ2;
    label = string(source);
    return
end
if isstruct(source) && isfield(source, 'fts1') && isfield(source, 'opts')
    BWA_AQ2 = source;
    sourceResult = source;
    return
end
error('refine_event_delay_map:BadSource', ...
    'Expected source carrying pathIn, a cluster-output root, or BWA_AQ2 directly.');
end

function event = localResolveEvent(BWA_AQ2, sourceResult, bundle, eventUID, eventID)
fileIDs = localFileIDs(BWA_AQ2);
nEvents = numel(fileIDs);
uidList = strings(nEvents, 1);
for i = 1:nEvents
    uidList(i) = sprintf('%s_E%d', fileIDs(i), localEventOrdinalWithinFile(fileIDs, i));
end

if strlength(eventUID) > 0
    idx = find(uidList == string(eventUID), 1, 'first');
elseif isfinite(eventID)
    idx = round(double(eventID));
else
    idx = localInferEventIndexFromBundle(sourceResult, bundle, uidList);
end
if isempty(idx) || idx < 1 || idx > nEvents
    error('refine_event_delay_map:MissingEvent', ...
        'Could not resolve event from eventUID/eventID/bundle.');
end

fileID = string(fileIDs(idx));
uniqueIDs = unique(cellstr(string(fileIDs)), 'stable');
fileRow = find(strcmp(uniqueIDs, char(fileID)), 1, 'first');
sz = double(BWA_AQ2.opts.sz(fileRow, :));
if numel(sz) < 4
    error('refine_event_delay_map:BadImageSize', 'opts.sz row must have at least 4 values.');
end

event = struct();
event.event_id = double(idx);
event.event_uid = uidList(idx);
event.file_id = fileID;
event.file_row = double(fileRow);
event.t0 = double(BWA_AQ2.fts1.loc.t0(idx));
event.t1 = double(BWA_AQ2.fts1.loc.t1(idx));
event.xSpaTemp = BWA_AQ2.fts1.loc.xSpaTemp{idx};
event.spatial_size = [round(sz(1)) round(sz(2)) max(1, round(sz(3))) round(sz(4))];
end

function idx = localInferEventIndexFromBundle(sourceResult, bundle, uidList)
idx = [];
uid = "";
if isfield(bundle, 'audit') && isstruct(bundle.audit) && isfield(bundle.audit, 'event_uid')
    uid = string(bundle.audit.event_uid);
elseif isstruct(sourceResult) && isfield(sourceResult, 'event_uid')
    uid = string(sourceResult.event_uid);
end
if strlength(uid) > 0
    idx = find(uidList == uid, 1, 'first');
end
end

function fileIDs = localFileIDs(BWA_AQ2)
if isfield(BWA_AQ2.fts1, 'basic') && isstruct(BWA_AQ2.fts1.basic) && isfield(BWA_AQ2.fts1.basic, 'fileID')
    fileIDs = string(BWA_AQ2.fts1.basic.fileID);
elseif isfield(BWA_AQ2.fts1, 'fileID')
    fileIDs = string(BWA_AQ2.fts1.fileID);
else
    error('refine_event_delay_map:MissingFileIDs', 'BWA_AQ2 is missing fts1 file IDs.');
end
end

function ord = localEventOrdinalWithinFile(fileIDs, idx)
ord = nnz(fileIDs(1:idx) == fileIDs(idx));
end

function movie = localResolveMovie(loadedStruct, BWA_AQ2, sourceResult, fileRow, spatialSize)
movie = struct('stack', [], 'frame_range', [], 'alignment', struct('method', "none", 'note', ""));
candidateValues = {};

if isstruct(loadedStruct)
    candidateValues = [candidateValues, localCollectMovieCandidates(loadedStruct)];
end
candidateValues = [candidateValues, localCollectMovieCandidates(BWA_AQ2)];
if isstruct(sourceResult)
    candidateValues = [candidateValues, localCollectMovieCandidates(sourceResult)];
end

for i = 1:numel(candidateValues)
    [tmpMovie, frameRange, alignInfo] = localNormalizeMovieCandidate(candidateValues{i}, fileRow, spatialSize, sourceResult);
    if ~isempty(tmpMovie)
        movie.stack = tmpMovie;
        movie.frame_range = frameRange;
        movie.alignment = alignInfo;
        return
    end
end

error('refine_event_delay_map:MissingRawMovie', ...
    ['Could not resolve raw movie stack from source. Expected one of datOrg, dat, dF, df, or data ', ...
     'to be reachable from loadClusterOutputs/BWA_AQ2/source result.']);
end

function out = localCollectMovieCandidates(S)
out = {};
if ~isstruct(S)
    return
end
fieldNames = {'datOrg', 'dat', 'dF', 'df', 'data'};
for i = 1:numel(fieldNames)
    if isfield(S, fieldNames{i})
        out{end+1} = S.(fieldNames{i}); %#ok<AGROW>
    end
end
if isfield(S, 'res') && isstruct(S.res)
    for i = 1:numel(fieldNames)
        if isfield(S.res, fieldNames{i})
            out{end+1} = S.res.(fieldNames{i}); %#ok<AGROW>
        end
    end
end
if isfield(S, 'BWA_AQ2') && isstruct(S.BWA_AQ2)
    out = [out, localCollectMovieCandidates(S.BWA_AQ2)];
end
end

function [movie, frameRange, alignInfo] = localNormalizeMovieCandidate(cand, fileRow, spatialSize, sourceResult)
movie = [];
frameRange = [];
alignInfo = struct('method', "none", 'note', "");
sx = spatialSize(1);
sy = spatialSize(2);
sz3 = spatialSize(3);
st = spatialSize(4);
sourceFrameRange = localSourceMovieFrameRange(sourceResult);

if iscell(cand)
    if fileRow < 1 || fileRow > numel(cand)
        return
    end
    cand = cand{fileRow};
end
if ~isnumeric(cand)
    return
end

sz = size(cand);
if ndims(cand) == 4
    if isequal(sz(1:4), [sx sy sz3 st])
        movie = squeeze(cand(:,:,1,:));
        frameRange = 1:st;
        alignInfo.method = "exact_4d";
        return
    end
    if isequal(sz(1:4), [sx sy st sz3])
        movie = reshape(cand, [sx sy st]);
        frameRange = 1:st;
        alignInfo.method = "exact_4d_permuted";
        return
    end
    if numel(sz) >= 4 && sz(3) == st
        movie = squeeze(cand(:,:,1:1,1:st));
        movie = reshape(movie, size(cand,1), size(cand,2), st);
        frameRange = 1:st;
        alignInfo.method = "collapsed_4d";
        return
    end
end
if ndims(cand) == 3
    if isequal(sz, [sx sy st])
        movie = cand;
        frameRange = 1:st;
        alignInfo.method = "exact_3d";
        return
    end
    if ~isempty(sourceFrameRange) && sz(1) == sx && sz(2) == sy && sz(3) == numel(sourceFrameRange)
        movie = cand;
        frameRange = sourceFrameRange;
        alignInfo.method = "bounded_frame_window_3d";
        alignInfo.note = sprintf('raw movie contains absolute frames %d:%d only', ...
            min(sourceFrameRange), max(sourceFrameRange));
        return
    end
    if isequal(sz, [sx sy sz3]) && st == 1
        movie = cand;
        frameRange = 1:sz3;
        alignInfo.method = "exact_single_frame";
        return
    end
    if sz(1) >= sx && sz(2) >= sy && sz(3) == st
        dRow = sz(1) - sx;
        dCol = sz(2) - sy;
        if dRow <= 20 && dCol <= 20
            r0 = floor(dRow / 2) + 1;
            c0 = floor(dCol / 2) + 1;
            rows = r0:(r0 + sx - 1);
            cols = c0:(c0 + sy - 1);
            movie = cand(rows, cols, :);
            frameRange = 1:st;
            alignInfo.method = "center_crop_3d";
            alignInfo.note = sprintf('center-cropped raw movie from %dx%dx%d to %dx%dx%d', ...
                sz(1), sz(2), sz(3), sx, sy, st);
            return
        end
    end
end
if ismatrix(cand) && numel(cand) == sx * sy * st
    movie = reshape(cand, [sx sy st]);
    frameRange = 1:st;
    alignInfo.method = "reshaped_matrix";
end
end

function frameRange = localSourceMovieFrameRange(sourceResult)
frameRange = [];
if isstruct(sourceResult) && isfield(sourceResult, 'raw_movie_frame_range')
    frameRange = round(double(sourceResult.raw_movie_frame_range(:)'));
    frameRange = frameRange(isfinite(frameRange));
end
end

function dtSec = localInferFrameDtSec(delaySec)
vals = unique(sort(delaySec(isfinite(delaySec))));
if numel(vals) < 2
    dtSec = 1.0;
    return
end
d = diff(vals);
d = d(d > 0);
if isempty(d)
    dtSec = 1.0;
else
    dtSec = median(d);
end
end

function spaRatio = localInferSpaRatio(delaySec, baseMask)
spaRatio = 6;
if ~any(baseMask(:))
    return
end
Q = round(double(delaySec) .* 1e6);
Q(~baseMask | ~isfinite(delaySec)) = intmin('int32');
lengths = [];
for r = 1:size(Q, 1)
    row = Q(r, :);
    lengths = [lengths; localRunLengths(row)]; %#ok<AGROW>
end
for c = 1:size(Q, 2)
    col = Q(:, c)';
    lengths = [lengths; localRunLengths(col)]; %#ok<AGROW>
end
lengths = lengths(lengths >= 2 & lengths <= 16);
if ~isempty(lengths)
    spaRatio = round(median(lengths));
end
end

function residual = localRankResidual(score, maxShift)
score = double(score(:));
n = numel(score);
if n < 2 || maxShift <= 0
    residual = zeros(size(score));
    return
end
[~, ord] = sort(score, 'ascend');
rankPos = zeros(size(score));
rankPos(ord) = 1:n;
rankCentered = (rankPos - mean(1:n)) ./ max(1, (n - 1) / 2);
residual = rankCentered .* maxShift;
end

function out = localRunLengths(v)
out = zeros(0, 1);
if isempty(v)
    return
end
startIdx = 1;
for i = 2:numel(v) + 1
    if i > numel(v) || v(i) ~= v(startIdx)
        if v(startIdx) ~= intmin('int32')
            out(end+1, 1) = i - startIdx; %#ok<AGROW>
        end
        startIdx = i;
    end
end
end

function [labelMap, stats] = localBuildApproxBlocks(baseMask, coarseFrameLocal, spaRatio)
sx = size(baseMask, 1);
sy = size(baseMask, 2);
gridRow = ceil((1:sx)' ./ spaRatio);
gridCol = ceil((1:sy) ./ spaRatio);
gridRows = max(gridRow);
gridCols = max(gridCol);
labelMap = zeros(size(baseMask));
stats = repmat(struct('label', 0, 'pixel_idx', [], 'coarse_frame_local', NaN), 0, 1);
nextLabel = 0;
for gr = 1:gridRows
    rr = (gridRow == gr);
    for gc = 1:gridCols
        cc = (gridCol == gc);
        tileMask = false(size(baseMask));
        tileMask(rr, cc) = true;
        pix = find(tileMask & baseMask & isfinite(coarseFrameLocal));
        if isempty(pix)
            continue
        end
        nextLabel = nextLabel + 1;
        labelMap(pix) = nextLabel;
        stats(nextLabel).label = nextLabel;
        stats(nextLabel).pixel_idx = pix;
        stats(nextLabel).coarse_frame_local = median(coarseFrameLocal(pix), 'omitnan');
    end
end
end

function residualSmooth = localSmoothResidual(residualMap, confidence, validMask, domainMask, sigmaPx)
residualSmooth = residualMap;
w = zeros(size(residualMap));
w(validMask) = max(confidence(validMask), 0.05);
num = imgaussfilt(localNanToZero(residualMap) .* w, sigmaPx, 'FilterDomain', 'spatial');
den = imgaussfilt(w, sigmaPx, 'FilterDomain', 'spatial');
keep = domainMask & den > 0;
residualSmooth(keep) = num(keep) ./ den(keep);
end

function residualOut = localMedianResidual(residualIn, validMask, domainMask)
residualOut = residualIn;
work = residualIn;
work(~validMask) = NaN;
pad = nan(size(work, 1) + 2, size(work, 2) + 2);
pad(2:end-1, 2:end-1) = work;
for dr = -1:1
    for dc = -1:1
        nbr = pad((2:end-1) + dr, (2:end-1) + dc);
        if dr == -1 && dc == -1
            stack = nbr;
        else
            stack(:,:,end+1) = nbr; %#ok<AGROW>
        end
    end
end
med = median(stack, 3, 'omitnan');
replace = domainMask & ~validMask & isfinite(med);
residualOut(replace) = med(replace);
smoothable = domainMask & validMask & isfinite(med);
residualOut(smoothable) = 0.5 * residualOut(smoothable) + 0.5 * med(smoothable);
end

function delayFrameLocal = localPreserveBlockAnchors(delayFrameLocal, coarseFrameLocal, blockStats)
for i = 1:numel(blockStats)
    pix = blockStats(i).pixel_idx;
    if isempty(pix)
        continue
    end
    valid = pix(isfinite(delayFrameLocal(pix)) & isfinite(coarseFrameLocal(pix)));
    if isempty(valid)
        continue
    end
    shift = median(delayFrameLocal(valid), 'omitnan') - median(coarseFrameLocal(valid), 'omitnan');
    delayFrameLocal(valid) = delayFrameLocal(valid) - shift;
end
end

function [delayOut, imputedMask] = localImputeDelay(delayIn, seedMask, domainMask)
delayOut = delayIn;
imputedMask = false(size(delayIn));
need = domainMask & ~isfinite(delayOut);
if ~any(need(:)) || ~any(seedMask(:))
    return
end
[~, nnIdx] = bwdist(seedMask);
delayOut(need) = delayOut(nnIdx(need));
imputedMask(need) = true;
end

function confOut = localImputeConfidence(confIn, validMask, domainMask)
confOut = confIn;
need = domainMask & ~validMask;
if ~any(need(:)) || ~any(validMask(:))
    confOut(need) = 0.25;
    return
end
[~, nnIdx] = bwdist(validMask);
confOut(need) = confOut(nnIdx(need));
confOut(need) = min(confOut(need), 0.60);
end

function delayOut = localDenoiseDelayMap(delayIn, coarseDelaySec, confidence, validRawMask, domainMask, frameDtSec)
delayOut = delayIn;
if ~any(domainMask(:))
    return
end

coarseWork = double(coarseDelaySec);
coarseWork(~domainMask) = NaN;
refinedWork = double(delayIn);
refinedWork(~domainMask) = NaN;

refinedMed = medfilt2(localNanToZero(refinedWork), [5 5], 'symmetric');
coarseMed = medfilt2(localNanToZero(coarseWork), [5 5], 'symmetric');
refinedSmooth = localGaussianWithinMask(refinedMed, domainMask, 1.4);
coarseSmooth = localGaussianWithinMask(coarseMed, domainMask, 1.0);

supportWeight = zeros(size(confidence));
supportWeight(validRawMask) = max(confidence(validRawMask), 0.15);
supportWeight = localGaussianWithinMask(supportWeight, domainMask, 1.2);
supportWeight(domainMask) = min(0.55, supportWeight(domainMask));

localResidual = refinedSmooth - coarseSmooth;
maxResidualSec = max(2 * frameDtSec, 4.0);
localResidual(domainMask) = max(-maxResidualSec, min(maxResidualSec, localResidual(domainMask)));

blend = zeros(size(delayIn));
blend(domainMask) = 0.12 + 0.88 * supportWeight(domainMask);
delayOut(domainMask) = coarseSmooth(domainMask) + blend(domainMask) .* localResidual(domainMask);
delayOut(domainMask) = max(0, delayOut(domainMask));
end

function out = localGaussianWithinMask(in, mask, sigmaPx)
w = double(mask & isfinite(in));
num = imgaussfilt(localNanToZero(in) .* w, sigmaPx, 'FilterDomain', 'spatial');
den = imgaussfilt(w, sigmaPx, 'FilterDomain', 'spatial');
out = NaN(size(in));
keep = mask & den > 0;
out(keep) = num(keep) ./ den(keep);
end

function x = localNanToZero(x)
x(~isfinite(x)) = 0;
end

function value = localGetFieldOr(S, pathExpr, defaultValue)
value = defaultValue;
if ~isstruct(S)
    return
end
parts = strsplit(char(pathExpr), '.');
curr = S;
for i = 1:numel(parts)
    if ~isstruct(curr) || ~isfield(curr, parts{i})
        return
    end
    curr = curr.(parts{i});
end
value = curr;
end
