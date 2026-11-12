function out = estimate_pixel_onset_time(trace, coarseFrameLocal, varargin)
%ESTIMATE_PIXEL_ONSET_TIME Conservative one-pixel onset estimate from a raw trace.
%   OUT = ESTIMATE_PIXEL_ONSET_TIME(TRACE, COARSEFRAMELOCAL, OPTS)
%   returns onset_frame_local, confidence, and peak_snr. The estimate is
%   intentionally anchored to a small window around the coarse AQuA timing.
%   The timing signal is a local baseline-normalized dF/F-like trace plus
%   a positive derivative term measured in the same normalized units.

if nargin >= 3 && isstruct(varargin{1})
    opts = varargin{1};
else
    p = inputParser;
    p.StructExpand = false;
    addParameter(p, 'anchorWindowFrames', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
    addParameter(p, 'baselineFrames', 4, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
    addParameter(p, 'smoothSpanFrames', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
    addParameter(p, 'onsetThresholdFrac', 0.35, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x <= 1);
    addParameter(p, 'onsetMinSnr', 1.5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
    addParameter(p, 'derivativeWeight', 0.50, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
    parse(p, varargin{:});
    opts = p.Results;
end

trace = double(trace(:)');
n = numel(trace);
coarseFrameLocal = round(double(coarseFrameLocal));

out = struct('onset_frame_local', NaN, 'confidence', 0, 'peak_snr', 0, ...
    'peak_value', NaN, 'baseline', NaN, 'noise', NaN, ...
    'best_frame_local', NaN, 'best_score', NaN, 'z_peak', NaN, ...
    'signal_type', "dff_like_local");
if isempty(trace) || all(~isfinite(trace))
    return
end

trace = fillmissing(trace, 'nearest');
traceSm = smoothdata(trace, 'gaussian', max(1, round(double(opts.smoothSpanFrames))));
windowStart = max(1, coarseFrameLocal - round(double(opts.anchorWindowFrames)));
windowEnd = min(n, coarseFrameLocal + round(double(opts.anchorWindowFrames)));
if windowStart > windowEnd
    return
end

baselineEnd = max(1, min(windowStart, coarseFrameLocal) - 1);
baselineStart = max(1, baselineEnd - round(double(opts.baselineFrames)) + 1);
if baselineEnd < baselineStart
    baselineIdx = 1:max(1, min(n, coarseFrameLocal));
else
    baselineIdx = baselineStart:baselineEnd;
end
baseline = median(traceSm(baselineIdx), 'omitnan');
noise = 1.4826 * mad(traceSm(baselineIdx), 1);
noise = max(noise, eps);
baselineScale = max(abs(baseline), 2 * noise);

amp = traceSm - baseline;
dAmp = [0 diff(amp)];
dffLike = amp ./ baselineScale;
dffDeriv = max(0, dAmp) ./ baselineScale;
zAmp = amp ./ noise;
score = dffLike + double(opts.derivativeWeight) * dffDeriv;
roi = windowStart:windowEnd;
peakValue = max(dffLike(roi), [], 'omitnan');
peakSnr = max(zAmp(roi), [], 'omitnan');
peakSnr = double(peakSnr(1));
[bestScore, relIdx] = max(score(roi), [], 'omitnan');
bestFrame = roi(relIdx);

out.peak_value = peakValue;
out.baseline = baseline;
out.noise = noise;
out.peak_snr = peakSnr;
out.best_frame_local = double(bestFrame);
out.best_score = double(bestScore);
out.z_peak = max(zAmp(roi), [], 'omitnan');

if ~isfinite(bestScore) || bestScore <= 0
    return
end

threshold = max(double(opts.onsetThresholdFrac) * peakValue, 0.02);
candidate = roi(dffLike(roi) >= threshold & score(roi) > 0);
if isempty(candidate)
    onset = bestFrame;
else
    onset = candidate(1);
end

distPenalty = abs(onset - coarseFrameLocal) / max(1, double(opts.anchorWindowFrames));
confidence = max(peakSnr, 0) / (max(peakSnr, 0) + 2);
confidence = max(confidence, min(0.55, max(bestScore, 0) / 6));
confidence = confidence * exp(-0.5 * distPenalty);

out.onset_frame_local = double(onset);
out.confidence = max(0, min(1, confidence));
end
