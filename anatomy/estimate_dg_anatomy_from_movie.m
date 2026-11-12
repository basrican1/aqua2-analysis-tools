function anatomy = estimate_dg_anatomy_from_movie(moviePath, eventMaps, metadata, opts)
% Estimate DG layer geometry from unmasked two-photon fluorescence.

if nargin < 2
    eventMaps = [];
end
if nargin < 3 || isempty(metadata)
    metadata = struct('label', 'dataset');
end
if nargin < 4 || isempty(opts)
    opts = struct();
end
opts = apply_defaults(opts, struct( ...
    'upperPercentile', 85, ...
    'displayLowPercentile', 0.15, ...
    'displayHighPercentile', 99.85, ...
    'displayCompression', 7, ...
    'minProcessLength', 28, ...
    'maxProcessCandidates', 28, ...
    'maxProposals', 4));

movie = read_movie(moviePath);
[ny, nx, nt] = size(movie);
opts.imageScale = min(max(min(ny, nx) / 512, 0.65), 1.75);
opts.minProcessLength = max(12, round(opts.minProcessLength * opts.imageScale));

projections = struct();
projections.median = median(movie, 3, 'omitnan');
projections.upper = prctile(movie, opts.upperPercentile, 3);
projections.maximum = max(movie, [], 3, 'omitnan');
projections.variability = std(movie, 0, 3, 'omitnan');
projections.mean = mean(movie, 3, 'omitnan');

projections.medianDisplay = compress_projection(projections.median, opts);
projections.upperDisplay = compress_projection(projections.upper, opts);
projections.maximumDisplay = compress_projection(projections.maximum, opts);
projections.variabilityDisplay = compress_projection(projections.variability, opts);

medianNorm = robust_normalize(projections.median, 0.5, 99.7);
upperNorm = robust_normalize(projections.upper, 0.5, 99.7);
maximumNorm = robust_normalize(projections.maximum, 0.5, 99.8);
variabilityNorm = robust_normalize(projections.variability, 1, 99.7);
structural = 0.34 .* medianNorm + 0.36 .* upperNorm + ...
    0.18 .* maximumNorm + 0.12 .* variabilityNorm;
structural = imgaussfilt(structural, 0.8, 'Padding', 'replicate');
structural = robust_normalize(structural, 0, 99.8);
gclBandEvidence = detect_gcl_band(structural, opts.imageScale);

[ridgeStrength, orientation, orientationCoherence] = process_orientation(structural);
processCandidates = detect_process_candidates(structural, ridgeStrength, ...
    orientation, orientationCoherence, opts);
[sgzBandEvidence, sgzCandidatePoints, blades] = fit_sgz_from_soma_chain( ...
    structural, ridgeStrength, processCandidates);
if isempty(blades)
    [sgzBandEvidence, blades] = fit_sgz_from_band(structural, ridgeStrength, ...
        processCandidates, [ny nx]);
end
if isempty(blades)
    blades = fit_sgz_blades(processCandidates, [ny nx]);
end
fluorescenceOnlyBlades = blades;
roiEvidence = load_roi_evidence(metadata, structural);
opts.morphologyScale = estimate_morphology_scale(roiEvidence, opts.imageScale);
roiCandidates = roi_process_candidates(roiEvidence, structural, opts.morphologyScale);
sgzRefinement = [];
if isfield(metadata, 'sgzRefinement')
    sgzRefinement = metadata.sgzRefinement;
end
roiBlades = fit_sgz_blades(roiCandidates, [ny nx], opts.morphologyScale, ...
    sgzRefinement);
if ~isempty(roiBlades)
    blades = reconcile_roi_and_fluorescence_blades(roiBlades, fluorescenceOnlyBlades);
    estimationMode = 'ROI-assisted, fluorescence-checked';
    if ~isempty(sgzRefinement)
        estimationMode = [estimationMode ', user-refined'];
    end
else
    estimationMode = 'fluorescence-only';
end
blades = prune_blade_hypotheses(blades, [ny nx], opts.morphologyScale);
blades = apply_gcl_support(blades, gclBandEvidence, opts.morphologyScale);
[apexEstimate, apexVisible, apexConfidence] = shared_apex_summary(blades, [ny nx]);
[layerLabels, signedDistance, anatomyConfidence, nearestRadial] = ...
    infer_layers([ny nx], blades, orientation, orientationCoherence);

activityScore = build_activity_score(eventMaps, variabilityNorm, structural);
[wholeCellProposals, distalProposals] = build_proposals(activityScore, ...
    structural, ridgeStrength, blades, nearestRadial, opts.maxProposals);

anatomy = struct();
anatomy.metadata = metadata;
anatomy.metadata.movieSize = [ny nx nt];
anatomy.metadata.structuralSource = ...
    '0.34 median + 0.36 p85 + 0.18 max + 0.12 temporal SD';
anatomy.projections = projections;
anatomy.structural = structural;
anatomy.gclBandEvidence = gclBandEvidence;
anatomy.ridgeStrength = ridgeStrength;
anatomy.orientation = orientation;
anatomy.orientationCoherence = orientationCoherence;
anatomy.sgzBandEvidence = sgzBandEvidence;
anatomy.sgzCandidatePoints = sgzCandidatePoints;
anatomy.processCandidates = processCandidates;
anatomy.blades = blades;
anatomy.apexEstimate = apexEstimate;
anatomy.apexVisible = apexVisible;
anatomy.apexConfidence = apexConfidence;
anatomy.fluorescenceOnlyBlades = fluorescenceOnlyBlades;
anatomy.roiEvidence = roiEvidence;
anatomy.roiProcessCandidates = roiCandidates;
anatomy.sgzRefinement = sgzRefinement;
anatomy.estimationMode = estimationMode;
anatomy.layerLabels = layerLabels;
anatomy.signedDistance = signedDistance;
anatomy.anatomyConfidence = anatomyConfidence;
anatomy.nearestRadial = nearestRadial;
anatomy.activityScore = activityScore;
anatomy.wholeCellProposals = wholeCellProposals;
anatomy.distalProposals = distalProposals;
anatomy.opts = opts;
end

function movie = read_movie(path)
movie = tiffreadVolume(path);
movie = squeeze(single(movie));
if ndims(movie) ~= 3
    error('Expected a Y x X x T movie at %s; got %s.', path, mat2str(size(movie)));
end
movie(~isfinite(movie)) = 0;
end

function out = compress_projection(img, opts)
img = single(img);
vals = img(isfinite(img));
if isempty(vals)
    out = zeros(size(img), 'single');
    return
end
lo = prctile(vals, opts.displayLowPercentile);
hi = prctile(vals, opts.displayHighPercentile);
if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
    lo = min(vals);
    hi = max(vals);
end
x = min(max((img - lo) ./ max(hi - lo, eps('single')), 0), 1);
gain = opts.displayCompression;
out = asinh(gain .* x) ./ asinh(gain);
end

function out = robust_normalize(img, loPct, hiPct)
img = single(img);
vals = img(isfinite(img));
if isempty(vals)
    out = zeros(size(img), 'single');
    return
end
lo = prctile(vals, loPct);
hi = prctile(vals, hiPct);
if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
    lo = min(vals);
    hi = max(vals);
end
out = min(max((img - lo) ./ max(hi - lo, eps('single')), 0), 1);
out(~isfinite(out)) = 0;
end

function evidence = detect_gcl_band(structural, scale)
% Favor broad dark bands with locally granular texture; keep this cue weak.
fineSigma = max(2.5, 4 * scale);
broadSigma = max(12, 24 * scale);
fine = imgaussfilt(structural, fineSigma, 'Padding', 'replicate');
broad = imgaussfilt(structural, broadSigma, 'Padding', 'replicate');
darkValley = robust_normalize(max(broad - fine, 0), 0, 99.5);

window = max(5, 2 * floor(4 * scale) + 1);
texture = stdfilt(structural, true(window));
texture = imgaussfilt(texture, max(3, 7 * scale), 'Padding', 'replicate');
texture = robust_normalize(texture, 2, 99.5);
evidence = robust_normalize(0.72 .* darkValley + 0.28 .* texture, 0, 99.5);
end

function [ridge, theta, coherence] = process_orientation(img)
ridge = fibermetric(img, [2 4 6 9], 'ObjectPolarity', 'bright');
ridge = robust_normalize(ridge, 0, 99.7);

sm = imgaussfilt(img, 1.4, 'Padding', 'replicate');
[gx, gy] = imgradientxy(sm, 'sobel');
jxx = imgaussfilt(gx .^ 2, 4, 'Padding', 'replicate');
jyy = imgaussfilt(gy .^ 2, 4, 'Padding', 'replicate');
jxy = imgaussfilt(gx .* gy, 4, 'Padding', 'replicate');
delta = sqrt((jxx - jyy) .^ 2 + 4 .* jxy .^ 2);
coherence = delta ./ (jxx + jyy + eps('single'));
coherence = min(max(single(coherence), 0), 1);

gradientAxis = 0.5 .* atan2(2 .* jxy, jxx - jyy);
theta = single(gradientAxis + pi / 2);
theta(ridge < 0.03) = nan;
end

function candidates = detect_process_candidates(structural, ridge, theta, coherence, opts)
supportVals = ridge(ridge > 0 & isfinite(ridge));
if isempty(supportVals)
    candidates = empty_process_candidates();
    return
end
ridgeFloor = max(0.08, prctile(supportVals, 62));
mask = ridge >= ridgeFloor & structural >= prctile(structural(:), 48);
mask = imclose(mask, strel('disk', 2));
mask = bwareaopen(mask, 18, 8);
skel = bwskel(mask, 'MinBranchLength', 7);

[H, T, R] = hough(skel, 'Theta', -90:0.5:89.5);
nPeaks = min(60, max(10, round(nnz(skel) / 45)));
peaks = houghpeaks(H, nPeaks, 'Threshold', ceil(0.17 * max(H(:))));
lines = houghlines(skel, T, R, peaks, 'FillGap', 16, ...
    'MinLength', opts.minProcessLength);

endpoints = bwmorph(skel, 'endpoints');
branchpoints = bwmorph(skel, 'branchpoints');
candidates = empty_process_candidates();
for i = 1:numel(lines)
    p1 = double(lines(i).point1);
    p2 = double(lines(i).point2);
    len = norm(p2 - p1);
    if len < opts.minProcessLength
        continue
    end
    borderMargin = 7;
    onBorder1 = p1(1) <= borderMargin || p1(1) >= size(structural, 2) - borderMargin || ...
        p1(2) <= borderMargin || p1(2) >= size(structural, 1) - borderMargin;
    onBorder2 = p2(1) <= borderMargin || p2(1) >= size(structural, 2) - borderMargin || ...
        p2(2) <= borderMargin || p2(2) >= size(structural, 1) - borderMargin;
    if onBorder1 && onBorder2
        continue
    end
    q1 = endpoint_complexity(p1, structural, endpoints, branchpoints);
    q2 = endpoint_complexity(p2, structural, endpoints, branchpoints);
    if q1 <= q2
        base = p1;
        distal = p2;
    else
        base = p2;
        distal = p1;
    end
    radial = (distal - base) ./ max(len, eps);
    lineSamples = sample_line(base, distal, size(structural));
    ridgeEvidence = mean(ridge(lineSamples), 'omitnan');
    coherenceEvidence = mean(coherence(lineSamples), 'omitnan');
    branchAsymmetry = abs(q2 - q1) ./ max(q1 + q2, eps);
    confidence = min(1, 0.32 .* min(len / 80, 1) + ...
        0.34 .* ridgeEvidence + 0.20 .* coherenceEvidence + ...
        0.14 .* branchAsymmetry);
    mid = 0.5 .* (base + distal);
    localTheta = theta(round_clip(mid(2), size(theta, 1)), ...
        round_clip(mid(1), size(theta, 2)));
    if isfinite(localTheta)
        localAxis = [cos(localTheta), sin(localTheta)];
        alignment = abs(dot(radial, localAxis));
        confidence = confidence .* (0.65 + 0.35 .* alignment);
    end
    candidate = struct('base', base, 'distal', distal, 'radial', radial, ...
        'length', len, 'confidence', confidence, 'ridgeEvidence', ridgeEvidence, ...
        'branchAsymmetry', branchAsymmetry);
    if ~duplicates_candidate(candidate, candidates)
        candidates(end + 1) = candidate; %#ok<AGROW>
    end
end

if numel(candidates) > opts.maxProcessCandidates
    [~, order] = sort([candidates.confidence] .* sqrt([candidates.length]), 'descend');
    candidates = candidates(order(1:opts.maxProcessCandidates));
end
end

function [bandEvidence, blades] = fit_sgz_from_band(structural, ridge, candidates, imgSize)
% Detect the dense, coarse SGZ soma band before interpreting radial processes.
density = imgaussfilt(structural, 5.5, 'Padding', 'replicate');
coarseFiber = fibermetric(density, [10 16 24 32], 'ObjectPolarity', 'bright');
density = robust_normalize(density, 0, 99.5);
coarseFiber = robust_normalize(coarseFiber, 0, 99.5);
bandEvidence = robust_normalize(0.62 .* coarseFiber + 0.38 .* density, 0, 99.7);

positive = bandEvidence(bandEvidence > 0 & isfinite(bandEvidence));
blades = empty_blades();
if isempty(positive)
    return
end
floorVal = max(0.16, prctile(positive, 73));
mask = bandEvidence >= floorVal & density >= prctile(density(:), 52);
mask = imclose(mask, strel('disk', 7));
mask = imopen(mask, strel('disk', 2));
mask = bwareaopen(mask, 90, 8);

cc = bwconncomp(mask, 8);
if cc.NumObjects == 0
    return
end
componentScore = zeros(cc.NumObjects, 1);
for i = 1:cc.NumObjects
    pix = cc.PixelIdxList{i};
    componentScore(i) = sum(bandEvidence(pix), 'omitnan') .* log1p(numel(pix));
end
[componentScore, order] = sort(componentScore, 'descend');
if isempty(order)
    return
end
keepCount = 1;
if numel(order) >= 2 && componentScore(2) >= 0.48 * componentScore(1)
    keepCount = 2;
end

for i = 1:keepCount
    pix = cc.PixelIdxList{order(i)};
    blade = fit_band_component(pix, bandEvidence, structural, ridge, candidates, imgSize);
    if ~isempty(blade)
        blades(end + 1) = blade; %#ok<AGROW>
    end
end
end

function blade = fit_band_component(pix, bandEvidence, structural, ridge, candidates, imgSize)
blade = [];
[yy, xx] = ind2sub(imgSize, pix);
weights = double(bandEvidence(pix));
weights = weights .^ 1.6;
weights = weights ./ max(sum(weights), eps);
center = [sum(double(xx) .* weights), sum(double(yy) .* weights)];
coords = [double(xx) - center(1), double(yy) - center(2)];
covar = coords' * (coords .* weights);
[vec, val] = eig(covar);
[~, majorIdx] = max(diag(val));
tangent = vec(:, majorIdx)';
if tangent(1) < 0
    tangent = -tangent;
end
normal = [-tangent(2), tangent(1)];
t = coords * tangent';
r = coords * normal';

tLo = weighted_percentile(t, weights, 3);
tHi = weighted_percentile(t, weights, 97);
span = tHi - tLo;
if span < 65
    return
end
keep = t >= tLo & t <= tHi;
tFit = t(keep);
rFit = r(keep);
wFit = weights(keep);
degree = min(2, max(1, floor(numel(tFit) / 80)));
X = ones(numel(tFit), degree + 1);
for p = 1:degree
    X(:, p + 1) = tFit .^ p;
end
beta = robust_weighted_fit(X, rFit, wFit);
residual = rFit - X * beta;
residualPx = sqrt(sum(wFit .* residual .^ 2) ./ max(sum(wFit), eps));

tSample = linspace(tLo - 28, tHi + 28, 240)';
R = ones(numel(tSample), degree + 1);
for p = 1:degree
    R(:, p + 1) = tSample .^ p;
end
rSample = R * beta;
curveXY = center + tSample .* tangent + rSample .* normal;
if degree == 1
    drdt = repmat(beta(2), size(tSample));
else
    drdt = beta(2) + 2 .* beta(3) .* tSample;
end
curveTangent = tangent + drdt .* normal;
curveTangent = curveTangent ./ max(vecnorm(curveTangent, 2, 2), eps);
curveNormal = [-curveTangent(:, 2), curveTangent(:, 1)];

[positiveSupport, negativeSupport] = side_support(curveXY, curveNormal, structural, ridge);
if negativeSupport > positiveSupport
    curveNormal = -curveNormal;
end
sideAsymmetry = abs(positiveSupport - negativeSupport) ./ ...
    max(positiveSupport + negativeSupport, eps);

processAgreement = 0;
if ~isempty(candidates)
    agreements = zeros(numel(candidates), 1);
    candidateWeights = zeros(numel(candidates), 1);
    for c = 1:numel(candidates)
        [d2, k] = min(sum((curveXY - candidates(c).base) .^ 2, 2));
        if d2 <= 90 ^ 2
            agreements(c) = abs(dot(candidates(c).radial, curveNormal(k, :)));
            candidateWeights(c) = candidates(c).confidence;
        end
    end
    if sum(candidateWeights) > 0
        processAgreement = sum(agreements .* candidateWeights) ./ sum(candidateWeights);
    end
end

spanTerm = min(1, span / 240);
residualTerm = exp(-residualPx / 34);
confidence = (0.42 .* spanTerm + 0.28 .* residualTerm + ...
    0.18 .* sideAsymmetry + 0.12 .* processAgreement);
avgRadial = mean_unit(curveNormal);
blade = struct('curve', [curveXY curveNormal], 'radial', avgRadial, ...
    'tangent', tangent, 'confidence', min(max(confidence, 0), 1), ...
    'residualPx', residualPx, 'candidateIndices', []);
end

function [positiveSupport, negativeSupport] = side_support(curveXY, curveNormal, structural, ridge)
evidence = 0.58 .* structural + 0.42 .* ridge;
positiveSupport = 0;
negativeSupport = 0;
offsets = 22:8:118;
indices = round(linspace(8, size(curveXY, 1) - 8, 70));
indices = indices(indices >= 1 & indices <= size(curveXY, 1));
for k = indices
    for d = offsets
        p = curveXY(k, :) + d .* curveNormal(k, :);
        n = curveXY(k, :) - d .* curveNormal(k, :);
        positiveSupport = positiveSupport + sample_image(evidence, p);
        negativeSupport = negativeSupport + sample_image(evidence, n);
    end
end
end

function value = sample_image(img, point)
x = round(point(1));
y = round(point(2));
if x < 1 || x > size(img, 2) || y < 1 || y > size(img, 1)
    value = 0;
else
    value = double(img(y, x));
end
end

function q = weighted_percentile(vals, weights, percentile)
vals = vals(:);
weights = weights(:);
keep = isfinite(vals) & isfinite(weights) & weights > 0;
vals = vals(keep);
weights = weights(keep);
if isempty(vals)
    q = nan;
    return
end
[vals, order] = sort(vals);
weights = weights(order);
cdf = cumsum(weights) ./ sum(weights);
q = vals(find(cdf >= percentile / 100, 1, 'first'));
end

function candidates = empty_process_candidates()
candidates = struct('base', {}, 'distal', {}, 'radial', {}, 'length', {}, ...
    'confidence', {}, 'ridgeEvidence', {}, 'branchAsymmetry', {});
end

function [somaEvidence, points, blades] = fit_sgz_from_soma_chain( ...
        structural, ridge, processCandidates)
% Find a long chain of blob-like putative somata without using manual ROIs.
blob = imgaussfilt(structural, 2.2, 'Padding', 'replicate');
background = imgaussfilt(structural, 10, 'Padding', 'replicate');
somaEvidence = robust_normalize(max(blob - background, 0), 0, 99.7);
regional = imextendedmax(blob, 0.018);
regional = regional & somaEvidence >= prctile(somaEvidence(:), 62) & ...
    structural >= prctile(structural(:), 52);
regional = bwareaopen(regional, 2, 8);
stats = regionprops(regional, somaEvidence, 'WeightedCentroid', 'MeanIntensity', 'Area');
points = zeros(numel(stats), 3);
for i = 1:numel(stats)
    points(i, 1:2) = stats(i).WeightedCentroid;
    points(i, 3) = stats(i).MeanIntensity .* sqrt(stats(i).Area);
end
points = points(all(isfinite(points), 2), :);
if size(points, 1) > 260
    [~, order] = sort(points(:, 3), 'descend');
    points = points(order(1:260), :);
end
if size(points, 1) < 5
    blades = empty_blades();
    return
end
points(:, 3) = points(:, 3) ./ max(points(:, 3), eps);

[blade, inliers, bestScore] = ransac_soma_chain(points, structural, ridge, ...
    processCandidates, 900, []);
blades = empty_blades();
if ~isempty(blade)
    blades(end + 1) = blade;
end

remaining = find(~inliers);
if numel(remaining) >= 5
    [secondBlade, ~, secondScore] = ransac_soma_chain( ...
        points, structural, ridge, processCandidates, 500, remaining);
    if ~isempty(secondBlade) && secondScore >= 0.58 * bestScore && ...
            curves_are_distinct(blade, secondBlade)
        blades(end + 1) = secondBlade;
    end
end
end

function [blade, inlierMask, bestScore] = ransac_soma_chain(points, structural, ...
        ridge, processCandidates, nTrials, allowedIndices)
blade = [];
inlierMask = false(size(points, 1), 1);
bestScore = -inf;
if isempty(allowedIndices)
    allowedIndices = 1:size(points, 1);
end
if numel(allowedIndices) < 3
    return
end
stream = RandStream('mt19937ar', 'Seed', 1907 + numel(allowedIndices));
for trial = 1:nTrials
    selected = allowedIndices(randperm(stream, numel(allowedIndices), 3));
    sample = points(selected, 1:2);
    sampleCenter = mean(sample, 1);
    sampleCoords = sample - sampleCenter;
    [~, ~, v] = svd(sampleCoords, 0);
    tangent = v(:, 1)';
    normal = [-tangent(2), tangent(1)];
    tSample = sampleCoords * tangent';
    if range(tSample) < 65
        continue
    end
    rSample = sampleCoords * normal';
    X = [ones(3, 1), tSample, tSample .^ 2];
    if rcond(X) < 1e-9
        continue
    end
    beta = X \ rSample;
    allCoords = points(:, 1:2) - sampleCenter;
    tAll = allCoords * tangent';
    rAll = allCoords * normal';
    residual = abs(rAll - (beta(1) + beta(2) .* tAll + beta(3) .* tAll .^ 2));
    candidateInliers = residual <= 13;
    candidateInliers(setdiff(1:size(points, 1), allowedIndices)) = false;
    if nnz(candidateInliers) < 5
        continue
    end
    tIn = tAll(candidateInliers);
    span = prctile(tIn, 95) - prctile(tIn, 5);
    if span < 90
        continue
    end
    occupancy = chain_occupancy(tIn, 9);
    curvaturePenalty = exp(-abs(beta(3)) * span * 2.5);
    weightedSupport = sum(points(candidateInliers, 3));
    score = weightedSupport * (0.35 + 0.65 * min(span / 330, 1)) * ...
        (0.45 + 0.55 * occupancy) * curvaturePenalty;
    if score > bestScore
        bestScore = score;
        inlierMask = candidateInliers;
    end
end

if ~any(inlierMask)
    return
end
blade = refine_soma_chain(points(inlierMask, :), structural, ridge, ...
    processCandidates, bestScore);
if isempty(blade)
    inlierMask(:) = false;
    bestScore = -inf;
end
end

function occupancy = chain_occupancy(t, nBins)
edges = linspace(min(t), max(t) + eps, nBins + 1);
counts = histcounts(t, edges);
occupancy = nnz(counts > 0) ./ nBins;
end

function blade = refine_soma_chain(points, structural, ridge, processCandidates, rawScore)
blade = [];
weights = points(:, 3);
weights = weights ./ max(sum(weights), eps);
center = sum(points(:, 1:2) .* weights, 1);
coords = points(:, 1:2) - center;
covar = coords' * (coords .* weights);
[vec, val] = eig(covar);
[~, majorIdx] = max(diag(val));
tangent = vec(:, majorIdx)';
if tangent(1) < 0
    tangent = -tangent;
end
normal = [-tangent(2), tangent(1)];
t = coords * tangent';
r = coords * normal';
X = [ones(numel(t), 1), t, t .^ 2];
beta = robust_weighted_fit(X, r, weights);
residual = r - X * beta;
residualPx = sqrt(sum(weights .* residual .^ 2));
tLo = weighted_percentile(t, weights, 2) - 25;
tHi = weighted_percentile(t, weights, 98) + 25;
span = tHi - tLo;
if span < 90
    return
end
tSample = linspace(tLo, tHi, 260)';
rSample = beta(1) + beta(2) .* tSample + beta(3) .* tSample .^ 2;
curveXY = center + tSample .* tangent + rSample .* normal;
drdt = beta(2) + 2 .* beta(3) .* tSample;
curveTangent = tangent + drdt .* normal;
curveTangent = curveTangent ./ max(vecnorm(curveTangent, 2, 2), eps);
curveNormal = [-curveTangent(:, 2), curveTangent(:, 1)];

[positiveSupport, negativeSupport] = side_support(curveXY, curveNormal, structural, ridge);
if negativeSupport > positiveSupport
    curveNormal = -curveNormal;
end
sideAsymmetry = abs(positiveSupport - negativeSupport) ./ ...
    max(positiveSupport + negativeSupport, eps);
perpendicularSupport = process_perpendicular_support(curveXY, curveNormal, processCandidates);

spanTerm = min(1, span / 320);
pointTerm = min(1, size(points, 1) / 16);
residualTerm = exp(-residualPx / 18);
scoreTerm = min(1, rawScore / 8);
confidence = 0.28 .* spanTerm + 0.20 .* pointTerm + 0.18 .* residualTerm + ...
    0.12 .* sideAsymmetry + 0.14 .* perpendicularSupport + 0.08 .* scoreTerm;
blade = struct('curve', [curveXY curveNormal], ...
    'radial', mean_unit(curveNormal), 'tangent', tangent, ...
    'confidence', min(max(confidence, 0), 1), ...
    'residualPx', residualPx, 'candidateIndices', []);
end

function support = process_perpendicular_support(curveXY, curveNormal, candidates)
if isempty(candidates)
    support = 0;
    return
end
values = zeros(numel(candidates), 1);
weights = zeros(numel(candidates), 1);
for i = 1:numel(candidates)
    midpoint = 0.5 .* (candidates(i).base + candidates(i).distal);
    [d2, k] = min(sum((curveXY - midpoint) .^ 2, 2));
    if d2 <= 75 ^ 2
        values(i) = abs(dot(candidates(i).radial, curveNormal(k, :)));
        weights(i) = candidates(i).confidence;
    end
end
if sum(weights) <= 0
    support = 0;
else
    support = sum(values .* weights) ./ sum(weights);
end
end

function tf = curves_are_distinct(a, b)
if isempty(a) || isempty(b)
    tf = false;
    return
end
distances = zeros(size(a.curve, 1), 1);
for i = 1:size(a.curve, 1)
    distances(i) = sqrt(min(sum((b.curve(:, 1:2) - a.curve(i, 1:2)) .^ 2, 2)));
end
tf = median(distances) > 35;
end

function roi = load_roi_evidence(metadata, structural)
roi = struct('cellBoundaries', {{}}, 'somaBoundaries', {{}}, ...
    'cellMasks', {{}}, 'somaMasks', {{}}, 'pairs', zeros(0, 2), ...
    'isAvailable', false);
if ~isfield(metadata, 'cellRoiPath') || ~isfield(metadata, 'somaRoiPath') || ...
        ~isfile(metadata.cellRoiPath) || ~isfile(metadata.somaRoiPath)
    return
end
roi.cellBoundaries = read_aqua_roi_boundaries(metadata.cellRoiPath);
roi.somaBoundaries = read_aqua_roi_boundaries(metadata.somaRoiPath);
if isempty(roi.cellBoundaries) || isempty(roi.somaBoundaries)
    return
end
roi.cellMasks = boundaries_to_masks(roi.cellBoundaries, size(structural));
roi.somaMasks = boundaries_to_masks(roi.somaBoundaries, size(structural));
cellCenters = mask_centers(roi.cellMasks);
somaCenters = mask_centers(roi.somaMasks);
roi.pairs = pair_cells_and_somas(roi.cellBoundaries, cellCenters, somaCenters);
roi.isAvailable = ~isempty(roi.pairs);
end

function scale = estimate_morphology_scale(roi, imageScale)
scale = imageScale;
if ~roi.isAvailable || isempty(roi.somaMasks)
    return
end
areas = cellfun(@nnz, roi.somaMasks);
areas = areas(areas >= 12);
if isempty(areas)
    return
end
equivalentDiameter = 2 .* sqrt(median(areas) ./ pi);
scale = min(max(equivalentDiameter ./ 18, 0.65), 1.8);
scale = 0.75 .* scale + 0.25 .* imageScale;
end

function boundaries = read_aqua_roi_boundaries(path)
S = load(path, 'bd0');
boundaries = {};
if ~isfield(S, 'bd0')
    return
end
for i = 1:numel(S.bd0)
    record = S.bd0{i};
    coords = [];
    if iscell(record) && ~isempty(record) && iscell(record{1}) && ...
            ~isempty(record{1}) && isnumeric(record{1}{1})
        coords = double(record{1}{1});
    elseif isnumeric(record) && size(record, 2) == 2
        coords = double(record);
    end
    if size(coords, 2) == 2 && size(coords, 1) >= 3
        % AQuA bd0 stores boundary coordinates as [row, column].
        boundaries{end + 1} = coords(:, [2 1]); %#ok<AGROW>
    end
end
end

function masks = boundaries_to_masks(boundaries, imgSize)
masks = cell(size(boundaries));
for i = 1:numel(boundaries)
    b = boundaries{i};
    masks{i} = poly2mask(b(:, 1), b(:, 2), imgSize(1), imgSize(2));
end
end

function centers = mask_centers(masks)
centers = nan(numel(masks), 2);
for i = 1:numel(masks)
    [yy, xx] = find(masks{i});
    if ~isempty(xx)
        centers(i, :) = [mean(xx), mean(yy)];
    end
end
end

function pairs = pair_cells_and_somas(cellBoundaries, cellCenters, somaCenters)
pairs = zeros(0, 2);
used = false(size(somaCenters, 1), 1);
for c = 1:size(cellCenters, 1)
    if any(~isfinite(cellCenters(c, :)))
        continue
    end
    inside = false(size(somaCenters, 1), 1);
    boundary = cellBoundaries{c};
    for s = 1:size(somaCenters, 1)
        if ~used(s) && all(isfinite(somaCenters(s, :)))
            inside(s) = inpolygon(somaCenters(s, 1), somaCenters(s, 2), ...
                boundary(:, 1), boundary(:, 2));
        end
    end
    choices = find(inside);
    if isempty(choices)
        choices = find(~used & all(isfinite(somaCenters), 2));
    end
    if isempty(choices)
        continue
    end
    [~, k] = min(vecnorm(somaCenters(choices, :) - cellCenters(c, :), 2, 2));
    s = choices(k);
    pairs(end + 1, :) = [c s]; %#ok<AGROW>
    used(s) = true;
end
end

function candidates = roi_process_candidates(roi, structural, scale)
candidates = empty_process_candidates();
if ~roi.isAvailable
    return
end
for i = 1:size(roi.pairs, 1)
    cellIdx = roi.pairs(i, 1);
    somaIdx = roi.pairs(i, 2);
    cellMask = roi.cellMasks{cellIdx};
    somaMask = roi.somaMasks{somaIdx};
    [sy, sx] = find(somaMask);
    [cy, cx] = find(cellMask);
    if numel(sx) < 3 || numel(cx) < 20
        continue
    end
    soma = [mean(sx), mean(sy)];
    coords = [double(cx), double(cy)];
    rel = coords - soma;
    weights = double(structural(cellMask));
    weights = 0.35 + 0.65 .* weights ./ max(max(weights), eps);
    covar = (rel' * (rel .* weights)) ./ max(sum(weights), eps);
    [vec, val] = eig(covar);
    [eigenvalues, order] = sort(diag(val), 'descend');
    major = vec(:, order(1))';
    cellCentroidDirection = sum(rel .* weights, 1) ./ max(sum(weights), eps);
    if norm(cellCentroidDirection) > 1 && dot(major, cellCentroidDirection) < 0
        major = -major;
    end
    projection = rel * major';
    positiveExtent = prctile(projection, 98);
    negativeExtent = abs(prctile(projection, 2));
    if negativeExtent > positiveExtent
        major = -major;
        projection = -projection;
        positiveExtent = prctile(projection, 98);
        negativeExtent = abs(prctile(projection, 2));
    end
    elongation = 1 - min(eigenvalues(2) / max(eigenvalues(1), eps), 1);
    polarity = positiveExtent ./ max(positiveExtent + negativeExtent, eps);
    somaSignal = mean(structural(somaMask), 'omitnan');
    ringRadius = max(5, round(12 * scale));
    ring = imdilate(somaMask, strel('disk', ringRadius)) & ~somaMask;
    ringSignal = mean(structural(ring), 'omitnan');
    signalSupport = min(max(0.5 + 0.5 .* (somaSignal - ringSignal), 0.35), 1);
    confidence = min(0.88, (0.30 + 0.38 .* elongation + 0.20 .* polarity + ...
        0.12 .* signalSupport));
    lengthValue = max(positiveExtent, 25 * scale);
    distal = soma + lengthValue .* major;
    candidates(end + 1) = struct('base', soma, 'distal', distal, ...
        'radial', major, 'length', lengthValue, 'confidence', confidence, ...
        'ridgeEvidence', elongation, 'branchAsymmetry', polarity); %#ok<AGROW>
end
end

function blades = reconcile_roi_and_fluorescence_blades(roiBlades, fluorescenceBlades)
blades = roiBlades;
if isempty(fluorescenceBlades)
    for i = 1:numel(blades)
        blades(i).confidence = 0.82 .* blades(i).confidence;
    end
    return
end
for i = 1:numel(blades)
    distances = inf(numel(fluorescenceBlades), 1);
    directionAgreement = zeros(numel(fluorescenceBlades), 1);
    for j = 1:numel(fluorescenceBlades)
        distances(j) = median_curve_distance(blades(i).curve, fluorescenceBlades(j).curve);
        directionAgreement(j) = abs(dot(blades(i).radial, fluorescenceBlades(j).radial));
    end
    [distance, j] = min(distances);
    agreement = exp(-distance / 75) .* directionAgreement(j);
    blades(i).confidence = blades(i).confidence .* (0.62 + 0.38 .* agreement);
end
end

function distance = median_curve_distance(a, b)
sampleIdx = round(linspace(1, size(a, 1), min(80, size(a, 1))));
values = zeros(numel(sampleIdx), 1);
for i = 1:numel(sampleIdx)
    p = a(sampleIdx(i), 1:2);
    values(i) = sqrt(min(sum((b(:, 1:2) - p) .^ 2, 2)));
end
distance = median(values);
end

function blades = prune_blade_hypotheses(blades, imgSize, scale)
if numel(blades) < 2
    return
end
[~, order] = sort([blades.confidence], 'descend');
blades = blades(order);
keep = true(size(blades));
distanceFloor = max(24 * scale, 0.045 * hypot(imgSize(1), imgSize(2)));
for i = 1:numel(blades)
    if ~keep(i)
        continue
    end
    for j = i + 1:numel(blades)
        if ~keep(j)
            continue
        end
        distance = 0.5 .* (median_curve_distance(blades(i).curve, blades(j).curve) + ...
            median_curve_distance(blades(j).curve, blades(i).curve));
        tangentAgreement = abs(dot(blades(i).tangent, blades(j).tangent));
        endpointDistance = blade_endpoint_distance(blades(i), blades(j));
        if all(endpointDistance < 60 * scale, 'all') && ...
                all(tangentAgreement > cosd(20), 'all')
            blades(i) = merge_collinear_blades(blades(i), blades(j), imgSize, scale);
            keep(j) = false;
        elseif all(distance < distanceFloor, 'all') && ...
                all(tangentAgreement > cosd(15), 'all')
            keep(j) = false;
        end
    end
end
blades = blades(keep);
end

function distance = blade_endpoint_distance(a, b)
endpointsA = [a.curve(1, 1:2); a.curve(end, 1:2)];
endpointsB = [b.curve(1, 1:2); b.curve(end, 1:2)];
distance = inf;
for i = 1:2
    distance = min(distance, min(vecnorm(endpointsB - endpointsA(i, :), 2, 2)));
end
end

function blade = merge_collinear_blades(a, b, imgSize, scale)
source = [a b];
template = struct('base', [], 'distal', [], 'radial', [], 'length', [], ...
    'confidence', [], 'ridgeEvidence', [], 'branchAsymmetry', []);
candidates = repmat(template, 1, 14);
candidateIndex = 0;
for i = 1:2
    curve = source(i).curve;
    sampleRows = unique(round(linspace(1, size(curve, 1), 7)));
    for k = sampleRows
        candidateIndex = candidateIndex + 1;
        radial = curve(k, 3:4);
        candidates(candidateIndex) = struct('base', curve(k, 1:2), ...
            'distal', curve(k, 1:2) + 60 * scale .* radial, ...
            'radial', radial, 'length', 60 * scale, ...
            'confidence', source(i).confidence, 'ridgeEvidence', 1, ...
            'branchAsymmetry', 1);
    end
end
candidates = candidates(1:candidateIndex);
blade = fit_one_blade(candidates, imgSize, scale);
blade.confidence = max([a.confidence b.confidence]) .* 0.96;
end

function blades = apply_gcl_support(blades, gclEvidence, scale)
if isempty(blades)
    return
end
offsets = scale .* [10 20 32];
for b = 1:numel(blades)
    curve = blades(b).curve;
    sampleRows = round(linspace(1, size(curve, 1), min(90, size(curve, 1))));
    values = zeros(numel(sampleRows), numel(offsets));
    for i = 1:numel(sampleRows)
        k = sampleRows(i);
        for j = 1:numel(offsets)
            point = curve(k, 1:2) + offsets(j) .* curve(k, 3:4);
            values(i, j) = sample_image(gclEvidence, point);
        end
    end
    support = median(values(:), 'omitnan');
    if ~isfinite(support)
        support = 0;
    end
    blades(b).confidence = blades(b).confidence .* (0.90 + 0.10 .* support);
end
end

function q = endpoint_complexity(point, structural, endpoints, branchpoints)
[xx, yy] = meshgrid(1:size(structural, 2), 1:size(structural, 1));
d2 = (xx - point(1)) .^ 2 + (yy - point(2)) .^ 2;
near = d2 <= 28 ^ 2;
inner = d2 <= 13 ^ 2;
mass = sum(structural(near), 'omitnan') ./ max(nnz(near), 1);
tipCount = nnz(endpoints & near);
branchCount = nnz(branchpoints & near);
core = sum(structural(inner), 'omitnan') ./ max(nnz(inner), 1);
q = double(0.45 .* mass + 0.20 .* core + 0.18 .* min(tipCount / 5, 1) + ...
    0.17 .* min(branchCount / 3, 1));
end

function idx = sample_line(p1, p2, imgSize)
n = max(2, ceil(norm(p2 - p1)));
x = round(linspace(p1(1), p2(1), n));
y = round(linspace(p1(2), p2(2), n));
x = min(max(x, 1), imgSize(2));
y = min(max(y, 1), imgSize(1));
idx = unique(sub2ind(imgSize, y, x));
end

function tf = duplicates_candidate(candidate, candidates)
tf = false;
for i = 1:numel(candidates)
    centerA = 0.5 .* (candidate.base + candidate.distal);
    centerB = 0.5 .* (candidates(i).base + candidates(i).distal);
    angleAgreement = abs(dot(candidate.radial, candidates(i).radial));
    if norm(centerA - centerB) < 22 && angleAgreement > cosd(18)
        tf = true;
        return
    end
end
end

function blades = fit_sgz_blades(candidates, imgSize, scale, refinement)
if nargin < 3 || isempty(scale)
    scale = min(max(min(imgSize) / 512, 0.65), 1.75);
end
if nargin < 4
    refinement = [];
end
blades = empty_blades();
if isempty(candidates)
    return
end
confidenceFloor = 0.16;
keep = [candidates.confidence] >= confidenceFloor;
candidates = candidates(keep);
if isempty(candidates)
    return
end

[groups, jointV] = spatial_soma_chain_groups(candidates, imgSize, scale);
if jointV.valid && ~isempty(refinement)
    [groups, jointV] = apply_sgz_anchor_refinement( ...
        candidates, groups, jointV, refinement, imgSize, scale);
end

for i = 1:numel(groups)
    if jointV.valid
        blade = fit_apex_constrained_blade(candidates(groups{i}), imgSize, ...
            scale, jointV.apex, jointV.tangents(i, :));
    else
        blade = fit_one_blade(candidates(groups{i}), imgSize, scale);
    end
    if ~isempty(blade)
        blade.candidateIndices = groups{i};
        blades(end + 1) = blade; %#ok<AGROW>
    end
end
blades = prune_blade_hypotheses(blades, imgSize, scale);
end

function [groups, model] = apply_sgz_anchor_refinement( ...
        candidates, groups, model, refinement, imgSize, scale)
required = {'referencePoints', 'includePoint', 'excludePoint'};
if numel(groups) ~= 2 || ~all(isfield(refinement, required))
    return
end
bases = vertcat(candidates.base);
reference = refinement.referencePoints;
groupDistance = inf(1, 2);
for i = 1:2
    distances = zeros(size(reference, 1), 1);
    for j = 1:size(reference, 1)
        distances(j) = min(vecnorm(bases(groups{i}, :) - reference(j, :), 2, 2));
    end
    groupDistance(i) = median(distances);
end
[~, target] = min(groupDistance);
other = 3 - target;
[~, includeIndex] = min(vecnorm(bases - refinement.includePoint, 2, 2));
[~, excludeIndex] = min(vecnorm(bases - refinement.excludePoint, 2, 2));
if norm(bases(includeIndex, :) - refinement.includePoint) > 18 * scale || ...
        norm(bases(excludeIndex, :) - refinement.excludePoint) > 18 * scale
    return
end
groups{target} = unique([setdiff(groups{target}, excludeIndex), includeIndex]);
groups{other} = setdiff(groups{other}, includeIndex);
if any(cellfun(@numel, groups) < 2)
    return
end
refinedModel = shared_apex_from_fixed_groups(candidates, groups, imgSize, scale);
if refinedModel.valid
    model = refinedModel;
end
end

function model = shared_apex_from_fixed_groups(candidates, groups, imgSize, ~)
model = struct('valid', false, 'apex', [nan nan], ...
    'tangents', zeros(0, 2), 'score', -inf);
bases = vertcat(candidates.base);
confidence = [candidates.confidence]';
centers = zeros(2, 2);
tangents = zeros(2, 2);
for i = 1:2
    idx = groups{i};
    weights = confidence(idx) ./ max(sum(confidence(idx)), eps);
    centers(i, :) = sum(bases(idx, :) .* weights, 1);
    coords = bases(idx, :) - centers(i, :);
    covariance = coords' * (coords .* weights);
    [vectors, values] = eig(covariance);
    [~, major] = max(diag(values));
    tangents(i, :) = vectors(:, major)';
end
system = [tangents(1, :)', -tangents(2, :)'];
if rcond(system) < 0.08
    return
end
parameters = system \ (centers(2, :) - centers(1, :))';
apex = centers(1, :) + parameters(1) .* tangents(1, :);
margin = 0.80 * hypot(imgSize(1), imgSize(2));
if apex(1) < 1 - margin || apex(1) > imgSize(2) + margin || ...
        apex(2) < 1 - margin || apex(2) > imgSize(1) + margin
    return
end
for i = 1:2
    if mean((bases(groups{i}, :) - apex) * tangents(i, :)') < 0
        tangents(i, :) = -tangents(i, :);
    end
end
opening = acosd(min(max(dot(tangents(1, :), tangents(2, :)), -1), 1));
acuteOpening = min(opening, 180 - opening);
if acuteOpening < 20 || acuteOpening > 82
    return
end
model = struct('valid', true, 'apex', apex, 'tangents', tangents, ...
    'score', 1);
end

function [groups, jointV] = spatial_soma_chain_groups(candidates, imgSize, scale)
groups = {1:numel(candidates)};
jointV = struct('valid', false, 'apex', [nan nan], ...
    'tangents', zeros(0, 2), 'score', -inf);
if numel(candidates) < 4
    return
end
[firstGroup, firstScore, firstTangent] = best_soma_line_group( ...
    candidates, 1:numel(candidates), scale);
if numel(firstGroup) < 2
    return
end
remaining = setdiff(1:numel(candidates), firstGroup);
bases = vertcat(candidates.base);
confidence = [candidates.confidence]';
groupWeights = confidence(firstGroup) ./ max(sum(confidence(firstGroup)), eps);
groupCenter = sum(bases(firstGroup, :) .* groupWeights, 1);
firstNormal = [-firstTangent(2), firstTangent(1)];
nearFirstChain = abs((bases - groupCenter) * firstNormal') <= 50 * scale;
remaining = remaining(~nearFirstChain(remaining));
if numel(remaining) < 2
    groups = {firstGroup};
    return
end
[secondGroup, secondScore, secondTangent] = best_soma_line_group( ...
    candidates, remaining, scale);
angle = acosd(min(max(abs(dot(firstTangent, secondTangent)), -1), 1));
if numel(secondGroup) >= 2 && secondScore >= 0.14 * firstScore && angle >= 20
    [jointGroups, jointV] = best_shared_apex_groups(candidates, imgSize, scale);
    if jointV.valid
        groups = jointGroups;
    else
        groups = {firstGroup};
    end
else
    groups = {1:numel(candidates)};
end
end

function [bestGroups, bestModel] = best_shared_apex_groups(candidates, imgSize, scale)
bestGroups = {};
bestModel = struct('valid', false, 'apex', [nan nan], ...
    'tangents', zeros(0, 2), 'score', -inf);
bases = vertcat(candidates.base);
radials = vertcat(candidates.radial);
confidence = [candidates.confidence]';
n = size(bases, 1);
distanceLimit = 25 * scale;
backwardTolerance = 30 * scale;
margin = 0.80 * hypot(imgSize(1), imgSize(2));

hypotheses = repmat(struct('point', [0 0], 'tangent', [1 0]), nchoosek(n, 2), 1);
hypothesisCount = 0;
for i = 1:n - 1
    for j = i + 1:n
        delta = bases(j, :) - bases(i, :);
        if norm(delta) < 32 * scale
            continue
        end
        hypothesisCount = hypothesisCount + 1;
        hypotheses(hypothesisCount) = struct('point', bases(i, :), ...
            'tangent', delta ./ norm(delta));
    end
end
hypotheses = hypotheses(1:hypothesisCount);

for a = 1:numel(hypotheses) - 1
    for b = a + 1:numel(hypotheses)
        lineA = hypotheses(a);
        lineB = hypotheses(b);
        angle = acosd(min(max(abs(dot(lineA.tangent, lineB.tangent)), 0), 1));
        if angle < 20 || angle > 82
            continue
        end
        system = [lineA.tangent(:), -lineB.tangent(:)];
        if rcond(system) < 0.08
            continue
        end
        parameters = system \ (lineB.point - lineA.point)';
        apex = lineA.point + parameters(1) .* lineA.tangent;
        if apex(1) < 1 - margin || apex(1) > imgSize(2) + margin || ...
                apex(2) < 1 - margin || apex(2) > imgSize(1) + margin
            continue
        end
        for signA = [-1 1]
            for signB = [-1 1]
                tangents = [signA .* lineA.tangent; signB .* lineB.tangent];
                [groups, score] = score_shared_apex_model(apex, tangents, ...
                    bases, radials, confidence, distanceLimit, ...
                    backwardTolerance, scale, imgSize);
                if score > bestModel.score
                    bestGroups = groups;
                    bestModel = struct('valid', true, 'apex', apex, ...
                        'tangents', tangents, 'score', score);
                end
            end
        end
    end
end

if ~bestModel.valid
    return
end
support = cellfun(@numel, bestGroups);
spans = zeros(1, 2);
for i = 1:2
    spans(i) = range((bases(bestGroups{i}, :) - bestModel.apex) * ...
        bestModel.tangents(i, :)');
end
if any(support < 2) || any(spans < 32 * scale)
    bestGroups = {};
    bestModel.valid = false;
end
end

function [groups, score] = score_shared_apex_model(apex, tangents, bases, ...
        radials, confidence, distanceLimit, backwardTolerance, scale, imgSize)
n = size(bases, 1);
cost = inf(n, 2);
perpendicular = inf(n, 2);
forward = zeros(n, 2);
axisQuality = zeros(n, 2);
for i = 1:2
    tangent = tangents(i, :);
    normal = [-tangent(2), tangent(1)];
    delta = bases - apex;
    forward(:, i) = delta * tangent';
    perpendicular(:, i) = abs(delta * normal');
    behind = max(-forward(:, i), 0);
    cost(:, i) = perpendicular(:, i) + 1.8 .* behind;
    axisQuality(:, i) = abs(radials * normal');
end
[~, assignment] = min(cost, [], 2);
groups = cell(1, 2);
score = 0;
for i = 1:2
    inlier = assignment == i & perpendicular(:, i) <= distanceLimit & ...
        forward(:, i) >= -backwardTolerance;
    groups{i} = find(inlier)';
    if numel(groups{i}) < 2
        score = -inf;
        return
    end
    t = forward(groups{i}, i);
    span = range(t);
    if span < 32 * scale
        score = -inf;
        return
    end
    % Soma positions define the blade. Process orientation is only a weak
    % tie-breaker because basal branches can point along the SGZ.
    support = sum(confidence(groups{i}) .* ...
        (0.85 + 0.15 .* axisQuality(groups{i}, i)) .* ...
        exp(-0.5 .* (perpendicular(groups{i}, i) ./ distanceLimit) .^ 2));
    occupancy = chain_occupancy(t, min(6, max(2, numel(t))));
    score = score + support .* (0.35 + 0.65 .* min(span / (220 * scale), 1)) .* ...
        (0.55 + 0.45 .* occupancy);
end

angle = acosd(min(max(dot(tangents(1, :), tangents(2, :)), -1), 1));
acuteAngle = min(angle, 180 - angle);
angleTerm = exp(-0.5 .* ((acuteAngle - 48) ./ 28) .^ 2);
outsideX = max([1 - apex(1), 0, apex(1) - imgSize(2)]);
outsideY = max([1 - apex(2), 0, apex(2) - imgSize(1)]);
outsideDistance = hypot(outsideX, outsideY);
apexTerm = exp(-outsideDistance / (0.55 .* hypot(imgSize(1), imgSize(2))));
balance = min(numel(groups{1}), numel(groups{2})) / ...
    max(numel(groups{1}), numel(groups{2}));
score = score .* (0.58 + 0.22 .* angleTerm + 0.12 .* apexTerm + 0.08 .* balance);
end

function [bestGroup, bestScore, bestTangent] = best_soma_line_group( ...
        candidates, allowed, scale)
bestGroup = [];
bestScore = -inf;
bestTangent = [1 0];
bases = vertcat(candidates.base);
radials = vertcat(candidates.radial);
confidence = [candidates.confidence]';
distanceLimit = 22 * scale;
for a = 1:numel(allowed) - 1
    for b = a + 1:numel(allowed)
        i = allowed(a);
        j = allowed(b);
        delta = bases(j, :) - bases(i, :);
        spanSeed = norm(delta);
        if spanSeed < 28 * scale
            continue
        end
        tangent = delta ./ spanSeed;
        normal = [-tangent(2), tangent(1)];
        residual = abs((bases - bases(i, :)) * normal');
        inliers = allowed(residual(allowed) <= distanceLimit);
        if numel(inliers) < 2
            continue
        end
        t = (bases(inliers, :) - bases(i, :)) * tangent';
        span = range(t);
        if span < 28 * scale
            continue
        end
        axisQuality = abs(radials(inliers, :) * normal');
        weightedSupport = sum(confidence(inliers) .* (0.45 + 0.55 .* axisQuality));
        occupancy = chain_occupancy(t, min(6, max(2, numel(inliers))));
        score = weightedSupport .* (0.35 + 0.65 .* min(span / (230 * scale), 1)) .* ...
            (0.55 + 0.45 .* occupancy);
        score = double(score(1));
        if score > bestScore
            bestGroup = inliers;
            bestScore = score;
            bestTangent = tangent;
        end
    end
end
end

function blades = empty_blades()
blades = struct('curve', {}, 'radial', {}, 'tangent', {}, 'confidence', {}, ...
    'residualPx', {}, 'candidateIndices', {});
end

function blade = fit_one_blade(candidates, imgSize, scale)
bases = vertcat(candidates.base);
radials = vertcat(candidates.radial);
weights = [candidates.confidence]' .* sqrt([candidates.length]');
weights = weights ./ max(sum(weights), eps);

% The soma chain defines the blade tangent. Process axes only choose side.
if size(bases, 1) >= 2
    provisionalCenter = sum(bases .* weights, 1);
    provisionalCoords = bases - provisionalCenter;
    covar = provisionalCoords' * (provisionalCoords .* weights);
    [vec, val] = eig(covar);
    [~, majorIdx] = max(diag(val));
    tangent = vec(:, majorIdx)';
    if range(provisionalCoords * tangent') < 18 * scale
        rawRadial = mean_unit(radials, weights);
        tangent = [-rawRadial(2), rawRadial(1)];
    end
else
    rawRadial = mean_unit(radials, weights);
    tangent = [-rawRadial(2), rawRadial(1)];
end
if tangent(1) < 0
    tangent = -tangent;
end
normal = [-tangent(2), tangent(1)];
normalProjection = radials * normal';
sideVote = sum(weights .* normalProjection);
if abs(sideVote) < 0.05
    sideVote = dot(mean_unit(radials, weights), normal);
end
if sideVote < 0
    normal = -normal;
    normalProjection = -normalProjection;
end
avgRadial = normal;
axisQuality = abs(normalProjection);
fitWeights = weights .* (0.18 + 0.82 .* axisQuality .^ 2);
fitWeights = fitWeights ./ max(sum(fitWeights), eps);

center = sum(bases .* fitWeights, 1);
rel = bases - center;
t = rel * tangent';
r = rel * avgRadial';

if numel(t) >= 5
    degree = 2;
elseif numel(t) >= 2
    degree = 1;
else
    degree = 0;
end
X = ones(numel(t), degree + 1);
for p = 1:degree
    X(:, p + 1) = t .^ p;
end
beta = robust_weighted_fit(X, r, fitWeights);
residual = r - X * beta;
residualPx = sqrt(sum(fitWeights .* residual .^ 2));

if isscalar(t)
    halfSpan = min(165 * scale, 0.42 .* hypot(imgSize(1), imgSize(2)));
    tSample = linspace(-halfSpan, halfSpan, 180)';
else
    tMin = min(t) - 45 * scale;
    tMax = max(t) + 45 * scale;
    if tMax - tMin < 150 * scale
        mid = 0.5 .* (tMin + tMax);
        tMin = mid - 75 * scale;
        tMax = mid + 75 * scale;
    end
    tSample = linspace(tMin, tMax, 220)';
end
R = ones(numel(tSample), degree + 1);
for p = 1:degree
    R(:, p + 1) = tSample .^ p;
end
rSample = R * beta;
curveXY = center + tSample .* tangent + rSample .* avgRadial;

if degree == 0
    drdt = zeros(size(tSample));
elseif degree == 1
    drdt = repmat(beta(2), size(tSample));
else
    drdt = beta(2) + 2 .* beta(3) .* tSample;
end
curveTangent = tangent + drdt .* avgRadial;
curveTangent = curveTangent ./ max(vecnorm(curveTangent, 2, 2), eps);
curveRadial = [-curveTangent(:, 2), curveTangent(:, 1)];
flip = curveRadial * avgRadial' < 0;
curveRadial(flip, :) = -curveRadial(flip, :);

support = min(1, numel(candidates) / 3);
directionAgreement = sum(weights .* axisQuality);
baseSpan = max(range(t), 0);
spanSupport = min(1, baseSpan / max(150 * scale, eps));
confidence = mean([candidates.confidence]) .* (0.38 + 0.62 .* support) .* ...
    (0.35 + 0.65 .* directionAgreement) .* (0.55 + 0.45 .* spanSupport) .* ...
    exp(-residualPx / max(24 * scale, eps));
blade = struct('curve', [curveXY curveRadial], 'radial', avgRadial, ...
    'tangent', tangent, 'confidence', min(max(confidence, 0), 1), ...
    'residualPx', residualPx, 'candidateIndices', 1:numel(candidates));
end

function blade = fit_apex_constrained_blade(candidates, ~, scale, apex, tangent)
bases = vertcat(candidates.base);
radials = vertcat(candidates.radial);
weights = [candidates.confidence]' .* sqrt([candidates.length]');
weights = weights ./ max(sum(weights), eps);
tangent = tangent ./ max(norm(tangent), eps);
if sum(weights .* ((bases - apex) * tangent')) < 0
    tangent = -tangent;
end
normal = [-tangent(2), tangent(1)];
normalProjection = radials * normal';
if sum(weights .* normalProjection) < 0
    normal = -normal;
    normalProjection = -normalProjection;
end
axisQuality = abs(normalProjection);
fitWeights = weights .* (0.18 + 0.82 .* axisQuality .^ 2);
fitWeights = fitWeights ./ max(sum(fitWeights), eps);

delta = bases - apex;
t = delta * tangent';
r = delta * normal';
degree = 1 + double(numel(t) >= 5);
X = t;
if degree == 2
    X = [X t .^ 2];
end
beta = robust_weighted_fit(X, r, fitWeights);
residual = r - X * beta;
residualPx = sqrt(sum(fitWeights .* residual .^ 2));

tMax = max(t) + 45 * scale;
if tMax < 150 * scale
    tMax = 150 * scale;
end
tSample = linspace(0, tMax, 240)';
R = tSample;
if degree == 2
    R = [R tSample .^ 2];
end
rSample = R * beta;
curveXY = apex + tSample .* tangent + rSample .* normal;
if degree == 1
    drdt = repmat(beta(1), size(tSample));
else
    drdt = beta(1) + 2 .* beta(2) .* tSample;
end
curveTangent = tangent + drdt .* normal;
curveTangent = curveTangent ./ max(vecnorm(curveTangent, 2, 2), eps);
curveRadial = [-curveTangent(:, 2), curveTangent(:, 1)];
flip = curveRadial * normal' < 0;
curveRadial(flip, :) = -curveRadial(flip, :);

support = min(1, numel(candidates) / 3);
directionAgreement = sum(weights .* axisQuality);
spanSupport = min(1, range(t) / max(150 * scale, eps));
confidence = mean([candidates.confidence]) .* (0.38 + 0.62 .* support) .* ...
    (0.35 + 0.65 .* directionAgreement) .* (0.55 + 0.45 .* spanSupport) .* ...
    exp(-residualPx / max(24 * scale, eps));
blade = struct('curve', [curveXY curveRadial], 'radial', normal, ...
    'tangent', tangent, 'confidence', min(max(confidence, 0), 1), ...
    'residualPx', residualPx, 'candidateIndices', 1:numel(candidates));
end

function [apex, visible, confidence] = shared_apex_summary(blades, imgSize)
apex = [nan nan];
visible = false;
confidence = 0;
if numel(blades) ~= 2
    return
end
endpoints = cell(1, 2);
for i = 1:2
    endpoints{i} = [blades(i).curve(1, 1:2); blades(i).curve(end, 1:2)];
end
distance = inf;
pair = [1 1];
for i = 1:2
    for j = 1:2
        thisDistance = norm(endpoints{1}(i, :) - endpoints{2}(j, :));
        if thisDistance < distance
            distance = thisDistance;
            pair = [i j];
        end
    end
end
if distance > 4
    return
end
apex = 0.5 .* (endpoints{1}(pair(1), :) + endpoints{2}(pair(2), :));
visible = apex(1) >= 1 && apex(1) <= imgSize(2) && ...
    apex(2) >= 1 && apex(2) <= imgSize(1);
confidence = min([blades.confidence]) .* exp(-distance / 4);
end

function beta = robust_weighted_fit(X, y, baseWeights)
y = y(:);
baseWeights = baseWeights(:);
if numel(baseWeights) ~= size(X, 1)
    error('Regression weights have %d values for %d observations.', ...
        numel(baseWeights), size(X, 1));
end
w = baseWeights ./ max(sum(baseWeights), eps);
beta = weighted_linear_fit(X, y, w);
if size(X, 1) <= size(X, 2) + 1
    return
end
for iter = 1:5
    residual = y - X * beta;
    scale = 1.4826 .* median(abs(residual - median(residual))) + eps;
    u = residual ./ (4.685 .* scale);
    robustWeights = ((abs(u) < 1) .* (1 - u .^ 2) .^ 2);
    combined = w .* robustWeights(:);
    if nnz(combined > 0) < size(X, 2)
        break
    end
    beta = weighted_linear_fit(X, y, combined);
end
end

function beta = weighted_linear_fit(X, y, weights)
sqrtWeights = sqrt(max(weights(:), 0));
weightedX = X .* sqrtWeights;
weightedY = y(:) .* sqrtWeights;
beta = pinv(weightedX) * weightedY;
end

function v = mean_unit(vectors, weights)
if nargin < 2 || isempty(weights)
    weights = ones(size(vectors, 1), 1) ./ size(vectors, 1);
end
v = sum(vectors .* weights(:), 1);
v = v ./ max(norm(v), eps);
end

function [labels, signedDistance, confidence, radialMap] = ...
        infer_layers(imgSize, blades, theta, coherence)
labels = zeros(imgSize, 'uint8');
signedDistance = nan(imgSize, 'single');
confidence = zeros(imgSize, 'single');
radialMap = nan([imgSize 2], 'single');
if isempty(blades)
    return
end

[xx, yy] = meshgrid(single(1:imgSize(2)), single(1:imgSize(1)));
bestD2 = inf(imgSize, 'single');
bestBladeConfidence = zeros(imgSize, 'single');
for b = 1:numel(blades)
    curve = blades(b).curve;
    for k = 1:3:size(curve, 1)
        dx = xx - single(curve(k, 1));
        dy = yy - single(curve(k, 2));
        d2 = dx .^ 2 + dy .^ 2;
        update = d2 < bestD2;
        if ~any(update(:))
            continue
        end
        radial = single(curve(k, 3:4));
        thisSigned = dx .* radial(1) + dy .* radial(2);
        signedDistance(update) = thisSigned(update);
        rx = radialMap(:, :, 1);
        ry = radialMap(:, :, 2);
        rx(update) = radial(1);
        ry(update) = radial(2);
        radialMap(:, :, 1) = rx;
        radialMap(:, :, 2) = ry;
        bestD2(update) = d2(update);
        bestBladeConfidence(update) = single(blades(b).confidence);
    end
end

valid = sqrt(bestD2) <= 185;
labels(valid & signedDistance < -14) = 1;                 % hilus
labels(valid & signedDistance >= -14 & signedDistance <= 8) = 2; % SGZ
labels(valid & signedDistance > 8 & signedDistance <= 46) = 3;   % GCL
labels(valid & signedDistance > 46) = 4;                 % molecular layer

localAxisX = cos(theta);
localAxisY = sin(theta);
orientationAgreement = abs(localAxisX .* radialMap(:, :, 1) + ...
    localAxisY .* radialMap(:, :, 2));
orientationAgreement(~isfinite(orientationAgreement)) = 0;
distanceTerm = exp(-sqrt(bestD2) ./ 125);
confidence = bestBladeConfidence .* distanceTerm .* ...
    (0.45 + 0.35 .* orientationAgreement + 0.20 .* coherence);
confidence(~valid) = 0;
end

function score = build_activity_score(eventMaps, variability, structural)
if isempty(eventMaps)
    score = 0.68 .* variability + 0.32 .* structural;
else
    cum = robust_normalize(eventMaps.cumDff, 0, 99.5);
    peak = robust_normalize(eventMaps.largestEventPeak, 0, 99.5);
    duration = robust_normalize(eventMaps.activeDurationSec, 0, 99.5);
    cum = match_image_size(cum, size(structural));
    peak = match_image_size(peak, size(structural));
    duration = match_image_size(duration, size(structural));
    score = 0.38 .* cum + 0.24 .* peak + 0.14 .* duration + ...
        0.14 .* variability + 0.10 .* structural;
end
score = imgaussfilt(score, 1, 'Padding', 'replicate');
score = robust_normalize(score, 0, 99.7);
end

function img = match_image_size(img, targetSize)
if ~isequal(size(img, 1), targetSize(1)) || ~isequal(size(img, 2), targetSize(2))
    img = imresize(img, targetSize(1:2), 'bilinear');
end
end

function [whole, distal] = build_proposals(score, structural, ridge, blades, radialMap, maxCount)
whole = empty_proposals();
distal = empty_proposals();
work = 0.70 .* score + 0.18 .* structural + 0.12 .* ridge;
positive = work(work > 0 & isfinite(work));
if isempty(positive)
    return
end
floorVal = prctile(positive, 93);
suppressed = work;
for i = 1:(maxCount * 4)
    [peak, idx] = max(suppressed(:));
    if ~isfinite(peak) || peak < floorVal
        break
    end
    [y, x] = ind2sub(size(work), idx);
    center = [x y];
    radial = squeeze(radialMap(y, x, :))';
    if any(~isfinite(radial)) || norm(radial) < 0.5
        radial = nearest_blade_radial(center, blades);
    end
    if any(~isfinite(radial)) || norm(radial) < 0.5
        radial = [0 -1];
    end
    radial = radial ./ norm(radial);

    [boundaryPoint, boundaryDistance] = nearest_boundary(center, blades);
    if isfinite(boundaryDistance) && boundaryDistance < 150
        wholeCenter = 0.44 .* boundaryPoint + 0.56 .* center;
        wholeCenter = wholeCenter + 18 .* radial;
    else
        wholeCenter = center;
    end
    wholeCandidate = make_proposal(wholeCenter, radial, 82, 54, peak, 'whole-cell');
    distalCandidate = make_proposal(center, radial, 58, 44, peak, 'distal');
    if ~proposal_overlaps(wholeCandidate, whole, 112)
        whole(end + 1) = wholeCandidate; %#ok<AGROW>
    end
    if ~proposal_overlaps(distalCandidate, distal, 88)
        distal(end + 1) = distalCandidate; %#ok<AGROW>
    end
    [xx, yy] = meshgrid(1:size(work, 2), 1:size(work, 1));
    suppressed((xx - x) .^ 2 + (yy - y) .^ 2 <= 68 ^ 2) = 0;
    if numel(whole) >= maxCount && numel(distal) >= maxCount
        break
    end
end
whole = whole(1:min(maxCount, numel(whole)));
distal = distal(1:min(maxCount, numel(distal)));
end

function proposals = empty_proposals()
proposals = struct('center', {}, 'radial', {}, 'poly', {}, 'score', {}, 'type', {});
end

function proposal = make_proposal(center, radial, halfLong, halfShort, score, type)
tangent = [-radial(2), radial(1)];
corners = [-halfLong -halfShort; halfLong -halfShort; ...
    halfLong halfShort; -halfLong halfShort];
poly = center + corners(:, 1) .* radial + corners(:, 2) .* tangent;
proposal = struct('center', center, 'radial', radial, 'poly', poly, ...
    'score', score, 'type', type);
end

function tf = proposal_overlaps(candidate, selected, minDistance)
tf = false;
for i = 1:numel(selected)
    if norm(candidate.center - selected(i).center) < minDistance
        tf = true;
        return
    end
    candidateShape = polyshape(candidate.poly(:, 1), candidate.poly(:, 2), ...
        'Simplify', false);
    selectedShape = polyshape(selected(i).poly(:, 1), selected(i).poly(:, 2), ...
        'Simplify', false);
    overlapArea = area(intersect(candidateShape, selectedShape));
    smallerArea = min(area(candidateShape), area(selectedShape));
    if overlapArea > 0.08 * smallerArea
        tf = true;
        return
    end
end
end

function radial = nearest_blade_radial(center, blades)
radial = [nan nan];
best = inf;
for b = 1:numel(blades)
    curve = blades(b).curve;
    [d2, idx] = min(sum((curve(:, 1:2) - center) .^ 2, 2));
    if d2 < best
        best = d2;
        radial = curve(idx, 3:4);
    end
end
end

function [point, distance] = nearest_boundary(center, blades)
point = [nan nan];
distance = inf;
for b = 1:numel(blades)
    curve = blades(b).curve;
    [d2, idx] = min(sum((curve(:, 1:2) - center) .^ 2, 2));
    if d2 < distance ^ 2
        distance = sqrt(d2);
        point = curve(idx, 1:2);
    end
end
end

function value = round_clip(value, limit)
value = min(max(round(value), 1), limit);
end

function opts = apply_defaults(opts, defaults)
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
end
