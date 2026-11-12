function out = dominantPremergeBranchPrototype(source, eventID, varargin)
%DOMINANTPREMERGEBRANCHPROTOTYPE Inspect pre-merge source branches for one event.
%   OUT = dominantPremergeBranchPrototype(RESULT, EVENTID)
%   OUT = dominantPremergeBranchPrototype(PATHIN, EVENTID)
%   OUT = dominantPremergeBranchPrototype(BWA_AQ2, EVENTID)
%
%   This isolates branch-origin selection from wavefront path tracking.

BWA_setupPaths(fileparts(mfilename('fullpath')));

p = inputParser;
p.StructExpand = false;
addParameter(p, 'showFigure', true, @islogical);
addParameter(p, 'figureName', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'minComponentPixels', 10, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'maxComponentsPerFrame', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'minAreaFractionOfLargest', 0.2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
addParameter(p, 'mergePersistenceFrames', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
parse(p, varargin{:});

[BWA_AQ2, sourceLabel] = localResolveSource(source);
fts = BWA_AQ2.fts1;
[imageSz, spatialRes] = localEventImagingMeta(BWA_AQ2, eventID);
geomAll = CalcEventFrameGeometry(fts.loc.xSpaTemp, fts.loc.t0, fts.loc.t1, imageSz, 'eventIDs', eventID);
G = geomAll{eventID};
W = CalcWavefrontMetrics(G, spatialRes);

frames = G.frames(:);
branchOriginIdx = find(frames == W.branchOriginFrame, 1, 'first');
branchPeakIdx = find(frames == W.branchPeakFrame, 1, 'first');

nCompRaw = zeros(numel(frames), 1);
nComp = zeros(numel(frames), 1);
for k = 1:numel(frames)
    compsRaw = localFrameComponents(logical(G.frameMasks{k}), 1);
    compsKeep = localLimitComponents( ...
        localFrameComponents(logical(G.frameMasks{k}), p.Results.minComponentPixels), ...
        p.Results.maxComponentsPerFrame, p.Results.minAreaFractionOfLargest);
    nCompRaw(k) = numel(compsRaw);
    nComp(k) = numel(compsKeep);
end

mergeFrameIdx = NaN;
for k = max(branchOriginIdx + 1, 2):branchPeakIdx
    if nComp(k - 1) > 1 && localHasPersistentMerge(nComp, k, p.Results.mergePersistenceFrames)
        mergeFrameIdx = k;
        break
    end
end

premergeIdx = NaN;
premergeComps = struct('mask', {}, 'pts', {}, 'centroid', {}, 'rep', {}, 'area', {});
dominantIdx = NaN;
localIdx = NaN;
if isfinite(mergeFrameIdx) && mergeFrameIdx > 1
    premergeIdx = mergeFrameIdx - 1;
    premergeComps = localLimitComponents( ...
        localFrameComponents(logical(G.frameMasks{premergeIdx}), p.Results.minComponentPixels), ...
        p.Results.maxComponentsPerFrame, p.Results.minAreaFractionOfLargest);
    if ~isempty(premergeComps)
        [~, dominantIdx] = max([premergeComps.area]);
        localIdx = localComponentIndexForPoint(premergeComps, W.branchOriginXY);
    end
end

statusText = "no premerge split detected";
if isfinite(mergeFrameIdx) && isempty(premergeComps)
    statusText = "merge detected, but no premerge components";
elseif isfinite(mergeFrameIdx) && ~isempty(premergeComps) && numel(premergeComps) == 1
    statusText = "merge frame found, but only one premerge component";
elseif isfinite(mergeFrameIdx) && numel(premergeComps) >= 2
    statusText = sprintf('premerge components=%d | dominant=%s | local=%s', ...
        numel(premergeComps), localIdxLabel(dominantIdx), localIdxLabel(localIdx));
end

hFig = [];
hAx = [];
hAxComp = [];
if p.Results.showFigure
    figName = char(string(p.Results.figureName));
    if isempty(strtrim(figName))
        figName = sprintf('DominantPremergeBranch_%d', eventID);
    end
    hFig = figure('Name', figName, 'Color', 'w');
    hAx = axes('Parent', hFig, 'Position', [0.08 0.12 0.58 0.78]);
    hold(hAx, 'on');
    axis(hAx, 'equal');
    set(hAx, 'YDir', 'reverse');
    xlabel(hAx, 'x (um)');
    ylabel(hAx, 'y (um)');
    hAxComp = axes('Parent', hFig, 'Position', [0.72 0.18 0.22 0.68]);
    hold(hAxComp, 'on');
    xlabel(hAxComp, 'frame');
    ylabel(hAxComp, 'n comp');
    title(hAxComp, sprintf('Components (min %d px, top %d, >= %.0f%% largest)', ...
        round(p.Results.minComponentPixels), round(p.Results.maxComponentsPerFrame), ...
        100 * p.Results.minAreaFractionOfLargest));
    plot(hAxComp, frames, nCompRaw, '-o', 'Color', [0.7 0.7 0.7], ...
        'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 3);
    plot(hAxComp, frames, nComp, '-ok', 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerSize', 4);
    xline(hAxComp, W.branchOriginFrame, '--y', 'Origin');
    xline(hAxComp, W.branchPeakFrame, '--r', 'Peak');
    if isfinite(localFrameValue(frames, premergeIdx))
        xline(hAxComp, localFrameValue(frames, premergeIdx), ':c', 'Premerge');
    end
    if isfinite(localFrameValue(frames, mergeFrameIdx))
        xline(hAxComp, localFrameValue(frames, mergeFrameIdx), ':m', 'Merge');
    end
    grid(hAxComp, 'on');

    bg = false(imageSz(1), imageSz(2));
    for k = 1:numel(G.frameMasks)
        bg = bg | logical(G.frameMasks{k});
    end
    [r, c] = find(bg);
    if ~isempty(r)
        plot(hAx, c * spatialRes, r * spatialRes, '.', 'Color', [0.82 0.82 0.82], 'MarkerSize', 6);
    end

    title(hAx, sprintf('Event %d | origin=%d | premerge=%s | peak=%d', ...
        eventID, round(W.branchOriginFrame), localFrameLabel(frames, premergeIdx), round(W.branchPeakFrame)));

    if isfinite(premergeIdx) && ~isempty(premergeComps)
        colors = lines(max(numel(premergeComps), 2));
        for i = 1:numel(premergeComps)
            pts = premergeComps(i).pts * spatialRes;
            plot(hAx, pts(:,1), pts(:,2), '.', 'Color', colors(i,:), 'MarkerSize', 8);
            rep = premergeComps(i).rep * spatialRes;
            label = sprintf('%d | area=%d', i, premergeComps(i).area);
            text(hAx, rep(1), rep(2), label, 'Color', colors(i,:), 'FontSize', 8, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        end
    end

    plot(hAx, W.branchOriginXY(1) * spatialRes, W.branchOriginXY(2) * spatialRes, 'ks', ...
        'MarkerFaceColor', 'y', 'MarkerSize', 8, 'LineWidth', 1.25);
    plot(hAx, W.branchPeakPointXY(1) * spatialRes, W.branchPeakPointXY(2) * spatialRes, 'rp', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 11, 'LineWidth', 1.25);

    if isfinite(dominantIdx)
        rep = premergeComps(dominantIdx).rep * spatialRes;
        plot(hAx, rep(1), rep(2), 'co', 'MarkerSize', 14, 'LineWidth', 2);
    end
    if isfinite(localIdx)
        rep = premergeComps(localIdx).rep * spatialRes;
        plot(hAx, rep(1), rep(2), 'mo', 'MarkerSize', 14, 'LineWidth', 2);
    end
    xl = xlim(hAx);
    yl = ylim(hAx);
    text(hAx, xl(1) + 0.02 * diff(xl), yl(1) + 0.04 * diff(yl), char(statusText), ...
        'Color', [1 1 1], 'FontSize', 9, 'FontWeight', 'bold', ...
        'BackgroundColor', [0 0 0], 'Margin', 4, 'VerticalAlignment', 'top');
end

out = struct();
out.figure = hFig;
out.axes = hAx;
out.componentAxes = hAxComp;
out.sourceLabel = sourceLabel;
out.geometry = G;
out.wavefront = W;
out.frames = frames;
out.nComponentsByFrameRaw = nCompRaw;
out.nComponentsByFrame = nComp;
out.minComponentPixels = p.Results.minComponentPixels;
out.maxComponentsPerFrame = p.Results.maxComponentsPerFrame;
out.minAreaFractionOfLargest = p.Results.minAreaFractionOfLargest;
out.mergePersistenceFrames = p.Results.mergePersistenceFrames;
out.mergeFrameIdx = mergeFrameIdx;
out.mergeFrame = localFrameValue(frames, mergeFrameIdx);
out.premergeFrameIdx = premergeIdx;
out.premergeFrame = localFrameValue(frames, premergeIdx);
out.lastSplitFrameIdx = localLastSplitFrameIdx(nComp, branchOriginIdx, branchPeakIdx);
out.lastSplitFrame = localFrameValue(frames, out.lastSplitFrameIdx);
out.premergeComponents = premergeComps;
out.premergeComponentAreas = reshape([premergeComps.area], 1, []);
out.premergeComponentReps = localComponentRepArray(premergeComps);
out.dominantComponentIdx = dominantIdx;
out.localOriginComponentIdx = localIdx;
out.statusText = statusText;
end

function txt = localFrameLabel(frames, idx)
if isempty(idx) || ~isfinite(idx)
    txt = 'none';
else
    txt = sprintf('%d', round(frames(idx)));
end
end

function val = localFrameValue(frames, idx)
val = NaN;
if ~isempty(idx) && isfinite(idx)
    val = frames(idx);
end
end

function idx = localLastSplitFrameIdx(nComp, startIdx, endIdx)
idx = NaN;
for k = min(endIdx, numel(nComp)):-1:max(startIdx, 1)
    if nComp(k) > 1
        idx = k;
        return
    end
end
end

function tf = localHasPersistentMerge(nComp, startIdx, persistenceFrames)
stopIdx = min(numel(nComp), startIdx + persistenceFrames - 1);
tf = startIdx <= numel(nComp) && all(nComp(startIdx:stopIdx) == 1);
end

function txt = localIdxLabel(idx)
if isempty(idx) || ~isfinite(idx)
    txt = "none";
else
    txt = string(idx);
end
end

function comps = localFrameComponents(mask, minPixels)
comps = struct('mask', {}, 'pts', {}, 'centroid', {}, 'rep', {}, 'area', {});
if isempty(mask) || ~any(mask(:))
    return
end
cc = bwconncomp(mask, 8);
for i = 1:cc.NumObjects
    if numel(cc.PixelIdxList{i}) < minPixels
        continue
    end
    compMask = false(size(mask));
    compMask(cc.PixelIdxList{i}) = true;
    [r, c] = find(compMask);
    ctr = [mean(c, 'omitnan'), mean(r, 'omitnan')];
    rep = localNearestPointToTarget([c, r], ctr);
    comps(end+1,1) = struct( ... %#ok<AGROW>
        'mask', compMask, ...
        'pts', [c, r], ...
        'centroid', ctr, ...
        'rep', rep, ...
        'area', numel(cc.PixelIdxList{i}));
end
end

function comps = localLimitComponents(comps, maxCount, minAreaFractionOfLargest)
if isempty(comps)
    return
end
[~, order] = sort([comps.area], 'descend');
comps = comps(order);
largestArea = comps(1).area;
minArea = max(1, minAreaFractionOfLargest * largestArea);
keep = [comps.area] >= minArea;
comps = comps(keep);
if isempty(comps)
    return
end
comps = comps(1:min(numel(comps), maxCount));
end

function reps = localComponentRepArray(comps)
if isempty(comps)
    reps = zeros(0, 2);
    return
end
reps = vertcat(comps.rep);
end

function idx = localComponentIndexForPoint(comps, pt)
idx = NaN;
best = inf;
for i = 1:numel(comps)
    pix = round(comps(i).pts);
    if any(pix(:,1) == round(pt(1)) & pix(:,2) == round(pt(2)))
        idx = i;
        return
    end
    d = min(hypot(comps(i).pts(:,1) - pt(1), comps(i).pts(:,2) - pt(2)));
    if d < best
        best = d;
        idx = i;
    end
end
end

function pt = localNearestPointToTarget(pts, target)
pt = [NaN NaN];
if isempty(pts)
    return
end
[~, ibest] = min(hypot(pts(:,1) - target(1), pts(:,2) - target(2)));
pt = pts(ibest, :);
end

function [BWA_AQ2, label] = localResolveSource(source)
label = "";
if isstruct(source) && isfield(source, 'derived') && isstruct(source.derived) && ...
        isfield(source.derived, 'context') && isstruct(source.derived.context) && ...
        isfield(source.derived.context, 'BWA_AQ2')
    cand = source.derived.context.BWA_AQ2;
    if isstruct(cand) && isfield(cand, 'fts1') && isfield(cand, 'opts')
        BWA_AQ2 = cand;
        if isfield(source, 'pathIn')
            label = string(source.pathIn);
        end
        return
    end
end
if isstruct(source) && isfield(source, 'pathIn') && strlength(string(source.pathIn)) > 0
    loaded = loadClusterOutputs(char(string(source.pathIn)), 'attachLegacyMeta', true, 'attachEventMeta', true);
    BWA_AQ2 = loaded.BWA_AQ2;
    label = string(source.pathIn);
    return
end
if ischar(source) || isstring(source)
    loaded = loadClusterOutputs(char(string(source)), 'attachLegacyMeta', true, 'attachEventMeta', true);
    BWA_AQ2 = loaded.BWA_AQ2;
    label = string(source);
    return
end
if isstruct(source) && isfield(source, 'fts1') && isfield(source, 'opts')
    BWA_AQ2 = source;
    return
end
error('dominantPremergeBranchPrototype:BadSource', 'Unsupported source type.');
end

function [imageSz, spatialRes] = localEventImagingMeta(BWA_AQ2, eventID)
opts = BWA_AQ2.opts;
imageSz = double(opts.sz(1,:));
spatialRes = double(opts.spatialRes(1));
fts = BWA_AQ2.fts1;
if isfield(fts, 'basic') && isfield(fts.basic, 'fileID') && ~isempty(fts.basic.fileID)
    fileIDs = string(fts.basic.fileID);
elseif isfield(fts, 'fileID') && ~isempty(fts.fileID)
    fileIDs = string(fts.fileID);
else
    fileIDs = strings(numel(fts.loc.t0), 1);
end
rowIDs = unique(fileIDs, 'stable');
if eventID <= numel(fileIDs)
    row = find(rowIDs == fileIDs(eventID), 1, 'first');
    if ~isempty(row)
        if size(opts.sz, 1) >= row
            imageSz = double(opts.sz(row,:));
        end
        if numel(opts.spatialRes) >= row && isfinite(opts.spatialRes(row))
            spatialRes = double(opts.spatialRes(row));
        end
    end
end
if numel(imageSz) < 3
    imageSz(3) = 1;
end
end
