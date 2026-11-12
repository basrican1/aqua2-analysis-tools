function out = canonicalEventCellPartitionMetrics(xSpaTempOrCache, t0Frame, t1Frame, spatialSize, cellMaskOrPix, coverageThreshold)
%CANONICALEVENTCELLPARTITIONMETRICS Compute cell-resolved occupancy metrics.
% These metrics treat an event-cell row as a partition of one parent event,
% not as an independent event. Timing is reported in absolute frame numbers.

out = struct( ...
    'cell_first_overlap_frame', NaN, ...
    'cell_entry_frame_cov_thresh', NaN, ...
    'cell_last_overlap_frame', NaN, ...
    'cell_overlap_duration_frames', NaN, ...
    'cell_peak_overlap_npix', NaN, ...
    'cell_peak_overlap_frac_of_cell', NaN, ...
    'cell_peak_overlap_frac_of_event_frame', NaN);

if nargin < 6 || ~isfinite(coverageThreshold)
    coverageThreshold = 0.10;
end
coverageThreshold = max(0, double(coverageThreshold));

cellPix = forcePixOrMask(cellMaskOrPix);
if isempty(cellPix) || ~isfinite(t0Frame) || ~isfinite(t1Frame)
    return
end

cache = prepareFrameCache(xSpaTempOrCache, t0Frame, t1Frame, spatialSize);
if isempty(cache.frames)
    return
end

firstF = NaN;
entryThreshF = NaN;
lastF = NaN;
peakN = 0;
peakFracCell = NaN;
peakFracFrame = NaN;
cellN = numel(cellPix);

for kk = 1:numel(cache.frames)
    ff = cache.frames(kk);
    framePix = cache.framePix{kk};
    if isempty(framePix)
        continue
    end
    ov = numel(intersect(framePix, cellPix));
    if ov <= 0
        continue
    end
    if ~isfinite(firstF)
        firstF = ff;
    end
    fracCell = ov / max(cellN, 1);
    if ~isfinite(entryThreshF) && fracCell >= coverageThreshold
        entryThreshF = ff;
    end
    lastF = ff;
    if ov > peakN
        peakN = ov;
        peakFracCell = fracCell;
        peakFracFrame = ov / max(numel(framePix), 1);
    end
end

out.cell_first_overlap_frame = firstF;
out.cell_entry_frame_cov_thresh = entryThreshF;
out.cell_last_overlap_frame = lastF;
if isfinite(firstF) && isfinite(lastF)
    out.cell_overlap_duration_frames = lastF - firstF + 1;
end
if peakN > 0
    out.cell_peak_overlap_npix = peakN;
    out.cell_peak_overlap_frac_of_cell = peakFracCell;
    out.cell_peak_overlap_frac_of_event_frame = peakFracFrame;
end
end

function cellPix = forcePixOrMask(cellMaskOrPix)
cellPix = [];
if isnumeric(cellMaskOrPix) || islogical(cellMaskOrPix)
    if islogical(cellMaskOrPix) || ~isvector(cellMaskOrPix)
        cellPix = find(cellMaskOrPix > 0);
    else
        cellPix = double(cellMaskOrPix(:));
    end
end
cellPix = unique(cellPix(:));
cellPix = cellPix(isfinite(cellPix) & cellPix > 0);
end

function cache = prepareFrameCache(xSpaTempOrCache, t0Frame, t1Frame, spatialSize)
cache = struct('frames', [], 'framePix', {{}});
if isstruct(xSpaTempOrCache) && isfield(xSpaTempOrCache, 'frames') && isfield(xSpaTempOrCache, 'framePix')
    cache.frames = double(xSpaTempOrCache.frames(:)');
    cache.framePix = xSpaTempOrCache.framePix;
    return
end
cache = canonicalEventFrameCache(xSpaTempOrCache, t0Frame, t1Frame, spatialSize);
end
