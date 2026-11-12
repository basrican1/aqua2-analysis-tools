function cache = canonicalEventFrameCache(xSpaTempEvent, t0Frame, t1Frame, spatialSize)
%CANONICALEVENTFRAMECACHE Build reusable per-event framewise spatial pixels.

cache = struct('frames', [], 'framePix', {{}});
if isempty(xSpaTempEvent) || ~isnumeric(xSpaTempEvent) || ~isfinite(t0Frame) || ~isfinite(t1Frame)
    return
end

spatialSize = double(spatialSize(:)');
if isempty(spatialSize)
    spatialSize = [NaN NaN NaN NaN];
end
if numel(spatialSize) < 4
    spatialSize(end+1:4) = NaN;
end
sx = round(spatialSize(1));
sy = round(spatialSize(2));
sz = round(spatialSize(3));
st = round(spatialSize(4));
if any(~isfinite([sx sy sz st])) || any([sx sy sz st] <= 0)
    return
end

frames = round(t0Frame):round(t1Frame);
if isempty(frames)
    return
end

vox = double(xSpaTempEvent(:));
try
    [ih, iw, il, it] = ind2sub([sx sy sz st], vox);
catch
    return
end

cache.frames = frames;
cache.framePix = cell(numel(frames), 1);
for kk = 1:numel(frames)
    ff = frames(kk);
    keep = (it == ff);
    if ~any(keep)
        cache.framePix{kk} = [];
        continue
    end
    pix = unique(sub2ind([sx sy sz], ih(keep), iw(keep), il(keep)));
    cache.framePix{kk} = pix(pix > 0 & isfinite(pix));
end
end
