function cellToSoma = resolveSingleSomaPerCellStandalone(inclu, cellMasks, lmkMasks, cellCenters, lmkCenters)
%RESOLVESINGLESOMAPERCELLSTANDALONE Enforce a single soma landmark per cell.
% Selection rule for cells with multiple candidate landmarks:
%   1) maximize mask overlap with the cell mask
%   2) minimize center-to-center distance to the cell center
%   3) smallest landmark index as deterministic tie-break

    if ~isnumeric(inclu) || isempty(inclu)
        cellToSoma = NaN(0,1);
        return
    end

    nCells = size(inclu,1);
    cellToSoma = NaN(nCells,1);
    for cc = 1:nCells
        cand = find(inclu(cc,:) > 0);
        if isempty(cand)
            continue
        end
        if numel(cand) == 1
            cellToSoma(cc) = cand(1);
            continue
        end
        cellToSoma(cc) = pickMostCentralLandmark(cand, cc, cellMasks, lmkMasks, cellCenters, lmkCenters);
    end
end

function best = pickMostCentralLandmark(cand, cellIdx, cellMasks, lmkMasks, cellCenters, lmkCenters)
    best = cand(1);
    scoreOverlap = -inf(numel(cand),1);
    scoreDist = inf(numel(cand),1);

    cMask = {};
    if iscell(cellMasks) && numel(cellMasks) >= cellIdx
        cMask = cellMasks(cellIdx);
    end

    cCenter = getCenterRow(cellCenters, cellIdx);
    if any(isnan(cCenter)) && ~isempty(cMask)
        cCenter = maskCenter(cMask{1});
    end

    for ii = 1:numel(cand)
        ll = cand(ii);
        lm = [];
        if iscell(lmkMasks) && numel(lmkMasks) >= ll
            lm = lmkMasks{ll};
        end

        if ~isempty(cMask) && ~isempty(lm) && (isnumeric(lm) || islogical(lm))
            cm = cMask{1};
            if isnumeric(cm) || islogical(cm)
                scoreOverlap(ii) = overlapCountSafe(cm, lm);
            end
        end

        lCenter = getCenterRow(lmkCenters, ll);
        if any(isnan(lCenter)) && ~isempty(lm)
            lCenter = maskCenter(lm);
        end
        if all(isfinite(cCenter(1:2))) && all(isfinite(lCenter(1:2)))
            scoreDist(ii) = hypot(cCenter(1)-lCenter(1), cCenter(2)-lCenter(2));
        end
    end

    [~, ord] = sortrows([ -scoreOverlap, scoreDist, cand(:) ], [1 2 3]);
    best = cand(ord(1));
end

function n = overlapCountSafe(a, b)
    n = -inf;
    if ~(isnumeric(a) || islogical(a)) || ~(isnumeric(b) || islogical(b))
        return
    end
    if isequal(size(a), size(b))
        n = nnz((a > 0) & (b > 0));
        return
    end
    ai = find(a > 0);
    bi = find(b > 0);
    if isempty(ai) || isempty(bi)
        n = 0;
        return
    end
    n = numel(intersect(ai(:), bi(:)));
end

function c = getCenterRow(C, idx)
    c = [NaN NaN NaN];
    if isnumeric(C) && size(C,1) >= idx
        row = C(idx,:);
        c(1:min(3,numel(row))) = row(1:min(3,numel(row)));
    end
end

function c = maskCenter(m)
    c = [NaN NaN NaN];
    if isempty(m) || ~(isnumeric(m) || islogical(m))
        return
    end
    nd = ndims(m);
    if nd <= 2
        [r, col] = find(m > 0);
        if ~isempty(r)
            c = [mean(r), mean(col), 1];
        end
    else
        [r, col, z] = ind2sub(size(m), find(m > 0));
        if ~isempty(r)
            c = [mean(r), mean(col), mean(z)];
        end
    end
end
