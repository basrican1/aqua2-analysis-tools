function outStruct = mergeStructs(toMergeStruct, treeloc, outStruct, newTreeLoc, dim)
%MERGESTRUCTS Merge data at treeloc into outStruct(newTreeLoc) along dim.
%   This local implementation is type-safe for numeric/cell/string payloads and
%   avoids padarray-based assumptions that fail on nonnumeric canonical fields.

    if nargin < 4 || isempty(newTreeLoc)
        newTreeLoc = treeloc;
    end
    if nargin < 5 || isempty(dim)
        dim = 2;
    end

    currentSize = FindFieldSize(outStruct, newTreeLoc);
    if isempty(currentSize)
        data = FindFieldValue(toMergeStruct, treeloc);
        outStruct = PutFieldValue(data, newTreeLoc, outStruct, 0);
        return
    end

    data1 = FindFieldValue(outStruct, newTreeLoc);
    data2 = FindFieldValue(toMergeStruct, treeloc);
    mergeData = mergeValueRecursive(data1, data2, dim);
    outStruct = PutFieldValue(mergeData, newTreeLoc, outStruct, 0);
end

function out = mergeValueRecursive(a, b, dim)
    if isempty(a)
        out = b;
        return
    end
    if isempty(b)
        out = a;
        return
    end

    if isstruct(a) || isstruct(b)
        out = mergeStructPayload(a, b, dim);
        return
    end

    if isnumeric(a) && isnumeric(b)
        out = concatNumericWithPad(a, b, dim);
        return
    end

    if islogical(a) && islogical(b)
        out = concatLogicalWithPad(a, b, dim);
        return
    end

    if iscell(a) && iscell(b)
        out = concatCellWithPad(a, b, dim);
        return
    end

    if isstring(a) && isstring(b)
        out = concatStringWithPad(a, b, dim);
        return
    end

    % Fallback for mixed/unsupported types.
    out = concatCellWithPad(wrapCell(a), wrapCell(b), dim);
end

function out = mergeStructPayload(a, b, dim)
    if ~isstruct(a)
        out = b;
        return
    end
    if ~isstruct(b)
        out = a;
        return
    end

    fa = fieldnames(a);
    fb = fieldnames(b);
    allf = unique([fa; fb], 'stable');
    out = struct();
    for ii = 1:numel(allf)
        f = allf{ii};
        hasA = isfield(a, f);
        hasB = isfield(b, f);
        if hasA && hasB
            out.(f) = mergeValueRecursive(a.(f), b.(f), dim);
        elseif hasA
            out.(f) = a.(f);
        else
            out.(f) = b.(f);
        end
    end
end

function out = concatNumericWithPad(a, b, dim)
    [a2, b2] = padPairToMatch(a, b, dim, @nanLike);
    out = cat(dim, a2, b2);
end

function out = concatLogicalWithPad(a, b, dim)
    [a2, b2] = padPairToMatch(a, b, dim, @falseLike);
    out = cat(dim, a2, b2);
end

function out = concatCellWithPad(a, b, dim)
    [a2, b2] = padPairToMatch(a, b, dim, @emptyCellLike);
    out = cat(dim, a2, b2);
end

function out = concatStringWithPad(a, b, dim)
    [a2, b2] = padPairToMatch(a, b, dim, @missingStringLike);
    out = cat(dim, a2, b2);
end

function [a2, b2] = padPairToMatch(a, b, dim, fillFactory)
    sa = size(a);
    sb = size(b);
    nDims = max([numel(sa), numel(sb), dim]);
    sa(end+1:nDims) = 1;
    sb(end+1:nDims) = 1;

    ta = sa;
    tb = sb;
    for kk = 1:nDims
        if kk ~= dim
            t = max(sa(kk), sb(kk));
            ta(kk) = t;
            tb(kk) = t;
        end
    end

    a2 = padToSize(a, ta, fillFactory);
    b2 = padToSize(b, tb, fillFactory);
end

function out = padToSize(in, targetSz, fillFactory)
    out = in;
    s = size(in);
    nDims = max(numel(s), numel(targetSz));
    s(end+1:nDims) = 1;
    targetSz(end+1:nDims) = 1;
    if isequal(s, targetSz)
        return
    end

    filler = fillFactory(targetSz, in);
    idx = repmat({':'}, 1, nDims);
    for kk = 1:nDims
        idx{kk} = 1:s(kk);
    end
    filler(idx{:}) = in;
    out = filler;
end

function v = nanLike(sz, in)
    if isa(in, 'single')
        v = nan(sz, 'single');
    else
        v = nan(sz);
    end
end

function v = falseLike(sz, ~)
    v = false(sz);
end

function v = emptyCellLike(sz, ~)
    v = cell(sz);
end

function v = missingStringLike(sz, ~)
    v = strings(sz);
    v(:) = missing;
end

function c = wrapCell(x)
    c = {x};
end
