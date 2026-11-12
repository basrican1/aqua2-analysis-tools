function subs = numSubplotsForN(n)
% Return [nRows nCols] for n subplot panels.
% Uses existing numSubplots() if available, otherwise a local fallback.

    if nargin < 1 || isempty(n) || n <= 0
        subs = [1 1];
        return;
    end

    if exist('numSubplots', 'file') == 2
        subs = numSubplots(n);
        return;
    end

    nRows = floor(sqrt(n));
    if nRows < 1
        nRows = 1;
    end
    nCols = ceil(n / nRows);

    while (nRows - 1) * nCols >= n
        nRows = nRows - 1;
    end

    subs = [nRows nCols];
end
