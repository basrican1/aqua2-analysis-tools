function out = summarize_CorrelationStats_top(source, varargin)
%SUMMARIZE_CORRELATIONSTATS_TOP Return top positive and negative pairs.
%
% Usage:
%   out = summarize_CorrelationStats_top(corrOut)
%   out = summarize_CorrelationStats_top(corrOut, 'topK', 10)
%   out = summarize_CorrelationStats_top(result, 'unitLevel', 'slice', ...
%       'metricFamilies', {'shape','synchrony'}, 'epochName', 'ExperimentalEpoch1')
%
% Name-value:
%   'topK'        : number of rows to return for positive/negative sets (default 10)
%   'minAbsRho'   : minimum absolute rho to keep (default 0)
%   'maxQValue'   : maximum q-value to keep (default Inf)
%   'maxPValue'   : maximum p-value to keep (default Inf)
%   'minN'        : minimum N to keep in the summary table (default 3)
%   plus any CorrelationStats name-value args when source is not already a CorrelationStats output
%
% Output:
%   out.all       : filtered pair table
%   out.positive  : top positive rows by q, p, |rho|
%   out.negative  : top negative rows by q, p, |rho|
%   out.topK      : requested topK
%
% Notes:
%   - Uses corrOut.pairTable as the source of truth.
%   - Positive/negative subsets are sorted by q_value, then p_value, then abs_rho.

if isstruct(source) && isfield(source, 'source') && strcmp(source.source, 'CorrelationStats')
    corrOut = source;
    corrArgs = {};
else
    corrOut = [];
    corrArgs = varargin;
end

p = inputParser;
p.StructExpand = false;
addParameter(p, 'topK', 10, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'minAbsRho', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'maxQValue', inf, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'maxPValue', inf, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'minN', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 2);
p.KeepUnmatched = true;
parse(p, varargin{:});

if isempty(corrOut)
    corrOut = CorrelationStats(source, corrArgs{:});
end

if ~isfield(corrOut, 'pairTable') || ~istable(corrOut.pairTable)
    error('summarize_CorrelationStats_top:MissingPairTable', ...
        'CorrelationStats output does not contain a usable pairTable.');
end

T = corrOut.pairTable;
if isempty(T)
    out = struct('all', T, 'positive', T, 'negative', T, 'topK', double(p.Results.topK));
    return
end

countCol = '';
if ismember('n_obs', T.Properties.VariableNames)
    countCol = 'n_obs';
elseif ismember('n', T.Properties.VariableNames)
    countCol = 'n';
end
req = {'rho','abs_rho','p_value','q_value'};
if strlength(string(countCol)) == 0 || ~all(ismember(req, T.Properties.VariableNames))
    miss = req(~ismember(req, T.Properties.VariableNames));
    if strlength(string(countCol)) == 0
        miss = [miss {'n_obs'}];
    end
    error('summarize_CorrelationStats_top:BadPairTable', ...
        'pairTable is missing required columns: %s', strjoin(miss, ', '));
end

nv = double(T.(countCol));
keep = isfinite(T.rho) & isfinite(T.abs_rho) & isfinite(T.p_value) & isfinite(T.q_value) & isfinite(nv) & ...
    nv >= double(p.Results.minN) & T.abs_rho >= double(p.Results.minAbsRho) & ...
    T.q_value <= double(p.Results.maxQValue) & T.p_value <= double(p.Results.maxPValue);
T = T(keep, :);

Tpos = T(T.rho > 0, :);
Tneg = T(T.rho < 0, :);
if ~isempty(Tpos)
    Tpos = sortrows(Tpos, {'q_value','p_value','abs_rho'}, {'ascend','ascend','descend'});
    Tpos = Tpos(1:min(height(Tpos), double(p.Results.topK)), :);
end
if ~isempty(Tneg)
    Tneg = sortrows(Tneg, {'q_value','p_value','abs_rho'}, {'ascend','ascend','descend'});
    Tneg = Tneg(1:min(height(Tneg), double(p.Results.topK)), :);
end

out = struct();
out.all = T;
out.positive = Tpos;
out.negative = Tneg;
out.topK = double(p.Results.topK);
out.minAbsRho = double(p.Results.minAbsRho);
out.maxQValue = double(p.Results.maxQValue);
out.maxPValue = double(p.Results.maxPValue);
out.minN = double(p.Results.minN);
out.unitLevel = corrOut.unitLevel;
if isfield(corrOut, 'metricFamilies')
    out.metricFamilies = corrOut.metricFamilies;
else
    out.metricFamilies = {};
end
out.epochName = corrOut.epochName;
out.dataView = corrOut.dataView;
out.method = corrOut.method;
end
