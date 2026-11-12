function out = computeMetricBySliceCanonical(events, files, metricCol, epochSpec, varargin)
% Compute table-native summaries for one event metric by slice/file/animal/etc.
% This function does not generate figures.
%
% Inputs:
%   events: canonical events table
%   files: canonical files table
%   metricCol: event column name (ex: 'area_px2', 'duration_frames', 'dff_max2')
%   epochSpec: struct from frameTimeConversionsByEpoch()
%
% Options:
%   'groupBy' grouping columns from files table (default {'animal_id','slice_num','file_id'})
%   'metricFilter' function handle, default keeps finite values
%
% Outputs:
%   out.eventValues
%   out.summaryByGroupEpoch
%   out.epochTable

p = inputParser;
p.StructExpand = false;
addRequired(p, 'events', @istable);
addRequired(p, 'files', @istable);
addRequired(p, 'metricCol', @(x) ischar(x) || isstring(x));
addRequired(p, 'epochSpec', @isstruct);
addParameter(p, 'groupBy', {'animal_id','slice_num','file_id'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'metricFilter', [], @(f) isempty(f) || isa(f, 'function_handle'));
parse(p, events, files, metricCol, epochSpec, varargin{:});

metricColRequested = char(string(metricCol));
metricCol = metricColRequested;
if ~ismember(metricCol, events.Properties.VariableNames)
    error('Metric column not found in events: %s', metricCol);
end

grpCols = cellstr(string(p.Results.groupBy(:)));
grpCols = grpCols(:)';
A = assignEpochsToEvents(events, epochSpec, 'method', 'auto');
E = A.events;
out = struct();
out.n_events_input = height(E);
out.n_events_excluded_non_epoch = sum(E.epoch_idx <= 0);
E = E(E.epoch_idx > 0, :);
out.n_events_assigned_epoch = height(E);

v = E.(metricCol);
metricFilter = p.Results.metricFilter;
if isempty(metricFilter)
    metricFilter = defaultMetricFilter(metricCol);
end
keep = metricFilter(v);
keep = logical(keep(:));
if numel(keep) ~= height(E)
    error('metricFilter must return logical vector with one value per row.');
end
out.n_events_excluded_invalid_metric = sum(~keep);
E = E(keep, :);

needCols = unique([{'file_id'}, grpCols], 'stable');
needCols = needCols(ismember(needCols, files.Properties.VariableNames));
if ~ismember('file_id', needCols)
    needCols = [{'file_id'}, needCols];
end
eventCols = {'event_uid','file_id','epoch_idx','epoch_name', metricCol};
if ismember('primary_cell_uid', E.Properties.VariableNames)
    eventCols{end+1} = 'primary_cell_uid';
end
T = outerjoin(E(:, eventCols), files(:, needCols), ...
    'Keys', 'file_id', 'Type', 'left', 'MergeKeys', true);

for i = 1:numel(grpCols)
    c = grpCols{i};
    if ~ismember(c, T.Properties.VariableNames)
        T.(c) = repmat("", height(T), 1);
    end
end

sumCols = unique([grpCols, {'epoch_idx','epoch_name'}], 'stable');
[G, grp] = groupTableRows(T, sumCols);
vals = T.(metricCol);
if isempty(G)
    meanVals = zeros(0,1);
    medianVals = zeros(0,1);
    stdVals = zeros(0,1);
    nEvents = zeros(0,1);
else
    nGrp = max(G);
    meanVals = NaN(nGrp,1);
    medianVals = NaN(nGrp,1);
    stdVals = NaN(nGrp,1);
    nEvents = zeros(nGrp,1);
    for gg = 1:nGrp
        x = vals(G == gg);
        x = x(isfinite(x));
        nEvents(gg) = numel(x);
        if isempty(x)
            continue
        end
        meanVals(gg) = mean(x);
        medianVals(gg) = median(x);
        stdVals(gg) = std(x);
    end
end
S = grp;
S.(['mean_' metricColRequested]) = meanVals;
S.(['median_' metricColRequested]) = medianVals;
S.(['std_' metricColRequested]) = stdVals;
S.n_events = nEvents;
S.sem = stdVals ./ sqrt(max(nEvents, 1));

out.metricCol = metricColRequested;
out.metricColResolved = metricCol;
out.eventValues = T;
out.summaryByGroupEpoch = S;
out.epochTable = A.epochTable;
out.groupBy = grpCols;
out.options = p.Results;
end

function f = defaultMetricFilter(metricCol)
switch lower(char(string(metricCol)))
    case 'decay_tau_frames'
        % Physically valid decay constants must be positive and finite.
        f = @(v) isfinite(v) & v > 0;
    otherwise
        f = @(v) isfinite(v);
end
end

function [G, grpTbl] = groupTableRows(T, cols)
n = height(T);
if n == 0
    G = zeros(0,1);
    grpTbl = T(:, cols);
    return
end

key = strings(n,1);
for ii = 1:numel(cols)
    c = cols{ii};
    key = key + "|" + toKeyStrings(T.(c));
end
[~, ~, G] = unique(key, 'stable');
G = full(double(G(:)));

firstIdx = accumarray(G, (1:n)', [], @min);
grpTbl = T(firstIdx, cols);
end

function s = toKeyStrings(v)
if isstring(v)
    s = v(:);
elseif iscell(v)
    s = strings(numel(v),1);
    for i = 1:numel(v)
        vi = v{i};
        if ischar(vi)
            s(i) = string(vi);
        elseif isstring(vi)
            s(i) = vi;
        elseif isempty(vi)
            s(i) = "";
        else
            s(i) = string(vi);
        end
    end
elseif ischar(v)
    s = repmat(string(v), size(v,1), 1);
elseif isnumeric(v) || islogical(v)
    s = string(v(:));
elseif iscategorical(v)
    s = string(v(:));
else
    s = string(v);
    s = s(:);
end
end
