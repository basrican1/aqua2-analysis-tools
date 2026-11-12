function out = assignEpochsToEvents(events, epochSpec, varargin)
% Assign epoch index/name to each event row in canonical events table.
%
% Inputs:
%   events: canonical events table (requires t0_frame and/or t0_sec)
%   epochSpec: struct from frameTimeConversionsByEpoch() or struct with fields:
%              - frameRange [nEpoch x 2] (optional)
%              - timeRange  [nEpoch x 2] in minutes (optional)
%              - epochNames {nEpoch x 1} (optional)
%
% Options:
%   'method': 'auto' (default), 'frame', or 'time'.
%
% Output struct:
%   out.events: input table + epoch_idx, epoch_name
%   out.epochTable: epoch metadata table

p = inputParser;
p.StructExpand = false;
addRequired(p, 'events', @istable);
addRequired(p, 'epochSpec', @isstruct);
addParameter(p, 'method', 'auto', @(x) ischar(x) || isstring(x));
parse(p, events, epochSpec, varargin{:});

method = lower(char(string(p.Results.method)));
if ~ismember(method, {'auto','frame','time'})
    error('method must be auto, frame, or time');
end

if isfield(epochSpec, 'assignmentTable') && istable(epochSpec.assignmentTable) && ~isempty(epochSpec.assignmentTable)
    out = assignRecordingRelative(events, epochSpec);
    return
end

if ~isfield(epochSpec, 'frameRange'), epochSpec.frameRange = []; end
if ~isfield(epochSpec, 'timeRange'), epochSpec.timeRange = []; end
if ~isfield(epochSpec, 'epochNames'), epochSpec.epochNames = {}; end

nFrame = size(epochSpec.frameRange, 1);
nTime = size(epochSpec.timeRange, 1);
nEpoch = max(nFrame, nTime);
if nEpoch == 0
    error('epochSpec must include frameRange or timeRange');
end

names = epochSpec.epochNames;
if isempty(names)
    names = arrayfun(@(k) sprintf('Epoch%d', k), 1:nEpoch, 'UniformOutput', false)';
else
    names = names(:);
    if numel(names) ~= nEpoch
        names = arrayfun(@(k) sprintf('Epoch%d', k), 1:nEpoch, 'UniformOutput', false)';
    end
end

useFrame = strcmp(method, 'frame') || (strcmp(method, 'auto') && nFrame > 0 && ismember('t0_frame', events.Properties.VariableNames));
useTime = strcmp(method, 'time') || (strcmp(method, 'auto') && ~useFrame);

if useFrame
    if isempty(epochSpec.frameRange) || ~ismember('t0_frame', events.Properties.VariableNames)
        error('Frame method requires epochSpec.frameRange and events.t0_frame');
    end
    x = events.t0_frame;
    ranges = epochSpec.frameRange;
elseif useTime
    if isempty(epochSpec.timeRange)
        error('Time method requires epochSpec.timeRange');
    end
    if ismember('t0_sec', events.Properties.VariableNames)
        x = events.t0_sec ./ 60; % convert to minutes
    elseif ismember('t0_frame', events.Properties.VariableNames)
        error('events.t0_sec missing for time method; provide frame method instead.');
    else
        error('events table missing t0_sec/t0_frame fields.');
    end
    ranges = epochSpec.timeRange;
else
    error('Could not determine epoch assignment method.');
end

x = forceNumericCol(x);
epochIdx = zeros(height(events), 1);
for e = 1:size(ranges,1)
    lo = ranges(e,1);
    hi = ranges(e,2);
    in = x >= lo & x <= hi;
    epochIdx(in) = e;
end

epochName = repmat({''}, height(events), 1);
inEpoch = epochIdx > 0;
epochName(inEpoch) = names(epochIdx(inEpoch));

eventsOut = events;
eventsOut.epoch_idx = epochIdx;
eventsOut.epoch_name = string(epochName);

epochTable = table((1:nEpoch)', string(names(:)), 'VariableNames', {'epoch_idx','epoch_name'});
if nFrame > 0
    epochTable.frame_start = epochSpec.frameRange(:,1);
    epochTable.frame_end = epochSpec.frameRange(:,2);
end
if nTime > 0
    epochTable.time_start_min = epochSpec.timeRange(:,1);
    epochTable.time_end_min = epochSpec.timeRange(:,2);
end

out = struct('events', eventsOut, 'epochTable', epochTable, 'method', method);
end

function out = assignRecordingRelative(events, epochSpec)
A = epochSpec.assignmentTable;
needEvents = {'file_id','t0_sec'};
needAssignments = {'file_id','epoch_idx','epoch_name','start_sec','end_sec','duration_min'};
if ~all(ismember(needEvents, events.Properties.VariableNames))
    error('assignEpochsToEvents:RecordingRelativeEvents', ...
        'Recording-relative assignment requires events.file_id and events.t0_sec.');
end
if ~all(ismember(needAssignments, A.Properties.VariableNames))
    error('assignEpochsToEvents:RecordingRelativeSpec', 'assignmentTable is missing required columns.');
end

eventsOut = events;
eventsOut.epoch_idx = zeros(height(events),1);
eventsOut.epoch_name = strings(height(events),1);
eventsOut.epoch_role = strings(height(events),1);
eventsOut.epoch_duration_min = nan(height(events),1);
t0 = forceNumericCol(events.t0_sec);
for aa = 1:height(A)
    rows = string(events.file_id) == string(A.file_id(aa)) & ...
        t0 >= double(A.start_sec(aa)) & t0 < double(A.end_sec(aa));
    if any(eventsOut.epoch_idx(rows) ~= 0)
        error('assignEpochsToEvents:OverlappingAssignments', ...
            'An event matched multiple recording-relative windows.');
    end
    eventsOut.epoch_idx(rows) = double(A.epoch_idx(aa));
    eventsOut.epoch_name(rows) = string(A.epoch_name(aa));
    if ismember('epoch_role', A.Properties.VariableNames)
        eventsOut.epoch_role(rows) = string(A.epoch_role(aa));
    end
    eventsOut.epoch_duration_min(rows) = double(A.duration_min(aa));
end

if isfield(epochSpec, 'epochTable') && istable(epochSpec.epochTable)
    epochTable = epochSpec.epochTable;
else
    epochTable = unique(A(:, {'epoch_idx','epoch_name'}), 'rows', 'stable');
end
out = struct('events', eventsOut, 'epochTable', epochTable, ...
    'assignmentTable', A, 'method', 'recording_relative');
end

function x = forceNumericCol(x)
if iscell(x)
    out = NaN(numel(x),1);
    for i = 1:numel(x)
        xi = x{i};
        if isnumeric(xi) || islogical(xi)
            out(i) = double(xi);
        else
            out(i) = str2double(string(xi));
        end
    end
    x = out;
elseif isstring(x)
    x = str2double(x);
elseif ischar(x)
    x = str2double(x);
else
    x = double(x);
end
x = x(:);
end
