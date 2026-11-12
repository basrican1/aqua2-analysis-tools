function [epochSpec, filesOut, meta] = buildRecordingRelativeEpochSpec(files, recordingManifest, epochManifest, windowManifest, varargin)
%BUILDRECORDINGRELATIVEEPOCHSPEC Expand relative windows for each recording.
% Events are assigned by onset with half-open windows: start <= t0 < end.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'files', @istable);
addRequired(p, 'recordingManifest', @(x) istable(x) || ischar(x) || isstring(x));
addRequired(p, 'epochManifest', @(x) istable(x) || ischar(x) || isstring(x));
addRequired(p, 'windowManifest', @(x) istable(x) || ischar(x) || isstring(x));
addParameter(p, 'includeRoles', {'baseline','response'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'analysisSpecPath', '', @(x) ischar(x) || isstring(x));
parse(p, files, recordingManifest, epochManifest, windowManifest, varargin{:});

R = readManifest(recordingManifest);
E = readManifest(epochManifest);
W = readManifest(windowManifest);
requireColumns(files, {'file_id'}, 'canonical files');
requireColumns(R, {'recording_id','canonical_file_id','stim_pattern'}, 'recording manifest');
requireColumns(E, {'recording_id','epoch_id','start_sec','end_sec'}, 'epoch manifest');
requireColumns(W, {'stim_pattern','window_id','start_rel_sec','end_rel_sec','role','label'}, 'window manifest');

assertUnique(string(files.file_id), 'canonical files.file_id');
assertUnique(string(R.canonical_file_id), 'recording manifest.canonical_file_id');
[matched, rIndex] = ismember(string(files.file_id), string(R.canonical_file_id));
if ~all(matched)
    error('RecordingRelativeEpochs:MissingRecordingMetadata', ...
        'Recording metadata is missing for canonical file_id(s): %s', ...
        strjoin(string(files.file_id(~matched)), ', '));
end

filesOut = files;
rNames = string(R.Properties.VariableNames);
for ii = 1:numel(rNames)
    name = rNames(ii);
    if name == "canonical_file_id", continue; end
    incoming = R.(char(name))(rIndex);
    if ismember(name, string(filesOut.Properties.VariableNames))
        validateExistingColumn(filesOut.(char(name)), incoming, char(name));
    else
        filesOut.(char(name)) = incoming;
    end
end

includeRoles = lower(string(p.Results.includeRoles(:)));
W = W(ismember(lower(string(W.role)), includeRoles), :);
if isempty(W)
    error('RecordingRelativeEpochs:NoWindows', 'No windows remain after includeRoles filtering.');
end
if any(double(W.end_rel_sec) <= double(W.start_rel_sec))
    error('RecordingRelativeEpochs:BadWindow', 'Every relative window must have end_rel_sec > start_rel_sec.');
end

windowIDs = unique(string(W.window_id), 'stable');
nEpoch = numel(windowIDs);
epochTable = table((1:nEpoch)', windowIDs, strings(nEpoch,1), strings(nEpoch,1), nan(nEpoch,1), ...
    'VariableNames', {'epoch_idx','epoch_name','epoch_role','epoch_label','duration_min'});
for ee = 1:nEpoch
    rows = string(W.window_id) == windowIDs(ee);
    roles = unique(string(W.role(rows)));
    labels = unique(string(W.label(rows)));
    if numel(roles) ~= 1 || numel(labels) ~= 1
        error('RecordingRelativeEpochs:WindowSemantics', ...
            'window_id %s must have one role and one label across patterns.', windowIDs(ee));
    end
    epochTable.epoch_role(ee) = roles;
    epochTable.epoch_label(ee) = labels;
    durations = unique((double(W.end_rel_sec(rows)) - double(W.start_rel_sec(rows))) ./ 60);
    if isscalar(durations), epochTable.duration_min(ee) = durations; end
end

nRows = 0;
for rr = 1:height(R)
    nRows = nRows + nnz(string(W.stim_pattern) == string(R.stim_pattern(rr)) | lower(string(W.stim_pattern)) == "all");
end
A = table('Size',[nRows 11], ...
    'VariableTypes', {'string','string','string','double','string','string','string','double','double','double','double'}, ...
    'VariableNames', {'file_id','recording_id','stim_pattern','epoch_idx','epoch_name','epoch_role','epoch_label','start_sec','end_sec','duration_min','anchor_sec'});

row = 0;
for rr = 1:height(R)
    recordingID = string(R.recording_id(rr));
    stimPattern = string(R.stim_pattern(rr));
    stimRows = string(E.recording_id) == recordingID & lower(string(E.epoch_id)) == "stimulation";
    if nnz(stimRows) ~= 1
        error('RecordingRelativeEpochs:StimulusAnchor', ...
            'Expected exactly one stimulation row for recording %s; found %d.', recordingID, nnz(stimRows));
    end
    anchorSec = double(E.end_sec(stimRows));
    if ~isfinite(anchorSec)
        error('RecordingRelativeEpochs:StimulusAnchor', 'Non-finite stimulation end for recording %s.', recordingID);
    end
    fileHit = string(filesOut.file_id) == string(R.canonical_file_id(rr));
    if nnz(fileHit) ~= 1
        error('RecordingRelativeEpochs:FileMapping', 'Expected one canonical file for recording %s.', recordingID);
    end
    windows = W(string(W.stim_pattern) == stimPattern | lower(string(W.stim_pattern)) == "all", :);
    if isempty(windows)
        error('RecordingRelativeEpochs:PatternWindows', 'No windows are declared for stimulation pattern %s.', stimPattern);
    end
    for ww = 1:height(windows)
        row = row + 1;
        epochName = string(windows.window_id(ww));
        epochIdx = find(windowIDs == epochName, 1);
        startSec = anchorSec + double(windows.start_rel_sec(ww));
        endSec = anchorSec + double(windows.end_rel_sec(ww));
        A(row,:) = {string(filesOut.file_id(fileHit)), recordingID, stimPattern, epochIdx, epochName, ...
            string(windows.role(ww)), string(windows.label(ww)), startSec, endSec, (endSec-startSec)/60, anchorSec};
    end
end

validateNoOverlaps(A);
analysisSpecPath = char(string(p.Results.analysisSpecPath));
if ~isempty(strtrim(analysisSpecPath)), validateFactorLevels(filesOut, analysisSpecPath); end

epochSpec = struct();
epochSpec.mode = 'recording_relative_times';
epochSpec.assignmentTable = A;
epochSpec.epochTable = epochTable;
epochSpec.epochNames = cellstr(epochTable.epoch_name);
epochSpec.frameRange = [];
epochSpec.timeRange = [];
epochSpec.boundaryPolicy = 'event_onset_half_open';
meta = struct('mode',epochSpec.mode, 'nFiles',height(filesOut), 'nAssignments',height(A), ...
    'includeRoles',{cellstr(includeRoles)}, 'boundaryPolicy',epochSpec.boundaryPolicy, ...
    'recordingManifest',sourceLabel(recordingManifest), 'epochManifest',sourceLabel(epochManifest), ...
    'windowManifest',sourceLabel(windowManifest), 'analysisSpec',analysisSpecPath);
end

function T = readManifest(value)
if istable(value), T = value; else, T = readtable(char(string(value)), 'TextType','string', 'VariableNamingRule','preserve'); end
end
function requireColumns(T,names,label)
missing=string(names(~ismember(names,T.Properties.VariableNames))); if ~isempty(missing), error('RecordingRelativeEpochs:MissingColumns','%s is missing: %s',label,strjoin(missing,', ')); end
end
function assertUnique(values,label), if numel(unique(values))~=numel(values), error('RecordingRelativeEpochs:DuplicateID','%s contains duplicate IDs.',label); end, end
function validateExistingColumn(current,incoming,name)
a=string(current); b=string(incoming); conflict=a~="" & b~="" & a~=b; if any(conflict), error('RecordingRelativeEpochs:MetadataConflict','Canonical and recording metadata disagree in %s for %d row(s).',name,nnz(conflict)); end
end
function validateNoOverlaps(A)
files=unique(A.file_id,'stable'); for ii=1:numel(files), T=sortrows(A(A.file_id==files(ii),:),'start_sec'); if height(T)>1 && any(T.start_sec(2:end)<T.end_sec(1:end-1)), error('RecordingRelativeEpochs:OverlappingWindows','Overlapping windows for file_id %s.',files(ii)); end, end
end
function validateFactorLevels(files,pathIn)
spec=jsondecode(fileread(pathIn)); if ~isfield(spec,'factors'),return;end; defs=spec.factors;
for ii=1:numel(defs)
    if iscell(defs), def=defs{ii}; else, def=defs(ii); end
    id=string(def.id); if ~isfield(def,'levels')||~ismember(id,string(files.Properties.VariableNames)),continue;end
    declared=string(def.levels(:)); observed=unique(string(files.(char(id))),'stable'); missing=observed(~ismember(observed,declared));
    if ~isempty(missing),error('RecordingRelativeEpochs:FactorLevel','Observed %s level(s) absent from analysis specification: %s',id,strjoin(missing,', '));end
end
end
function label=sourceLabel(value), if istable(value),label='<table>';else,label=char(string(value));end, end
