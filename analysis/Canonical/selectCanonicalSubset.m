function S = selectCanonicalSubset(files, cells, events, event_cell_map, varargin)
% Return filtered canonical tables for downstream compute/plot steps.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'files', @istable);
addRequired(p, 'cells', @istable);
addRequired(p, 'events', @istable);
addRequired(p, 'event_cell_map', @istable);
addParameter(p, 'file_ids', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'animal_ids', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'slice_nums', {}, @(x) isnumeric(x) || iscell(x) || isstring(x));
addParameter(p, 'genotypes', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'treatments', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'qc_mask_files', [], @(x) isempty(x) || islogical(x));
parse(p, files, cells, events, event_cell_map, varargin{:});

F = files;
if ~isempty(p.Results.qc_mask_files)
    q = p.Results.qc_mask_files(:);
    if numel(q)==height(F)
        F = F(q,:);
    end
end

F = filterEq(F, 'file_id', p.Results.file_ids);
F = filterEq(F, 'animal_id', p.Results.animal_ids);
F = filterEq(F, 'slice_num', p.Results.slice_nums);
F = filterEq(F, 'genotype', p.Results.genotypes);
F = filterEq(F, 'treatment', p.Results.treatments);

keepFile = canonCellstr(F.file_id);

E = events(ismember(canonCellstr(events.file_id), keepFile), :);
C = cells(ismember(canonCellstr(cells.file_id), keepFile), :);

keepEvent = canonCellstr(E.event_uid);
keepCell = canonCellstr(C.cell_uid);
M = event_cell_map(ismember(canonCellstr(event_cell_map.event_uid), keepEvent) & ...
                   ismember(canonCellstr(event_cell_map.cell_uid), keepCell), :);

S = struct('files', F, 'cells', C, 'events', E, 'event_cell_map', M);
end

function T = filterEq(T, col, vals)
if isempty(vals) || ~ismember(col, T.Properties.VariableNames)
    return
end
lhs = canonCellstr(T.(col));
rhs = canonCellstr(vals);
T = T(ismember(lhs, rhs), :);
end

function c = canonCellstr(v)
if isempty(v)
    c = {};
elseif iscell(v)
    c = cellfun(@asCharVal, v, 'UniformOutput', false);
elseif isstring(v)
    c = cellstr(v(:));
elseif ischar(v)
    c = {v};
elseif isnumeric(v) || islogical(v)
    c = cellstr(string(v(:)));
else
    c = cellstr(string(v));
end
end

function c = asCharVal(v)
if ischar(v), c = v;
elseif isstring(v), c = char(v);
else, c = char(string(v));
end
end
