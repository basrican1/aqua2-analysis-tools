function T = build_metric_qc_report(canonicalOut, derivedStandard)
% Build a compact per-metric QC report table from canonical and derived outputs.

rows = repmat(rowTemplate(), 0, 1);

% Canonical frequency
if isstruct(canonicalOut) && isfield(canonicalOut, 'frequency') && isstruct(canonicalOut.frequency)
    f = canonicalOut.frequency;
    r = rowTemplate();
    r.source = "canonical";
    r.metric = "events_per_min";
    r.n_input = getFieldOrNaN(f, 'n_events_input');
    r.n_excluded_non_epoch = getFieldOrNaN(f, 'n_events_excluded_non_epoch');
    r.n_excluded_invalid = NaN;
    r.n_negative = NaN;
    r.negative_dropped = false;
    r.notes = "";
    rows(end+1,1) = r; %#ok<AGROW>
end

% Canonical direct metrics
if isstruct(canonicalOut) && isfield(canonicalOut, 'metrics') && isstruct(canonicalOut.metrics)
    mNames = fieldnames(canonicalOut.metrics);
    for i = 1:numel(mNames)
        m = mNames{i};
        mo = canonicalOut.metrics.(m);
        r = rowTemplate();
        r.source = "canonical";
        r.metric = string(m);
        r.n_input = getFieldOrNaN(mo, 'n_events_input');
        r.n_excluded_non_epoch = getFieldOrNaN(mo, 'n_events_excluded_non_epoch');
        r.n_excluded_invalid = getFieldOrNaN(mo, 'n_events_excluded_invalid_metric');
        r.n_negative = NaN;
        r.negative_dropped = false;
        r.notes = "";
        rows(end+1,1) = r; %#ok<AGROW>
    end
end

% Derived BA-style metrics
if isstruct(derivedStandard)
    rows = addDerivedRows(rows, derivedStandard);
end

if isempty(rows)
    T = table();
else
    T = struct2table(rows);
end
end

function rows = addDerivedRows(rows, S)
if isfield(S, 'baviolin') && isstruct(S.baviolin)
    fn = fieldnames(S.baviolin);
    for i = 1:numel(fn)
        rows = addCompRow(rows, string(fn{i}), S.baviolin.(fn{i}));
    end
end
if isfield(S, 'centroid') && isstruct(S.centroid) && isfield(S.centroid, 'distComp')
    rows = addCompRow(rows, "centroid_dist", S.centroid.distComp);
end
if isfield(S, 'propagation') && isstruct(S.propagation)
    if isfield(S.propagation, 'distance')
        rows = addCompRow(rows, "propagation_distance", S.propagation.distance);
    end
    if isfield(S.propagation, 'angle')
        rows = addCompRow(rows, "propagation_angle", S.propagation.angle);
    end
end
if isfield(S, 'footprint')
    rows = addCompRow(rows, "footprint", S.footprint);
end
end

function rows = addCompRow(rows, metricName, comp)
if ~isstruct(comp)
    return
end
r = rowTemplate();
r.source = "derived";
r.metric = string(metricName);
r.n_input = getFieldOrNaN(comp, 'n_events_input');
r.n_excluded_non_epoch = getFieldOrNaN(comp, 'n_non_epoch');
r.n_excluded_invalid = getFieldOrNaN(comp, 'n_nonfinite_input');
r.n_negative = getFieldOrNaN(comp, 'n_negative_values');
r.negative_dropped = logical(getFieldOrDefault(comp, 'negative_values_dropped', false));
r.notes = "";
rows(end+1,1) = r; %#ok<AGROW>
end

function r = rowTemplate()
r = struct( ...
    'source', "", ...
    'metric', "", ...
    'n_input', NaN, ...
    'n_excluded_non_epoch', NaN, ...
    'n_excluded_invalid', NaN, ...
    'n_negative', NaN, ...
    'negative_dropped', false, ...
    'notes', "");
end

function v = getFieldOrNaN(S, f)
if isstruct(S) && isfield(S, f) && ~isempty(S.(f))
    v = double(S.(f));
else
    v = NaN;
end
end

function v = getFieldOrDefault(S, f, d)
if isstruct(S) && isfield(S, f) && ~isempty(S.(f))
    v = S.(f);
else
    v = d;
end
end

