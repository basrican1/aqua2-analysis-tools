function out = addCalibratedCanonicalMetrics(canonicalOut, epochSpec)
% Add calibrated canonical summary metrics as local companions.
% Raw canonical fields remain unchanged.

out = canonicalOut;
if ~isstruct(canonicalOut) || ~isfield(canonicalOut, 'subset') || ~isstruct(canonicalOut.subset)
    return
end
if ~isfield(canonicalOut.subset, 'events') || ~istable(canonicalOut.subset.events) || ...
        ~isfield(canonicalOut.subset, 'files') || ~istable(canonicalOut.subset.files)
    return
end
if ~isfield(canonicalOut, 'metrics') || ~isstruct(canonicalOut.metrics)
    out.metrics = struct();
end
E = canonicalOut.subset.events;
F = canonicalOut.subset.files;
if isempty(E) || isempty(F)
    return
end
fileKeep = intersect({'file_id','frame_rate_s','spatial_res'}, F.Properties.VariableNames);
if ismember('file_id', fileKeep)
    E = outerjoin(E, unique(F(:, fileKeep), 'rows', 'stable'), 'Keys', {'file_id'}, 'Type', 'left', 'MergeKeys', true);
end
pairs = { ...
    'area_px2', 'area_um2', 'spatial_res', 'area'; ...
    'duration_frames', 'duration_sec', 'frame_rate_s', 'frame'; ...
    'rise19_frames', 'rise19_sec', 'frame_rate_s', 'frame'; ...
    'fall91_frames', 'fall91_sec', 'frame_rate_s', 'frame'; ...
    'decay_tau_frames', 'decay_tau_sec', 'frame_rate_s', 'frame'; ...
    'dff_auc', 'dff_auc_sec', 'frame_rate_s', 'frame'; ...
    'perimeter_px', 'perimeter_um', 'spatial_res', 'distance'; ...
    'max_prop_speed', 'max_prop_speed_um_per_sec', '', 'speed'; ...
    'propagation_distance_px', 'propagation_distance_um', 'spatial_res', 'distance'};
for i = 1:size(pairs,1)
    src = pairs{i,1};
    dst = pairs{i,2};
    calCol = pairs{i,3};
    mode = pairs{i,4};
    if ~ismember(src, E.Properties.VariableNames)
        continue
    end
    if ~strcmp(mode, 'speed') && ~ismember(calCol, E.Properties.VariableNames)
        continue
    end
    x = double(E.(src));
    if strcmp(mode, 'speed')
        if ~all(ismember({'spatial_res','frame_rate_s'}, E.Properties.VariableNames))
            continue
        end
        c = [double(E.spatial_res), double(E.frame_rate_s)];
        y = NaN(size(x));
        ok = isfinite(x) & isfinite(c(:,1)) & isfinite(c(:,2)) & c(:,2) > 0;
    else
        c = double(E.(calCol));
        y = NaN(size(x));
        ok = isfinite(x) & isfinite(c);
    end
    switch mode
        case 'frame'
            y(ok) = x(ok) .* c(ok);
        case 'distance'
            y(ok) = x(ok) .* c(ok);
        case 'area'
            y(ok) = x(ok) .* (c(ok).^2);
        case 'speed'
            y(ok) = x(ok) ./ c(ok,2);
    end
    E.(dst) = y;
    E2 = E;
    try
        mOut = computeMetricBySliceCanonical(E2, F, dst, epochSpec, 'groupBy', {'animal_id','slice_num','file_id'});
        mOut.titleName = dst;
        out.metrics.(dst) = mOut;
    catch
    end
end
out.subset.events = E;
end
