function transform = bwaPreferredTransformForMetricID(metricID, varargin)
%BWAPREFERREDTRANSFORMFORMETRICID Recommend normalization transform per metric.
%
% Policy:
%   - Canonical positive magnitude metrics -> log2_fold
%   - Cell / event / recruitment / network metrics -> delta
%   - Bounded or signed metrics (fractions, circularity, motion, angles,
%     net-distance-change) -> delta

metricID = lower(char(string(metricID)));
unitLevel = "";
if nargin > 1 && ~isempty(varargin{1})
    unitLevel = lower(string(varargin{1}));
elseif startsWith(metricID, "cell.")
    unitLevel = "cell";
elseif startsWith(metricID, "slice.")
    unitLevel = "slice";
elseif startsWith(metricID, "event.")
    unitLevel = "event";
else
    unitLevel = "canonical";
end

% Non-canonical metrics default to delta.
if unitLevel ~= "canonical"
    transform = 'delta';
    return
end

% Canonical metrics that should not be folded.
deltaCanonical = [ ...
    "circularity", ...
    "propagation_angle_deg", ...
    "motion_index"];
if any(strcmp(metricID, deltaCanonical)) || contains(metricID, "angle")
    transform = 'delta';
    return
end

% Positive magnitude canonical metrics.
log2Canonical = [ ...
    "frequency", ...
    "area_px2","area_um2", ...
    "duration_frames","duration_sec", ...
    "rise19_frames","rise19_sec", ...
    "fall91_frames","fall91_sec", ...
    "decay_tau_frames","decay_tau_sec", ...
    "dff_auc","dff_auc_sec", ...
    "dff_max2", ...
    "perimeter_px","perimeter_um", ...
    "propagation_distance_px","propagation_distance_um", ...
    "max_prop_speed","max_prop_speed_um_per_sec"];
if any(strcmp(metricID, log2Canonical))
    transform = 'log2_fold';
    return
end

% Safe fallback.
transform = 'delta';
end
