function family = bwaResolveMetricFamily(metricID, registry)
%BWARESOLVEMETRICFAMILY Resolve one or more metric IDs/labels to biological families.

metricID = string(metricID(:));
family = strings(size(metricID));

for i = 1:numel(metricID)
    row = table();
    if nargin >= 2 && ~isempty(registry)
        row = bwaLookupMetricVocabulary(metricID(i), registry);
    else
        row = bwaLookupMetricVocabulary(metricID(i));
    end
    if ~isempty(row) && ismember('family', row.Properties.VariableNames) && strlength(string(row.family(1))) > 0
        family(i) = string(row.family(1));
        continue
    end
    family(i) = resolveFamilyHeuristic(metricID(i));
end
end

function family = resolveFamilyHeuristic(metricID)
id = lower(char(string(metricID)));
id = strrep(id, '-', '_');
id = strrep(id, ' ', '_');

if contains(id, 'coincidence_rate') || contains(id, 'burst_fraction') || contains(id, 'pairwise') || contains(id, 'network')
    family = "synchrony";
elseif contains(id, 'freq') || contains(id, 'event_rate') || contains(id, 'frac_active_cells')
    family = "activity";
elseif contains(id, 'perimeter') || contains(id, 'circularity') || contains(id, 'area') || contains(id, 'footprint')
    family = "shape";
elseif contains(id, 'duration') || contains(id, 'rise19') || contains(id, 'fall91') || contains(id, 'decay_tau') || contains(id, 'decaytau')
    family = "kinetics";
elseif contains(id, 'dff_auc') || contains(id, 'dffauc') || contains(id, 'dff_max2') || contains(id, 'dffmax2') || contains(id, 'amplitude')
    family = "intensity";
elseif contains(id, 'propagation') || contains(id, 'propdist') || contains(id, 'propangle') || contains(id, 'max_prop_speed') || contains(id, 'maxpropspeed') || contains(id, 'centroid')
    family = "propagation";
elseif contains(id, 'outside_origin') || contains(id, 'outside_cell') || contains(id, 'within_origin') || contains(id, 'contained')
    family = "containment";
elseif contains(id, 'multicell') || contains(id, 'cells_per') || contains(id, 'neighbor') || contains(id, 'recruit')
    family = "recruitment";
elseif contains(id, 'soma') || contains(id, 'toward') || contains(id, 'away') || contains(id, 'traverse') || contains(id, 'motion_index')
    family = "soma";
else
    family = "";
end
end
