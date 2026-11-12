function key = bwaMetricConceptKey(metricID)
%BWAMETRICCONCEPTKEY Resolve a metric or legacy panel key to a concept key for duplicate suppression.

id = lower(char(string(metricID)));
id = strrep(id, '-', '_');

switch id
    case {'area_um2','area_px2','event.area_um2','basic_area','cell.mean_event_area_um2'}
        key = "area";
    case {'duration_sec','duration_frames','event.duration_sec','curve_duration','cell.mean_event_duration_sec'}
        key = "duration";
    case {'rise19_sec','rise19_frames','event.rise19_sec','curve_rise19','cell.mean_event_rise19_sec'}
        key = "rise19";
    case {'fall91_sec','fall91_frames','event.fall91_sec','curve_fall91','cell.mean_event_fall91_sec'}
        key = "fall91";
    case {'decay_tau_sec','decay_tau_frames','event.decay_tau_sec','curve_decaytau','cell.mean_event_decay_tau_sec'}
        key = "decay_tau";
    case {'dff_auc_sec','dff_auc','event.signal_auc_dff_sec','curve_dffauc','cell.mean_event_signal_auc_dff_sec'}
        key = "dff_auc";
    case {'dff_max2','event.peak_amplitude_dff','curve_dffmax2','cell.mean_event_peak_amplitude_dff'}
        key = "dff_max2";
    case {'max_prop_speed_um_per_sec','max_prop_speed','event.max_propagation_speed_um_per_sec','propagation_maxpropspeed','cell.mean_event_max_propagation_speed_um_per_sec'}
        key = "max_prop_speed";
    case {'perimeter_um','perimeter_px','event.perimeter_um','cell.mean_perimeter_um'}
        key = "perimeter";
    case {'circularity','event.circularity','cell.mean_circularity'}
        key = "circularity";
    case {'propagation_distance_um','propagation_distance_px','event.propagation_distance_um','propagation_distance','propdist','cell.mean_propagation_distance_um'}
        key = "propagation_distance";
    case {'propagation_angle','propangle','event.propagation_angle_deg'}
        key = "propagation_angle";
    case {'event.frac_within_origin_cell','event.frac_outside_origin_cell','cell.mean_frac_outside_origin_cell'}
        key = "containment_mean";
    case {'cell.frac_escape_events_gt_threshold'}
        key = "containment_frac";
    case {'event.has_neighbor_recruitment','cell.frac_events_with_neighbor_recruitment','slice.frac_primary_events_with_neighbor_recruitment'}
        key = "neighbor_recruitment_frac";
    case {'event.mean_recruited_neighbor_entry_delay_sec','event.mean_neighbor_entry_delay_sec','cell.mean_neighbor_entry_delay_sec','slice.mean_neighbor_entry_delay_sec'}
        key = "neighbor_delay";
    case {'event.peak_recruited_neighbor_occupancy','event.cell_peak_overlap_frac_of_cell','cell.mean_neighbor_peak_overlap_frac_of_cell','slice.mean_neighbor_peak_overlap_frac_of_cell'}
        key = "neighbor_peak";
    case {'event.n_recruited_neighbor_partitions','event.n_neighbor_partitions','cell.mean_neighbor_partitions_per_event','slice.mean_neighbor_partitions_per_event'}
        key = "neighbor_partitions";
    otherwise
        key = string(id);
end
end
