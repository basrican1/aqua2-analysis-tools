function T = bwaCollectDerivedSummaryEntries(S, includeOverlappingDerived)
keys = strings(0,1);
titles = strings(0,1);
comps = cell(0,1);
overlapKeys = {'basic_area','curve_duration','curve_rise19','curve_fall91', ...
    'curve_decayTau','curve_dffAUC','curve_dffMax2','propagation_maxPropSpeed'};

if isfield(S, 'baviolin') && isstruct(S.baviolin)
    fn = fieldnames(S.baviolin);
    for i = 1:numel(fn)
        if ~includeOverlappingDerived && any(strcmp(fn{i}, overlapKeys))
            continue
        end
        c = S.baviolin.(fn{i});
        if isstruct(c) && isfield(c, 'dffAll') && isfield(c, 'epMembIdx')
            keys(end+1,1) = string(['derived_' fn{i}]); %#ok<AGROW>
            titles(end+1,1) = string(strrep(fn{i}, '_', '.')); %#ok<AGROW>
            comps{end+1,1} = c; %#ok<AGROW>
        end
    end
end

if isfield(S, 'centroid') && isstruct(S.centroid) && isfield(S.centroid, 'distComp')
    c = S.centroid.distComp;
    if isstruct(c) && isfield(c, 'dffAll') && isfield(c, 'epMembIdx')
        keys(end+1,1) = "derived_centroid_dist"; %#ok<AGROW>
        titles(end+1,1) = "centroid.distTot"; %#ok<AGROW>
        comps{end+1,1} = c; %#ok<AGROW>
    end
end
if isfield(S, 'propagation') && isstruct(S.propagation)
    if isfield(S.propagation, 'distance')
        c = S.propagation.distance;
        if isstruct(c) && isfield(c, 'dffAll') && isfield(c, 'epMembIdx')
            keys(end+1,1) = "derived_prop_distance"; %#ok<AGROW>
            titles(end+1,1) = "propagation.propdist"; %#ok<AGROW>
            comps{end+1,1} = c; %#ok<AGROW>
        end
    end
end
if isfield(S, 'footprint') && isstruct(S.footprint) && isfield(S.footprint, 'dffAll')
    keys(end+1,1) = "derived_footprint"; %#ok<AGROW>
    titles(end+1,1) = "propagation.pCell"; %#ok<AGROW>
    comps{end+1,1} = S.footprint; %#ok<AGROW>
end

T = table(keys, titles, comps, 'VariableNames', {'key','title','comp'});
end
