function [order, familiesOrdered, ranks] = bwaOrderSpiderMetricsByFamily(metricLabels)
%BWAORDERSPIDERMETRICSBYFAMILY Order spider metrics by a stable family sequence.

labels = string(metricLabels(:));
families = string(bwaResolveMetricFamily(labels));
families = families(:);
if numel(families) ~= numel(labels)
    families = repmat("", numel(labels), 1);
end

familyOrder = [ ...
    "activity"
    "intensity"
    "kinetics"
    "shape"
    "soma"
    "propagation"
    "containment"
    "recruitment"
    "synchrony"
    "uncategorized"
    ""];

ranks = zeros(numel(labels), 1);
for i = 1:numel(labels)
    hit = find(familyOrder == families(i), 1, 'first');
    if isempty(hit)
        ranks(i) = numel(familyOrder) + 1;
    else
        ranks(i) = hit;
    end
end

tbl = table((1:numel(labels))', ranks, families, labels, ...
    'VariableNames', {'idx','rank','family','label'});
tbl = sortrows(tbl, {'rank','family','label','idx'});
order = tbl.idx;
familiesOrdered = families(order);
end
