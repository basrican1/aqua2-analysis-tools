function reg = bwaEnrichMetricRegistry(reg)
%BWAENRICHMETRICREGISTRY Add family/lineage metadata to the local metric registry.

if ~istable(reg) || isempty(reg) || ~ismember('metric_id', reg.Properties.VariableNames)
    return
end

catalog = bwaMetricVocabularyCatalog();
reg.metric_id = string(reg.metric_id);

[lia, loc] = ismember(reg.metric_id, string(catalog.metric_id));

reg = ensureStringColumn(reg, 'public_metric_id');
reg = ensureStringColumn(reg, 'family');
reg = ensureStringColumn(reg, 'native_unit');
reg = ensureStringColumn(reg, 'display_unit');
reg = ensureStringColumn(reg, 'aggregation_path');
reg = ensureStringColumn(reg, 'aggregation_rule');
reg = ensureStringColumn(reg, 'vocabulary_notes');
reg = ensureStringColumn(reg, 'storage_status');
reg = ensureStringColumn(reg, 'storage_level');
reg = ensureStringColumn(reg, 'metric_class');
reg = ensureStringColumn(reg, 'default_slice_summary');
reg = ensureStringColumn(reg, 'alternate_slice_summary');
reg = ensureStringColumn(reg, 'storage_notes');

reg.public_metric_id(~lia) = reg.metric_id(~lia);
reg.family(~lia) = "";
reg.native_unit(~lia) = string(reg.unit_level(~lia));
reg.display_unit(~lia) = string(reg.unit_level(~lia));
reg.aggregation_path(~lia) = string(reg.unit_level(~lia));
reg.aggregation_rule(~lia) = "";
reg.vocabulary_notes(~lia) = "";
reg.storage_status(~lia) = "stored";
reg.storage_level(~lia) = string(reg.unit_level(~lia));
reg.metric_class(~lia) = "native";
reg.default_slice_summary(~lia) = "none";
reg.alternate_slice_summary(~lia) = "none";
reg.storage_notes(~lia) = "";

hitIdx = find(lia);
for i = 1:numel(hitIdx)
    r = hitIdx(i);
    c = loc(r);
    reg.public_metric_id(r) = string(catalog.public_metric_id(c));
    reg.family(r) = string(catalog.family(c));
    reg.native_unit(r) = string(catalog.native_unit(c));
    reg.display_unit(r) = string(catalog.display_unit(c));
    reg.aggregation_path(r) = string(catalog.aggregation_path(c));
    reg.aggregation_rule(r) = string(catalog.aggregation_rule(c));
    reg.vocabulary_notes(r) = string(catalog.notes(c));
    reg.storage_status(r) = string(catalog.storage_status(c));
    reg.storage_level(r) = string(catalog.storage_level(c));
    reg.metric_class(r) = string(catalog.metric_class(c));
    reg.default_slice_summary(r) = string(catalog.default_slice_summary(c));
    reg.alternate_slice_summary(r) = string(catalog.alternate_slice_summary(c));
    reg.storage_notes(r) = string(catalog.storage_notes(c));
end
end

function T = ensureStringColumn(T, name)
if ~ismember(name, T.Properties.VariableNames)
    T.(name) = strings(height(T), 1);
else
    T.(name) = string(T.(name));
end
end
