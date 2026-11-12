function polData = filter_spider_metric_families(polData, varargin)
%FILTER_SPIDER_METRIC_FAMILIES Keep only selected spider metric families.

p = inputParser;
p.StructExpand = false;
addParameter(p, 'metricFamilies', 'all', @(x) ischar(x) || isstring(x) || iscell(x));
parse(p, varargin{:});

families = lower(string(p.Results.metricFamilies));
if ~isstruct(polData) || ~all(isfield(polData, {'labels','valuesA','valuesC'}))
    return
end

labels = string(polData.labels(:));
legacyDrop = lower(labels) == "propdist";
if isempty(families) || any(families == "all")
    keep = ~legacyDrop;
    if ~all(keep)
        polData = subsetSpiderPolars(polData, labels, keep);
    end
    polData = orderSpiderPolars(polData);
    return
end
metricFamilies = bwaResolveMetricFamily(labels);
keep = ismember(metricFamilies, families);
if ~any(keep)
    return
end

% Drop legacy propagation-distance spider spoke. Keep the canonical
% propagation_distance_um / cell.mean_propagation_distance_um paths instead.
keep = keep & ~legacyDrop;
if ~any(keep)
    return
end

polData = subsetSpiderPolars(polData, labels, keep);
polData = orderSpiderPolars(polData);

keptFamilies = metricFamilies(keep);
synchronyCols = keptFamilies == "synchrony";
if any(synchronyCols)
    if isnumeric(polData.valuesA) && size(polData.valuesA,2) == numel(keptFamilies)
        polData.valuesA(:, synchronyCols) = NaN;
    end
    if isfield(polData, 'sigMaskA') && isnumeric(polData.sigMaskA) && size(polData.sigMaskA,2) == numel(keptFamilies)
        polData.sigMaskA(:, synchronyCols) = false;
    end
end
end

function polData = orderSpiderPolars(polData)
if ~isstruct(polData) || ~isfield(polData, 'labels') || isempty(polData.labels)
    return
end
labels = string(polData.labels(:));
[order, ~] = bwaOrderSpiderMetricsByFamily(labels);
polData = subsetSpiderPolars(polData, labels, order);
end

function polData = subsetSpiderPolars(polData, labels, keep)
if islogical(keep)
    idx = find(keep);
else
    idx = keep(:);
end
polData.labels = cellstr(labels(keep));
if isnumeric(polData.valuesA) && size(polData.valuesA,2) == numel(labels)
    polData.valuesA = polData.valuesA(:, idx);
end
if isnumeric(polData.valuesC) && size(polData.valuesC,2) == numel(labels)
    polData.valuesC = polData.valuesC(:, idx);
end
for f = ["entityValuesA","entityValuesC","entityValuesAByFile","entityValuesCByFile"]
    if isfield(polData, char(f)) && iscell(polData.(char(f))) && numel(polData.(char(f))) == numel(labels)
        polData.(char(f)) = polData.(char(f))(idx);
    end
end
for f = ["sigMaskA","sigMaskC","panelLimitsLog2"]
    if isfield(polData, char(f)) && isnumeric(polData.(char(f))) && size(polData.(char(f)),2) == numel(labels)
        polData.(char(f)) = polData.(char(f))(:, idx);
    end
end
end
