function out = export_local_analysis_spider_mean_set(source, varargin)
%EXPORT_LOCAL_ANALYSIS_SPIDER_MEAN_SET Export a compact E-branch spider-mean set.
%
% Usage:
%   out = export_local_analysis_spider_mean_set(result)
%   out = export_local_analysis_spider_mean_set('.../analysis_result.mat')
%
% Default panel layout is chosen to cover event-relevant metric families
% without overloading any single spider:
%   1. core_event          : activity, intensity, kinetics
%   2. shape_and_motion    : shape, soma, propagation, containment
%   3. multicell_spread    : recruitment
%
% Synchrony is omitted in the default E-branch export because synchrony
% metrics are slice/network-level by construction rather than event-level.

result = load_local_analysis_result(source);

p = inputParser;
p.StructExpand = false;
addParameter(p, 'comparisonLevel', 'by_event', @(x) any(strcmpi(string(x), ["by_slice","by_event","by_branch"])));
addParameter(p, 'exportDir', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'subdirName', 'spider_mean_exports_E', @(x) ischar(x) || isstring(x));
addParameter(p, 'closeFigures', true, @islogical);
addParameter(p, 'exportFormat', 'png', @(x) any(strcmpi(string(x), ["png","pdf"])));
addParameter(p, 'showEntityPoints', true, @islogical);
addParameter(p, 'entityPointSize', 18, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'entityPointAlpha', 0.35, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
addParameter(p, 'panelSpecs', defaultPanelSpecs(), @(x) isstruct(x) && ~isempty(x));
parse(p, varargin{:});

exportRoot = resolveSpiderMeanExportDir(source, result, string(p.Results.exportDir), string(p.Results.subdirName));
if exist(exportRoot, 'dir') ~= 7
    mkdir(exportRoot);
end

panelSpecs = p.Results.panelSpecs;
allStats = table();
manifest = table();
metricGlossary = table();
exports = repmat(struct( ...
    'panelID', "", ...
    'title', "", ...
    'metricFamilies', strings(0,1), ...
    'figurePath', "", ...
    'statsFigurePath', "", ...
    'statsCsvPath', "", ...
    'metricGlossaryPath', "", ...
    'nMetrics', 0), 0, 1);

for i = 1:numel(panelSpecs)
    spec = panelSpecs(i);
    outMean = plot_local_analysis_spider_mean(result, ...
        'branchMode', 'E', ...
        'comparisonLevel', p.Results.comparisonLevel, ...
        'metricFamilies', cellstr(string(spec.metricFamilies(:))'), ...
        'showEntityPoints', p.Results.showEntityPoints, ...
        'entityPointSize', p.Results.entityPointSize, ...
        'entityPointAlpha', p.Results.entityPointAlpha, ...
        'showStatsFigure', true, ...
        'figureName', char(spec.figureName), ...
        'statsFigureName', char(spec.statsFigureName));

    if isempty(outMean.statsSummary) || ~istable(outMean.statsSummary) || height(outMean.statsSummary) == 0
        if p.Results.closeFigures
            closeIfGraphics(outMean);
        end
        continue
    end

    panelID = string(spec.id);
    familyLabel = strjoin(string(spec.metricFamilies(:))', ", ");
    statsPanel = outMean.statsSummary;
    statsPanel.panel_id = repmat(panelID, height(statsPanel), 1);
    statsPanel.panel_title = repmat(string(spec.title), height(statsPanel), 1);
    statsPanel.panel_families = repmat(familyLabel, height(statsPanel), 1);
    statsPanel = movevars(statsPanel, {'panel_id','panel_title','panel_families'}, 'Before', 1);
    allStats = [allStats; statsPanel]; %#ok<AGROW>

    figPath = fullfile(exportRoot, sprintf('%s.%s', spec.fileStem, lower(char(string(p.Results.exportFormat)))));
    statsFigPath = fullfile(exportRoot, sprintf('%s_stats.%s', spec.fileStem, lower(char(string(p.Results.exportFormat)))));
    statsCsvPath = fullfile(exportRoot, sprintf('%s_stats.csv', spec.fileStem));
    glossaryCsvPath = fullfile(exportRoot, sprintf('%s_metric_glossary.csv', spec.fileStem));

    exportgraphics(outMean.figure, figPath, 'Resolution', 300);
    if isfield(outMean, 'statsFigure') && ~isempty(outMean.statsFigure) && isgraphics(outMean.statsFigure)
        exportgraphics(outMean.statsFigure, statsFigPath, 'Resolution', 300);
    else
        statsFigPath = "";
    end
    writetable(statsPanel, statsCsvPath);
    glossaryPanel = buildMetricExplanationTable(statsPanel);
    writetable(glossaryPanel, glossaryCsvPath);
    metricGlossary = [metricGlossary; glossaryPanel]; %#ok<AGROW>

    exports(end+1,1) = struct( ... %#ok<AGROW>
        'panelID', panelID, ...
        'title', string(spec.title), ...
        'metricFamilies', string(spec.metricFamilies(:)), ...
        'figurePath', string(figPath), ...
        'statsFigurePath', string(statsFigPath), ...
        'statsCsvPath', string(statsCsvPath), ...
        'metricGlossaryPath', string(glossaryCsvPath), ...
        'nMetrics', height(statsPanel));

    manifest = [manifest; table(panelID, string(spec.title), familyLabel, height(statsPanel), string(figPath), string(statsCsvPath), string(glossaryCsvPath), ...
        'VariableNames', {'panel_id','panel_title','panel_families','n_metrics','figure_path','stats_csv_path','metric_glossary_csv_path'})]; %#ok<AGROW>

    if p.Results.closeFigures
        closeIfGraphics(outMean);
    end
end

combinedCsv = fullfile(exportRoot, 'spider_mean_E_combined_stats.csv');
manifestCsv = fullfile(exportRoot, 'spider_mean_E_manifest.csv');
glossaryCsv = fullfile(exportRoot, 'spider_mean_E_metric_glossary.csv');
glossaryMd = fullfile(exportRoot, 'spider_mean_E_metric_glossary.md');
if ~isempty(allStats)
    writetable(allStats, combinedCsv);
else
    combinedCsv = "";
end
if ~isempty(manifest)
    writetable(manifest, manifestCsv);
else
    manifestCsv = "";
end
if ~isempty(metricGlossary)
    metricGlossary = unique(metricGlossary, 'rows', 'stable');
    writetable(metricGlossary, glossaryCsv);
    writeHumanReadableGlossary(metricGlossary, glossaryMd);
else
    glossaryCsv = "";
    glossaryMd = "";
end

out = struct();
out.exportDir = exportRoot;
out.combinedStats = allStats;
out.manifest = manifest;
out.metricGlossary = metricGlossary;
out.exports = exports;
out.combinedStatsCsv = string(combinedCsv);
out.manifestCsv = string(manifestCsv);
out.metricGlossaryCsv = string(glossaryCsv);
out.metricGlossaryMarkdown = string(glossaryMd);
out.note = "Synchrony omitted from default E-branch spider_mean export because it is not event-level.";
end

function specs = defaultPanelSpecs()
specs = struct( ...
    'id', {'core_event','shape_and_motion','multicell_spread'}, ...
    'title', {'Core Event', 'Shape And Motion', 'Multicell Spread'}, ...
    'metricFamilies', { ...
        ["activity";"intensity";"kinetics"], ...
        ["shape";"soma";"propagation";"containment"], ...
        ["recruitment"]}, ...
    'figureName', {'BAPol_Mean_E_CoreEvent','BAPol_Mean_E_ShapeMotion','BAPol_Mean_E_Multicell'}, ...
    'statsFigureName', {'BAPol_Mean_E_CoreEvent_Stats','BAPol_Mean_E_ShapeMotion_Stats','BAPol_Mean_E_Multicell_Stats'}, ...
    'fileStem', {'spider_mean_E_core_event','spider_mean_E_shape_motion','spider_mean_E_multicell'});
end

function exportDir = resolveSpiderMeanExportDir(source, result, override, subdirName)
if strlength(override) > 0
    root = string(override);
elseif isstruct(result) && isfield(result, 'exportDir') && strlength(string(result.exportDir)) > 0
    root = string(result.exportDir);
elseif (ischar(source) || (isstring(source) && isscalar(source))) && exist(char(string(source)), 'file') == 2
    root = string(fileparts(char(string(source))));
elseif (ischar(source) || (isstring(source) && isscalar(source))) && exist(char(string(source)), 'dir') == 7
    root = string(source);
elseif isstruct(result) && isfield(result, 'pathIn') && strlength(string(result.pathIn)) > 0
    root = string(fullfile(char(string(result.pathIn)), 'derived', 'local_outputs', ['manual_' datestr(now, 'yyyymmdd_HHMMSS')]));
else
    error('export_local_analysis_spider_mean_set:NoExportDir', ...
        'Could not determine an export directory from the supplied source/result.');
end
exportDir = char(fullfile(char(root), char(subdirName)));
end

function closeIfGraphics(outMean)
if isstruct(outMean) && isfield(outMean, 'figure') && ~isempty(outMean.figure) && isgraphics(outMean.figure)
    close(outMean.figure);
end
if isstruct(outMean) && isfield(outMean, 'statsFigure') && ~isempty(outMean.statsFigure) && isgraphics(outMean.statsFigure)
    close(outMean.statsFigure);
end
end

function T = buildMetricExplanationTable(statsPanel)
metrics = string(statsPanel.metric(:));
families = string(statsPanel.metric_family(:));
panelID = string(statsPanel.panel_id(:));
panelTitle = string(statsPanel.panel_title(:));

T = table(metrics, families, panelID, panelTitle, ...
    'VariableNames', {'metric','metric_family','panel_id','panel_title'});

publicMetric = repmat("", height(T), 1);
displayUnit = repmat("", height(T), 1);
nativeUnit = repmat("", height(T), 1);
aggregationPath = repmat("", height(T), 1);
aggregationRule = repmat("", height(T), 1);
notes = repmat("", height(T), 1);

for i = 1:height(T)
    row = bwaLookupMetricVocabulary(T.metric(i));
    if isempty(row)
        continue
    end
    publicMetric(i) = firstStringOrEmpty(row, 'public_metric_id');
    displayUnit(i) = firstStringOrEmpty(row, 'display_unit');
    nativeUnit(i) = firstStringOrEmpty(row, 'native_unit');
    aggregationPath(i) = firstStringOrEmpty(row, 'aggregation_path');
    aggregationRule(i) = firstStringOrEmpty(row, 'aggregation_rule');
    notes(i) = firstStringOrEmpty(row, 'vocabulary_notes');
end

T.public_metric_id = publicMetric;
T.display_unit = displayUnit;
T.native_unit = nativeUnit;
T.aggregation_path = aggregationPath;
T.aggregation_rule = aggregationRule;
T.explanation = notes;
end

function out = firstStringOrEmpty(T, varName)
out = "";
if istable(T) && ismember(varName, T.Properties.VariableNames) && ~isempty(T.(varName))
    out = string(T.(varName)(1));
end
end

function writeHumanReadableGlossary(T, outFile)
if isempty(T) || ~istable(T)
    return
end

lines = strings(0,1);
lines(end+1) = "# Spider Mean E-Branch Metric Glossary";
lines(end+1) = "";
lines(end+1) = "Human-readable explanations for the metrics used in the exported E-branch spider_mean panels.";
lines(end+1) = "";
lines(end+1) = "Metrics are grouped and ordered to match the spider plots.";
lines(end+1) = "";

families = unique(string(T.metric_family), 'stable');
for i = 1:numel(families)
    fam = families(i);
    rows = T(string(T.metric_family) == fam, :);
    lines(end+1) = "## " + familyTitle(fam); %#ok<AGROW>
    lines(end+1) = ""; %#ok<AGROW>
    for j = 1:height(rows)
        metric = string(rows.metric(j));
        desc = explainSpiderMetric(metric);
        units = string(rows.display_unit(j));
        if strlength(units) == 0
            units = string(rows.native_unit(j));
        end
        aggRule = string(rows.aggregation_rule(j));
        if strlength(aggRule) == 0
            aggRule = inferAggregationRule(metric);
        end
        lines(end+1) = "### `" + metric + "`"; %#ok<AGROW>
        lines(end+1) = ""; %#ok<AGROW>
        lines(end+1) = desc.what; %#ok<AGROW>
        lines(end+1) = ""; %#ok<AGROW>
        if strlength(units) > 0
            lines(end+1) = "- Units: `" + units + "`"; %#ok<AGROW>
        end
        if strlength(aggRule) > 0
            lines(end+1) = "- Summary rule: " + aggRule; %#ok<AGROW>
        end
        if strlength(desc.interpretation) > 0
            lines(end+1) = "- Interpretation: " + desc.interpretation; %#ok<AGROW>
        end
        lines(end+1) = ""; %#ok<AGROW>
    end
end

fid = fopen(outFile, 'w');
if fid < 0
    error('export_local_analysis_spider_mean_set:OpenGlossaryFailed', ...
        'Could not open glossary markdown for writing: %s', outFile);
end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, strjoin(lines, newline), 'char');
end

function title = familyTitle(fam)
switch string(fam)
    case "activity"
        title = "Activity";
    case "intensity"
        title = "Intensity";
    case "kinetics"
        title = "Kinetics";
    case "shape"
        title = "Shape";
    case "soma"
        title = "Soma Relationship";
    case "propagation"
        title = "Propagation";
    case "containment"
        title = "Containment";
    case "recruitment"
        title = "Recruitment";
    otherwise
        title = char(string(fam));
end
end

function desc = explainSpiderMetric(metric)
metric = string(metric);
desc = struct('what', "Metric description unavailable.", 'interpretation', "");
switch metric
    case "dffAUC"
        desc.what = "Integrated event intensity over time. This summarizes the total fluorescence response across the event rather than just its peak.";
        desc.interpretation = "Higher values mean larger overall calcium signal load during the event.";
    case "dffMax2"
        desc.what = "Peak event amplitude. This captures the strongest fluorescence excursion reached by the event.";
        desc.interpretation = "Higher values mean stronger peak calcium elevation.";
    case "decayTau"
        desc.what = "Approximate decay time constant of the event signal after the peak.";
        desc.interpretation = "Higher values mean slower recovery or more prolonged signal decay.";
    case "duration"
        desc.what = "Total event duration from onset to end.";
        desc.interpretation = "Higher values mean events persist longer in time.";
    case "fall91"
        desc.what = "Late falling-phase time, measuring how long the signal takes to decline through the main decay range.";
        desc.interpretation = "Higher values mean a slower falling phase.";
    case "rise19"
        desc.what = "Early rising-phase time, measuring how long the signal takes to rise through the main activation range.";
        desc.interpretation = "Higher values mean slower event onset.";
    case "area"
        desc.what = "Spatial size of the event footprint.";
        desc.interpretation = "Higher values mean larger events.";
    case "circularity"
        desc.what = "Compactness of the event shape, comparing how circle-like the footprint is.";
        desc.interpretation = "Higher values mean more compact and round shapes; lower values mean more irregular or elongated shapes.";
    case "footprint"
        desc.what = "Total unique area occupied by the event over its lifetime.";
        desc.interpretation = "Higher values mean the event visited more spatial territory.";
    case "perimeter_um"
        desc.what = "Boundary length of the event footprint.";
        desc.interpretation = "Higher values mean a larger or more complex outline.";
    case "frac_away_events"
        desc.what = "Fraction of events whose motion is oriented away from the soma.";
        desc.interpretation = "Higher values mean outward-directed event motion is more common.";
    case "frac_net_away_soma_events"
        desc.what = "Fraction of events with a net displacement away from the soma over the full event.";
        desc.interpretation = "Higher values mean more events end farther from the soma than where they started.";
    case "frac_net_toward_soma_events"
        desc.what = "Fraction of events with a net displacement toward the soma over the full event.";
        desc.interpretation = "Higher values mean more events finish closer to the soma.";
    case "frac_soma_events"
        desc.what = "Fraction of events that involve the soma region itself.";
        desc.interpretation = "Higher values mean soma-involving events are more common.";
    case "frac_toward_events"
        desc.what = "Fraction of events with movement oriented toward the soma.";
        desc.interpretation = "Higher values mean inward-directed motion is more common.";
    case "frac_traverse_events"
        desc.what = "Fraction of events that traverse past the soma-relative axis or effectively cross through the local soma-centered geometry.";
        desc.interpretation = "Higher values mean more events span across the soma-centered spatial field rather than staying one-sided.";
    case "soma_end_dist_um"
        desc.what = "Distance from the soma to the event endpoint.";
        desc.interpretation = "Higher values mean events terminate farther from the soma.";
    case "soma_net_distance_change_um"
        desc.what = "Net change in distance from the soma between event start and end.";
        desc.interpretation = "Positive shifts indicate outward net motion; negative shifts indicate inward net motion.";
    case "soma_start_dist_um"
        desc.what = "Distance from the soma to the event starting location.";
        desc.interpretation = "Higher values mean events initiate farther from the soma.";
    case "centroid"
        desc.what = "Total centroid travel distance across the event.";
        desc.interpretation = "Higher values mean the event center of mass moves more over time.";
    case "maxPropSpeed"
        desc.what = "Maximum observed propagation speed during the event.";
        desc.interpretation = "Higher values mean faster local spread at the peak propagation moment.";
    case "propagation_distance_um"
        desc.what = "Overall event propagation reach or displacement magnitude.";
        desc.interpretation = "Higher values mean the event spreads farther across space.";
    case "propangle"
        desc.what = "Preferred propagation direction angle.";
        desc.interpretation = "This is directional rather than simply larger/smaller; changes indicate a shift in the dominant direction of spread.";
    case "soma_centroid_dist_um"
        desc.what = "Characteristic distance between the event centroid trajectory and the soma.";
        desc.interpretation = "Higher values mean centroid motion occurs farther from the soma.";
    case "frac_outside_origin_cell"
        desc.what = "Fraction of the event footprint that lies outside the originating cell.";
        desc.interpretation = "Higher values mean weaker containment within the origin cell and greater escape into surrounding territory.";
    case "cells_per_multicell_event"
        desc.what = "Average number of participating cells for events that recruit more than one cell.";
        desc.interpretation = "Higher values mean multicell events recruit broader cell ensembles.";
    case "frac_events_with_neighbor_recruitment"
        desc.what = "Fraction of events that recruit at least one neighboring cell.";
        desc.interpretation = "Higher values mean intercellular spread is more common.";
    case "frac_multicell_events"
        desc.what = "Fraction of all events that involve more than one cell.";
        desc.interpretation = "Higher values mean multicell propagation is more common.";
    case "neighbor_entry_delay_sec"
        desc.what = "Delay between origin-cell event onset and entry into a neighboring cell.";
        desc.interpretation = "Higher values mean slower recruitment of neighboring cells.";
    case "neighbor_partitions_per_event"
        desc.what = "Number of distinct recruited-neighbor partitions or recruited cell compartments per event.";
        desc.interpretation = "Higher values mean broader or more segmented multicell recruitment.";
    case "neighbor_peak_overlap_frac_of_cell"
        desc.what = "Peak fraction of a recruited neighboring cell occupied or overlapped by the event.";
        desc.interpretation = "Higher values mean deeper penetration into recruited neighboring cells.";
end
end

function rule = inferAggregationRule(metric)
switch string(metric)
    case {"dffAUC","dffMax2","decayTau","duration","fall91","rise19","area","circularity","footprint","perimeter_um", ...
            "frac_away_events","frac_net_away_soma_events","frac_net_toward_soma_events","frac_soma_events","frac_toward_events","frac_traverse_events", ...
            "soma_end_dist_um","soma_net_distance_change_um","soma_start_dist_um","centroid","maxPropSpeed","propagation_distance_um","propangle","soma_centroid_dist_um", ...
            "frac_outside_origin_cell","neighbor_entry_delay_sec","neighbor_partitions_per_event","neighbor_peak_overlap_frac_of_cell"}
        rule = "Event-level summary, pooled across post events for the E branch.";
    case {"cells_per_multicell_event","frac_events_with_neighbor_recruitment","frac_multicell_events"}
        rule = "Event-level recruitment summary pooled across post events for the E branch.";
    otherwise
        rule = "";
end
end
