function result = run_synthetic_workflow(outDir)
% Artificial arithmetic/geometry demonstration; no biological data or inference.
% Optional output directory must not already contain synthetic_result.json.
setup_aqua2_tools();
spatialSize = [6 6 1];
eventPixels = sub2ind(spatialSize(1:2), [2 2 3], [2 3 3]);
cellMask = false(6); cellMask(2,2:3) = true;
result.geometry = canonicalEventMetrics(eventPixels, spatialSize, cellMask, [], []);
assert(result.geometry.projected_area_npix == 3);
assert(abs(result.geometry.cell_overlap_frac_of_event - 2/3) < 1e-12);
canonical.frequency.epochTable = table([1;2], string(["Baseline";"Epoch2"]), ...
    'VariableNames', {'epoch_idx','epoch_name'});
canonical.frequency.ratesByCellEpoch = table(string(["c1";"c1";"c2";"c2"]), ...
    [1;2;1;2], string(["Baseline";"Epoch2";"Baseline";"Epoch2"]), [2;4;10;5], ...
    'VariableNames', {'cell_uid','epoch_idx','epoch_name','rate_events_per_min'});
[normOut, meta] = normalizeCanonicalOutput(canonical, ...
    'mode', 'by_cell_baseline_mean', 'transform', 'log2_fold');
actual = normOut.frequency.ratesByCellEpoch.rate_events_per_min;
assert(max(abs(actual - [0;1;0;-1])) < 1e-12);
result.synthetic = true;
result.scope = 'Artificial geometry and within-cell log2-fold calculation; no experimental results or method approval.';
result.normalized_values = actual;
result.normalization = meta;
if nargin > 0 && ~isempty(outDir)
    if ~isfolder(outDir), mkdir(outDir); end
    target = fullfile(outDir, 'synthetic_result.json');
    assert(~isfile(target), 'Refusing to overwrite an existing output.');
    fid = fopen(target, 'w'); assert(fid >= 0, 'Cannot open output.');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '%s', jsonencode(result, 'PrettyPrint', true));
end
end
