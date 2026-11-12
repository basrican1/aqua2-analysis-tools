function render_dg_anatomy_overview(anatomy, outPath)
% Render structural evidence and the inferred DG anatomy.

fig = figure('Color', 'w', 'Position', [20 20 1800 1320], 'Visible', 'off');
tiledlayout(fig, 3, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

show_gray(anatomy.projections.maximumDisplay, 'Unmasked maximum projection');
show_gray(anatomy.projections.upperDisplay, 'Upper-percentile projection');
show_gray(anatomy.structural, 'Structural composite');
show_gray(anatomy.gclBandEvidence, 'Low-frequency GCL-band evidence');

plot_orientation(anatomy);
plot_process_candidates(anatomy);
plot_clean_anatomy(anatomy);
plot_confidence(anatomy);

plot_activity_anatomy(anatomy);
plot_proposals(anatomy, anatomy.wholeCellProposals, 'Whole-cell proposals');
plot_proposals(anatomy, anatomy.distalProposals, 'Distal-process proposals');
plot_metadata(anatomy);

sgtitle(strrep(anatomy.metadata.label, '_', ' '), ...
    'Interpreter', 'none', 'FontWeight', 'bold');
exportgraphics(fig, outPath, 'Resolution', 220);
close(fig);
end

function show_gray(img, ttl)
nexttile;
imagesc(img, [0 1]);
axis image off
colormap(gca, gray(256));
title(ttl, 'FontWeight', 'normal');
end

function plot_orientation(anatomy)
nexttile;
imagesc(anatomy.projections.upperDisplay, [0 1]);
axis image off
colormap(gca, gray(256));
title('Local process orientation', 'FontWeight', 'normal');
hold on
step = 24;
[xx, yy] = meshgrid(1:step:size(anatomy.structural, 2), ...
    1:step:size(anatomy.structural, 1));
idx = sub2ind(size(anatomy.structural), yy, xx);
theta = anatomy.orientation(idx);
support = anatomy.ridgeStrength(idx) .* anatomy.orientationCoherence(idx);
keep = isfinite(theta) & support > prctile(support(isfinite(support)), 45);
u = 9 .* cos(theta);
v = 9 .* sin(theta);
quiver(xx(keep) - u(keep), yy(keep) - v(keep), 2 .* u(keep), 2 .* v(keep), 0, ...
    'Color', [0 0.95 0.85], 'LineWidth', 0.7, 'MaxHeadSize', 0);
hold off
end

function plot_process_candidates(anatomy)
nexttile;
imagesc(anatomy.projections.maximumDisplay, [0 1]);
axis image off
colormap(gca, gray(256));
title('Process, cell-boundary, and soma evidence', 'FontWeight', 'normal');
hold on
if isfield(anatomy, 'sgzCandidatePoints') && ~isempty(anatomy.sgzCandidatePoints)
    scatter(anatomy.sgzCandidatePoints(:, 1), anatomy.sgzCandidatePoints(:, 2), ...
        7, [0 1 1], 'o');
end
for i = 1:numel(anatomy.processCandidates)
    c = anatomy.processCandidates(i);
    color = confidence_color(c.confidence);
    plot([c.base(1), c.distal(1)], [c.base(2), c.distal(2)], '-', ...
        'Color', color, 'LineWidth', 1.2);
    plot(c.base(1), c.base(2), 'o', 'Color', color, ...
        'MarkerFaceColor', color, 'MarkerSize', 4);
end
if isfield(anatomy, 'roiEvidence') && anatomy.roiEvidence.isAvailable
    for i = 1:numel(anatomy.roiEvidence.cellBoundaries)
        b = anatomy.roiEvidence.cellBoundaries{i};
        plot([b(:, 1); b(1, 1)], [b(:, 2); b(1, 2)], '-', ...
            'Color', [0.75 0.20 0.95], 'LineWidth', 0.8);
    end
    for i = 1:numel(anatomy.roiEvidence.somaBoundaries)
        b = anatomy.roiEvidence.somaBoundaries{i};
        plot([b(:, 1); b(1, 1)], [b(:, 2); b(1, 2)], '-', ...
            'Color', [0.20 1.00 0.35], 'LineWidth', 1.0);
    end
    for i = 1:numel(anatomy.roiProcessCandidates)
        c = anatomy.roiProcessCandidates(i);
        plot([c.base(1), c.distal(1)], [c.base(2), c.distal(2)], '-', ...
            'Color', [1.00 0.80 0.05], 'LineWidth', 1.5);
        plot(c.base(1), c.base(2), 'o', 'Color', [0.20 1.00 0.35], ...
            'MarkerFaceColor', [0.20 1.00 0.35], 'MarkerSize', 4);
    end
end
hold off
end

function plot_clean_anatomy(anatomy)
nexttile;
imagesc(anatomy.projections.maximumDisplay, [0 1]);
axis image off
colormap(gca, gray(256));
title('Inferred DG anatomy', 'FontWeight', 'normal');
hold on
overlay_layers(anatomy.layerLabels, 0.13);
plot_blades(anatomy.blades, true);
plot_shared_apex(anatomy);
hold off
end

function plot_confidence(anatomy)
nexttile;
imagesc(anatomy.anatomyConfidence, [0 1]);
axis image off
colormap(gca, parula(256));
title('Anatomy confidence', 'FontWeight', 'normal');
hold on
plot_blades(anatomy.blades, false);
plot_shared_apex(anatomy);
plot_blade_confidence_labels(anatomy.blades);
hold off
colorbar;
end

function plot_activity_anatomy(anatomy)
nexttile;
imagesc(anatomy.activityScore, [0 1]);
axis image off
colormap(gca, hot(256));
title('Activity evidence with SGZ estimate', 'FontWeight', 'normal');
hold on
plot_blades(anatomy.blades, false);
plot_shared_apex(anatomy);
hold off
colorbar;
end

function plot_proposals(anatomy, proposals, ttl)
nexttile;
imagesc(anatomy.projections.maximumDisplay, [0 1]);
axis image off
colormap(gca, gray(256));
title(ttl, 'FontWeight', 'normal');
hold on
plot_blades(anatomy.blades, false);
plot_shared_apex(anatomy);
for i = 1:numel(proposals)
    p = proposals(i).poly;
    plot([p(:, 1); p(1, 1)], [p(:, 2); p(1, 2)], '-', ...
        'Color', [1 0.82 0.05], 'LineWidth', 1.4);
    text(proposals(i).center(1), proposals(i).center(2), sprintf('%d', i), ...
        'Color', [1 0.82 0.05], 'FontWeight', 'bold', 'FontSize', 9);
end
hold off
end

function plot_metadata(anatomy)
nexttile;
axis off
meanConfidence = 0;
if ~isempty(anatomy.blades)
    meanConfidence = mean([anatomy.blades.confidence]);
end
apexStatus = 'not applicable';
if isfield(anatomy, 'apexEstimate') && all(isfinite(anatomy.apexEstimate))
    if anatomy.apexVisible
        apexStatus = sprintf('visible at (%.0f, %.0f)', anatomy.apexEstimate);
    else
        apexStatus = sprintf('projected at (%.0f, %.0f)', anatomy.apexEstimate);
    end
end
lines = { ...
    sprintf('Animal: %s', anatomy.metadata.animal), ...
    sprintf('Slice: %s', anatomy.metadata.slice), ...
    sprintf('Condition: %s', anatomy.metadata.condition), ...
    sprintf('Stimulation: %s', anatomy.metadata.stimulation), ...
    sprintf('Tamoxifen: %d dpi', anatomy.metadata.tamoxDpi), ...
    sprintf('Movie: %d x %d x %d', anatomy.metadata.movieSize), ...
    sprintf('Process candidates: %d', numel(anatomy.processCandidates)), ...
    sprintf('Blade hypotheses: %d', numel(anatomy.blades)), ...
    sprintf('Shared apex: %s', apexStatus), ...
    sprintf('Mean blade confidence: %.2f', meanConfidence), ...
    sprintf('Mode: %s', anatomy.estimationMode), ...
    '', ...
    'Dotted cyan: inferred hilus/SGZ-GCL boundary', ...
    'Orange guides: predicted radial direction', ...
    'Purple / green: cell / soma annotations', ...
    'Layer colors: hilus / SGZ / GCL / molecular'};
text(0.03, 0.96, strjoin(lines, newline), 'VerticalAlignment', 'top', ...
    'FontName', 'Menlo', 'FontSize', 10, 'Interpreter', 'none');
end

function overlay_layers(labels, alphaValue)
colors = single([ ...
    0.12 0.28 0.95; ...
    0.00 0.95 0.95; ...
    0.20 0.88 0.25; ...
    1.00 0.52 0.05]);
rgb = zeros([size(labels), 3], 'single');
alpha = zeros(size(labels), 'single');
for k = 1:4
    mask = labels == k;
    for ch = 1:3
        plane = rgb(:, :, ch);
        plane(mask) = colors(k, ch);
        rgb(:, :, ch) = plane;
    end
    alpha(mask) = alphaValue;
end
h = image(rgb);
set(h, 'AlphaData', alpha);
end

function plot_blades(blades, showRadials)
for b = 1:numel(blades)
    curve = blades(b).curve;
    alphaColor = confidence_color(blades(b).confidence);
    plot(curve(:, 1), curve(:, 2), ':', 'Color', [0 1 1], 'LineWidth', 2.0);
    if showRadials
        indices = round(linspace(20, size(curve, 1) - 20, 5));
        indices = indices(indices >= 1 & indices <= size(curve, 1));
        for k = indices
            p0 = curve(k, 1:2) + 8 .* curve(k, 3:4);
            p1 = curve(k, 1:2) + 68 .* curve(k, 3:4);
            plot([p0(1), p1(1)], [p0(2), p1(2)], '-', ...
                'Color', alphaColor, 'LineWidth', 0.9);
        end
    end
end
end

function plot_shared_apex(anatomy)
if ~isfield(anatomy, 'apexEstimate') || ~anatomy.apexVisible
    return
end
plot(anatomy.apexEstimate(1), anatomy.apexEstimate(2), 'o', ...
    'Color', [0 1 1], 'MarkerFaceColor', [0 1 1], 'MarkerSize', 5);
end

function plot_blade_confidence_labels(blades)
xl = xlim(gca);
yl = ylim(gca);
for i = 1:numel(blades)
    curve = blades(i).curve;
    visible = curve(:, 1) >= min(xl) & curve(:, 1) <= max(xl) & ...
        curve(:, 2) >= min(yl) & curve(:, 2) <= max(yl);
    indices = find(visible);
    if isempty(indices)
        continue
    end
    k = indices(round(0.62 * numel(indices)));
    text(curve(k, 1), curve(k, 2), sprintf(' %.2f ', blades(i).confidence), ...
        'Color', [1 1 1], 'BackgroundColor', [0.08 0.08 0.08], ...
        'FontSize', 8, 'FontWeight', 'bold', 'Margin', 1);
end
end

function color = confidence_color(value)
value = min(max(value, 0), 1);
color = (1 - value) .* [0.95 0.25 0.20] + value .* [1.00 0.78 0.05];
end
