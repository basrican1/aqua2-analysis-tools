# AQUA2 analysis tools

MATLAB extensions for turning calcium-imaging event outputs into linked scientific measurements, quality-controlled summaries and readable figures. This curated source release includes canonical ingestion, event/cell geometry, soma-distance and partition measurements, normalization, epoch summaries, group statistics, correlation analysis, multidimensional profiles, anatomical estimation and exploratory spatial diagnostics.

Explore the [scientific pipeline](https://basrican1.github.io/explore-research-software-aqua2.html), [cells/anatomy/compartments and measurements](https://basrican1.github.io/explore-research-software-aqua2-metrics.html), [59-metric inventory and examples](https://basrican1.github.io/explore-research-software-aqua2-metrics-data.html), and [statistical outputs](https://basrican1.github.io/explore-research-software-aqua2-summary-data.html).

## Start with implemented source

| Area | Entry points | Scientific role |
|---|---|---|
| Ingest | [ingestAQuA2Canonical](ingest/ingestAQuA2Canonical.m), [ingestCanonicalFromBWA](ingest/ingestCanonicalFromBWA.m) | Linked file, cell, event and event-to-cell representations; original event identities retained |
| Measurements | [canonicalEventMetrics](ingest/canonicalEventMetrics.m), [canonicalEventCellPartitionMetrics](ingest/canonicalEventCellPartitionMetrics.m), [canonicalSomaDistanceMetrics](ingest/canonicalSomaDistanceMetrics.m) | Event geometry, partition involvement and soma-referenced quantities |
| QC | [applyQCToCanonical](analysis/QC/applyQCToCanonical.m), [createQCTables](analysis/QC/createQCTables.m) | Explicit review decisions and filtered views |
| Normalization | [normalizeCanonicalOutput](analysis/Canonical/normalizeCanonicalOutput.m), [baseline computation](analysis/Utils/bwaComputeBaselineValue.m) | Cell/global baselines, log2-fold, delta and percent-change transforms |
| Statistics | [runGroupComparison](analysis/Stats/runGroupComparison.m), [epoch vs baseline](analysis/Stats/summarizeMetricEpochVsBaseline.m) | Configurable comparisons, units and pairing; implemented tests do not establish suitability for every experiment |
| Integrated output | [legacy multipanel renderer](analysis/Workflow/plot_local_analysis_legacy_panels.m), [summary renderer](analysis/Workflow/plot_local_analysis_summary.m) | Current rendering family corresponding to historical dashboard demonstrations; historical plots were not regenerated |
| Profiles | [mean spider profiles](analysis/Workflow/plot_local_analysis_spider_mean.m), [slice spider](analysis/Plots/plotBAviolinSubplotSliceSpider.m) | Multidimensional measurements with explicit scale and aggregation choices |
| Relationships | [CorrelationStats](analysis/CorrelationStats/CorrelationStats.m), [heatmap](analysis/CorrelationStats/plot_CorrelationStats_heatmap.m), [paired display](analysis/CorrelationStats/plot_CorrelationStats_pair.m) | Unit-specific correlation tables and visual inspection |
| DG anatomy | [estimate_dg_anatomy_from_movie](anatomy/estimate_dg_anatomy_from_movie.m), [render_dg_anatomy_overview](anatomy/render_dg_anatomy_overview.m) | Diagnostic DG/GCL hypotheses, local direction and soft ROI evidence; no validated segmentation claim |
| Spatial prototypes | [event/cell graph](analysis/Workflow/MulticellGraph/build_event_cell_adjacency_graph.m), [refined pixel onset](experimental/refined-timing/estimate_pixel_onset_time.m) | Exploratory event geometry/timing; technical arrows are not proof of biological propagation |
| Computing | [preset reader](hpc/bwa_aqua2_read_preset.m), [OME/ROI validator](hpc/bwa_aqua2_validate.m), [Slurm template](hpc/submit_aqua2.sl) | Explicit preset parsing and user-configured infrastructure |

## Installation and artificial example

Use MATLAB R2024b or a compatible version. From this repository root:

```matlab
setup_aqua2_tools
result = run_synthetic_workflow
results = run_public_tests
```

The [synthetic workflow](examples/run_synthetic_workflow.m) exercises actual geometry and within-cell baseline normalization on artificial values. It reads no research dataset and makes no biological inference. `run_public_tests` selects existing tests that create their own small fixtures. Other historical tests are included as source and may require optional dependencies. See [verification scope](docs/VERIFICATION.md).

The source is a curated collection, not a promise that every workflow is turnkey. [Dependencies and execution notes](docs/DEPENDENCIES.md) distinguish source included here from required external components. Some plotting paths require helpers such as violin/spider plotting utilities; the repository does not redistribute them. Real AQuA2 ingest requires compatible upstream output structures. Real anatomical estimation requires image-processing functionality and a user-supplied registered TIFF; its outputs remain hypotheses requiring human review.

## Contributions and upstream boundaries

The original [AQuA2 software](https://github.com/yu-lab-vt/AQuA2) and its authors own the upstream event-detection methods; this repository presents the canonical interpretation, metrics, post-analysis, presentation and prototype extensions separately. Upstream AQuA2 code and lab fork history are not bundled. Some functions consume upstream AQuA2 structures and depend on that independently installed software.

The legacy batch-engine wrapper is withheld pending upstream derivation/redistribution review. The included Slurm launcher calls a separately installed function selected by the user; it does not ship or silently select an engine. No trained cell/compartment segmentation model is included. Current manual annotation/editor implementations are not all part of this source selection; the Explorer distinguishes those demonstrations from the released subset.

Existing legacy `CNO1`, `CNO2` and `CNO_post` programmatic epoch keys remain where needed for compatibility. They are configurable labels, not stimulation prescriptions. Their presence does not identify an unpublished experimental combination; do not equate administration periods with stimulation without actual source metadata. Comment examples use generic labels. Two study-bound tests and private regression references are not distributed.

See [NOTICES](NOTICES.md). No new open-source license is assigned. Publication permission and scientific qualification are distinct.
