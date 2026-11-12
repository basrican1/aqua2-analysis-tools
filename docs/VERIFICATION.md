# Verification scope

MATLAB R2024b Update 6: **35 selected synthetic-input tests passed, zero failed/incomplete**. The actual artificial geometry/normalization workflow also completed. The selected suite covers epoch assignment, frequency and per-slice metric summaries, QC creation/application, group comparisons, subset selection, correlations, metric vocabulary, heatmap rendering and paired correlation plots.

The initial larger test attempt exposed the absent external `spider_plot` dependency. Rendering tests for mean/slice spider aggregation are therefore optional in `run_public_tests`: they run when that helper is installed, and the runner explicitly reports their absence otherwise. No helper was downloaded or substituted to make those tests pass. The source remains available and the existing historical scientific spider figures remain distinct evidence.

Other source tests are present but were not all executed for this release; some require optional shared modules. Real-data regression tests/private fixtures and a study-specific timing fixture are excluded. No upstream AQuA2 detection, raw research run, empirical figure regeneration, new scientific inference or model validation was performed.

Anatomical estimation and refined-propagation modules are source-inspected diagnostic/prototype implementations. No native runtime claim is made for those modules in this release. Technical plotting and successful arithmetic tests do not demonstrate biological propagation, establish cell identity or validate statistical independence.
