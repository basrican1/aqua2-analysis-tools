# Dependencies and execution scope

The setup function adds only this export's source folders. It never searches a private home directory, installs software, downloads data or contacts a cluster.

- MATLAB R2024b is the verification target. Statistical routines can require Statistics and Machine Learning Toolbox; anatomical/image operations can require Image Processing Toolbox. See function calls for exact path-specific dependencies.
- AQuA2 upstream structures are inputs to ingest; install the upstream software separately under its own terms. No upstream engine or GPL source is redistributed here.
- ImageJ ROI ZIP parsing, Bio-Formats, violin/spider plotting helpers and related plotting utilities are optional external dependencies, not bundled. Paths that call them require separate compatible installations.
- SharedAnalysisMetadata integration is an optional dependency. No private shared metadata or experiment manifests ship here.
- `bwa_aqua2_validate` requires explicit `BIOFORMATS_ROOT` and `AQUA2_HELPERS_ROOT`; the latter must provide the source-specific OME calibration helper. `bwa_aqua2_read_preset` is an independent CSV reader.
- Slurm is optional. Run `sbatch hpc/submit_aqua2.sl` from the repository root only after configuring your own resources, MATLAB path, external runner, input and preset. The launcher does not infer scientific settings or overwrite private environment configuration.
- Refined timing, wavefront, corridor and anatomical code is exploratory. Parameter defaults are software defaults, not approved biological models. Preserve source data, units and human review.

The source tree preserves original function names and relative organization where possible. Curated components may still refer to optional helpers absent from this release. This is source access with qualified examples, not a current end-to-end replication of the historical empirical figures.
