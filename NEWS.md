# karioCaS 0.99.2

## New Features

* **Group overlay plots are now the default output** of `taxa_retention()`
  (001), `reads_per_taxa()` (003) and `optimize_CS()` (006). Instead of one set
  of PDFs per sample, each function draws a single figure per biological group
  in which every sample is a faint line and the group mean (+/-SD band) is
  highlighted, faceted by Domain. `optimize_CS()` additionally marks each
  domain's median Primary Stability Index as a dashed reference line. This
  drastically reduces the number of generated PDFs and makes group-level trends
  obvious at a glance.
* Groups are inferred from sample names by stripping trailing digits
  (e.g. SAMPLE33, SAMPLE34 -> group SAMPLE; CONTROL01, TREATED01 -> CONTROL,
  TREATED).
* New `detail_samples` argument on those three functions restores detailed
  per-sample panels on demand: `NULL` (default) draws only the group overlay,
  `"all"` renders every sample, and a comma-separated string such as
  `"SAMPLE33, SAMPLE45"` renders just those. Detailed PDFs are written to a
  `per_sample/` subfolder.

## Bug Fixes

* **Confidence Score = 1.0 is now handled correctly.** Filenames using the
  natural `CS10` (or `CS100`, or the decimal `CS1.0`) for maximum stringency are
  now parsed as 1.0 (100%) instead of being silently misread as 0.1. CS values
  are stored internally on a single canonical integer-percent scale (0–100).
* `retrieve_selected_taxa()` and `heatmaps_karioCaS()` now accept a CS supplied
  either as a Kraken fraction (e.g. `1.0`) or a percentage (e.g. `40`); a `1.0`
  request previously produced an empty result.
* Silenced spurious `max()`/`-Inf` warnings emitted when a domain has no taxa at
  a given Confidence Score (common at high stringency).
* Removed a leftover hard-coded genus name from the `taxa_resolution()` audit log.
* `import_karioCaS()` now aggregates duplicate taxonomy rows explicitly
  (`values_fn = sum`), avoiding the cryptic `pivot_wider()` "Can't convert
  `fill` <double> to <list>" error on repeated taxonomy entries.

# karioCaS 0.99.1

## New Features

* Initial Bioconductor submission.
* `import_karioCaS()`: Imports Kraken2 MPA-style reports into a `TreeSummarizedExperiment` object.
* `taxa_retention()`: Evaluates taxonomic retention across Confidence Scores.
* `upset_kariocas()`: Generates UpSet plots to identify persistent vs transient taxa.
* `reads_per_taxa()`: Saturation analysis of reads per taxon.
* `taxa_resolution()`: Parent-to-child taxonomic resolution analysis.
* `heatmaps_karioCaS()`: Relative abundance heatmaps with extinction patterns.
* `optimize_CS()`: Multi-strategy Stability Index engine to find optimal Confidence Scores.
* `retrieve_selected_taxa()`: Generates the final high-confidence biological mosaic.