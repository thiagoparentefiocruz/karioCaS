# karioCaS 0.99.2

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