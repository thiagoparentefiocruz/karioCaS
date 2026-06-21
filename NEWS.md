# karioCaS 0.99.10

## Changes

* `upset_kariocas()` now takes a single `tax_level` argument (default
  `"Species"`), consistent with the other functions' rank flags, and draws one
  UpSet plot per sample and domain at that rank. Previously it produced a fixed
  set of three ranks (Species/Genus/Family) per sample/domain; the new behaviour
  reduces output clutter and lets you pick any rank.

# karioCaS 0.99.9

## Changes

* `retrieve_selected_taxa()` now **always retains all taxonomic ranks** in the
  mosaic (an MPA profile naturally has every level). The `tax_level` argument no
  longer filters the output; it now only selects which rank's optimization audit
  the `"auto"`/`"secondary"` thresholds are read from (`SI_Audit_<tax_level>` /
  `Reads_Audit_<tax_level>`; `NULL` -> `"Species"`). This also makes
  `taxa_resolution()` on the final mosaic reliable without any special setup.

# karioCaS 0.99.8

## Changes

* Applied `styler` (4-space indentation, Bioconductor style) across the package
  for consistent formatting.
* Added funding (`fnd`) roles to `Authors@R`: Fiocruz, IOC-Fiocruz, and CAPES.

# karioCaS 0.99.7

## Changes

* **`taxa_resolution()` no longer draws a plot for every Confidence Score.** It
  gains a `CS` argument: by default (`CS = NULL`) it analyses the **final mosaic**
  from `retrieve_selected_taxa()` (`1000_final_selection/`) - one figure per
  sample - importing and parsing the mosaic `.tsv` directly. Passing a numeric
  `CS` (fraction or percent) analyses the imported data at that single score
  instead of looping over all of them. Output filenames now end in
  `_Final_Mosaic` or `_CS<nn>`. For a meaningful mosaic resolution, build the
  mosaic with `retrieve_selected_taxa(tax_level = NULL)` so parent ranks are kept.

# karioCaS 0.99.6

## Changes

* **`retrieve_selected_taxa()` can now choose the minimum reads automatically.**
  The `reads_min_*` arguments accept `"auto"` / `"secondary"` (in addition to a
  manual number), pulling the optimal minimum reads from the `Reads_Audit`
  written by `reads_per_taxa()`, looked up at each domain's resolved Confidence
  Score. The final mosaic therefore combines both data-driven thresholds - the
  optimal CS and the optimal min-reads - per domain, closing the workflow loop.

# karioCaS 0.99.5

## Changes

* **`reads_per_taxa()` now reports an optimal minimum-reads threshold.** Using
  the same elbow engine as the optimal CS (default `"kneedle"`, on the log read
  axis), it finds the knee of each domain's saturation curve - the read count
  above which the stable taxa core persists and below which the rare/background
  tail is shed. The group overlay marks each domain's median optimal reads with a
  dashed line, and per-sample values are written to
  `Reads_Audit_<rank>.tsv`/`.rds`. Together with the optimal CS (Step 001) this
  gives two quantitative thresholds for excluding background false positives.
* **The `"Rare_Taxa"` view and the `x_max_*` arguments are removed.** That linear
  zoom of the 1-10 read region duplicated the low-count end of the saturation
  curve, which already covers it on the log axis. `reads_per_taxa()` now produces
  a single saturation plot per CS and gains a `method=` argument.

# karioCaS 0.99.4

## Changes

* **`optimize_CS()` has been merged into `taxa_retention()` and removed.** The
  two functions produced a near-identical group overlay; `taxa_retention()` now
  computes the Stability Index in the same step, marks each domain's median
  optimal CS on its overlay, writes the `SI_Audit_<rank>.tsv`/`.rds` tables, and
  returns the audit data frame invisibly. It gains `method=` (default
  `"kneedle"`) and `manual_toll=` arguments.
* The SI audit now lives in `001_taxa_retention/`. `retrieve_selected_taxa()`
  reads it from there, with a backward-compatible fallback to the old
  `006_optimize_CS/` location for existing projects.

# karioCaS 0.99.3

## Changes

* **`optimize_CS()` now defaults to a new `"kneedle"` method** (parameter-free
  elbow detection) and adds a `"postcliff"` method. The previous default,
  `"dynamic"`, picked the *first* CS whose step-wise taxa loss fell within tail
  noise; on domains whose retention curve starts with a plateau (Archaea,
  Eukaryota, Viruses in typical data) this stopped *before* the main drop,
  giving an inconsistent optimum (e.g. CS10 for those domains but CS60 for
  Bacteria). `"kneedle"` locates the inflection between the steep
  noise-removal phase and the stable signal floor consistently across domains.
  `"dynamic"`, `"segmented"` and `"manual"` remain available via `method =`.
* Under `"kneedle"`, the Secondary Stability Index is the more conservative
  post-cliff floor, so `retrieve_selected_taxa(CS_* = "secondary")` yields a
  stricter threshold than the Primary.

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