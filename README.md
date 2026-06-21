# karioCaS 🧬

[![DOI](https://zenodo.org/badge/1141170581.svg)](https://doi.org/10.5281/zenodo.19392342)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

**Kraken Confidence Scores for Reliable Domain-Specific Microbiota Inference and Discovery**

> Turn the noisy output of k-mer classifiers into a defensible, high-confidence microbial profile — with thresholds chosen by the data, not by guesswork.

## ✨ Why karioCaS?

k-mer classifiers such as Kraken2 are blazing fast but notoriously **hyper-inflate false positives** — thousands of taxa end up with a handful of reads each. The usual remedy, a single global Confidence Score (CS) for the whole sample, is a blunt instrument: **Bacteria, Archaea, Eukaryota and Viruses** differ wildly in database representation, genome size and k-mer uniqueness. **One size does not fit all.**

karioCaS makes the cleanup **objective, reproducible, and domain-specific** by deriving *two* data-driven thresholds for every domain:

* 🎯 **Optimal Confidence Score** — the exact mathematical inflection point where extra stringency stops removing statistical noise and starts deleting true biological signal (the *Biological Dark Matter*).
* 🔢 **Optimal minimum number of reads** — the saturation-curve elbow below which taxa are indistinguishable from background.

It then assembles a **biological mosaic**: the taxa that survive each domain's own optimal CS *and* read cutoff — a refined, high-confidence profile ready for downstream ecological or clinical analysis. Every step also produces **publication-ready, per-group figures**.

* **Objective thresholds** — no more eyeballing a cutoff; the numbers come from the decay/saturation curves.
* **Domain-aware** — each domain is analyzed separately, so rare domains aren't drowned out by Bacteria.
* **Group-aware** — samples from the same biological group are overlaid on one figure, and **core** vs **unique** taxa are surfaced (the expected pathogen / false-positive signature).
* **Bioconductor-native** — your data lives in a `TreeSummarizedExperiment`.

## 🚀 Installation

```r
# install.packages("devtools")
devtools::install_github("thiagoparentefiocruz/karioCaS")
library(karioCaS)
```
*(karioCaS is currently under preparation for Bioconductor submission.)*

## 📁 Folder Architecture

karioCaS reads one input folder and writes everything else for you. Place your Kraken2 / MetaPhlAn-style reports inside a folder named **`000_mpa_original`** at your project root — one file per sample × Confidence Score:

```text
YOUR_PROJECT_DIR/
 └── 000_mpa_original/
      ├── SAMPLE01_CS00.mpa      # CS 0.0  (most permissive)
      ├── SAMPLE01_CS50.mpa      # CS 0.5
      └── SAMPLE01_CS100.mpa     # CS 1.0  (most stringent)
```

**File naming rules**

1. **Sample name:** anything you like, but **no underscores** in the sample name itself.
2. **`_CS` separator:** must be present.
3. **Confidence Score:** written as a **percentage** (recommended) — `CS00` = 0.0, `CS40` = 0.4, `CS90` = 0.9, `CS100` = 1.0 — or as an **explicit decimal** (`CS0.0`, `CS0.05`, `CS0.9`, `CS1.0`).

> A CS of 0.1 means at least 10% of a read's k-mers must map to the same genome; CS 0.0 keeps a read if even a single k-mer maps (maximum sensitivity, maximum noise), while CS 1.0 requires every k-mer to agree (maximum stringency). Samples are grouped automatically by name prefix (`SAMPLE01`, `SAMPLE02` → group `SAMPLE`; `CONTROL01`, `TREATED01` → `CONTROL`, `TREATED`).

Function arguments such as `confidence_score=` and `CS_B=` accept either a Kraken fraction in `(0, 1]` (e.g. `1.0` for maximum stringency) or a percentage `> 1` (e.g. `40`).

## 🛠️ The Workflow

Outputs are organized into numbered subfolders that follow the analysis logic — **harmonize → quantify the optima → decide → describe → interpret**:

**1. Harmonize → `001_imported_matrix`**
* `import_karioCaS()` — standardizes the multi-stringency reports into one Bioconductor `TreeSummarizedExperiment`, the dataset every other function reads.

**2. Quantify the optima → `002_taxa_retention`, `003_reads_saturation`**
* `taxa_retention()` — taxa-retention decay curves **and** the optimal CS per domain (Stability Index). Marks each domain's median optimal CS on the group plot and writes the `SI_Audit` tables. Strategies (`method=`): **Kneedle** (default, parameter-free elbow), **Post-cliff** (conservative), **Segmented** (broken-stick), **Dynamic / Manual**.
* `reads_per_taxa()` — saturation analysis on a log read axis **and** the optimal minimum reads per domain (same elbow engine). Writes the `Reads_Audit` tables.

**3. Decide → `004_final_mosaic`**
* `retrieve_selected_taxa()` — builds the biological mosaic using **both** thresholds per domain (optimal CS and optimal min-reads, looked up together). `CS_*` and `reads_min_*` each accept `"auto"`, `"secondary"`, or a manual value. Writes `mpa/` and `tsv/` profiles, keeping **all** taxonomic ranks.

**4. Describe → `005_taxa_intersections_across_CS`, `006_relative_abundance_across_CS`**
* `upset_kariocas()` — UpSet of taxon persistence **across Confidence Scores** (per sample/domain).
* `heatmaps_karioCaS()` — relative-abundance heatmaps showing taxa extinction patterns across CS.

**5. Interpret → `007_taxa_resolution`, `008_taxa_intersections_across_samples`**
* `taxa_resolution()` — Parent-to-Child resolution of the final mosaic.
* `group_upset()` — cross-sample UpSet within each biological group; separates **core** taxa (in every sample) from **unique**/rare taxa, and writes a Core/Shared/Unique membership table.

> **Group overlays by default.** `taxa_retention()` and `reads_per_taxa()` draw one figure **per biological group** (each sample a faint line, the group mean ± SD highlighted) instead of a flood of per-sample PDFs. Add `detail_samples=` to also export per-sample panels.

## 📖 Quick Example

Every call below shows its full set of flags; the commented `#` lines are the defaults — uncomment to customize.

```r
library(karioCaS)
proj_dir <- "path/to/your_project"

# 1. Harmonize the multi-stringency reports into a TSE
import_karioCaS(project_dir = proj_dir)

# 2. Optimal Confidence Score per domain (+ retention overlay)
taxa_retention(
  project_dir = proj_dir
  # tax_level      = "Species",   # rank used for the optimal-CS curve
  # method         = "kneedle",   # "kneedle" | "postcliff" | "segmented" | "dynamic" | "manual"
  # manual_toll    = 1.0,         # acceptable step-loss %, only for method = "manual"
  # detail_samples = NULL         # NULL = group only | "all" | "SAMPLE01, SAMPLE02"
)

# 3. Optimal minimum reads per domain (+ saturation overlay)
reads_per_taxa(
  project_dir = proj_dir
  # analysis_level = "Species",
  # method         = "kneedle",   # "kneedle" | "postcliff" | "segmented"
  # detail_samples = NULL
)

# 4. Build the final biological mosaic.
#    CS_* and reads_min_* each take "auto" | "secondary" | a number.
retrieve_selected_taxa(
  project_dir = proj_dir,
  CS_B = "auto", reads_min_B = "auto",   # Bacteria  : fully data-driven
  CS_A = "auto", reads_min_A = "auto",   # Archaea   : fully data-driven
  CS_E = 40,     reads_min_E = 10,       # Eukaryota : manual overrides
  CS_V = 0,      reads_min_V = 0         # Viruses   : keep everything
  # tax_level = NULL                     # which rank's audit "auto" reads (NULL = Species)
)

# --- Descriptive diagnostics (optional, run any time after step 1) ---
upset_kariocas(project_dir = proj_dir)          # tax_level = "Species"
heatmaps_karioCaS(project_dir = proj_dir)        # analysis_rank = "Genus", confidence_score = NULL, top_n = 20

# --- Biological insight from the final mosaic ---
taxa_resolution(
  project_dir = proj_dir
  # parent_level = "Genus", child_level = "Species",
  # CS           = NULL,    # NULL = final mosaic | a number = a single CS
  # top_n        = 10
)

group_upset(
  project_dir = proj_dir
  # tax_level = "Species",
  # CS        = NULL        # NULL = final mosaic | a number = a single CS
)
```
