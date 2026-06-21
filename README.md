# karioCaS 🧬

[![DOI](https://zenodo.org/badge/1141170581.svg)](https://doi.org/10.5281/zenodo.19392342)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

**Kraken Confidence Scores for Reliable Domain-Specific Microbiota Inference and Discovery**

**karioCaS** is an R package designed to enhance the reliability and clarity of metagenomic data analysis originating from Kraken2 (or any k-mer based classifier). It prevents the hyper-inflation of false-positive assignments by providing a robust analytical framework to evaluate taxonomic stability across multiple stringency levels.

**karioCaS** offers two distinct approaches to deep-diving into metagenomic data:

**1. Domain-Specific Analysis:** Based on the principle that "one size does **not** fit all," all analyses and visualizations are performed individually by Biological Domain (**Archaea, Bacteria, Eukaryota, and Viruses**).

**2. Confidence Score Exploration:** We leverage a crucial yet underutilized Kraken parameter—the **Confidence Score (CS)**. 

With **karioCaS**, users can objectively identify the optimal Confidence Score unique to each Domain and appropriate to their project goals, finding the exact mathematical inflection point to separate statistical noise from true biological signal (Biological Dark Matter).

## 🚀 Installation

You can install the development version of **karioCaS** directly from GitHub:

```r
install.packages("devtools")
devtools::install_github("thiagoparentefiocruz/karioCaS")
```
*(Note: karioCaS is currently under preparation for Bioconductor submission).*

## 📁 Strict Folder Architecture

**IMPORTANT:** The package requires a specific directory structure to function correctly. You must place your raw Metaphlan-style outputs inside a base directory, specifically within a folder named `000_mpa_original`.

```text
YOUR_PROJECT_DIR/
 └── 000_mpa_original/
      ├── SAMPLE01_CS00.mpa
      ├── SAMPLE01_CS05.mpa
      └── SAMPLE01_CS10.mpa
```

**File naming rules:**
1. **SAMPLE:** Any name (DO NOT use underscores `_` in the sample name itself).
2. **_CS:** **Must** be present.
3. **XX:** The Confidence Score. karioCaS accepts three equivalent notations:
   * **Percentage** (recommended): `CS00` = 0.0, `CS40` = 0.4, `CS90` = 0.9, `CS100` = 1.0.
   * **Legacy 0.1-step shorthand:** the zero-padded `CS00`–`CS10`, where `CS09` = 0.9 and `CS10` = 1.0. *(Note: `CS05` resolves to 0.5, not 0.05.)*
   * **Explicit decimal:** `CS0.0`, `CS0.05`, `CS0.9`, `CS1.0`.

*Quick note on CS values: A CS = 0.1 means at least 10% of the k-mers from a read must map identically to a genome for assignment. A CS = 0.0 assigns a read if even a single k-mer maps (highest sensitivity, highest noise); a CS = 1.0 requires every k-mer to agree (maximum stringency).*

*Function arguments such as `confidence_score=` and `CS_B=` accept either a Kraken fraction in `(0, 1]` (e.g. `1.0` for maximum stringency) or a percentage `> 1` (e.g. `40`).*

## 🛠️ Key Features & Workflow

The package follows a logical, step-by-step workflow for metagenomic validation, automatically organizing outputs into dynamically generated subfolders:

**1. Data Harmonization (Step 000)**
* `import_karioCaS()`: Standardizes raw outputs from multiple stringencies into a cohesive Bioconductor `TreeSummarizedExperiment` (TSE) object.

**2. Visual Exploration & Objective Thresholding (Steps 001 - 005)**
* `taxa_retention()`: Quantifies the percentage of taxa retained as stringency increases **and** computes the objective optimal CS (Stability Index) per domain in the same step. Its group plot marks each domain's median optimal CS with a dashed line, and it writes the `SI_Audit_<rank>` tables consumed by Step 1000. Strategies (`method=`):
  * **Kneedle (default):** Parameter-free elbow detection — finds the inflection between the steep noise-removal phase and the stable signal floor, consistently across domains.
  * **Post-cliff:** A more conservative threshold deeper in the plateau (first stable CS after the steepest drop).
  * **Segmented:** Broken-stick regression for regime shifts (ideal for ecology/dark matter).
  * **Dynamic / Manual:** First CS within tail noise / expert-defined loss tolls.
* `taxa_resolution()`: Evaluates taxonomic depth (Parent-to-Child resolution). By default it analyzes the **final mosaic** (`retrieve_selected_taxa()` output) — one figure per sample; pass `CS=` to instead analyze the imported data at a single Confidence Score.
* `group_upset()`: Cross-sample UpSet **within each biological group** (inferred from name prefixes) — separates **core** taxa (in every sample) from **unique**/rare taxa (in one or a few samples), the expected pattern for pathogens and false positives. Writes a membership TSV (presence matrix + Core/Shared/Unique category). Default source is the final mosaic; `CS=` compares at a single Confidence Score.
* `reads_per_taxa()`: Saturation analysis on a log read axis. Computes the **optimal minimum reads** per domain (the saturation-curve elbow, via the same engine as the optimal CS), marks each domain's median on the group plot, and writes `Reads_Audit_<rank>` tables — a quantitative threshold for excluding low-abundance background/false-positive taxa.
* `upset_kariocas()`: Identifies "transient" vs. "persistent" taxa.
* `heatmaps_karioCaS()`: Detailed abundance heatmaps showing taxa extinction patterns.

> **Group overlays by default.** `taxa_retention()` and `reads_per_taxa()` draw a single figure **per biological group** instead of one set of PDFs per sample: every sample is a faint line and the group mean (± SD) is highlighted, faceted by Domain. Groups are inferred from sample names by stripping trailing digits (e.g. `SAMPLE33`, `SAMPLE34` → group `SAMPLE`; `CONTROL01`, `TREATED01` → `CONTROL`, `TREATED`). Use `detail_samples=` to also render detailed per-sample panels into a `per_sample/` subfolder: `NULL` (default) = group only, `"all"` = every sample, or a comma-separated list such as `"SAMPLE33, SAMPLE45"`.

**3. The Ultimate Biological Mosaic (Step 1000)**
* `retrieve_selected_taxa()`: Extracts surviving taxa per domain using **both** data-driven thresholds — the optimal CS (`SI_Audit` from `taxa_retention()`) and the optimal minimum reads (`Reads_Audit` from `reads_per_taxa()`), looked up at the resolved CS. `CS_*` and `reads_min_*` each accept `"auto"`, `"secondary"`, or a manual value. Generates a highly refined, high-confidence `.mpa` file, scrubbed of statistical noise and ready for downstream analysis.

## 📖 Quick Example

```r
library(karioCaS)

proj_dir <- "path/to/your_project"

# 1. Harmonize Data
import_karioCaS(project_dir = proj_dir)

# 2. Visual Exploration + Objective Thresholding

    # Taxa retention overlay + optimal CS (Stability Index) per group — one
    # figure + the SI_Audit tables. Default method = Kneedle elbow.
      taxa_retention(project_dir = proj_dir)

    # ...drill into specific samples (detailed PDFs in per_sample/)
      taxa_retention(project_dir = proj_dir, detail_samples = "SAMPLE33, SAMPLE45")

    # Evaluate the taxa "persistance" over increasing Confidence Scores 
      upset_kariocas(project_dir = proj_dir)

    # Evaluate NGS Read Retention (group overlay; detail_samples = "all" for every sample)
      reads_per_taxa(project_dir = proj_dir)

    # Parent-to-Child taxa resolution at a chosen Confidence Score
      taxa_resolution(project_dir = proj_dir, CS = 40)

    # Evaluate taxa extinction patterns
      heatmaps_karioCaS(project_dir = proj_dir)

# 3. Retrieve the Final Biological Mosaic
# Both CS_* and reads_min_* accept "auto"/"secondary" (pulled from the SI and
# Reads audits) or a manual value. Here: fully data-driven for Bacteria/Archaea,
# manual overrides for Eukaryota/Viruses.
retrieve_selected_taxa(
  project_dir = proj_dir, 
  tax_level = "Species",
  CS_B = "auto", reads_min_B = "auto",
  CS_A = "auto", reads_min_A = "auto",
  CS_E = 40,     reads_min_E = 10,
  CS_V = 0,      reads_min_V = 0
)

# 4. Inspect Parent-to-Child resolution of the FINAL mosaic (default source)
#    (the mosaic always keeps all ranks, so resolution works out of the box)
taxa_resolution(project_dir = proj_dir)

# 5. Core vs unique taxa across samples within each biological group
group_upset(project_dir = proj_dir)
```
