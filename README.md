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
3. **XX:** The Confidence Score (e.g., `00` for 0.0, `05` for 0.05, `90` for 0.9).

*Quick note on CS values: A CS = 0.1 means at least 10% of the k-mers from a read must map identically to a genome for assignment. A CS = 0.0 assigns a read if even a single k-mer maps (highest sensitivity, highest noise).*

## 🛠️ Key Features & Workflow

The package follows a logical, step-by-step workflow for metagenomic validation, automatically organizing outputs into dynamically generated subfolders:

**1. Data Harmonization (Step 000)**
* `import_karioCaS()`: Standardizes raw outputs from multiple stringencies into a cohesive Bioconductor `TreeSummarizedExperiment` (TSE) object.

**2. Visual Exploration (Steps 001 - 005)**
* `taxa_retention()`: Quantifies the percentage of taxa and reads that remain as confidence stringency increases.
* `taxa_resolution()`: Evaluates taxonomic depth and Parent-to-Child resolution across different scores.
* `reads_per_taxa()`: Visualizes read distribution efficiency across ranks.
* `upset_kariocas()`: Identifies "transient" vs. "persistent" taxa.
* `heatmaps_karioCaS()`: Detailed abundance heatmaps showing taxa extinction patterns.

**3. Objective Thresholding (Step 006)**
* `optimize_CS()`: The core mathematical engine. Replaces subjective thresholding with a Multi-Strategy Stability Index (SI):
  * **Dynamic:** Adapts to the basal noise of the sample (ideal for pathogen focus).
  * **Segmented:** Uses a Broken-Stick regression model to find regime shifts (ideal for ecology/dark matter).
  * **Manual:** Allows expert-defined acceptable loss tolls.

**4. The Ultimate Biological Mosaic (Step 1000)**
* `retrieve_selected_taxa()`: Seamlessly integrates with the `optimize_CS` audit file to extract surviving taxa based on domain-specific mathematical thresholds. Generates a highly refined, high-confidence `.mpa` file, scrubbed of statistical noise and ready for downstream analysis.

## 📖 Quick Example

```r
library(karioCaS)

proj_dir <- "path/to/your_project"

# 1. Harmonize Data
import_karioCaS(project_dir = proj_dir)

# 2. Visual Exploration 

    # Evaluate Taxonomic Retention
      taxa_retention(project_dir = proj_dir)

    # Evaluate the taxa "persistance" over increasing Confidence Scores 
      upset_kariocas(project_dir = proj_dir)

    # Evaluate NGS Read Retention
      reads_per_taxa(project_dir = proj_dir)

    # Evaluate Parent-to-Child taxa resolution
      taxa_resolution(project_dir = proj_dir)

    # Evaluate taxa extinction patterns
      heatmaps_karioCaS(project_dir = proj_dir)

# 3. Optimize Confidence Scores (Objective Mathematics)
# Calculates the Stability Index automatically
optimize_CS(project_dir = proj_dir, tax_level = "Species", method = "dynamic")

# 4. Retrieve the Final Biological Mosaic
# Uses mathematically optimal thresholds for Bacteria and Archaea, 
# and manual overrides for Eukaryota and Viruses.
retrieve_selected_taxa(
  project_dir = proj_dir, 
  tax_level = "Species",
  CS_B = "auto", 
  CS_A = "auto", 
  CS_E = 40, 
  CS_V = 0
)
```
