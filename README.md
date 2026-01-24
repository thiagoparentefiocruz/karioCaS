# karioCaS 🧬

**Kraken Confidence Scores for Reliable Domain-Specific Microbiota
Inference and Discovery**

**karioCaS** is an R package designed to enhance the reliability and
clarity of metagenomic data analysis originating from Kraken2.
**karioCaS** offers two distinct approaches to deep-diving into
metagenomic data:

**1. Domain-Specific Analysis:** Based on the principle that "one size
does **not** fit all," all analyses and visualizations are performed
individually by Biological Domain (**Archaea, Bacteria, Eukaryota, and
Viruses**).

**2. Confidence Score Exploration:** We leverage a crucial yet
underutilized Kraken parameter—the **Confidence Score (CS)**.

The comparative analysis of taxa retention across different ranks is
astonishingly informative. With **karioCaS**, users can identify
Confidence Scores unique to each Domain and appropriate to the project
goals, finding the fine balance between false positives and false
negatives. True-positive calling strength can be further increased by
setting a minimal read threshold for each Domain and taxonomic rank.

Ultimately, **karioCaS** provides a chimeric `.mpa` file listing all
taxa ranks confidently detected in a sample, providing a robust
framework to evaluate taxonomic stability across multiple stringency
levels.

## 🚀 Installation

Since this is currently a private repository, install it using a
Personal Access Token (PAT):

r install.packages("devtools")
devtools::install_github("thiagoparentefiocruz/karioCaS", auth_token =
"YOUR_PERSONAL_ACCESS_TOKEN")

**Strict Folder Architecture**

IMPORTANT: The package requires a specific directory structure to
function correctly.

YOUR_BASE_DIR/000_mpa_original/ SAMPLE_CSXX/\
\
Folder name rules:\
1. SAMPLE: Any name (DO NOT use underscores "\_")\
\
2. \_CS: **Must** be present\
\
3. XX: The Confidence Score (e.g., 01 for 0.1, 09 for 0.9)

*Quick note on CS values: A CS = 0.1 means at least 10% of the kmers
from a read must map identically to a genome for assignment. A CS = 0.0
assigns a read if even a single kmer maps.*

**Key Features**

The package follows a logical workflow for metagenomic validation:

1.  Data Harmonization

**import_mpa_data()**: Standardizes raw Kraken2/MPA outputs from the
same sample run with different Confidence Scores.

Note: We are currently optimizing data importation to offer deeper
functionalities in future versions.

2.  Stability Metrics

**taxa_retention()**: Quantifies the percentage of the sample and reads
that remain as confidence stringency increases.

**taxa_resolution()**: Evaluates the taxonomic depth achieved across
different scores.

3.  Advanced Visualization

**heatmaps_karioCaS()**: Detailed abundance heatmaps with built-in
support for confidence clustering.

**upset_karioCaS_never_are()**: Identifies "transient" vs. "persistent"
taxa using UpSet plots.

**reads_per_taxa()**: Visualizes read distribution efficiency across
ranks.

4.  The ultimate Output

**retrieve_selected_taxa()**: saves the list of all biological ranks
confidently identified in the sample, filtering for Domain, Taxonomic
Rank, Confidence Score and Minimal Number of Reads. Files are saved as
'.tsv', for easy conference/audit, an also as '.mpa', ready for further
downstream analysis.

***Quick Example***

library(karioCaS)

**1. Import MPA (Metaphlan-style) data from Kraken2**

This harmonizes the different confidence score outputs into one object

`my_data <- import_mpa_data("your_project_name")`

**2. Evaluate Taxonomic Retention**

Visualize how many taxa "survive" as you increase stringency

`taxa_retention(my_data)`

**3. Assess Taxonomic Resolution**

Identify the lowest common ancestor (LCA) stability

`taxa_resolution(my_data)`

**4. Generate the Heatmap**

Compare abundances across samples and confidence levels

`heatmaps_karioCaS(my_data)`
