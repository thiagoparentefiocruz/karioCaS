# karioCaS validation (paper benchmark)

This branch holds the benchmark used to validate karioCaS against **ground truth**.
It is intentionally kept off `main` so the submitted Bioconductor package stays
package-only.

## What we are testing

1. **Domain-specific beats global** — per-domain optimal CS + min-reads gives higher
   F1 than a single global threshold, especially for under-represented domains
   (Archaea / Eukaryota / Viruses).
2. **The elbow finds the right threshold** — karioCaS's chosen CS/min-reads land
   near the *oracle* optimum (the CS/cutoff that truly maximises F1, computable
   because we know the answer).

## Pipeline (run order)

| Step | Script | What it does |
|---|---|---|
| 1 | `select_genomes.sh` | Download a multi-domain genome set (Bacteria/Archaea/Eukaryota/Viruses) + human, from NCBI RefSeq via the `datasets` CLI. |
| 2 | `build_community.py` | Turn the genomes into CAMISIM inputs (genome_to_id, metadata) + abundance designs (even / log-normal / staggered, with a rare noise tail + human contaminant). Emits `truth_species.tsv` / `truth_genus.tsv`. |
| 3 | `camisim_config.ini` + `run_camisim.sbatch` | Simulate Illumina paired reads with known composition (CAMISIM → ART). |
| 4 | `kraken_cs_sweep.sbatch` | Kraken2 CS sweep (0.0→1.0) on the **complete RefSeq** DB → karioCaS-ready `SAMPLE_CS<x>.mpa` in each project's `000_mpa_original/`. |
| 5 | run karioCaS | `import_karioCaS()` → `taxa_retention()` → `reads_per_taxa()` (produces the SI_Audit + Reads_Audit). |
| 6 | `score_karioCaS.R` | Per-domain Precision/Recall/F1 across the CS×cutoff grid; oracle-vs-karioCaS; baselines. Writes result tables + figures. |

## Mock design (controls the experiment)

* ~30 Bacteria, 5 Archaea, 5 Eukaryota, 5 Viruses — all present in RefSeq, so every
  miss is a true stringency cost (FN), not a database gap.
* Abundance profiles: **even** and **log-normal/staggered** (the staggered one
  stresses the rare-taxon / min-reads claim).
* **False-positive sources** baked in: a low-abundance noise tail (taxa at a few
  reads) and **human** contamination. Near-neighbour mis-assignment arises
  naturally at low CS because RefSeq is comprehensive.
* Replicate across 2–3 sequencing depths. Name samples cleanly (`MOCKeven01`,
  `MOCKstag01`, …) so karioCaS group parsing works (group = name minus trailing
  digits).

## Scoring definitions (see `score_karioCaS.R`)

* **Truth** = the species/genus names of the genomes we put in (per domain).
* **Detected** at a given (CS, read-cutoff) = taxa at that rank with count ≥ cutoff.
* **TP** = detected ∩ truth, **FP** = detected − truth, **FN** = truth − detected.
* **Precision** TP/(TP+FP), **Recall** TP/(TP+FN), **F1**, reported per domain.
* **Oracle** = the (CS, cutoff) maximising F1. **karioCaS pick** = Primary_SI from
  the SI_Audit (CS) and Reads_Audit (cutoff). We report the gap.

## Baselines compared
1. Raw Kraken2, CS 0, no filter (naive).
2. Fixed global CS = 0.1, no read filter.
3. Single **global** elbow-CS (ablation: isolates the domain-specific gain).
4. Read filter only (≥10), no CS optimisation.
5. KrakenUniq unique-k-mer filter (most direct competitor).

## Phase 2 — dark matter / discovery (real data)
Run karioCaS on a **published** environmental shotgun dataset (with SRA accessions
so reviewers can reproduce) to show the retention-curve behaviour on genuine
"biological dark matter". Optionally complement with the maintainer's own
environmental data as a real-world case study.
