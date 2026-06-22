#!/usr/bin/env Rscript
# Score karioCaS against the mock ground truth.
#
# Computes per-domain Precision/Recall/F1 across the CS x read-cutoff grid, then:
#   (1) the ORACLE optimum (CS, cutoff maximising F1, per sample/domain), and
#   (2) karioCaS's PICK (Primary SI from SI_Audit + Reads_Audit), and the gap;
#   (3) baselines (raw, fixed-CS, global-elbow, read-filter-only).
#
# Usage (adapt paths):
#   Rscript validation/score_karioCaS.R \
#     --project   /path/kariocas_project \
#     --truth_sp  camisim_in/truth_species.tsv \
#     --truth_ge  camisim_in/truth_genus.tsv \
#     --out       validation_results
#
# Detection model (matches karioCaS): from the mpa-style reports, a species is
# "detected" at a given (CS, cutoff) if its s__ line has count >= cutoff; a genus
# if its g__ (non-s__) line does. Truth = the taxa we put in the mock.

suppressMessages({
    library(dplyr); library(tidyr); library(readr); library(ggplot2); library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) {
    i <- which(args == flag)
    if (length(i) == 1 && i < length(args)) args[i + 1] else default
}
project  <- getopt("--project", ".")
mpa_dir  <- getopt("--mpa_dir", file.path(project, "000_mpa_original"))
truth_sp <- getopt("--truth_sp", "camisim_in/truth_species.tsv")
truth_ge <- getopt("--truth_ge", "camisim_in/truth_genus.tsv")
out_dir  <- getopt("--out", "validation_results")
cutoffs  <- as.numeric(strsplit(getopt("--cutoffs", "1,2,3,5,10,20,50,100"), ",")[[1]])
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

norm_name <- function(x) str_squish(tolower(x))

# ---- parse mpa reports -> long detection table -----------------------------
parse_mpa <- function(dir) {
    files <- list.files(dir, pattern = "_CS[0-9.]+\\.(mpa|txt|report)$", full.names = TRUE)
    if (!length(files)) stop("No _CS*.mpa files in ", dir)
    do.call(rbind, lapply(files, function(f) {
        b <- basename(f)
        m <- str_match(b, "^(.*)_CS([0-9.]+)\\.(mpa|txt|report)$")
        samp <- m[2]; tok <- as.numeric(m[3])
        cs <- if (grepl("\\.", m[3]) || (tok > 0 && tok < 1)) tok * 100 else tok
        d <- tryCatch(read.delim(f, header = FALSE, sep = "\t",
                                 quote = "", stringsAsFactors = FALSE),
                      error = function(e) NULL)
        if (is.null(d) || ncol(d) < 2) return(NULL)
        lin <- d[[1]]; cnt <- as.numeric(d[[2]])
        domain  <- str_match(lin, "d__([^|]+)")[, 2]
        species <- str_match(lin, "\\|s__([^|]+)")[, 2]
        genus   <- str_match(lin, "\\|g__([^|]+)")[, 2]
        sp <- !is.na(species)
        ge <- !is.na(genus) & is.na(species)
        rbind(
            data.frame(sample = samp, cs = cs, domain = domain[sp],
                       rank = "Species", taxon = species[sp], count = cnt[sp],
                       stringsAsFactors = FALSE),
            data.frame(sample = samp, cs = cs, domain = domain[ge],
                       rank = "Genus", taxon = genus[ge], count = cnt[ge],
                       stringsAsFactors = FALSE)
        )
    }))
}

prf <- function(detected, truth) {
    detected <- unique(detected); truth <- unique(truth)
    tp <- length(intersect(detected, truth))
    fp <- length(setdiff(detected, truth))
    fn <- length(setdiff(truth, detected))
    prec <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
    rec  <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
    f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0)
        2 * prec * rec / (prec + rec) else 0
    c(TP = tp, FP = fp, FN = fn, Precision = prec, Recall = rec, F1 = f1)
}

det <- parse_mpa(mpa_dir)
det$taxon_n <- norm_name(det$taxon)

tsp <- read_tsv(truth_sp, show_col_types = FALSE) |>
    transmute(sample, domain, rank = "Species", taxon_n = norm_name(species))
tge <- read_tsv(truth_ge, show_col_types = FALSE) |>
    transmute(sample, domain, rank = "Genus", taxon_n = norm_name(genus))
truth <- bind_rows(tsp, tge)

# ---- full grid: per sample x rank x domain x CS x cutoff -------------------
grid <- expand_grid(
    distinct(det, sample, rank, domain, cs),
    cutoff = cutoffs
)
score_row <- function(sample, rank, domain, cs, cutoff) {
    d <- det$taxon_n[det$sample == sample & det$rank == rank &
                     det$domain == domain & det$cs == cs & det$count >= cutoff]
    t <- truth$taxon_n[truth$sample == sample & truth$rank == rank &
                       truth$domain == domain]
    if (length(t) == 0) return(NULL)
    as.list(prf(d, t))
}
res <- grid |>
    rowwise() |>
    mutate(m = list(score_row(sample, rank, domain, cs, cutoff))) |>
    filter(!is.null(m)) |>
    unnest_wider(m) |>
    ungroup()
write_tsv(res, file.path(out_dir, "metrics_grid.tsv"))

# ---- oracle (max-F1 over CS x cutoff) per sample/rank/domain ---------------
oracle <- res |>
    group_by(sample, rank, domain) |>
    slice_max(F1, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(sample, rank, domain, oracle_cs = cs, oracle_cutoff = cutoff,
              oracle_F1 = F1, oracle_Prec = Precision, oracle_Rec = Recall)

# ---- karioCaS pick (Primary SI: CS from SI_Audit, cutoff from Reads_Audit) --
read_audit <- function(folder, prefix) {
    f <- list.files(file.path(project, folder), pattern = paste0("^", prefix, ".*\\.rds$"),
                    full.names = TRUE)
    if (!length(f)) return(NULL)
    readRDS(f[1])
}
si <- read_audit("002_taxa_retention", "SI_Audit")
ra <- read_audit("003_reads_saturation", "Reads_Audit")

pick <- NULL
if (!is.null(si)) {
    si_pick <- si |> filter(.data$SI_Type == "Primary_SI") |>
        transmute(sample = Sample, domain = Domain, pick_cs = CS)
    pick <- si_pick
    if (!is.null(ra)) {
        ra_pick <- ra |> filter(.data$SI_Type == "Primary_SI") |>
            transmute(sample = Sample, domain = Domain, pick_cs = CS, pick_cutoff = Cutoff)
        pick <- left_join(si_pick, ra_pick, by = c("sample", "domain", "pick_cs"))
    }
    if (is.null(pick$pick_cutoff)) pick$pick_cutoff <- 1
    pick$pick_cutoff[is.na(pick$pick_cutoff)] <- 1
    # F1 at karioCaS's chosen (CS, cutoff) for Species (extend to Genus as needed)
    pick <- pick |> rowwise() |>
        mutate(rank = "Species",
               kcs = list(score_row(sample, "Species", domain, pick_cs, pick_cutoff))) |>
        filter(!is.null(kcs)) |> unnest_wider(kcs) |> ungroup() |>
        transmute(sample, rank, domain, pick_cs, pick_cutoff,
                  kcs_F1 = F1, kcs_Prec = Precision, kcs_Rec = Recall)
}

# ---- baselines (Species) ---------------------------------------------------
baseline <- function(label, cs_val, cutoff_val) {
    res |> filter(rank == "Species", cs == cs_val, cutoff == cutoff_val) |>
        group_by(domain) |>
        summarise(method = label, F1 = mean(F1, na.rm = TRUE),
                  Precision = mean(Precision, na.rm = TRUE),
                  Recall = mean(Recall, na.rm = TRUE), .groups = "drop")
}
bl <- bind_rows(
    baseline("raw_CS0_nofilter",  0,  1),
    baseline("fixedCS10_nofilter", 10, 1),
    baseline("CS0_reads10",        0, 10)
)
write_tsv(bl, file.path(out_dir, "baselines_species.tsv"))

# ---- compare oracle vs karioCaS (Species) ---------------------------------
if (!is.null(pick)) {
    cmp <- oracle |> filter(rank == "Species") |>
        left_join(pick, by = c("sample", "rank", "domain"))
    write_tsv(cmp, file.path(out_dir, "oracle_vs_kariocas_species.tsv"))
    summ <- cmp |> group_by(domain) |>
        summarise(median_oracle_cs = median(oracle_cs),
                  median_pick_cs = median(pick_cs, na.rm = TRUE),
                  mean_oracle_F1 = mean(oracle_F1, na.rm = TRUE),
                  mean_kariocas_F1 = mean(kcs_F1, na.rm = TRUE),
                  mean_F1_gap = mean(oracle_F1 - kcs_F1, na.rm = TRUE),
                  .groups = "drop")
    write_tsv(summ, file.path(out_dir, "summary_species.tsv"))
    print(summ)

    # F1 vs CS per domain (mean over samples), marking karioCaS median pick
    fdat <- res |> filter(rank == "Species", cutoff == 1) |>
        group_by(domain, cs) |> summarise(F1 = mean(F1, na.rm = TRUE), .groups = "drop")
    vmark <- cmp |> group_by(domain) |>
        summarise(pick = median(pick_cs, na.rm = TRUE),
                  orc = median(oracle_cs), .groups = "drop")
    p <- ggplot(fdat, aes(cs, F1)) +
        geom_line(linewidth = 1) + geom_point() +
        geom_vline(data = vmark, aes(xintercept = pick), linetype = "dashed", color = "red") +
        geom_vline(data = vmark, aes(xintercept = orc), linetype = "dotted", color = "blue") +
        facet_wrap(~domain, scales = "free_y") +
        labs(title = "F1 vs Confidence Score (Species)",
             subtitle = "red dashed = karioCaS median pick; blue dotted = oracle",
             x = "Confidence Score (%)", y = "F1") + theme_bw()
    ggsave(file.path(out_dir, "f1_vs_cs_species.pdf"), p, width = 10, height = 7)
}

cat("\nDone. Outputs in ", out_dir, "/\n", sep = "")
