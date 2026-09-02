# =============================================================================
# phase5_01_combine_args.R
# Phase 5 — Step 1: Load and combine ARG results from RGI, Abricate, AMRFinder
#
# INPUT:
#   amr_results/rgi/SAMPLE/rgi_output.txt
#   amr_results/abricate/SAMPLE/abricate_resfinder.txt
#   amr_results/abricate/SAMPLE/abricate_card.txt
#   amr_results/amrfinder/SAMPLE/amrfinder_results.txt
#
# OUTPUT:
#   phase5/args_combined.tsv
#   Columns: sample | contig | start | end | strand | gene | drug_class |
#            mechanism | tool
# =============================================================================

library(tidyverse)

BASE    <- "/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
OUTDIR  <- file.path(BASE, "phase5")
dir.create(OUTDIR, showWarnings = FALSE)

samples <- readLines(file.path(BASE, "r_analysis/decontam/sample_list_real.txt"))
cat("Samples to process:", length(samples), "\n\n")

# ── Helper: clean contig name (remove trailing hit info) ─────────────────────
clean_contig <- function(x) str_extract(x, "^[^\\s#]+")

# ── 1. RGI ────────────────────────────────────────────────────────────────────
# Uses Strict + Perfect hits only (Loose = too many false positives)
load_rgi <- function(sample) {
  f <- file.path(BASE, "amr_results/rgi", sample, "rgi_output.txt")
  if (!file.exists(f)) return(NULL)
  tryCatch({
    read_tsv(f, show_col_types = FALSE) %>%
      filter(Cut_Off %in% c("Perfect", "Strict")) %>%
      transmute(
        sample     = sample,
        contig     = clean_contig(Contig),
        start      = as.integer(Start),
        end        = as.integer(Stop),
        strand     = Orientation,
        gene       = `Best_Hit_ARO`,
        drug_class = `Drug Class`,
        mechanism  = `Resistance Mechanism`,
        identity   = as.numeric(`Percentage Length of Reference Sequence`),
        tool       = "RGI"
      )
  }, error = function(e) { message("  RGI error for ", sample, ": ", e$message); NULL })
}

# ── 2. Abricate ───────────────────────────────────────────────────────────────
# Run for resfinder and card databases; filter on identity + coverage
load_abricate <- function(sample, db) {
  f <- file.path(BASE, "amr_results/abricate", sample, paste0("abricate_", db, ".txt"))
  if (!file.exists(f)) return(NULL)
  tryCatch({
    d <- read_tsv(f, show_col_types = FALSE, comment = "#")
    if (nrow(d) == 0) return(NULL)
    d %>%
      filter(`%IDENTITY` >= 90, `%COVERAGE` >= 80) %>%
      transmute(
        sample     = sample,
        contig     = clean_contig(SEQUENCE),
        start      = as.integer(START),
        end        = as.integer(END),
        strand     = STRAND,
        gene       = GENE,
        drug_class = RESISTANCE,
        mechanism  = NA_character_,
        identity   = as.numeric(`%IDENTITY`),
        tool       = paste0("Abricate_", db)
      )
  }, error = function(e) { message("  Abricate error for ", sample, "/", db, ": ", e$message); NULL })
}

# ── 3. AMRFinder ─────────────────────────────────────────────────────────────
load_amrfinder <- function(sample) {
  f <- file.path(BASE, "amr_results/amrfinder", sample, "amrfinder_results.txt")
  if (!file.exists(f)) return(NULL)
  tryCatch({
    d <- read_tsv(f, show_col_types = FALSE)
    if (nrow(d) == 0) return(NULL)
    d %>%
      transmute(
        sample     = sample,
        contig     = clean_contig(`Contig id`),
        start      = as.integer(Start),
        end        = as.integer(Stop),
        strand     = Strand,
        gene       = `Element symbol`,
        drug_class = `Element name`,
        mechanism  = NA_character_,
        identity   = NA_real_,
        tool       = "AMRFinder"
      ) %>%
      filter(!is.na(contig), !is.na(start))
  }, error = function(e) { message("  AMRFinder error for ", sample, ": ", e$message); NULL })
}

# ── Combine all ───────────────────────────────────────────────────────────────
cat("Loading RGI results...\n")
rgi_all  <- map_dfr(samples, load_rgi)
cat("  ", nrow(rgi_all), "Strict/Perfect RGI hits\n")

cat("Loading Abricate (ResFinder)...\n")
ab_res   <- map_dfr(samples, load_abricate, db = "resfinder")
cat("  ", nrow(ab_res), "hits\n")

cat("Loading Abricate (CARD)...\n")
ab_card  <- map_dfr(samples, load_abricate, db = "card")
cat("  ", nrow(ab_card), "hits\n")

cat("Loading AMRFinder...\n")
amrf_all <- map_dfr(samples, load_amrfinder)
cat("  ", nrow(amrf_all), "hits\n")

args_combined <- bind_rows(rgi_all, ab_res, ab_card, amrf_all) %>%
  filter(!is.na(contig), !is.na(start), !is.na(end)) %>%
  # Deduplicate: same sample + contig + gene + approximate position → keep RGI > Abricate > AMRFinder
  mutate(tool_rank = case_when(
    tool == "RGI"              ~ 1,
    str_starts(tool, "Abricate") ~ 2,
    tool == "AMRFinder"        ~ 3,
    TRUE                       ~ 4
  )) %>%
  arrange(sample, contig, gene, tool_rank) %>%
  group_by(sample, contig, gene) %>%
  slice(1) %>%
  ungroup() %>%
  select(-tool_rank)

# Save
outfile <- file.path(OUTDIR, "args_combined.tsv")
write_tsv(args_combined, outfile)

cat("\n======================================\n")
cat("ARGs combined and saved\n")
cat("  Total ARG hits:      ", nrow(args_combined), "\n")
cat("  Unique genes:        ", n_distinct(args_combined$gene), "\n")
cat("  Unique drug classes: ", n_distinct(args_combined$drug_class), "\n")
cat("  Samples with ARGs:   ", n_distinct(args_combined$sample), "\n")
cat("  By tool:\n")
print(count(args_combined, tool) %>% arrange(desc(n)))
cat("  Output:", outfile, "\n")
cat("======================================\n")
