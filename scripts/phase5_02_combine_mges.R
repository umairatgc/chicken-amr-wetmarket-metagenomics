# =============================================================================
# phase5_02_combine_mges.R
# Phase 5 — Step 2: Load and combine MGE results from all Phase 3 tools
#
# INPUT:
#   mobile_elements/SAMPLE/isescan/isescan/*.tsv       (IS elements)
#   mobile_elements/SAMPLE/integronfinder/Results_*/consensus.integrons
#   mobile_elements/SAMPLE/genomad/consensus_summary/consensus_plasmid_summary.tsv
#   mobile_elements/SAMPLE/mobsuite/contig_report.txt  (plasmids)
#
# OUTPUT:
#   phase5/mges_combined.tsv
#   Columns: sample | contig | start | end | mge_type | mge_name | tool
# =============================================================================

library(tidyverse)

BASE   <- "/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
OUTDIR <- file.path(BASE, "phase5")
dir.create(OUTDIR, showWarnings = FALSE)

samples <- readLines(file.path(BASE, "r_analysis/decontam/sample_list_real.txt"))

clean_contig <- function(x) str_extract(x, "^[^\\s#]+")

# ── 1. ISEScan — Insertion Sequences ─────────────────────────────────────────
# Columns: seqID family cluster isBegin isEnd isLen ...
load_isescan <- function(sample) {
  dir <- file.path(BASE, "mobile_elements", sample, "isescan/isescan")
  if (!dir.exists(dir)) return(NULL)
  files <- list.files(dir, pattern = "\\.tsv$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  tryCatch({
    map_dfr(files, ~read_tsv(.x, show_col_types = FALSE)) %>%
      transmute(
        sample   = sample,
        contig   = clean_contig(seqID),
        start    = as.integer(isBegin),
        end      = as.integer(isEnd),
        mge_type = "IS_element",
        mge_name = family,
        tool     = "ISEScan"
      ) %>%
      filter(!is.na(start), !is.na(end))
  }, error = function(e) { message("ISEScan error ", sample, ": ", e$message); NULL })
}

# ── 2. IntegronFinder — Integrons ────────────────────────────────────────────
# Columns: ID_integron ID_replicon element pos_beg pos_end strand evalue type_elt annotation model type
load_integrons <- function(sample) {
  f <- file.path(BASE, "mobile_elements", sample,
                 "integronfinder/Results_Integron_Finder_consensus/consensus.integrons")
  if (!file.exists(f)) return(NULL)
  tryCatch({
    d <- read_tsv(f, show_col_types = FALSE, comment = "#")
    if (is.null(d) || nrow(d) == 0) return(NULL)
    # Filter to integron-level rows (not individual gene cassettes)
    d %>%
      filter(type_elt == "integron" | type == "complete" | type == "In0") %>%
      transmute(
        sample   = sample,
        contig   = clean_contig(ID_replicon),
        start    = as.integer(pos_beg),
        end      = as.integer(pos_end),
        mge_type = "Integron",
        mge_name = paste0("Integron_class", type),
        tool     = "IntegronFinder"
      ) %>%
      filter(!is.na(start), !is.na(end), !is.na(contig))
  }, error = function(e) { NULL })  # No integrons found is normal
}

# ── 3. geNomad — Plasmids ────────────────────────────────────────────────────
# Columns: seq_name length topology n_genes plasmid_score fdr n_hallmarks ...
load_genomad <- function(sample, score_threshold = 0.7) {
  f <- file.path(BASE, "mobile_elements", sample,
                 "genomad/consensus_summary/consensus_plasmid_summary.tsv")
  if (!file.exists(f)) return(NULL)
  tryCatch({
    d <- read_tsv(f, show_col_types = FALSE)
    if (nrow(d) == 0) return(NULL)
    d %>%
      filter(plasmid_score >= score_threshold) %>%
      transmute(
        sample   = sample,
        contig   = clean_contig(seq_name),
        start    = 1L,
        end      = as.integer(length),
        mge_type = "Plasmid",
        mge_name = paste0("Plasmid_geNomad"),
        tool     = "geNomad"
      )
  }, error = function(e) { message("geNomad error ", sample, ": ", e$message); NULL })
}

# ── 4. MOBsuite — Plasmids ───────────────────────────────────────────────────
# Columns: sample_id molecule_type primary_cluster_id contig_id size predicted_mobility ...
load_mobsuite <- function(sample) {
  f <- file.path(BASE, "mobile_elements", sample, "mobsuite/contig_report.txt")
  if (!file.exists(f)) return(NULL)
  tryCatch({
    d <- read_tsv(f, show_col_types = FALSE)
    if (nrow(d) == 0) return(NULL)
    d %>%
      filter(molecule_type == "plasmid") %>%
      transmute(
        sample   = sample,
        contig   = clean_contig(contig_id),
        start    = 1L,
        end      = as.integer(size),
        mge_type = "Plasmid",
        mge_name = paste0("Plasmid_", primary_cluster_id),
        tool     = "MOBsuite"
      ) %>%
      filter(!is.na(contig))
  }, error = function(e) { message("MOBsuite error ", sample, ": ", e$message); NULL })
}

# ── Combine all ───────────────────────────────────────────────────────────────
cat("Loading ISEScan results...\n")
is_all   <- map_dfr(samples, load_isescan);   cat(" ", nrow(is_all),   "IS elements\n")

cat("Loading IntegronFinder results...\n")
int_all  <- map_dfr(samples, load_integrons); cat(" ", nrow(int_all),  "integrons\n")

cat("Loading geNomad results (plasmid score ≥0.7)...\n")
gen_all  <- map_dfr(samples, load_genomad);   cat(" ", nrow(gen_all),  "plasmid contigs\n")

cat("Loading MOBsuite results...\n")
mob_all  <- map_dfr(samples, load_mobsuite);  cat(" ", nrow(mob_all),  "plasmid contigs\n")

mges_combined <- bind_rows(is_all, int_all, gen_all, mob_all) %>%
  filter(!is.na(contig), !is.na(start), !is.na(end))

outfile <- file.path(OUTDIR, "mges_combined.tsv")
write_tsv(mges_combined, outfile)

cat("\n======================================\n")
cat("MGEs combined and saved\n")
cat("  Total MGE records:   ", nrow(mges_combined), "\n")
cat("  Samples with MGEs:   ", n_distinct(mges_combined$sample), "\n")
cat("  By type:\n")
print(count(mges_combined, mge_type, tool) %>% arrange(mge_type, desc(n)))
cat("  Output:", outfile, "\n")
cat("======================================\n")
