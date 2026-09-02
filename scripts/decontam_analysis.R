# =============================================================================
# decontam_analysis.R
# Step 4.1 — Remove contaminant species using negative controls
#
# WHAT THIS SCRIPT DOES:
#   Identifies and removes species that are more abundant in negative controls
#   (AE-buffer, EB-neg) than in real samples. These are reagent contaminants,
#   not genuine microbiome members. Removing them before diversity analysis
#   prevents inflated richness and biased community comparisons.
#
# METHOD:
#   Prevalence method (isContaminant, threshold = 0.5)
#   A species is flagged as a contaminant if its prevalence score > 0.5,
#   meaning it is more prevalent in negative controls than in true samples.
#   We use prevalence rather than frequency because negative controls have
#   NA DNA concentration values in the metadata.
#
# INPUT:
#   combined_bracken_all_samples.txt  — 63 samples x 9482 species count matrix
#   samples.txt                       — metadata with negative control labels
#
# OUTPUT:
#   decontam_results_all_species.csv  — full results table (all species)
#   contaminant_species_list.txt      — list of flagged contaminants
#   decontam_prevalence_plot.pdf      — visualisation of prevalence scores
#   decontam_cleaned_matrix.txt       — count matrix with contaminants removed
#   ps_all_samples.rds                — phyloseq object (all 63 samples, clean)
#   decontam_summary.txt              — summary statistics
#
# Run AFTER: run_krakentools.sh
# Run BEFORE: run_phyloseq.sh and run_vegan.sh
# =============================================================================

suppressPackageStartupMessages({
  library(decontam)
  library(phyloseq)
  library(ggplot2)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE        <- "/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
BRACKEN     <- file.path(BASE, "kraken2_output", "combined_bracken_all_samples.txt")
META        <- file.path(BASE, "samples.txt")
OUT_DIR     <- file.path(BASE, "r_analysis", "decontam")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("======================================================\n")
cat("decontam — Contamination Removal\n")
cat("Start time:", format(Sys.time()), "\n")
cat("======================================================\n\n")

# ── Load Bracken combined matrix ──────────────────────────────────────────────
cat("Loading Bracken combined matrix...\n")
bracken <- read.table(BRACKEN, sep = "\t", header = TRUE,
                      quote = "", check.names = FALSE)
cat("  Raw dimensions:", nrow(bracken), "species x", ncol(bracken), "columns\n")

# Extract raw count columns (_num suffix) — ignore _frac columns
count_cols <- grep("_num$", colnames(bracken), value = TRUE)
otu_mat    <- as.matrix(bracken[, count_cols])
colnames(otu_mat) <- gsub("_num$", "", colnames(otu_mat))  # strip _num suffix
rownames(otu_mat) <- bracken$name
cat("  OTU table:", nrow(otu_mat), "species x", ncol(otu_mat), "samples\n")

# ── Load and align metadata ───────────────────────────────────────────────────
cat("\nLoading metadata...\n")
meta <- read.table(META, sep = "\t", header = TRUE,
                   quote = "", check.names = FALSE, comment.char = "")
rownames(meta) <- meta[["Sample ID"]]
cat("  Metadata rows:", nrow(meta), "\n")

# Keep only samples present in both files
common <- intersect(colnames(otu_mat), rownames(meta))
cat("  Matched samples:", length(common), "\n")
if (length(common) < ncol(otu_mat)) {
  cat("  WARNING: some Bracken samples not found in metadata:\n")
  cat("   ", paste(setdiff(colnames(otu_mat), rownames(meta)), collapse = ", "), "\n")
}
otu_mat <- otu_mat[, common]
meta    <- meta[common, ]

# ── Build phyloseq object (all 63 samples) ───────────────────────────────────
cat("\nBuilding phyloseq object...\n")
OTU  <- otu_table(otu_mat, taxa_are_rows = TRUE)
SAMP <- sample_data(meta)
ps   <- phyloseq(OTU, SAMP)
cat("  phyloseq:", ntaxa(ps), "taxa,", nsamples(ps), "samples\n")

# ── Run decontam ──────────────────────────────────────────────────────────────
cat("\nRunning decontam (prevalence method, threshold = 0.5)...\n")

# Mark negative controls: AE-buffer and EB-neg have sample == "negative_control"
sample_data(ps)$is_neg <- sample_data(ps)[["sample"]] == "negative_control"
n_neg  <- sum(sample_data(ps)$is_neg)
n_real <- sum(!sample_data(ps)$is_neg)
cat("  Negative controls:", n_neg, "(AE-buffer, EB-neg)\n")
cat("  True samples:     ", n_real, "\n")

# isContaminant: prevalence method
# A species scores high (closer to 1) if it appears more in controls than samples
# Threshold 0.5 means: flag as contaminant if it is more common in controls
contamdf <- isContaminant(ps, method = "prevalence",
                          neg = "is_neg", threshold = 0.5)
n_contam <- sum(contamdf$contaminant, na.rm = TRUE)
n_clean  <- sum(!contamdf$contaminant, na.rm = TRUE)

cat("\n  Species tested:          ", nrow(contamdf), "\n")
cat("  Contaminants identified: ", n_contam, "\n")
cat("  Clean species retained:  ", n_clean, "\n")

# ── Save full results table ───────────────────────────────────────────────────
write.csv(contamdf,
          file = file.path(OUT_DIR, "decontam_results_all_species.csv"),
          quote = FALSE)

# List contaminant species names
contam_species <- rownames(contamdf)[contamdf$contaminant & !is.na(contamdf$contaminant)]
writeLines(contam_species,
           con = file.path(OUT_DIR, "contaminant_species_list.txt"))

if (length(contam_species) > 0) {
  cat("\nContaminant species flagged:\n")
  cat(paste(" -", contam_species, collapse = "\n"), "\n")
} else {
  cat("\n  No contaminants identified at threshold = 0.5\n")
  cat("  (This can happen with only 2 negative controls — the signal is limited)\n")
  cat("  Consider running with threshold = 0.1 for a more sensitive check\n")
}

# ── Prevalence score plot ─────────────────────────────────────────────────────
cat("\nGenerating prevalence score plot...\n")
contamdf$status <- ifelse(
  is.na(contamdf$contaminant), "Untested",
  ifelse(contamdf$contaminant, "Contaminant", "Non-contaminant")
)
contamdf$prev_score <- contamdf$p   # prevalence score column from decontam

p_plot <- ggplot(contamdf[!is.na(contamdf$p), ], aes(x = p, fill = status)) +
  geom_histogram(bins = 40, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.5, colour = "black", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~status, scales = "free_y") +
  scale_fill_manual(values = c("Contaminant"     = "#ef4444",
                               "Non-contaminant" = "#22c55e",
                               "Untested"        = "#94a3b8")) +
  labs(title    = "decontam -- Prevalence Score Distribution",
       subtitle  = paste0("Threshold = 0.5 | Prevalence method | ",
                          n_neg, " negative controls | ",
                          n_contam, " contaminants removed"),
       x = "Prevalence Score  (higher = more prevalent in controls)",
       y = "Number of Species") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "decontam_prevalence_plot.pdf"),
       p_plot, width = 12, height = 5)

# ── Remove contaminants and export cleaned matrix ─────────────────────────────
cat("\nRemoving contaminants and exporting cleaned matrix...\n")
keep    <- !contamdf$contaminant | is.na(contamdf$contaminant)
ps_clean <- prune_taxa(keep, ps)
cat("  Cleaned phyloseq:", ntaxa(ps_clean), "taxa,", nsamples(ps_clean), "samples\n")

# Export as tab-delimited count matrix
clean_otu_df       <- as.data.frame(otu_table(ps_clean))
clean_otu_df       <- cbind(name = rownames(clean_otu_df), clean_otu_df)
write.table(clean_otu_df,
            file      = file.path(OUT_DIR, "decontam_cleaned_matrix.txt"),
            sep       = "\t",
            row.names = FALSE,
            quote     = FALSE)

# Save phyloseq object (all 63 samples including controls) for downstream use
saveRDS(ps_clean, file = file.path(OUT_DIR, "ps_all_samples.rds"))

# ── Summary report ────────────────────────────────────────────────────────────
sink(file.path(OUT_DIR, "decontam_summary.txt"))
cat("decontam Summary Report\n")
cat("========================\n")
cat("Date:               ", format(Sys.time()), "\n")
cat("Input species:      ", nrow(bracken), "\n")
cat("Input samples:      ", ncol(otu_mat), "\n")
cat("Negative controls:  ", n_neg, "(AE-buffer, EB-neg)\n")
cat("True samples:       ", n_real, "\n")
cat("Contaminants found: ", n_contam, "\n")
cat("Species retained:   ", n_clean, "\n\n")
if (length(contam_species) > 0) {
  cat("Contaminant species:\n")
  cat(paste(" -", contam_species, collapse = "\n"), "\n")
} else {
  cat("No contaminants identified at threshold = 0.5\n")
}
sink()

cat("\n======================================================\n")
cat("decontam complete!\n")
cat("Output directory:", OUT_DIR, "\n")
cat("Files created:\n")
cat("  decontam_results_all_species.csv\n")
cat("  contaminant_species_list.txt\n")
cat("  decontam_prevalence_plot.pdf\n")
cat("  decontam_cleaned_matrix.txt\n")
cat("  ps_all_samples.rds\n")
cat("  decontam_summary.txt\n")
cat("End time:", format(Sys.time()), "\n")
cat("======================================================\n")
