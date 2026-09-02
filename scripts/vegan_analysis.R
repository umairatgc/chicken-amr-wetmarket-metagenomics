# =============================================================================
# vegan_analysis.R
# Step 4.3 — PERMANOVA, Betadispersion, and NMDS
#
# WHAT THIS SCRIPT DOES:
#   Tests whether microbial community composition differs between groups using
#   multivariate statistics on Bray-Curtis dissimilarities. Complements the
#   phyloseq alpha/beta diversity results with formal hypothesis tests.
#
# ANALYSES:
#   1. PERMANOVA (adonis2)  — does community composition differ between groups?
#   2. Betadispersion       — are differences due to location or spread?
#   3. NMDS ordination      — visual 2D representation of community distances
#
# GROUPING VARIABLES TESTED:
#   - Environment    (Butcher vs Farm)
#   - Sample Type    (fine-grained: Flies, Meat, Drainage, …)
#   - Source 1       (Shop 1-6, farms)
#   - mcr-1, tetX4, fosA3, tmex, NDM, OXA-48, KPC  (binary PCR results)
#     → skipped for any gene where either group has < 3 samples
#
# INPUT:
#   ps_real_samples.rds  — phyloseq object from phyloseq_analysis.R
#                          (61 real samples, negative controls removed,
#                           contaminants removed, rarefied)
#
# OUTPUT (r_analysis/vegan/):
#   permanova_results.csv       — table of all PERMANOVA tests
#   betadisper_results.csv      — table of all betadispersion tests
#   nmds_environment.pdf        — NMDS coloured by Butcher vs Farm
#   nmds_sampletype.pdf         — NMDS coloured by Sample Type
#   nmds_mcr1.pdf               — NMDS coloured by mcr-1 PCR result
#   nmds_all_pcr.pdf            — multi-panel NMDS for all 7 PCR genes
#   vegan_summary.txt           — plain text summary of all results
#
# Run AFTER: run_phyloseq.sh (needs ps_real_samples.rds)
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    <- "/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
PS_IN   <- file.path(BASE, "r_analysis", "phyloseq", "ps_real_samples.rds")
OUT_DIR <- file.path(BASE, "r_analysis", "vegan")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("======================================================\n")
cat("vegan — PERMANOVA, Betadispersion & NMDS\n")
cat("Start time:", format(Sys.time()), "\n")
cat("======================================================\n\n")

# ── Load rarefied phyloseq object from phyloseq_analysis.R ────────────────────
cat("Loading ps_real_samples.rds ...\n")
if (!file.exists(PS_IN)) {
  stop("ps_real_samples.rds not found at: ", PS_IN,
       "\nRun phyloseq_analysis.R first (sbatch run_phyloseq.sh)")
}
ps <- readRDS(PS_IN)
cat("  phyloseq:", ntaxa(ps), "taxa,", nsamples(ps), "samples\n\n")

# ── Extract OTU table (samples as rows) and metadata ─────────────────────────
otu  <- as.data.frame(t(otu_table(ps)))    # samples × taxa
meta <- as.data.frame(sample_data(ps))

# ── Compute Bray-Curtis distance matrix ──────────────────────────────────────
cat("Computing Bray-Curtis dissimilarities...\n")
bc_dist <- vegdist(otu, method = "bray")
cat("  Distance matrix:", nsamples(ps), "x", nsamples(ps), "\n\n")

# ── Colour palettes ───────────────────────────────────────────────────────────
env_colours <- c(Butcher = "#f97316", Farm = "#22c55e")
pcr_colours <- c(Negative = "#94a3b8", Positive = "#ef4444")

# Palette for Sample Type (up to 18 types)
sampletype_pal <- c(
  "#f97316","#fb923c","#fdba74","#fde68a",   # warm orange/yellow (Butcher)
  "#ea580c","#c2410c","#b45309","#92400e",   # darker butcher
  "#22c55e","#4ade80","#86efac","#bbf7d0",   # green (Farm)
  "#16a34a","#15803d","#166534","#14532d",   # darker farm
  "#0ea5e9","#7c3aed"                        # extras
)

# ── Helper: run PERMANOVA + betadispersion for one grouping variable ──────────
# Returns list: permanova_row, betadisper_row
run_tests <- function(group_var, label, min_n = 3) {

  grp <- meta[[group_var]]

  # Check for at least 2 levels, each with >= min_n samples
  tbl <- table(grp)
  tbl <- tbl[tbl >= min_n]
  if (length(tbl) < 2) {
    msg <- paste0("  SKIP (", label, "): fewer than 2 groups with >= ",
                  min_n, " samples — skipping")
    cat(msg, "\n")
    return(NULL)
  }

  # Subset if some levels were dropped (unlikely after binary recoding, but safe)
  keep <- grp %in% names(tbl)
  otu_sub  <- otu[keep, ]
  grp_sub  <- grp[keep]
  dist_sub <- as.dist(as.matrix(bc_dist)[keep, keep])

  # ── PERMANOVA ──────────────────────────────────────────────────────────────
  set.seed(42)
  perm <- adonis2(dist_sub ~ grp_sub, permutations = 999)
  R2  <- round(perm$R2[1], 4)
  pval <- round(perm$`Pr(>F)`[1], 4)
  F_stat <- round(perm$F[1], 3)
  n_used <- sum(keep)

  cat(sprintf("  PERMANOVA  %-22s  R2=%.4f  F=%.3f  p=%s  (n=%d)\n",
              label, R2, F_stat,
              ifelse(pval < 0.001, "<0.001", pval), n_used))

  perm_row <- data.frame(
    Variable     = label,
    n_samples    = n_used,
    n_groups     = length(tbl),
    R2           = R2,
    F_stat       = F_stat,
    p_value      = pval,
    Significance = ifelse(pval <= 0.001, "***",
                   ifelse(pval <= 0.01,  "**",
                   ifelse(pval <= 0.05,  "*", "ns")))
  )

  # ── Betadispersion ────────────────────────────────────────────────────────
  bd   <- betadisper(dist_sub, grp_sub)
  bd_p <- round(permutest(bd, permutations = 999)$tab$`Pr(>F)`[1], 4)
  bd_F <- round(permutest(bd, permutations = 999)$tab$F[1], 3)

  cat(sprintf("  Betadisper %-22s  F=%.3f  p=%s\n",
              label, bd_F,
              ifelse(bd_p < 0.001, "<0.001", bd_p)))

  bd_row <- data.frame(
    Variable     = label,
    n_samples    = n_used,
    n_groups     = length(tbl),
    F_stat       = bd_F,
    p_value      = bd_p,
    Interpretation = ifelse(bd_p > 0.05,
      "Homogeneous dispersion (PERMANOVA result is reliable)",
      "Heterogeneous dispersion (interpret PERMANOVA with caution)")
  )

  list(perm_row = perm_row, bd_row = bd_row)
}

# ── Define all tests to run ───────────────────────────────────────────────────
pcr_short  <- c("mcr1", "tetX4", "fosA3", "tmex", "NDM", "OXA48", "KPC")
pcr_labels <- c("mcr-1", "tetX4", "fosA3", "tmex", "NDM", "OXA-48", "KPC")

test_vars <- list(
  list(col = "Environment",   label = "Environment"),
  list(col = "Sample Type",   label = "Sample Type"),
  list(col = "Source 1",      label = "Source 1")
)
for (i in seq_along(pcr_short)) {
  col <- paste0(pcr_short[i], "_pos")
  test_vars <- c(test_vars, list(list(col = col, label = pcr_labels[i])))
}

# ── Run all tests ─────────────────────────────────────────────────────────────
cat("Running PERMANOVA and betadispersion tests...\n\n")
perm_rows <- list()
bd_rows   <- list()

for (tv in test_vars) {
  if (!tv$col %in% colnames(meta)) {
    cat("  WARNING: column '", tv$col, "' not found in metadata — skipping\n", sep = "")
    next
  }
  res <- run_tests(tv$col, tv$label)
  if (!is.null(res)) {
    perm_rows <- c(perm_rows, list(res$perm_row))
    bd_rows   <- c(bd_rows,   list(res$bd_row))
  }
}

# ── Save results tables ───────────────────────────────────────────────────────
perm_df <- do.call(rbind, perm_rows)
bd_df   <- do.call(rbind, bd_rows)

write.csv(perm_df, file = file.path(OUT_DIR, "permanova_results.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(bd_df,   file = file.path(OUT_DIR, "betadisper_results.csv"),
          row.names = FALSE, quote = FALSE)

cat("\n  Saved: permanova_results.csv\n")
cat("  Saved: betadisper_results.csv\n\n")

# ── NMDS ordination ───────────────────────────────────────────────────────────
cat("Running NMDS (k=2, trymax=100)...\n")
set.seed(42)
nmds <- metaMDS(otu, distance = "bray", k = 2, trymax = 100, trace = FALSE)
cat(sprintf("  Stress = %.4f", nmds$stress))
if (nmds$stress < 0.1) {
  cat(" — Excellent representation\n")
} else if (nmds$stress < 0.2) {
  cat(" — Good representation\n")
} else {
  cat(" — Fair (consider k=3 for publication)\n")
}

# Build score data frame with metadata
scores_df <- as.data.frame(scores(nmds, display = "sites"))
scores_df$SampleID <- rownames(scores_df)
scores_df <- merge(scores_df, meta, by.x = "SampleID", by.y = "row.names",
                   all.x = TRUE)

stress_lab <- sprintf("Stress = %.4f", nmds$stress)

# ── Helper: NMDS plot function ────────────────────────────────────────────────
make_nmds <- function(df, colour_col, title, palette, shape_col = NULL,
                      file_out, w = 8, h = 6) {

  # remove rows where colour column is NA
  df <- df[!is.na(df[[colour_col]]), ]

  aes_args <- if (!is.null(shape_col) && shape_col %in% colnames(df)) {
    aes(x = NMDS1, y = NMDS2,
        colour = .data[[colour_col]],
        shape  = .data[[shape_col]])
  } else {
    aes(x = NMDS1, y = NMDS2, colour = .data[[colour_col]])
  }

  p <- ggplot(df, aes_args) +
    geom_point(size = 3.5, alpha = 0.85) +
    stat_ellipse(aes(colour = .data[[colour_col]]), level = 0.95,
                 linetype = "dashed", linewidth = 0.6) +
    annotate("text", x = Inf, y = Inf, label = stress_lab,
             hjust = 1.1, vjust = 1.5, size = 3.5, colour = "grey40") +
    labs(title   = title,
         x = "NMDS1", y = "NMDS2",
         colour = colour_col) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")

  # Apply named palette if provided and length matches
  if (is.character(palette) && !is.null(names(palette))) {
    p <- p + scale_colour_manual(values = palette)
  } else if (length(unique(df[[colour_col]])) <= length(palette)) {
    p <- p + scale_colour_manual(values = palette)
  }

  ggsave(file_out, p, width = w, height = h)
  cat("  Saved:", basename(file_out), "\n")
}

cat("\nGenerating NMDS plots...\n")

# 1. Environment (Butcher vs Farm)
# shape_col omitted — 16 sample types exceeds ggplot's 6-shape limit,
# which drops points. Colour alone distinguishes Butcher vs Farm clearly.
make_nmds(
  df         = scores_df,
  colour_col = "Environment",
  title      = "NMDS -- Bray-Curtis | Butcher vs Farm",
  palette    = env_colours,
  file_out   = file.path(OUT_DIR, "nmds_environment.pdf")
)

# 2. Sample Type
make_nmds(
  df         = scores_df,
  colour_col = "Sample Type",
  title      = "NMDS -- Bray-Curtis | Sample Type",
  palette    = sampletype_pal,
  file_out   = file.path(OUT_DIR, "nmds_sampletype.pdf"),
  w = 10, h = 7
)

# 3. mcr-1 (standalone highlighted plot)
if ("mcr1_pos" %in% colnames(scores_df)) {
  make_nmds(
    df         = scores_df,
    colour_col = "mcr1_pos",
    title      = "NMDS -- Bray-Curtis | mcr-1 PCR Result",
    palette    = pcr_colours,
    shape_col  = "Environment",
    file_out   = file.path(OUT_DIR, "nmds_mcr1.pdf")
  )
}

# 4. Multi-panel: all 7 PCR genes
cat("  Building multi-panel PCR NMDS...\n")
pcr_panels <- list()
for (i in seq_along(pcr_short)) {
  col <- paste0(pcr_short[i], "_pos")
  if (!col %in% colnames(scores_df)) next
  tmp <- scores_df[!is.na(scores_df[[col]]), ]
  if (nrow(tmp) < 6) next   # skip if too few to plot meaningfully
  grp_counts <- table(tmp[[col]])
  if (length(grp_counts) < 2) next

  p <- ggplot(tmp, aes(x = NMDS1, y = NMDS2, colour = .data[[col]])) +
    geom_point(size = 2.5, alpha = 0.85) +
    stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.5) +
    scale_colour_manual(values = pcr_colours,
                        labels = c(Negative = "Neg", Positive = "Pos")) +
    annotate("text", x = Inf, y = Inf, label = stress_lab,
             hjust = 1.1, vjust = 1.5, size = 2.8, colour = "grey40") +
    labs(title  = pcr_labels[i],
         colour = NULL,
         x = "NMDS1", y = "NMDS2") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "bottom",
          legend.margin = margin(0, 0, 0, 0))

  pcr_panels[[pcr_labels[i]]] <- p
}

if (length(pcr_panels) > 0) {
  # Arrange panels using gridExtra (2 columns)
  library(gridExtra)
  ncols  <- 2
  nrows  <- ceiling(length(pcr_panels) / ncols)
  pdf(file.path(OUT_DIR, "nmds_all_pcr.pdf"),
      width = 5 * ncols, height = 5 * nrows)
  grid.arrange(grobs = pcr_panels, ncol = ncols,
               top = grid::textGrob(
                 "NMDS -- Bray-Curtis | PCR Gene Positivity",
                 gp = grid::gpar(fontface = "bold", fontsize = 14)))
  dev.off()
  cat("  Saved: nmds_all_pcr.pdf\n")
}

# ── Plain-text summary ─────────────────────────────────────────────────────────
cat("\nWriting summary report...\n")
sink(file.path(OUT_DIR, "vegan_summary.txt"))
cat("vegan Analysis Summary\n")
cat("=======================\n")
cat("Date:            ", format(Sys.time()), "\n")
cat("Input:           ps_real_samples.rds (rarefied, contaminants removed)\n")
cat("Samples:         ", nsamples(ps), "\n")
cat("Taxa:            ", ntaxa(ps), "\n")
cat("NMDS stress:     ", round(nmds$stress, 4), "\n\n")

cat("PERMANOVA Results (adonis2, Bray-Curtis, 999 permutations)\n")
cat("------------------------------------------------------------\n")
print(perm_df, row.names = FALSE)

cat("\n\nBetadispersion Results (permutest, 999 permutations)\n")
cat("------------------------------------------------------\n")
print(bd_df, row.names = FALSE)

cat("\n\nInterpretation Guide\n")
cat("--------------------\n")
cat("PERMANOVA R2      : proportion of variance explained by grouping (0–1)\n")
cat("PERMANOVA p-value : *** p<=0.001, ** p<=0.01, * p<=0.05, ns p>0.05\n")
cat("Betadispersion    : if p>0.05, groups have homogeneous spread —\n")
cat("                    PERMANOVA differences reflect true community shift\n")
cat("                    if p<=0.05, one group is more variable — interpret\n")
cat("                    PERMANOVA result with caution\n")
cat("NMDS Stress       : <0.1 Excellent, <0.2 Good, <0.3 Fair\n")
sink()

cat("  Saved: vegan_summary.txt\n\n")

# ── Final summary to console ──────────────────────────────────────────────────
cat("======================================================\n")
cat("vegan complete!\n")
cat("Output directory:", OUT_DIR, "\n")
cat("Files created:\n")
cat("  permanova_results.csv\n")
cat("  betadisper_results.csv\n")
cat("  nmds_environment.pdf\n")
cat("  nmds_sampletype.pdf\n")
cat("  nmds_mcr1.pdf\n")
cat("  nmds_all_pcr.pdf\n")
cat("  vegan_summary.txt\n")
cat("End time:", format(Sys.time()), "\n")
cat("======================================================\n")
