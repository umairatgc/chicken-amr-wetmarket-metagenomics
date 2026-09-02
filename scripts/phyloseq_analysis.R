# =============================================================================
# phyloseq_analysis.R
# Step 4.2 — Alpha and beta diversity analysis
#
# WHAT THIS SCRIPT DOES:
#   Analyses microbial community diversity across 61 real samples (negative
#   controls removed). Produces alpha diversity plots (richness and evenness
#   per sample), beta diversity ordination (community similarity between
#   samples), and relative abundance bar charts.
#
# GROUPING VARIABLES USED:
#   - Environment: "Butcher" vs "Farm" (derived from Sample Type)
#   - Sample Type: Butcher Flies, Butcher Knife, Farm Litter, etc.
#   - PCR results: mcr-1, tetX4, fosA3, tmex, NDM, OXA-48, KPC
#     (recoded as binary: positive = 1, negative/absent = 0)
#
# ALPHA DIVERSITY METRICS:
#   Observed species, Shannon (evenness + richness), Chao1 (richness estimator),
#   Simpson (evenness-weighted)
#
# BETA DIVERSITY:
#   Bray-Curtis distance, PCoA ordination
#
# INPUT:
#   r_analysis/decontam/ps_all_samples.rds  — cleaned phyloseq object
#   samples.txt                              — metadata
#
# OUTPUT (all in r_analysis/phyloseq/):
#   alpha_by_environment.pdf         — Shannon/Chao1/Simpson: Butcher vs Farm
#   alpha_by_sample_type.pdf         — all Sample Types
#   alpha_by_pcr_*.pdf               — one plot per PCR gene
#   alpha_stats.csv                  — Wilcoxon/Kruskal-Wallis test results
#   pcoa_environment.pdf             — PCoA coloured by Environment
#   pcoa_sample_type.pdf             — PCoA coloured by Sample Type
#   pcoa_mcr1.pdf                    — PCoA coloured by mcr-1 status
#   relative_abundance_top20.pdf     — stacked bar chart (top 20 species)
#   ps_real_samples.rds              — phyloseq object (61 samples, rarefied)
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(dplyr)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    <- "/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
PS_FILE <- file.path(BASE, "r_analysis", "decontam", "ps_all_samples.rds")
META    <- file.path(BASE, "samples.txt")
OUT_DIR <- file.path(BASE, "r_analysis", "phyloseq")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("======================================================\n")
cat("phyloseq — Diversity Analysis\n")
cat("Start time:", format(Sys.time()), "\n")
cat("======================================================\n\n")

# ── Load phyloseq object and metadata ────────────────────────────────────────
cat("Loading cleaned phyloseq object...\n")
ps_all <- readRDS(PS_FILE)
cat("  All samples:", ntaxa(ps_all), "taxa,", nsamples(ps_all), "samples\n")

# ── Remove negative controls — keep only 61 real samples ─────────────────────
cat("Removing negative controls...\n")
ps <- subset_samples(ps_all, sample != "negative_control")
cat("  Real samples:", nsamples(ps), "\n")

# ── Reload metadata and prepare variables ─────────────────────────────────────
cat("Preparing metadata variables...\n")
meta <- read.table(META, sep = "\t", header = TRUE,
                   quote = "", check.names = FALSE, comment.char = "")
rownames(meta) <- meta[["Sample ID"]]

# Keep only real samples
meta <- meta[meta[["sample"]] == "sample", ]

# Create broad Environment variable (Butcher vs Farm)
meta$Environment <- ifelse(
  grepl("^Butcher", meta[["Sample Type"]]), "Butcher",
  ifelse(grepl("^Farm", meta[["Sample Type"]]), "Farm", "Other")
)

# Recode PCR columns to binary (any positive result = 1, blank/NA = 0)
pcr_genes <- c("mcr-1 C2", "tetX4 C2", "fosA3 C2", "tmex multi C2",
               "NDM C2", "OXA-48 C2", "KPC C2")
pcr_short  <- c("mcr1", "tetX4", "fosA3", "tmex", "NDM", "OXA48", "KPC")

for (i in seq_along(pcr_genes)) {
  val <- meta[[pcr_genes[i]]]
  meta[[paste0(pcr_short[i], "_pos")]] <- factor(
    ifelse(is.na(val) | val == "" | val == "NA", "Negative", "Positive"),
    levels = c("Negative", "Positive")
  )
}

# Update sample data in phyloseq with enriched metadata
common_samples <- intersect(sample_names(ps), rownames(meta))
ps <- prune_samples(common_samples, ps)
sample_data(ps) <- sample_data(meta[common_samples, ])

cat("  Environment counts:\n")
print(table(meta$Environment))
cat("  PCR gene positivity:\n")
for (i in seq_along(pcr_short)) {
  tbl <- table(meta[[paste0(pcr_short[i], "_pos")]])
  cat("   ", pcr_genes[i], "— Positive:", tbl["Positive"],
      "Negative:", tbl["Negative"], "\n")
}

# ── Filter low-abundance species ──────────────────────────────────────────────
# Remove species with total counts < 2 across all samples (likely noise)
cat("\nFiltering low-abundance species...\n")
ps <- prune_taxa(taxa_sums(ps) >= 2, ps)
cat("  Species after filter:", ntaxa(ps), "\n")

# ── Rarefy to even depth ──────────────────────────────────────────────────────
# Required for fair comparison of alpha diversity between samples with
# different sequencing depths
cat("\nRarefying to even sequencing depth...\n")
lib_sizes <- sample_sums(ps)
min_lib   <- min(lib_sizes)
cat("  Library sizes — min:", min_lib, " max:", max(lib_sizes),
    " median:", round(median(lib_sizes)), "\n")
cat("  Rarefying to:", min_lib, "reads per sample\n")
set.seed(42)
ps_rare <- rarefy_even_depth(ps, sample.size = min_lib,
                              replace = FALSE, verbose = FALSE)
cat("  After rarefaction:", ntaxa(ps_rare), "species,", nsamples(ps_rare), "samples\n")

# Save rarefied phyloseq object
saveRDS(ps_rare, file = file.path(OUT_DIR, "ps_real_samples.rds"))

# ── Helper: alpha diversity plot ──────────────────────────────────────────────
plot_alpha <- function(ps_obj, group_var, fill_colours = NULL, title = "") {
  alpha_df <- estimate_richness(ps_obj,
                                measures = c("Observed", "Shannon", "Chao1", "Simpson"))
  alpha_df$SampleID   <- rownames(alpha_df)
  alpha_df[[group_var]] <- as.character(sample_data(ps_obj)[[group_var]])

  # Reshape to long format
  alpha_long <- reshape(alpha_df,
                        varying   = c("Observed", "Shannon", "Chao1", "Simpson"),
                        v.names   = "Value",
                        timevar   = "Metric",
                        times     = c("Observed", "Shannon", "Chao1", "Simpson"),
                        direction = "long")
  alpha_long$Metric <- factor(alpha_long$Metric,
                               levels = c("Observed", "Chao1", "Shannon", "Simpson"))

  p <- ggplot(alpha_long, aes(x = .data[[group_var]],
                               y = Value, fill = .data[[group_var]])) +
    geom_boxplot(outlier.shape = 21, outlier.size = 2, linewidth = 0.5, alpha = 0.85) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, colour = "black") +
    facet_wrap(~Metric, scales = "free_y", ncol = 4) +
    labs(title = title, x = NULL, y = "Diversity Value") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x     = element_text(angle = 35, hjust = 1, size = 10),
          plot.title      = element_text(face = "bold"),
          strip.background = element_rect(fill = "#1e293b"),
          strip.text       = element_text(colour = "white", face = "bold"))
  if (!is.null(fill_colours)) {
    p <- p + scale_fill_manual(values = fill_colours)
  }
  return(p)
}

# ── Alpha diversity: by Environment ──────────────────────────────────────────
cat("\nPlotting alpha diversity by Environment...\n")
env_colours <- c("Butcher" = "#f97316", "Farm" = "#22c55e")
p_alpha_env <- plot_alpha(ps_rare, "Environment", env_colours,
                          "Alpha Diversity -- Butcher vs Farm Environments")
ggsave(file.path(OUT_DIR, "alpha_by_environment.pdf"),
       p_alpha_env, width = 14, height = 5)

# ── Alpha diversity: by Sample Type ──────────────────────────────────────────
cat("Plotting alpha diversity by Sample Type...\n")
p_alpha_type <- plot_alpha(ps_rare, "Sample Type", title =
                             "Alpha Diversity -- By Sample Type")
ggsave(file.path(OUT_DIR, "alpha_by_sample_type.pdf"),
       p_alpha_type, width = 18, height = 6)

# ── Alpha diversity: by each PCR gene ────────────────────────────────────────
cat("Plotting alpha diversity by PCR gene status...\n")
pcr_colours <- c("Negative" = "#94a3b8", "Positive" = "#ef4444")
for (i in seq_along(pcr_short)) {
  col_name <- paste0(pcr_short[i], "_pos")
  tbl <- table(sample_data(ps_rare)[[col_name]])
  # Only plot if both groups have >= 3 samples
  if (length(tbl) == 2 && all(tbl >= 3)) {
    p_pcr <- plot_alpha(ps_rare, col_name, pcr_colours,
                        paste0("Alpha Diversity -- ", pcr_genes[i],
                               " Status (Positive=", tbl["Positive"],
                               ", Negative=", tbl["Negative"], ")"))
    ggsave(file.path(OUT_DIR, paste0("alpha_by_", pcr_short[i], ".pdf")),
           p_pcr, width = 14, height = 5)
  } else {
    cat("  Skipping", pcr_genes[i], "— insufficient group sizes:",
        paste(names(tbl), tbl, sep = "=", collapse = ", "), "\n")
  }
}

# ── Alpha diversity statistics ────────────────────────────────────────────────
cat("Running statistical tests for alpha diversity...\n")
alpha_df <- estimate_richness(ps_rare, measures = c("Observed", "Shannon", "Chao1", "Simpson"))
alpha_df$Environment  <- sample_data(ps_rare)$Environment
alpha_df$Sample.Type  <- sample_data(ps_rare)[["Sample Type"]]

stats_list <- list()
metrics    <- c("Observed", "Shannon", "Chao1", "Simpson")

# Wilcoxon test: Butcher vs Farm
for (m in metrics) {
  butcher <- alpha_df[[m]][alpha_df$Environment == "Butcher"]
  farm    <- alpha_df[[m]][alpha_df$Environment == "Farm"]
  w       <- wilcox.test(butcher, farm)
  stats_list[[paste0(m, "_Env")]] <- data.frame(
    Metric = m, Comparison = "Butcher_vs_Farm",
    Test = "Wilcoxon", P_value = round(w$p.value, 4),
    Median_A = round(median(butcher, na.rm = TRUE), 3),
    Median_B = round(median(farm, na.rm = TRUE), 3)
  )
}

# Kruskal-Wallis: across Sample Types
for (m in metrics) {
  kw <- kruskal.test(alpha_df[[m]] ~ factor(alpha_df$Sample.Type))
  stats_list[[paste0(m, "_SampleType")]] <- data.frame(
    Metric = m, Comparison = "SampleType_KW",
    Test = "Kruskal-Wallis", P_value = round(kw$p.value, 4),
    Median_A = NA, Median_B = NA
  )
}

stats_df <- do.call(rbind, stats_list)
write.csv(stats_df, file = file.path(OUT_DIR, "alpha_stats.csv"),
          row.names = FALSE, quote = FALSE)
cat("  Alpha diversity tests saved\n")

# ── Beta diversity: PCoA on Bray-Curtis distances ─────────────────────────────
cat("\nComputing Bray-Curtis distances and PCoA...\n")
ps_rel  <- transform_sample_counts(ps_rare, function(x) x / sum(x))

# phyloseq's plot_ordination cannot find columns with spaces in their names.
# Add a no-space copy of "Sample Type" so colour mapping works correctly.
sample_data(ps_rel)$SampleType <- sample_data(ps_rel)[["Sample Type"]]
ord_bc  <- ordinate(ps_rel, method = "PCoA", distance = "bray")
var_exp <- round(ord_bc$values$Relative_eig[1:2] * 100, 1)

# PCoA by Environment
p_pcoa_env <- plot_ordination(ps_rel, ord_bc, color = "Environment",
                               shape = "Environment") +
  geom_point(size = 3, alpha = 0.85) +
  scale_colour_manual(values = env_colours) +
  scale_shape_manual(values  = c("Butcher" = 16, "Farm" = 17)) +
  stat_ellipse(type = "t", linewidth = 0.7, level = 0.95) +
  labs(title    = "PCoA -- Bray-Curtis Dissimilarity",
       subtitle  = "Coloured by Environment (Butcher vs Farm)",
       x = paste0("PC1 [", var_exp[1], "%]"),
       y = paste0("PC2 [", var_exp[2], "%]")) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "pcoa_environment.pdf"),
       p_pcoa_env, width = 8, height = 6)

# PCoA by Sample Type
p_pcoa_type <- plot_ordination(ps_rel, ord_bc, color = "SampleType") +
  geom_point(size = 3, alpha = 0.85) +
  labs(title    = "PCoA -- Bray-Curtis Dissimilarity",
       subtitle  = "Coloured by Sample Type",
       colour   = "Sample Type",
       x = paste0("PC1 [", var_exp[1], "%]"),
       y = paste0("PC2 [", var_exp[2], "%]")) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "pcoa_sample_type.pdf"),
       p_pcoa_type, width = 10, height = 6)

# PCoA by mcr-1 status
p_pcoa_mcr1 <- plot_ordination(ps_rel, ord_bc, color = "mcr1_pos",
                                 shape = "mcr1_pos") +
  geom_point(size = 3, alpha = 0.85) +
  scale_colour_manual(values = pcr_colours, name = "mcr-1 status") +
  scale_shape_manual(values  = c("Negative" = 1, "Positive" = 16),
                     name = "mcr-1 status") +
  stat_ellipse(type = "t", linewidth = 0.7, level = 0.95) +
  labs(title    = "PCoA -- Bray-Curtis Dissimilarity",
       subtitle  = "Coloured by mcr-1 PCR status",
       x = paste0("PC1 [", var_exp[1], "%]"),
       y = paste0("PC2 [", var_exp[2], "%]")) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "pcoa_mcr1.pdf"),
       p_pcoa_mcr1, width = 8, height = 6)

# ── Relative abundance bar chart — Top 20 species ────────────────────────────
cat("Generating relative abundance bar chart (top 20 species)...\n")
ps_top20 <- prune_taxa(names(sort(taxa_sums(ps_rel), decreasing = TRUE)[1:20]), ps_rel)

# Melt to long format
otu_melt <- psmelt(ps_top20)
otu_melt$Environment <- factor(otu_melt$Environment, levels = c("Butcher", "Farm"))

# Sort samples by Environment then sample name
otu_melt <- otu_melt[order(otu_melt$Environment, otu_melt$Sample), ]
otu_melt$Sample <- factor(otu_melt$Sample,
                           levels = unique(otu_melt$Sample))

p_relabund <- ggplot(otu_melt, aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~Environment, scales = "free_x", space = "free_x") +
  labs(title = "Relative Abundance -- Top 20 Species",
       x = NULL, y = "Relative Abundance",
       fill = "Species") +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        legend.position = "right",
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        plot.title   = element_text(face = "bold"),
        strip.background = element_rect(fill = "#1e293b"),
        strip.text       = element_text(colour = "white", face = "bold"))
ggsave(file.path(OUT_DIR, "relative_abundance_top20.pdf"),
       p_relabund, width = 20, height = 8)

cat("\n======================================================\n")
cat("phyloseq analysis complete!\n")
cat("Output directory:", OUT_DIR, "\n")
cat("Files created:\n")
cat("  alpha_by_environment.pdf\n")
cat("  alpha_by_sample_type.pdf\n")
cat("  alpha_by_[pcr_gene].pdf  (one per gene with sufficient groups)\n")
cat("  alpha_stats.csv\n")
cat("  pcoa_environment.pdf\n")
cat("  pcoa_sample_type.pdf\n")
cat("  pcoa_mcr1.pdf\n")
cat("  relative_abundance_top20.pdf\n")
cat("  ps_real_samples.rds\n")
cat("End time:", format(Sys.time()), "\n")
cat("======================================================\n")
