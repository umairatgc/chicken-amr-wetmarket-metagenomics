# =============================================================================
# amr_analysis.R
# AMR Gene Analysis — One Health Metagenomics Project
#
# Analyses results from 5 AMR detection approaches:
#   1. AMRFinderPlus (assembly-based)
#   2. RGI / CARD (assembly-based)
#   3. Abricate — 7 databases (assembly-based)
#   4. ResFinder — assembly-based
#   5. ResFinder — read-based
#
# USAGE:
#   Set RESULTS_DIR below to your local results folder, then:
#   source("amr_analysis.R")  in RStudio / Spyder
#
# REQUIRED PACKAGES (install once):
#   install.packages(c("tidyverse","readxl","writexl","pheatmap","ggplot2",
#                      "RColorBrewer","patchwork","scales","janitor"))
# =============================================================================

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(scales)
library(janitor)

# ── 0. CONFIGURE PATHS ────────────────────────────────────────────────────────
# Set this to the folder where you downloaded all AMR results from ARC
RESULTS_DIR <- "~/OneDrive-Nexus365/A_Files/1_DPhil/6_Metagenome-Analysis/1_Analysis-Files"

AMRFINDER_DIR     <- file.path(RESULTS_DIR, "amrfinder_results")
RGI_DIR           <- file.path(RESULTS_DIR, "rgi_results")
ABRICATE_DIR      <- file.path(RESULTS_DIR, "abricate_results")
RESFINDER_ASM_DIR <- file.path(RESULTS_DIR, "resfinder_assembly_results")
RESFINDER_READ_DIR<- file.path(RESULTS_DIR, "resfinder_results")
OUTPUT_DIR        <- file.path(RESULTS_DIR, "R_outputs")

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ── 1. LOAD AMRFinderPlus RESULTS ─────────────────────────────────────────────

load_amrfinder <- function(dir) {
  files <- list.files(dir, pattern = "amrfinder_results\\.txt$",
                      recursive = TRUE, full.names = TRUE)

  if (length(files) == 0) stop("No AMRFinderPlus results found in: ", dir)

  map_dfr(files, function(f) {
    sample <- basename(dirname(f))
    df <- read_tsv(f, show_col_types = FALSE, comment = "#")
    df$sample <- sample
    df
  }) %>%
    clean_names() %>%
    mutate(sample = as.character(sample))
}

amrfinder <- load_amrfinder(AMRFINDER_DIR)

# Filter to AMR genes only (exclude STRESS, VIRULENCE if not needed)
amrfinder_amr <- amrfinder %>%
  filter(element_type == "AMR")

cat("AMRFinderPlus: ", nrow(amrfinder_amr), "AMR hits across",
    n_distinct(amrfinder_amr$sample), "samples\n")

# ── 2. LOAD RGI / CARD RESULTS ────────────────────────────────────────────────

load_rgi <- function(dir) {
  files <- list.files(dir, pattern = "rgi_output\\.txt$",
                      recursive = TRUE, full.names = TRUE)

  if (length(files) == 0) stop("No RGI results found in: ", dir)

  map_dfr(files, function(f) {
    sample <- basename(dirname(f))
    df <- read_tsv(f, show_col_types = FALSE)
    df$sample <- sample
    df
  }) %>%
    clean_names() %>%
    mutate(sample = as.character(sample))
}

rgi <- load_rgi(RGI_DIR)

# Split by confidence level
rgi_strict <- rgi %>% filter(cut_off %in% c("Perfect", "Strict"))
rgi_loose  <- rgi %>% filter(cut_off == "Loose")

cat("RGI: ", nrow(rgi), "total hits |",
    nrow(rgi_strict), "Strict/Perfect |",
    nrow(rgi_loose), "Loose across",
    n_distinct(rgi$sample), "samples\n")

# ── 3. LOAD ABRICATE RESULTS ──────────────────────────────────────────────────

ABRICATE_DBS <- c("ncbi", "card", "resfinder", "vfdb",
                  "plasmidfinder", "ecoh", "argannot")

load_abricate <- function(dir, databases = ABRICATE_DBS) {
  all_results <- list()

  for (db in databases) {
    files <- list.files(dir, pattern = paste0("abricate_", db, "\\.txt$"),
                        recursive = TRUE, full.names = TRUE)

    if (length(files) == 0) {
      warning("No Abricate results for database: ", db)
      next
    }

    df <- map_dfr(files, function(f) {
      sample <- basename(dirname(f))
      d <- read_tsv(f, show_col_types = FALSE, comment = "#")
      d$sample <- sample
      d
    }) %>%
      clean_names() %>%
      mutate(sample = as.character(sample),
             database_name = db)

    all_results[[db]] <- df
  }

  bind_rows(all_results)
}

abricate <- load_abricate(ABRICATE_DIR)

cat("Abricate: ", nrow(abricate), "total hits across",
    n_distinct(abricate$sample), "samples and",
    n_distinct(abricate$database_name), "databases\n")

# Separate virulence from AMR
abricate_amr <- abricate %>%
  filter(database_name %in% c("ncbi", "card", "resfinder", "argannot"))

abricate_virulence <- abricate %>%
  filter(database_name == "vfdb")

abricate_plasmid <- abricate %>%
  filter(database_name == "plasmidfinder")

# ── 4. LOAD RESFINDER RESULTS (ASSEMBLY-BASED) ────────────────────────────────

load_resfinder <- function(dir, mode = "assembly") {
  files <- list.files(dir, pattern = "ResFinder_results_tab\\.txt$",
                      recursive = TRUE, full.names = TRUE)

  if (length(files) == 0) stop("No ResFinder results found in: ", dir)

  map_dfr(files, function(f) {
    sample <- basename(dirname(f))

    # ResFinder tab file can be empty (no hits) — handle gracefully
    df <- tryCatch(
      read_tsv(f, show_col_types = FALSE),
      error = function(e) tibble()
    )

    if (nrow(df) == 0) return(tibble(sample = sample))

    df$sample  <- sample
    df$mode    <- mode
    df
  }) %>%
    clean_names() %>%
    mutate(sample = as.character(sample))
}

resfinder_asm  <- load_resfinder(RESFINDER_ASM_DIR,  mode = "assembly")
resfinder_read <- load_resfinder(RESFINDER_READ_DIR, mode = "reads")

cat("ResFinder assembly: ", nrow(resfinder_asm),  "hits across",
    n_distinct(resfinder_asm$sample),  "samples\n")
cat("ResFinder reads:    ", nrow(resfinder_read), "hits across",
    n_distinct(resfinder_read$sample), "samples\n")

# Combined ResFinder
resfinder_all <- bind_rows(resfinder_asm, resfinder_read)

# ── 5. PER-SAMPLE SUMMARIES ───────────────────────────────────────────────────

# AMRFinderPlus summary
amrfinder_summary <- amrfinder_amr %>%
  group_by(sample) %>%
  summarise(
    amrfinder_hits        = n(),
    amrfinder_unique_genes = n_distinct(gene_symbol),
    amrfinder_classes     = n_distinct(class),
    .groups = "drop"
  )

# RGI summary
rgi_summary <- rgi %>%
  group_by(sample) %>%
  summarise(
    rgi_total_hits   = n(),
    rgi_strict_hits  = sum(cut_off %in% c("Perfect", "Strict")),
    rgi_loose_hits   = sum(cut_off == "Loose"),
    rgi_unique_genes = n_distinct(aro_term, na.rm = TRUE),
    rgi_drug_classes = n_distinct(drug_class, na.rm = TRUE),
    .groups = "drop"
  )

# Abricate summary (AMR databases only)
abricate_summary <- abricate_amr %>%
  group_by(sample) %>%
  summarise(
    abricate_hits        = n(),
    abricate_unique_genes = n_distinct(gene),
    abricate_databases   = n_distinct(database_name),
    .groups = "drop"
  )

# ResFinder summary
resfinder_summary <- resfinder_all %>%
  filter(!is.na(resistance_gene)) %>%
  group_by(sample, mode) %>%
  summarise(
    hits         = n(),
    unique_genes = n_distinct(resistance_gene, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = mode,
              values_from = c(hits, unique_genes),
              names_prefix = "resfinder_")

# Master summary table
all_samples <- tibble(
  sample = unique(c(amrfinder_amr$sample, rgi$sample,
                    abricate_amr$sample, resfinder_all$sample))
)

master_summary <- all_samples %>%
  left_join(amrfinder_summary,  by = "sample") %>%
  left_join(rgi_summary,        by = "sample") %>%
  left_join(abricate_summary,   by = "sample") %>%
  left_join(resfinder_summary,  by = "sample") %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))

cat("\nMaster summary table: ", nrow(master_summary),
    "samples with", ncol(master_summary), "columns\n")

# Save master summary
write_csv(master_summary,
          file.path(OUTPUT_DIR, "amr_master_summary.csv"))

# ── 6. GENE PRESENCE/ABSENCE MATRIX ──────────────────────────────────────────

# Using RGI Strict/Perfect hits for clean heatmap
make_presence_matrix <- function(df, gene_col, sample_col,
                                 min_samples = 2) {
  df %>%
    select(all_of(c(sample_col, gene_col))) %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = all_of(gene_col),
                values_from = present,
                values_fill = 0) %>%
    column_to_rownames(sample_col) %>%
    # Keep only genes found in >= min_samples
    select(where(~ sum(.) >= min_samples))
}

# RGI strict/perfect presence matrix
rgi_matrix <- make_presence_matrix(
  rgi_strict, "aro_term", "sample", min_samples = 2
)

cat("RGI gene matrix: ", nrow(rgi_matrix), "samples ×",
    ncol(rgi_matrix), "genes (found in ≥2 samples)\n")

# AMRFinderPlus presence matrix
amr_matrix <- make_presence_matrix(
  amrfinder_amr, "gene_symbol", "sample", min_samples = 2
)

# ResFinder (read-based) presence matrix
resfinder_read_clean <- resfinder_read %>%
  filter(!is.na(resistance_gene))

resfinder_matrix <- make_presence_matrix(
  resfinder_read_clean, "resistance_gene", "sample", min_samples = 2
)

# ── 7. VISUALISATIONS ─────────────────────────────────────────────────────────

# ── 7a. Heatmap: RGI gene presence/absence ────────────────────────────────────
if (ncol(rgi_matrix) > 0 && nrow(rgi_matrix) > 1) {

  pdf(file.path(OUTPUT_DIR, "heatmap_rgi_genes.pdf"),
      width = max(12, ncol(rgi_matrix) * 0.3 + 4),
      height = max(10, nrow(rgi_matrix) * 0.25 + 4))

  pheatmap(
    t(rgi_matrix),
    color          = c("white", "#1a7abf"),
    border_color   = "grey90",
    cluster_rows   = TRUE,
    cluster_cols   = TRUE,
    show_rownames  = TRUE,
    show_colnames  = TRUE,
    fontsize_row   = 7,
    fontsize_col   = 7,
    main           = "RGI/CARD — Gene Presence/Absence (Strict + Perfect)",
    legend_breaks  = c(0, 1),
    legend_labels  = c("Absent", "Present"),
    angle_col      = 45
  )

  dev.off()
  cat("Saved: heatmap_rgi_genes.pdf\n")
}

# ── 7b. Heatmap: AMRFinderPlus gene presence/absence ─────────────────────────
if (ncol(amr_matrix) > 0 && nrow(amr_matrix) > 1) {

  pdf(file.path(OUTPUT_DIR, "heatmap_amrfinder_genes.pdf"),
      width  = max(12, ncol(amr_matrix) * 0.3 + 4),
      height = max(10, nrow(amr_matrix) * 0.25 + 4))

  pheatmap(
    t(amr_matrix),
    color          = c("white", "#22863a"),
    border_color   = "grey90",
    cluster_rows   = TRUE,
    cluster_cols   = TRUE,
    show_rownames  = TRUE,
    show_colnames  = TRUE,
    fontsize_row   = 7,
    fontsize_col   = 7,
    main           = "AMRFinderPlus — Gene Presence/Absence (AMR genes only)",
    legend_breaks  = c(0, 1),
    legend_labels  = c("Absent", "Present"),
    angle_col      = 45
  )

  dev.off()
  cat("Saved: heatmap_amrfinder_genes.pdf\n")
}

# ── 7c. Bar chart: AMR gene count per sample (RGI strict) ────────────────────
p_rgi_bar <- rgi_summary %>%
  pivot_longer(c(rgi_strict_hits, rgi_loose_hits),
               names_to = "confidence", values_to = "hits") %>%
  mutate(confidence = recode(confidence,
                             rgi_strict_hits = "Strict/Perfect",
                             rgi_loose_hits  = "Loose")) %>%
  ggplot(aes(x = reorder(sample, rgi_total_hits), y = hits,
             fill = confidence)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Strict/Perfect" = "#1a7abf",
                                "Loose"          = "#a8d4f0")) +
  coord_flip() +
  labs(title   = "RGI/CARD — AMR Hits per Sample",
       subtitle = "Strict/Perfect and Loose confidence hits",
       x = "Sample", y = "Number of Hits", fill = "Confidence") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(OUTPUT_DIR, "barplot_rgi_per_sample.pdf"),
       p_rgi_bar, width = 10, height = 14)
cat("Saved: barplot_rgi_per_sample.pdf\n")

# ── 7d. Drug class distribution (RGI Strict) ─────────────────────────────────
rgi_drug_counts <- rgi_strict %>%
  filter(!is.na(drug_class)) %>%
  # Split multi-class entries (separated by semicolon)
  mutate(drug_class = str_split(drug_class, ";\\s*")) %>%
  unnest(drug_class) %>%
  mutate(drug_class = str_trim(drug_class)) %>%
  count(drug_class, sort = TRUE) %>%
  slice_head(n = 20)

p_drug_class <- ggplot(rgi_drug_counts,
                       aes(x = reorder(drug_class, n), y = n)) +
  geom_bar(stat = "identity", fill = "#1a7abf") +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(title = "RGI/CARD — Top 20 Drug Classes (Strict/Perfect hits)",
       x = "Drug Class", y = "Number of Hits") +
  theme_minimal(base_size = 11) +
  expand_limits(y = max(rgi_drug_counts$n) * 1.1)

ggsave(file.path(OUTPUT_DIR, "barplot_rgi_drug_classes.pdf"),
       p_drug_class, width = 10, height = 8)
cat("Saved: barplot_rgi_drug_classes.pdf\n")

# ── 7e. Abricate: hits per database ──────────────────────────────────────────
abricate_db_summary <- abricate %>%
  count(database_name, sort = TRUE)

p_abricate_db <- ggplot(abricate_db_summary,
                        aes(x = reorder(database_name, n),
                            y = n, fill = database_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = comma(n)), hjust = -0.2, size = 3.5) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Abricate — Total Hits by Database (all 61 samples)",
       x = "Database", y = "Total Hits") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none") +
  expand_limits(y = max(abricate_db_summary$n) * 1.12)

ggsave(file.path(OUTPUT_DIR, "barplot_abricate_by_database.pdf"),
       p_abricate_db, width = 8, height = 5)
cat("Saved: barplot_abricate_by_database.pdf\n")

# ── 7f. ResFinder: read-based vs assembly comparison (shared samples) ─────────
resfinder_compare <- resfinder_all %>%
  filter(!is.na(resistance_gene)) %>%
  group_by(sample, mode) %>%
  summarise(hits = n(), .groups = "drop") %>%
  pivot_wider(names_from = mode, values_from = hits, values_fill = 0)

if ("assembly" %in% colnames(resfinder_compare) &&
    "reads"    %in% colnames(resfinder_compare)) {

  p_compare <- ggplot(resfinder_compare,
                      aes(x = assembly, y = reads, label = sample)) +
    geom_point(colour = "#1a7abf", alpha = 0.7, size = 2.5) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey50") +
    labs(title    = "ResFinder: Assembly-based vs Read-based hits per sample",
         subtitle = "Points above the line = read-based detects more; below = assembly detects more",
         x = "Assembly-based hits",
         y = "Read-based hits") +
    theme_minimal(base_size = 11)

  ggsave(file.path(OUTPUT_DIR, "scatter_resfinder_asm_vs_read.pdf"),
         p_compare, width = 7, height = 6)
  cat("Saved: scatter_resfinder_asm_vs_read.pdf\n")
}

# ── 7g. AMRFinderPlus: resistance class breakdown ─────────────────────────────
amrfinder_class <- amrfinder_amr %>%
  filter(!is.na(class)) %>%
  count(class, sort = TRUE) %>%
  slice_head(n = 20)

p_amrfinder_class <- ggplot(amrfinder_class,
                             aes(x = reorder(class, n), y = n)) +
  geom_bar(stat = "identity", fill = "#22863a") +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(title = "AMRFinderPlus — Top 20 Resistance Classes",
       x = "Resistance Class", y = "Total Hits") +
  theme_minimal(base_size = 11) +
  expand_limits(y = max(amrfinder_class$n) * 1.12)

ggsave(file.path(OUTPUT_DIR, "barplot_amrfinder_classes.pdf"),
       p_amrfinder_class, width = 10, height = 8)
cat("Saved: barplot_amrfinder_classes.pdf\n")

# ── 8. TOOL AGREEMENT — GENE-LEVEL OVERLAP ───────────────────────────────────
# Which resistance gene families are detected by multiple tools?

# Standardise gene names as best as possible for cross-tool comparison
rgi_genes       <- rgi_strict %>%
                     distinct(sample, gene = aro_term) %>%
                     mutate(tool = "RGI_strict")

amrfinder_genes <- amrfinder_amr %>%
                     distinct(sample, gene = gene_symbol) %>%
                     mutate(tool = "AMRFinderPlus")

resfinder_genes <- resfinder_read %>%
                     filter(!is.na(resistance_gene)) %>%
                     distinct(sample, gene = resistance_gene) %>%
                     mutate(tool = "ResFinder_reads")

abricate_amr_genes <- abricate_amr %>%
                        distinct(sample, gene) %>%
                        mutate(tool = "Abricate_AMR")

all_genes_long <- bind_rows(
  rgi_genes, amrfinder_genes, resfinder_genes, abricate_amr_genes
)

# Count how many unique tools detect each gene (across all samples)
gene_tool_overlap <- all_genes_long %>%
  group_by(gene) %>%
  summarise(
    n_tools   = n_distinct(tool),
    tools     = paste(sort(unique(tool)), collapse = " | "),
    n_samples = n_distinct(sample),
    .groups   = "drop"
  ) %>%
  arrange(desc(n_tools), desc(n_samples))

write_csv(gene_tool_overlap,
          file.path(OUTPUT_DIR, "gene_tool_overlap.csv"))

cat("\nTop genes detected by multiple tools:\n")
print(head(gene_tool_overlap %>% filter(n_tools > 1), 20))

# ── 9. EXPORT ALL CLEANED DATA ───────────────────────────────────────────────
write_csv(amrfinder,       file.path(OUTPUT_DIR, "amrfinder_all.csv"))
write_csv(rgi,             file.path(OUTPUT_DIR, "rgi_all.csv"))
write_csv(abricate,        file.path(OUTPUT_DIR, "abricate_all.csv"))
write_csv(resfinder_asm,   file.path(OUTPUT_DIR, "resfinder_assembly.csv"))
write_csv(resfinder_read,  file.path(OUTPUT_DIR, "resfinder_reads.csv"))
write_csv(gene_tool_overlap, file.path(OUTPUT_DIR, "gene_tool_overlap.csv"))

cat("\n")
cat("=============================================================\n")
cat("Analysis complete. Outputs saved to:\n")
cat(" ", OUTPUT_DIR, "\n")
cat("Files created:\n")
cat("  amr_master_summary.csv       — per-sample hit counts all tools\n")
cat("  amrfinder_all.csv            — all AMRFinderPlus results\n")
cat("  rgi_all.csv                  — all RGI results\n")
cat("  abricate_all.csv             — all Abricate results (7 DBs)\n")
cat("  resfinder_assembly.csv       — ResFinder assembly-based\n")
cat("  resfinder_reads.csv          — ResFinder read-based\n")
cat("  gene_tool_overlap.csv        — cross-tool gene agreement\n")
cat("  heatmap_rgi_genes.pdf        — RGI gene presence/absence heatmap\n")
cat("  heatmap_amrfinder_genes.pdf  — AMRFinderPlus gene heatmap\n")
cat("  barplot_rgi_per_sample.pdf   — RGI hits per sample\n")
cat("  barplot_rgi_drug_classes.pdf — RGI drug class distribution\n")
cat("  barplot_abricate_by_database.pdf — Abricate hits by database\n")
cat("  barplot_amrfinder_classes.pdf    — AMRFinderPlus class breakdown\n")
cat("  scatter_resfinder_asm_vs_read.pdf — assembly vs read comparison\n")
cat("=============================================================\n")
