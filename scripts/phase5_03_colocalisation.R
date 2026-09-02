# =============================================================================
# phase5_03_colocalisation.R
# Phase 5 — Step 3: Co-localisation analysis
#
# For each ARG hit: find MGEs on the same contig within 10 kb.
# Join with MAG taxonomy (GTDB-Tk) and sample metadata.
#
# INPUT:
#   phase5/args_combined.tsv
#   phase5/mges_combined.tsv
#   phase5/contig_mag_map.tsv
#   sample_metadata.csv  (built from Excel file)
#
# OUTPUT:
#   phase5/colocalisation_master.tsv   — one row per ARG hit, with MGE context
#   phase5/coloc_summary.tsv           — per-sample summary statistics
# =============================================================================

library(tidyverse)

BASE      <- "/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
OUTDIR    <- file.path(BASE, "phase5")
THRESHOLD <- 10000L   # 10 kb co-localisation threshold

# ── Load data ────────────────────────────────────────────────────────────────
cat("Loading ARGs...\n")
args <- read_tsv(file.path(OUTDIR, "args_combined.tsv"), show_col_types = FALSE)

cat("Loading MGEs...\n")
mges <- read_tsv(file.path(OUTDIR, "mges_combined.tsv"), show_col_types = FALSE)

cat("Loading contig → MAG map...\n")
mag_map <- read_tsv(file.path(OUTDIR, "contig_mag_map.tsv"), show_col_types = FALSE)

cat("Loading sample metadata...\n")
# Adjust path to wherever sample_metadata.csv was saved
metadata <- read_csv(file.path(OUTDIR, "sample_metadata.csv"), show_col_types = FALSE)

cat("  ARG hits:   ", nrow(args), "\n")
cat("  MGE records:", nrow(mges), "\n")
cat("  MAG contigs:", nrow(mag_map), "\n\n")

# ── Step 1: Find co-localised ARG–MGE pairs ───────────────────────────────────
# Join on sample + contig, then calculate distance between features
coloc_pairs <- args %>%
  inner_join(
    mges %>% rename(mge_start = start, mge_end = end),
    by = c("sample", "contig"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    # Distance = 0 if overlapping, otherwise gap between features
    distance = pmax(0L,
      pmax(start, mge_start) - pmin(end, mge_end)
    ),
    colocalised = distance <= THRESHOLD
  ) %>%
  # For each ARG, keep the closest MGE of each type
  group_by(sample, contig, gene, start, end, mge_type) %>%
  slice_min(distance, n = 1, with_ties = FALSE) %>%
  ungroup()

# ── Step 2: Summarise per ARG — closest MGE overall ─────────────────────────
closest_mge <- coloc_pairs %>%
  group_by(sample, contig, gene, start, end) %>%
  slice_min(distance, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(sample, contig, gene, start, end,
         closest_mge_type = mge_type,
         closest_mge_name = mge_name,
         closest_mge_distance = distance,
         colocalised,
         mge_tool = tool.y)

# ── Step 3: Join back with all ARG hits ──────────────────────────────────────
# ARGs with no MGE on same contig get NA for MGE columns + colocalised = FALSE
args_with_mge <- args %>%
  left_join(closest_mge,
            by = c("sample", "contig", "gene", "start", "end")) %>%
  mutate(
    colocalised           = replace_na(colocalised, FALSE),
    mge_on_same_contig    = !is.na(closest_mge_type)
  )

# ── Step 4: Add MAG taxonomy ─────────────────────────────────────────────────
args_taxon <- args_with_mge %>%
  left_join(mag_map %>% select(sample, contig, bin, completeness,
                                contamination, taxonomy),
            by = c("sample", "contig")) %>%
  mutate(
    in_quality_mag = !is.na(bin),
    # Parse taxonomy components
    phylum  = str_extract(taxonomy, "(?<=p__)[^;]+"),
    class   = str_extract(taxonomy, "(?<=c__)[^;]+"),
    order   = str_extract(taxonomy, "(?<=o__)[^;]+"),
    family  = str_extract(taxonomy, "(?<=f__)[^;]+"),
    genus   = str_extract(taxonomy, "(?<=g__)[^;]+"),
    species = str_extract(taxonomy, "(?<=s__)[^;]+")
  )

# ── Step 5: Add sample metadata ──────────────────────────────────────────────
# Rename CSV columns to match expected names in script
metadata_clean <- metadata %>%
  rename(
    sample_id    = `Sample ID`,
    sample_type  = `Sample Type`,
    sample_type2 = `Sample Type 2`,
    source       = `Source 1`,
    ont_run      = `ONT Run`,
    mcr1         = `mcr-1 C2`,
    tetX4        = `tetX4 C2`,
    fosA3        = `fosA3 C2`,
    NDM          = `NDM C2`,
    OXA48        = `OXA-48 C2`,
    KPC          = `KPC C2`,
    pos_count    = `Pos Count`
  )

# Join ARG hits with cleaned metadata using sample ID
args_final <- args_taxon %>%
  left_join(metadata_clean %>%
              select(sample_id, sample_type, sample_type2, source,
                     ont_run, mcr1, tetX4, fosA3, NDM, OXA48, KPC, pos_count),
            by = c("sample" = "sample_id")) %>%
  mutate(
    # Classify samples into broad environment categories
    environment = case_when(
      str_detect(sample_type, "Butcher") ~ "Butcher shop",
      str_detect(sample_type, "Farm")    ~ "Farm",
      TRUE                                ~ "Other"
    ),
    # Extract day number for longitudinal Farm Litter samples (e.g. D1, D14, D35)
    litter_day = as.integer(str_extract(sample_type2, "(?<=D)\\d+"))
  )

# ── Save master table ─────────────────────────────────────────────────────────
write_tsv(args_final, file.path(OUTDIR, "colocalisation_master.tsv"))

# ── Per-sample summary ────────────────────────────────────────────────────────
sample_summary <- args_final %>%
  group_by(sample, sample_type, environment, source) %>%
  summarise(
    total_arg_hits      = n(),
    unique_genes        = n_distinct(gene),
    unique_drug_classes = n_distinct(drug_class, na.rm = TRUE),
    args_in_mag         = sum(in_quality_mag),
    args_with_mge_contig = sum(mge_on_same_contig),
    args_colocalised_10kb = sum(colocalised),
    pct_colocalised     = round(mean(colocalised) * 100, 1),
    n_is_elements       = sum(closest_mge_type == "IS_element", na.rm = TRUE),
    n_integrons         = sum(closest_mge_type == "Integron", na.rm = TRUE),
    n_plasmids          = sum(closest_mge_type == "Plasmid", na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(sample_summary, file.path(OUTDIR, "coloc_summary.tsv"))

# ── Print summary statistics ──────────────────────────────────────────────────
cat("\n==========================================\n")
cat("CO-LOCALISATION ANALYSIS COMPLETE\n")
cat("==========================================\n")
cat("Total ARG hits analysed:         ", nrow(args_final), "\n")
cat("ARGs on MAG contigs:             ", sum(args_final$in_quality_mag), "\n")
cat("ARGs with MGE on same contig:    ", sum(args_final$mge_on_same_contig), "\n")
cat("ARGs co-localised ≤10 kb:        ", sum(args_final$colocalised), "\n")
cat("Overall % co-localised:          ",
    round(mean(args_final$colocalised) * 100, 1), "%\n\n")

cat("Co-localisation by MGE type:\n")
print(args_final %>%
  filter(colocalised) %>%
  count(closest_mge_type) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  arrange(desc(n)))

cat("\nCo-localisation by environment:\n")
print(args_final %>%
  group_by(environment) %>%
  summarise(
    arg_hits    = n(),
    colocalised = sum(colocalised),
    pct         = round(mean(colocalised) * 100, 1),
    .groups = "drop"
  ))

cat("\nOutputs:\n")
cat("  ", file.path(OUTDIR, "colocalisation_master.tsv"), "\n")
cat("  ", file.path(OUTDIR, "coloc_summary.tsv"), "\n")
cat("==========================================\n")
