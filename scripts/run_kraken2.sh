#!/bin/bash
#SBATCH --job-name=kraken2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --array=0-62                    # 63 jobs: one per sample/control
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/kraken2_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/kraken2_%A_%a.err

# =============================================================================
# run_kraken2.sh
# Step 3.1 — Taxonomic classification using Kraken2 + Bracken
#
# What does this script do?
#   Takes each sample's host-depleted reads (bacteria only, no chicken/human)
#   and identifies which organisms are present and in what proportions.
#   Both Kraken2 and Bracken run in the same script because Bracken needs
#   Kraken2's output as its input — they always run together.
#
# STEP 1 — Kraken2 (classification):
#   Compares every read against the PlusPF database (79.8 GB of reference
#   genomes from bacteria, archaea, viruses, fungi, protozoa, human, plasmids).
#   For each read it asks "which organism does this DNA come from?" and
#   assigns it to the best matching species. Produces two output files:
#     - report file  : summary of how many reads assigned to each species
#     - output file  : individual assignment for every single read
#
# STEP 2 — Bracken (abundance correction):
#   Kraken2 sometimes assigns reads to a genus or family level when it cannot
#   confidently pick a species (e.g. reads shared between E. coli strains).
#   Bracken takes Kraken2's report and mathematically redistributes these
#   ambiguous reads back down to the species level, giving more accurate
#   abundance estimates. This is the file used for all downstream analysis.
#
# Why 128 GB memory?
#   Kraken2 loads the entire database into RAM for fast searching.
#   The PlusPF database is ~103 GB after extraction, so 128 GB RAM is needed
#   to hold the database plus working memory for the analysis.
#
# Input:  host_depleted/SAMPLE_depleted.fastq.gz
# Output: kraken2_output/SAMPLE_kraken2.report    (Kraken2 report)
#         kraken2_output/SAMPLE_kraken2.output    (per-read assignments)
#         kraken2_output/SAMPLE_bracken.report    (Bracken corrected report)
#         kraken2_output/SAMPLE_bracken.output    (Bracken output)
#
# Run AFTER:
#   1. download_kraken2_db.sh has completed (database ready)
#   2. run_host_depletion.sh has completed (63 depleted files exist)
#
# Verify before submitting:
#   ls host_depleted/ | wc -l         → should be 63
#   ls databases/kraken2/hash.k2d     → database file must exist
# =============================================================================

# ── Load conda and activate taxonomy environment ──────────────────────────────
# The taxonomy environment contains both Kraken2 and Bracken
module load Anaconda3
eval "$(conda shell.bash hook)"   # CRITICAL: enables conda activate in SLURM
conda activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/taxonomy

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
INDIR=$BASE/host_depleted           # input: host-depleted reads
OUTDIR=$BASE/kraken2_output         # output: Kraken2 and Bracken results
DB=$BASE/databases/kraken2          # Kraken2 PlusPF database location

# Create output directory if it does not already exist
mkdir -p "$OUTDIR"

# ── Sample list (61 samples + 2 negative controls) ───────────────────────────
# Order must match all previous scripts exactly
SAMPLES=(
    PK-UF-1002 PK-UF-1009 PK-UF-1015 PK-UF-1022 PK-UF-1025
    PK-UF-1028 PK-UF-1044 PK-UF-1064 PK-UF-1069 PK-UF-1070
    PK-UF-1101 PK-UF-1139 PK-UF-1140 PK-UF-1154 PK-UF-1155
    PK-UF-1162 PK-UF-1180 PK-U3-32EW113 PK-U3-32EE122 PK-U3-42ER128
    PK-U3-42EM130 PK-U3-42EW142 PK-U3-2121P14 PK-U3-2127P15
    PK-U3-1211P40 PK-U3-1235P51 PK-U3-4203P78 PK-U3-4215P84
    PK-U3-1101P2 PK-U3-4105P26 PK-UF-1010 PK-UF-1016 PK-UF-1020
    PK-UF-1034 PK-UF-1040 PK-UF-1055 PK-UF-1059 PK-UF-1075
    PK-UF-1078 PK-UF-1111 PK-UF-1131 PK-UF-1143 PK-UF-1147
    PK-UF-1158 PK-UF-1160 PK-UF-1181 PK-U3-12ES96 PK-U3-12EF97
    PK-U3-12EW99 PK-U3-32EF109 PK-U3-32EE125 PK-U3-1105P3
    PK-U3-2101P9 PK-U3-3101P17 PK-U3-3103P18 PK-U3-3127P23
    PK-U3-1215P42 PK-U3-1221P45 PK-U3-2235P62 PK-U3-3207P67 PK-U3-4211P83
    EB-neg        # index 61 — negative extraction blank control
    AE-buffer     # index 62 — negative AE buffer control
)

# ── Pick this job's sample using the array task index ────────────────────────
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

INPUT="$INDIR/${SAMPLE}_depleted.fastq.gz"

# Output file paths — one set per sample
K2_REPORT="$OUTDIR/${SAMPLE}_kraken2.report"    # Kraken2 summary report
K2_OUTPUT="$OUTDIR/${SAMPLE}_kraken2.output"    # per-read assignments
BR_REPORT="$OUTDIR/${SAMPLE}_bracken.report"    # Bracken corrected report
BR_OUTPUT="$OUTDIR/${SAMPLE}_bracken.output"    # Bracken species abundances

echo "======================================================"
echo "Array task ID : $SLURM_ARRAY_TASK_ID"
echo "Sample        : $SAMPLE"
echo "Input         : $INPUT"
echo "Database      : $DB"
echo "Start time    : $(date)"
echo "======================================================"

# ── Check input file exists ───────────────────────────────────────────────────
if [ ! -f "$INPUT" ]; then
    echo "WARNING: $INPUT not found — skipping $SAMPLE"
    echo "         Check that host depletion completed for this sample"
    exit 0
fi

# Count reads going into classification
READS_IN=$(gunzip -c "$INPUT" | wc -l)
READS_IN=$((READS_IN / 4))
echo "Reads entering classification: $READS_IN"

# ==============================================================================
# STEP 1 — Run Kraken2
#
# --db $DB
#   Path to the PlusPF database. Kraken2 loads this into RAM (~103 GB).
#
# --threads 16
#   Use 16 CPU threads for faster processing.
#
# --report $K2_REPORT
#   The summary report — shows how many reads were assigned to each
#   species, genus, family etc. This is the file Bracken reads as input.
#   Format: percentage | reads covered | reads assigned | rank | taxID | name
#
# --output $K2_OUTPUT
#   Per-read assignments — one line per read showing exactly which species
#   (or "unclassified") each read was assigned to. Large file.
#
# --gzip-compressed
#   Tells Kraken2 the input file is gzip compressed (.fastq.gz)
#   Without this flag Kraken2 cannot read the compressed file.
#
# --minimum-hit-groups 3
#   A read must match at least 3 groups of k-mers in the database before
#   being classified. Reduces false positive classifications.
#   Recommended setting for long ONT reads.
# ==============================================================================
echo ""
echo "------------------------------------------------------"
echo "STEP 1: Running Kraken2 classification..."
echo "------------------------------------------------------"

kraken2 \
    --db "$DB" \
    --threads 16 \
    --report "$K2_REPORT" \
    --output "$K2_OUTPUT" \
    --gzip-compressed \
    --minimum-hit-groups 3 \
    "$INPUT"

echo "Kraken2 complete for $SAMPLE"

# ==============================================================================
# STEP 2 — Run Bracken
#
# Bracken reads Kraken2's report and redistributes reads assigned at higher
# taxonomic levels (genus, family) back down to species level for more
# accurate abundance estimates.
#
# -d $DB
#   Same database path as Kraken2. Bracken uses the database files to
#   understand the taxonomic relationships between organisms.
#
# -i $K2_REPORT
#   Input: the Kraken2 report file produced in Step 1.
#
# -o $BR_OUTPUT
#   Bracken's main output — a table of species with corrected read counts
#   and abundance percentages. This is what you use for downstream analysis.
#
# -w $BR_REPORT
#   A new Kraken2-format report with Bracken's corrected read estimates.
#   Useful for visualisation with Krona.
#
# -r 100
#   Read length parameter. Bracken databases are pre-built for specific
#   read lengths (100, 150, 200 bp). For ONT long reads we use 100 as
#   the minimum available — the classification is still accurate for
#   longer reads.
#
# -l S
#   Taxonomic level to estimate abundance at. S = Species level.
#   Other options: G (genus), F (family), P (phylum), D (domain).
#   Species level gives the most biologically meaningful results.
#
# -t 10
#   Minimum number of reads a species must have in the Kraken2 report
#   before Bracken will estimate its abundance. Filters out species
#   supported by very few reads which are likely noise.
# ==============================================================================
echo ""
echo "------------------------------------------------------"
echo "STEP 2: Running Bracken abundance estimation..."
echo "------------------------------------------------------"

bracken \
    -d "$DB" \
    -i "$K2_REPORT" \
    -o "$BR_OUTPUT" \
    -w "$BR_REPORT" \
    -r 100 \
    -l S \
    -t 10

echo "Bracken complete for $SAMPLE"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "======================================================"
echo "Done: $SAMPLE"
echo "  Kraken2 report : $K2_REPORT"
echo "  Bracken output : $BR_OUTPUT"
echo "  End time       : $(date)"
echo "======================================================"

conda deactivate
