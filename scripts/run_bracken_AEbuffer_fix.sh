#!/bin/bash
#SBATCH --job-name=bracken_AEbuffer_fix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/bracken_AEbuffer_fix_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/bracken_AEbuffer_fix_%j.err

# =============================================================================
# run_bracken_AEbuffer_fix.sh
# One-off fix — Re-run Bracken for the AE-buffer negative control sample
#
# WHY DOES THIS SCRIPT EXIST?
#   During the main Bracken array job (run_kraken2.sh), AE-buffer failed with:
#
#     "Error: no reads found. Please check your Kraken report"
#
#   This happened because Bracken's default minimum read threshold (-t 10)
#   rejected AE-buffer, which only had 26 sequences (23 classified).
#   That is expected for a negative extraction control (buffer only, no sample).
#
#   However, AE-buffer MUST be present in the combined Bracken matrix because:
#     - decontam uses negative controls to identify contaminant species
#     - Without AE-buffer in the matrix, decontam has no reference and fails
#     - phyloseq alpha/beta diversity should also include the control
#
# THE FIX:
#   Re-run Bracken for AE-buffer only, with -t 0 (threshold = 0).
#   This removes the minimum read requirement and forces Bracken to write
#   output even when read counts are very low. The resulting file will have
#   low counts — that is correct and expected for a negative control.
#
# WHAT CHANGES vs. THE ORIGINAL RUN:
#   Two parameters differ from the original Bracken command:
#     -t 10  →  -t 0   forces output despite low read count
#     -w               explicitly sets report filename to match all other samples
#                      (without -w, Bracken auto-names it using the input filename,
#                      producing AE-buffer_kraken2_bracken_species.report instead
#                      of the expected AE-buffer_bracken.report)
#
# Input:  kraken2_output/AE-buffer_kraken2.report   (already exists)
# Output: kraken2_output/AE-buffer_bracken.output   (was missing, now created)
#         kraken2_output/AE-buffer_bracken.report   (consistent naming with all other samples)
#
# Run AFTER: confirming AE-buffer_kraken2.report exists
# Run BEFORE: run_krakentools.sh (which needs all 63 bracken files)
#
# Verify after: ls kraken2_output/*_bracken.output | wc -l  → should be 63
# =============================================================================

# ── Load conda and activate taxonomy environment ──────────────────────────────
module load Anaconda3
eval "$(conda shell.bash hook)"   # CRITICAL: enables conda activate in SLURM
conda activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/taxonomy

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
DB=$BASE/databases/kraken2        # Kraken2 PlusPF database — same as run_kraken2.sh
INDIR=$BASE/kraken2_output

SAMPLE="AE-buffer"
REPORT="$INDIR/${SAMPLE}_kraken2.report"
OUTPUT="$INDIR/${SAMPLE}_bracken.output"
BREPORT="$INDIR/${SAMPLE}_bracken.report"   # explicit report name to match all other samples

echo "======================================================"
echo "Bracken — AE-buffer negative control fix"
echo "Start time: $(date)"
echo "Sample: $SAMPLE"
echo "Database: $DB"
echo "======================================================"

# ── Confirm the database exists ───────────────────────────────────────────────
if [ ! -d "$DB" ]; then
    echo "ERROR: Kraken2 database not found: $DB"
    exit 1
fi

# ── Confirm the Kraken2 report exists ────────────────────────────────────────
if [ ! -f "$REPORT" ]; then
    echo "ERROR: Kraken2 report not found: $REPORT"
    echo "       Check that run_kraken2.sh completed for this sample."
    exit 1
fi

echo ""
echo "Kraken2 report found: $REPORT"
echo "Read count in report:"
awk 'NR==1{print "  Total reads: " $2}' "$REPORT"
echo ""

# ==============================================================================
# Run Bracken with threshold = 0 and explicit report filename
#
# -d : path to the Kraken2 database (same database used in run_kraken2.sh)
# -i : input Kraken2 report for this sample
# -o : output Bracken file (will be picked up by combine_bracken_outputs.py)
# -w : output report file — explicitly named to match all other samples
#      Without -w, Bracken auto-generates the report name from the input
#      filename, producing AE-buffer_kraken2_bracken_species.report which
#      is inconsistent with all other samples that have _bracken.report
# -r : read length (150 bp, same as the original array job)
# -l : taxonomic level (S = species, same as the original array job)
# -t : minimum read threshold set to 0
#      Default is 10 — AE-buffer only has 23 classified reads, so it fails.
#      Setting -t 0 removes this restriction and forces output to be written.
#      Low counts in the output are CORRECT for a negative control.
# ==============================================================================
echo "------------------------------------------------------"
echo "Running Bracken with -t 0 (no minimum read threshold)"
echo "------------------------------------------------------"

bracken \
    -d "$DB" \
    -i "$REPORT" \
    -o "$OUTPUT" \
    -w "$BREPORT" \
    -r 150 \
    -l S \
    -t 0

# ── Verify output ─────────────────────────────────────────────────────────────
if [ -f "$OUTPUT" ]; then
    NLINES=$(wc -l < "$OUTPUT")
    echo ""
    echo "======================================================"
    echo "Bracken complete for AE-buffer!"
    echo "End time: $(date)"
    echo ""
    echo "Output files created:"
    echo "  $OUTPUT"
    echo "  $BREPORT"
    echo ""
    echo "Lines in bracken output (species + header): $NLINES"
    echo ""
    echo "Preview:"
    head -3 "$OUTPUT"
    echo ""
    echo "Total bracken files now available:"
    ls "$INDIR"/*_bracken.output | wc -l
    echo "(should be 63 — all samples including AE-buffer)"
    echo ""
    echo "Next step: sbatch run_krakentools.sh"
    echo "======================================================"
else
    echo "ERROR: Bracken output was not created — check the error log"
    exit 1
fi

conda deactivate
