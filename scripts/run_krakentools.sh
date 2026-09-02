#!/bin/bash
#SBATCH --job-name=krakentools
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/krakentools_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/krakentools_%j.err

# =============================================================================
# run_krakentools.sh
# Step 3.2 (Part 1) — Merge all 63 Bracken outputs into one combined matrix
#
# WHAT IS KRAKENTOOLS?
#   KrakenTools is a collection of Python helper scripts written by the same
#   team that built Kraken2. Kraken2 only classifies reads — it does not
#   combine or reshape results for downstream analysis. KrakenTools fills
#   that gap. This script uses one specific KrakenTools script:
#
#     combine_bracken_outputs.py
#
#   That is the only command that runs here. The job name "krakentools" and
#   the filename "run_krakentools.sh" both refer to this tool.
#
# WHAT DOES combine_bracken_outputs.py DO?
#   After Kraken2 + Bracken ran for all 63 samples, you have 63 separate
#   Bracken output files — one per sample. Each file contains a table of
#   species with their read counts for THAT sample only.
#
#   combine_bracken_outputs.py merges all 63 files into ONE matrix where:
#     - Each ROW    is a species
#     - Each COLUMN is a sample
#     - Each CELL   is the read count for that species in that sample
#
#   This combined matrix is the input for ALL downstream R analysis:
#     - decontam  : uses read counts to identify contaminant species
#     - phyloseq  : uses the matrix to calculate alpha and beta diversity
#     - vegan     : uses the matrix for PERMANOVA and ordination
#
# HOW DOES IT FIND THE SAMPLE NAMES?
#   The script automatically scans the kraken2_output/ folder for all files
#   ending in _bracken.output — no hardcoded sample list needed.
#   It strips the _bracken.output suffix from each filename to get the
#   sample name. This means it works even if some samples failed — missing
#   files are simply not included rather than causing an error.
#
# Input:  kraken2_output/*_bracken.output    (63 individual Bracken files)
# Output: kraken2_output/combined_bracken_all_samples.txt   (one merged matrix)
#
# Run AFTER: all 63 run_kraken2.sh array jobs have completed
# Verify:  ls kraken2_output/*_bracken.output | wc -l   → should be 63
# =============================================================================

# ── Load conda and activate taxonomy environment ──────────────────────────────
# KrakenTools (combine_bracken_outputs.py) is installed in the taxonomy env
module load Anaconda3
eval "$(conda shell.bash hook)"   # CRITICAL: enables conda activate in SLURM
conda activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/taxonomy

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
INDIR=$BASE/kraken2_output        # folder containing all bracken output files
OUTPUT=$INDIR/combined_bracken_all_samples.txt  # merged output matrix

echo "======================================================"
echo "KrakenTools — Combining all Bracken outputs"
echo "Start time: $(date)"
echo "Scanning: $INDIR"
echo "======================================================"

# ── Find all Bracken output files ─────────────────────────────────────────────
# Automatically collects all files ending in _bracken.output
FILES=("$INDIR"/*_bracken.output)

# Check that files were found
if [ ${#FILES[@]} -eq 0 ]; then
    echo "ERROR: No *_bracken.output files found in $INDIR"
    echo "       Check that run_kraken2.sh has completed successfully"
    exit 1
fi

echo "Found ${#FILES[@]} Bracken output files"

# ── Build sample names list ────────────────────────────────────────────────────
# combine_bracken_outputs.py needs sample names as a comma-separated list
# We extract names by stripping the directory path and _bracken.output suffix
# Example: /path/to/kraken2_output/PK-UF-1002_bracken.output → PK-UF-1002
#
# basename strips the directory path
# The second argument to basename strips the suffix
# tr '\n' ',' joins names with commas
# sed 's/,$//' removes the trailing comma from the last name
NAMES=$(for f in "${FILES[@]}"; do
    basename "$f" _bracken.output
done | tr '\n' ',' | sed 's/,$//')

echo ""
echo "Sample names extracted: $(echo "$NAMES" | tr ',' '\n' | wc -l) samples"
echo "First few: $(echo "$NAMES" | cut -d',' -f1-3)..."

# ==============================================================================
# Run combine_bracken_outputs.py
#
# --files  : list of all bracken output files to combine
#            The shell expands *_bracken.output into the full file list
#
# --names  : comma-separated sample names — these become the column headers
#            in the output matrix. Order must match the order of --files.
#            We extract names from filenames so the order matches automatically.
#
# -o       : output file path for the combined matrix
#
# Output format:
#   name              | taxid | taxonomy_lvl | PK-UF-1002_num | PK-UF-1002_frac | ...
#   Escherichia coli  | 562   | S            | 1200           | 0.0045          | ...
#
#   _num  = raw read count for that species in that sample
#   _frac = fraction of classified reads assigned to that species
# ==============================================================================
echo ""
echo "------------------------------------------------------"
echo "Running combine_bracken_outputs.py..."
echo "------------------------------------------------------"

combine_bracken_outputs.py \
    --files "${FILES[@]}" \
    --names "$NAMES" \
    -o "$OUTPUT"

# ── Verify output ─────────────────────────────────────────────────────────────
if [ -f "$OUTPUT" ]; then
    NROWS=$(wc -l < "$OUTPUT")
    NCOLS=$(head -1 "$OUTPUT" | tr '\t' '\n' | wc -l)
    echo ""
    echo "======================================================"
    echo "KrakenTools complete!"
    echo "End time: $(date)"
    echo ""
    echo "Output: $OUTPUT"
    echo "  Rows (species + header): $NROWS"
    echo "  Columns (name+taxid+level + 2 per sample): $NCOLS"
    echo ""
    echo "Preview (first 3 rows, first 6 columns):"
    cut -f1-6 "$OUTPUT" | head -3
    echo ""
    echo "Next steps:"
    echo "  1. Run run_krona.sh for interactive visualisation"
    echo "  2. Use combined_bracken_all_samples.txt as input for R analysis"
    echo "     (decontam, phyloseq, vegan)"
    echo "======================================================"
else
    echo "ERROR: Output file was not created — check the error log"
    exit 1
fi

conda deactivate
