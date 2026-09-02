#!/bin/bash
#SBATCH --job-name=krona
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/krona_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/krona_%j.err

# =============================================================================
# run_krona.sh
# Step 3.2 (Part 2) — Interactive taxonomic visualisation using Krona
#
# What does Krona produce?
#   A single interactive HTML file that visualises the taxonomic composition
#   of ALL 63 samples together. When you open it in your browser you see:
#     - A multi-level pie chart (kingdom → phylum → class → order → family
#       → genus → species) where you can click any slice to zoom in
#     - A dropdown menu to switch between samples and compare them
#     - Percentage labels showing how much of each sample is each organism
#
#   This is the fastest way to get an overview of what organisms are present
#   and whether anything looks unusual before running formal statistics in R.
#
# How does it work — two steps:
#
# STEP 1 — kreport2krona.py (KrakenTools)
#   Krona cannot read Kraken2 reports directly — the formats are different.
#   kreport2krona.py converts each Kraken2 report into Krona's input format.
#   This runs in a loop, once per sample, producing one _krona.txt per sample.
#
# STEP 2 — ktImportText (Krona)
#   Reads all 63 _krona.txt files at once and combines them into one HTML.
#   Each sample appears as a separate dataset you can select in the browser.
#   The "name" you give each file (sample ID) becomes the label in the dropdown.
#
# Input:  kraken2_output/*_kraken2.report    (63 Kraken2 reports)
# Output: krona_output/*_krona.txt           (intermediate conversion files)
#         krona_output/krona_all_samples.html (final interactive HTML — open this)
#
# Run AFTER: all 63 run_kraken2.sh array jobs have completed
# Verify:  ls kraken2_output/*_kraken2.report | wc -l   → should be 63
# =============================================================================

# ── Load conda and activate taxonomy environment ──────────────────────────────
# Both KrakenTools (kreport2krona.py) and Krona (ktImportText) are in taxonomy env
module load Anaconda3
eval "$(conda shell.bash hook)"   # CRITICAL: enables conda activate in SLURM
conda activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/taxonomy

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
K2DIR=$BASE/kraken2_output        # folder containing Kraken2 reports
OUTDIR=$BASE/krona_output         # folder for Krona intermediate and final files

# Create output directory if it does not already exist
mkdir -p "$OUTDIR"

echo "======================================================"
echo "Krona — Interactive taxonomic visualisation"
echo "Start time: $(date)"
echo "Reading Kraken2 reports from: $K2DIR"
echo "Output directory: $OUTDIR"
echo "======================================================"

# ── Check Kraken2 reports exist ───────────────────────────────────────────────
REPORTS=("$K2DIR"/*_kraken2.report)
if [ ${#REPORTS[@]} -eq 0 ]; then
    echo "ERROR: No *_kraken2.report files found in $K2DIR"
    echo "       Check that run_kraken2.sh has completed successfully"
    exit 1
fi
echo "Found ${#REPORTS[@]} Kraken2 reports"

# ==============================================================================
# STEP 1 — Convert each Kraken2 report to Krona format
#
# kreport2krona.py reads a Kraken2 report and outputs a two-column text file
# that Krona understands:
#   Column 1: read count
#   Column 2: full taxonomic path (e.g. Bacteria;Proteobacteria;...;E.coli)
#
# -r : input Kraken2 report file
# -o : output Krona-format text file
#
# This loop runs once per sample. Each sample gets its own _krona.txt file.
# ==============================================================================
echo ""
echo "------------------------------------------------------"
echo "STEP 1: Converting Kraken2 reports to Krona format..."
echo "------------------------------------------------------"

KRONA_FILES=()   # will collect all _krona.txt paths for Step 2

for REPORT in "${REPORTS[@]}"; do
    # Extract sample name by stripping path and _kraken2.report suffix
    # Example: /path/PK-UF-1002_kraken2.report → PK-UF-1002
    SAMPLE=$(basename "$REPORT" _kraken2.report)
    KRONA_TXT="$OUTDIR/${SAMPLE}_krona.txt"

    echo "  Converting: $SAMPLE"

    kreport2krona.py \
        -r "$REPORT" \
        -o "$KRONA_TXT"

    # Add this file to the list for Step 2
    # Format: filepath,samplename — Krona uses the name as the dropdown label
    KRONA_FILES+=("$KRONA_TXT,$SAMPLE")
done

echo "Conversion complete — ${#REPORTS[@]} files converted"

# ==============================================================================
# STEP 2 — Build interactive HTML with ktImportText
#
# ktImportText reads all the Krona-format text files and combines them into
# one interactive HTML file. Each file becomes one dataset (one pie chart)
# accessible via the dropdown menu in the browser.
#
# Input format for each file: filepath,SampleLabel
#   The SampleLabel (after the comma) is what appears in the dropdown menu.
#   We use the sample ID (e.g. PK-UF-1002) as the label.
#
# -o : output HTML file path
#
# The "${KRONA_FILES[@]}" expands to all the filepath,name pairs we built
# in the loop above — one per sample.
# ==============================================================================
echo ""
echo "------------------------------------------------------"
echo "STEP 2: Building interactive Krona HTML..."
echo "------------------------------------------------------"

HTML_OUTPUT="$OUTDIR/krona_all_samples.html"

ktImportText \
    "${KRONA_FILES[@]}" \
    -o "$HTML_OUTPUT"

# ── Verify output and report ──────────────────────────────────────────────────
if [ -f "$HTML_OUTPUT" ]; then
    echo ""
    echo "======================================================"
    echo "Krona complete!"
    echo "End time: $(date)"
    echo ""
    echo "Interactive HTML: $HTML_OUTPUT"
    echo "File size: $(du -sh "$HTML_OUTPUT" | cut -f1)"
    echo ""
    echo "Download and open in your browser:"
    echo "  scp YOUR_USERNAME@arc-login.arc.ox.ac.uk:$HTML_OUTPUT ~/Desktop/"
    echo ""
    echo "Intermediate files (can delete after checking HTML):"
    echo "  $OUTDIR/*_krona.txt"
    echo "======================================================"
else
    echo "ERROR: Krona HTML was not created — check the error log"
    exit 1
fi

conda deactivate
