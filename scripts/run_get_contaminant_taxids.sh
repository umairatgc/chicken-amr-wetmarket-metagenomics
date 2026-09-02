#!/bin/bash
#SBATCH --job-name=get_taxids
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/get_taxids_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/get_taxids_%j.err

# =============================================================================
# run_get_contaminant_taxids.sh
# Step 5.1 — Look up NCBI taxids for contaminant species from decontam output
#
# WHAT THIS SCRIPT DOES:
#   Reads contaminant species names from a plain text file (one name per line)
#   and searches the Kraken2 report files to find the NCBI taxonomy ID for
#   each species. Saves the taxids to a file used by run_extract_kraken.sh.
#
# INPUTS (set the three variables below):
#   SPECIES_FILE  — plain text file, one species name per line
#                   (the contaminant_species_list.txt from decontam)
#   REPORT_DIR    — folder containing *_kraken2.report files
#   OUT_DIR       — where to save the taxid output files
#
# HOW TAXID LOOKUP WORKS:
#   Kraken2 report format (tab-separated):
#   [%]  [reads_clade]  [reads_direct]  [rank]  [taxid]  [name]
#   This script greps each species name across all reports and extracts
#   column 5 (the taxid). Uses the first match found across all reports.
#
# OUTPUT:
#   contaminant_taxids.txt        — one taxid per line (used by run_extract_kraken.sh)
#   contaminant_taxids_named.txt  — taxid + species name (human readable check)
#
# Run AFTER:  run_decontam.sh  (needs contaminant_species_list.txt)
# Run BEFORE: run_extract_kraken.sh
# =============================================================================

# ── Configure these three paths ───────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1

SPECIES_FILE=$BASE/r_analysis/decontam/contaminant_species_list.txt
REPORT_DIR=$BASE/kraken2_output
OUT_DIR=$BASE/r_analysis/decontam

# ── Derived output paths ──────────────────────────────────────────────────────
# Saved alongside contaminant_species_list.txt — all decontam outputs in one place
OUT_IDS=$OUT_DIR/contaminant_taxids.txt
OUT_NAMED=$OUT_DIR/contaminant_taxids_named.txt

mkdir -p "$OUT_DIR"

echo "======================================================"
echo "get_contaminant_taxids"
echo "Start time: $(date)"
echo "======================================================"
echo ""
echo "Species file:  $SPECIES_FILE"
echo "Report dir:    $REPORT_DIR"
echo "Output dir:    $OUT_DIR"
echo ""

# ── Validate inputs ───────────────────────────────────────────────────────────
if [ ! -f "$SPECIES_FILE" ]; then
    echo "ERROR: Species file not found: $SPECIES_FILE"
    echo "       Run run_decontam.sh first to generate contaminant_species_list.txt"
    exit 1
fi

N_REPORTS=$(ls "$REPORT_DIR"/*_kraken2.report 2>/dev/null | wc -l)
if [ "$N_REPORTS" -eq 0 ]; then
    echo "ERROR: No Kraken2 report files found in $REPORT_DIR"
    echo "       Expected files named *_kraken2.report"
    exit 1
fi

N_SPECIES=$(grep -c . "$SPECIES_FILE")
echo "Species to look up: $N_SPECIES"
echo "Kraken2 reports:    $N_REPORTS"
echo ""

# ── Look up taxid for each species ────────────────────────────────────────────
> "$OUT_IDS"
> "$OUT_NAMED"

FOUND=0
NOT_FOUND=()

while IFS= read -r species || [ -n "$species" ]; do

    # Skip empty lines
    [ -z "$species" ] && continue

    # Grep all reports for this species name, extract taxid (column 5)
    # -h suppresses filename prefix
    # awk filters to lines where col 5 is a number (taxid)
    # sort -u removes duplicates, head -1 takes first result
    taxid=$(grep -rh "$species" "$REPORT_DIR"/*_kraken2.report 2>/dev/null \
            | awk 'NF>=5 && $5~/^[0-9]+$/ {print $5}' \
            | sort -u \
            | head -1)

    if [ -n "$taxid" ]; then
        echo "$taxid" >> "$OUT_IDS"
        echo "$taxid    $species" >> "$OUT_NAMED"
        printf "  FOUND   taxid=%-10s %s\n" "$taxid" "$species"
        FOUND=$((FOUND + 1))
    else
        NOT_FOUND+=("$species")
        printf "  MISSING %-10s %s\n" "" "$species"
    fi

done < "$SPECIES_FILE"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "======================================================"
echo "Taxid lookup complete"
echo "  Found:   $FOUND / $N_SPECIES species"
echo "  Missing: ${#NOT_FOUND[@]} species"
echo ""

if [ ${#NOT_FOUND[@]} -gt 0 ]; then
    echo "  Species not found in any Kraken2 report:"
    for sp in "${NOT_FOUND[@]}"; do
        echo "    - $sp"
    done
    echo ""
    echo "  NOTE: Missing species were never detected in any of your samples."
    echo "  This is fine — if Kraken2 never saw them they cannot be in your reads."
    echo "  They are simply not added to the taxid list."
    echo ""
fi

echo "  Output files:"
echo "    $OUT_IDS"
echo "    $OUT_NAMED"
echo ""
echo "  Full taxid list:"
cat "$OUT_NAMED"
echo ""
echo "Next step:"
echo "  sbatch run_extract_kraken.sh"
echo "End time: $(date)"
echo "======================================================"
