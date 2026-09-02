#!/bin/bash
#SBATCH --job-name=extract_kraken
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/extract_kraken_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/extract_kraken_%A_%a.err

# =============================================================================
# run_extract_kraken.sh
# Step 5.2 — Remove contaminant reads from filtered FASTQ files
#
# WHAT THIS SCRIPT DOES:
#   For each of the 61 real samples (negative controls excluded), uses
#   extract_kraken_reads.py to remove reads that Kraken2 classified as one
#   of the 39 decontam-identified contaminant species.
#   Outputs a clean FASTQ.gz file per sample ready for:
#     - Assembly (Flye/metaFlye)
#     - Contig-based AMR detection (AMRFinderPlus)
#     - Read-based AMR detection (CARD/RGI)
#
# HOW IT WORKS:
#   Kraken2 classified every read in your filtered FASTQ files and recorded
#   the taxonomy ID (taxid) of each classified read in the _kraken2.output
#   files. extract_kraken_reads.py reads those classifications and uses
#   --exclude to keep every read EXCEPT those classified to the contaminant
#   taxids. Unclassified reads (taxid=0) are always kept.
#
# INPUT PER SAMPLE:
#   filtered/SampleName_filtered.fastq.gz    — chopper-filtered ONT reads
#   kraken2_output/SampleName_kraken2.output — per-read Kraken2 classifications
#   scripts/contaminant_taxids.txt           — taxids from get_contaminant_taxids.sh
#
# OUTPUT PER SAMPLE:
#   clean_reads/SampleName_clean.fastq.gz    — contaminant-free reads
#
# RUN AFTER:  bash get_contaminant_taxids.sh
# RUN BEFORE: run_assembly.sh  and  run_amr.sh
#
# SUBMISSION:
#   sbatch run_extract_kraken.sh
#   (runs as 61-task array — one task per real sample, all in parallel)
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/taxonomy

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
FILTERED=$BASE/filtered
KRAKEN_OUT=$BASE/kraken2_output
CLEAN=$BASE/clean_reads
TAXID_FILE=$BASE/r_analysis/decontam/contaminant_taxids.txt
SAMPLE_LIST=$BASE/r_analysis/decontam/sample_list_real.txt

mkdir -p "$CLEAN"

# ── Build sample list once (task 1 only, if file doesn't exist) ───────────────
# Excludes negative controls AE-buffer and EB-neg
if [ ! -f "$SAMPLE_LIST" ]; then
    ls "$FILTERED"/*_filtered.fastq.gz \
        | xargs -n1 basename \
        | sed 's/_filtered\.fastq\.gz//' \
        | grep -v "^AE-buffer$\|^EB-neg$" \
        | sort > "$SAMPLE_LIST"
    echo "Sample list created: $(wc -l < $SAMPLE_LIST) real samples"
fi

# ── Get this task's sample name ───────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

if [ -z "$SAMPLE" ]; then
    echo "ERROR: No sample found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# ── Input/output paths for this sample ───────────────────────────────────────
IN_FASTQ="$FILTERED/${SAMPLE}_filtered.fastq.gz"
IN_KRAKEN="$KRAKEN_OUT/${SAMPLE}_kraken2.output"
IN_REPORT="$KRAKEN_OUT/${SAMPLE}_kraken2.report"
OUT_FASTQ="$CLEAN/${SAMPLE}_clean.fastq"
OUT_GZ="$CLEAN/${SAMPLE}_clean.fastq.gz"

echo "======================================================"
echo "extract_kraken_reads — Remove Contaminant Reads"
echo "Array task:  $SLURM_ARRAY_TASK_ID / 61"
echo "Sample:      $SAMPLE"
echo "Start time:  $(date)"
echo "======================================================"

# ── Validate inputs ───────────────────────────────────────────────────────────
if [ ! -f "$IN_FASTQ" ]; then
    echo "ERROR: Filtered FASTQ not found: $IN_FASTQ"
    exit 1
fi

if [ ! -f "$IN_KRAKEN" ]; then
    echo "ERROR: Kraken2 output not found: $IN_KRAKEN"
    exit 1
fi

if [ ! -f "$IN_REPORT" ]; then
    echo "ERROR: Kraken2 report not found: $IN_REPORT"
    exit 1
fi

if [ ! -f "$TAXID_FILE" ]; then
    echo "ERROR: Taxid file not found: $TAXID_FILE"
    echo "       Run bash get_contaminant_taxids.sh first"
    exit 1
fi

# ── Load taxids into a space-separated list ───────────────────────────────────
TAXIDS=$(cat "$TAXID_FILE" | tr '\n' ' ' | sed 's/ $//')
N_TAXIDS=$(wc -l < "$TAXID_FILE")

echo ""
echo "Contaminant taxids loaded: $N_TAXIDS"
echo "Taxids: $TAXIDS"
echo ""

# ── Count input reads ─────────────────────────────────────────────────────────
IN_READS=$(zcat "$IN_FASTQ" | awk 'NR%4==1' | wc -l)
echo "Input reads:  $IN_READS"

# ── Run extract_kraken_reads.py ───────────────────────────────────────────────
echo ""
echo "Running extract_kraken_reads.py (--exclude mode)..."
echo "  Input FASTQ:   $IN_FASTQ"
echo "  Kraken output: $IN_KRAKEN"
echo "  Kraken report: $IN_REPORT"
echo "  Output FASTQ:  $OUT_FASTQ"
echo ""

extract_kraken_reads.py \
    -k "$IN_KRAKEN" \
    -s "$IN_FASTQ" \
    -o "$OUT_FASTQ" \
    --report "$IN_REPORT" \
    --exclude \
    --include-children \
    --fastq-output \
    --taxid $TAXIDS

if [ $? -ne 0 ]; then
    echo "ERROR: extract_kraken_reads.py failed for $SAMPLE"
    exit 1
fi

# ── Compress output ───────────────────────────────────────────────
echo "Compressing output..."
gzip -c "$OUT_FASTQ" > "${OUT_GZ}.tmp"

if [ $? -ne 0 ]; then
    echo "ERROR: gzip compression failed for $SAMPLE"
    rm -f "${OUT_GZ}.tmp"
    exit 1
fi

# Verify the gzip file is complete and not truncated
echo "Verifying gzip integrity..."
gzip -t "${OUT_GZ}.tmp"
if [ $? -ne 0 ]; then
    echo "ERROR: gzip output is corrupted/truncated for $SAMPLE"
    rm -f "${OUT_GZ}.tmp"
    exit 1
fi

# Only move to final name after verified complete
mv "${OUT_GZ}.tmp" "$OUT_GZ"
rm -f "$OUT_FASTQ"

if [ ! -f "$OUT_GZ" ]; then
    echo "ERROR: Final gz file missing after move — $OUT_GZ"
    exit 1
fi

# ── Count output reads and report ─────────────────────────────────────────────
OUT_READS=$(zcat "$OUT_GZ" | awk 'NR%4==1' | wc -l)
REMOVED=$((IN_READS - OUT_READS))
PCT_REMOVED=$(awk "BEGIN {printf \"%.2f\", ($REMOVED/$IN_READS)*100}")
PCT_KEPT=$(awk "BEGIN {printf \"%.2f\", ($OUT_READS/$IN_READS)*100}")

echo ""
echo "======================================================"
echo "extract_kraken_reads complete for: $SAMPLE"
echo "End time: $(date)"
echo ""
echo "  Input reads:    $IN_READS"
echo "  Reads removed:  $REMOVED  ($PCT_REMOVED% — contaminant reads)"
echo "  Reads kept:     $OUT_READS  ($PCT_KEPT%)"
echo ""
echo "  Output: $OUT_GZ"
echo "======================================================"

conda deactivate
