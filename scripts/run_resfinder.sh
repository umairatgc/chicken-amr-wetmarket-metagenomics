#!/bin/bash
#SBATCH --job-name=resfinder
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/resfinder_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/resfinder_%A_%a.err

# =============================================================================
# run_resfinder.sh
# Step 7 -- Per-sample read-based AMR detection with ResFinder v4.7.2
#
# WHAT THIS SCRIPT DOES:
#   Screens decontaminated ONT reads (clean_reads/) directly against the CGE
#   ResFinder acquired resistance gene database using the KMA aligner.
#   KMA is designed for variable-length reads and is suitable for ONT long reads.
#   One job per sample, 61 jobs in parallel.
#
# WHY READ-BASED (NOT ASSEMBLY-BASED):
#   Read-based analysis detects AMR genes present at low coverage that may
#   be lost or fragmented during assembly. ResFinder + KMA is the recommended
#   approach for ONT metagenomics read-level AMR screening.
#
# PARAMETERS:
#   --acquired  : detect acquired resistance genes (horizontal gene transfer)
#   -l 0.60     : minimum gene coverage 60%
#   -t 0.90     : minimum nucleotide identity 90% (suitable for SUP R10.4.1 reads)
#
# NOTE ON POINTFINDER:
#   Chromosomal point mutation detection (--point) is NOT run here because
#   it requires species specification (--species). For isolate follow-up,
#   re-run with: --point --species "Escherichia coli" (or relevant species).
#
# INPUT PER SAMPLE:
#   clean_reads/SampleName_clean.fastq.gz   -- decontaminated ONT reads
#
# OUTPUT PER SAMPLE:
#   amr_results/resfinder/SampleName/ResFinder_results_tab.txt  -- main hits table
#   amr_results/resfinder/SampleName/ResFinder_results.txt      -- detailed results
#   amr_results/resfinder/SampleName/pheno_table.txt            -- phenotype predictions
#
# SAFE RESTART:
#   If ResFinder_results_tab.txt already exists and is non-empty, the task skips.
#   Delete the sample folder to force re-run.
#
# DATABASES USED:
#   ResFinder DB : databases/resfinder/resfinder_db/  (CGE, updated 2026-01-26)
#
# KMA INDEXING:
#   KMA requires binary index files (.length.b etc.) built from the .fsa files.
#   A fresh git clone does NOT include these. This script auto-builds them on
#   first run -- task 1 builds, all other tasks wait. No manual action needed.
#
# RUN AFTER:  run_extract_kraken.sh  (clean reads must exist)
# RUN BEFORE: downstream AMR summary / R analysis
#
# SUBMISSION:
#   sbatch run_resfinder.sh
#   Monitor:  squeue -u $USER
#   Progress: tail -f logs/resfinder_JOBID_TASKID.log
#   Count done: ls amr_results/resfinder/*/ResFinder_results_tab.txt 2>/dev/null | wc -l
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
CLEAN=$BASE/clean_reads
SAMPLE_LIST=$BASE/r_analysis/decontam/sample_list_real.txt
DB_RES=$BASE/databases/resfinder/resfinder_db
OUT_BASE=$BASE/amr_results/resfinder
CONDA_ENV=$BASE/envs/resfinder

# ── Validate conda environment exists ─────────────────────────────────────────
if [ ! -d "$CONDA_ENV" ]; then
    echo "ERROR: Conda environment not found: $CONDA_ENV"
    exit 1
fi

conda activate "$CONDA_ENV"

# ── Validate ResFinder is available ───────────────────────────────────────────
if ! python -m resfinder --help &> /dev/null; then
    echo "ERROR: ResFinder not found in conda environment: $CONDA_ENV"
    exit 1
fi

RESFINDER_VERSION=$(python -m resfinder --version 2>&1 | head -1)

# ── Validate databases exist ───────────────────────────────────────────────────
if [ ! -d "$DB_RES" ]; then
    echo "ERROR: ResFinder database not found: $DB_RES"
    exit 1
fi

# ── Build KMA database indices if missing ─────────────────────────────────────
# A fresh git clone of resfinder_db only contains .fsa text files.
# KMA requires pre-built binary index files (.length.b, .seq.b, .name etc.)
# Task 1 builds the indices; all other tasks wait until they are ready.
INDEX_CHECK="$DB_RES/aminoglycoside.length.b"

if [ ! -f "$INDEX_CHECK" ]; then
    if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
        echo "KMA index files not found -- building database indices (task 1 only)..."
        echo "This only runs once after a fresh database clone."
        for fsa in "$DB_RES"/*.fsa; do
            base=$(basename "$fsa" .fsa)
            echo "  Indexing ${base}.fsa..."
            kma index -i "$fsa" -o "$DB_RES/$base"
        done
        echo "KMA database indexing complete."
    else
        echo "Task $SLURM_ARRAY_TASK_ID waiting for task 1 to build KMA indices..."
        WAIT=0
        while [ ! -f "$INDEX_CHECK" ]; do
            sleep 15
            WAIT=$((WAIT + 15))
            echo "  Still waiting... ${WAIT}s elapsed"
            if [ "$WAIT" -gt 600 ]; then
                echo "ERROR: Timed out waiting for KMA indices after 600s"
                exit 1
            fi
        done
        echo "KMA indices ready -- proceeding."
    fi
else
    echo "KMA database indices found -- skipping indexing."
fi

# ── Validate sample list exists ───────────────────────────────────────────────
if [ ! -f "$SAMPLE_LIST" ]; then
    echo "ERROR: Sample list not found: $SAMPLE_LIST"
    exit 1
fi

mkdir -p "$OUT_BASE"

# ── Get this task's sample name ───────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

if [ -z "$SAMPLE" ]; then
    echo "ERROR: No sample found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# ── Input/output paths ────────────────────────────────────────────────────────
IN_FASTQ="$CLEAN/${SAMPLE}_clean.fastq.gz"
OUT_DIR="$OUT_BASE/$SAMPLE"
RESULT_FILE="$OUT_DIR/ResFinder_results_tab.txt"

echo "======================================================"
echo "ResFinder read-based AMR detection"
echo "Array task:          $SLURM_ARRAY_TASK_ID / 61"
echo "Sample:              $SAMPLE"
echo "ResFinder version:   $RESFINDER_VERSION"
echo "Identity threshold:  90%"
echo "Coverage threshold:  60%"
echo "CPUs:                $SLURM_CPUS_PER_TASK"
echo "Start time:          $(date)"
echo "======================================================"

# ── Validate input FASTQ ──────────────────────────────────────────────────────
if [ ! -f "$IN_FASTQ" ]; then
    echo "ERROR: Clean FASTQ not found: $IN_FASTQ"
    echo "       Run run_extract_kraken.sh first"
    exit 1
fi

if [ ! -s "$IN_FASTQ" ]; then
    echo "ERROR: Clean FASTQ is empty: $IN_FASTQ"
    exit 1
fi

echo "Input validated: $IN_FASTQ"
echo "File size:       $(du -h "$IN_FASTQ" | cut -f1)"

# ── Skip if already complete ───────────────────────────────────────────────────
if [ -f "$RESULT_FILE" ] && [ -s "$RESULT_FILE" ]; then
    N=$(tail -n +2 "$RESULT_FILE" | wc -l)
    echo ""
    echo "ResFinder results already exist with $N hits: $RESULT_FILE"
    echo "Skipping $SAMPLE -- delete $OUT_DIR to re-run"
    exit 0
fi

mkdir -p "$OUT_DIR"

# ── Run ResFinder ──────────────────────────────────────────────────────────────
echo ""
echo "Running ResFinder..."
echo "  Input:    $IN_FASTQ"
echo "  Database: $DB_RES"
echo "  Output:   $OUT_DIR"
echo ""

python -m resfinder \
    -ifq "$IN_FASTQ" \
    --acquired \
    -db_res "$DB_RES" \
    -l 0.60 \
    -t 0.90 \
    -o "$OUT_DIR"

RESFINDER_EXIT=$?

if [ $RESFINDER_EXIT -ne 0 ]; then
    echo ""
    echo "ERROR: ResFinder failed for $SAMPLE (exit code $RESFINDER_EXIT)"
    exit 1
fi

# ── Validate output ───────────────────────────────────────────────────────────
if [ ! -f "$RESULT_FILE" ]; then
    echo "ERROR: ResFinder_results_tab.txt not found after run"
    echo "       ResFinder may have found no hits or encountered an error"
    exit 1
fi

# ── Final summary ─────────────────────────────────────────────────────────────
N_HITS=$(tail -n +2 "$RESULT_FILE" | grep -v "^$" | wc -l)

echo ""
echo "======================================================"
echo "ResFinder complete: $SAMPLE"
echo "End time:  $(date)"
echo ""
echo "  AMR genes detected: $N_HITS"
echo "  Results:  $RESULT_FILE"
echo ""
if [ "$N_HITS" -gt 0 ]; then
    echo "  Top hits:"
    tail -n +2 "$RESULT_FILE" | awk -F'\t' '{print "    " $1 " | " $2 " | Identity: " $4}' | head -10
fi
echo "======================================================"

conda deactivate
