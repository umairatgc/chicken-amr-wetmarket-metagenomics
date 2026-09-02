#!/bin/bash
#SBATCH --job-name=amrfinder
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/amrfinder_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/amrfinder_%A_%a.err

# =============================================================================
# run_amrfinder.sh
# Step 7b-i -- Per-sample assembly-based AMR detection with AMRFinderPlus 4.2.7
#
# WHAT THIS SCRIPT DOES:
#   Screens Medaka-polished assemblies (consensus.fasta) against the NCBI
#   AMRFinderPlus database (NDARO 2026-05-15.1). Detects acquired resistance
#   genes and point mutations in assembled contigs.
#
# WHY ASSEMBLY-BASED (NOT READ-BASED):
#   AMRFinderPlus requires assembled sequences as input -- it has no read
#   mapping mode. Using polished assemblies gives complete gene context
#   including flanking sequence and mobile element associations.
#
# NO --organism FLAG:
#   Samples are metagenomes with mixed bacterial communities. Specifying a
#   single organism would be incorrect. Running without --organism detects
#   all acquired resistance genes across all species present.
#
# INPUT PER SAMPLE:
#   polished/SampleName/consensus.fasta     -- Medaka-polished assembly
#
# OUTPUT PER SAMPLE:
#   amr_results/amrfinder/SampleName/amrfinder_results.txt  -- main hits table
#
# SAFE RESTART:
#   If amrfinder_results.txt already exists and is non-empty, the task skips.
#   Delete the sample folder to force re-run.
#
# DATABASE:
#   NCBI NDARO 2026-05-15.1
#   Location: envs/amrfinder/share/amrfinderplus/data/
#
# RUN AFTER:  run_medaka.sh   (polished assemblies must exist)
# RUN BEFORE: downstream AMR summary / R analysis
#
# SUBMISSION:
#   sbatch scripts/run_amrfinder.sh
#   Monitor:  squeue -u $USER
#   Progress: ls amr_results/amrfinder/*/amrfinder_results.txt | wc -l
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
SAMPLE_LIST=$BASE/r_analysis/decontam/sample_list_real.txt
POLISHED=$BASE/polished
OUT_BASE=$BASE/amr_results/amrfinder
CONDA_ENV=$BASE/envs/amrfinder

# ── Validate conda environment ────────────────────────────────────────────────
if [ ! -d "$CONDA_ENV" ]; then
    echo "ERROR: AMRFinderPlus environment not found: $CONDA_ENV"
    echo "       Run: sbatch scripts/install_amrfinder.sh"
    exit 1
fi

conda activate "$CONDA_ENV"

# ── Validate AMRFinderPlus is available ───────────────────────────────────────
if ! command -v amrfinder &> /dev/null; then
    echo "ERROR: amrfinder not found in environment: $CONDA_ENV"
    exit 1
fi

AMRFINDER_VERSION=$(amrfinder --version 2>/dev/null)

# ── Validate database exists ──────────────────────────────────────────────────
DB_DIR=$(ls -d "$CONDA_ENV/share/amrfinderplus/data/"*/ 2>/dev/null | tail -1)
if [ -z "$DB_DIR" ]; then
    echo "ERROR: AMRFinderPlus database not found"
    echo "       Run on login node: amrfinder --update"
    exit 1
fi
DB_VERSION=$(basename "$DB_DIR")

# ── Validate sample list ──────────────────────────────────────────────────────
if [ ! -f "$SAMPLE_LIST" ]; then
    echo "ERROR: Sample list not found: $SAMPLE_LIST"
    exit 1
fi

mkdir -p "$OUT_BASE"

# ── Get this task's sample ────────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

if [ -z "$SAMPLE" ]; then
    echo "ERROR: No sample found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# ── Input / output paths ──────────────────────────────────────────────────────
IN_FASTA="$POLISHED/$SAMPLE/consensus.fasta"
OUT_DIR="$OUT_BASE/$SAMPLE"
RESULT_FILE="$OUT_DIR/amrfinder_results.txt"

echo "======================================================"
echo "AMRFinderPlus assembly-based AMR detection"
echo "Array task:   $SLURM_ARRAY_TASK_ID / 61"
echo "Sample:       $SAMPLE"
echo "Version:      $AMRFINDER_VERSION"
echo "Database:     $DB_VERSION"
echo "CPUs:         $SLURM_CPUS_PER_TASK"
echo "Start time:   $(date)"
echo "======================================================"

# ── Validate input ────────────────────────────────────────────────────────────
if [ ! -f "$IN_FASTA" ]; then
    echo "ERROR: Polished assembly not found: $IN_FASTA"
    echo "       Run run_medaka.sh first"
    exit 1
fi

if [ ! -s "$IN_FASTA" ]; then
    echo "ERROR: Assembly file is empty: $IN_FASTA"
    exit 1
fi

N_CONTIGS=$(grep -c "^>" "$IN_FASTA" || true)
echo "Input contigs: $N_CONTIGS"

# ── Skip if already complete ──────────────────────────────────────────────────
if [ -f "$RESULT_FILE" ] && [ -s "$RESULT_FILE" ]; then
    N=$(tail -n +2 "$RESULT_FILE" | wc -l)
    echo ""
    echo "AMRFinder results already exist with $N hits: $RESULT_FILE"
    echo "Skipping $SAMPLE -- delete $OUT_DIR to re-run"
    exit 0
fi

mkdir -p "$OUT_DIR"

# ── Run AMRFinderPlus ─────────────────────────────────────────────────────────
echo ""
echo "Running AMRFinderPlus..."
echo "  Input:  $IN_FASTA"
echo "  Output: $RESULT_FILE"
echo ""

amrfinder \
    --nucleotide "$IN_FASTA" \
    --output "$RESULT_FILE" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --plus

AMRFINDER_EXIT=$?

if [ $AMRFINDER_EXIT -ne 0 ]; then
    echo ""
    echo "ERROR: AMRFinderPlus failed for $SAMPLE (exit code $AMRFINDER_EXIT)"
    exit 1
fi

# ── Validate output ───────────────────────────────────────────────────────────
if [ ! -f "$RESULT_FILE" ]; then
    echo "ERROR: amrfinder_results.txt not found after run"
    exit 1
fi

# ── Final summary ─────────────────────────────────────────────────────────────
N_HITS=$(tail -n +2 "$RESULT_FILE" | grep -v "^$" | wc -l)

echo ""
echo "======================================================"
echo "AMRFinderPlus complete: $SAMPLE"
echo "End time:  $(date)"
echo ""
echo "  AMR genes detected: $N_HITS"
echo "  Results:  $RESULT_FILE"
echo ""
if [ "$N_HITS" -gt 0 ]; then
    echo "  Top hits:"
    tail -n +2 "$RESULT_FILE" | awk -F'\t' '{print "    " $6 " | " $9 " | Identity: " $12}' | head -10
fi
echo "======================================================"

conda deactivate
