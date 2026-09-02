#!/bin/bash
#SBATCH --job-name=rgi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/rgi_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/rgi_%A_%a.err

# =============================================================================
# run_rgi.sh
# Step 7b-ii -- Per-sample assembly-based AMR detection with RGI 6.0.5 / CARD 4.0.1
#
# WHAT THIS SCRIPT DOES:
#   Screens Medaka-polished assemblies (consensus.fasta) against the CARD
#   database using RGI (Resistance Gene Identifier). Detects resistance genes
#   using BLAST, diamond, and the CARD protein homolog and variant models.
#
# RGI ALIGNMENT MODEL:
#   --input_type contig  -- assembly-based mode (not read-based)
#   --alignment_tool DIAMOND -- faster than BLAST for large metagenomes
#   --include_loose      -- include loose hits (lower confidence, more sensitive)
#
# INPUT PER SAMPLE:
#   polished/SampleName/consensus.fasta     -- Medaka-polished assembly
#
# OUTPUT PER SAMPLE:
#   amr_results/rgi/SampleName/rgi_output.txt   -- main hits table
#   amr_results/rgi/SampleName/rgi_output.json  -- full JSON results
#
# SAFE RESTART:
#   If rgi_output.txt already exists and is non-empty, the task skips.
#   Delete the sample folder to force re-run.
#
# DATABASE:
#   CARD 4.0.1 (loaded globally into envs/rgi with rgi load)
#
# RUN AFTER:  run_medaka.sh   (polished assemblies must exist)
# RUN BEFORE: downstream AMR summary / R analysis
#
# SUBMISSION:
#   sbatch scripts/run_rgi.sh
#   Monitor:  squeue -u $USER
#   Progress: ls amr_results/rgi/*/rgi_output.txt | wc -l
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
SAMPLE_LIST=$BASE/r_analysis/decontam/sample_list_real.txt
POLISHED=$BASE/polished
OUT_BASE=$BASE/amr_results/rgi
CONDA_ENV=$BASE/envs/rgi

# ── Validate conda environment ────────────────────────────────────────────────
if [ ! -d "$CONDA_ENV" ]; then
    echo "ERROR: RGI environment not found: $CONDA_ENV"
    echo "       Run: sbatch scripts/install_rgi.sh"
    exit 1
fi

conda activate "$CONDA_ENV"

# ── Validate RGI is available ─────────────────────────────────────────────────
if ! command -v rgi &> /dev/null; then
    echo "ERROR: rgi not found in environment: $CONDA_ENV"
    exit 1
fi

RGI_VERSION=$(rgi --version 2>&1 | head -1)
CARD_VERSION=$(rgi database --version 2>/dev/null)

# ── Validate CARD database is loaded ─────────────────────────────────────────
if [ -z "$CARD_VERSION" ]; then
    echo "ERROR: CARD database not loaded into RGI environment"
    echo "       Run: conda activate envs/rgi && rgi load --card_json databases/card/card.json"
    exit 1
fi

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
RESULT_FILE="$OUT_DIR/rgi_output.txt"

echo "======================================================"
echo "RGI assembly-based AMR detection"
echo "Array task:   $SLURM_ARRAY_TASK_ID / 61"
echo "Sample:       $SAMPLE"
echo "RGI:          $RGI_VERSION"
echo "CARD:         $CARD_VERSION"
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
    echo "RGI results already exist with $N hits: $RESULT_FILE"
    echo "Skipping $SAMPLE -- delete $OUT_DIR to re-run"
    exit 0
fi

mkdir -p "$OUT_DIR"

# ── Run RGI ───────────────────────────────────────────────────────────────────
echo ""
echo "Running RGI..."
echo "  Input:  $IN_FASTA"
echo "  Output: $OUT_DIR/rgi_output"
echo ""

rgi main \
    --input_sequence "$IN_FASTA" \
    --output_file "$OUT_DIR/rgi_output" \
    --input_type contig \
    --alignment_tool DIAMOND \
    --num_threads "$SLURM_CPUS_PER_TASK" \
    --include_loose \
    --clean

RGI_EXIT=$?

if [ $RGI_EXIT -ne 0 ]; then
    echo ""
    echo "ERROR: RGI failed for $SAMPLE (exit code $RGI_EXIT)"
    exit 1
fi

# ── Validate output ───────────────────────────────────────────────────────────
if [ ! -f "$RESULT_FILE" ]; then
    echo "ERROR: rgi_output.txt not found after run"
    exit 1
fi

# ── Final summary ─────────────────────────────────────────────────────────────
N_HITS=$(tail -n +2 "$RESULT_FILE" | grep -v "^$" | wc -l)
N_STRICT=$(tail -n +2 "$RESULT_FILE" | awk -F'\t' '$17=="Strict"' | wc -l)
N_LOOSE=$(tail -n +2 "$RESULT_FILE" | awk -F'\t' '$17=="Loose"' | wc -l)

echo ""
echo "======================================================"
echo "RGI complete: $SAMPLE"
echo "End time:  $(date)"
echo ""
echo "  Total hits:   $N_HITS"
echo "  Strict hits:  $N_STRICT"
echo "  Loose hits:   $N_LOOSE"
echo "  Results:      $RESULT_FILE"
echo ""
if [ "$N_HITS" -gt 0 ]; then
    echo "  Top strict hits:"
    tail -n +2 "$RESULT_FILE" | awk -F'\t' '$17=="Strict" {print "    " $9 " | " $16 " | Identity: " $10}' | head -10
fi
echo "======================================================"

conda deactivate
