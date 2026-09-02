#!/bin/bash
#SBATCH --job-name=medaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-61
#SBATCH --clusters=htc
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/medaka_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/medaka_%A_%a.err

# =============================================================================
# run_medaka.sh
# Step 6b -- Per-sample assembly polishing with Medaka (GPU-accelerated)
#
# WHAT THIS SCRIPT DOES:
#   Polishes Flye metagenomic assemblies using the original ONT reads and a
#   neural network model trained on R10.4.1 SUP data. Corrects homopolymer
#   errors and other systematic errors remaining after Flye assembly.
#   Runs mini_align (CPU) then medaka consensus (GPU) then medaka stitch (CPU).
#
# MODEL:
#   Chemistry:    dna_r10.4.1_e8.2_400bps_sup@v5.2.0
#   Medaka model: r1041_e82_400bps_sup_v5.0.0
#
# CUDNN FIX (built-in):
#   Some HTC GPU nodes reject non-contiguous tensors in the GRU forward pass,
#   causing: RuntimeError: cuDNN error: CUDNN_STATUS_NOT_SUPPORTED
#   This script auto-patches medaka's gru.py to call .contiguous() before GRU.
#   The patch is idempotent -- safe to apply multiple times.
#
# INPUT PER SAMPLE:
#   assemblies/SampleName/assembly.fasta     -- Flye assembly
#   clean_reads/SampleName_clean.fastq.gz    -- decontaminated ONT reads
#
# OUTPUT PER SAMPLE:
#   polished/SampleName/consensus.fasta      -- Medaka-polished assembly
#
# SAFE RESTART:
#   If consensus.fasta already exists and is non-empty, the task skips.
#   Delete the sample folder under polished/ to force re-run.
#
# CLUSTER:
#   Medaka requires a GPU -- must run on HTC cluster (not ARC).
#   Submit from htc-login: sbatch scripts/run_medaka.sh
#
# RUN AFTER:  run_assembly.sh   (Flye assemblies must exist)
# RUN BEFORE: run_amrfinder.sh, run_rgi.sh, run_abricate.sh
#
# SUBMISSION:
#   sbatch scripts/run_medaka.sh
#   Monitor:    squeue --clusters=htc -u $USER
#   Progress:   tail -f logs/medaka_JOBID_TASKID.log
#   Count done: ls polished/*/consensus.fasta 2>/dev/null | wc -l
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
SAMPLE_LIST=$BASE/r_analysis/decontam/sample_list_real.txt
ASSEMBLY_DIR=$BASE/assemblies
CLEAN=$BASE/clean_reads
OUT_BASE=$BASE/polished
CONDA_ENV=$BASE/envs/medaka
MODEL=r1041_e82_400bps_sup_v5.0.0

# ── Validate conda environment ────────────────────────────────────────────────
if [ ! -d "$CONDA_ENV" ]; then
    echo "ERROR: Medaka environment not found: $CONDA_ENV"
    echo "       Run: sbatch scripts/install_medaka.sh"
    exit 1
fi

conda activate "$CONDA_ENV"

# ── Validate medaka is available ──────────────────────────────────────────────
if ! command -v medaka &> /dev/null; then
    echo "ERROR: medaka not found in environment: $CONDA_ENV"
    exit 1
fi

MEDAKA_VERSION=$(medaka --version 2>&1 | head -1)

# ── Apply cuDNN GRU patch ─────────────────────────────────────────────────────
# Some HTC GPU nodes (older architectures) fail with:
#   RuntimeError: cuDNN error: CUDNN_STATUS_NOT_SUPPORTED
#   This error may appear if you passed in a non-contiguous input.
#
# Root cause: medaka's GRU layer receives a non-contiguous tensor.
# cuDNN on certain GPUs requires contiguous memory layout.
#
# Fix: patch gru.py to call .contiguous() before passing to GRU.
# This patch is idempotent -- if already applied, sed makes no change.
echo ""
echo "Applying cuDNN GRU patch..."
GRU_FILE=$(python -c "import os, medaka; print(os.path.join(os.path.dirname(medaka.__file__), 'architectures', 'gru.py'))")

if [ ! -f "$GRU_FILE" ]; then
    echo "ERROR: gru.py not found at: $GRU_FILE"
    exit 1
fi

# Check if patch already applied
if grep -q "x.contiguous()" "$GRU_FILE"; then
    echo "cuDNN patch already applied -- skipping."
else
    sed -i 's/x = self\.gru(x)\[0\]/x = self.gru(x.contiguous())[0]/' "$GRU_FILE"
    if grep -q "x.contiguous()" "$GRU_FILE"; then
        echo "cuDNN patch applied successfully."
    else
        echo "ERROR: cuDNN patch failed to apply to $GRU_FILE"
        exit 1
    fi
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
IN_ASSEMBLY="$ASSEMBLY_DIR/$SAMPLE/assembly.fasta"
IN_READS="$CLEAN/${SAMPLE}_clean.fastq.gz"
OUT_DIR="$OUT_BASE/$SAMPLE"
RESULT_FILE="$OUT_DIR/consensus.fasta"

echo "======================================================"
echo "Medaka assembly polishing"
echo "Array task:   $SLURM_ARRAY_TASK_ID / 61"
echo "Sample:       $SAMPLE"
echo "Medaka:       $MEDAKA_VERSION"
echo "Model:        $MODEL"
echo "CPUs:         $SLURM_CPUS_PER_TASK"
echo "GPU:          $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'unavailable')"
echo "Start time:   $(date)"
echo "======================================================"

# ── Validate inputs ───────────────────────────────────────────────────────────
if [ ! -f "$IN_ASSEMBLY" ]; then
    echo "ERROR: Flye assembly not found: $IN_ASSEMBLY"
    echo "       Run run_assembly.sh first"
    exit 1
fi

if [ ! -s "$IN_ASSEMBLY" ]; then
    echo "ERROR: Assembly file is empty: $IN_ASSEMBLY"
    exit 1
fi

if [ ! -f "$IN_READS" ]; then
    echo "ERROR: Clean reads not found: $IN_READS"
    echo "       Run run_extract_kraken.sh first"
    exit 1
fi

if [ ! -s "$IN_READS" ]; then
    echo "ERROR: Clean reads file is empty: $IN_READS"
    exit 1
fi

ASSEMBLY_SIZE=$(grep -c "^>" "$IN_ASSEMBLY" || true)
echo "Assembly contigs: $ASSEMBLY_SIZE"
echo "Reads:            $(du -h "$IN_READS" | cut -f1)"

# ── Skip if already complete ──────────────────────────────────────────────────
if [ -f "$RESULT_FILE" ] && [ -s "$RESULT_FILE" ]; then
    N=$(grep -c "^>" "$RESULT_FILE" || true)
    echo ""
    echo "Polished assembly already exists with $N contigs: $RESULT_FILE"
    echo "Skipping $SAMPLE -- delete $OUT_DIR to re-run"
    exit 0
fi

mkdir -p "$OUT_DIR"

# ── Run Medaka ────────────────────────────────────────────────────────────────
# medaka_consensus wrapper runs:
#   1. mini_align       -- maps reads to assembly (CPU)
#   2. medaka consensus -- neural network polishing (GPU)
#   3. medaka stitch    -- stitches polished segments back together (CPU)
echo ""
echo "Running Medaka polishing..."
echo "  Assembly: $IN_ASSEMBLY"
echo "  Reads:    $IN_READS"
echo "  Model:    $MODEL"
echo "  Output:   $OUT_DIR"
echo ""

medaka_consensus \
    -i "$IN_READS" \
    -d "$IN_ASSEMBLY" \
    -o "$OUT_DIR" \
    -m "$MODEL" \
    -t "$SLURM_CPUS_PER_TASK"

MEDAKA_EXIT=$?

if [ $MEDAKA_EXIT -ne 0 ]; then
    echo ""
    echo "ERROR: Medaka failed for $SAMPLE (exit code $MEDAKA_EXIT)"
    exit 1
fi

# ── Validate output ───────────────────────────────────────────────────────────
if [ ! -f "$RESULT_FILE" ]; then
    echo "ERROR: consensus.fasta not found after Medaka run"
    exit 1
fi

if [ ! -s "$RESULT_FILE" ]; then
    echo "ERROR: consensus.fasta is empty"
    exit 1
fi

# ── Final summary ─────────────────────────────────────────────────────────────
N_CONTIGS=$(grep -c "^>" "$RESULT_FILE" || true)
FILESIZE=$(du -h "$RESULT_FILE" | cut -f1)

echo ""
echo "======================================================"
echo "Medaka complete: $SAMPLE"
echo "End time:  $(date)"
echo ""
echo "  Input contigs:   $ASSEMBLY_SIZE"
echo "  Output contigs:  $N_CONTIGS"
echo "  Output size:     $FILESIZE"
echo "  Result:          $RESULT_FILE"
echo "======================================================"

conda deactivate
