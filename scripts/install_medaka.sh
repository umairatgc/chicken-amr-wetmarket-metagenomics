#!/bin/bash
#SBATCH --job-name=install_medaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --clusters=htc
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/install_medaka_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/install_medaka_%j.err

# =============================================================================
# install_medaka.sh
# Install Medaka environment for GPU-accelerated assembly polishing
#
# ENVIRONMENT CREATED:
#   envs/medaka  — Medaka (latest stable) + minimap2  (Python 3.10)
#
# WHAT MEDAKA DOES:
#   Polishes Flye assemblies using the original reads and a neural network
#   model trained on ONT data. Corrects remaining errors in contigs,
#   especially homopolymer regions. Runs on GPU for speed.
#
# MODEL SELECTION:
#   Your chemistry: dna_r10.4.1_e8.2_400bps_sup@v5.2.0
#   Medaka model:   r1041_e82_400bps_sup_v5.0.0
#   This script verifies the model is available after installation.
#
# CLUSTER:
#   Medaka requires a GPU — runs on the HTC cluster, not ARC.
#   Submitted with --clusters=htc and --gres=gpu:1
#
# NO DATABASE NEEDED:
#   Medaka uses neural network models bundled with the installation.
#   No separate database download required.
#
# RUN AFTER:  run_assembly.sh   (assemblies must exist)
# RUN BEFORE: run_medaka.sh     (polishing script)
#
# SUBMISSION:
#   sbatch scripts/install_medaka.sh
#   Monitor: tail -f logs/install_medaka_JOBID.log
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"

BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
MEDAKA_ENV=$BASE/envs/medaka

mkdir -p "$BASE/logs"

echo "=============================================="
echo "Medaka Environment Installation"
echo "Start time:  $(date)"
echo "Node:        $(hostname)"
echo "GPU:         $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'GPU info unavailable')"
echo "=============================================="

# ── Remove existing medaka env if present ─────────────────────────
if [ -d "$MEDAKA_ENV" ]; then
    echo ""
    echo "Removing existing medaka environment..."
    conda env remove --prefix "$MEDAKA_ENV" -y
    echo "Old medaka environment removed."
fi

# ── Create environment ────────────────────────────────────────────
echo ""
echo "Creating medaka environment (Python 3.10)..."
echo "  Medaka (latest stable)"
echo "  minimap2 (bundled dependency)"
echo ""

conda create --prefix "$MEDAKA_ENV" \
    -c conda-forge -c bioconda \
    python=3.10 \
    medaka \
    -y

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create medaka environment"
    exit 1
fi

conda activate "$MEDAKA_ENV"

# ── Verify Medaka ─────────────────────────────────────────────────
echo ""
echo "Verifying medaka environment..."

MEDAKA_VERSION=$(medaka --version 2>&1 | head -1)
if [ -z "$MEDAKA_VERSION" ]; then
    echo "ERROR: Medaka not working in environment"
    exit 1
fi
echo "Medaka:   OK  ($MEDAKA_VERSION)"

MINIMAP_VERSION=$(minimap2 --version 2>&1)
echo "minimap2: OK  ($MINIMAP_VERSION)"

# ── Check GPU is visible ──────────────────────────────────────────
echo ""
echo "Checking GPU availability..."
if nvidia-smi &> /dev/null; then
    nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
    echo "GPU: accessible"
else
    echo "WARNING: nvidia-smi not available -- GPU may not be accessible"
fi

# ── Verify correct Medaka model is available ─────────────────────
# Model must match basecall chemistry:
# dna_r10.4.1_e8.2_400bps_sup@v5.2.0 → r1041_e82_400bps_sup_v5.0.0
echo ""
echo "Checking Medaka models for R10.4.1 SUP chemistry..."
medaka tools list_models 2>/dev/null | grep r1041 || echo "No r1041 models found"

TARGET_MODEL="r1041_e82_400bps_sup_v5.0.0"
if medaka tools list_models 2>/dev/null | grep -q "$TARGET_MODEL"; then
    echo ""
    echo "Target model confirmed available: $TARGET_MODEL"
else
    echo ""
    echo "WARNING: Target model $TARGET_MODEL not found"
    echo "         Available r1041 models listed above"
    echo "         Update MODEL variable in run_medaka.sh accordingly"
fi

conda deactivate

# ── Final summary ─────────────────────────────────────────────────
echo ""
echo "=============================================="
echo "Installation complete: $(date)"
echo ""
echo "  envs/medaka — Medaka + minimap2  ✓"
echo "  No database download required    ✓"
echo ""
echo "Next step: sbatch scripts/run_medaka.sh"
echo "=============================================="
