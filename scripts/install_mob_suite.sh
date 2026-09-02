#!/bin/bash
#SBATCH --job-name=install_mob_suite
#SBATCH --clusters=arc
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/install_mob_suite.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/install_mob_suite.err

# =============================================================================
# install_mob_suite.sh
# Creates a new conda environment and installs MOBsuite
#
# USAGE:
#   sbatch scripts/install_mobsuite.sh
#
# When complete, verify with:
#   conda run -p "$ENV_DIR" mob_recon --version
# =============================================================================

ENV_DIR="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/mob_suite"

echo "========================================"
echo "  Installing MOBsuite"
echo "  Environment: $ENV_DIR"
echo "  Started: $(date)"
echo "========================================"

# ── Load Anaconda ─────────────────────────────────────────────────────────────
module purge
module load Anaconda3

# ── Remove existing environment if broken ─────────────────────────────────────
if [ -d "$ENV_DIR" ]; then
    echo "Removing existing environment..."
    conda env remove -p "$ENV_DIR" -y
fi

# ── Create environment and install MOBsuite ───────────────────────────────────
echo ""
echo "Creating conda environment and installing MOBsuite..."
conda create -y \
    -p "$ENV_DIR" \
    -c conda-forge \
    -c bioconda \
    mob_suite

# ── Verify installation ───────────────────────────────────────────────────────
echo ""
echo "Verifying installation..."
conda run -p "$ENV_DIR" mob_recon --version

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "  MOBsuite installed successfully!"
    echo "  Completed: $(date)"
    echo "========================================"
else
    echo ""
    echo "ERROR: Installation may have failed — check the log"
    exit 1
fi
