#!/bin/bash
#SBATCH --job-name=install_amrfinder
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/install_amrfinder_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/install_amrfinder_%j.err

# =============================================================================
# install_amrfinder.sh
# Install AMRFinderPlus environment from scratch
#
# ENVIRONMENT CREATED:
#   envs/amrfinder  — AMRFinderPlus 4.2.7  (Python 3.11)
#
# DATABASE:
#   AMRFinderPlus database must be downloaded ON THE LOGIN NODE after this
#   job finishes. Compute nodes have no internet access.
#   ⚠️  Do NOT use --database flag with --update (causes error)
#   Database saves automatically to: envs/amrfinder/share/amrfinderplus/data/
#
# WHY SEPARATE ENVIRONMENT:
#   AMRFinderPlus and ResFinder are kept in separate environments so that
#   issues with one tool do not require reinstalling the other.
#
# WORKFLOW:
#   Step 1:  sbatch scripts/install_amrfinder.sh    (this script)
#   Step 2:  After job finishes, ON LOGIN NODE:
#              module load Anaconda3
#              eval "$(conda shell.bash hook)"
#              conda activate envs/amrfinder
#              amrfinder --update
#              conda deactivate
#
# SUBMISSION:
#   sbatch scripts/install_amrfinder.sh
#   Monitor: tail -f logs/install_amrfinder_JOBID.log
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"

BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
AMRFINDER_ENV=$BASE/envs/amrfinder

mkdir -p "$BASE/logs"

echo "=============================================="
echo "AMRFinderPlus Environment Installation"
echo "Start time: $(date)"
echo "=============================================="

# ── Remove existing amrfinder env if present ─────────────────────
if [ -d "$AMRFINDER_ENV" ]; then
    echo ""
    echo "Removing existing amrfinder environment..."
    conda env remove --prefix "$AMRFINDER_ENV" -y
    echo "Old amrfinder environment removed."
fi

# ── Create environment ────────────────────────────────────────────
echo ""
echo "Creating amrfinder environment (Python 3.11)..."
echo "  AMRFinderPlus 4.2.7"
echo ""

conda create --prefix "$AMRFINDER_ENV" \
    -c conda-forge -c bioconda \
    python=3.11 \
    ncbi-amrfinderplus=4.2.7 \
    -y

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create amrfinder environment"
    exit 1
fi

conda activate "$AMRFINDER_ENV"

# ── Verify AMRFinderPlus ──────────────────────────────────────────
echo ""
echo "Verifying amrfinder environment..."

AMRFINDER_VERSION=$(amrfinder --version 2>/dev/null)
if [ -z "$AMRFINDER_VERSION" ]; then
    echo "ERROR: amrfinder not working in environment"
    exit 1
fi
echo "AMRFinderPlus: OK  (version: $AMRFINDER_VERSION)"

conda deactivate

# ── Final summary ─────────────────────────────────────────────────
echo ""
echo "=============================================="
echo "Installation complete: $(date)"
echo ""
echo "  envs/amrfinder — AMRFinderPlus 4.2.7  ✓"
echo ""
echo "ACTION REQUIRED — run on LOGIN NODE after this job:"
echo "  module load Anaconda3"
echo "  eval \"\$(conda shell.bash hook)\""
echo "  conda activate $AMRFINDER_ENV"
echo "  amrfinder --update"
echo "  conda deactivate"
echo ""
echo "  Database will save to:"
echo "  $AMRFINDER_ENV/share/amrfinderplus/data/"
echo "=============================================="
