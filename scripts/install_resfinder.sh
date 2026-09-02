#!/bin/bash
#SBATCH --job-name=install_resfinder
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/install_resfinder_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/install_resfinder_%j.err

# =============================================================================
# install_resfinder.sh
# Install ResFinder environment from scratch
#
# ENVIRONMENT CREATED:
#   envs/resfinder  — ResFinder 4.7.2 + KMA aligner  (Python 3.11)
#
# DATABASE:
#   ResFinder DB  : databases/resfinder/resfinder_db/  (CGE, updated 2026-01-26)
#   PointFinder DB: databases/resfinder/pointfinder_db/ (CGE, updated 2025-11-05)
#
# KMA INDEXING:
#   KMA requires binary index files (.length.b etc.) built from the .fsa files.
#   A fresh git clone does NOT include these binary files — only text .fsa files.
#   This script builds the KMA indices automatically during installation.
#   This is a one-time step. All future ResFinder jobs reuse these indices.
#
# WHY SEPARATE ENVIRONMENT:
#   AMRFinderPlus and ResFinder are kept in separate environments so that
#   issues with one tool do not require reinstalling the other.
#
# WORKFLOW:
#   Step 1:  Ensure databases exist:
#              databases/resfinder/resfinder_db/   (git clone from CGE)
#              databases/resfinder/pointfinder_db/ (git clone from CGE)
#   Step 2:  sbatch scripts/install_resfinder.sh   (this script)
#   No login-node step needed — no internet required for ResFinder.
#
# SUBMISSION:
#   sbatch scripts/install_resfinder.sh
#   Monitor: tail -f logs/install_resfinder_JOBID.log
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"

BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
RESFINDER_ENV=$BASE/envs/resfinder
DB_RES=$BASE/databases/resfinder/resfinder_db
POINTFINDER_DB=$BASE/databases/resfinder/pointfinder_db

mkdir -p "$BASE/logs"

echo "=============================================="
echo "ResFinder Environment Installation"
echo "Start time: $(date)"
echo "=============================================="

# ── Validate databases exist before starting ──────────────────────
echo ""
echo "Checking databases..."

if [ ! -d "$DB_RES" ]; then
    echo "ERROR: ResFinder database not found: $DB_RES"
    echo "       Clone it first:"
    echo "       cd databases/resfinder"
    echo "       git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git"
    exit 1
fi
echo "ResFinder DB:   OK  ($DB_RES)"

if [ ! -d "$POINTFINDER_DB" ]; then
    echo "WARNING: PointFinder database not found: $POINTFINDER_DB"
    echo "         Point mutation detection will not be available."
    echo "         Clone it with:"
    echo "         cd databases/resfinder"
    echo "         git clone https://bitbucket.org/genomicepidemiology/pointfinder_db.git"
else
    echo "PointFinder DB: OK  ($POINTFINDER_DB)"
fi

# ── Remove existing resfinder env if present ─────────────────────
if [ -d "$RESFINDER_ENV" ]; then
    echo ""
    echo "Removing existing resfinder environment..."
    conda env remove --prefix "$RESFINDER_ENV" -y
    echo "Old resfinder environment removed."
fi

# ── Create environment ────────────────────────────────────────────
echo ""
echo "Creating resfinder environment (Python 3.11)..."
echo "  ResFinder 4.7.2"
echo "  KMA aligner (bundled with ResFinder)"
echo ""

conda create --prefix "$RESFINDER_ENV" \
    -c conda-forge -c bioconda \
    python=3.11 \
    resfinder=4.7.2 \
    -y

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create resfinder environment"
    exit 1
fi

conda activate "$RESFINDER_ENV"

# ── Verify ResFinder ──────────────────────────────────────────────
echo ""
echo "Verifying resfinder environment..."

RESFINDER_VERSION=$(python -m resfinder --version 2>&1 | head -1)
if [ -z "$RESFINDER_VERSION" ]; then
    echo "ERROR: ResFinder not working in environment"
    exit 1
fi
echo "ResFinder: OK  (version: $RESFINDER_VERSION)"

if ! command -v kma &> /dev/null; then
    echo "ERROR: KMA aligner not found in environment"
    exit 1
fi
KMA_VERSION=$(kma -v 2>&1 | head -1)
echo "KMA:       OK  (version: $KMA_VERSION)"

# ── Build KMA indices for ResFinder database ──────────────────────
# KMA searches binary index files, not raw .fsa text files.
# These must be built once from the .fsa files after cloning the database.
echo ""
echo "Building KMA indices for ResFinder database..."
echo "  Source: $DB_RES"
echo "  (One-time step — indices are reused by all future ResFinder jobs)"
echo ""

KMA_ERRORS=0
for fsa in "$DB_RES"/*.fsa; do
    base=$(basename "$fsa" .fsa)
    echo "  Indexing ${base}.fsa..."
    kma index -i "$fsa" -o "$DB_RES/$base"
    if [ $? -ne 0 ]; then
        echo "  ERROR: Failed to index ${base}.fsa"
        KMA_ERRORS=$((KMA_ERRORS + 1))
    fi
done

if [ $KMA_ERRORS -gt 0 ]; then
    echo "ERROR: $KMA_ERRORS KMA index builds failed"
    exit 1
fi

INDEX_COUNT=$(ls "$DB_RES"/*.length.b 2>/dev/null | wc -l)
echo ""
echo "KMA indexing complete.  Index files created: $INDEX_COUNT"

conda deactivate

# ── Final summary ─────────────────────────────────────────────────
echo ""
echo "=============================================="
echo "Installation complete: $(date)"
echo ""
echo "  envs/resfinder — ResFinder 4.7.2 + KMA  ✓"
echo "  KMA indices built for ResFinder database  ✓"
echo ""
echo "No further action required — ResFinder jobs can now run."
echo ""
echo "To update the ResFinder database in the future:"
echo "  cd $DB_RES"
echo "  git pull origin master"
echo "  Then re-run this script to rebuild KMA indices."
echo "=============================================="
