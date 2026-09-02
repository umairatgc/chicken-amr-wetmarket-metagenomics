#!/bin/bash
#SBATCH --job-name=kma_index
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/kma_index_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/kma_index_%j.err

# =============================================================================
# build_kma_index.sh
# Build KMA binary index files for the ResFinder database
#
# WHAT THIS SCRIPT DOES:
#   KMA (the aligner used by ResFinder) cannot search raw .fsa text files.
#   It requires pre-built binary index files (.length.b, .seq.b, .name etc.)
#   A fresh git clone of resfinder_db only contains .fsa files -- not indices.
#   This script builds those indices once. All future ResFinder jobs reuse them.
#
# WHEN TO RUN:
#   - After cloning resfinder_db for the first time
#   - After pulling database updates (git pull)
#   - If ResFinder fails with: KMA Error: No such file or directory *.length.b
#
# INPUT:
#   databases/resfinder/resfinder_db/*.fsa   -- gene database files per drug class
#
# OUTPUT:
#   databases/resfinder/resfinder_db/*.length.b  -- KMA index files (one set per .fsa)
#
# SUBMISSION:
#   sbatch scripts/build_kma_index.sh
#   Monitor: tail -f logs/kma_index_JOBID.log
#   Check:   ls databases/resfinder/resfinder_db/*.length.b | wc -l
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"

BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
RESFINDER_DB=$BASE/databases/resfinder/resfinder_db
CONDA_ENV=$BASE/envs/resfinder

mkdir -p "$BASE/logs"

echo "======================================================"
echo "KMA database index builder"
echo "Start time: $(date)"
echo "Database:   $RESFINDER_DB"
echo "======================================================"

# ── Validate conda environment ────────────────────────────────────
if [ ! -d "$CONDA_ENV" ]; then
    echo "ERROR: ResFinder environment not found: $CONDA_ENV"
    exit 1
fi

conda activate "$CONDA_ENV"

if ! command -v kma &> /dev/null; then
    echo "ERROR: KMA not found in environment: $CONDA_ENV"
    exit 1
fi

KMA_VERSION=$(kma -v 2>&1 | head -1)
echo "KMA version: $KMA_VERSION"

# ── Validate database directory ───────────────────────────────────
if [ ! -d "$RESFINDER_DB" ]; then
    echo "ERROR: ResFinder database not found: $RESFINDER_DB"
    echo "       Clone it first from:"
    echo "       https://bitbucket.org/genomicepidemiology/resfinder_db.git"
    exit 1
fi

FSA_COUNT=$(ls "$RESFINDER_DB"/*.fsa 2>/dev/null | wc -l)
if [ "$FSA_COUNT" -eq 0 ]; then
    echo "ERROR: No .fsa files found in $RESFINDER_DB"
    exit 1
fi

echo "FSA files found: $FSA_COUNT"
echo ""
echo "Building KMA indices..."
echo ""

# ── Build KMA index for each .fsa file ───────────────────────────
KMA_ERRORS=0
for fsa in "$RESFINDER_DB"/*.fsa; do
    base=$(basename "$fsa" .fsa)
    echo "  Indexing ${base}.fsa..."
    kma index -i "$fsa" -o "$RESFINDER_DB/$base"
    if [ $? -ne 0 ]; then
        echo "  ERROR: Failed to index ${base}.fsa"
        KMA_ERRORS=$((KMA_ERRORS + 1))
    else
        echo "  Done: ${base}"
    fi
done

echo ""

if [ $KMA_ERRORS -gt 0 ]; then
    echo "ERROR: $KMA_ERRORS index builds failed -- check database files"
    exit 1
fi

# ── Verify indices were created ───────────────────────────────────
INDEX_COUNT=$(ls "$RESFINDER_DB"/*.length.b 2>/dev/null | wc -l)

echo "======================================================"
echo "KMA indexing complete: $(date)"
echo ""
echo "  FSA files indexed:   $FSA_COUNT"
echo "  Index files created: $INDEX_COUNT"
echo ""
echo "  ResFinder jobs can now run successfully."
echo "======================================================"

conda deactivate
