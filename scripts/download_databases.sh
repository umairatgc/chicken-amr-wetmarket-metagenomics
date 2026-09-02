#!/bin/bash
#SBATCH --job-name=download_dbs
#SBATCH --clusters=arc
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/download_dbs.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/download_dbs.err

# =============================================================================
# download_databases.sh
# Downloads all required databases for the MGE/MAG pipeline phases.
#
# DATABASES DOWNLOADED:
#   1. geNomad    ~2 GB   — required for Phase 3 (plasmid/virus classification)
#   2. CheckM2    ~3 GB   — required for Phase 4 (MAG quality assessment)
#   3. Bakta      ~60 GB  — required for Phase 4 (MAG gene annotation)
#
# ALREADY AVAILABLE — no download needed:
#   GTDB-Tk DB  — databases/gtdbtk/  (check version below)
#   Kraken2 DB  — databases/kraken2/
#   CARD DB     — bundled with RGI env
#   ResFinder   — databases/resfinder/
#
# NOTE: Bakta is large (~60 GB). If time limit is hit, re-submit this script
#       — it will skip already-completed downloads.
#
# USAGE:
#   sbatch scripts/download_databases.sh
# =============================================================================

module purge
module load Anaconda3

BASE_DIR="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
DB_DIR="${BASE_DIR}/databases"
ENVS="${BASE_DIR}/envs"

# =============================================================================
# CHECK GTDB-Tk DATABASE VERSION (informational)
# =============================================================================
echo "========================================"
echo "  GTDB-Tk database version check"
echo "========================================"
GTDBTK_DB="${DB_DIR}/gtdbtk"
if [ -f "${GTDBTK_DB}/metadata/genome_metadata.tsv" ]; then
    echo "  Checking genome_metadata.tsv..."
    head -1 "${GTDBTK_DB}/metadata/genome_metadata.tsv"
elif [ -f "${GTDBTK_DB}/metadata/metadata.tsv" ]; then
    echo "  Checking metadata.tsv..."
    head -1 "${GTDBTK_DB}/metadata/metadata.tsv"
else
    echo "  Metadata file not found — listing taxonomy folder:"
    ls "${GTDBTK_DB}/taxonomy/" 2>/dev/null || echo "  No taxonomy folder"
fi
echo ""

# =============================================================================
# 1. geNomad database (~2 GB)
# =============================================================================
echo "========================================"
echo "  1. Downloading geNomad database (~2 GB)"
echo "     Started: $(date)"
echo "========================================"

GENOMAD_DB="${DB_DIR}/genomad"

if [ -d "${GENOMAD_DB}/genomad_db" ]; then
    echo "  geNomad database already exists — skipping"
    ls -lh "${GENOMAD_DB}/genomad_db/" | head -5
else
    mkdir -p "${GENOMAD_DB}"
    conda run -p "${ENVS}/genomad" \
        genomad download-database "${GENOMAD_DB}"

    if [ $? -eq 0 ]; then
        echo "  ✅ geNomad database downloaded successfully"
    else
        echo "  ❌ geNomad database download FAILED"
    fi
fi

echo ""

# =============================================================================
# 2. CheckM2 database (~3 GB)
# =============================================================================
echo "========================================"
echo "  2. Downloading CheckM2 database (~3 GB)"
echo "     Started: $(date)"
echo "========================================"

CHECKM2_DB="${DB_DIR}/checkm2"

if [ -f "${CHECKM2_DB}/uniref100.KO.1.dmnd" ]; then
    echo "  CheckM2 database already exists — skipping"
    ls -lh "${CHECKM2_DB}/" | head -5
else
    mkdir -p "${CHECKM2_DB}"
    conda run -p "${ENVS}/checkm2" \
        checkm2 database --download --path "${CHECKM2_DB}/"

    if [ $? -eq 0 ]; then
        echo "  ✅ CheckM2 database downloaded successfully"
    else
        echo "  ❌ CheckM2 database download FAILED"
    fi
fi

echo ""

# =============================================================================
# 3. Bakta database (~60 GB — takes longest)
# =============================================================================
echo "========================================"
echo "  3. Downloading Bakta full database (~60 GB)"
echo "     Started: $(date)"
echo "  WARNING: This is large — may take 1-2 hours"
echo "========================================"

BAKTA_DB="${DB_DIR}/bakta"

if [ -d "${BAKTA_DB}" ] && [ "$(ls -A ${BAKTA_DB} 2>/dev/null)" ]; then
    echo "  Bakta database directory exists and is non-empty:"
    ls -lh "${BAKTA_DB}/" | head -10
    echo ""
    echo "  If download was incomplete, delete ${BAKTA_DB} and re-submit."
    echo "  Skipping to avoid overwriting partial data."
else
    mkdir -p "${BAKTA_DB}"
    conda run -p "${ENVS}/bakta" \
        bakta_db download --type full \
        --output "${BAKTA_DB}/"

    if [ $? -eq 0 ]; then
        echo "  ✅ Bakta database downloaded successfully"
        echo "  Size: $(du -sh ${BAKTA_DB}/ | cut -f1)"
    else
        echo "  ❌ Bakta database download FAILED or incomplete"
        echo "  Check the log and re-submit if needed"
    fi
fi

echo ""

# =============================================================================
# SUMMARY
# =============================================================================
echo "========================================"
echo "  Download summary: $(date)"
echo "========================================"
echo ""
echo "  Database locations:"
echo "    geNomad : ${DB_DIR}/genomad/genomad_db"
echo "    CheckM2 : ${DB_DIR}/checkm2"
echo "    Bakta   : ${DB_DIR}/bakta"
echo "    GTDB-Tk : ${DB_DIR}/gtdbtk  (pre-existing)"
echo ""
echo "  Set this variable before running GTDB-Tk jobs:"
echo "    export GTDBTK_DATA_PATH=${DB_DIR}/gtdbtk"
echo ""
du -sh "${DB_DIR}"/*/  2>/dev/null
echo "========================================"
