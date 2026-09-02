#!/bin/bash
#SBATCH --job-name=amr_databases
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/amr_databases_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/amr_databases_%j.err

# =============================================================================
# download_amr_databases.sh
# Downloads the latest databases for all four AMR tools:
#   1. AMRFinderPlus  — NCBI curated AMR gene database (v2026-05-15.1)
#   2. CARD / RGI     — Comprehensive Antibiotic Resistance Database (v4.0.1)
#   3. ResFinder      — DTU resistance gene database (latest git)
#   4. ABRicate       — 12 built-in databases (argannot, card, vfdb etc.)
#
# FIXES APPLIED:
#   - module load Anaconda3 (no version number) — loads full conda shell
#     function so conda activate/deactivate work correctly
#   - conda activate/deactivate instead of source activate/deactivate
#   - AMRFinderPlus: no --database flag with --update (causes error)
#   - ABRicate in separate environment due to RGI dependency conflicts
#   - Skip re-download if databases already exist (safe to rerun)
# =============================================================================

# ── Load conda ────────────────────────────────────────────────────────────────
# IMPORTANT: Use 'module load Anaconda3' WITHOUT version number.
# The versioned module (Anaconda3/2024.06-1) loads a stripped conda that
# breaks conda activate and conda deactivate. The unversioned module loads
# the full conda shell function correctly.
module load Anaconda3

# ── Paths ────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
DB_DIR=$BASE/databases
AMR_ENV=$BASE/envs/amr
ABRICATE_ENV=$BASE/envs/abricate

# Create database directories (no folder for amrfinderplus — it saves
# automatically inside the conda environment)
mkdir -p "$DB_DIR/card"
mkdir -p "$DB_DIR/resfinder"

echo "======================================================"
echo "AMR Database Download"
echo "Start time: $(date)"
echo "Saving databases to: $DB_DIR"
echo "======================================================"


# ==============================================================================
# STEP 1 — AMRFinderPlus database
# Latest version: 2026-05-15.1
# Downloads from: https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/
# Saves to: $AMR_ENV/share/amrfinderplus/data/2026-05-15.1/
#
# IMPORTANT: Do NOT use --database flag with --update — AMRFinderPlus does
#            not allow both flags together and will throw an error.
#            It saves the database to its own location inside the conda env.
# ==============================================================================
echo ""
echo "STEP 1: Downloading AMRFinderPlus database..."

conda activate "$AMR_ENV"

# --update automatically fetches the latest database version from NCBI
amrfinder --update

echo "AMRFinderPlus database downloaded."
echo "Location: $AMR_ENV/share/amrfinderplus/data/"

conda deactivate


# ==============================================================================
# STEP 2 — CARD database (for RGI)
# Latest version: 4.0.1 (released 2025-05-29)
# Downloads from: https://card.mcmaster.ca/latest/data
# Saves to: $DB_DIR/card/
#
# After downloading, rgi load indexes card.json so RGI can search it fast.
# If card.json already exists from a previous run, skips download and
# just re-runs rgi load to ensure the index is up to date.
# ==============================================================================
echo ""
echo "STEP 2: Downloading CARD database for RGI..."

conda activate "$AMR_ENV"

cd "$DB_DIR/card"

if [ ! -f "$DB_DIR/card/card.json" ]; then
    echo "Downloading CARD data bundle from McMaster University..."
    # -L follows redirects, -O saves with original filename
    wget -L -O card-data.tar.bz2 https://card.mcmaster.ca/latest/data
    # Extract the archive — produces card.json and supporting files
    tar -xjf card-data.tar.bz2
else
    echo "card.json already exists — skipping download"
fi

# Load card.json into RGI — indexes it for fast resistance gene searching
# --local tells RGI to use the current directory as the database location
rgi load --card_json "$DB_DIR/card/card.json" --local

echo "CARD database downloaded and loaded into RGI."
echo "Location: $DB_DIR/card/"

conda deactivate


# ==============================================================================
# STEP 3 — ResFinder database
# Latest version: current git HEAD (DTU Food Institute, Bitbucket)
# Downloads from: https://bitbucket.org/genomicepidemiology/resfinder_db.git
# Saves to: $DB_DIR/resfinder/resfinder_db/
#
# The database is a git repository of FASTA files (one per antibiotic class).
# If already cloned from a previous run, pulls latest updates instead.
# --depth 1 clones only the latest commit — faster and less disk space.
# ==============================================================================
echo ""
echo "STEP 3: Downloading ResFinder database..."

conda activate "$AMR_ENV"

if [ -d "$DB_DIR/resfinder/resfinder_db" ]; then
    echo "ResFinder already cloned — pulling latest updates..."
    cd "$DB_DIR/resfinder/resfinder_db"
    git pull
else
    echo "Cloning ResFinder database from Bitbucket..."
    cd "$DB_DIR/resfinder"
    git clone --depth 1 https://bitbucket.org/genomicepidemiology/resfinder_db.git
fi

echo "ResFinder database downloaded."
echo "Location: $DB_DIR/resfinder/resfinder_db/"

conda deactivate


# ==============================================================================
# STEP 4 — ABRicate databases
# Latest version: ABRicate 1.4.0 with 12 built-in databases
# Databases: argannot, bacmet2, card, ecoh, ecoli_vf, megares, ncbi,
#            plasmidfinder, resfinder, upec_expec_vf, vfdb, victors
#
# ABRicate is in a SEPARATE conda environment from RGI because ABRicate 1.4.0
# has dependency conflicts with RGI (openssl and perl version clashes).
# --setupdb downloads and indexes all 12 databases at once.
# ==============================================================================
echo ""
echo "STEP 4: Setting up ABRicate databases..."

conda activate "$ABRICATE_ENV"

# Confirm correct version before running
echo "ABRicate version: $(abricate --version)"

# Download and build BLAST indices for all 12 built-in databases
abricate --setupdb

# List all databases to confirm everything is installed correctly
echo ""
echo "ABRicate databases installed:"
abricate --list

conda deactivate


# ==============================================================================
# SUMMARY
# ==============================================================================
echo ""
echo "======================================================"
echo "All AMR databases downloaded successfully!"
echo "End time: $(date)"
echo ""
echo "Database locations:"
echo "  AMRFinderPlus : $AMR_ENV/share/amrfinderplus/data/"
echo "  CARD (RGI)    : $DB_DIR/card/"
echo "  ResFinder     : $DB_DIR/resfinder/resfinder_db/"
echo "  ABRicate      : inside $ABRICATE_ENV (managed by abricate)"
echo "======================================================"
