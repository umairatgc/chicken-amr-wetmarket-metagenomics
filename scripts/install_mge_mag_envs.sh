#!/bin/bash
#SBATCH --job-name=install_mge_mag
#SBATCH --clusters=arc
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/install_all_envs.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/install_all_envs.err

# =============================================================================
# install_all_envs.sh
# Installs all remaining conda environments for the AMR pipeline.
#
# Environments created:
#   Phase 2: hamronization  — hAMRonization (harmonise AMR outputs)
#   Phase 3: genomad        — geNomad (plasmid + virus)
#             integronfinder — IntegronFinder 2 (integrons)
#             isescan        — ISEScan (insertion sequences)
#   Phase 4: metabat2       — MetaBAT2 + minimap2 (MAG binning)
#             checkm2        — CheckM2 (bin QC)
#             gtdbtk         — GTDB-Tk (MAG taxonomy)
#             bakta          — Bakta (gene annotation)
#
# NOTE: MOBsuite is handled separately by install_mobsuite.sh
#
# LARGE DATABASES (download separately after install — see bottom of script):
#   geNomad DB   ~2 GB   — genomad download-database
#   CheckM2 DB   ~3 GB   — checkm2 database --download
#   GTDB-Tk DB   ~110 GB — check if already on ARC shared storage first
#   Bakta DB     ~60 GB  — bakta_db --type full --output ...
#
# USAGE:
#   sbatch scripts/install_all_envs.sh
# =============================================================================

module purge
module load Anaconda3

ENVS="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs"

install_env() {
  local name=$1; shift
  local path="${ENVS}/${name}"
  echo ""
  echo "========================================"
  echo "  Installing: ${name}"
  echo "  Path: ${path}"
  echo "  Started: $(date)"
  echo "========================================"
  [ -d "$path" ] && conda env remove -p "$path" -y
  conda create -y -p "$path" "$@"
  if [ $? -eq 0 ]; then
    echo "  ✅ ${name} installed successfully"
  else
    echo "  ❌ ${name} FAILED — check log"
  fi
}

# =============================================================================
# PHASE 2 — hAMRonization
# =============================================================================
install_env hamronization \
    -c conda-forge -c bioconda \
    hamronization

# =============================================================================
# PHASE 3 — Mobile Genetic Elements
# =============================================================================
install_env genomad \
    -c conda-forge -c bioconda \
    genomad

install_env integronfinder \
    -c conda-forge -c bioconda \
    integron_finder

install_env isescan \
    -c conda-forge -c bioconda \
    isescan

# =============================================================================
# PHASE 4 — MAG Binning & Taxonomy
# =============================================================================
# MetaBAT2 + minimap2 (minimap2 needed to generate coverage for binning)
install_env metabat2 \
    -c conda-forge -c bioconda \
    metabat2 minimap2 samtools

# CheckM2
install_env checkm2 \
    -c conda-forge -c bioconda \
    checkm2

# GTDB-Tk
install_env gtdbtk \
    -c conda-forge -c bioconda \
    gtdbtk

# Bakta
install_env bakta \
    -c conda-forge -c bioconda \
    bakta

# =============================================================================
# SUMMARY
# =============================================================================
echo ""
echo "========================================"
echo "  All installations complete: $(date)"
echo "  Verifying versions..."
echo "========================================"
echo ""

verify() {
  local name=$1; local cmd=$2
  result=$(conda run -p "${ENVS}/${name}" $cmd 2>&1 | head -1)
  echo "  ${name}: ${result}"
}

verify hamronization   "hamronize --version"
verify genomad         "genomad --version"
verify integronfinder  "integron_finder --version"
verify isescan         "isescan.py --version 2>&1 || isescan --version"
verify metabat2        "metabat2"
verify checkm2         "checkm2 --version"
verify gtdbtk          "gtdbtk --version"
verify bakta           "bakta --version"

echo ""
echo "========================================"
echo "  NEXT STEPS — Download databases"
echo "========================================"
echo ""
echo "  1. geNomad database (~2 GB):"
echo "     conda run -p ${ENVS}/genomad \\"
echo "       genomad download-database databases/"
echo ""
echo "  2. CheckM2 database (~3 GB):"
echo "     conda run -p ${ENVS}/checkm2 \\"
echo "       checkm2 database --download --path databases/checkm2/"
echo ""
echo "  3. GTDB-Tk database — ALREADY AVAILABLE at:"
echo "     /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/databases/gtdbtk"
echo "     Set environment variable before running GTDB-Tk jobs:"
echo "     export GTDBTK_DATA_PATH=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/databases/gtdbtk"
echo ""
echo "  4. Bakta database (~60 GB) — NOT YET downloaded:"
echo "     conda run -p ${ENVS}/bakta \\"
echo "       bakta_db --type full \\"
echo "       --output /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/databases/bakta/"
echo "========================================"
