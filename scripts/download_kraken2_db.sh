#!/bin/bash
#SBATCH --job-name=kraken2_db
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/kraken2_db_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/kraken2_db_%j.err

# =============================================================================
# download_kraken2_db.sh
# Downloads the Kraken2 PlusPF database for taxonomic classification.
#
# What is this database?
#   Kraken2 works by comparing your reads against a large collection of
#   reference genomes. The database contains the reference sequences from:
#     - Bacteria        ← what you are studying
#     - Archaea         ← single-celled organisms without a nucleus
#     - Viruses         ← to detect any viral sequences
#     - Fungi           ← to detect fungal contamination
#     - Protozoa        ← single-celled parasites
#     - Human (hg38)    ← additional human reference
#     - Plasmids        ← mobile genetic elements carrying resistance genes
#   "PlusPF" = the standard Kraken2 database PLUS Protozoa and Fungi
#
# Why is it so large (~140 GB)?
#   The database contains millions of reference genomes. Each organism's
#   DNA is broken into short k-mers (short DNA words) and stored in a
#   hash table for ultra-fast lookup. The size reflects the breadth of
#   microbial diversity represented.
#
# Why 64 GB memory?
#   Building and loading the Kraken2 database requires holding large
#   portions of it in RAM. 64 GB is the recommended minimum.
#
# What files are produced?
#   databases/kraken2/hash.k2d    — the main k-mer hash table (~140 GB)
#   databases/kraken2/opts.k2d    — database options and parameters
#   databases/kraken2/taxo.k2d    — taxonomy information
#
# Runtime: ~6-12 hours depending on download speed and ARC load.
# Submit before leaving for the day and check results tomorrow.
# =============================================================================

# ── Load conda and activate taxonomy environment ──────────────────────────────
# The taxonomy environment contains Kraken2 and Bracken
module load Anaconda3
eval "$(conda shell.bash hook)"   # CRITICAL: enables conda activate in SLURM
conda activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/taxonomy

# ── Set database directory path ──────────────────────────────────────────────
DB=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/databases/kraken2

# Create the directory if it does not already exist
mkdir -p "$DB"

echo "======================================================"
echo "Kraken2 PlusPF Database Download"
echo "Start time: $(date)"
echo "Saving to: $DB"
echo "Estimated size: ~79.8 GB download / ~103.4 GB after extraction"
echo "======================================================"

# ==============================================================================
# STEP 1 — Download the pre-built Kraken2 PlusPF database
# We download the pre-built version from the Kraken2 authors' server rather
# than building it from scratch. Building from scratch would take days and
# require internet access to NCBI. The pre-built version is identical in
# content and ready to use immediately after extraction.
#
# wget flags:
#   --no-verbose    : suppresses most output but shows errors
#   --show-progress : shows a progress bar for the download
#   -O              : save with this specific filename
# ==============================================================================
echo ""
echo "------------------------------------------------------"
echo "STEP 1: Downloading Kraken2 PlusPF database..."
echo "This will take several hours. Started: $(date)"
echo "------------------------------------------------------"

wget --no-verbose --show-progress \
    https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20260226.tar.gz \
    -O "$DB/k2_pluspf.tar.gz"

echo "Download complete: $(date)"
echo "Downloaded file size: $(du -sh "$DB/k2_pluspf.tar.gz" | cut -f1)"

# ==============================================================================
# STEP 2 — Extract the database archive
# The downloaded file is a compressed archive (.tar.gz).
# tar -xvzf extracts it:
#   -x : extract files from archive
#   -v : verbose — shows each file as it is extracted
#   -z : decompress using gzip
#   -f : the archive filename follows
#   -C : extract into this directory
# ==============================================================================
echo ""
echo "------------------------------------------------------"
echo "STEP 2: Extracting database archive..."
echo "Started: $(date)"
echo "------------------------------------------------------"

tar -xvzf "$DB/k2_pluspf.tar.gz" -C "$DB"

echo "Extraction complete: $(date)"

# ==============================================================================
# STEP 3 — Remove the compressed archive to save disk space
# After extraction the .tar.gz file is no longer needed.
# The extracted files (hash.k2d, opts.k2d, taxo.k2d) are what Kraken2 uses.
# Keeping the archive would waste ~140 GB of additional disk space.
# ==============================================================================
echo ""
echo "------------------------------------------------------"
echo "STEP 3: Removing compressed archive to save disk space..."
echo "------------------------------------------------------"

rm "$DB/k2_pluspf.tar.gz"
echo "Archive removed."

# ==============================================================================
# STEP 4 — Verify the database is complete and usable
# kraken2-inspect reads the database and prints a summary of its contents.
# If this command runs without error, the database is valid and ready to use.
# ==============================================================================
echo ""
echo "------------------------------------------------------"
echo "STEP 4: Verifying database integrity..."
echo "------------------------------------------------------"

kraken2-inspect \
    --db "$DB" \
    --report-zero-counts 2>/dev/null | head -20

echo ""
echo "Database verification complete."

# ==============================================================================
# SUMMARY
# ==============================================================================
echo ""
echo "======================================================"
echo "Kraken2 database download complete!"
echo "End time: $(date)"
echo ""
echo "Database location: $DB"
echo "Database contents:"
ls -lh "$DB"
echo ""
echo "Next steps:"
echo "  1. Verify host depletion is complete (63 files in host_depleted/)"
echo "  2. Submit run_kraken2.sh to classify all 63 samples"
echo "======================================================"

conda deactivate
