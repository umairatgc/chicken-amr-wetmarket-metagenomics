#!/bin/bash
#SBATCH --job-name=host_genomes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/host_genomes_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/host_genomes_%j.err

# =============================================================================
# download_host_genomes.sh
# Downloads and prepares reference genomes for host depletion (Step 2.4).
#
# What is host depletion?
#   Your DNA samples contain DNA from the bacteria you want to study, BUT also
#   DNA from the host organism (chicken/human). If you don't remove host DNA,
#   it will appear in your results and confuse the analysis. This script
#   downloads the reference genomes so minimap2 can identify and remove them.
#
# Which hosts?
#   1. Gallus gallus (chicken) GRCg7b — Ensembl release 111
#      Used because samples may contain chicken-origin DNA
#   2. Homo sapiens (human) GRCh38 — Ensembl release 111
#      Used to remove any human DNA from researcher handling
#
# Why concatenate both genomes?
#   minimap2 can screen against both hosts in a single pass if they are
#   combined into one file — faster and simpler than running twice.
#
# What is a minimap2 index (.mmi)?
#   minimap2 needs to search through ~4 billion base pairs of host DNA for
#   every read. An index pre-processes the genome into a format that allows
#   very fast searching — like the index at the back of a textbook vs reading
#   every page. Building the index takes ~30-60 minutes but only needs to be
#   done ONCE. After that, every host depletion job uses this index directly.
#
# Output files:
#   databases/host_genomes/gallus_gallus.fa.gz    — chicken genome
#   databases/host_genomes/homo_sapiens.fa.gz     — human genome
#   databases/host_genomes/host_combined.fa.gz    — both concatenated
#   databases/host_genomes/host_combined.mmi      — minimap2 index (the one you use)
#
# Run BEFORE submitting run_host_depletion.sh
# =============================================================================

# ── Load conda and activate the qc environment ───────────────────────────────
# The qc environment contains minimap2 which is needed to build the index
module load Anaconda3
eval "$(conda shell.bash hook)"   # required for conda activate to work in SLURM
conda activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/qc

# ── Set database directory path ──────────────────────────────────────────────
# All host genome files will be stored here
DB=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/databases/host_genomes

# Create the directory if it does not already exist
mkdir -p "$DB"

echo "======================================================"
echo "Host Genome Download and Indexing"
echo "Start time: $(date)"
echo "Saving to: $DB"
echo "======================================================"


# ==============================================================================
# STEP 1 — Download Gallus gallus (chicken) reference genome
# Source: Ensembl release 111, assembly GRCg7b
# File size: ~1.1 GB compressed
# -q           : quiet mode (suppresses verbose wget output)
# --show-progress : still shows a progress bar so you can see it is running
# -O           : save the file with this specific name
# ==============================================================================
echo ""
echo "------------------------------------------------------"
echo "STEP 1: Downloading Gallus gallus (chicken) genome..."
echo "------------------------------------------------------"

wget -q --show-progress \
    https://ftp.ensembl.org/pub/release-111/fasta/gallus_gallus/dna/Gallus_gallus.GRCg7b.dna.toplevel.fa.gz \
    -O "$DB/gallus_gallus.fa.gz"

echo "Gallus gallus download complete."
echo "File size: $(du -sh "$DB/gallus_gallus.fa.gz" | cut -f1)"


# ==============================================================================
# STEP 2 — Download Homo sapiens (human) reference genome
# Source: Ensembl release 111, assembly GRCh38
# File size: ~3.1 GB compressed
# ==============================================================================
echo ""
echo "------------------------------------------------------"
echo "STEP 2: Downloading Homo sapiens (human) genome..."
echo "------------------------------------------------------"

wget -q --show-progress \
    https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
    -O "$DB/homo_sapiens.fa.gz"

echo "Homo sapiens download complete."
echo "File size: $(du -sh "$DB/homo_sapiens.fa.gz" | cut -f1)"


# ==============================================================================
# STEP 3 — Concatenate both genomes into one combined reference file
# cat joins files together end-to-end.
# The > symbol saves the output to a new file instead of printing to screen.
# Result: one file containing all chicken chromosomes + all human chromosomes.
# minimap2 will align reads against this combined file in one pass.
# ==============================================================================
echo ""
echo "------------------------------------------------------"
echo "STEP 3: Combining both genomes into one file..."
echo "------------------------------------------------------"

cat "$DB/gallus_gallus.fa.gz" "$DB/homo_sapiens.fa.gz" > "$DB/host_combined.fa.gz"

echo "Combined genome created."
echo "File size: $(du -sh "$DB/host_combined.fa.gz" | cut -f1)"


# ==============================================================================
# STEP 4 — Build minimap2 index from the combined reference
# This is the most time-consuming step (~30-60 minutes).
# The index only needs to be built ONCE — all future jobs use the .mmi file.
#
# -x map-ont  : use the ONT (Oxford Nanopore) preset — tells minimap2 to
#               expect long reads with ONT-typical error rates
# -d          : output the index to this file (.mmi format)
# -t 8        : use 8 CPU threads to speed up index building
# ==============================================================================
echo ""
echo "------------------------------------------------------"
echo "STEP 4: Building minimap2 index (30-60 minutes)..."
echo "------------------------------------------------------"

minimap2 \
    -x map-ont \
    -d "$DB/host_combined.mmi" \
    "$DB/host_combined.fa.gz" \
    -t 8

echo "minimap2 index built successfully."
echo "Index size: $(du -sh "$DB/host_combined.mmi" | cut -f1)"


# ==============================================================================
# SUMMARY
# ==============================================================================
echo ""
echo "======================================================"
echo "Host genome setup complete!"
echo "End time: $(date)"
echo ""
echo "Files created:"
echo "  Chicken genome : $DB/gallus_gallus.fa.gz"
echo "  Human genome   : $DB/homo_sapiens.fa.gz"
echo "  Combined       : $DB/host_combined.fa.gz"
echo "  Index (use this): $DB/host_combined.mmi"
echo ""
echo "Next step: submit run_host_depletion.sh"
echo "======================================================"

conda deactivate
