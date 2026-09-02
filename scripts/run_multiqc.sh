#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/multiqc_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/multiqc_%j.err

# =============================================================================
# run_multiqc.sh
# Step 2.1 (Part 2) — Aggregate all NanoPlot QC reports into one summary.
#
# What does MultiQC do?
#   When NanoPlot ran, it produced a separate HTML report for each of the
#   63 samples — stored in 63 subfolders inside the qc/ directory. That means
#   you would need to open 63 separate files to get an overview of all samples.
#
#   MultiQC solves this by automatically scanning the qc/ folder, finding all
#   63 NanoPlot reports, reading the statistics from each one, and combining
#   them into a single interactive HTML file. You open one file and see all
#   63 samples side by side — making it easy to compare and spot outliers.
#
# How does MultiQC find the NanoPlot reports?
#   You simply point MultiQC at the qc/ folder. It then:
#     1. Scans every subfolder inside qc/ automatically
#     2. Recognises NanoPlot output files by their filename patterns
#        (e.g. *_NanoStats.txt, *_NanoPlot-report.html)
#     3. Reads the statistics from each file
#     4. Combines everything into one summary report
#   You do not need to list the samples individually — MultiQC finds them all.
#
# Input:  qc/                        (63 subfolders, one per sample)
# Output: multiqc_report/multiqc_report.html  (one combined HTML file)
#
# Run AFTER run_qc.sh and run_qc_cont.sh have both completed.
# Verify: ls qc/ | wc -l  → should be 63 before running this
# =============================================================================

# ── Load conda and activate qc environment ───────────────────────────────────
# MultiQC is installed in the qc environment
module load Anaconda3
eval "$(conda shell.bash hook)"   # CRITICAL: enables conda activate in SLURM
conda activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/qc

# ── Set base path ─────────────────────────────────────────────────────────────
# BASE is the root of your project folder — used to build all other paths
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1

echo "======================================================"
echo "MultiQC — Aggregating all 63 NanoPlot reports"
echo "Start time: $(date)"
echo "Scanning: $BASE/qc/"
echo "======================================================"

# ── Run MultiQC ───────────────────────────────────────────────────────────────
# MultiQC arguments explained:
#
#   $BASE/qc/
#     The folder to scan. MultiQC searches this folder and ALL subfolders
#     inside it recursively. It automatically finds and reads every NanoPlot
#     output file it encounters. This is how it collects all 63 reports
#     without you needing to list them individually.
#
#   --outdir $BASE/multiqc_report
#     Where to save the output HTML file. The folder is created automatically
#     if it does not already exist.
#
#   --title "ONT Metagenomics QC — All 63 Samples"
#     The title displayed at the top of the HTML report — helps identify
#     the report when you open it later.
#
#   --filename multiqc_report.html
#     The name of the output HTML file. Open this in your browser to view
#     all 63 samples aggregated into one interactive summary.

multiqc $BASE/qc/ \
    --outdir $BASE/multiqc_report \
    --title "ONT Metagenomics QC — All 63 Samples" \
    --filename multiqc_report.html

echo ""
echo "======================================================"
echo "MultiQC complete!"
echo "End time: $(date)"
echo ""
echo "Output saved to: $BASE/multiqc_report/multiqc_report.html"
echo ""
echo "Download and open in your browser:"
echo "  scp -r YOUR_USERNAME@arc-login.arc.ox.ac.uk:$BASE/multiqc_report ~/Desktop/"
echo "======================================================"

conda deactivate
