#!/bin/bash
#SBATCH --job-name=minimap2_coverage
#SBATCH --clusters=arc
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/minimap2_coverage_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/minimap2_coverage_%A_%a.err

# =============================================================================
# run_minimap2_coverage.sh
# Phase 4a — Map clean reads (coverage for binning, not host depletion) back to polished assembly to generate coverage depth
#
# WHY THIS STEP IS NEEDED:
#   MetaBAT2 bins contigs using two signals: (1) tetranucleotide composition
#   and (2) coverage depth. Coverage depth = how many reads landed on each
#   contig. Contigs from the same organism will have similar depth because
#   that organism is present at the same abundance throughout the sample.
#   This script generates the depth file that MetaBAT2 requires.
#
# WORKFLOW PER SAMPLE:
#   1. minimap2  — aligns clean reads to polished assembly → SAM file
#   2. samtools sort — sorts alignments by genomic position → BAM file
#   3. samtools index — indexes the BAM for fast access
#   4. jgi_summarize_bam_contig_depths — calculates per-contig depth → depth.txt
#   5. SAM file deleted after BAM is created (large intermediate file)
#
# INPUT PER SAMPLE:
#   clean_reads/SampleName_clean.fastq.gz   — quality-filtered ONT reads
#   polished/SampleName/consensus.fasta     — polished assembly (reference)
#
# OUTPUT PER SAMPLE:
#   bins/SampleName/coverage/SampleName_sorted.bam      — sorted alignment
#   bins/SampleName/coverage/SampleName_sorted.bam.bai  — BAM index
#   bins/SampleName/coverage/SampleName_depth.txt       — per-contig depth
#
# CONDA ENV: metabat2 (contains minimap2, samtools, jgi_summarize_bam_contig_depths)
#
# RUN BEFORE: run_metabat2.sh
# RUN AFTER:  Phase 3 complete (polished/SampleName/consensus.fasta must exist)
#
# SUBMISSION:
#   sbatch scripts/run_minimap2.sh
#   Monitor:  squeue -u YOUR_USERNAME
#   Check done: ls bins/*/coverage/*_depth.txt 2>/dev/null | wc -l
# =============================================================================

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
CLEAN_DIR="${BASE_DIR}/clean_reads"
POLISHED_DIR="${BASE_DIR}/polished"
BINS_DIR="${BASE_DIR}/bins"
SAMPLE_LIST="${BASE_DIR}/r_analysis/decontam/sample_list_real.txt"
ENV_PATH="${BASE_DIR}/envs/metabat2"

# ── Get sample name from list using SLURM array task ID ───────────────────────
# The sample list has one sample name per line.
# Task 1 → line 1 → first sample. Task 61 → line 61 → last sample.
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

# ── Build input/output paths ───────────────────────────────────────────────────
READS="${CLEAN_DIR}/${SAMPLE_NAME}_clean.fastq.gz"
ASSEMBLY="${POLISHED_DIR}/${SAMPLE_NAME}/consensus.fasta"
SAMPLE_OUT="${BINS_DIR}/${SAMPLE_NAME}/coverage"

echo "========================================"
echo "  minimap2 + samtools — Coverage Depth"
echo "  Array task:  ${SLURM_ARRAY_TASK_ID}"
echo "  Sample:      ${SAMPLE_NAME}"
echo "  Started:     $(date)"
echo "========================================"

# ── Guard: check sample name was found ────────────────────────────────────────
if [[ -z "${SAMPLE_NAME}" ]]; then
    echo "ERROR: No sample found at line ${SLURM_ARRAY_TASK_ID} in ${SAMPLE_LIST}"
    exit 1
fi

# ── Guard: check reads file exists ────────────────────────────────────────────
if [[ ! -f "${READS}" ]]; then
    echo "ERROR: Clean reads not found: ${READS}"
    exit 1
fi

# ── Guard: check assembly exists ──────────────────────────────────────────────
if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "ERROR: Polished assembly not found: ${ASSEMBLY}"
    exit 1
fi

# ── Guard: check conda environment exists ─────────────────────────────────────
if [[ ! -d "${ENV_PATH}" ]]; then
    echo "ERROR: Conda environment not found: ${ENV_PATH}"
    exit 1
fi

# ── Create output directory ────────────────────────────────────────────────────
mkdir -p "${SAMPLE_OUT}"

# ── Define output file paths ───────────────────────────────────────────────────
# NOTE: No SAM file — minimap2 is piped directly into samtools sort.
# This avoids writing a huge intermediate SAM file and prevents the
# "truncated file" error caused by conda run stdout buffering with > redirect.
BAM="${SAMPLE_OUT}/${SAMPLE_NAME}_sorted.bam"    # final sorted alignment
DEPTH="${SAMPLE_OUT}/${SAMPLE_NAME}_depth.txt"   # per-contig depth for MetaBAT2

# ── Steps 1+2: minimap2 piped directly into samtools sort ────────────────────
# Running both tools inside a single 'conda run bash -c' ensures the pipe is
# fully in-process — no file redirect race condition, no intermediate SAM file.
#
# minimap2 flags:
#   -a         : output SAM format (piped to samtools, not written to disk)
#   -x map-ont : Nanopore long-read preset
#   -t 8       : CPU threads (minimap2 uses -t, not --threads)
#
# samtools sort flags:
#   -@ 8 : sort threads
#   -o   : output BAM file
#   -    : read from stdin (the pipe from minimap2)
echo ""
echo "  [1+2/4] minimap2 → samtools sort (piped, no intermediate SAM)..."
echo "          Reference: $(basename ${ASSEMBLY})"
echo "          Reads:     $(basename ${READS})"

conda run -p "${ENV_PATH}" bash -c "
    minimap2 -a -x map-ont -t 8 '${ASSEMBLY}' '${READS}' | \
    samtools sort -@ 8 -o '${BAM}' -
"

if [[ $? -ne 0 ]]; then
    echo "  ERROR: minimap2 | samtools sort failed for ${SAMPLE_NAME}"
    rm -f "${BAM}"
    exit 1
fi
echo "  [1+2/4] minimap2 + samtools sort complete (no SAM written to disk)"

# ── Step 3: samtools index — index the BAM file ────────────────────────────────
# Creates a .bai index file alongside the BAM.
# Required for random access into the BAM (jgi_summarize needs this).
echo ""
echo "  [3/4] samtools index — indexing BAM..."


conda run -p "${ENV_PATH}" \
    samtools index "${BAM}"

if [[ $? -ne 0 ]]; then
    echo "  ERROR: samtools index failed for ${SAMPLE_NAME}"
    exit 1
fi
echo "  [3/4] samtools index complete"

# ── Step 4: jgi_summarize_bam_contig_depths — calculate per-contig depth ──────
# This tool is part of MetaBAT2. It reads the sorted BAM and for every contig
# calculates: mean coverage depth, length, and variance.
# The output depth.txt is the exact format MetaBAT2 expects as input.
echo ""
echo "  [4/4] jgi_summarize_bam_contig_depths — calculating depth..."

conda run -p "${ENV_PATH}" \
    jgi_summarize_bam_contig_depths \
        --outputDepth "${DEPTH}" \
        "${BAM}"

if [[ $? -ne 0 ]]; then
    echo "  ERROR: jgi_summarize_bam_contig_depths failed for ${SAMPLE_NAME}"
    exit 1
fi
echo "  [4/4] Depth calculation complete"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
CONTIG_COUNT=$(tail -n +2 "${DEPTH}" 2>/dev/null | wc -l | tr -d ' ')
MEAN_DEPTH=$(tail -n +2 "${DEPTH}" 2>/dev/null | awk '{sum+=$3; count++} END {if(count>0) printf "%.1f", sum/count; else print "?"}')

echo "  ✅ ${SAMPLE_NAME} complete: $(date)"
echo "  Contigs in depth file: ${CONTIG_COUNT}"
echo "  Mean coverage depth:   ${MEAN_DEPTH}x"
echo "  Depth file:            ${DEPTH}"
echo ""
echo "  Next step: sbatch scripts/run_metabat2.sh"
