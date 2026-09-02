#!/bin/bash
#SBATCH --job-name=kraken2_assembly
#SBATCH --clusters=arc
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/kraken2_assembly_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/kraken2_assembly_%A_%a.err

# =============================================================================
# kraken2_assembly.sh
# Classifies polished metagenomic assemblies (contigs) with Kraken2.
# Produces per-contig species assignments for AMR gene annotation.
#
# BEFORE SUBMITTING:
#   1. Set KRAKEN2_DB to the correct database path (see below)
#   2. Run: sbatch kraken2_assembly.sh
#
# OUTPUTS (per sample):
#   kraken2_assembly/SampleName/kraken2_output.txt  — per-contig: C/U, contig, taxid
#   kraken2_assembly/SampleName/kraken2_report.txt  — summary report
# =============================================================================

KRAKEN2_DB="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/databases/kraken2"

# ── Project paths ──────────────────────────────────────────────────────────────
BASE_DIR="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
POLISHED_DIR="${BASE_DIR}/polished"
OUT_DIR="${BASE_DIR}/kraken2_assembly"
mkdir -p "${OUT_DIR}"

# ── Get sample list ────────────────────────────────────────────────────────────
SAMPLES=($(ls -d ${POLISHED_DIR}/*/))
SAMPLE_PATH="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"
SAMPLE_NAME=$(basename "${SAMPLE_PATH}")
ASSEMBLY="${SAMPLE_PATH}consensus.fasta"

echo "================================="
echo "Array task:  ${SLURM_ARRAY_TASK_ID}"
echo "Sample:      ${SAMPLE_NAME}"
echo "Assembly:    ${ASSEMBLY}"
echo "================================="

if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "ERROR: Assembly not found: ${ASSEMBLY}"
    exit 1
fi

# ── Output directory for this sample ──────────────────────────────────────────
SAMPLE_OUT="${OUT_DIR}/${SAMPLE_NAME}"
mkdir -p "${SAMPLE_OUT}"

# ── Load Kraken2 ──────────────────────────────────────────────────────────────
module purge
module load Kraken2

# ── Run Kraken2 on polished assembly ──────────────────────────────────────────
kraken2 \
    --db "${KRAKEN2_DB}" \
    --threads 8 \
    --output  "${SAMPLE_OUT}/kraken2_output.txt" \
    --report  "${SAMPLE_OUT}/kraken2_report.txt" \
    --use-names \
    "${ASSEMBLY}"

echo ""
echo "Done: ${SAMPLE_NAME}"

# ── Quick check ───────────────────────────────────────────────────────────────
TOTAL=$(wc -l < "${SAMPLE_OUT}/kraken2_output.txt")
CLASSIFIED=$(grep -c "^C" "${SAMPLE_OUT}/kraken2_output.txt" || true)
echo "Contigs classified: ${CLASSIFIED} / ${TOTAL}"
