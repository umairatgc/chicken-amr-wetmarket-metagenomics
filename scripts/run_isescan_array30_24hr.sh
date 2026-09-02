#!/bin/bash
#SBATCH --job-name=isescan
#SBATCH --clusters=arc
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --array=30
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/isescan_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/isescan_%A_%a.err

# =============================================================================
# run_isescan.sh
# Identifies insertion sequences (IS elements) in polished assemblies using ISEScan.
#
# INPUT:  polished/SampleName/consensus.fasta
# OUTPUT: mobile_elements/SampleName/isescan/
#           <sample>.fasta.tsv    — IS element summary table
#           <sample>.fasta.gff    — GFF3 annotation of IS elements
#           <sample>.fasta.is.fna — IS element nucleotide sequences
#
# USAGE:
#   sbatch scripts/run_isescan.sh
# =============================================================================

BASE_DIR="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
POLISHED_DIR="${BASE_DIR}/polished"
OUT_DIR="${BASE_DIR}/mobile_elements"
ENVS="${BASE_DIR}/envs"

# ── Get sample ─────────────────────────────────────────────────────────────────
SAMPLES=($(ls -d ${POLISHED_DIR}/*/))
SAMPLE_PATH="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"
SAMPLE_NAME=$(basename "${SAMPLE_PATH}")
ASSEMBLY="${SAMPLE_PATH}consensus.fasta"

echo "========================================"
echo "  ISEScan"
echo "  Array task:  ${SLURM_ARRAY_TASK_ID}"
echo "  Sample:      ${SAMPLE_NAME}"
echo "  Started:     $(date)"
echo "========================================"

if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "ERROR: Assembly not found: ${ASSEMBLY}"
    exit 1
fi

SAMPLE_OUT="${OUT_DIR}/${SAMPLE_NAME}/isescan"
mkdir -p "${SAMPLE_OUT}"

# ISEScan writes output AND .list temp files relative to the current working directory
# cd to output dir so all files land there, not in the submission directory
cd "${SAMPLE_OUT}"

# Copy assembly to output dir so results land here cleanly
cp "${ASSEMBLY}" "${SAMPLE_OUT}/${SAMPLE_NAME}.fasta"

conda run -p "${ENVS}/isescan" \
    isescan.py \
        --seqfile "${SAMPLE_OUT}/${SAMPLE_NAME}.fasta" \
        --output  "${SAMPLE_OUT}" \
        --nthread 32

if [ $? -eq 0 ]; then
    echo ""
    echo "  ✅ ${SAMPLE_NAME} complete: $(date)"
    IS_COUNT=$(grep -v "^#" "${SAMPLE_OUT}/${SAMPLE_NAME}.fasta.tsv" 2>/dev/null | wc -l || echo 0)
    echo "  IS elements found: ${IS_COUNT}"
    # Remove the copied assembly to save space
    rm -f "${SAMPLE_OUT}/${SAMPLE_NAME}.fasta"
else
    echo "  ❌ ${SAMPLE_NAME} FAILED — check log"
    rm -f "${SAMPLE_OUT}/${SAMPLE_NAME}.fasta"
    exit 1
fi
