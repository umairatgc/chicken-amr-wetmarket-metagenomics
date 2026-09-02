#!/bin/bash
#SBATCH --job-name=genomad
#SBATCH --clusters=arc
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/genomad_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/genomad_%A_%a.err

# =============================================================================
# run_genomad.sh
# Classifies polished assemblies into plasmid / virus / chromosome using geNomad.
#
# INPUT:  polished/SampleName/consensus.fasta
# OUTPUT: mobile_elements/SampleName/genomad/
#           <sample>_summary/
#             <sample>_plasmid_summary.tsv   — plasmid contigs + scores
#             <sample>_virus_summary.tsv     — virus contigs + scores
#           <sample>_aggregated_classification/
#             <sample>_aggregated_classification.tsv  — per-contig classification
#
# USAGE:
#   sbatch scripts/run_genomad.sh
# =============================================================================

BASE_DIR="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
POLISHED_DIR="${BASE_DIR}/polished"
OUT_DIR="${BASE_DIR}/mobile_elements"
DB="${BASE_DIR}/databases/genomad/genomad_db"
ENVS="${BASE_DIR}/envs"

# ── Get sample ─────────────────────────────────────────────────────────────────
SAMPLES=($(ls -d ${POLISHED_DIR}/*/))
SAMPLE_PATH="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"
SAMPLE_NAME=$(basename "${SAMPLE_PATH}")
ASSEMBLY="${SAMPLE_PATH}consensus.fasta"

echo "========================================"
echo "  geNomad"
echo "  Array task:  ${SLURM_ARRAY_TASK_ID}"
echo "  Sample:      ${SAMPLE_NAME}"
echo "  Started:     $(date)"
echo "========================================"

if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "ERROR: Assembly not found: ${ASSEMBLY}"
    exit 1
fi

SAMPLE_OUT="${OUT_DIR}/${SAMPLE_NAME}/genomad"
mkdir -p "${SAMPLE_OUT}"

conda run -p "${ENVS}/genomad" \
    genomad end-to-end \
        --threads 8 \
        --cleanup \
        "${ASSEMBLY}" \
        "${SAMPLE_OUT}" \
        "${DB}"

if [ $? -eq 0 ]; then
    echo ""
    echo "  ✅ ${SAMPLE_NAME} complete: $(date)"
    # Quick summary (tail -n +2 skips header row)
    PLASMIDS=$(tail -n +2 "${SAMPLE_OUT}/${SAMPLE_NAME}_summary/${SAMPLE_NAME}_plasmid_summary.tsv" 2>/dev/null | wc -l)
    VIRUSES=$(tail  -n +2 "${SAMPLE_OUT}/${SAMPLE_NAME}_summary/${SAMPLE_NAME}_virus_summary.tsv"   2>/dev/null | wc -l)
    echo "  Plasmid contigs: ${PLASMIDS}   Viral contigs: ${VIRUSES}"
else
    echo "  ❌ ${SAMPLE_NAME} FAILED — check log"
    exit 1
fi
