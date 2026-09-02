#!/bin/bash
#SBATCH --job-name=integronfinder
#SBATCH --clusters=arc
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/integronfinder_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/integronfinder_%A_%a.err

# =============================================================================
# run_integronfinder.sh
# Detects integrons in polished assemblies using IntegronFinder 2.
#
# INPUT:  polished/SampleName/consensus.fasta
# OUTPUT: mobile_elements/SampleName/integronfinder/
#           Results_consensus/
#             consensus.integrons   — integron summary table
#             consensus.gbk         — GenBank annotation
#           (output named after input file, not sample name)
#
# USAGE:
#   sbatch scripts/run_integronfinder.sh
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
echo "  IntegronFinder 2"
echo "  Array task:  ${SLURM_ARRAY_TASK_ID}"
echo "  Sample:      ${SAMPLE_NAME}"
echo "  Started:     $(date)"
echo "========================================"

if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "ERROR: Assembly not found: ${ASSEMBLY}"
    exit 1
fi

SAMPLE_OUT="${OUT_DIR}/${SAMPLE_NAME}/integronfinder"
mkdir -p "${SAMPLE_OUT}"

conda run -p "${ENVS}/integronfinder" \
    integron_finder \
        --cpu 8 \
        --outdir "${SAMPLE_OUT}" \
        --linear \
        --calin-threshold 2 \
        "${ASSEMBLY}"

# --linear          : treats contigs as linear (replaces --metagenome in this version)
#                     disables circularity assumption — appropriate for assembled contigs
# --calin-threshold 2 : minimum attC sites to report CALIN elements (default 2)
# NOTE: no --prefix flag in IntegronFinder 2; output folder named after input file (consensus)

if [ $? -eq 0 ]; then
    echo ""
    echo "  ✅ ${SAMPLE_NAME} complete: $(date)"
    INTEGRONS=$(grep -v "^#" "${SAMPLE_OUT}/Results_Integron_Finder_consensus/consensus.integrons" 2>/dev/null | wc -l || echo 0)
    echo "  Integron entries found: ${INTEGRONS}"
else
    echo "  ❌ ${SAMPLE_NAME} FAILED — check log"
    exit 1
fi
