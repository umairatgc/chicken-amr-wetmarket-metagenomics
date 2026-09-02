#!/bin/bash
#SBATCH --job-name=mobsuite
#SBATCH --clusters=arc
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/mobsuite_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/mobsuite_%A_%a.err

# =============================================================================
# run_mobsuite.sh
# Classifies contigs as chromosome or plasmid and types plasmids using MOBsuite.
#
# INPUT:  polished/SampleName/consensus.fasta
# OUTPUT: mobile_elements/SampleName/mobsuite/
#           contig_report.txt    — per-contig: chromosome / plasmid / ambiguous
#           plasmid_report.txt   — per-plasmid: replicon type, mobility, GC%
#           chromosome.fasta     — chromosomal contigs
#           plasmid_*.fasta      — individual plasmid sequences
#
# KEY COLUMNS in contig_report.txt:
#   molecule_type      : chromosome | plasmid | ambiguous
#   primary_cluster_id : plasmid cluster ID
#   rep_type           : replicon type (e.g. IncF, IncI, Col)
#   mob_type           : mobility class (conjugative / mobilizable / non-mobilizable)
#
# USAGE:
#   sbatch scripts/run_mobsuite.sh
# =============================================================================

BASE_DIR="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
POLISHED_DIR="${BASE_DIR}/polished"
OUT_DIR="${BASE_DIR}/mobile_elements"
ENV_PATH="${BASE_DIR}/envs/mob_suite"

# ── Get sample ────────────────────────────────────────────────────────────────
SAMPLES=($(ls -d ${POLISHED_DIR}/*/))
SAMPLE_PATH="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"
SAMPLE_NAME=$(basename "${SAMPLE_PATH}")
ASSEMBLY="${SAMPLE_PATH}consensus.fasta"

echo "========================================"
echo "  MOBsuite — mob_recon"
echo "  Array task:  ${SLURM_ARRAY_TASK_ID}"
echo "  Sample:      ${SAMPLE_NAME}"
echo "  Started:     $(date)"
echo "========================================"

# ── Input guard ───────────────────────────────────────────────────────────────
if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "ERROR: Assembly not found: ${ASSEMBLY}"
    exit 1
fi

if [[ ! -d "${ENV_PATH}" ]]; then
    echo "ERROR: Conda environment not found: ${ENV_PATH}"
    exit 1
fi

# ── Output directory ──────────────────────────────────────────────────────────
SAMPLE_OUT="${OUT_DIR}/${SAMPLE_NAME}/mobsuite"
mkdir -p "${SAMPLE_OUT}"

# ── Run MOBsuite ──────────────────────────────────────────────────────────────
# Prepend conda env bin/ to PATH so mob_recon can find blastn/makeblastdb/tblastn.
# We invoke Python directly (not conda run) to bypass the broken shebang that was
# left behind when envs/mobsuite was renamed to envs/mob_suite.
export PATH="${ENV_PATH}/bin:${PATH}"

"${ENV_PATH}/bin/python3" "${ENV_PATH}/bin/mob_recon" \
    --infile      "${ASSEMBLY}" \
    --outdir      "${SAMPLE_OUT}" \
    --num_threads 8 \
    --force

# --force : overwrite existing output (safe for re-runs)

EXIT_CODE=$?
if [ ${EXIT_CODE} -eq 0 ]; then
    echo ""
    echo "  ✅ ${SAMPLE_NAME} complete: $(date)"
    # Count by molecule_type column (col 2), skipping header
    PLASMIDS=$(tail -n +2 "${SAMPLE_OUT}/contig_report.txt" 2>/dev/null | awk -F'\t' '$2=="plasmid"' | wc -l || echo 0)
    CHROMOS=$(tail  -n +2 "${SAMPLE_OUT}/contig_report.txt" 2>/dev/null | awk -F'\t' '$2=="chromosome"' | wc -l || echo 0)
    PLASMID_TYPES=$(tail -n +2 "${SAMPLE_OUT}/plasmid_report.txt" 2>/dev/null | wc -l || echo 0)
    echo "  Chromosome contigs: ${CHROMOS}   Plasmid contigs: ${PLASMIDS}   Plasmid types: ${PLASMID_TYPES}"
else
    echo "  ❌ ${SAMPLE_NAME} FAILED — check log"
    exit 1
fi
