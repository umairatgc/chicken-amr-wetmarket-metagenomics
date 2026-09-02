#!/bin/bash
#SBATCH --job-name=metabat2
#SBATCH --clusters=arc
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/metabat2_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/metabat2_%A_%a.err

# =============================================================================
# run_metabat2.sh
# Phase 4b — Bin contigs into Metagenome-Assembled Genomes (MAGs)
#
# WHY THIS STEP IS NEEDED:
#   The polished assembly (consensus.fasta) contains contigs from ALL organisms
#   in the sample mixed together. MetaBAT2 groups contigs that likely came from
#   the same organism into "bins" (MAGs) using two signals:
#     1. Tetranucleotide composition — same organism = similar DNA signature
#     2. Coverage depth — same organism = similar abundance across the sample
#   The depth file from run_minimap2.sh provides signal 2.
#
# INPUT PER SAMPLE:
#   polished/SampleName/consensus.fasta               — the mixed assembly
#   bins/SampleName/coverage/SampleName_depth.txt     — per-contig depth
#
# OUTPUT PER SAMPLE:
#   bins/SampleName/metabat2/bin.1.fa   — MAG 1 (one organism's contigs)
#   bins/SampleName/metabat2/bin.2.fa   — MAG 2
#   bins/SampleName/metabat2/bin.N.fa   — MAG N ...
#   bins/SampleName/metabat2/unbinned.fa — contigs too short/ambiguous to bin
#
# CONDA ENV: metabat2
#
# RUN AFTER:  run_minimap2.sh  (depth file must exist)
# RUN BEFORE: run_checkm2.sh   (assess quality of these bins)
#
# SUBMISSION:
#   sbatch scripts/run_metabat2.sh
#   Check done: ls bins/*/metabat2/bin.*.fa 2>/dev/null | wc -l
# =============================================================================

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
POLISHED_DIR="${BASE_DIR}/polished"
BINS_DIR="${BASE_DIR}/bins"
SAMPLE_LIST="${BASE_DIR}/r_analysis/decontam/sample_list_real.txt"
ENV_PATH="${BASE_DIR}/envs/metabat2"

# ── Get sample name ────────────────────────────────────────────────────────────
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

# ── Build input/output paths ───────────────────────────────────────────────────
ASSEMBLY="${POLISHED_DIR}/${SAMPLE_NAME}/consensus.fasta"
DEPTH="${BINS_DIR}/${SAMPLE_NAME}/coverage/${SAMPLE_NAME}_depth.txt"
SAMPLE_OUT="${BINS_DIR}/${SAMPLE_NAME}/metabat2"

echo "========================================"
echo "  MetaBAT2 — Contig Binning"
echo "  Array task:  ${SLURM_ARRAY_TASK_ID}"
echo "  Sample:      ${SAMPLE_NAME}"
echo "  Started:     $(date)"
echo "========================================"

# ── Guards ────────────────────────────────────────────────────────────────────
if [[ -z "${SAMPLE_NAME}" ]]; then
    echo "ERROR: No sample found at line ${SLURM_ARRAY_TASK_ID} in ${SAMPLE_LIST}"
    exit 1
fi

if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "ERROR: Assembly not found: ${ASSEMBLY}"
    exit 1
fi

if [[ ! -f "${DEPTH}" ]]; then
    echo "ERROR: Depth file not found: ${DEPTH}"
    echo "       Run run_minimap2.sh first"
    exit 1
fi

if [[ ! -d "${ENV_PATH}" ]]; then
    echo "ERROR: Conda environment not found: ${ENV_PATH}"
    exit 1
fi

# ── Create output directory ────────────────────────────────────────────────────
mkdir -p "${SAMPLE_OUT}"

# ── Run MetaBAT2 ───────────────────────────────────────────────────────────────
# -i  : input assembly FASTA (all contigs mixed)
# -a  : depth file from jgi_summarize_bam_contig_depths (coverage signal)
# -o  : output prefix — MetaBAT2 will create bin.1.fa, bin.2.fa etc.
# -m  : minimum contig length to bin (default 2500 bp)
#        Shorter contigs do not have enough signal for reliable binning.
#        Contigs below this threshold go to unbinned.fa
# --numThreads : CPU threads to use
# --unbinned   : write unbinned contigs to a separate file (useful to know what was skipped)
echo ""
echo "  Running MetaBAT2..."
echo "  Assembly: ${ASSEMBLY}"
echo "  Depth:    ${DEPTH}"
echo "  Output:   ${SAMPLE_OUT}"

conda run -p "${ENV_PATH}" \
    metabat2 \
        -i  "${ASSEMBLY}" \
        -a  "${DEPTH}" \
        -o  "${SAMPLE_OUT}/bin" \
        -m  2500 \
        --numThreads 8 \
        --unbinned

EXIT_CODE=$?

if [[ ${EXIT_CODE} -ne 0 ]]; then
    echo ""
    echo "  ❌ ${SAMPLE_NAME} FAILED — check log"
    exit 1
fi

# ── Summary ───────────────────────────────────────────────────────────────────
BIN_COUNT=$(ls "${SAMPLE_OUT}"/bin.*.fa 2>/dev/null | wc -l)
TOTAL_CONTIGS=$(grep -c "^>" "${ASSEMBLY}" 2>/dev/null || echo "?")

# Count contigs that ended up in bins (binned) vs total
BINNED_CONTIGS=0
for BIN in "${SAMPLE_OUT}"/bin.*.fa; do
    if [[ -f "${BIN}" ]]; then
        COUNT=$(grep -c "^>" "${BIN}" 2>/dev/null || echo 0)
        BINNED_CONTIGS=$((BINNED_CONTIGS + COUNT))
    fi
done

echo ""
echo "  ✅ ${SAMPLE_NAME} complete: $(date)"
echo "  Total contigs in assembly: ${TOTAL_CONTIGS}"
echo "  MAG bins produced:         ${BIN_COUNT}"
echo "  Contigs placed in bins:    ${BINNED_CONTIGS}"
echo ""
echo "  Next step: sbatch scripts/run_checkm2.sh"
