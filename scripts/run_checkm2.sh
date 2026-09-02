#!/bin/bash
#SBATCH --job-name=checkm2
#SBATCH --clusters=arc
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/checkm2_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/checkm2_%A_%a.err

# =============================================================================
# run_checkm2.sh
# Phase 4c — Assess MAG quality: completeness and contamination
#
# WHY THIS STEP IS NEEDED:
#   MetaBAT2 bins contigs but has no way to know how good each bin is.
#   A bin might be:
#     - Missing large chunks of the genome (low completeness)
#     - Contaminated with contigs from another organism (high contamination)
#   CheckM2 uses machine-learning with universal marker genes to estimate:
#     - Completeness  (%) — what fraction of expected genes are present
#     - Contamination (%) — how much foreign DNA is in the bin
#
#   Standard quality thresholds (MIMAG standards):
#     High quality MAG:   completeness ≥ 90%,  contamination < 5%
#     Medium quality MAG: completeness ≥ 50%,  contamination < 10%
#     Low quality MAG:    completeness < 50%   (not used for downstream analysis)
#
# INPUT PER SAMPLE:
#   bins/SampleName/metabat2/          — folder containing bin.*.fa files
#
# OUTPUT PER SAMPLE:
#   bins/SampleName/checkm2/quality_report.tsv   — completeness & contamination
#   bins/SampleName/checkm2/diamond_output/       — intermediate DIAMOND results
#
# CONDA ENV: checkm2
# DATABASE:  databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd
#            (DIAMOND database of universal marker genes)
#            ⚠️  Update CHECKM2_DB path below if your database is elsewhere
#
# RUN AFTER:  run_metabat2.sh
# RUN BEFORE: run_gtdbtk.sh (only high-quality bins proceed to taxonomy)
#
# SUBMISSION:
#   sbatch scripts/run_checkm2.sh
#   Check done: ls bins/*/checkm2/quality_report.tsv 2>/dev/null | wc -l
# =============================================================================

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
BINS_DIR="${BASE_DIR}/bins"
SAMPLE_LIST="${BASE_DIR}/r_analysis/decontam/sample_list_real.txt"
ENV_PATH="${BASE_DIR}/envs/checkm2"

# Database confirmed at: databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd
CHECKM2_DB="${BASE_DIR}/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd"

# ── Get sample name ────────────────────────────────────────────────────────────
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

# ── Build input/output paths ───────────────────────────────────────────────────
BIN_DIR="${BINS_DIR}/${SAMPLE_NAME}/metabat2"
SAMPLE_OUT="${BINS_DIR}/${SAMPLE_NAME}/checkm2"

echo "========================================"
echo "  CheckM2 — MAG Quality Assessment"
echo "  Array task:  ${SLURM_ARRAY_TASK_ID}"
echo "  Sample:      ${SAMPLE_NAME}"
echo "  Started:     $(date)"
echo "========================================"

# ── Guards ────────────────────────────────────────────────────────────────────
if [[ -z "${SAMPLE_NAME}" ]]; then
    echo "ERROR: No sample found at line ${SLURM_ARRAY_TASK_ID} in ${SAMPLE_LIST}"
    exit 1
fi

if [[ ! -d "${BIN_DIR}" ]]; then
    echo "ERROR: MetaBAT2 bin directory not found: ${BIN_DIR}"
    echo "       Run run_metabat2.sh first"
    exit 1
fi

# Count bins — skip sample if no bins were produced
BIN_COUNT=$(ls "${BIN_DIR}"/bin.*.fa 2>/dev/null | wc -l)
if [[ ${BIN_COUNT} -eq 0 ]]; then
    echo "WARNING: No bins found in ${BIN_DIR}"
    echo "         MetaBAT2 may have produced no bins for this sample (low coverage or very fragmented assembly)"
    echo "         Skipping CheckM2 for ${SAMPLE_NAME}"
    exit 0
fi
echo "  Bins to assess: ${BIN_COUNT}"

if [[ ! -f "${CHECKM2_DB}" ]]; then
    echo "ERROR: CheckM2 database not found: ${CHECKM2_DB}"
    echo "       Download with: checkm2 database --download"
    exit 1
fi

if [[ ! -d "${ENV_PATH}" ]]; then
    echo "ERROR: Conda environment not found: ${ENV_PATH}"
    exit 1
fi

# ── Create output directory ────────────────────────────────────────────────────
mkdir -p "${SAMPLE_OUT}"

# ── Run CheckM2 ───────────────────────────────────────────────────────────────
# predict         : the CheckM2 subcommand for quality prediction
# --input         : directory containing the bin FASTA files
# --output-directory : where to write quality_report.tsv and intermediate files
# --threads       : CPU threads (CheckM2 runs DIAMOND which is highly parallelisable)
# --database_path : path to the DIAMOND marker gene database
# --extension     : file extension of the bin files (fa, not fasta or fna)
# --force         : overwrite previous output if re-running
echo ""
echo "  Running CheckM2 on ${BIN_COUNT} bins..."

conda run -p "${ENV_PATH}" \
    checkm2 predict \
        --input              "${BIN_DIR}" \
        --output-directory   "${SAMPLE_OUT}" \
        --threads            16 \
        --database_path      "${CHECKM2_DB}" \
        --extension          fa \
        --force

EXIT_CODE=$?

if [[ ${EXIT_CODE} -ne 0 ]]; then
    echo ""
    echo "  ❌ ${SAMPLE_NAME} FAILED — check log"
    exit 1
fi

# ── Summary — count bins passing quality thresholds ───────────────────────────
QUALITY_REPORT="${SAMPLE_OUT}/quality_report.tsv"

if [[ -s "${QUALITY_REPORT}" ]]; then
    # High quality: completeness >= 90 AND contamination < 5
    HQ=$(tail -n +2 "${QUALITY_REPORT}" | awk -F'\t' '$2 >= 90 && $3 < 5' | wc -l)
    # Medium quality: completeness >= 50 AND contamination < 10
    MQ=$(tail -n +2 "${QUALITY_REPORT}" | awk -F'\t' '$2 >= 50 && $3 < 10 && !($2 >= 90 && $3 < 5)' | wc -l)
    # Low quality: everything else
    LQ=$(tail -n +2 "${QUALITY_REPORT}" | awk -F'\t' '$2 < 50 || $3 >= 10' | wc -l)

    echo ""
    echo "  ✅ ${SAMPLE_NAME} complete: $(date)"
    echo "  Total bins assessed:  ${BIN_COUNT}"
    echo "  High quality  (≥90% complete, <5% contam):  ${HQ}"
    echo "  Medium quality (≥50% complete, <10% contam): ${MQ}"
    echo "  Low quality (<50% complete or ≥10% contam):  ${LQ}"
    echo ""
    echo "  Quality report: ${QUALITY_REPORT}"
    echo "  Next step: sbatch scripts/run_gtdbtk.sh"
else
    echo "  ❌ ${SAMPLE_NAME} FAILED — quality_report.tsv not produced"
    exit 1
fi
