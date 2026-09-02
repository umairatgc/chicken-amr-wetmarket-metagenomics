#!/bin/bash
#SBATCH --job-name=gtdbtk
#SBATCH --clusters=arc
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/gtdbtk_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/gtdbtk_%A_%a.err

# =============================================================================
# run_gtdbtk.sh
# Phase 4d — Assign taxonomy to medium/high quality MAGs using GTDB-Tk
#
# WHY THIS STEP IS NEEDED:
#   MetaBAT2 bins have no names — they are just bin.1.fa, bin.2.fa etc.
#   GTDB-Tk compares each bin's genome sequence against the Genome Taxonomy
#   Database (GTDB) to assign a species/genus/family label.
#   This tells you WHICH BACTERIUM each MAG represents.
#
# IMPORTANT — MEMORY:
#   GTDB-Tk loads the entire reference database into RAM during the run.
#   This requires ~200GB of RAM. The medium partition is used here.
#   ⚠️  Do NOT run on the short partition or the job will fail with OOM.
#
# QUALITY FILTER:
#   Only medium and high quality bins are passed to GTDB-Tk.
#   Low quality bins (completeness < 50% OR contamination >= 10%) are skipped.
#   This is read automatically from the CheckM2 quality_report.tsv.
#
# INPUT PER SAMPLE:
#   bins/SampleName/metabat2/bin.*.fa                — MAG bin files
#   bins/SampleName/checkm2/quality_report.tsv       — CheckM2 results (for filtering)
#
# OUTPUT PER SAMPLE:
#   taxonomy_mags/SampleName/gtdbtk.bac120.summary.tsv  — bacterial taxonomy
#   taxonomy_mags/SampleName/gtdbtk.ar53.summary.tsv    — archaeal taxonomy
#   taxonomy_mags/SampleName/classify/                  — detailed placement files
#
# CONDA ENV: gtdbtk
# DATABASE:  databases/gtdbtk/  (GTDB reference data ~80GB)
#            ⚠️  Update GTDBTK_DATA_PATH below if your database is elsewhere
#
# RUN AFTER:  run_checkm2.sh
# RUN BEFORE: run_bakta.sh (annotation of quality bins)
#
# SUBMISSION:
#   sbatch scripts/run_gtdbtk.sh
#   Check done: ls taxonomy_mags/*/gtdbtk.bac120.summary.tsv 2>/dev/null | wc -l
# =============================================================================

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
BINS_DIR="${BASE_DIR}/bins"
TAX_DIR="${BASE_DIR}/taxonomy_mags"
SAMPLE_LIST="${BASE_DIR}/r_analysis/decontam/sample_list_real.txt"
ENV_PATH="${BASE_DIR}/envs/gtdbtk"

# Database confirmed at: databases/gtdbtk/ (contains markers, masks, metadata, msa, pplacer etc.)
export GTDBTK_DATA_PATH="${BASE_DIR}/databases/gtdbtk"

# Quality thresholds — bins below these values are skipped
MIN_COMPLETENESS=50    # minimum completeness % to pass to GTDB-Tk
MAX_CONTAMINATION=10   # maximum contamination % to pass to GTDB-Tk

# ── Get sample name ────────────────────────────────────────────────────────────
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

# ── Build input/output paths ───────────────────────────────────────────────────
BIN_DIR="${BINS_DIR}/${SAMPLE_NAME}/metabat2"
QUALITY_REPORT="${BINS_DIR}/${SAMPLE_NAME}/checkm2/quality_report.tsv"
SAMPLE_OUT="${TAX_DIR}/${SAMPLE_NAME}"
FILTERED_BINS_DIR="${SAMPLE_OUT}/quality_bins"   # temp dir for quality-filtered bins

echo "========================================"
echo "  GTDB-Tk — MAG Taxonomy Assignment"
echo "  Array task:  ${SLURM_ARRAY_TASK_ID}"
echo "  Sample:      ${SAMPLE_NAME}"
echo "  Started:     $(date)"
echo "  Memory:      200G (GTDB-Tk loads full database into RAM)"
echo "========================================"

# ── Guards ────────────────────────────────────────────────────────────────────
if [[ -z "${SAMPLE_NAME}" ]]; then
    echo "ERROR: No sample found at line ${SLURM_ARRAY_TASK_ID} in ${SAMPLE_LIST}"
    exit 1
fi

if [[ ! -f "${QUALITY_REPORT}" ]]; then
    echo "ERROR: CheckM2 quality report not found: ${QUALITY_REPORT}"
    echo "       Run run_checkm2.sh first"
    exit 1
fi

if [[ ! -d "${ENV_PATH}" ]]; then
    echo "ERROR: Conda environment not found: ${ENV_PATH}"
    exit 1
fi

if [[ ! -d "${GTDBTK_DATA_PATH}" ]]; then
    echo "ERROR: GTDB-Tk database not found: ${GTDBTK_DATA_PATH}"
    echo "       Download with: gtdbtk download_db --db_version r220 --download_path ${GTDBTK_DATA_PATH}"
    exit 1
fi

# ── Create output and filtered bins directories ────────────────────────────────
mkdir -p "${SAMPLE_OUT}"
# Remove and recreate filtered_bins dir on every run to avoid stale symlinks
# from a previous failed/partial run pointing at bins that no longer qualify
rm -rf "${FILTERED_BINS_DIR}"
mkdir -p "${FILTERED_BINS_DIR}"

# ── Filter bins by quality — only pass medium/high quality bins ────────────────
# Read CheckM2 quality_report.tsv (tab-separated, header on line 1)
# Column 1 = bin name, Column 2 = completeness, Column 3 = contamination
echo ""
echo "  Filtering bins by quality (completeness ≥${MIN_COMPLETENESS}%, contamination <${MAX_CONTAMINATION}%)..."

QUALITY_BINS=()
while IFS=$'\t' read -r BIN_NAME COMPLETENESS CONTAMINATION REST; do
    # Skip header line
    [[ "${BIN_NAME}" == "Name" ]] && continue

    # Check quality thresholds using awk for float comparison
    PASSES=$(echo "${COMPLETENESS} ${CONTAMINATION}" | \
             awk -v minc="${MIN_COMPLETENESS}" -v maxcon="${MAX_CONTAMINATION}" \
             '{if ($1 >= minc && $2 < maxcon) print "yes"; else print "no"}')

    if [[ "${PASSES}" == "yes" ]]; then
        BIN_FILE="${BIN_DIR}/${BIN_NAME}.fa"
        if [[ -f "${BIN_FILE}" ]]; then
            # Symlink the quality bin into filtered_bins dir for GTDB-Tk input
            ln -sf "${BIN_FILE}" "${FILTERED_BINS_DIR}/${BIN_NAME}.fa"
            QUALITY_BINS+=("${BIN_NAME}")
            echo "    ✅ ${BIN_NAME}  completeness=${COMPLETENESS}%  contamination=${CONTAMINATION}%"
        fi
    else
        echo "    ⏭  ${BIN_NAME}  completeness=${COMPLETENESS}%  contamination=${CONTAMINATION}% — below threshold, skipped"
    fi
done < "${QUALITY_REPORT}"

QUALITY_BIN_COUNT=${#QUALITY_BINS[@]}
echo ""
echo "  Quality bins passing filter: ${QUALITY_BIN_COUNT}"

# If no bins passed quality filter, exit gracefully (not an error)
if [[ ${QUALITY_BIN_COUNT} -eq 0 ]]; then
    echo "  WARNING: No bins passed quality filter for ${SAMPLE_NAME}"
    echo "           This sample will not have taxonomy assignments"
    echo "           This is not an error — some samples genuinely yield no quality MAGs"
    exit 0
fi

# ── Run GTDB-Tk classify_wf ───────────────────────────────────────────────────
# classify_wf : the main workflow — runs marker gene identification,
#               placement into reference tree, and taxonomy assignment
# --genome_dir : directory containing the quality-filtered bin FASTA files
# --out_dir    : output directory
# --cpus       : number of CPU threads
# --extension  : file extension of the bin files
echo ""
echo "  Running GTDB-Tk on ${QUALITY_BIN_COUNT} quality bins..."

conda run -p "${ENV_PATH}" \
    gtdbtk classify_wf \
        --genome_dir        "${FILTERED_BINS_DIR}" \
        --out_dir           "${SAMPLE_OUT}" \
        --cpus              16 \
        --extension         fa

EXIT_CODE=$?

if [[ ${EXIT_CODE} -ne 0 ]]; then
    echo ""
    echo "  ❌ ${SAMPLE_NAME} FAILED — check log"
    exit 1
fi

# ── Summary ───────────────────────────────────────────────────────────────────
BAC_SUMMARY="${SAMPLE_OUT}/classify/gtdbtk.bac120.summary.tsv"
ARC_SUMMARY="${SAMPLE_OUT}/classify/gtdbtk.ar53.summary.tsv"

BAC_COUNT=0
ARC_COUNT=0
[[ -s "${BAC_SUMMARY}" ]] && BAC_COUNT=$(tail -n +2 "${BAC_SUMMARY}" | wc -l)
[[ -s "${ARC_SUMMARY}" ]] && ARC_COUNT=$(tail -n +2 "${ARC_SUMMARY}" | wc -l)

# ── Validate output was actually produced ─────────────────────────────────────
if [[ ${BAC_COUNT} -eq 0 && ${ARC_COUNT} -eq 0 ]]; then
    echo ""
    echo "  ⚠️  WARNING: GTDB-Tk produced no taxonomy assignments for ${SAMPLE_NAME}"
    echo "     ${QUALITY_BIN_COUNT} bins were submitted but neither bac120 nor ar53 summary has any rows."
    echo "     This may indicate a pplacer failure — check the log for pplacer errors."
    # Don't exit 1 here — zero assignments can legitimately happen with unusual organisms,
    # but the warning will be visible in the log for morning review.
fi

echo ""
echo "  ✅ ${SAMPLE_NAME} complete: $(date)"
echo "  Quality bins processed:   ${QUALITY_BIN_COUNT}"
echo "  Bacterial MAGs assigned:  ${BAC_COUNT}"
echo "  Archaeal MAGs assigned:   ${ARC_COUNT}"
echo ""
echo "  Results: ${SAMPLE_OUT}"
echo "  Next step: run_bakta.sh (can run in parallel with this)"
