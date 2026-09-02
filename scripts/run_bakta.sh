#!/bin/bash
#SBATCH --job-name=bakta
#SBATCH --clusters=arc
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/bakta_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/bakta_%A_%a.err

# =============================================================================
# run_bakta.sh
# Phase 4e — Annotate quality MAGs with Bakta
#
# WHY THIS STEP IS NEEDED:
#   After binning and quality filtering, we have MAG genomes but no information
#   about what genes they contain. Bakta annotates each MAG by identifying:
#     - Coding sequences (CDS) and their protein functions
#     - rRNA genes (16S, 23S, 5S)
#     - tRNA and tmRNA genes
#     - CRISPR arrays
#     - Signal peptides and transmembrane helices
#   This is needed for Phase 5: to know which ARGs are on which MAG and in
#   what genomic context they sit (e.g. next to a transposase = on a MGE).
#
# HOW THIS SCRIPT DIFFERS FROM THE OTHERS:
#   Other scripts process one file per sample (one assembly, one depth file).
#   Bakta must run on each BIN separately — one sample may have 5 bins,
#   another may have 20. This script loops over all quality bins within a sample.
#   One SLURM task = one sample = loops over all its quality bins.
#
# QUALITY FILTER:
#   Only bins passing CheckM2 thresholds are annotated:
#     Completeness >= 50%  AND  Contamination < 10%
#   Annotating low-quality bins wastes compute time and produces unreliable results.
#
# INPUT PER SAMPLE:
#   bins/SampleName/metabat2/bin.*.fa               — MAG bin files
#   bins/SampleName/checkm2/quality_report.tsv      — CheckM2 quality filter
#
# OUTPUT PER BIN:
#   annotation/SampleName/bakta_chromosome/BinName/BinName.gff3  — genome annotation (GFF3)
#   annotation/SampleName/bakta_chromosome/BinName/BinName.gbff  — GenBank flat file
#   annotation/SampleName/bakta_chromosome/BinName/BinName.faa   — annotated protein sequences
#   annotation/SampleName/bakta_chromosome/BinName/BinName.tsv   — tab-delimited feature table
#
# OUTPUT DIRECTORY STRUCTURE:
#   annotation/SampleName/
#     bakta_chromosome/   — chromosomal MAG annotations (this script)
#       bin.1/
#       bin.2/
#       ...
#     bakta_plasmid/      — plasmid annotations (separate script)
#
# CONDA ENV: bakta
# DATABASE:  databases/bakta/db/  (Bakta reference database ~30GB)
#            ⚠️  Update BAKTA_DB path below if your database is elsewhere
#
# RUN AFTER:  run_checkm2.sh (quality report must exist)
# RUN BEFORE: Phase 5 — ARG-MGE co-localisation
#
# SUBMISSION:
#   sbatch scripts/run_bakta.sh
#   Check done: ls annotation/*/bakta_chromosome/*/bin.*.gff3 2>/dev/null | wc -l
# =============================================================================

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
BINS_DIR="${BASE_DIR}/bins"
ANNOT_DIR="${BASE_DIR}/annotation"
SAMPLE_LIST="${BASE_DIR}/r_analysis/decontam/sample_list_real.txt"
ENV_PATH="${BASE_DIR}/envs/bakta"

# Database confirmed at: databases/bakta/db/
BAKTA_DB="${BASE_DIR}/databases/bakta/db"

# Quality thresholds — same as used in run_gtdbtk.sh
MIN_COMPLETENESS=50
MAX_CONTAMINATION=10

# ── Get sample name ────────────────────────────────────────────────────────────
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

# ── Build input paths ──────────────────────────────────────────────────────────
BIN_DIR="${BINS_DIR}/${SAMPLE_NAME}/metabat2"
QUALITY_REPORT="${BINS_DIR}/${SAMPLE_NAME}/checkm2/quality_report.tsv"
SAMPLE_ANNOT="${ANNOT_DIR}/${SAMPLE_NAME}"

# Chromosomal bin annotations go into bakta_chromosome/ — parallel to bakta_plasmid/
CHROM_ANNOT="${SAMPLE_ANNOT}/bakta_chromosome"

echo "========================================"
echo "  Bakta — MAG Genome Annotation"
echo "  Array task:  ${SLURM_ARRAY_TASK_ID}"
echo "  Sample:      ${SAMPLE_NAME}"
echo "  Started:     $(date)"
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

if [[ ! -d "${BAKTA_DB}" ]]; then
    echo "ERROR: Bakta database not found: ${BAKTA_DB}"
    echo "       Download with: bakta_db download --output ${BASE_DIR}/databases/bakta"
    exit 1
fi

if [[ ! -d "${ENV_PATH}" ]]; then
    echo "ERROR: Conda environment not found: ${ENV_PATH}"
    exit 1
fi

# ── Create sample annotation directory ────────────────────────────────────────
mkdir -p "${CHROM_ANNOT}"

# ── Filter bins by quality and annotate each passing bin ──────────────────────
echo ""
echo "  Filtering bins by quality (completeness ≥${MIN_COMPLETENESS}%, contamination <${MAX_CONTAMINATION}%)..."

PASS_COUNT=0
FAIL_COUNT=0
ANNOT_SUCCESS=0
ANNOT_FAIL=0

while IFS=$'\t' read -r BIN_NAME COMPLETENESS CONTAMINATION REST; do
    # Skip header line
    [[ "${BIN_NAME}" == "Name" ]] && continue

    # Check quality thresholds
    PASSES=$(echo "${COMPLETENESS} ${CONTAMINATION}" | \
             awk -v minc="${MIN_COMPLETENESS}" -v maxcon="${MAX_CONTAMINATION}" \
             '{if ($1 >= minc && $2 < maxcon) print "yes"; else print "no"}')

    if [[ "${PASSES}" != "yes" ]]; then
        echo "  ⏭  Skipping ${BIN_NAME} (completeness=${COMPLETENESS}%, contamination=${CONTAMINATION}%) — below threshold"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    BIN_FILE="${BIN_DIR}/${BIN_NAME}.fa"
    if [[ ! -f "${BIN_FILE}" ]]; then
        echo "  WARNING: Bin file not found: ${BIN_FILE} — skipping"
        continue
    fi

    PASS_COUNT=$((PASS_COUNT + 1))
    BIN_OUT="${CHROM_ANNOT}/${BIN_NAME}"
    mkdir -p "${BIN_OUT}"

    echo ""
    echo "  Annotating: ${BIN_NAME}"
    echo "    Completeness=${COMPLETENESS}%  Contamination=${CONTAMINATION}%"
    echo "    Output: ${BIN_OUT}"

    # ── Run Bakta annotation ───────────────────────────────────────────────────
    # --db              : path to Bakta reference database
    # --output          : output directory for this bin
    # --prefix          : prefix for all output files (uses bin name)
    # --threads         : CPU threads
    # --meta            : metagenome mode — relaxes gene prediction thresholds
    #                     suitable for potentially incomplete MAG sequences
    # --keep-contig-headers : preserve original contig names from the bin FASTA.
    #                     Without this, Bakta renames contigs sequentially
    #                     (contig_1, contig_2 …), which breaks genomic context
    #                     maps — GFF3 seqid will not match AMRFinder contig
    #                     names and no flanking genes appear in plots.
    #                     Bug identified and fixed 12-Jun-2026.
    # --skip-plot       : skip genome visualisation plot (saves time)
    # --force           : overwrite previous output if re-running this bin
    conda run -p "${ENV_PATH}" \
        bakta \
            --db                  "${BAKTA_DB}" \
            --output              "${BIN_OUT}" \
            --prefix              "${BIN_NAME}" \
            --threads             16 \
            --meta \
            --keep-contig-headers \
            --skip-plot \
            --force \
            "${BIN_FILE}"

    BIN_EXIT=$?

    if [[ ${BIN_EXIT} -eq 0 ]]; then
        # Count annotated features from the TSV output
        FEATURE_COUNT=$(tail -n +2 "${BIN_OUT}/${BIN_NAME}.tsv" 2>/dev/null | wc -l | tr -d ' ')
        echo "    ✅ ${BIN_NAME} annotated — features: ${FEATURE_COUNT}"
        ANNOT_SUCCESS=$((ANNOT_SUCCESS + 1))
    else
        echo "    ❌ ${BIN_NAME} FAILED — check log for details"
        ANNOT_FAIL=$((ANNOT_FAIL + 1))
    fi

done < "${QUALITY_REPORT}"

# ── Final summary ──────────────────────────────────────────────────────────────
echo ""
echo "========================================"
echo "  ✅ ${SAMPLE_NAME} complete: $(date)"
echo "  Bins below quality threshold (skipped): ${FAIL_COUNT}"
echo "  Bins passed quality filter:              ${PASS_COUNT}"
echo "  Annotations successful:                  ${ANNOT_SUCCESS}"
echo "  Annotations failed:                      ${ANNOT_FAIL}"
echo ""
echo "  Output directory: ${CHROM_ANNOT}"
echo "  Key output files per bin:"
echo "    BinName.gff3  — genome annotation (use for Phase 5 co-localisation)"
echo "    BinName.faa   — protein sequences"
echo "    BinName.tsv   — tab-delimited feature table"
echo "========================================"

# Exit with failure if any annotation failed
if [[ ${ANNOT_FAIL} -gt 0 ]]; then
    echo "  WARNING: ${ANNOT_FAIL} bin(s) failed annotation — review logs above"
    exit 1
fi
