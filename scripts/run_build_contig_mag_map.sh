#!/bin/bash
#SBATCH --job-name=contig_mag_map
#SBATCH --clusters=arc
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/contig_mag_map.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/contig_mag_map.err

# =============================================================================
# run_build_contig_mag_map.sh
# Phase 5 — Step 0: Build contig → MAG mapping table
#
# For every quality MAG bin (passing CheckM2 ≥50% / <10%), extract all contig
# names and join with:
#   - CheckM2 completeness + contamination
#   - GTDB-Tk taxonomy (bacterial + archaeal)
#
# OUTPUT:
#   phase5/contig_mag_map.tsv
#   Columns: sample | contig | bin | completeness | contamination | taxonomy
#
# This file is the backbone of Phase 5 — every ARG and MGE result is joined
# against this table to assign organism identity to each contig.
#
# SUBMISSION:
#   sbatch scripts/run_build_contig_mag_map.sh
#   Or run interactively: bash scripts/run_build_contig_mag_map.sh
# =============================================================================

BASE="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
SAMPLE_LIST="${BASE}/r_analysis/decontam/sample_list_real.txt"
OUTFILE="${BASE}/phase5/contig_mag_map.tsv"

MIN_COMPLETENESS=50
MAX_CONTAMINATION=10

mkdir -p "${BASE}/phase5"
mkdir -p "${BASE}/logs"

echo "========================================"
echo "  Building contig → MAG mapping table"
echo "  Started: $(date)"
echo "========================================"

# Write header
echo -e "sample\tcontig\tbin\tcompleteness\tcontamination\ttaxonomy" > "${OUTFILE}"

TOTAL_CONTIGS=0
TOTAL_BINS=0
SAMPLES_PROCESSED=0

while IFS= read -r SAMPLE; do
    CHECKM2="${BASE}/bins/${SAMPLE}/checkm2/quality_report.tsv"
    GTDBTK_BAC="${BASE}/taxonomy_mags/${SAMPLE}/classify/gtdbtk.bac120.summary.tsv"
    GTDBTK_ARC="${BASE}/taxonomy_mags/${SAMPLE}/classify/gtdbtk.ar53.summary.tsv"

    if [[ ! -f "${CHECKM2}" ]]; then
        echo "  WARNING: No CheckM2 report for ${SAMPLE} — skipping"
        continue
    fi

    SAMPLE_CONTIGS=0
    SAMPLE_BINS=0

    while IFS=$'\t' read -r BIN_NAME COMP CONT REST; do
        # Skip header
        [[ "${BIN_NAME}" == "Name" ]] && continue

        # Quality filter
        PASSES=$(echo "${COMP} ${CONT}" | awk \
            -v minc="${MIN_COMPLETENESS}" -v maxcon="${MAX_CONTAMINATION}" \
            '{if ($1+0 >= minc && $2+0 < maxcon) print "yes"; else print "no"}')
        [[ "${PASSES}" != "yes" ]] && continue

        BIN_FASTA="${BASE}/bins/${SAMPLE}/metabat2/${BIN_NAME}.fa"
        [[ ! -f "${BIN_FASTA}" ]] && continue

        # Look up taxonomy — check bac120 first, then ar53
        TAXONOMY="unclassified"
        if [[ -f "${GTDBTK_BAC}" ]]; then
            TAX=$(grep "^${BIN_NAME}"$'\t' "${GTDBTK_BAC}" 2>/dev/null | cut -f2 | head -1)
            [[ -n "${TAX}" && "${TAX}" != "N/A" ]] && TAXONOMY="${TAX}"
        fi
        if [[ "${TAXONOMY}" == "unclassified" && -f "${GTDBTK_ARC}" ]]; then
            TAX=$(grep "^${BIN_NAME}"$'\t' "${GTDBTK_ARC}" 2>/dev/null | cut -f2 | head -1)
            [[ -n "${TAX}" && "${TAX}" != "N/A" ]] && TAXONOMY="${TAX}"
        fi

        # Extract contig names from bin FASTA and write one row per contig
        while IFS= read -r CONTIG; do
            echo -e "${SAMPLE}\t${CONTIG}\t${BIN_NAME}\t${COMP}\t${CONT}\t${TAXONOMY}"
            SAMPLE_CONTIGS=$((SAMPLE_CONTIGS + 1))
        done < <(grep "^>" "${BIN_FASTA}" | sed 's/^>//' | awk '{print $1}')

        SAMPLE_BINS=$((SAMPLE_BINS + 1))

    done < "${CHECKM2}"

    TOTAL_CONTIGS=$((TOTAL_CONTIGS + SAMPLE_CONTIGS))
    TOTAL_BINS=$((TOTAL_BINS + SAMPLE_BINS))
    SAMPLES_PROCESSED=$((SAMPLES_PROCESSED + 1))
    echo "  ✅ ${SAMPLE}: ${SAMPLE_BINS} bins, ${SAMPLE_CONTIGS} contigs"

done < "${SAMPLE_LIST}" >> "${OUTFILE}"

echo ""
echo "========================================"
echo "  ✅ Complete: $(date)"
echo "  Samples processed:  ${SAMPLES_PROCESSED}"
echo "  Quality bins total: ${TOTAL_BINS}"
echo "  Contigs mapped:     ${TOTAL_CONTIGS}"
echo "  Output: ${OUTFILE}"
echo "========================================"
