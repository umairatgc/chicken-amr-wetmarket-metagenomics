#!/bin/bash
BASE_DIR="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
SAMPLE_LIST="${BASE_DIR}/r_analysis/decontam/sample_list_real.txt"
BINS_DIR="${BASE_DIR}/bins"
LOG_DIR="${BASE_DIR}/logs"

PASS=0; FAIL=0

while IFS= read -r SAMPLE; do
    DEPTH="${BINS_DIR}/${SAMPLE}/coverage/${SAMPLE}_depth.txt"
    BAM="${BINS_DIR}/${SAMPLE}/coverage/${SAMPLE}_sorted.bam"
    BAI="${BINS_DIR}/${SAMPLE}/coverage/${SAMPLE}_sorted.bam.bai"

    ISSUES=()
    [[ ! -s "${DEPTH}" ]] && ISSUES+=("depth.txt missing/empty")
    [[ ! -s "${BAM}"   ]] && ISSUES+=("BAM missing/empty")
    [[ ! -f "${BAI}"   ]] && ISSUES+=("BAM index missing")

    if [[ ${#ISSUES[@]} -eq 0 ]]; then
        CONTIG_COUNT=$(tail -n +2 "${DEPTH}" | wc -l | tr -d ' ')
        echo "  ✅  ${SAMPLE}  (${CONTIG_COUNT} contigs)"
        PASS=$((PASS + 1))
    else
        echo "  ❌  ${SAMPLE}  — ${ISSUES[*]}"
        FAIL=$((FAIL + 1))
    fi
done < "${SAMPLE_LIST}"

echo ""
echo "========================================"
echo "  PASS: ${PASS} / 61"
echo "  FAIL: ${FAIL} / 61"
echo "========================================"
