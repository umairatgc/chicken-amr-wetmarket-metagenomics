#!/bin/bash
# =============================================================================
# run_bakta_plasmid_annot.sh
# SLURM array job — annotate all MOB-suite plasmid FASTAs with Bakta
#
# One array task per sample (61 samples total).
# For each sample, loops over all plasmid_*.fasta files produced by MOB-suite
# and runs Bakta on each one individually.
#
# Output structure (mirrors existing annotation/ layout):
#   annotation/<SAMPLE>/bakta_plasmid/<plasmid_ID>/
#       <plasmid_ID>.gff3   ← seqids match polished assembly contig names
#       <plasmid_ID>.faa
#       <plasmid_ID>.fna
#       ...
#
# Logs: logs/bakta_plasmid_<array_ID>.log / .err
#
# Submit: sbatch scripts/run_bakta_plasmid_annot.sh
# =============================================================================

#SBATCH --job-name=bakta_plasmid
#SBATCH --array=1-61
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/bakta_plasmid_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/bakta_plasmid_%a.err

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
BAKTA_DB=$BASE/databases/bakta/db
BAKTA_ENV=$BASE/envs/bakta

# ── All 61 samples with MOB-suite plasmid output ──────────────────────────────
SAMPLES=(
  PK-U3-1101P2  PK-U3-1105P3  PK-U3-1211P40 PK-U3-1215P42 PK-U3-1221P45
  PK-U3-1235P51 PK-U3-12EF97  PK-U3-12ES96  PK-U3-12EW99  PK-U3-2101P9
  PK-U3-2121P14 PK-U3-2127P15 PK-U3-2235P62 PK-U3-3101P17 PK-U3-3103P18
  PK-U3-3127P23 PK-U3-3207P67 PK-U3-32EE122 PK-U3-32EE125 PK-U3-32EF109
  PK-U3-32EW113 PK-U3-4105P26 PK-U3-4203P78 PK-U3-4211P83 PK-U3-4215P84
  PK-U3-42EM130 PK-U3-42ER128 PK-U3-42EW142 PK-UF-1002    PK-UF-1009
  PK-UF-1010    PK-UF-1015    PK-UF-1016    PK-UF-1020    PK-UF-1022
  PK-UF-1025    PK-UF-1028    PK-UF-1034    PK-UF-1040    PK-UF-1044
  PK-UF-1055    PK-UF-1059    PK-UF-1064    PK-UF-1069    PK-UF-1070
  PK-UF-1075    PK-UF-1078    PK-UF-1101    PK-UF-1111    PK-UF-1131
  PK-UF-1139    PK-UF-1140    PK-UF-1143    PK-UF-1147    PK-UF-1154
  PK-UF-1155    PK-UF-1158    PK-UF-1160    PK-UF-1162    PK-UF-1180
  PK-UF-1181
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# ── Per-sample I/O directories ────────────────────────────────────────────────
MOBSUITE_DIR=$BASE/mobile_elements/$SAMPLE/mobsuite   # MOB-suite plasmid FASTAs
ANNOT_OUT=$BASE/annotation/$SAMPLE/bakta_plasmid       # output lives here

echo "========================================================"
echo "Sample     : $SAMPLE"
echo "Array task : $SLURM_ARRAY_TASK_ID"
echo "Start      : $(date)"
echo "Node       : $(hostname)"
echo "========================================================"

# Create output root for this sample
mkdir -p "$ANNOT_OUT"

# Count plasmid FASTAs available for this sample
N_PLASMIDS=$(ls "$MOBSUITE_DIR"/plasmid_*.fasta 2>/dev/null | wc -l)

if [ "$N_PLASMIDS" -eq 0 ]; then
    echo "ERROR: No plasmid FASTAs found in $MOBSUITE_DIR — exiting."
    exit 1
fi

echo "Plasmids to annotate: $N_PLASMIDS"
echo ""

# ── Main loop — one Bakta run per plasmid FASTA ───────────────────────────────
DONE=0
SKIPPED=0
FAILED=0

for fasta in "$MOBSUITE_DIR"/plasmid_*.fasta; do

    # Derive plasmid name from filename (e.g. plasmid_AA372)
    plasmid_name=$(basename "$fasta" .fasta)
    out_dir=$ANNOT_OUT/$plasmid_name

    # Skip if GFF3 already exists (resume-safe)
    if [ -f "$out_dir/${plasmid_name}.gff3" ]; then
        echo "[SKIP] $plasmid_name — already annotated"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    echo "[START] $plasmid_name  ($(date +%H:%M:%S))"

    # Run Bakta in metagenome mode:
    #   --meta              : metagenome mode (no genus-specific tuning, faster)
    #   --keep-contig-headers: preserve original contig IDs (contig_XXXX) in GFF3
    #   --skip-plot         : skip circular genome plot (saves time)
    #   --threads           : use all allocated CPUs
    conda run -p "$BAKTA_ENV" bakta \
        --db "$BAKTA_DB" \
        --output "$out_dir" \
        --prefix "$plasmid_name" \
        --threads "$SLURM_CPUS_PER_TASK" \
        --meta \
        --keep-contig-headers \
        --skip-plot \
        "$fasta"

    EXIT_CODE=$?

    if [ $EXIT_CODE -eq 0 ]; then
        echo "[DONE]  $plasmid_name  ($(date +%H:%M:%S))"
        DONE=$((DONE + 1))
    else
        echo "[FAIL]  $plasmid_name — exit code $EXIT_CODE"
        FAILED=$((FAILED + 1))
        # Remove incomplete output so it can be retried
        rm -rf "$out_dir"
    fi

done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "========================================================"
echo "Finished: $SAMPLE  $(date)"
echo "  Annotated : $DONE"
echo "  Skipped   : $SKIPPED  (already done)"
echo "  Failed    : $FAILED"
echo "========================================================"
