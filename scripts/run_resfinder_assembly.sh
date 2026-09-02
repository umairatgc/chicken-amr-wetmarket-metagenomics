#!/bin/bash
#SBATCH --job-name=resfinder_asm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/resfinder_asm_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/resfinder_asm_%A_%a.err

# =============================================================================
# run_resfinder_assembly.sh
# Step 7b-iv -- Per-sample assembly-based AMR detection with ResFinder 4.7.2
#
# WHAT THIS SCRIPT DOES:
#   Screens Medaka-polished assemblies (consensus.fasta) against the CGE
#   ResFinder database. Complements the read-based ResFinder analysis by
#   working on assembled contigs, providing full gene context and enabling
#   detection of complete resistance genes with flanking sequences.
#
# KEY DIFFERENCE FROM READ-BASED:
#   Read-based (run_resfinder.sh):    -ifq reads.fastq.gz   (FASTQ input)
#   Assembly-based (this script):     -ifa assembly.fasta   (FASTA input)
#
# PARAMETERS:
#   --acquired  : detect acquired resistance genes
#   -l 0.60     : minimum gene coverage 60%
#   -t 0.90     : minimum nucleotide identity 90%
#
# INPUT PER SAMPLE:
#   polished/SampleName/consensus.fasta      -- Medaka-polished assembly
#
# OUTPUT PER SAMPLE:
#   amr_results/resfinder_assembly/SampleName/ResFinder_results_tab.txt
#   amr_results/resfinder_assembly/SampleName/ResFinder_results.txt
#   amr_results/resfinder_assembly/SampleName/pheno_table.txt
#
# SAFE RESTART:
#   If ResFinder_results_tab.txt already exists and is non-empty, the task skips.
#   Delete the sample folder to force re-run.
#
# DATABASES USED:
#   ResFinder DB : databases/resfinder/resfinder_db/  (CGE, updated 2026-01-26)
#
# RUN AFTER:  run_medaka.sh   (polished assemblies must exist)
# RUN BEFORE: downstream AMR summary / R analysis
#
# SUBMISSION:
#   sbatch scripts/run_resfinder_assembly.sh
#   Monitor:  squeue -u $USER
#   Progress: ls amr_results/resfinder_assembly/*/ResFinder_results_tab.txt | wc -l
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
SAMPLE_LIST=$BASE/r_analysis/decontam/sample_list_real.txt
POLISHED=$BASE/polished
DB_RES=$BASE/databases/resfinder/resfinder_db
OUT_BASE=$BASE/amr_results/resfinder_assembly
CONDA_ENV=$BASE/envs/resfinder

# ── Validate conda environment ────────────────────────────────────────────────
if [ ! -d "$CONDA_ENV" ]; then
    echo "ERROR: ResFinder environment not found: $CONDA_ENV"
    echo "       Run: sbatch scripts/install_resfinder.sh"
    exit 1
fi

conda activate "$CONDA_ENV"

# ── Validate ResFinder is available ──────────────────────────────────────────
if ! python -m resfinder --help &> /dev/null; then
    echo "ERROR: ResFinder not found in conda environment: $CONDA_ENV"
    exit 1
fi

RESFINDER_VERSION=$(python -m resfinder --version 2>&1 | head -1)

# ── Validate database exists ──────────────────────────────────────────────────
if [ ! -d "$DB_RES" ]; then
    echo "ERROR: ResFinder database not found: $DB_RES"
    exit 1
fi

# ── Check KMA indices exist ───────────────────────────────────────────────────
INDEX_CHECK="$DB_RES/aminoglycoside.length.b"
if [ ! -f "$INDEX_CHECK" ]; then
    echo "ERROR: KMA index files not found in $DB_RES"
    echo "       Run: sbatch scripts/build_kma_index.sh"
    exit 1
fi

# ── Validate sample list ──────────────────────────────────────────────────────
if [ ! -f "$SAMPLE_LIST" ]; then
    echo "ERROR: Sample list not found: $SAMPLE_LIST"
    exit 1
fi

mkdir -p "$OUT_BASE"

# ── Get this task's sample ────────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

if [ -z "$SAMPLE" ]; then
    echo "ERROR: No sample found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# ── Input / output paths ──────────────────────────────────────────────────────
IN_FASTA="$POLISHED/$SAMPLE/consensus.fasta"
OUT_DIR="$OUT_BASE/$SAMPLE"
RESULT_FILE="$OUT_DIR/ResFinder_results_tab.txt"

echo "======================================================"
echo "ResFinder assembly-based AMR detection"
echo "Array task:          $SLURM_ARRAY_TASK_ID / 61"
echo "Sample:              $SAMPLE"
echo "ResFinder version:   $RESFINDER_VERSION"
echo "Identity threshold:  90%"
echo "Coverage threshold:  60%"
echo "Input type:          FASTA assembly (-ifa)"
echo "CPUs:                $SLURM_CPUS_PER_TASK"
echo "Start time:          $(date)"
echo "======================================================"

# ── Validate input ────────────────────────────────────────────────────────────
if [ ! -f "$IN_FASTA" ]; then
    echo "ERROR: Polished assembly not found: $IN_FASTA"
    echo "       Run run_medaka.sh first"
    exit 1
fi

if [ ! -s "$IN_FASTA" ]; then
    echo "ERROR: Assembly file is empty: $IN_FASTA"
    exit 1
fi

N_CONTIGS=$(grep -c "^>" "$IN_FASTA" || true)
echo "Input contigs: $N_CONTIGS"
echo "Assembly size: $(du -h "$IN_FASTA" | cut -f1)"

# ── Skip if already complete ──────────────────────────────────────────────────
if [ -f "$RESULT_FILE" ] && [ -s "$RESULT_FILE" ]; then
    N=$(tail -n +2 "$RESULT_FILE" | wc -l)
    echo ""
    echo "ResFinder results already exist with $N hits: $RESULT_FILE"
    echo "Skipping $SAMPLE -- delete $OUT_DIR to re-run"
    exit 0
fi

mkdir -p "$OUT_DIR"

# ── Run ResFinder ─────────────────────────────────────────────────────────────
# -ifa : FASTA input (assembly-based) -- different from -ifq (FASTQ/read-based)
echo ""
echo "Running ResFinder (assembly mode)..."
echo "  Input:    $IN_FASTA"
echo "  Database: $DB_RES"
echo "  Output:   $OUT_DIR"
echo ""

python -m resfinder \
    -ifa "$IN_FASTA" \
    --acquired \
    -db_res "$DB_RES" \
    -l 0.60 \
    -t 0.90 \
    -o "$OUT_DIR"

RESFINDER_EXIT=$?

if [ $RESFINDER_EXIT -ne 0 ]; then
    echo ""
    echo "ERROR: ResFinder failed for $SAMPLE (exit code $RESFINDER_EXIT)"
    exit 1
fi

# ── Validate output ───────────────────────────────────────────────────────────
if [ ! -f "$RESULT_FILE" ]; then
    echo "ERROR: ResFinder_results_tab.txt not found after run"
    exit 1
fi

# ── Final summary ─────────────────────────────────────────────────────────────
N_HITS=$(tail -n +2 "$RESULT_FILE" | grep -v "^$" | wc -l)

echo ""
echo "======================================================"
echo "ResFinder assembly complete: $SAMPLE"
echo "End time:  $(date)"
echo ""
echo "  AMR genes detected: $N_HITS"
echo "  Results:  $RESULT_FILE"
echo ""
if [ "$N_HITS" -gt 0 ]; then
    echo "  Top hits:"
    tail -n +2 "$RESULT_FILE" | awk -F'\t' '{print "    " $1 " | " $2 " | Identity: " $4}' | head -10
fi
echo "======================================================"

conda deactivate
