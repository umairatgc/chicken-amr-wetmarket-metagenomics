#!/bin/bash
#SBATCH --job-name=abricate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/abricate_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/abricate_%A_%a.err

# =============================================================================
# run_abricate.sh
# Step 7b-iii -- Per-sample assembly-based AMR/virulence screening with Abricate 1.0.4
#
# WHAT THIS SCRIPT DOES:
#   Screens Medaka-polished assemblies against 7 databases covering AMR genes,
#   virulence factors, and plasmid replicons. Each database is run separately
#   and results are saved per database. A combined summary is also generated.
#
# DATABASES SCREENED:
#   ncbi        -- NCBI AMRFinderPlus gene database
#   card        -- CARD (Comprehensive Antibiotic Resistance Database)
#   resfinder   -- ResFinder acquired resistance genes
#   vfdb        -- Virulence Factor Database
#   plasmidfinder -- Plasmid replicon typing
#   ecoh        -- E. coli O/H serotyping antigens
#   argannot    -- ARG-ANNOT resistance gene database
#
# PARAMETERS:
#   --minid 80  -- minimum nucleotide identity 80%
#   --mincov 60 -- minimum gene coverage 60%
#
# INPUT PER SAMPLE:
#   polished/SampleName/consensus.fasta     -- Medaka-polished assembly
#
# OUTPUT PER SAMPLE:
#   amr_results/abricate/SampleName/abricate_ncbi.txt
#   amr_results/abricate/SampleName/abricate_card.txt
#   amr_results/abricate/SampleName/abricate_resfinder.txt
#   amr_results/abricate/SampleName/abricate_vfdb.txt
#   amr_results/abricate/SampleName/abricate_plasmidfinder.txt
#   amr_results/abricate/SampleName/abricate_ecoh.txt
#   amr_results/abricate/SampleName/abricate_argannot.txt
#   amr_results/abricate/SampleName/abricate_summary.txt  -- combined
#
# SAFE RESTART:
#   If abricate_summary.txt already exists and is non-empty, the task skips.
#   Delete the sample folder to force re-run.
#
# RUN AFTER:  run_medaka.sh   (polished assemblies must exist)
# RUN BEFORE: downstream AMR summary / R analysis
#
# SUBMISSION:
#   sbatch scripts/run_abricate.sh
#   Monitor:  squeue -u $USER
#   Progress: ls amr_results/abricate/*/abricate_summary.txt | wc -l
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
SAMPLE_LIST=$BASE/r_analysis/decontam/sample_list_real.txt
POLISHED=$BASE/polished
OUT_BASE=$BASE/amr_results/abricate
CONDA_ENV=$BASE/envs/abricate

# ── Abricate parameters ───────────────────────────────────────────────────────
MIN_ID=80
MIN_COV=60
DATABASES="ncbi card resfinder vfdb plasmidfinder ecoh argannot"

# ── Validate conda environment ────────────────────────────────────────────────
if [ ! -d "$CONDA_ENV" ]; then
    echo "ERROR: Abricate environment not found: $CONDA_ENV"
    exit 1
fi

conda activate "$CONDA_ENV"

# ── Validate Abricate is available ───────────────────────────────────────────
if ! command -v abricate &> /dev/null; then
    echo "ERROR: abricate not found in environment: $CONDA_ENV"
    exit 1
fi

ABRICATE_VERSION=$(abricate --version 2>&1 | head -1)

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
SUMMARY_FILE="$OUT_DIR/abricate_summary.txt"

echo "======================================================"
echo "Abricate assembly-based AMR/virulence screening"
echo "Array task:   $SLURM_ARRAY_TASK_ID / 61"
echo "Sample:       $SAMPLE"
echo "Version:      $ABRICATE_VERSION"
echo "Databases:    $DATABASES"
echo "Min identity: ${MIN_ID}%"
echo "Min coverage: ${MIN_COV}%"
echo "Start time:   $(date)"
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

# ── Skip if already complete ──────────────────────────────────────────────────
if [ -f "$SUMMARY_FILE" ] && [ -s "$SUMMARY_FILE" ]; then
    N=$(tail -n +2 "$SUMMARY_FILE" | wc -l)
    echo ""
    echo "Abricate results already exist: $SUMMARY_FILE"
    echo "Skipping $SAMPLE -- delete $OUT_DIR to re-run"
    exit 0
fi

mkdir -p "$OUT_DIR"

# ── Run Abricate against each database ───────────────────────────────────────
ABRICATE_ERRORS=0

for DB in $DATABASES; do
    OUT_FILE="$OUT_DIR/abricate_${DB}.txt"
    echo ""
    echo "Screening against database: $DB"

    abricate \
        --db "$DB" \
        --minid "$MIN_ID" \
        --mincov "$MIN_COV" \
        --threads "$SLURM_CPUS_PER_TASK" \
        "$IN_FASTA" > "$OUT_FILE"

    EXIT_CODE=$?

    if [ $EXIT_CODE -ne 0 ] || [ ! -s "$OUT_FILE" ]; then
        echo "  WARNING: Abricate failed or returned no output for database: $DB"
        ABRICATE_ERRORS=$((ABRICATE_ERRORS + 1))
    else
        N_HITS=$(tail -n +2 "$OUT_FILE" | grep -v "^$" | wc -l)
        echo "  Hits: $N_HITS"
    fi
done

# ── Generate combined summary ─────────────────────────────────────────────────
echo ""
echo "Generating combined summary..."

abricate --summary \
    "$OUT_DIR"/abricate_*.txt > "$SUMMARY_FILE"

if [ $? -ne 0 ]; then
    echo "WARNING: Summary generation failed -- individual results still available"
fi

# ── Final summary ─────────────────────────────────────────────────────────────
echo ""
echo "======================================================"
echo "Abricate complete: $SAMPLE"
echo "End time:  $(date)"
echo ""
for DB in $DATABASES; do
    OUT_FILE="$OUT_DIR/abricate_${DB}.txt"
    if [ -f "$OUT_FILE" ]; then
        N=$(tail -n +2 "$OUT_FILE" | grep -v "^$" | wc -l)
        printf "  %-15s %s hits\n" "$DB:" "$N"
    fi
done
echo ""
if [ $ABRICATE_ERRORS -gt 0 ]; then
    echo "  WARNING: $ABRICATE_ERRORS database(s) had errors -- check logs"
fi
echo "  Summary:  $SUMMARY_FILE"
echo "======================================================"

conda deactivate
