#!/bin/bash
#SBATCH --job-name=flye_assembly
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --array=1-61
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/assembly_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/assembly_%A_%a.err

# =============================================================================
# run_assembly.sh
# Step 6 -- Per-sample metagenome assembly with Flye
#
# WHAT THIS SCRIPT DOES:
#   Assembles decontaminated ONT reads (clean_reads/) into contigs using
#   Flye in metagenome mode. One assembly per sample, 61 jobs in parallel.
#
# FLYE MODE CONFIRMED:
#   basecall_model_version_id = dna_r10.4.1_e8.2_400bps_sup@v5.2.0
#   R10.4.1 flow cell + Dorado SUP = --nano-hq is correct.
#   DO NOT change to --nano-raw unless you re-sequence with older chemistry.
#
# INPUT PER SAMPLE:
#   clean_reads/SampleName_clean.fastq.gz   -- decontaminated FASTQ reads
#
# OUTPUT PER SAMPLE:
#   assemblies/SampleName/assembly.fasta      -- final contigs
#   assemblies/SampleName/assembly_info.txt   -- per-contig stats from Flye
#   assemblies/SampleName/assembly_stats.txt  -- summary (N50, contigs, size)
#   assemblies/SampleName/flye.log            -- full Flye log
#
# SAFE RESTART:
#   If assembly.fasta already exists and is non-empty, the task skips.
#   Delete the sample folder to force re-assembly.
#
# RUN AFTER:  run_extract_kraken.sh  (confirmed FASTQ format, @ headers)
# RUN BEFORE: run_amr_assembly.sh    (contig-based AMR with AMRFinderPlus)
#
# SUBMISSION:
#   sbatch run_assembly.sh
#   Monitor:  squeue -u $USER
#   Progress: tail -f logs/assembly_JOBID_TASKID.log
#   Count done: ls assemblies/*/assembly.fasta 2>/dev/null | wc -l
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
CLEAN=$BASE/clean_reads
ASSEMBLIES=$BASE/assemblies
SAMPLE_LIST=$BASE/r_analysis/decontam/sample_list_real.txt
CONDA_ENV=$BASE/envs/assembly

# ── Validate conda environment exists ─────────────────────────────────────────
if [ ! -d "$CONDA_ENV" ]; then
    echo "ERROR: Conda environment not found: $CONDA_ENV"
    exit 1
fi

conda activate "$CONDA_ENV"

# ── Validate Flye is available ────────────────────────────────────────────────
if ! command -v flye &> /dev/null; then
    echo "ERROR: flye not found in conda environment: $CONDA_ENV"
    exit 1
fi

FLYE_VERSION=$(flye --version 2>&1)

# ── Validate sample list exists ───────────────────────────────────────────────
if [ ! -f "$SAMPLE_LIST" ]; then
    echo "ERROR: Sample list not found: $SAMPLE_LIST"
    echo "       Run run_extract_kraken.sh first to generate it"
    exit 1
fi

mkdir -p "$ASSEMBLIES"

# ── Get this task's sample name ───────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

if [ -z "$SAMPLE" ]; then
    echo "ERROR: No sample found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# ── Input/output paths ────────────────────────────────────────────────────────
IN_FASTQ="$CLEAN/${SAMPLE}_clean.fastq.gz"
OUT_DIR="$ASSEMBLIES/$SAMPLE"
FINAL_CONTIGS="$OUT_DIR/assembly.fasta"
STATS_FILE="$OUT_DIR/assembly_stats.txt"

echo "======================================================"
echo "Flye metagenome assembly"
echo "Array task:   $SLURM_ARRAY_TASK_ID / 61"
echo "Sample:       $SAMPLE"
echo "Flye version: $FLYE_VERSION"
echo "Read mode:    --nano-hq  (R10.4.1 Dorado SUP confirmed)"
echo "CPUs:         $SLURM_CPUS_PER_TASK"
echo "Memory:       64G"
echo "Start time:   $(date)"
echo "======================================================"

# ── Validate input FASTQ ──────────────────────────────────────────────────────
if [ ! -f "$IN_FASTQ" ]; then
    echo "ERROR: Clean FASTQ not found: $IN_FASTQ"
    echo "       Run run_extract_kraken.sh first"
    exit 1
fi

# Check file is non-empty
if [ ! -s "$IN_FASTQ" ]; then
    echo "ERROR: Clean FASTQ is empty: $IN_FASTQ"
    exit 1
fi

# Check FASTQ format -- first character must be @ not >
FIRST_CHAR=$(zcat "$IN_FASTQ" | head -1 | cut -c1)
if [ "$FIRST_CHAR" != "@" ]; then
    echo "ERROR: Input file is not FASTQ format (first char = '$FIRST_CHAR', expected '@')"
    echo "       File may be FASTA -- re-run run_extract_kraken.sh with --fastq-output"
    exit 1
fi

echo "Input validated: $IN_FASTQ"
echo "Format check:    FASTQ confirmed (@)"
echo "File size:       $(du -h "$IN_FASTQ" | cut -f1)"

# ── Skip if assembly already complete ─────────────────────────────────────────
if [ -f "$FINAL_CONTIGS" ] && [ -s "$FINAL_CONTIGS" ]; then
    N=$(grep -c "^>" "$FINAL_CONTIGS")
    echo ""
    echo "Assembly already exists with $N contigs: $FINAL_CONTIGS"
    echo "Skipping $SAMPLE -- delete $OUT_DIR to re-run"
    exit 0
fi

# ── Create output directory ───────────────────────────────────────────────────
mkdir -p "$OUT_DIR"

# ── Run Flye ──────────────────────────────────────────────────────────────────
echo ""
echo "Running Flye..."
echo "  Input:  $IN_FASTQ"
echo "  Output: $OUT_DIR"
echo ""

flye \
    --nano-hq "$IN_FASTQ" \
    --meta \
    --out-dir "$OUT_DIR" \
    --threads "$SLURM_CPUS_PER_TASK"

FLYE_EXIT=$?

if [ $FLYE_EXIT -ne 0 ]; then
    echo ""
    echo "ERROR: Flye failed for $SAMPLE (exit code $FLYE_EXIT)"
    echo "       Check Flye log: $OUT_DIR/flye.log"
    exit 1
fi

# ── Validate output ───────────────────────────────────────────────────────────
if [ ! -f "$FINAL_CONTIGS" ]; then
    echo "ERROR: assembly.fasta not found after Flye completed"
    echo "       Flye produced no contigs -- sample may have too few reads"
    echo "       Check: $OUT_DIR/flye.log"
    exit 1
fi

if [ ! -s "$FINAL_CONTIGS" ]; then
    echo "ERROR: assembly.fasta is empty"
    echo "       Check: $OUT_DIR/flye.log"
    exit 1
fi

# ── Calculate assembly statistics ─────────────────────────────────────────────
N_CONTIGS=$(grep -c "^>" "$FINAL_CONTIGS")
TOTAL_BP=$(grep -v "^>" "$FINAL_CONTIGS" | tr -d '\n' | wc -c)
TOTAL_MB=$(awk "BEGIN {printf \"%.2f\", $TOTAL_BP/1000000}")

# N50 calculation
N50=$(awk '/^>/{if(seq) print length(seq); seq=""; next} {seq=seq$0}
           END{if(seq) print length(seq)}' "$FINAL_CONTIGS" \
    | sort -rn \
    | awk -v total="$TOTAL_BP" '
        BEGIN{cumsum=0; n50=0}
        { cumsum+=$1; if(cumsum >= total/2 && n50==0) n50=$1 }
        END{print n50}')

LONGEST=$(awk '/^>/{if(seq) print length(seq); seq=""; next} {seq=seq$0}
               END{if(seq) print length(seq)}' "$FINAL_CONTIGS" \
         | sort -rn | head -1)

# ── Write stats file ──────────────────────────────────────────────────────────
cat > "$STATS_FILE" << STATS
Sample:        $SAMPLE
Assembly date: $(date)
Flye version:  $FLYE_VERSION
Read mode:     --nano-hq --meta

Contigs:       $N_CONTIGS
Total size:    ${TOTAL_MB} Mb (${TOTAL_BP} bp)
N50:           $N50 bp
Longest:       $LONGEST bp

Input:  $IN_FASTQ
Output: $FINAL_CONTIGS
STATS

# ── Final summary ─────────────────────────────────────────────────────────────
echo ""
echo "======================================================"
echo "Assembly complete: $SAMPLE"
echo "End time: $(date)"
echo ""
echo "  Contigs:    $N_CONTIGS"
echo "  Total size: ${TOTAL_MB} Mb"
echo "  N50:        $N50 bp"
echo "  Longest:    $LONGEST bp"
echo ""
echo "  Output: $FINAL_CONTIGS"
echo "  Stats:  $STATS_FILE"
echo "======================================================"

conda deactivate
