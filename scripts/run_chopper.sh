#!/bin/bash
#SBATCH --job-name=chopper_filter
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-62
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/chopper_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/chopper_%A_%a.err

# =============================================================================
# run_chopper.sh
# Quality and length filtering of adapter-trimmed ONT reads using Chopper.
#
# Thresholds:
#   -q 15   →  minimum Phred quality score 15 (96.8% accuracy per base)
#   -l 1000 →  minimum read length 1000 bp
#
# Input:  trimmed/*.fastq.gz  (Porechop_ABI output)
# Output: filtered/*.fastq.gz
#
# IMPORTANT: Run AFTER run_porechop.sh has completed.
#            Verify all 63 trimmed files exist before submitting this job.
#
# FIX: eval "$(conda shell.bash hook)" is REQUIRED in SLURM batch scripts.
#      Without it, conda activate fails in non-interactive shells.
# =============================================================================

# ── Load conda and initialise shell functions ────────────────────────────────
module load Anaconda3
eval "$(conda shell.bash hook)"   # CRITICAL: enables conda activate in SLURM

# ── Activate the qc environment (contains Chopper) ──────────────────────────
conda activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/qc

# ── Paths ────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
TRIMMED=$BASE/trimmed        # input: adapter-trimmed files from Porechop
OUTDIR=$BASE/filtered        # output: quality and length filtered files

mkdir -p "$OUTDIR"

# ── Sample list (61 samples + 2 controls) ────────────────────────────────────
# Order must match exactly with run_porechop.sh
# SLURM_ARRAY_TASK_ID picks the sample by its position in this list
SAMPLES=(
    PK-UF-1002 PK-UF-1009 PK-UF-1015 PK-UF-1022 PK-UF-1025
    PK-UF-1028 PK-UF-1044 PK-UF-1064 PK-UF-1069 PK-UF-1070
    PK-UF-1101 PK-UF-1139 PK-UF-1140 PK-UF-1154 PK-UF-1155
    PK-UF-1162 PK-UF-1180 PK-U3-32EW113 PK-U3-32EE122 PK-U3-42ER128
    PK-U3-42EM130 PK-U3-42EW142 PK-U3-2121P14 PK-U3-2127P15
    PK-U3-1211P40 PK-U3-1235P51 PK-U3-4203P78 PK-U3-4215P84
    PK-U3-1101P2 PK-U3-4105P26 PK-UF-1010 PK-UF-1016 PK-UF-1020
    PK-UF-1034 PK-UF-1040 PK-UF-1055 PK-UF-1059 PK-UF-1075
    PK-UF-1078 PK-UF-1111 PK-UF-1131 PK-UF-1143 PK-UF-1147
    PK-UF-1158 PK-UF-1160 PK-UF-1181 PK-U3-12ES96 PK-U3-12EF97
    PK-U3-12EW99 PK-U3-32EF109 PK-U3-32EE125 PK-U3-1105P3
    PK-U3-2101P9 PK-U3-3101P17 PK-U3-3103P18 PK-U3-3127P23
    PK-U3-1215P42 PK-U3-1221P45 PK-U3-2235P62 PK-U3-3207P67 PK-U3-4211P83
    EB-neg                          # negative extraction blank control
    AE-buffer                       # negative AE buffer control
)

# ── Pick this job's sample using the array task index ────────────────────────
# SLURM_ARRAY_TASK_ID is automatically set by SLURM for each job in the array
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

INPUT="$TRIMMED/${SAMPLE}_trimmed.fastq.gz"
OUTPUT="$OUTDIR/${SAMPLE}_filtered.fastq.gz"

echo "======================================================"
echo "Array task ID : $SLURM_ARRAY_TASK_ID"
echo "Sample        : $SAMPLE"
echo "Input         : $INPUT"
echo "Output        : $OUTPUT"
echo "Start time    : $(date)"
echo "======================================================"

# ── Check input exists ───────────────────────────────────────────────────────
if [ ! -f "$INPUT" ]; then
    echo "WARNING: $INPUT not found — skipping $SAMPLE"
    echo "         (Check that Porechop completed successfully for this sample)"
    exit 0
fi

# ── Run Chopper ──────────────────────────────────────────────────────────────
# Chopper reads from stdin (pipe), so we decompress with gunzip -c
# then pipe through chopper, then recompress with gzip
#
# -q 15     : discard reads with mean Phred quality below Q15
# -l 1000   : discard reads shorter than 1000 bp
# --threads : number of CPU threads to use

gunzip -c "$INPUT" | \
    chopper \
        -q 15 \
        -l 1000 \
        --threads 8 | \
    gzip > "$OUTPUT"

# ── Verify output was created ────────────────────────────────────────────────
if [ -f "$OUTPUT" ]; then
    # Count reads in filtered output (each read = 4 lines in FASTQ)
    READ_COUNT=$(gunzip -c "$OUTPUT" | wc -l)
    READ_COUNT=$((READ_COUNT / 4))
    echo ""
    echo "Done: $SAMPLE"
    echo "  Filtered reads: $READ_COUNT"
    echo "  Output file   : $OUTPUT"
    echo "  File size     : $(du -sh "$OUTPUT" | cut -f1)"
else
    echo "ERROR: Output file not created for $SAMPLE"
    exit 1
fi

echo "End time: $(date)"
echo "======================================================"

conda deactivate
