#!/bin/bash
#SBATCH --job-name=porechop_trim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-62                    # 61 samples + 2 controls = 63 total, indexed 0 to 62
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/porechop_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/porechop_%A_%a.err

# =============================================================================
# run_porechop.sh
# Step 2.3 — Adapter Trimming using Porechop_ABI v0.5.1
# SLURM Array Job version
#
# How array jobs work:
#   SLURM launches 63 separate jobs simultaneously — one per sample/control.
#   Each job gets a unique index ($SLURM_ARRAY_TASK_ID) from 0 to 62, which
#   is used to pick the corresponding sample from the SAMPLES array below.
#
#   %A in the log filename = the main array job ID
#   %a in the log filename = the individual task index (0, 1, 2 ... 62)
#   So each sample gets its own separate log and err file.
#
# Samples: 61 patient samples + 2 controls
#   EB-neg    — negative extraction blank control (should have very few reads)
#   AE-buffer — negative buffer control (should have very few reads)
#   Controls are processed through the same pipeline as samples to detect
#   any contamination that occurred during DNA extraction or library prep.
#
# What it does:
#   Removes adapter sequences ligated to DNA fragments during ONT library
#   preparation. --discard_middle removes chimeric reads where two fragments
#   accidentally joined together during sequencing.
#
# Input:  basecalled/SAMPLE.fastq.gz
# Output: trimmed/SAMPLE_trimmed.fastq.gz
#
# Run this BEFORE run_chopper.sh
# =============================================================================

# ── Load conda ────────────────────────────────────────────────────────────────
module load Anaconda3
eval "$(conda shell.bash hook)" 
conda activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/qc

# ── Paths ────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
READS=$BASE/basecalled       # input: raw basecalled fastq files
OUTDIR=$BASE/trimmed         # output: adapter-trimmed fastq files

mkdir -p "$OUTDIR"

# ── Sample list (61 samples + 2 controls) ────────────────────────────────────
# Order must not change — SLURM_ARRAY_TASK_ID picks sample by position
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
# SLURM_ARRAY_TASK_ID is automatically set by SLURM (0 for first job, 1 for
# second, etc.). We use it to select the corresponding sample from the list.
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

INPUT="$READS/${SAMPLE}.fastq.gz"
OUTPUT="$OUTDIR/${SAMPLE}_trimmed.fastq.gz"

echo "======================================================"
echo "Array task ID : $SLURM_ARRAY_TASK_ID"
echo "Sample        : $SAMPLE"
echo "Input         : $INPUT"
echo "Output        : $OUTPUT"
echo "Start time    : $(date)"
echo "======================================================"

# ── Skip if input file does not exist ────────────────────────────────────────
if [ ! -f "$INPUT" ]; then
    echo "WARNING: $INPUT not found — skipping $SAMPLE"
    exit 0
fi

# ── Run Porechop_ABI ──────────────────────────────────────────────────────────
echo "Trimming adapters: $SAMPLE"

porechop_abi \
    --input "$INPUT" \
    --output "$OUTPUT" \
    --threads 8 \
    --discard_middle

echo ""
echo "Done: $SAMPLE → $OUTPUT"
echo "End time: $(date)"
echo "======================================================"

conda deactivate
