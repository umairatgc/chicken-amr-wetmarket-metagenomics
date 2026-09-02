#!/bin/bash
#SBATCH --job-name=host_depletion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-62                    # 63 jobs: one per sample/control
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/host_depletion_%A_%a.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/host_depletion_%A_%a.err

# =============================================================================
# run_host_depletion.sh
# Step 2.4 — Host Depletion using minimap2 + samtools
#
# What does this script do?
#   Your filtered reads (from Chopper) contain a mix of:
#     - Bacterial DNA  ← what you WANT to study
#     - Chicken DNA    ← host contamination from the sample source
#     - Human DNA      ← contamination from researcher handling
#   This script removes the chicken and human reads, keeping only the
#   bacterial reads for downstream assembly and analysis.
#
# How does it work?
#   Step 1 — minimap2 aligns every read against the combined host genome
#             (chicken + human). It marks each read as either:
#             MAPPED   = matches host DNA → discard this read
#             UNMAPPED = does not match host → keep this read (it is bacterial)
#
#   Step 2 — samtools extracts only the UNMAPPED reads using the -f 4 flag
#             (-f 4 means "keep reads where flag 4 is set" = unmapped reads)
#
#   Step 3 — samtools converts the output back to FASTQ format for the
#             next pipeline step (assembly)
#
# Input:  filtered/SAMPLE_filtered.fastq.gz    (Chopper output)
# Output: host_depleted/SAMPLE_depleted.fastq.gz
#
# Run AFTER run_chopper.sh has completed and all 63 filtered files exist.
# Verify: ls filtered/ | wc -l  → should be 63
#
# Note on controls:
#   EB-neg and AE-buffer (negative controls) are processed through host
#   depletion too. They should have very few reads — if many reads survive
#   host depletion in the controls, it indicates contamination.
# =============================================================================

# ── Load conda and activate the qc environment ───────────────────────────────
# The qc environment contains both minimap2 and samtools
module load Anaconda3
eval "$(conda shell.bash hook)"   # CRITICAL: enables conda activate in SLURM
conda activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/qc

# ── Paths ────────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
INDIR=$BASE/filtered          # input: Chopper quality-filtered reads
OUTDIR=$BASE/host_depleted    # output: reads with host DNA removed
INDEX=$BASE/databases/host_genomes/host_combined.mmi  # minimap2 index

# Create output directory if it does not exist
mkdir -p "$OUTDIR"

# ── Sample list (61 samples + 2 negative controls) ───────────────────────────
# Order must match run_porechop_abi.sh and run_chopper.sh exactly
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
    EB-neg        # index 61 — negative extraction blank control
    AE-buffer     # index 62 — negative AE buffer control
)

# ── Pick this job's sample using the array task index ────────────────────────
# SLURM automatically sets SLURM_ARRAY_TASK_ID to 0, 1, 2 ... 62
# We use it to select the matching sample from the list above
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

INPUT="$INDIR/${SAMPLE}_filtered.fastq.gz"
OUTPUT="$OUTDIR/${SAMPLE}_depleted.fastq.gz"

echo "======================================================"
echo "Array task ID : $SLURM_ARRAY_TASK_ID"
echo "Sample        : $SAMPLE"
echo "Input         : $INPUT"
echo "Output        : $OUTPUT"
echo "Start time    : $(date)"
echo "======================================================"

# ── Check input file exists ───────────────────────────────────────────────────
# If Chopper did not produce output for this sample, skip rather than crash
if [ ! -f "$INPUT" ]; then
    echo "WARNING: $INPUT not found — skipping $SAMPLE"
    echo "         Check that Chopper completed successfully for this sample"
    exit 0
fi

# Count reads going in so we can compare with reads coming out
READS_IN=$(gunzip -c "$INPUT" | wc -l)
READS_IN=$((READS_IN / 4))   # FASTQ format = 4 lines per read
echo "Reads entering host depletion: $READS_IN"

# ==============================================================================
# STEP 1 + 2 — Align reads to host genome, then keep only unmapped reads
#
# This is a pipeline of three commands chained with pipes (|):
#
# minimap2:
#   -a              : output in SAM format (standard alignment format)
#   -x map-ont      : use ONT preset (optimised for long noisy reads)
#   -t 8            : use 8 CPU threads
#   --secondary=no  : do not report secondary alignments (one result per read)
#   $INDEX          : the pre-built host genome index file (.mmi)
#   $INPUT          : your sample's filtered reads
#
# samtools view:
#   -b              : output in BAM format (compressed binary version of SAM)
#   -f 4            : keep ONLY reads with flag 4 set = unmapped reads
#                     (flag 4 = "this read did not align to the reference")
#   -@  8           : use 8 threads for compression
#
# samtools fastq:
#   converts the surviving unmapped reads back to FASTQ format
#   -@ 8            : use 8 threads
#   | gzip          : compress the output
# ==============================================================================
echo ""
echo "Running minimap2 alignment and samtools filtering..."

minimap2 \
    -a \
    -x map-ont \
    -t 8 \
    --secondary=no \
    "$INDEX" \
    "$INPUT" | \
samtools view \
    -b \
    -f 4 \
    -@ 8 | \
samtools fastq \
    -@ 8 | \
gzip > "$OUTPUT"

# ── Verify output and report how many reads were removed ─────────────────────
if [ -f "$OUTPUT" ]; then
    # Count reads that survived host depletion
    READS_OUT=$(gunzip -c "$OUTPUT" | wc -l)
    READS_OUT=$((READS_OUT / 4))

    # Calculate how many host reads were removed
    READS_REMOVED=$((READS_IN - READS_OUT))

    # Calculate percentage of reads kept
    if [ "$READS_IN" -gt 0 ]; then
        PCT_KEPT=$(awk "BEGIN {printf \"%.1f\", ($READS_OUT/$READS_IN)*100}")
    else
        PCT_KEPT="0.0"
    fi

    echo ""
    echo "======================================================"
    echo "Done: $SAMPLE"
    echo "  Reads in (filtered)    : $READS_IN"
    echo "  Reads removed (host)   : $READS_REMOVED"
    echo "  Reads out (bacterial)  : $READS_OUT  ($PCT_KEPT% kept)"
    echo "  Output file            : $OUTPUT"
    echo "  File size              : $(du -sh "$OUTPUT" | cut -f1)"
    echo "End time: $(date)"
    echo "======================================================"
else
    echo "ERROR: Output file not created for $SAMPLE"
    exit 1
fi

conda deactivate
