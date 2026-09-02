#!/bin/bash
#SBATCH --job-name=nanoplot_qc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/qc_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/qc_%j.err

# ── Load conda and activate qc environment ──────────────────────────────────
module load Anaconda3/2024.06-1
source activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/qc

# ── Paths ───────────────────────────────────────────────────────────────────
BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1
READS=$BASE/basecalled
OUTDIR=$BASE/qc

mkdir -p "$OUTDIR"
mkdir -p "$BASE/logs"

# ── Sample list ─────────────────────────────────────────────────────────────
SAMPLES=(
    PK-UF-1002 PK-UF-1009 PK-UF-1015 PK-UF-1022 PK-UF-1025
    PK-UF-1028 PK-UF-1044 PK-UF-1064 PK-UF-1069 PK-UF-1070
    PK-UF-1101 PK-UF-1139 PK-UF-1140 PK-UF-1154 PK-UF-1155
    PK-UF-1162 PK-UF-1180 PK-U3-32EW113 PK-U3-32EE122
    PK-U3-42ER128 PK-U3-42EM130 PK-U3-42EW142 PK-U3-2121P14
    PK-U3-2127P15 PK-U3-1211P40 PK-U3-1235P51 PK-U3-4203P78
    PK-U3-4215P84 PK-U3-1101P2 PK-U3-4105P26 PK-UF-1010
    PK-UF-1016 PK-UF-1020 PK-UF-1034 PK-UF-1040 PK-UF-1055
    PK-UF-1059 PK-UF-1075 PK-UF-1078 PK-UF-1111 PK-UF-1131
    PK-UF-1143 PK-UF-1147 PK-UF-1158 PK-UF-1160 PK-UF-1181
    PK-U3-12ES96 PK-U3-12EF97 PK-U3-12EW99 PK-U3-32EF109
    PK-U3-32EE125 PK-U3-1105P3 PK-U3-2101P9 PK-U3-3101P17
    PK-U3-3103P18 PK-U3-3127P23 PK-U3-1215P42 PK-U3-1221P45
    PK-U3-2235P62 PK-U3-3207P67 PK-U3-4211P83
)

# ── Run NanoPlot for each sample ────────────────────────────────────────────
echo "Starting NanoPlot QC for ${#SAMPLES[@]} samples"
echo "Start time: $(date)"

for SAMPLE in "${SAMPLES[@]}"; do
    FASTQ="$READS/${SAMPLE}.fastq.gz"

    if [ ! -f "$FASTQ" ]; then
        echo "WARNING: $FASTQ not found — skipping $SAMPLE"
        continue
    fi

    echo "Processing: $SAMPLE"

    NanoPlot \
        --fastq "$FASTQ" \
        --outdir "$OUTDIR/$SAMPLE" \
        --prefix "${SAMPLE}_" \
        --threads 4 \
        --plots dot \
        --N50 \
        --loglength \
        --title "$SAMPLE"

    echo "Done: $SAMPLE"
done

echo "All QC complete!"
echo "End time: $(date)"
echo "Results saved to: $OUTDIR"
