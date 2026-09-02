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
    AE-buffer EB-neg
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
