#!/bin/bash
#SBATCH --job-name=vegan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/vegan_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/vegan_%j.err

# =============================================================================
# run_vegan.sh
# Step 4.3 — PERMANOVA, Betadispersion, and NMDS using vegan
#
# Run AFTER:  run_phyloseq.sh (needs ps_real_samples.rds)
# This is the final R analysis step in the pipeline.
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/r_env

BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1

echo "======================================================"
echo "vegan — PERMANOVA, Betadispersion & NMDS"
echo "Start time: $(date)"
echo "======================================================"

Rscript $BASE/scripts/vegan_analysis.R

if [ $? -eq 0 ]; then
    echo ""
    echo "======================================================"
    echo "vegan complete!"
    echo "End time: $(date)"
    echo "Output: $BASE/r_analysis/vegan/"
    echo ""
    echo "Pipeline complete! All R analysis steps finished:"
    echo "  r_analysis/decontam/   — contamination removal results"
    echo "  r_analysis/phyloseq/   — alpha & beta diversity results"
    echo "  r_analysis/vegan/      — PERMANOVA & NMDS results"
    echo "======================================================"
else
    echo "ERROR: vegan_analysis.R failed — check the log"
    exit 1
fi

conda deactivate
