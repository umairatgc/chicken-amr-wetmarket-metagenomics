#!/bin/bash
#SBATCH --job-name=phyloseq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/phyloseq_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/phyloseq_%j.err

# =============================================================================
# run_phyloseq.sh
# Step 4.2 — Alpha & Beta Diversity Analysis using phyloseq
#
# Run AFTER:  run_decontam.sh (needs ps_all_samples.rds)
# Run BEFORE: run_vegan.sh (which uses the Bray-Curtis distance matrix)
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/r_env

BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1

echo "======================================================"
echo "phyloseq — Alpha & Beta Diversity Analysis"
echo "Start time: $(date)"
echo "======================================================"

Rscript $BASE/scripts/phyloseq_analysis.R

if [ $? -eq 0 ]; then
    echo ""
    echo "======================================================"
    echo "phyloseq complete!"
    echo "End time: $(date)"
    echo "Output: $BASE/r_analysis/phyloseq/"
    echo ""
    echo "Next step:"
    echo "  sbatch run_vegan.sh"
    echo "======================================================"
else
    echo "ERROR: phyloseq_analysis.R failed — check the log"
    exit 1
fi

conda deactivate
