#!/bin/bash
#SBATCH --job-name=decontam
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/decontam_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/decontam_%j.err

# =============================================================================
# run_decontam.sh
# Step 4.1 — Remove contaminant species using negative controls
#
# Run AFTER:  run_krakentools.sh (needs combined_bracken_all_samples.txt)
# Run BEFORE: run_phyloseq.sh and run_vegan.sh
# =============================================================================

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate /data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/envs/r_env

BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1

echo "======================================================"
echo "decontam — Contamination Removal"
echo "Start time: $(date)"
echo "======================================================"

Rscript $BASE/scripts/decontam_analysis.R

if [ $? -eq 0 ]; then
    echo ""
    echo "======================================================"
    echo "decontam complete!"
    echo "End time: $(date)"
    echo "Output: $BASE/r_analysis/decontam/"
    echo ""
    echo "Next steps:"
    echo "  sbatch run_phyloseq.sh"
    echo "  sbatch run_vegan.sh"
    echo "======================================================"
else
    echo "ERROR: decontam_analysis.R failed — check the log"
    exit 1
fi

conda deactivate
