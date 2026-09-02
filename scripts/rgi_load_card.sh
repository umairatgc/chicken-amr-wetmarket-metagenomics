#!/bin/bash
#SBATCH --job-name=rgi_load_card
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/rgi_load_card_%j.out
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/rgi_load_card_%j.err

BASE=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1

module load Anaconda3
source /apps/system/easybuild/software/Anaconda3/2025.06-1/etc/profile.d/conda.sh
conda activate ${BASE}/envs/rgi

echo "========================================"
echo "Loading CARD database into RGI"
echo "CARD path: ${BASE}/databases/card/card.json"
echo "Started: $(date)"
echo "========================================"

rgi load \
    --card_json ${BASE}/databases/card/card.json

echo ""
echo "=== Verifying database ==="
rgi database --version

echo ""
echo "Done: $(date)"
