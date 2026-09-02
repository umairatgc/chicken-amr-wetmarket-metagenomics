#!/bin/bash
#SBATCH --job-name=phase5_05_combined_csvs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/phase5_05_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/phase5_05_%j.err

# =============================================================================
# Phase 5 — Step 5: Combined AMR + MGE + Taxonomy + Metadata CSVs
# Submit: sbatch scripts/run_phase5_05.sh
# Monitor: tail -f logs/phase5_05_JOBID.log
# =============================================================================

BASE="/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
ENV="$BASE/envs/r_env"
SCRIPT="$BASE/scripts/phase5_05_combined_csvs.py"

echo "========================================"
echo "  Job ID   : $SLURM_JOB_ID"
echo "  Node     : $SLURMD_NODENAME"
echo "  Started  : $(date)"
echo "  Script   : $SCRIPT"
echo "========================================"

# Create logs directory if it doesn't exist
mkdir -p "$BASE/logs"

# Activate conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV"

echo "Python   : $(which python3)"
echo "Conda env: $CONDA_DEFAULT_ENV"
echo ""

# Run the script
python3 "$SCRIPT"
EXIT_CODE=$?

echo ""
echo "========================================"
echo "  Finished : $(date)"
echo "  Exit code: $EXIT_CODE"
echo "========================================"

exit $EXIT_CODE
