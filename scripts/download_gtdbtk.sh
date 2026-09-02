#!/bin/bash
#SBATCH --job-name=gtdbtk_download
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/gtdbtk_download_%j.log
#SBATCH --error=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/logs/gtdbtk_download_%j.err

DB=/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1/databases/gtdbtk

mkdir -p "$DB"

echo "Downloading GTDB-Tk database..."
echo "Start time: $(date)"

wget https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz \
     -O "$DB/gtdbtk_r232_data.tar.gz"

echo "Extracting..."
tar -xvzf "$DB/gtdbtk_r232_data.tar.gz" -C "$DB" --strip 1 > /dev/null

echo "Cleaning up..."
rm "$DB/gtdbtk_r232_data.tar.gz"

echo "Done! End time: $(date)"
