#!/bin/bash
#SBATCH --partition=zen4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
. /home/tmandonnet/.cdsapi_env/bin/activate

for year in {2022..2024..1}
do
    python3 /home/tmandonnet/download_ERA5/download_ERA5.py "$year"

done
