#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64gb
source /home/tmandonnet/.dev_era5/bin/activate

python3 /home/tmandonnet/download_ERA5/ERA5_convert_tg.py