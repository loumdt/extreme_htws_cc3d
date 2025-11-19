#!/bin/bash
#SBATCH --partition=zen4
#SBATCH --time=12:00:00
#SBATCH --mem=32G

source /home/tmandonnet/.dev_era5/bin/activate
python3 /home/tmandonnet/extreme_htws_cc3d/run_analysis_loop.py