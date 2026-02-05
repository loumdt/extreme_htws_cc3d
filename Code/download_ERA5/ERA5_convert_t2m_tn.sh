#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64gb
work_dir="/path/to/dir"

source ${work_dir}/.dev_era5/bin/activate

python3 ${work_dir}/download_ERA5/ERA5_convert_tn.py