#!/bin/bash
#SBATCH --partition=zen4
#SBATCH --mem=32G
#SBATCH --time=24:00:00

work_dir="/path/to/dir"

. ${work_dir}/.cdsapi_env/bin/activate

## Loop
for ((year = 1975 ; year < 2025 ; year++ ))
do
	python3 ${work_dir}/download_ERA5/download_ERA5_t2m.py $year
done