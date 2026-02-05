#!/bin/bash
#SBATCH --partition=zen4
#SBATCH --mem=32G
#SBATCH --time=24:00:00

module load cdo/2.3.0

work_dir="/path/to/dir"
for ((year = 1975 ; year < 2025 ; year++ ))
do
	cdo daymean ${work_dir}/ERA5/t2m/ERA5_NorthAtlantic_hour_t2m_${year}010100-${year}123123.nc ${work_dir}/ERA5/t2m/ERA5_NorthAtlantic_day_t2m_tg_${year}.nc
done
