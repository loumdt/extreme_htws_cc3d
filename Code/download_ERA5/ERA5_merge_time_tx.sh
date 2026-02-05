#!/bin/bash
#SBATCH --partition=zen4
#SBATCH --mem=32G
#SBATCH --time=24:00:00

module load cdo/2.3.0

work_dir="/path/to/dir"

cdo mergetime ${work_dir}/ERA5/t2m/day/ERA5_Europe_day_tx*.nc ${work_dir}/ERA5/ERA5_t2m_tx_Europe_day_0.25deg_1975-2024.nc