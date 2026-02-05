#!/bin/bash
#SBATCH --partition=zen4
#SBATCH --time=4:00:00
#SBATCH --mem=64G

homedir="/path/to/homedir"
source ${homedir}/.gdal_env/bin/activate

datadir="/path/to/dir/GHS-POP"
datadir2=/path/to/datadir2"

year_start=1975
year_end=2030

for ((year=year_start;year<=year_end;year=year+5))
do
    echo "Year $year"
    FILE_IN_TIF="$datadir/GHS_POP_E${year}_GLOBE_R2023A_4326_30ss_V1_0.tif"
    echo "gdal_translate tif to netCDF"
    gdal_translate -of NetCDF "$FILE_IN_TIF" "$datadir/GHS_POP_R2023A_4326_30ss_${year}.nc"
    echo "cdo sellonlatbox to reduce memory load"
    cdo sellonlatbox,-13,46,29,72 "$datadir/GHS_POP_R2023A_4326_30ss_${year}.nc" "$datadir/GHS_POP_R2023A_4326_30ss_${year}_smallbox.nc"
    rm "$datadir/GHS_POP_R2023A_4326_30ss_${year}.nc"
    echo "python3 add_time_dim.py"
    python3 ${homedir}/GHS_POP_formatting/add_time_dim.py "$datadir/GHS_POP_R2023A_4326_30ss_${year}_smallbox.nc" $year
    rm "$datadir/GHS_POP_R2023A_4326_30ss_${year}_smallbox.nc"

done
echo "cdo mergetime"
cdo mergetime $datadir/GHS_POP_R2023A_4326_30ss_*_smallbox_time.nc $datadir/GHS_POP_R2023A_4326_30ss_${year_start}_${year_end}_merged.nc
rm $datadir/GHS_POP_R2023A_4326_30ss_*_smallbox_time.nc

echo "cdo inttime"
cdo inttime,${year_start}-01-01,00:00:00,1year $datadir/GHS_POP_R2023A_4326_30ss_${year_start}_${year_end}_merged.nc $datadir/GHS_POP_R2023A_4326_30ss_${year_start}_${year_end}_interpolated.nc
rm $datadir/GHS_POP_R2023A_4326_30ss_${year_start}_${year_end}_merged.nc

cdo -remapcon,${datadir2}/ERA5/t2m/ERA5_t2m_tx_Europe_day_0.25deg_1950-2021.nc "$datadir/GHS_POP_R2023A_4326_30ss_${year_start}_${year_end}_interpolated.nc" "$datadir/GHS_POP_R2023A_${year_start}_${year_end}_ERA5_Europe_grid.nc"

cdo gridarea ${datadir2}/ERA5/t2m/ERA5_t2m_tx_Europe_day_0.25deg_1950-2021.nc ${datadir2}/ERA5/ERA5_Europe_cellarea.nc

echo "Done"