#!/bin/bash
#SBATCH --partition=zen4
#SBATCH --time=2:00:00
#SBATCH --mem=32G

output_dir="/path/to/dir"
datadir="/path/to/dir"


source /home/username/.gdal_env/bin/activate

for ssp in {2,5}
do
    echo "SSP ${ssp}"
    for year in {2020..2100..10}
    do
        echo "Year $year"
        FILE_IN_TIF="${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_${year}.tif"
        echo "gdal_translate tif to netCDF"
        gdal_translate -of NetCDF "$FILE_IN_TIF" "${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_${year}.nc"
        echo "cdo sellonlatbox to reduce memory load"
        cdo sellonlatbox,-55,75,15,80 "${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_${year}.nc" "${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_${year}_smallbox.nc"
        rm "${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_${year}.nc"
        echo "python3 add_time_dim.py"
        python3 /home/tmandonnet/CORDEX/add_time_dim.py "${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_${year}_smallbox.nc" $year
        rm "${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_${year}_smallbox.nc"

    done
    echo "cdo mergetime"
    cdo mergetime ${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_*_smallbox_time.nc ${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_merged.nc
    rm ${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_*_smallbox_time.nc
    
    echo "cdo inttime"
    cdo inttime,2020-01-01,00:00:00,1year ${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_merged.nc ${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_2020_2100.nc
    rm ${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_merged.nc

    echo "cdo -remapbil Lambert-Conformal"
    cdo -remapbil,${datadir}/CORDEX/sftlf/CORDEX_EUR-11_land_area_fraction_Lambert_Conformal.nc "${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_2020_2100.nc" "${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_2020_2100_Lambert_Conformal.nc"
    echo "cdo -remapbil Rotated Latitude-Longitude"
    cdo -remapbil,${datadir}/CORDEX/sftlf/CORDEX_EUR-11_land_area_fraction_rotated_latitude_longitude.nc "${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_2020_2100.nc" "${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_2020_2100_rotated_latitude_longitude.nc"
    echo "cdo -remapbil Rotated Pole"
    cdo -remapbil,${datadir}/CORDEX/sftlf/CORDEX_EUR-11_land_area_fraction_rotated_pole.nc "${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_2020_2100.nc" "${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_2020_2100_rotated_pole.nc"
    
    rm "${output_dir}/FPOP/FPOP_SSP${ssp}/FPOP_SSP${ssp}_2020_2100.nc"
    echo "Done"
    done
done