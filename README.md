# Extreme heatwaves in Europe 1950-2021: analysis of the links between meteorology, population, and impacts

Disclaimer : these scripts would benefit a rewrite using xarray instead of netCDF4.

## Data

### E-OBS dataset
https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php#datafiles

### ERA5 dataset
https://cds.climate.copernicus.eu/cdsapp\#!/dataset/reanalysis-era5-single-levels

### EM-DAT dataset
https://public.emdat.be/data

### GHS-POP dataset
https://ghsl.jrc.ec.europa.eu/ghs_pop2019.php

## Data requirements
All data files uploaded here are necessary to run scripts. Users will also have to download data of the four datasets mentioned above (EM-DAT, ERA5 and/or E-OBS, GHS-POP), and preprocess data. ERA5 data should be formatted into daily data. Other preprocessing consists in running the scripts in Pop_formatting and GDIS_EM-DAT formatting. Users will have to check file names and locations to match with scripts. 

## Code
Set parameters in run_all_detection_overlap_analysis.py and run this script to carry out the entire heatwaves detection and overlap analysis process. It calls functions in detection_overlap_functions.py and analysis_classification_plot_functions.py.

The script run_all_detection_overlap_analysis_loop.py allows to run the entire process for several parameters combinations.