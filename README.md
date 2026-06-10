# Accounting for exposure in 3D spatiotemporally contiguous heatwaves in Europe 1975-2024
Lou Mandonnet, Aglaé Jézéquel, Fabio D'Andrea, Améline Vallet

Working Paper available on HAL: [https://hal.science/hal-05495839](https://hal.science/hal-05495839)

Users will can either download data of the three datasets (ERA5, GHS-POP, EM-DAT) and preprocess data or can download the preprocessed dataset directly on Zenodo, where output data is also available for one sensitivity configuration: [https://doi.org/10.5281/zenodo.20623473](https://doi.org/10.5281/zenodo.20623473). 
Users will also have to set correct folder locations in the scripts.

## Datasets

### ERA5 dataset
[https://cds.climate.copernicus.eu/cdsapp\#!/dataset/reanalysis-era5-single-levels](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download)

### GHS-POP dataset
[https://human-settlement.emergency.copernicus.eu/download.php?ds=pop](https://human-settlement.emergency.copernicus.eu/download.php?ds=pop)

### EM-DAT dataset
[https://public.emdat.be/](https://public.emdat.be/)

## Data preprocessing
### ERA5
Scripts within the download_ERA5 folder can be used to download and preprocess ERA5 data. ```[chosen_variable]``` refers to either tx (daily maximum), tg (daily average) or tn (daily minimum):
1. ```download_ERA5_t2m.sh```
2. ```ERA5_convert_[chosen_variable].sh```
3. ```ERA5_merge_time_[chosen_variable].sh```

### GHS-POP
Scripts within the download_ERA5 folder can be used to download and preprocess ERA5 data:
1. ```GHS_POP_download.sh```
2. ```GHS_POP_format.sh```

## Code
The user can set parameters in the ```run_analysis.py``` script, such as the chosen daily variable (tg, tx, or tn), the temperature variable, the intensity threshold, the duration threshold, the period of study.
Then, running this script will carry out the entire heatwaves detection and overlap analysis process. It calls functions defined in ```utils.py```. Some fine tuning is necessary to plot the correct figures in the ```validate_indices_vs_emdat_impacts``` function, specifically in lines 678-726 of ```utils.py```.

The script ```run_analysis_loop.py``` allows to run the entire process for several parameters combinations.

Some complementary analysis is carried out in ```undetected_and_trends.ipynb```, and should be run after all desired sensitivity configurations of ```run_analysis.py```.
