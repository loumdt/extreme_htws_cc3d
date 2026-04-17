#%%
from os.path import join,exists
from pathlib import Path #check existence and create folders

from utils import *
#%%
if __name__ == "__main__":
    temp_variable = 't2m' # 't2m', 'wbgt' or 'utci'
    daily_var = 'tx' # 'tg', 'tn' or 'tx' (mean, min, max), default value is 'tg'
    start_year = 1975 # beginning of the studied period, default 1975
    end_year = 2021 # end of the studied period, default 2021

    start_year_ref = 1975 #beginning of the seasonal_cycle period, default 1975
    end_year_ref = 2021 #end of the seasonal_cycle period, default 2021

    anomaly = True #If True, seasonal_cycle is removed from temperature data before computing percentiles. If False, absolute values are used. Default is True
    # I strongly recommend to keep anomaly True if using relative thresholds, cf https://doi.org/10.1038/s41467-024-46349-x
    nb_days = 4 #nb of days used as a heatwave duration threshold, default value is 4
    relative_threshold = True #If True, the threshold value should be a percentile and locally defined. If False, the threshold should be an absolute value. Default is True
    threshold_value = 95 # If relative_threshold is True, percentile used as a threshold for heatwave occurence, value in ]0;100[, default value is 95 (consistent with default value of relative_threshold). If absolute value, expected to be a value in °C.
    distrib_window_size = 15 #size (in days) of the temporal window that is used to compute the temperature distribution (on which is based the threshold) of each calendar day, default value is 15
    flex_time_span = 3 #In order to account for potential EM-DAT imprecisions, set a flexibility window of flex_time_span days, default value is 3
    dust_threshold=775
    nb_top_events=10 #number of top detected events to look for in the litterature
    connectivity = 26
    
    merge_method='weighted' #'weighted', 'equal' or 'strongest'
    dont_skip_figure=True # If True, plot figures in validate_indices_vs_emdat_impacts. Bubble plot needs tuning to avoid errors

    name_dict_threshold = {True : 'th', False : 'C'} #If relative threshold, value is a percentile; if absolute threshold, value is in °C

    read_directory = "Data"

    pop_file_path = join(read_directory,"GHS-POP","GHS_POP_R2023A_1975_2030_ERA5_Europe_grid.nc")
    emdat_file_path = join(read_directory,"EM-DAT","EMDAT_Europe_Turkey-1975-2021-heatwaves.xlsx")

    print("daily_var:", daily_var)
    print("nb_days:", nb_days)
    print("threshold_value:", threshold_value)
    print("connectivity:", connectivity)
    print(f"study period: {start_year}-{end_year}")
    print(f"baseline period: {start_year_ref}-{end_year_ref}")
    print(f"merge_method: {merge_method}")
    write_directory = join(read_directory,f"ERA5_{temp_variable}_{daily_var}_{threshold_value}{name_dict_threshold[relative_threshold]}_connec_{connectivity}_{nb_days}days_flex_{flex_time_span}d_{start_year}_{end_year}_ref_{start_year_ref}_{end_year_ref}{'_anomaly'*anomaly}_merge_{merge_method}")
    Path(join(write_directory,'figs')).mkdir(parents=True, exist_ok=True)
    temp_file_path = join(read_directory,"ERA5",temp_variable,f"ERA5_{temp_variable}_{daily_var}_Europe_day_0.25deg_1975-2024.nc")
    if distrib_window_size%2==0:
        raise ValueError('distrib_window_size is even. It has to be odd so the window can be centered on the computed day.')
    if relative_threshold==False and anomaly==True:
        raise ValueError("Using an absolute threshold can only work by working on absolute values, and not anomalies. The parameter 'anomaly' should be set to False.")
    if relative_threshold==False and threshold_value>=60 and threshold_value<100 :
        raise ValueError("It seems that the value given for the threshold is a percentile and not an absolute value, but the parameter 'relative_threshold' is set to False.")

    overwrite_files=False #If True, overwrite output files that already exists (may be relevant in case of code or data update)
    # if overwrite_file is True or if output file does not exist : call function ; else pass
    if (overwrite_files or exists(join(write_directory,f"seasonal_cycle.nc"))==False) :
        print("\n Running compute_seasonal_cycle... \n")
        compute_seasonal_cycle(write_directory=write_directory,temp_file_path=temp_file_path,start_year_ref=start_year_ref,end_year_ref=end_year_ref,temp_variable=temp_variable)

    if (overwrite_files or exists(join(write_directory,f"distrib_threshold_{threshold_value}.nc"))==False) and relative_threshold==True : # Only used for relative thresholds
        print("\n Running compute_distrib_percentile... \n")
        compute_distrib_percentile(write_directory=write_directory,temp_file_path=temp_file_path,start_year_ref=start_year_ref,end_year_ref=end_year_ref,temp_variable=temp_variable,threshold_value=threshold_value,distrib_window_size=distrib_window_size,anomaly=anomaly)

    if overwrite_files or exists(join(write_directory,"labels_cc3d.nc"))==False :
        print("\n Running cc3d_scan_heatwaves... \n")
        cc3d_scan_heatwaves(read_directory=read_directory,write_directory=write_directory,temp_file_path=temp_file_path,start_year=start_year,end_year=end_year,temp_variable=temp_variable,threshold_value=threshold_value,relative_threshold=relative_threshold,anomaly=anomaly,nb_days=nb_days,dust_threshold=dust_threshold,connectivity=connectivity)
        
    if overwrite_files or exists(join(write_directory,"Russo_HWMId.nc"))==False :
        print("\n Running compute_Russo_HWMId... \n")
        compute_Russo_HWMId(write_directory=write_directory,temp_file_path=temp_file_path,start_year=start_year,end_year=end_year,start_year_ref=start_year_ref,end_year_ref=end_year_ref,temp_variable=temp_variable)

    if overwrite_files or exists(join(write_directory,"df_htws.csv"))==False :
        print("\n Running create_heatwaves_indices_database... \n")
        create_heatwaves_indices_database(read_directory=read_directory,write_directory=write_directory,temp_file_path=temp_file_path,pop_file_path=pop_file_path,start_year=start_year,end_year=end_year,temp_variable=temp_variable,threshold_value=threshold_value,anomaly=anomaly)

    if overwrite_files or exists(join(write_directory,"df_htws_step2.csv"))==False :
        print("\n Running analyze_emdat_overlap... \n")
        analyze_emdat_overlap(read_directory=read_directory,write_directory=write_directory,emdat_file_path=emdat_file_path,flex_time_span=flex_time_span,start_year=start_year,end_year=end_year,merge_method=merge_method)

    if overwrite_files or exists(join(write_directory,'figs',"distrib_4idx.pdf"))==False :
        print("\n Running validate_indices_vs_emdat_impacts... \n")
        validate_indices_vs_emdat_impacts(read_directory=read_directory,write_directory=write_directory,emdat_file_path=emdat_file_path,pop_file_path=pop_file_path,start_year=start_year,end_year=end_year,temp_variable=temp_variable,daily_var=daily_var,start_year_ref=start_year_ref,end_year_ref=end_year_ref,anomaly=anomaly,nb_days=nb_days,threshold_value=threshold_value,relative_threshold=relative_threshold,flex_time_span=flex_time_span,connectivity=connectivity,dont_skip_figure=dont_skip_figure)

    if overwrite_files or exists(join(write_directory,f"df_htws_top_events_NOCHANGE.csv"))==False :
        print("\n Running analysis_top_detected_events... \n")
        analysis_top_detected_events(read_directory=read_directory,write_directory=write_directory)