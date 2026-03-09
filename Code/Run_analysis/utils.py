import numpy as np
import xarray as xr
from tqdm import tqdm #create a user-friendly feedback while script is running
from os.path import join
import pandas as pd #handle dataframes
import cc3d #connected components patterns
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as clrs
from scipy import stats
from sklearn import metrics
import cftime
from datetime import datetime
from ast import literal_eval
import subprocess
import json
from sklearn.linear_model import LinearRegression

def compute_seasonal_cycle(write_directory,temp_file_path,start_year_ref=1975,end_year_ref=2021,temp_variable='t2m') :
    '''This function computes a seasonal_cycle for each calendar day of the year.
    By default, the seasonal_cycle is computed over 1975-2021.
    This function can be used with several models and variables.'''

    # Load dataset
    ds = xr.open_dataset(join(temp_file_path), engine='netcdf4')
    da = getattr(ds, temp_variable)
    # Keep only years of interest
    da.sel(time=(da.time.dt.year>=start_year_ref) & (da.time.dt.year<=end_year_ref))

    # Drop Feb 29
    da = da.convert_calendar("noleap")
    # Group using dayofyear and sum to compute mean at the end
    seasonal_cycle = da.groupby(da.time.dt.dayofyear).mean(dim="time")

    # Export data to netcdf file
    seasonal_cycle.to_netcdf(join(write_directory,f"seasonal_cycle.nc"))

    # Clear resources
    da.close()
    seasonal_cycle.close()
    ds.close()
    return


def compute_distrib_percentile(write_directory,temp_file_path,start_year_ref=1975,end_year_ref=2021,temp_variable='t2m',threshold_value=95,distrib_window_size=15,anomaly=True) :
    '''This function computes, for every calendar day, the n-th (n is the threshold_value, default 95) percentile of the corresponding distribution of daily variable. 
    By default, the distribution is computed over 1975-2021.'''

    if distrib_window_size%2==0:
        raise ValueError('distrib_window_size is even. It has to be odd so the window can be centered on the computed day.')

    ds = xr.open_dataset(join(temp_file_path), engine="netcdf4")

    # Load seasonal_cycle file to create output data structure and compute anomaly
    seasonal_cycle = xr.open_dataarray(join(write_directory,f"seasonal_cycle.nc"), engine='netcdf4')
    
    # Initialize data array with the first file
    da = getattr(ds, temp_variable) # Iterate over files, except first one which has already been used in initialization
    # Drop 29 Feb and only keep the reference period
    da = da.convert_calendar("noleap")
    da = da.sel(time=(da.time.dt.year>=start_year_ref) & (da.time.dt.year<=end_year_ref))

    # Create threshold table by copying seasonal_cycle table, values will be updated later
    threshold = seasonal_cycle.copy()

    if anomaly :
        for year in tqdm(range(len(da.time)//365)) : # Iterate over the number of years
            da[year*365:(year+1)*365,:,:] = da[year*365:(year+1)*365,:,:] - seasonal_cycle.data # Compute anomaly
    else :
        seasonal_cycle.close()
    
    for day in tqdm(range(1,366)) : # Calendar days ranging from 1 to 365 (no leap years)
        if day - (distrib_window_size//2) <= 0  : # day is at the beginning of January, window overlapping with December
            distrib_start_bound = 365 + day - (distrib_window_size//2) # start bound of the distribution window
            distrib_end_bound = day + (distrib_window_size//2) # end bound of the distribution window
            mask = (da.time.dt.dayofyear >= distrib_start_bound) + (da.time.dt.dayofyear <= distrib_end_bound) # Create a boolean mask for days between distrib_start_bound and distrib_end_bound
        elif day + (distrib_window_size//2) > 365 : # day is at the end of December , window overlapping with January
            distrib_start_bound = day - (distrib_window_size//2)
            distrib_end_bound = day + (distrib_window_size//2) - 365
            mask = (da.time.dt.dayofyear >= distrib_start_bound) + (da.time.dt.dayofyear <= distrib_end_bound) # Create a boolean mask for days between distrib_start_bound and distrib_end_bound
        else :
            distrib_start_bound = day - (distrib_window_size//2)
            distrib_end_bound = day + (distrib_window_size//2)
            mask = (da.time.dt.dayofyear >= distrib_start_bound) & (da.time.dt.dayofyear <= distrib_end_bound) # Create a boolean mask for days between distrib_start_bound and distrib_end_bound
        # Apply the mask to select the corresponding days and compute percentile, take day-1 because day ranges from 1 to 365 but python indexing ranges from 0 to 364
        threshold.data[day-1,:,:] = np.percentile(da.sel(time=mask).data,threshold_value,axis=0)

    # Export data to netcdf file
    threshold.to_netcdf(join(write_directory,f"distrib_threshold_{threshold_value}.nc"))

    # Clear resources
    da.close()
    threshold.close()
    if anomaly :
        seasonal_cycle.close()
    ds.close()
    return

def cc3d_scan_heatwaves(read_directory,write_directory,temp_file_path,start_year=1975,end_year=2021,temp_variable='t2m',threshold_value=95,relative_threshold=True,anomaly=True,nb_days=4,dust_threshold=775,connectivity=26) :
    '''This function carries out a cc3d scan (https://pypi.org/project/connected-components-3d/) to detect heatwaves in the meteorological database (default ERA5, t2m, tx).
    The heatwaves point are labeled with a number corresponding to a heatwave identifier.
    Otherwise, values are set to -9999.'''
    # Set dust threshold to supress small amount of points. The threshold of 775 have been established by testing with different values and set to the point that the number of heatwaves did not change by more than 1% between two increment of 25.

    
    ds = xr.open_dataset(join(temp_file_path), engine="netcdf4")

    if anomaly :
        seasonal_cycle = xr.open_dataarray(join(write_directory,f"seasonal_cycle.nc"), engine='netcdf4')
        # Keep only JJA values
        mask = (seasonal_cycle.dayofyear>=152) & (seasonal_cycle.dayofyear<=243) # dayofyear ranges from 1 to 365 ; 152 is June 1st, 243 is August 31st
        seasonal_cycle = seasonal_cycle.sel(dayofyear=mask)

    if relative_threshold : # Load temperature threshold for reference period :
        threshold = xr.open_dataarray(join(write_directory,f"distrib_threshold_{threshold_value}.nc"), engine='netcdf4')
        # Keep only JJA values
        mask = (threshold.dayofyear>=152) & (threshold.dayofyear<=243) 
        threshold = threshold.sel(dayofyear=mask)
    else : # If absolute threshold, only need a scalar, not a 3D array
        threshold = threshold_value

    N_labels = 0 # Count the numbers of patterns
    da = getattr(ds, temp_variable)
    # Drop 29 Feb and correct day of year
    da = da.convert_calendar("noleap")
    # Keep only JJA values and years of interest
    da = da.sel(time=(da.time.dt.year>=start_year) & (da.time.dt.year<=end_year))
    da = da.sel(time=(da.time.dt.season=='JJA'))

    # Initialize output array 
    label = da.copy()

    print("Computing cc3d.connected_components labels and dusting...")

    for year in tqdm(range(end_year-start_year+1)) :# Iterate over the years
        da_year = da[year*92:(year+1)*92,:,:]# Select data for the given year
        if anomaly : # Substract seasonal_cycle to compute anomaly
            da_year = da_year - seasonal_cycle.data
        da_year = da_year*(da_year>threshold.data)-9999*(da_year<=threshold.data) # Set to -9999 the values that do not exceed threshold
        stack_temp = -9999*np.ones((92,np.shape(da_year)[1],np.shape(da_year)[2])) # Create a 3D array that will hold the temperature values when and where there are heatwaves
        stack_where = np.zeros((np.shape(da_year)[1],np.shape(da_year)[2])) # Create a 2D array that holds the number of consecutive hot days for each location, computed for each day
        for day in range(92):
            stack_where[:,:] = stack_where[:,:] + np.ones((np.shape(da_year)[1],np.shape(da_year)[2]))*(da_year[day,:,:]!=-9999) # Add one day to each potential heatwave location
            stack_where[:,:] = stack_where[:,:]*(da_year[day,:,:]!=-9999) # When not adding a day, have to set back the duration to zero
            if day>=nb_days-1 :
                stack_temp[day-(nb_days-1):day+1,:,:] = stack_temp[day-(nb_days-1):day+1,:,:]*(stack_temp[day-(nb_days-1):day+1,:,:]!=-9999)+da_year[day-(nb_days-1):day+1,:,:]*((stack_where>=nb_days)*stack_temp[day-(nb_days-1):day+1,:,:]==-9999)+(-9999*((stack_where<nb_days)*(stack_temp[day-(nb_days-1):day+1,:,:]==-9999))) #record the last four days for the corresponding scanning window, add new consecutive hot days and set not hot days to -9999

        # Compute connected components for the remaining values of stack_temp
        
        labels_in = cc3d.dust((stack_temp!=-9999),dust_threshold)
        labels_out, N_added = cc3d.connected_components(labels_in,connectivity=connectivity,return_N=True) # Return the table of labels and the number of added patterns
        
        # Record labels and add N_labels offset where labels are nonzero
        label.data[year*92:(year+1)*92,:,:] = labels_out + N_labels*(labels_out>0)

        #Remove sea heatwaves
        label.data[year*92:(year+1)*92,:,:] = remove_outside_heatwaves(read_directory=read_directory,labels=label[year*92:(year+1)*92,:,:])
        
        # Update N_labels
        N_labels += N_added

    labels_list = np.unique(label.data)
    labels_list = labels_list[np.where(labels_list!=0)]
    print(len(labels_list),"heatwaves detected")
    # Set DataArray to Dataset to set variable name
    label = label.to_dataset(name="label")
    # Save to netCDF 
    label.to_netcdf(join(write_directory,f"labels_cc3d.nc"))
    
    # Clear resources
    da.close()
    ds.close()
    label.close()
    threshold.close()
    if anomaly :
        seasonal_cycle.close()
    return

def remove_outside_heatwaves(read_directory,labels,dust_threshold=775,connectivity=26) :
    '''Remove sea heatwaves and heatwaves occurring outside continental Europe
    '''
    
    land_mask = xr.open_dataset(join(read_directory,'ERA5','Mask',f"Mask_Europe_land_only_ERA5_0.25deg.nc"),engine='netcdf4').mask

    # Remove sea points, North Africa points, and Middle East points that are South of Turkey
    labels = labels * (land_mask.data==0)

    # Remove small heatwaves with the cc3d dust function
    labels = labels * cc3d.dust((labels>0),dust_threshold,connectivity=connectivity)
    land_mask.close()
    return labels

def compute_Russo_HWMId(write_directory,temp_file_path,start_year=1975,end_year=2021,temp_variable='t2m',anomaly=True) :
    """Compute the Russo_HWMId index map.
    Based on HWMId defined by Russo et al (2015, https://dx.doi.org/10.1088/1748-9326/10/12/124003 )."""

    ds = xr.open_dataset(join(temp_file_path), engine='netcdf4')
    da = getattr(ds, temp_variable)
    # Drop Feb 29 and only keep JJA days and years of interest
    da = da.convert_calendar("noleap")
    da = da.sel(time=(da.time.dt.year>=start_year) & (da.time.dt.year<=end_year))
    da = da.sel(time=(da.time.dt.season=='JJA'))
    
    # if anomaly :
    #     # Keep only JJA values
    #     seasonal_cycle = xr.open_dataarray(join(write_directory,f"seasonal_cycle.nc"), engine='netcdf4')
    #     mask = (seasonal_cycle.dayofyear>=152) & (seasonal_cycle.dayofyear<=243) # dayofyear ranges from 1 to 365 ; 152 is June 1st, 243 is August 31st
    #     seasonal_cycle = seasonal_cycle.sel(dayofyear=mask)
    #     for year in tqdm(range(len(da.time)//92)) : # Iterate over the number of years
    #         da[year*92:(year+1)*92,:,:] = da[year*92:(year+1)*92,:,:] - seasonal_cycle.data # Compute anomaly
        
    temp_25p = np.percentile(da.groupby(da.time.dt.year).max(), 25, axis=0)
    temp_75p = np.percentile(da.groupby(da.time.dt.year).max(), 75, axis=0)
    
    Russo_HWMId = da.copy()
    Russo_HWMId.data = np.maximum((da - temp_25p)/(temp_75p - temp_25p), 0)
    
    Russo_HWMId.to_dataset(name="HWMId")
    Russo_HWMId.to_netcdf(join(write_directory,"Russo_HWMId.nc"))

    Russo_HWMId.close()
    da.close()
    ds.close()
    #if anomaly :
    #    seasonal_cycle.close()
    return


def create_heatwaves_indices_database(read_directory,write_directory,temp_file_path,pop_file_path,start_year=1975,end_year=2021,temp_variable='t2m',threshold_value=95,anomaly=True):
    '''This function is used to create the dataset of the indices of the detected heatwaves. The set of detected heatwaves depends on all the parameters.'''

    # Load temperature-related data
    ds_labels = xr.open_dataset(join(write_directory,f"labels_cc3d.nc"),engine='netcdf4')
    da_labels = ds_labels.label
    ds_temp = xr.open_dataset(join(temp_file_path), engine='netcdf4')
    da_temp = getattr(ds_temp,temp_variable)
    da_temp = da_temp.convert_calendar("noleap") # Drop 29 Feb
    da_threshold = xr.open_dataarray(join(write_directory,f"distrib_threshold_{threshold_value}.nc"),engine='netcdf4')
    da_HWMId = xr.open_dataarray(join(write_directory,"Russo_HWMId.nc"),engine='netcdf4')

    # Keep only JJA values; labels and HWMId are already JJA-only
    da_temp = da_temp.sel(time=(da_temp.time.dt.year>=start_year) & (da_temp.time.dt.year<=end_year))
    da_temp = da_temp.sel(time = (da_temp.time.dt.season=='JJA'))
    mask = (da_threshold.dayofyear>=152) & (da_threshold.dayofyear<=243) # dayofyear ranges from 1 to 365 ; 152 is June 1st, 243 is August 31st
    da_threshold = da_threshold.sel(dayofyear=mask)

    if anomaly :
        seasonal_cycle = xr.open_dataarray(join(write_directory,f"seasonal_cycle.nc"), engine='netcdf4')
        # Keep only JJA values
        seasonal_cycle = seasonal_cycle.sel(dayofyear=mask)
        for year in tqdm(range(len(da_temp.time)//92)) : # Iterate over the number of years
            da_temp[year*92:(year+1)*92,:,:] = da_temp[year*92:(year+1)*92,:,:] - seasonal_cycle.data # Compute anomaly

    # Load cell area
    ds_cell_area = xr.open_dataset(join(read_directory,"ERA5","ERA5_Europe_cellarea.nc"),engine='netcdf4') # Area of each grid cell, in m²
    da_cell_area = ds_cell_area.cell_area/1e6 # Load DataArray and convert to km²

    # Load population data
    ds_pop = xr.open_dataset(join(pop_file_path), engine='netcdf4') # Population count, need to convert to density
    da_pop = ds_pop.Band1*1000 # Data is in thousands of people and we want the real number of people

    da_pop_density = da_pop/da_cell_area # Population density in person/km²

    labels_list = np.unique(da_labels.data)
    labels_list = labels_list[np.where(labels_list!=0)] # Remove zero which corresponds to the absence of hot point, not a heatwave label ID
    df_htws = pd.DataFrame(columns=['Year','Start Date','End Date','Intensity','Spatial extent','Duration','Max','Temp_sum','HWMId_sum','HWMId_mean','Exposed_population','Total Exposed_population','Temp_pop','HWMId_pop','HWMId_pop_mean','EM-DAT heatwaves','EM-DAT Total Deaths'],index=[int(i) for i in labels_list],data=None)    

    for year in range(end_year-start_year+1) : # Compute temperature exceedance relative to the threshold, iterate over years (92 JJA days per year)
        da_temp[year*92:(year+1)*92,:,:] = da_temp[year*92:(year+1)*92,:,:] - da_threshold.data

    # Compute weights for latitude-weighted mean
    weights = np.cos(np.deg2rad(da_temp.lat))
    weights.name = "weights"

    for label in tqdm(labels_list) : # Iterate on heatwaves
        da_bool_htw = da_labels.where(da_labels==label, drop=True).fillna(0)>0 # Select days and grid points for the heatwave of interest and convert to bool array
        da_temp_htw = da_temp.where(da_labels==label, drop=True)
        da_HWMId_htw = da_HWMId.where(da_labels==label, drop=True)

        year_event = da_temp_htw.time.dt.year.data[0]
        df_htws.loc[label,'Year'] = year_event
        df_htws.loc[label,'Start Date'] = da_temp_htw.time.data[0]
        df_htws.loc[label,'End Date'] = da_temp_htw.time.data[-1]

        labels_bool_2D = np.max(da_labels==label,axis=0) # Squeeze heatwave labels on a boolean 2D-map to see maximum spatial extension
        da_pop_htw = da_pop.sel(time=(da_pop.time.dt.year==year_event))
        da_pop_htw = da_pop_htw.where(labels_bool_2D,drop=True)
        da_pop_density_htw = da_pop_density.sel(time=(da_pop.time.dt.year==year_event))
        da_pop_density_htw = da_pop_density_htw.where(labels_bool_2D,drop=True)
    
        df_htws.loc[label,'Intensity'] = da_temp_htw.weighted(weights).mean().data
        df_htws.loc[label,'Spatial extent'] = (da_cell_area*labels_bool_2D).sum().data
        df_htws.loc[label,'Duration'] = len(da_temp_htw.time)
        df_htws.loc[label,'Max'] = da_temp_htw.max().data
        df_htws.loc[label,'Temp_sum'] = da_temp_htw.weighted(weights).sum().data
        df_htws.loc[label,'HWMId_sum'] = da_HWMId_htw.weighted(weights).sum().data
        df_htws.loc[label,'HWMId_mean'] = da_HWMId_htw.weighted(weights).mean().data
        df_htws.loc[label,'Exposed_population'] = da_pop_htw.sum().data/1e3 # Compute population in thousands to avoid later memory issues with bootstrap and MK test 
        df_htws.loc[label,'Total Exposed_population'] = (da_bool_htw*da_pop_htw.data).sum().data/1e3 # Compute population in thousands to avoid later memory issues with bootstrap and MK test 
        df_htws.loc[label,'Temp_pop'] = (da_temp_htw*da_pop_density_htw.data).weighted(weights).sum().data
        df_htws.loc[label,'HWMId_pop'] = (da_HWMId_htw*da_pop_density_htw.data).weighted(weights).sum().data
        df_htws.loc[label,'HWMId_pop_mean'] = (da_HWMId_htw*da_pop_density_htw.data).weighted(weights).mean().data

    # Save dataframe 
    df_htws.to_csv(join(write_directory,f"df_htws.csv"))
    
    ds_labels.close()
    ds_temp.close()
    da_threshold.close()
    ds_cell_area.close()
    ds_pop.close()
    if anomaly :
        seasonal_cycle.close()
    return


def analyze_emdat_overlap(read_directory,write_directory,emdat_file_path,flex_time_span=3,start_year=1975,end_year=2021):
    '''This function is used to analyse the spatial and temporal overlap between EM-DAT heatwaves and the meteorological database heatwaves (default ERA5) detected with the CC3D scan.
    The detection threshold depends on the parameters used precedently, which is why all these parameters are required.
    This function can be used with several variables'''

    # Load EM-DAT dataset
    df_emdat = pd.read_excel(join(emdat_file_path),header=0, index_col=0) 
    df_emdat = df_emdat[(df_emdat['Start Year']>=start_year) & (df_emdat['Start Year']<=end_year)] # Only keep events of the studied period (default 1975-2021)
    df_emdat = df_emdat[(df_emdat['Start Month'].isin([6,7,8])) | (df_emdat['End Month'].isin([6,7,8]))] # Remove heatwaves occurring outside JJA period
    df_emdat = df_emdat[['ISO','Country','Subregion','Start Year','Start Month','Start Day','End Year','End Month','End Day','Total Deaths']]
    df_emdat['overlap CC3D'] = None
    # Load temperature-related data in case indices need to be recomputed (only occur if one EM-DAT heatwave matches several CC3D heatwaves)
    ds_labels = xr.open_dataset(join(write_directory,f"labels_cc3d.nc"),engine='netcdf4')
    da_labels = ds_labels.label

    # Load heatwaves indices database
    df_htws = pd.read_csv(join(write_directory,f"df_htws.csv"),header=0,index_col=0)

    df_htws['EM-DAT Total Deaths'] = [0]*len(df_htws)
    df_htws['EM-DAT heatwaves'] = [" "]*len(df_htws)

    # Link EM-DAT country names format to netCDF mask country names format
    country_dict = {'Albania':'Albania', 'Austria':'Austria', 'Belarus':'Belarus',
                    'Belgium':'Belgium', 'Bosnia and Herzegovina':'Bosnia_and_Herzegovina',
                    'Bulgaria':'Bulgaria', 'Croatia':'Croatia', 'Cyprus':'Cyprus', 
                    'Czechia':'Czechia', 'Denmark':'Denmark', 'Estonia':'Estonia', 
                    'Finland':'Finland', 'France':'France', 'Germany':'Germany', 'Greece':'Greece', 
                    'Hungary':'Hungary', 'Iceland':'Iceland', 'Ireland':'Ireland', 
                    'Italy':'Italy', 'Latvia':'Latvia', 'Lithuania':'Lithuania',
                    'Luxembourg':'Luxembourg', 'Montenegro':'Montenegro',
                    'North Macedonia':'Macedonia',
                    'Moldova':'Moldova', 'Netherlands (Kingdom of the)':'Netherlands', 'Norway':'Norway', 
                    'Poland':'Poland','Portugal':'Portugal', 'Romania':'Romania',
                    'Russian Federation':'Russia', 'Serbia':'Serbia', 
                    'Serbia Montenegro':'Serbia', # The corresponding heatwave happened in Serbia, cf 'Location' data of EM-DAT
                    'Slovakia':'Slovakia', 'Slovenia':'Slovenia', 'Spain':'Spain', 'Sweden':'Sweden',
                    'Switzerland':'Switzerland', 'Türkiye':'Turkey',
                    'United Kingdom of Great Britain and Northern Ireland':'United_Kingdom',
                    'Ukraine':'Ukraine','Yugoslavia':'Serbia'} # The corresponding heatwave happened in Serbia, cf 'Location' data of EM-DAT

    # Load cell area
    ds_cell_area = xr.open_dataset(join(read_directory,"ERA5","ERA5_Europe_cellarea.nc"),engine='netcdf4') # Area of each grid cell, in m²
    da_cell_area = ds_cell_area.cell_area/1e6 # Load DataArray and convert to km²

    for emdat_event in tqdm(df_emdat.index) :
        country = df_emdat.loc[emdat_event,'Country']
        mask_country = xr.open_dataset(join(read_directory,"ERA5","Mask",f"Mask_{country_dict[country]}_ERA5_0.25deg.nc"),engine='netcdf4').mask
        year_event = int(df_emdat.loc[emdat_event,'Start Year'])
        da_labels_event = da_labels.sel(time=(da_labels.time.dt.year==year_event)) # Keep only values matching the year of EM-DAT event

        if np.isnan(df_emdat.loc[emdat_event,'Start Day']) : # If Start Day is not recorded, start searching at 1st day of the month
            date_beg = pd.Timestamp(year_event, int(df_emdat.loc[emdat_event,'Start Month']), 1)
        else : # Otherwise, take recorded day and add flex_time_span to 
            date_beg = pd.Timestamp(year_event, int(df_emdat.loc[emdat_event,'Start Month']), int(df_emdat.loc[emdat_event,'Start Day'])) - pd.DateOffset(days=flex_time_span)

        if np.isnan(df_emdat.loc[emdat_event,'End Day']) : # If End Day is not recorded, keep searching until last day of the month
            date_end = pd.Timestamp(year_event, int(df_emdat.loc[emdat_event,'End Month']), 1)
            date_end = pd.offsets.MonthEnd().rollforward(date_end)
        else : # Otherwise, take recorded day and add flex_time_span
            date_end = pd.Timestamp(year_event, int(df_emdat.loc[emdat_event,'End Month']), int(df_emdat.loc[emdat_event,'End Day'])) + pd.DateOffset(days=flex_time_span)

        date_beg = cftime.DatetimeNoLeap(date_beg.year, date_beg.month, date_beg.day)
        date_end = cftime.DatetimeNoLeap(date_end.year, date_end.month, date_end.day)

        da_labels_event = da_labels.sel(time=slice(date_beg,date_end)) # Keep only values matching the period of EM-DAT event
        da_labels_event = da_labels_event*([mask_country==0]) # Keep only values matching the country of EM-DAT event, mask_country is 0 in the country, 1 elsewhere
        
        labels_list = np.unique(da_labels_event.data)
        labels_list = labels_list[np.where(labels_list!=0)] # Remove zero which corresponds to the absence of hot point, not a heatwave label ID
        summed_area_dict = {}
        total_summed_area = 0
        for label in labels_list :
            da_lab_htw = da_labels_event.where(da_labels_event==label,drop=False)
            da_lab_htw.data = (da_lab_htw.data>0)
            summed_area = (da_lab_htw*da_cell_area.data).sum().data
            summed_area_dict[label] = summed_area
            total_summed_area += summed_area

        for label in labels_list :
            htw_list = df_htws.loc[label,'EM-DAT heatwaves']
            df_htws.loc[label,'EM-DAT heatwaves'] = htw_list + ","*(len(htw_list)>1) + emdat_event
            if ~(df_emdat.loc[emdat_event,'Total Deaths'] is None or np.isnan(df_emdat.loc[emdat_event,'Total Deaths'])) :
                df_htws.loc[label,'EM-DAT Total Deaths'] += round(df_emdat.loc[emdat_event,'Total Deaths']*summed_area_dict[label]/total_summed_area)
        df_emdat.loc[emdat_event,'overlap CC3D'] = str(labels_list)
    # Save dataframe 
    df_htws.to_csv(join(write_directory,"df_htws_step2.csv"))
    df_emdat.to_csv(join(write_directory,"df_emdat_overlap.csv"))
    ds_labels.close()
    return

def r2_score_for_bootstrap(x,y):
    x = x.reshape(-1,1)
    model = LinearRegression().fit(x, y)
    return(model.score(x, y))

def correlation_for_bootstrap(x,y):
    correlation = stats.linregress(x, y, nan_policy='omit')
    return correlation.rvalue

def validate_indices_vs_emdat_impacts(read_directory,write_directory,emdat_file_path,pop_file_path,temp_variable='t2m',daily_var='tx',start_year=1975,end_year=2021,start_year_ref=1975,end_year_ref=2021,anomaly=True,nb_days=4,threshold_value=95,relative_threshold=True,flex_time_span=3,connectivity=26,dont_skip_figure=True):
    # Load heatwaves indices database
    df_htws = pd.read_csv(join(write_directory,f"df_htws_step2.csv"),header=0, index_col=0)
    # Load EM-DAT dataset
    df_emdat = pd.read_excel(join(emdat_file_path),header=0, index_col=0) 
    df_emdat = df_emdat[(df_emdat['Start Year']>=start_year) & (df_emdat['Start Year']<=end_year)] # Only keep events of the studied period (default 1975-2021)
    df_emdat = df_emdat[(df_emdat['Start Month'].isin([6,7,8])) | (df_emdat['End Month'].isin([6,7,8]))] # Remove heatwaves occurring outside JJA period

    land_mask = xr.open_dataset(join(read_directory,'ERA5','Mask',f"Mask_Europe_land_only_ERA5_0.25deg.nc"),engine='netcdf4').mask
    ds_pop = xr.open_dataset(join(pop_file_path), engine='netcdf4') # Population count, need to convert to density
    da_pop = ds_pop.Band1 # Data is in thousands of people
    pop_Europe = da_pop.where(land_mask==0).sum(dim=("lat","lon"))

    htw_indices = ['Intensity','Spatial extent','Duration','Max','Temp_sum','HWMId_sum','HWMId_mean','Exposed_population','Temp_pop','HWMId_pop','HWMId_pop_mean']#,'Total Exposed_population'
    shown_indices = ['Intensity', 'HWMId_sum','Exposed_population','HWMId_pop']

    df_scores = pd.DataFrame(index=htw_indices,columns=['R Pearson','p-value','significance','AUC ROC','Total score'],dtype=object)
    df_correlation = df_htws[df_htws['EM-DAT heatwaves'].map(len)>1]

    log_scale_dict = {'Intensity':False,'Spatial extent':True,'Duration':False,'Max':False,
    'Temp_sum':True,'HWMId_sum':True,'HWMId_mean':False,'Exposed_population':True,#'Total Exposed_population':True,
    'Temp_pop':True,'HWMId_pop':True,'HWMId_pop_mean':False}

    sns.set_theme(style="whitegrid")
    sns.set(font_scale=1.25)
    norm = clrs.LogNorm(vmin=1, vmax=df_htws['EM-DAT Total Deaths'].max())
    sm = plt.cm.ScalarMappable(cmap="YlOrRd", norm=norm)
    sm.set_array([])

    for index in tqdm(htw_indices):
        # Fill scores table
        if np.isnan(df_htws.loc[:,index]).any():
            df_scores.loc[index,'R2'] = str({'median':0,'q5':0,'q95':0})
        else:
            r2_bootstrap = stats.bootstrap(data=(df_correlation.loc[:,index],df_correlation.loc[:,'EM-DAT Total Deaths']),statistic=r2_score_for_bootstrap, n_resamples=1e4, paired=True, confidence_level=0.90, alternative='two-sided', method='percentile')
            df_scores.loc[index,'R2'] = str({'median':np.median(r2_bootstrap.bootstrap_distribution),'q5':r2_bootstrap.confidence_interval.low,'q95':r2_bootstrap.confidence_interval.high}) #np.round(correlation.rvalue,6)
        
        if np.isnan(df_htws.loc[:,index]).any():
            df_scores.loc[index,'AUC ROC'] = str({'median':0,'q5':0,'q95':0})
        else:
            roc_auc_bootstrap = stats.bootstrap(data=((df_htws['EM-DAT heatwaves'].map(len)>1).values,df_htws.loc[:,index]),statistic=metrics.roc_auc_score, n_resamples=1e4, paired=True, confidence_level=0.90, alternative='two-sided', method='percentile')
            df_scores.loc[index,'AUC ROC'] = str({'median':np.median(roc_auc_bootstrap.bootstrap_distribution),'q5':roc_auc_bootstrap.confidence_interval.low,'q95':roc_auc_bootstrap.confidence_interval.high}) #np.round(roc_auc,6)

        correlation_bootstrap = stats.bootstrap(data=(df_correlation.loc[:,index],df_correlation.loc[:,'EM-DAT Total Deaths']),statistic=correlation_for_bootstrap, n_resamples=1e4, paired=True, confidence_level=0.90, alternative='two-sided', method='percentile')
        df_scores.loc[index,'R Pearson'] = str({'median':np.median(correlation_bootstrap.bootstrap_distribution),'q5':correlation_bootstrap.confidence_interval.low,'q95':correlation_bootstrap.confidence_interval.high})#np.round(correlation.rvalue,6)

        
        df_scores.loc[index,'R2'] = str({'median':np.median(r2_bootstrap.bootstrap_distribution),'q5':r2_bootstrap.confidence_interval.low,'q95':r2_bootstrap.confidence_interval.high}) #np.round(correlation.rvalue,6)
        df_scores.loc[index,'Total score'] = np.round(np.median(r2_bootstrap.bootstrap_distribution)*np.median(roc_auc_bootstrap.bootstrap_distribution)*np.median(correlation_bootstrap.bootstrap_distribution),6)

        df_plot = df_htws[df_htws[index]>0]
        df_corrplot = df_correlation[df_correlation[index]>0]
        # Plot correlation between index values and impact values (only for overlapping heatwaves) and save correlation scores in df
        sns.set_theme(style="whitegrid")
        ax = sns.regplot(data=df_corrplot, x="EM-DAT Total Deaths", y=index, marker="x", line_kws=dict(color="r"))
        plt.savefig(join(write_directory,'figs',f"correlation_{index}.pdf"),dpi=1200)
        plt.close()

        # Plot distributions showing overlapping heatwaves
        sns.set_theme(style="whitegrid")
        ax = sns.displot(df_plot, x=index, kind="kde", log_scale=log_scale_dict[index], bw_adjust=0.5) # Plot KDE distribution
        sns.rugplot(df_plot,x=index) 
        sns.rugplot(df_corrplot,x=index,hue="EM-DAT Total Deaths",palette=sns.color_palette("YlOrBr", as_cmap=True),hue_norm=norm,height=0.1,legend=False)
        # Add colorbar for EM-DAT Total Deaths
        ax.figure.colorbar(sm,ax=ax.ax,label='Total Deaths')
        plt.savefig(join(write_directory,'figs',f"distrib_{index}.pdf"),dpi=1200)
        plt.close()
    df_scores.to_csv(join(write_directory,"df_scores.csv"))
    
    # Compute trends with Mann-Kendall trend test and Sen's slope (R script)
    df_htws["Exposed_population (relative)"] = None
    for i in df_htws.index : # Compute relative exposed population 
        df_htws.loc[i,"Exposed_population (relative)"] = np.round(100 * df_htws.loc[i,"Exposed_population"]/pop_Europe.sel(time=(pop_Europe.time.dt.year==df_htws.loc[i,"Year"])).data[0],4)
    
    # Export to csv for usage in R script
    df_htws.to_csv(join(write_directory,"df_htws_step3.csv"))

    # Compute and export dataframe for frequency trend
    df_frequency_export = pd.DataFrame(columns=["Heatwaves"],index=range(start_year,end_year+1),data=None)
    for i in df_frequency_export.index:
        try :
            df_frequency_export.loc[i,"Heatwaves"]=df_htws.groupby("Year").count().loc[i,"Start Date"]
        except: # Some years may be absent because there are no heatwaves
            df_frequency_export.loc[i,"Heatwaves"]=0
    # Export frequency to csv for usage in R script
    df_frequency_export.rename_axis(index="Year").to_csv(join(write_directory,"frequency.csv"))

    print("Calling mk_senslope_CI.R...")
    subprocess.call(f"Rscript /home/user/These/extreme_htws_cc3d/Code/Run_analysis/mk_senslope_CI.R {write_directory}", shell=True)
    print("Done")

    if dont_skip_figure :
        # Plot distribution and correlation for 4 shown_indices (4 subplots)
        # Correlations
        sns.set_theme(style="whitegrid")
        f, axs = plt.subplots(2, 2, figsize=(12, 6))
        for i in range(len(shown_indices)) :
            df_corrplot = df_correlation[df_correlation[shown_indices[i]]>0]
            sns.regplot(data=df_corrplot, x="EM-DAT Total Deaths", y=shown_indices[i], marker="x", line_kws=dict(color="r"), ax=axs[i//2,i%2])
        f.tight_layout()
        f.savefig(join(write_directory,'figs',"correlation_4idx.pdf"),dpi=1200)
        plt.close()
        
        # Distributions
        sns.set_theme(style="whitegrid")
        f, axs = plt.subplots(2, 2, figsize=(8, 6),sharey=True)
        for i in range(len(shown_indices)) :
            df_plot = df_htws[df_htws[shown_indices[i]]>0]
            df_corrplot = df_correlation[df_correlation[shown_indices[i]]>0]
            ax = sns.kdeplot(df_plot, x=shown_indices[i], log_scale=log_scale_dict[shown_indices[i]], bw_adjust=0.5, ax=axs[i//2,i%2]) # Plot KDE distribution
            sns.rugplot(df_plot,x=shown_indices[i], ax=axs[i//2,i%2],color='b') 
            sns.rugplot(df_corrplot,x=shown_indices[i],hue="EM-DAT Total Deaths",palette=sns.color_palette("YlOrBr", as_cmap=True),hue_norm=norm,height=0.25,lw=1,legend=False, ax=axs[i//2,i%2])
            plt.ylim([0,1.05])
            if i%2==1 :
                ax.set(ylabel=None)
        f.tight_layout()
        cbar = f.figure.colorbar(sm,ax=axs,label='Total Deaths')
        f.savefig(join(write_directory,'figs',"distrib_4idx.pdf"),dpi=1200)
        plt.close()

        df_htws["Overlap"] = None
        df_htws["edgecolors"] = None
        df_htws["linewidth"] = None
        drop_list = []
        for i in df_htws.index :
            df_htws.loc[i,"Overlap"] = len(df_htws.loc[i,"EM-DAT heatwaves"])>1
            if df_htws.loc[i,"HWMId_pop"]==0 :
                drop_list.append(i)
        df_htws_bbplot = df_htws.drop(drop_list) # Drop the heatwaves where HWMId_pop = 0 in order to plot with logartihmic color scale
        # Bubbleplot, only for HWMId_pop
        label_2003 = 159
        label_2010 = [220,222,225,226]
        handles, labels = sns.scatterplot(data=df_htws_bbplot, x="Year", y="Spatial extent",  size="Duration",
            sizes=(20, 200),color='black').get_legend_handles_labels()
        plt.close()
        norm = clrs.LogNorm(vmin=1, vmax=df_htws_bbplot['HWMId_pop'].max())
        sm = plt.cm.ScalarMappable(cmap="rocket_r", norm=norm)
        sm.set_array([])
        sns.set_theme(style="whitegrid")
        f = plt.figure(figsize=(8, 6))
        ax = sns.scatterplot(data=df_htws.drop(drop_list), x="Year", y="Spatial extent", hue="HWMId_pop", hue_norm=norm, size="Duration",
            sizes=(20, 200),palette=sns.color_palette("rocket_r", as_cmap=True),style="Overlap",markers=["o","v"])
        #ax = sns.scatterplot(data=df_htws_bbplot, x="Year", y="Spatial extent", hue="HWMId_pop", hue_norm=norm, size="Duration",
        #    sizes=(20, 200),palette=sns.color_palette("rocket_r", as_cmap=True))#,style="Overlap")
        #sns.scatterplot(data=df_correlation, x="Year", y="Spatial extent", hue="HWMId_pop", hue_norm=norm, size="Duration",
        #    sizes=(20, 200),palette=sns.color_palette("rocket_r", as_cmap=True),edgecolor='limegreen',linewidth=1)#,style="Overlap")
        ax.annotate('2003', 
                xy=(df_htws_bbplot.loc[label_2003,'Year'],df_htws_bbplot.loc[label_2003,'Spatial extent']), 
                xytext=(1995, 3.5e6),
                arrowprops=dict(facecolor='red', width=1.5,connectionstyle='arc3, rad=-0.3',alpha=0.5),
                fontsize=12,
                color='red', alpha = 0.7)
        ax.annotate('2010', 
                xy=(df_htws_bbplot.loc[label_2010[3],'Year'],df_htws_bbplot.loc[label_2010[3],'Spatial extent']), 
                xytext=(2013, 2.5e6),
                arrowprops=dict(facecolor='red', width=1.5,connectionstyle='arc3, rad=0.4',alpha=0.5),
                fontsize=12,
                color='red', alpha = 0.7)
        ax.annotate('', 
                xy=(df_htws_bbplot.loc[label_2010[2],'Year'],df_htws_bbplot.loc[label_2010[2],'Spatial extent']), 
                xytext=(2013, 2.5e6),
                arrowprops=dict(facecolor='red', width=1.5,connectionstyle='arc3, rad=0.4',alpha=0.5),
                fontsize=12,
                color='red', alpha = 0.7)
        ax.annotate('', 
                xy=(df_htws_bbplot.loc[label_2010[1],'Year'],df_htws_bbplot.loc[label_2010[1],'Spatial extent']), 
                xytext=(2013, 2.5e6),
                arrowprops=dict(facecolor='red', width=1.5,connectionstyle='arc3, rad=0.4',alpha=0.5),
                fontsize=12,
                color='red', alpha = 0.7)
        ax.annotate('', 
                xy=(df_htws_bbplot.loc[label_2010[0],'Year'],df_htws_bbplot.loc[label_2010[0],'Spatial extent']), 
                xytext=(2013, 2.5e6),
                arrowprops=dict(facecolor='red', width=1.5,connectionstyle='arc3, rad=0.4',alpha=0.5),
                fontsize=12,
                color='red', alpha = 0.7)

        ax.semilogy()
        plt.legend(handles, labels, loc='upper left',title='Duration')#loc=(0.1,0.625)
        plt.ylabel('Spatial extent (km²)')
        cbar = f.figure.colorbar(sm,ax=ax,label='HWMId_pop')
        plt.tight_layout()
        plt.savefig(join(write_directory,'figs',"bubble_plot_HWMId_pop.pdf"),dpi=1200)
        plt.close()

    df_best_scores = pd.read_csv(join(read_directory,f"summary_detection_overlap_sensitivity.csv"),header=0, index_col=0)

    df_scores = df_scores.loc[["Intensity","HWMId_sum","Exposed_population","HWMId_pop"]]

    get_index = df_best_scores[(df_best_scores['temp_variable']==temp_variable)&(df_best_scores['daily_var']==daily_var)&(df_best_scores['start_year']==start_year)&(df_best_scores['end_year']==end_year)&(df_best_scores['start_year_ref']==start_year_ref)&(df_best_scores['end_year_ref']==end_year_ref)&(df_best_scores['anomaly']==anomaly)&(df_best_scores['nb_days']==nb_days)&(df_best_scores['relative_threshold']==relative_threshold)&(df_best_scores['threshold_value']==threshold_value)&(df_best_scores['flex_time_span']==flex_time_span)&(df_best_scores['connectivity']==connectivity)].index.values[0]

    best_auc_roc_idx = "Intensity"
    best_pearson_R_idx = "Intensity"
    best_R2_idx	= "Intensity"
    total_score_best_index = "Intensity"

    # Record scores and best scores in sensitivity analysis table
    for index in ["Intensity","HWMId_sum","Exposed_population","HWMId_pop"] :
        pearson = eval(df_scores.loc[index,'R Pearson'])['median']
        df_best_scores.loc[get_index,f"pearson_R_{index}"] = pearson
        if pearson > eval(df_scores.loc[best_pearson_R_idx,'R Pearson'])['median']:
            best_pearson_R_idx = index

        roc_auc = eval(df_scores.loc[index,'AUC ROC'])['median']
        df_best_scores.loc[get_index,f"auc_roc_{index}"] = roc_auc
        if roc_auc > eval(df_scores.loc[best_auc_roc_idx,'AUC ROC'])['median']:
            best_auc_roc_idx = index

        R2 = eval(df_scores.loc[index,'R2'])['median']
        df_best_scores.loc[get_index,f"R2_{index}"] = R2
        if R2 > eval(df_scores.loc[best_R2_idx,'R2'])['median']:
            best_R2_idx = index
            
        total_score = df_scores.loc[index,'Total score']
        df_best_scores.loc[get_index,f"Total_score_{index}"] = total_score
        if total_score > df_scores.loc[total_score_best_index,'Total score']:
            total_score_best_index = index
    df_best_scores = df_best_scores.astype({"best_pearson_R_idx":str,"best_auc_roc_idx":str,"best_R2_idx":str,"total_score_best_index":str})
    df_best_scores.loc[get_index,"best_pearson_R_idx"] = best_pearson_R_idx
    df_best_scores.loc[get_index,"best_auc_roc_idx"] = best_auc_roc_idx
    df_best_scores.loc[get_index,"best_R2_idx"] = best_R2_idx
    df_best_scores.loc[get_index,"total_score_best_index"] = total_score_best_index

    # Load Mann-Kendall trend results
    try:
        with open(join(write_directory,"mk_result.json"), mode='r') as f:
            mk_result = json.load(f)

        df_mk_result = pd.json_normalize(mk_result)
        df_mk_result = df_mk_result.set_index('_row')

        
        for index in df_mk_result.index:
            df_best_scores = df_best_scores.astype({f"trend_{index}":str})
            df_best_scores.loc[get_index,f"trend_{index}"] = str({'median':df_mk_result.loc[index,'slope'],'q5':df_mk_result.loc[index,'LCL'],'q95':df_mk_result.loc[index,'UCL']})
    except: #In case unable to compute MK trend test, leave NaNs in sensitivity table
        pass

    df_best_scores.loc[get_index,"nb_detected_htws"] = len(df_htws)
    df_best_scores.loc[get_index,"nb_overlap_htws"] = len(df_correlation)

    # Record the number of matching EM-DAT heatwaves
    list_htws = []
    for htws in df_htws["EM-DAT heatwaves"] :
        htws = htws.replace(' ','')
        if len(htws)>0 :
            htws = htws.split(",")
            list_htws += htws
    list_htws = np.unique(list_htws)
    df_best_scores.loc[get_index,"nb_emdat_matching_htws"] = len(list_htws)
    df_best_scores = df_best_scores.astype({"emdat_matching_htws":str})
    df_best_scores.loc[get_index,"emdat_matching_htws"] = str(list_htws)

    df_best_scores.to_csv(join(read_directory,"summary_detection_overlap_sensitivity.csv"))

    return


def analysis_top_detected_events(read_directory,write_directory) :
    '''This function is used to search for the top detected heatwaves in an alternate impact database.'''
    
    # Link database country names format to netCDF mask country names format
    country_dict = {'Albania':'Albania', 'Austria':'Austria', 'Belarus':'Belarus',
                    'Belgium':'Belgium', 'Bosnia and Herzegovina':'Bosnia_and_Herzegovina',
                    'Bulgaria':'Bulgaria', 'Croatia':'Croatia', 'Cyprus':'Cyprus', 
                    'Czechia':'Czechia', 'Denmark':'Denmark', 'Estonia':'Estonia', 
                    'Finland':'Finland', 'France':'France', 'Germany':'Germany', 'Greece':'Greece', 
                    'Hungary':'Hungary', 'Iceland':'Iceland', 'Ireland':'Ireland', 
                    'Italy':'Italy', 'Latvia':'Latvia', 'Lithuania':'Lithuania',
                    'Luxembourg':'Luxembourg', 'Montenegro':'Montenegro',
                    'North Macedonia':'Macedonia',
                    'Moldova':'Moldova', 'Netherlands (Kingdom of the)':'Netherlands', 'Norway':'Norway', 
                    'Poland':'Poland','Portugal':'Portugal', 'Romania':'Romania',
                    'Russian Federation':'Russia', 'Serbia':'Serbia', 
                    'Serbia Montenegro':'Serbia', # The corresponding heatwave happened in Serbia, cf 'Location' data of EM-DAT
                    'Slovakia':'Slovakia', 'Slovenia':'Slovenia', 'Spain':'Spain', 'Sweden':'Sweden',
                    'Switzerland':'Switzerland', 'Türkiye':'Turkey', 'Turkey':'Turkey',
                    'United Kingdom of Great Britain and Northern Ireland':'United_Kingdom',
                    'Ukraine':'Ukraine','Yugoslavia':'Serbia', # The corresponding heatwave happened in Serbia, cf 'Location' data of EM-DAT
                    'England':'United_Kingdom','England and Wales':'United_Kingdom','Czech Republic':'Czechia'} 
    
    Hammond_to_mask_dict = {'Czech Republic':'Czechia','Turkey':'Turkey','England':'United_Kingdom','England and Wales':'United_Kingdom'} # Dictionary to convert Hammond country names to country names compatible with mask files
    
    # Load temperature-related data in case indices have to be recomputed (only occur if one EM-DAT heatwave matches several CC3D heatwaves)
    ds_labels = xr.open_dataset(join(write_directory,f"labels_cc3d.nc"),engine='netcdf4')
    da_labels = ds_labels.label

    # Load cell area
    ds_cell_area = xr.open_dataset(join(read_directory,"ERA5","ERA5_Europe_cellarea.nc"),engine='netcdf4') # Area of each grid cell, in m²
    da_cell_area = ds_cell_area.cell_area/1e6 # Load DataArray and convert to km² 

    # Load heatwaves indices database
    df_htws = pd.read_csv(join(write_directory,f"df_htws_step3.csv"),header=0, index_col=0)
    df_scores = pd.read_csv(join(write_directory,"df_scores.csv"),header=0, index_col=0)

    df_htws["Affected countries >10%"] = [" "]*len(df_htws)
    df_htws["Affected countries <10%"] = [" "]*len(df_htws)
    df_htws["Hammond Matching"] = [" "]*len(df_htws)
    df_htws["Hammond Countries"] = [" "]*len(df_htws)
    df_htws["Hammond Deaths"] = [0]*len(df_htws)
    df_htws["Other Deaths"] = None
    df_htws["Other refs"] = None

    best_scoring_index = df_scores['Total score'].idxmax()

    df_Hammond = pd.read_excel(join(read_directory,"EM-DAT","Lucy_Hammond_ETE_data_V2.xlsx"), header=0, index_col=0)
    df_Hammond = df_Hammond[df_Hammond['Country'].isin(country_dict.keys())]
    for idx in df_Hammond.index :
        if df_Hammond.loc[idx,'Country'] in Hammond_to_mask_dict.keys():
            df_Hammond.loc[idx,'Country'] = Hammond_to_mask_dict[df_Hammond.loc[idx,'Country']]
    
    for htw_label in tqdm(df_htws.index) :
        if len(df_htws.loc[htw_label,"EM-DAT heatwaves"])<2 : # Only check heatwaves that do not overlap with EM-DAT; default string is " " so if len()<2, it means the cell is empty
            year_event = int(df_htws.loc[htw_label,'Year'])
            da_labels_event = da_labels.sel(time=(da_labels.time.dt.year==year_event)) # Keep only values matching the year of EM-DAT event
            for country in np.unique(list(country_dict.values())) : # For each country, check if it is affected by the considered heatwave
                flag_country_for_search = False
                mask_country = xr.open_dataset(join(read_directory,"ERA5","Mask",f"Mask_{country}_ERA5_0.25deg.nc"),engine='netcdf4').mask
                mask_country['lat'] = da_labels.lat # Fix broken latitude
                mask_country['lon'] = da_labels.lon # Fix broken longitude
                label = da_labels_event.where(da_labels_event==htw_label,drop=False)
                label.data = ~np.isnan(label.data) # Convert label to bool, True for heatwave points, False elsewhere
                area_label = (label.where(mask_country==0,drop=False)*da_cell_area.data).sum()
                area_country = da_cell_area.where(mask_country==0,drop=True).sum()
                if area_label > area_country/10 :
                    country_list = df_htws.loc[htw_label,"Affected countries >10%"]
                    df_htws.loc[htw_label,"Affected countries >10%"] = country_list + ","*(len(country_list)>1) + country
                    flag_country_for_search = True
                elif area_label > 0 :
                    country_list = df_htws.loc[htw_label,"Affected countries <10%"]
                    df_htws.loc[htw_label,"Affected countries <10%"] = country_list + ","*(len(country_list)>1) + country
                    flag_country_for_search = True
                area_label = label.where(mask_country==0,drop=False)
                area_label = area_label.where(da_labels_event==htw_label,drop=True)
                start_date = area_label.time.data[0]
                end_date = area_label.time.data[-1]
                sub_df_Hammond = df_Hammond[(df_Hammond["Country"]==country)]
                if len(sub_df_Hammond)>0 and flag_country_for_search :
                    for idx in sub_df_Hammond.index :
                        start_date_hammmond = [int(i) for i in datetime.strftime(df_Hammond.loc[idx,"Start date"],format="%Y/%m/%d").split("/")] # Get Hammond Start date and convert it to [YYYY,MM,DD]
                        start_date_hammmond = cftime.datetime(*start_date_hammmond,calendar='noleap') # Convert it to cftime datetim to compare with xarray.DataArray time
                        end_date_hammond = [int(i) for i in datetime.strftime(df_Hammond.loc[idx,"End date"],format="%Y/%m/%d").split("/")] # Get Hammond End date and convert it to [YYYY,MM,DD]
                        end_date_hammond = cftime.datetime(*end_date_hammond,calendar='noleap') # Convert it to cftime datetime to compare with xarray.DataArray time

                        if (start_date>=start_date_hammmond and start_date<=end_date_hammond) or (end_date>=start_date_hammmond and end_date<=end_date_hammond) or (start_date<=start_date_hammmond and end_date>=end_date_hammond) : # Check if dates match
                            idx_list = df_htws.loc[htw_label,"Hammond Matching"]
                            df_htws.loc[htw_label,"Hammond Matching"] = idx_list + ","*(len(idx_list)>1) + str(idx)
                            country_list = df_htws.loc[htw_label,"Hammond Countries"]
                            df_htws.loc[htw_label,"Hammond Countries"] = country_list + ","*(len(country_list)>1) + country
                            deaths = df_htws.loc[htw_label,"Hammond Deaths"]
                            df_htws.loc[htw_label,"Hammond Deaths"] = deaths + df_Hammond.loc[idx,"Deaths"]
        
    df_htws = df_htws.sort_values(by=best_scoring_index,ascending=False)        
    df_htws.to_csv(join(write_directory,"df_htws_top_events_NOCHANGE.csv"))

    return