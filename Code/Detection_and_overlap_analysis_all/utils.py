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

def compute_climatology_smooth(read_directory,write_directory,temp_file,start_year_ref=1975,end_year_ref=2021,temp_variable='t2m',daily_var='tx',smooth_span=15) :
    '''This function computes a climatology for each calendar day of the year. The seasonal cycle is then smoothed with a 31-day window. 
    By default, the climatology is computed over 1975-2021.
    This function can be used with several models and variables.'''

    # Load dataset
    ds = xr.open_dataset(join(read_directory,'ERA5',temp_variable,temp_file), engine='netcdf4')
    da = getattr(ds, temp_variable)
    # Keep only years of interest
    da.sel(time=(da.time.dt.year>=start_year_ref) & (da.time.dt.year<=end_year_ref))

    # Drop Feb 29
    da = da.convert_calendar("noleap")
    # Group using dayofyear and sum to compute mean at the end
    climatology = da.groupby(da.time.dt.dayofyear).mean(dim="time")

    # Smoothing
    # Handle easily first and last day of year by taking thrice the entire dataset and working on the middle one (avoiding border effects)
    extended_temp=np.zeros((365*3,np.shape(climatology.data)[1],np.shape(climatology.data)[2]))
    extended_temp[0:365,:,:]=climatology.data[:,:,:]
    extended_temp[365:730,:,:]=climatology.data[:,:,:]
    extended_temp[730:,:,:]=climatology.data[:,:,:]

    print("Smoothing...")
    #For each day d, compute the mean of the values between day d-smooth_span and d+smooth_span to have a physically realistic seasonal cycle
    for i in tqdm(range(365,730)):
        climatology.data[i-365,:,:] = np.nanmean(extended_temp[i-smooth_span:i+smooth_span+1,:,:],axis=0)

    # Export data to netcdf file
    climatology.to_netcdf(join(write_directory,f"{temp_variable}_{daily_var}_climatology_{start_year_ref}_{end_year_ref}.nc"))

    # Clear resources
    da.close()
    climatology.close()
    ds.close()
    return


def compute_distrib_percentile(read_directory,write_directory,temp_file,start_year=1975,end_year=2021,start_year_ref=1975,end_year_ref=2021,temp_variable='t2m',daily_var='tx',threshold_value=95,distrib_window_size=15,anomaly=False) :
    '''This function computes, for every calendar day, the n-th (n is the threshold_value, default 95) percentile of the corresponding distribution of daily variable. 
    By default, the distribution is computed over 1975-2021.'''

    if distrib_window_size%2==0:
        raise ValueError('distrib_window_size is even. It has to be odd so the window can be centered on the computed day.')

    ds = xr.open_dataset(join(read_directory,'ERA5',temp_variable,temp_file), engine="netcdf4")

    # Load climatology file to create output data structure and compute anomaly
    climatology = xr.open_dataarray(join(write_directory,f"{temp_variable}_{daily_var}_climatology_{start_year_ref}_{end_year_ref}.nc"), engine='netcdf4')
    
    # Initialize data array with the first file
    da = getattr(ds, temp_variable) # Iterate over files, except first one which has already been used in initialization
    # Drop 29 Feb and only keep the reference period
    da = da.convert_calendar("noleap")
    da = da.sel(time=(da.time.dt.year>=start_year_ref) & (da.time.dt.year<=end_year_ref))

    # Create threshold table by copying climatology table, values will be updated later
    threshold = climatology.copy()

    if anomaly :
        for year in tqdm(range(len(da.time)//365)) : # Iterate over the number of years
            da[year*365:(year+1)*365,:,:] = da[year*365:(year+1)*365,:,:] - climatology.data # Compute anomaly
    else :
        climatology.close()
    
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
    threshold.to_netcdf(join(write_directory,f"{temp_variable}_{daily_var}_distrib_{'ano_'*anomaly}{start_year_ref}_{end_year_ref}_threshold_{threshold_value}_window_{distrib_window_size}d.nc"))

    # Clear resources
    da.close()
    threshold.close()
    if anomaly :
        climatology.close()
    ds.close()
    return

def cc3d_scan_heatwaves(read_directory,write_directory,temp_file,start_year=1975,end_year=2021,start_year_ref=1975,end_year_ref=2021,temp_variable='t2m',daily_var='tx',threshold_value=95,relative_threshold=True,distrib_window_size=15,anomaly=False,nb_days=4,dust_threshold=775) :
    '''This function carries out a cc3d scan (https://pypi.org/project/connected-components-3d/) to detect heatwaves in the meteorological database (default ERA5, t2m, tx).
    The heatwaves point are labeled with a number corresponding to a heatwave identifier.
    Otherwise, values are set to -9999.'''
    # Set dust threshold to supress small amount of points. The threshold of 775 have been established by testing with different values and set to the point that the number of heatwaves did not change by more than 1% between two increment of 25.

    
    ds = xr.open_dataset(join(read_directory,'ERA5',temp_variable,temp_file), engine="netcdf4")

    if anomaly :
        climatology = xr.open_dataarray(join(write_directory,f"{temp_variable}_{daily_var}_climatology_{start_year_ref}_{end_year_ref}.nc"), engine='netcdf4')
        # Keep only JJA values
        mask = (climatology.dayofyear>=152) & (climatology.dayofyear<=243) # dayofyear ranges from 1 to 365 ; 152 is June 1st, 243 is August 31st
        climatology = climatology.sel(dayofyear=mask)

    if relative_threshold : # Load temperature threshold for reference period :
        threshold = xr.open_dataarray(join(write_directory,f"{temp_variable}_{daily_var}_distrib_{'ano_'*anomaly}{start_year_ref}_{end_year_ref}_threshold_{threshold_value}_window_{distrib_window_size}d.nc"), engine='netcdf4')
        # Keep only JJA values
        mask = (threshold.dayofyear>=152) & (threshold.dayofyear<=243) 
        threshold = threshold.sel(dayofyear=mask)
    else : # If absolute threshold, only need a scalar, not a 3D array
        threshold = threshold_value

    connectivity = 26 # only 4,8 (2D) and 26, 18, and 6 (3D) are allowed
    N_labels = 0 # Count the numbers of patterns
    da = getattr(ds, temp_variable)
    # Drop 29 Feb and correct day of year
    da = da.convert_calendar("noleap")
    # Keep only JJA values and years of interest
    da.sel(time=(da.time.dt.year>=start_year) & (da.time.dt.year<=end_year))
    mask = (da.time.dt.season=='JJA')
    da = da.sel(time=mask)

    # Initialize output array 
    label = da.copy()

    print("Computing cc3d.connected_components labels and dusting...")

    for year in tqdm(range(end_year-start_year+1)) :# Iterate over the years
        da_year = da[year*92:(year+1)*92,:,:]# Select data for the given year
        if anomaly : # Substract climatology to compute anomaly
            da_year = da_year - climatology.data
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
        label.data[year*92:(year+1)*92,:,:] = remove_sea_heatwaves(read_directory="TODO",labels=label[year*92:(year+1)*92,:,:])
        
        # Update N_labels
        N_labels += N_added

    # Set DataArray to Dataset to set variable name
    label = label.to_dataset(name="label")
    # Save to netCDF 
    label.to_netcdf(join(write_directory,f"{temp_variable}_{daily_var}_labels_{'ano_'*anomaly}years_{start_year}_{end_year}_ref_{start_year_ref}_{end_year_ref}_threshold_{threshold_value}_{nb_days}d_window_{distrib_window_size}d.nc"))
    
    # Clear resources
    da.close()
    ds.close()
    label.close()
    threshold.close()
    if anomaly :
        climatology.close()
    print(N_labels,"heatwaves detected")
    return

def remove_outside_heatwaves(read_directory,labels,land_area_fraction_threshold=0.5,dust_threshold=775) :
    '''Remove sea heatwaves and heatwaves occurring outside continental Europe
    '''
    
    land_mask = xr.open_dataset(join(read_directory,'ERA5','Mask',f"Mask_Europe_land_only_ERA5_0.25deg.nc"),engine='netcdf4').mask

    # Remove sea points, North Africa points, and Middle East points that are South of Turkey
    labels = labels * (land_mask.data==0)

    # Remove small heatwaves with the cc3d dust function
    connectivity = 26 # only 4,8 (2D) and 26, 18, and 6 (3D) are allowed
    labels = labels * cc3d.dust((labels>0),dust_threshold)
    land_area_fraction.close()
    return labels

def compute_Russo_HWMId(read_directory,write_directory,temp_file,start_year=1975,end_year=2021,start_year_ref=1975,end_year_ref=2021,temp_variable='t2m',daily_var='tx',distrib_window_size=15,anomaly=False) :
    """Compute the Russo_HWMId index map.
    Based on HWMId defined by Russo et al (2015, https://dx.doi.org/10.1088/1748-9326/10/12/124003 )."""

    temp_file = f"ERA5_{temp_variable}_{daily_var}_Europe_day_0.25deg_{start_year}-{end_year}.nc"
    ds = xr.open_dataset(join(read_directory,temp_file), engine='netcdf4')
    da = getattr(ds, temp_variable)
    # Drop Feb 29 and only keep JJA days
    da = da.convert_calendar("noleap")
    da = da.sel(time=(da.time.dt.season=='JJA'))

    if anomaly :
        # Keep only JJA values
        climatology = xr.open_dataarray(join(write_directory,f"{temp_variable}_{daily_var}_climatology_{start_year_ref}_{end_year_ref}.nc"), engine='netcdf4')
        mask = (climatology.dayofyear>=152) & (climatology.dayofyear<=243) # dayofyear ranges from 1 to 365 ; 152 is June 1st, 243 is August 31st
        climatology = climatology.sel(dayofyear=mask)
        da = da - climatology
        

    temp_25p = np.percentile(da.groupby(da.time.dt.year).max(), 25, axis=0)
    temp_75p = np.percentile(da.groupby(da.time.dt.year).max(), 75, axis=0)
    
    Russo_HWMId = da.copy()
    Russo_HWMId.data = np.maximum((da - temp_25p)/(temp_75p - temp_25p), 0)
    Russo_HWMId.to_dataset(name="HWMId")

    Russo_HWMId.to_netcdf(join(write_directory,f"{temp_variable}_{daily_var}_Russo_HWMId_{'ano_'*anomaly}years_{start_year}_{end_year}_ref_{start_year_ref}_{end_year_ref}.nc"))

    Russo_HWMId.close()
    da.close()
    ds.close()
    if anomaly :
        climatology.close()
    return


def create_heatwaves_indices_database(read_directory,write_directory,start_year=1975,end_year=2021,start_year_ref=1975,end_year_ref=2021,temp_variable='t2m',daily_var='tx',distrib_window_size=15,anomaly=False,nb_days=4,relative_threshold=True):
    '''This function is used to create the dataset of the indices of the detected heatwaves. The set of detected heatwaves depends on all the parameters.'''

    # Load temperature-related data
    ds_labels = xr.open_dataset(join(read_directory,f"{temp_variable}_{daily_var}_labels_{'ano_'*anomaly}years_{start_year}_{end_year}_ref_{start_year_ref}_{end_year_ref}_threshold_{threshold_value}_{nb_days}d_window_{distrib_window_size}d.nc"),enine='netcdf4')
    da_labels = ds_labels.label
    ds_temp = xr.open_dataset(join(read_directory,f"ERA5_{temp_variable}_{daily_var}_Europe_day_0.25deg_{start_year}-{end_year}.nc"), engine='netcdf4')
    da_temp = getattr(ds_temp,temp_variable)
    da_temp = da_temp.convert_calendar("noleap") # Drop 29 Feb
    ds_threshold = xr.open_dataset(join(read_directory,f"{temp_variable}_{daily_var}_distrib_{'ano_'*anomaly}{start_year_ref}_{end_year_ref}_threshold_{threshold_value}_window_{distrib_window_size}d.nc"),engine='netcdf4')
    da_threshold = ds_threshold.threshold
    ds_HWMId = xr.open_dataset(join(read_directory,f"{temp_variable}_{daily_var}_Russo_HWMId_{'ano_'*anomaly}years_{start_year}_{end_year}_ref_{start_year_ref}_{end_year_ref}.nc"),engine='netcdf4')
    da_HWMId = ds_HWMId.HWMId

    # Keep only JJA values; labels and HWMId are already JJA-only
    da_temp = da_temp.sel(time = (da_temp.da.time.dt.season=='JJA'))
    mask = (da_threshold.dayofyear>=152) & (da_threshold.dayofyear<=243) # dayofyear ranges from 1 to 365 ; 152 is June 1st, 243 is August 31st
    da_threshold = da_threshold.sel(dayofyear=mask)

    if anomaly :
        climatology = xr.open_dataarray(join(write_directory,f"{temp_variable}_{daily_var}_climatology_{start_year_ref}_{end_year_ref}.nc"), engine='netcdf4')
        # Keep only JJA values
        climatology = climatology.sel(dayofyear=mask)
        da_temp = da_temp - climatology.data

    # Load cell area
    ds_cell_area = xr.open_dataset(join(read_directory,"ERA5_Europe_cellarea.nc"),engine='netcdf4') # Area of each grid cell, in m²
    da_cell_area = ds.cell_area/1e6 # Load DataArray and convert to km²

    # Load population data
    ds_pop = xr.open_dataset(join(TODO_read_directory_pop,"GHS_POP_R2023A_1975_2030_ERA5_Europe_grid.nc"), engine='netcdf4') # Population count, need to convert to density
    da_pop = ds_pop.Band1

    da_pop_density = da_pop/da_cell_area # Population density in person/km²

    labels_list = np.unique(da_labels.data)
    labels_list = labels_list[np.where(labels_list!=0)] # Remove zero which corresponds to the absence of hot point, not a heatwave label ID
    df_htw = pd.DataFrame(columns=['Year','Start Date','End Date','Mean','Spatial extent','Duration','Max','Temp_sum','HWMId_sum','Affected population','Total affected pop','Temp_pop','HWMId_pop','EM-DAT heatwaves','EM-DAT Total Deaths'],index=labels_list,data=None)    

    for year in range(end_year-start_year+1) : # Compute temperature exceedance relative to the threshold, iterate over years (92 JJA days per year)
        da_temp[year*92,(year+1)*92,:,:] = da_temp[year*92,(year+1)*92,:,:] - da_threshold

    output_dir = os.path.join("Output",database,f"{datavar}_{daily_var}", TODO,
                            f"{database}_{datavar}_{daily_var}_{name_dict_anomaly[anomaly]}_JJA_{nb_days}days_before_scan_{year_beg}_{year_end}_{threshold_value}{name_dict_threshold[relative_threshold]}_{distrib_window_size}days_window_climatology_{year_beg_climatology}_{year_end_climatology}")
    
    # Compute weights for latitude-weighted mean
    weights = np.cos(np.deg2rad(da_temp.lat))
    weights.name = "weights"

    for label in tqdm(labels_list) : # Iterate on heatwaves
        da_bool_htw = da_labels.where(da_labels==label, drop=True).fillna(0)>0 # Select days and grid points for the heatwave of interest and convert to bool array
        da_temp_htw = da_temp.where(da_labels==label, drop=True)
        da_HWMId_htw = da_HWMId.where(da_labels==label, drop=True)

        labels_bool_2D = np.max(labels==label,axis=0) # Squeeze heatwave labels on a boolean 2D-map to see maximum spatial extension
        da_pop_htw = da_pop.where(labels_bool_2D,drop=True)
        da_pop_density_htw = da_pop_density.where(labels_bool_2D,drop=True)

        df_htw.loc[label,'Year'] = da_temp_htw.time.dt.year.data[0]
        df_htw.loc[label,'Start Date'] = da_temp_htw.time.dt.date.data[0]
        df_htw.loc[label,'End Date'] = da_temp_htw.time.dt.date.data[-1]
    
        df_htw.loc[label,'Mean'] = da_temp_htw.weighted(weights).mean().data
        df_htw.loc[label,'Spatial extent'] = (da_cell_area*labels_bool_2D).sum().data
        df_htw.loc[label,'Duration'] = len(da_temp_htw.time)
        df_htw.loc[label,'Max'] = da_temp_htw.max().data
        df_htw.loc[label,'Temp_sum'] = da_temp_htw.weighted(weights).sum().data
        df_htw.loc[label,'HWMId_sum'] = da_HWMId_htw.weighted(weights).sum().data
        df_htw.loc[label,'Affected population'] = da_pop_htw.sum().data
        df_htw.loc[label,'Total affected population'] = (da_bool_htw*da_pop_htw).sum().data
        df_htw.loc[label,'Temp_pop'] = (da_temp_htw*da_pop_density_htw).weighted(weights).sum().data
        df_htw.loc[label,'HWMId_pop'] = (da_HWMId_htw*da_pop_density_htw).weighted(weights).sum().data

        df_htw.loc[label,'EM-DAT heatwaves'] = []
        df_htw.loc[label,'EM-DAT Total Deaths'] = 0

    # Save dataframe 
    df_htw.to_csv(os.path.join(output_dir,f"df_htws_{'anomaly_'*anomaly}.csv"))
    
    ds_labels.close()
    ds_temp.close()
    ds_threshold.close()
    ds_cell_area.close()
    ds_pop.close()
    if anomaly :
        climatology.close()
    return


def analyze_emdat_overlap(read_directory,write_directory,start_year=1975,end_year=2021,start_year_ref=1975,end_year_ref=2021,temp_variable='t2m',daily_var='tx',threshold_value=95,relative_threshold=True,distrib_window_size=15,anomaly=False,nb_days=4,flex_time_span=3):
    '''This function is used to analyse the spatial and temporal overlap between EM-DAT heatwaves and the meteorological database heatwaves (default ERA5) detected with the CC3D scan.
    The detection threshold depends on the parameters used precedently, which is why all these parameters are required.
    This function can be used with several variables'''

    # Load EM-DAT dataset
    df_emdat = pd.read_excel(join(TODO,"EM-DAT","EMDAT_Europe_Turkey-1975-2021-heatwaves.xlsx"),header=0, index_col=0) 
    df_emdat = df_emdat[(df_emdat['Year']>=year_beg) & (df_emdat['Year']<=year_end)] # Only keep events of the studied period (default 1975-2021)
    df_emdat = df_emdat[(df_emdat['Start Month'].isin([6,7,8])) | (df_emdat['End Month'].isin([6,7,8]))] # Remove heatwaves occurring outside JJA period

    # Load temperature-related data in case indices need to be recomputed (only occur if one EM-DAT heatwave matches several CC3D heatwaves)
    ds_labels = xr.open_dataset(join(read_directory,f"{temp_variable}_{daily_var}_labels_{'ano_'*anomaly}years_{start_year}_{end_year}_ref_{start_year_ref}_{end_year_ref}_threshold_{threshold_value}_{nb_days}d_window_{distrib_window_size}d.nc"),enine='netcdf4')
    da_labels = ds_labels.label

    # Load heatwaves indices database
    df_htws = pd.read_csv(os.path.join(TODO,f"df_htws_{'anomaly_'*anomaly}days.csv"),header=0,index_col=0)

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

    for emdat_event in tqdm(df_emdat.index) :
        country = df_emdat.loc[emdat_event,'Country']
        mask_country = xr.open_dataset(os.path.join(TODO,"Mask",f"Mask_{country_dict[country]}_ERA5_0.25deg.nc"),engine=',netcdf4').mask
        year_event = int(df_emdat.loc[emdat_event,'Year'])
        da_labels_event = da_labels.sel(time=(da_labels.time.dt.year==year_event)) # Keep only values matching the year of EM-DAT event

        if np.isnan(df_emdat.loc[emdat_event,'Start Day']) : # If Start Day is not recorded, start searching at 1st day of the month
            date_beg = pd.Timestamp(year_event, int(df_emdat.loc[emdat_event,'Start Month']), 1)
        else : # Otherwise, take recorded day and add flex_time_span to 
            date_beg = pd.Timestamp(year_event, int(df_emdat.loc[emdat_event,'Start Month']), int(df_emdat.loc[emdat_event,'Start Day'])) - pd.DateOffset(days=flex_time_span)

        if np.isnan(df_emdat.loc[emdat_event,'End Day']) : # If End Day is not recorded, keep searching until last day of the month
            date_end = pd.Timestamp(year_event, int(df_emdat.loc[emdat_event,'End Month']), 1)
            date_end = pd.offsets.MonthEnd().rollforward(d)
        else : # Otherwise, take recorded day and add flex_time_span
            date_end = Timestamp(year_event, int(df_emdat.loc[emdat_event,'End Month']), int(df_emdat.loc[emdat_event,'End Day'])) + pd.DateOffset(days=flex_time_span)

        da_labels_event = da_labels.sel(time=slice(date_beg,date_end)) # Keep only values matching the period of EM-DAT event
        da_labels_event = da_labels_event*(mask_country==0) # Keep only values matching the country of EM-DAT event, mask_country is 0 in the country, 1 elsewhere
        
        labels_list = np.unique(da_labels_event.data)
        labels_list = labels_list[np.where(labels_list!=0)] # Remove zero which corresponds to the absence of hot point, not a heatwave label ID
        for label in np.unique(da_labels_event.data) :
            df_htws.loc[label,'EM-DAT heatwaves'].append(emdat_event)
            df_htws.loc[label,'EM-DAT Total Deaths'] += int(df_emdat.loc[emdat_event,'Total Deaths'])

    # Save dataframe 
    df_htw.to_csv(os.path.join(TODO,f"df_htws_{'anomaly_'*anomaly}_flex_time_{flex_time_span}days.csv"))

    ds_labels.close()
    return

def validate_indices_vs_emdat_impacts(read_directory,write_directory,start_year=1975,end_year=2021,start_year_ref=1975,end_year_ref=2021,temp_variable='t2m',daily_var='tx',threshold_value=95,relative_threshold=True,distrib_window_size=15,anomaly=False,nb_days=4,flex_time_span=3) :
    # Load heatwaves indices database
    df_htws = pd.read_csv(os.path.join(TODO,f"df_htws_{'anomaly_'*anomaly}_flex_time_{flex_time_span}days.csv"),header=0, index_col=0)
    
    # Load EM-DAT dataset
    df_emdat = pd.read_excel(join(TODO,"EM-DAT","EMDAT_Europe_Turkey-1975-2021-heatwaves.xlsx"),header=0, index_col=0) 
    df_emdat = df_emdat[(df_emdat['Year']>=year_beg) & (df_emdat['Year']<=year_end)] # Only keep events of the studied period (default 1975-2021)
    df_emdat = df_emdat[(df_emdat['Start Month'].isin([6,7,8])) | (df_emdat['End Month'].isin([6,7,8]))] # Remove heatwaves occurring outside JJA period

    htw_indices = ['Mean','Spatial extent','Duration','Max','Temp_sum','HWMId_sum','Affected population','Total affected pop','Temp_pop','HWMId_pop']
    shown_indices = ['Mean', 'HWMId_sum','Affected population','HWMId_pop']
    # 'Year','Start Date','End Date','EM-DAT heatwaves','EM-DAT Total Deaths'
    df_scores = pd.DataFrame(index=htw_indices,columns=['R Pearson','p-value','AUC ROC'])
    df_correlation = df_htws[df_htws['EM-DAT heatwaves'].map(len)>0]

    log_scale_dict = {'Mean':False,'Spatial extent':False,'Duration':False,'Max':False,
    'Temp_sum':True,'HWMId_sum':True,'Affected population':False,'Total affected pop':False,'Temp_pop':True,'HWMId_pop':True}

    figure_directory = "TODO"

    sns.set_theme()
    norm = clrs.LogNorm(vmin=1, vmax=df_htws['EM-DAT Total Deaths'].max())
    sm = plt.cm.ScalarMappable(cmap="YlOrRd", norm=norm)
    sm.set_array([])

    for index in htw_indices :
        # Fill scores table
        correlation = stats.linregress(df_correlation.loc[:,index], df_correlation.loc[:,'EM-DAT Total Deaths'], nan_policy='omit')
        significance = ''+'*'*(correlation.pvalue<0.05)+'*'*(correlation.pvalue<0.01)+'*'*(correlation.pvalue<0.001)
        roc_auc = metrics.roc_auc_score((df_htws['EM-DAT heatwaves'].map(len)>0).values,df_htws.loc[:,index])
        df_scores.loc[index,(chosen_impact,'r_pearson')] = str(np.round(correlation.rvalue,6))+significance
        df_scores.loc[index,(chosen_impact,'roc_auc')] = np.round(roc_auc,6)

        # Plot correlation between index values and impact values (only for overlapping heatwaves) and save correlation scores in df
        ax = sns.regplot(data=df_emdat, x="EM-DAT Total Deaths", y=index, marker="x", line_kws=dict(color="r"))
        ax.savefig(join(figure_directory,f"correlation_{index}_{'anomaly_'*anomaly}_flex_time_{flex_time_span}days.pdf"),dpi=1200)
        plt.close()

        # Plot trends
        ax = None
        ax.savefig(join(figure_directory,f"trend_{index}_{'anomaly_'*anomaly}_flex_time_{flex_time_span}days.pdf"),dpi=1200)

        # Plot distributions showing overlapping heatwaves
        ax = sns.displot(df_emdat, x=index, kind="kde", log_scale=log_scale_dict[shown_indices[i]], bw_adjust=0.5) # Plot KDE distribution
        sns.rugplot(df_emdat,x=index) 
        sns.rugplot(df_correlation,x=index,hue="EM-DAT Total Deaths",palette=sns.color_palette("YlOrBr", as_cmap=True),hue_norm=norm,height=0.1,legend=False)
        # Add colorbar for EM-DAT Total Deaths
        ax.figure.colorbar(sm,ax=ax.ax,label='Total Deaths')
        ax.savefig(join(figure_directory,f"distrib_{index}_{'anomaly_'*anomaly}_flex_time_{flex_time_span}days.pdf"),dpi=1200)
        plt.close()

    # Plot trends, distribution, correlation for 4 shown_indices (4 subplots)
    # Correlations
    f, axs = plt.subplots(2, 2, figsize=(8, 4))
    for i in range(len(shown_indices)) :
        sns.regplot(data=df_correlation, x="EM-DAT Total Deaths", y=shown_indices[i], marker="x", line_kws=dict(color="r"), ax=axs[i//2,i%2])
    f.tight_layout()
    f.savefig(join(figure_directory,f"correlation_4idx_{'anomaly_'*anomaly}_flex_time_{flex_time_span}days.pdf"),dpi=1200)
    plt.close()
    # Trends
    f, axs = plt.subplots(2, 2, figsize=(8, 4))
    for i in range(len(shown_indices)) :
        pass
    f.tight_layout()
    f.savefig(join(figure_directory,f"trend_4idx_{'anomaly_'*anomaly}_flex_time_{flex_time_span}days.pdf"),dpi=1200)
    plt.close()
    # Distributions
    f, axs = plt.subplots(2, 2, figsize=(8, 4))
    for i in range(len(shown_indices)) :
        ax = sns.kdeplot(df_emdat, x=index, log_scale=log_scale_dict[shown_indices[i]], bw_adjust=0.5, ax=axs[i//2,i%2]) # Plot KDE distribution
        sns.rugplot(df_emdat,x=index, ax=axs[i//2,i%2],color='b') 
        sns.rugplot(df_correlation,x=index,hue="EM-DAT Total Deaths",palette=sns.color_palette("YlOrBr", as_cmap=True),hue_norm=norm,height=0.1,legend=False, ax=axs[i//2,i%2])
        if i%2==1 :
            ax.set(ylabel=None)
    f.tight_layout()
    f.figure.colorbar(sm,ax=axs,label='Total Deaths')
    f.savefig(join(figure_directory,f"distrib_{index}_{'anomaly_'*anomaly}_flex_time_{flex_time_span}days.pdf"),dpi=1200)
    plt.close()

    return