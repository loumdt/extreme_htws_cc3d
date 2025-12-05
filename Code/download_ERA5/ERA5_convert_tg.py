#!/usr/bin/env python3

#############################
## Yoann Robin             ##
## yoann.robin.k@gmail.com ##
#############################

import sys,os
import datetime as dt

import numpy  as np
import pandas as pd
import xarray as xr


if __name__ == "__main__":
	
	## Parameters
	pin  = os.path.join( "/data/tmandonnet" , "ERA5" , "t2m" )
	pout = os.path.join( "/data/tmandonnet" , "ERA5" , "t2m", "day" )
	
	## List files
	lfiles = os.listdir(pin)
	lfiles.sort()
	
	## Loop
	for f in lfiles:
		
		## Load hourly data
		idata = xr.open_dataset( os.path.join( pin , f ) )
		
		## Time axis
		year  = f.split("_")[-1][:4]
		timeh = idata.time
		timed = pd.date_range( f"{year}-01-01" , f"{year}-12-31" , freq = "d" )
		
		## Reformat the grid
		lat   = idata.latitude.values
		lon   = idata.longitude.values
		ilat  = np.argsort(lat)
		ilon  = np.argsort(lon)
		
		lat = lat[ilat]
		lon = lon[ilon]
		
		## Reformat data
		t2mh = xr.DataArray( idata.t2m.values[:,ilat,:][:,:,ilon] , dims = ["time","lat","lon"] , coords = [timeh,lat,lon] )
		t2md = t2mh.groupby("time.dayofyear").mean().rename( dayofyear = "time" ).assign_coords( time = timed )
		
		## Output
		odatad = xr.Dataset( { "t2m" : t2md } )
		
		## Add attributes
		odatad["time"].attrs["standard_name"] = "time"
		odatad["time"].attrs["long_name"]     = "Time axis"
		odatad["time"].attrs["axis"]          = "T"
		
		odatad["lat"].attrs["standard_name"] = "latitude"
		odatad["lat"].attrs["long_name"]     = "latitude coordinate"
		odatad["lat"].attrs["units"]         = "degrees_north"
		odatad["lat"].attrs["axis"]          = "Y"
		
		odatad["lon"].attrs["standard_name"] = "longitude"
		odatad["lon"].attrs["long_name"]     = "longitude coordinate"
		odatad["lon"].attrs["units"]         = "degrees_east"
		odatad["lon"].attrs["axis"]          = "X"
		
		odatad["t2m"].attrs["coordinates"]   = "lat lon"
		odatad["t2m"].attrs["standard_name"] = "2m_temperature"
		odatad["t2m"].attrs["long_name"]     = "t2m"
		odatad["t2m"].attrs["units"]         = "K"
		
		product = 'reanalysis-era5-single-levels'
		if year < "1959":
			product = 'reanalysis-era5-single-levels-preliminary-back-extension'
		
		odatad.attrs["title"]         = "ERA5"
		odatad.attrs["Conventions"]   = "CF-1.6"
		odatad.attrs["source"]        = "Climate Data Store"
		odatad.attrs["product"]       = product
		odatad.attrs["creation_date"] = str(dt.datetime.utcnow())[:19] + " (UTC)"
		
		## Encoding
		ny = lat.size
		nx = lon.size
		encodingd = { "time" : { "dtype" : "float32" , "zlib" : True , "complevel": 5 , "chunksizes" : (1,) , "units" : "days since 1850-01-01 00:00:00" } ,
					  "lon"  : { "dtype" : "float32" , "zlib" : True , "complevel": 5 , "chunksizes" : (nx,) } ,
					  "lat"  : { "dtype" : "float32" , "zlib" : True , "complevel": 5 , "chunksizes" : (ny,) } ,
					  "t2m"  : { "dtype" : "float32" , "zlib" : True , "complevel": 5 , "chunksizes" : (1,ny,nx) } }
		
		## And now in netcdf
		poutd = os.path.join( pout )
		if not os.path.isdir(poutd):
			os.makedirs(poutd)
		foutd = f"ERA5_Europe_day_t2m_{year}0101-{year}1231.nc"
		odatad.to_netcdf( os.path.join( poutd , foutd ) , encoding = encodingd )
	
		del idata
		del odatad
		del t2md
	
	print("Done")