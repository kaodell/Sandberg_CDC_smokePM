#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mk_cnty_rastor_kPM.py
    python script to make a rastor of county-assignments for each 
    WRF-chem grid cell for the krigged PM2.5 data 
    for John Sanberg's CDC projectv
    this code is adapted from the "average_var_bycounty.py" code 
    written by me for the GeoXO project
Created on Mon Jan 17 21:38:05 2022
@author: Katelyn O'Dell
"""
#%% user inputs
proj_folder =  '/Volumes/GoogleDrive/My Drive/Ongoing Projects/Sandberg_CDC/'

# shapefile with counties to average by
shp_fn = 'cb_2018_us_county_500k/cb_2018_us_county_500k'

# kPM netCDFfile with grid lat/lons
kPM_fn = 'krigedPM25_2018_v2.nc'

# kPM csv file from Bonne, to check this matches, it does!
# kPM_csv_fn = 'smokePM_2018.csv'

#%% load modules
from netCDF4 import Dataset
import pandas as pd
import shapefile
import matplotlib as mplt
import numpy as np
import datetime as dt
from ODell_udf_CDCprj import plt_map

#%% load data
shps_file = shapefile.Reader(proj_folder+shp_fn)
shps_shp = shps_file.shape(0)
shps_records = shps_file.records()
shps_shapes = shps_file.shapes()

nc_fid = Dataset(proj_folder + kPM_fn)
glat = nc_fid['lat'][:]
glon = nc_fid['lon'][:]
bk_pm = nc_fid['Background_PM25'][:]
tot_pm = nc_fid['PM25'][:]
hms = nc_fid['HMS_Smoke'][:]
nc_fid.close()

# kPM_bonne = pd.read_csv(proj_folder + kPM_csv_fn,header=2)

#%% calculate smoke PM to compare to csv version
smokePM = hms*(tot_pm-bk_pm)

#%% loop through counties and make grid
si = 0
fips = []
names = []
# check if grid cells are being assigned to two counties
grid_count = np.zeros(glon.shape)
# area mask
area_mask_full = np.zeros(glon.shape)
area_mask_full[:,:] = -99999
# loop through counties
for j in range(len(shps_records)):
        name = shps_records[j][5]
        fip = shps_records[j][4]
        fips.append(fip)
        names.append(name)
        print(100.0*(j/len(shps_records)))
        area_mask = np.zeros(glon.shape)
        area_shp = shps_shapes[j]
        # find grid cells in county
        for i in range(len(area_shp.parts)):
            i0 = area_shp.parts[i]
            if i < len(area_shp.parts)-1:
            		i1 = area_shp.parts[i+1] - 1
            else:
            		i1 = len(area_shp.points)
            seg = area_shp.points[i0:i1+1]
            mpath = mplt.path.Path(seg)
            points = np.array((glon.flatten(), glat.flatten())).T
            mask = mpath.contains_points(points).reshape(glon.shape)
            area_inds = np.where(mask==True)
            area_mask[area_inds] = 1
        full_area_inds = np.where(area_mask==1)
        area_mask_full[full_area_inds[0],full_area_inds[1]] = int(fip)
        grid_count[full_area_inds[0],full_area_inds[1]]+=1
        # looks like only 3 grid cells are close enough to be double-assigned! sweet.
        
fips = np.array(fips)
names = np.array(names)

#%% save to csv - use geopandas? see what method Bonne uses
# flatten lons and lats
lon_csv = glon.flatten()
lat_csv = glat.flatten()
fips_csv = area_mask_full.flatten()

# add to dataframe
smoke_df = pd.DataFrame({'lat':lat_csv,'lon':lon_csv,'fips':fips_csv})

# create datetime for 2018
a = dt.date(2018,1,1)
numdays = 365
day_strs = []
for x in range (0, smokePM.shape[0]):
    day_dt = a + dt.timedelta(days = x)
    day_strs.append(day_dt.strftime("%Y%m%d"))

# loop through days and add to csv
for di in range(smokePM.shape[0]):
    print(day_strs[di])
    dayPM = smokePM[di,:,:].flatten()
    day_df = pd.DataFrame({day_strs[di]:dayPM})
    smoke_df = pd.concat([smoke_df,day_df],axis=1)

# write to file
smoke_df.to_csv(proj_folder+'smokePM2018_wFIPs.csv',index=False)

#%% plot to check csv version looks ok
plt_map(glon,glat,area_mask_full,'jet','FIPS','fips gridded',clim=[0,50000])
plt_map(lon_csv,lat_csv,fips_csv,'jet','FIPS','fips csv',clim=[0,50000])

# california zoom on fips
cal_inds = np.where(np.logical_and(fips_csv>6000,fips_csv<6200))
plt_map(lon_csv[cal_inds],lat_csv[cal_inds],fips_csv[cal_inds],'jet','FIPS','fips csv',clim = [6000,6115])

# plot smoke PM in both versions
plt_map(glon,glat,smokePM[241,:,:],'magma','smokePM','smoke PM gridded',clim=[-1,20])
plt_map(lon_csv,lat_csv,smoke_df[day_strs[241]].values,'magma','smokePM','smoke PM csv',clim=[-1,20])
