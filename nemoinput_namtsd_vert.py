#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 10:30:51 2024

@author: x_tilin
"""


import netCDF4 as nc
import numpy as np
from scipy.interpolate import griddata
from scipy import interpolate
from datetime import datetime, timedelta
import time

def parse_args(data):
    import argparse

    parser = argparse.ArgumentParser(description="data root")
    parser.add_argument("--data-path", default=data, help="era5 root")
    args = parser.parse_args()

    return args
dataname=['data_1m_salinity_nomask.nc','data_1m_potential_temperature_nomask.nc']
var=['vosaline','votemper']
outputname=['newtsd_sss.nc','newtsd_sst.nc']
lon=['nav_lon','lon']
lat=['nav_lat','lat']


idx=0

data="/home/x_tilin/snic2021-23-400/users/x_tilin/run/nemo_4.2.2/cfgs/my_test/EXP00/"+dataname[idx]
args = parse_args(data)

dataset = nc.Dataset(args.data_path, mode='r')
dataA=dataset.variables[var[idx]][:].filled(np.nan)
dataA_lon=dataset.variables[lon[idx]][:]
dataA_lat=dataset.variables[lat[idx]][:]
dataA_depth=dataset.variables['depth'][:]



data="/home/x_tilin/snic2021-23-400/users/x_tilin/run/nemo_4.2.2/cfgs/my_test/EXP00/domain_cfg.nc"
args = parse_args(data)
dataset = nc.Dataset(args.data_path, mode='r')
dataB_depth=dataset.variables['gdept_0'][:]

x=dataA_depth.filled(np.nan)
new_x=dataB_depth[0,:,0,0].filled(np.nan)




new_matrix_zeros = np.zeros((dataA.shape[0],dataB_depth.shape[1],dataA.shape[2],dataA.shape[3]))
for i in range(0,dataA.shape[0]):
    for j in range(0,dataA.shape[2]):
        for k in range(0,dataA.shape[3]):   
            
                y = dataA[i,:,j,k]   
                f = interpolate.interp1d(x, y, kind='slinear',fill_value='extrapolate', bounds_error=False)
                new_matrix_zeros[i,:,j,k]=f(new_x)
        

        
if idx==0:
    save_path='/home/x_tilin/snic2021-23-400/users/x_tilin/run/nemo_4.2.2/cfgs/my_test/EXP00/'
    
    data_NC = nc.Dataset(save_path+outputname[idx], 'w', format='NETCDF4')
    
    ## define dimesions
    data_NC.createDimension('x',new_matrix_zeros.shape[3])
    data_NC.createDimension('y',new_matrix_zeros.shape[2])
    data_NC.createDimension('z',new_matrix_zeros.shape[1])
    data_NC.createDimension('time_counter',None)
    
    
    xlon=data_NC.createVariable("nav_lon", 'f', ('y','x'))
    xlon.long_name="longitude"
    xlon.units  = "degrees_east"
    
    
    ylat=data_NC.createVariable("nav_lat", 'f', ('y','x'))
    ylat.long_name="latitude"
    ylat.units  = "degrees_north"
    
    
    p_lvl=data_NC.createVariable("depth", 'i4', ('z'))
    p_lvl.units = "m"
    p_lvl.long_name="depth"
    
    t=data_NC.createVariable("time_counter", None)
    t.units = 'month'
    t.long_name="time"
    t.calendar="gregorian"
    
    value=data_NC.createVariable(var[idx], 'f', ("time_counter", "z", "y", "x"))
    value.units="PSU"
    value.long_name="salinity"
    #value.standard_name="eastward_wind"
    
    
    data_NC.variables['nav_lat'][:] = dataA_lat
    data_NC.variables['nav_lon'][:] = dataA_lon
    data_NC.variables['time_counter'][:] = 11
    data_NC.variables['depth'][:] = new_x
    data_NC.variables[var[idx]][:] = new_matrix_zeros
    
    
    data_NC.close() 
    
else:
    
    save_path='/home/x_tilin/snic2021-23-400/users/x_tilin/run/nemo_4.2.2/cfgs/my_test/EXP00/'

    data_NC = nc.Dataset(save_path+outputname[idx], 'w', format='NETCDF4')

    ## define dimesions
    data_NC.createDimension('x',new_matrix_zeros.shape[3])
    data_NC.createDimension('y',new_matrix_zeros.shape[2])
    data_NC.createDimension('z',new_matrix_zeros.shape[1])
    data_NC.createDimension('time_counter',None)


    xlon=data_NC.createVariable("lon", 'f', ('x'))
    xlon.long_name="longitude"
    xlon.units  = "degrees_east"


    ylat=data_NC.createVariable("lat", 'f', ('y'))
    ylat.long_name="latitude"
    ylat.units  = "degrees_north"


    p_lvl=data_NC.createVariable("depth", 'i4', ('z'))
    p_lvl.units = "m"
    p_lvl.long_name="depth"

    t=data_NC.createVariable("time_counter", None)
    t.units = 'month'
    t.long_name="time"
    t.calendar="gregorian"

    value=data_NC.createVariable(var[idx], 'f', ("time_counter", "z", "y", "x"))
    value.units="Celsius"
    value.long_name="Potential Temperature : WOA09 + PHC3 (up to 65N)"
    #value.standard_name="eastward_wind"


    data_NC.variables['lat'][:] = dataA_lat
    data_NC.variables['lon'][:] = dataA_lon
    data_NC.variables['time_counter'][:] = 11
    data_NC.variables['depth'][:] = new_x
    data_NC.variables[var[idx]][:] = new_matrix_zeros


    data_NC.close() 