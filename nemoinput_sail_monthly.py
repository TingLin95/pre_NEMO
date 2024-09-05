#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 12:20:19 2024

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


data="/home/x_tilin/snic2021-23-400/users/x_tilin/run/nemo_4.2.2/cfgs/my_test/EXP00/domain_cfg.nc"
args = parse_args(data)

dataset = nc.Dataset(args.data_path, mode='r')
dataB_depth=dataset.variables['gdept_0'][:]
dataB_lon=dataset.variables['nav_lon'][:]
dataB_lat=dataset.variables['nav_lat'][:]



new_matrix=np.ones((12,dataB_depth.shape[1],dataB_depth.shape[2],dataB_depth.shape[3]))



save_path='/home/x_tilin/snic2021-23-400/users/x_tilin/run/nemo_4.2.2/cfgs/my_test/EXP00/'

data_NC = nc.Dataset(save_path+'sali_ref_clim_monthly.nc', 'w', format='NETCDF4')

## define dimesions
data_NC.createDimension('x',new_matrix.shape[3])
data_NC.createDimension('y',new_matrix.shape[2])
data_NC.createDimension('z',new_matrix.shape[1])
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

t=data_NC.createVariable("time_counter", 'f4', ('time_counter',))
t.units = 'month'
t.long_name="time"
t.calendar="gregorian"

value=data_NC.createVariable("vosaline", 'i2', ("time_counter", "z", "y", "x"))
value.units="PSU"
value.long_name="salinity"
#value.standard_name="eastward_wind"


data_NC.variables['nav_lat'][:] = dataB_lat
data_NC.variables['nav_lon'][:] = dataB_lon
data_NC.variables['time_counter'][:] = np.linspace(0, 11, 12)
data_NC.variables['depth'][:] = dataB_depth[0,:,0,0]
data_NC.variables['vosaline'][:] = new_matrix


data_NC.close() 
