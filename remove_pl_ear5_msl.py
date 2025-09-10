#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 13:37:21 2024

@author: x_tilin
"""

import pandas as pd
import netCDF4 as nc
import numpy as np
from scipy.interpolate import griddata
from datetime import datetime, timedelta
import time
def nearest_neighbor_interpolation(matrix):
    """
    对矩阵中的 NaN 值进行最近邻插值。
    
    参数:
    matrix : np.ndarray
        包含 NaN 值的输入矩阵。
    
    返回:
    np.ndarray
        插值后的矩阵，NaN 值被替换为最近邻的有效值。
    """
    # 获取矩阵的行列索引
    x, y = np.indices(matrix.shape)

    # 找到有效的（非 NaN）和无效的（NaN）位置
    valid_mask = ~np.isnan(matrix)  # 有效数据的位置
    invalid_mask = np.isnan(matrix)  # NaN 值的位置

    # 使用有效位置的数据进行插值，将无效位置的 NaN 值替换
    matrix[invalid_mask] = griddata(
        (x[valid_mask], y[valid_mask]),  # 有效数据的坐标
        matrix[valid_mask],              # 有效数据的值
        (x[invalid_mask], y[invalid_mask]),  # 要插值的目标位置
        method='linear'                 # 最近邻插值
    )

    return matrix
def parse_args(data):
    import argparse

    parser = argparse.ArgumentParser(description="data root")

    parser.add_argument("--listnorth-path", default="/home/x_tilin/snic2021-23-400/users/x_tilin/run/project4/Climatology_Northern_Hemisphere.csv", help="PL list root") 
    parser.add_argument("--listsouth-path", default="/home/x_tilin/snic2021-23-400/users/x_tilin/run/project4/Climatology_Southern_Hemisphere.csv", help="PL list root") 
    parser.add_argument("--data-path", default=data, help="era5 root")
    parser.add_argument("--newdata-path", default="/home/x_tilin/snic2021-23-400/users/x_tilin/input_data/Boundary/NEMO/nopl3h/", help="removed PL era5 root")

    args = parser.parse_args()

    return args
def modify_matrix_x(matrix,lon):
    
    x=matrix.shape[1]
    mid=x//2
    shifted_matrix = np.hstack((matrix[:, mid:], matrix[:, :mid]))
    shifted_lon=np.hstack((lon[mid:]-360, lon[:mid]))
    
    return shifted_matrix,shifted_lon,mid


VARS = [['mean_sea_level_pressure','msl'], ['10m_u_component_of_wind','u10'], ['10m_v_component_of_wind','v10'], ['2m_dewpoint_temperature','d2m'],['2m_temperature','t2m'],
        ['mean_surface_downward_long_wave_radiation_flux','msdwlwrf'],['mean_surface_downward_short_wave_radiation_flux','msdwswrf'], 
        ['mean_total_precipitation_rate','mtpr'],
        ['mean_snowfall_rate','msr']]

VARS = [['mean_sea_level_pressure','msl']]



YEARS = [x for x in map(str, range(2018, 2021))]





for V in VARS:
    for Y in YEARS:
        print(Y)
        target_filename = 'era5-hourly-'+V[0]+'_y'+Y+'.nc'

        data="/home/x_tilin/snic2021-23-400/users/x_tilin/input_data/Boundary/NEMO/ERA5/2006-2020/out/"+target_filename
        args = parse_args(data)


        dataset = nc.Dataset(args.data_path, mode='r')
        # era5_time=dataset.variables['time'][:]
        era5_time=dataset.variables['valid_time'][:]
        era5_lon=dataset.variables['longitude'][:]
        era5_lat=dataset.variables['latitude'][:]
        
        
        npl_lon= pd.read_csv(args.listnorth_path,usecols=['lon'])
        npl_lon.loc[npl_lon['lon']<0,'lon']+=360
        npl_lat= pd.read_csv(args.listnorth_path,usecols=['lat'])
        npl_time= pd.read_csv(args.listnorth_path,usecols=['time'])
        npl_time_step= pd.read_csv(args.listnorth_path,usecols=['PL_time_step'])



        spl_lon= pd.read_csv(args.listsouth_path,usecols=['lon'])
        spl_lon.loc[spl_lon['lon']<0,'lon']+=360
        
        spl_lat= pd.read_csv(args.listsouth_path,usecols=['lat'])
        spl_time= pd.read_csv(args.listsouth_path,usecols=['time'])
        spl_time_step= pd.read_csv(args.listsouth_path,usecols=['PL_time_step'])


# =============================================================================
# creat a nc file
# =============================================================================
        era5_var=np.full((era5_time.shape[0],era5_lat.shape[0], era5_lon.shape[0]), np.nan, dtype=np.float32)
        # npl_center=np.full((10000,grid_length, grid_length), np.nan)
        # k=0
        
        
        data_NC = nc.Dataset(args.newdata_path+target_filename, 'w', format='NETCDF4')
        
        data_NC.createDimension('longitude',len(era5_lon))
        data_NC.createDimension('latitude', len(era5_lat))
        data_NC.createDimension('time', None)
        
        
        xlon=data_NC.createVariable("longitude", 'f', ('longitude'))
        xlon.long_name="longitude"
        xlon.units  = "degrees_east"
        
        ylat=data_NC.createVariable("latitude", 'f', ('latitude'))
        ylat.long_name="latitude"
        ylat.units  = "degrees_north"
        
        t=data_NC.createVariable("time", 'i4', ('time'))
        t.units = 'hour'
        t.long_name="time"
        t.calendar="gregorian"
        
        value=data_NC.createVariable(V[1], 'f', ("time", "latitude", "longitude"))
        # value.units="m/s"
        # value.long_name="10m wind speed"
        value.standard_name=V[1]
        
        
        data_NC.variables['latitude'][:] = era5_lat
        data_NC.variables['longitude'][:] = era5_lon
        data_NC.variables['time'][:] = era5_time-min(era5_time)
        data_NC.variables[V[1]][:] = era5_var
        
        data_NC.close() 
# =============================================================================
# track pmc in north hemisphere
# =============================================================================
        for i in range(0,era5_time.shape[0],1):
            print(i)
            
            # mid_time=datetime.timestamp(datetime.strptime('1900-01-01 00:00:00', '%Y-%m-%d %H:%M:%S'))+era5_time[i]*3600
            mid_time=datetime.timestamp(datetime.strptime('1970-01-01 00:00:00', '%Y-%m-%d %H:%M:%S'))+era5_time[i]
            time=datetime.fromtimestamp(mid_time)
            # u10=dataset.variables['u10'][i]
            # v10=dataset.variables['v10'][i]
            para = dataset.variables[V[1]][i]  
            
            #locate the PL in ERA5 data
            for j in range(0,npl_time.shape[0],1):
                if npl_time.time[j] == str(time):
                  
                    # print (j, time)
                    
                    npl_x=np.argmin(np.abs(era5_lon-npl_lon.lon[j]))
                    npl_y=np.argmin(np.abs(era5_lat-npl_lat.lat[j]))
                    
                    grid_xlength= round(650/(111.32*np.cos(np.radians(abs(npl_lat.lat[j])))/4)) #650km is mean radius of PLs
                    grid_ylength=round(650/(111.32/4))
                    if npl_x+grid_xlength >= (era5_lon.shape) or npl_x-grid_xlength < 0:
                        
                        shifted_x= np.argmin(np.abs(modify_matrix_x(para,era5_lon)[1]-npl_lon.lon[j]))    
        
                        # npl_center[k]=modify_matrix_x(para,era5_lon)[0][npl_y-grid_length:npl_y+grid_length,shifted_x-grid_length:shifted_x+grid_length]
                        midpara=modify_matrix_x(para,era5_lon)[0]                     
                        median_value=np.percentile(midpara[npl_y-grid_ylength:npl_y+grid_ylength,shifted_x-grid_xlength:shifted_x+grid_xlength],20)             
                        midpara[npl_y-grid_ylength:npl_y+grid_ylength,shifted_x-grid_xlength:shifted_x+grid_xlength][midpara[npl_y-grid_ylength:npl_y+grid_ylength,shifted_x-grid_xlength:shifted_x+grid_xlength]<median_value]=np.nan
                        
                        midpara[~np.isnan(midpara)] = 0
                        midpara_var=modify_matrix_x(para,era5_lon)[0] 
                        midpara_var=midpara_var+midpara
                        
                        
                        midpara_var=nearest_neighbor_interpolation(midpara_var)
                        mid=modify_matrix_x(para,era5_lon)[2]
                        para = np.hstack((midpara_var[:, mid:], midpara_var[:, :mid]))
         
                        
                    else:
                        midpara=para.copy()
                        
                        # npl_center[k]=para[npl_y-grid_length:npl_y+grid_length,npl_x-grid_length:npl_x+grid_length]
                        
                        median_value=np.percentile(midpara[npl_y-grid_ylength:npl_y+grid_ylength,npl_x-grid_xlength:npl_x+grid_xlength],20)
                       
                        midpara[npl_y-grid_ylength:npl_y+grid_ylength,npl_x-grid_xlength:npl_x+grid_xlength][midpara[npl_y-grid_ylength:npl_y+grid_ylength,npl_x-grid_xlength:npl_x+grid_xlength]<median_value]=np.nan
                        midpara[~np.isnan(midpara)] = 0
                        midpara_var=para+midpara
                        
                        para=nearest_neighbor_interpolation(midpara_var)
                        # para[npl_y-grid_length:npl_y+grid_length,npl_x-grid_length:npl_x+grid_length] = np.nan


            
# =============================================================================
# track pmc in south hemisphere
# =============================================================================

            #locate the PL in ERA5 data
            for j in range(0,spl_time.shape[0],1):
                # print(j)
                if spl_time.time[j] == str(time):
                  
                    # print (j, time)
                    
                    spl_x=np.argmin(np.abs(era5_lon-spl_lon.lon[j]))
                    spl_y=np.argmin(np.abs(era5_lat-spl_lat.lat[j]))
                    
                    
                    

                    grid_xlength= round(650/(111.32*np.cos(np.radians(abs(spl_lat.lat[j])))/4)) #650km is mean radius of PLs
                    grid_ylength=round(650/(111.32/4))
                    if spl_x+grid_xlength >= (era5_lon.shape) or spl_x-grid_xlength < 0:
                        
                        shifted_x= np.argmin(np.abs(modify_matrix_x(para,era5_lon)[1]-spl_lon.lon[j]))    
        
                        # spl_center[k]=modify_matrix_x(para,era5_lon)[0][spl_y-grid_length:spl_y+grid_length,shifted_x-grid_length:shifted_x+grid_length]
                        midpara=modify_matrix_x(para,era5_lon)[0]

                        
                        median_value=np.percentile(midpara[spl_y-grid_ylength:spl_y+grid_ylength,shifted_x-grid_xlength:shifted_x+grid_xlength],20)

                        
                        midpara[spl_y-grid_ylength:spl_y+grid_ylength,shifted_x-grid_xlength:shifted_x+grid_xlength][midpara[spl_y-grid_ylength:spl_y+grid_ylength,shifted_x-grid_xlength:shifted_x+grid_xlength]<median_value]=np.nan
                        
                       # midpara[spl_y-grid_length:spl_y+grid_length,shifted_x-grid_length:shifted_x+grid_length] = np.nan
                        
                        midpara[~np.isnan(midpara)] = 0
                        midpara_var=modify_matrix_x(para,era5_lon)[0] 
                        midpara_var=midpara_var+midpara
                        
                        
                        midpara_var=nearest_neighbor_interpolation(midpara_var)
                        mid=modify_matrix_x(para,era5_lon)[2]
                        para = np.hstack((midpara_var[:, mid:], midpara_var[:, :mid]))
         
                        
                    else:
                        midpara=para.copy()
              
                        # spl_center[k]=para[spl_y-grid_length:spl_y+grid_length,spl_x-grid_length:spl_x+grid_length]
                        
                        median_value=np.percentile(midpara[spl_y-grid_ylength:spl_y+grid_ylength,spl_x-grid_xlength:spl_x+grid_xlength],20)
                       
                        midpara[spl_y-grid_ylength:spl_y+grid_ylength,spl_x-grid_xlength:spl_x+grid_xlength][midpara[spl_y-grid_ylength:spl_y+grid_ylength,spl_x-grid_xlength:spl_x+grid_xlength]<median_value]=np.nan
                        midpara[~np.isnan(midpara)] = 0
                        midpara_var=para+midpara
                        
                        para=nearest_neighbor_interpolation(midpara_var)
                    
                    

            
            data_NC = nc.Dataset(args.newdata_path+target_filename, mode='r+', format='NETCDF4')
       
            data_NC.variables[V[1]][i] = para
            data_NC.close()
            
        



            

            
            

    
    


s





















