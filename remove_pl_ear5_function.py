#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 16:59:23 2024

@author: x_tilin
"""
import pandas as pd
import netCDF4 as nc
import numpy as np
import os
from save_nc import save_to_netcdf


from datetime import datetime, timedelta
import time

def parse_args(data):
    import argparse

    parser = argparse.ArgumentParser(description="data root")

    parser.add_argument("--list-path", default="/home/x_tilin/snic2021-23-400/users/x_tilin/run/project4/Polar-low-list_Stoll_2020.csv", help="PL list root") 
    parser.add_argument("--data-path", default=data, help="era5 root")
    parser.add_argument("--newdata-path", default="/home/x_tilin/snic2021-23-400/users/x_tilin/input_data/Boundary/WRF/project4/removed_pl_era5/", help="removed PL era5 root")

    args = parser.parse_args()

    return args

def modify_matrix_x(matrix,lon):
    
    x=matrix.shape[1]
    mid=x//2
    shifted_matrix = np.hstack((matrix[:, mid:], matrix[:, :mid]))
    shifted_lon=np.hstack((lon[mid:], lon[:mid]+360))
    
    return shifted_matrix,shifted_lon,mid

def inputpath(inputdata):
    # 使用 os.path.join 来安全地拼接路径
    base_path = "/home/x_tilin/snic2021-23-400/users/x_tilin/input_data/Boundary/WRF/project4/"
    data = os.path.join(base_path, inputdata)
    # 调用 parse_args 函数，并传入拼接好的路径
    args = parse_args(data)
    return data, args


def removed_pl(variable,data_path):


    dataset = nc.Dataset(data_path, mode='r')
    era5_time=dataset.variables['time'][:]
    era5_lon=dataset.variables['longitude'][:]
    era5_lat=dataset.variables['latitude'][:]

    era5_para=np.full((era5_time.shape[0],era5_lat.shape[0], era5_lon.shape[0]), np.nan)
    # pl_center=np.full((10000,grid_length, grid_length), np.nan)
    # k=0
    
    for i in range(0,era5_time.shape[0],1):
    # for i in range(0,450,1):    
        print(i)
        
        mid_time=datetime.timestamp(datetime.strptime('1900-01-01 00:00:00', '%Y-%m-%d %H:%M:%S'))+era5_time[i]*3600
        time=datetime.fromtimestamp(mid_time)

        para = dataset.variables[variable][i] ###
        
        #locate the PL in ERA5 data
        for j in range(0,pl_time.shape[0],1):
        
            if pl_time.time[j] == str(time):
                print (j, time)
                
                pl_x=np.argmin(np.abs(era5_lon-pl_lon.Longitude[j]))
                pl_y=np.argmin(np.abs(era5_lat-pl_lat.Latitude[j]))
                
                grid_length=10
               
                
                if era5_lon[pl_x+grid_length] > max(era5_lon) or era5_lon[pl_x-grid_length] < min(era5_lon):
                    
                    shifted_x= np.argmin(np.abs(modify_matrix_x(para,era5_lon)[1]-pl_lon.Longitude[j]))    
    
                    # pl_center[k]=modify_matrix_x(para,era5_lon)[0][pl_y-grid_length:pl_y+grid_length,shifted_x-grid_length:shifted_x+grid_length]
                    modify_matrix_x(para,era5_lon)[0][pl_y-grid_length:pl_y+grid_length,shifted_x-grid_length:shifted_x+grid_length] = np.nan
                    mid=modify_matrix_x(para,era5_lon)[2]
                    para = np.hstack((modify_matrix_x(para,era5_lon)[0][:, mid:], modify_matrix_x(para,era5_lon)[0][:, :mid]))
                    
                    
                else:
          
                    # pl_center[k]=para[pl_y-grid_length:pl_y+grid_length,pl_x-grid_length:pl_x+grid_length]
                    para[pl_y-grid_length:pl_y+grid_length,pl_x-grid_length:pl_x+grid_length] = np.nan
                    
            # k+=1
            
        era5_para[i]=para
        
    return era5_para

data,args=inputpath("ERA5-2000-heat.nc")   ####input

pl_lon= pd.read_csv(args.list_path,usecols=['Longitude'])
pl_lat= pd.read_csv(args.list_path,usecols=['Latitude'])
pl_time= pd.read_csv(args.list_path,usecols=['time'])


data_path=args.data_path
variable = 'sshf'    ####input
era5_para=removed_pl(variable,data_path)

file_path = args.newdata_path + 'year2000_sshf.nc'
save_to_netcdf(data_path, file_path, era5_para, variable_name='sshf')


data_path=args.data_path
variable = 'slhf'    ####input
era5_para=removed_pl(variable,data_path)

file_path = args.newdata_path + 'year2000_slhf.nc'
save_to_netcdf(data_path, file_path, era5_para, variable_name='slhf')
            
   
            

    
           
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            

    
    























