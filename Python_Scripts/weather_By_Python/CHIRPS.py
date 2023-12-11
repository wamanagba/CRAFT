#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 18:09:03 2023

@author: yacoub
"""

import os, sys
from osgeo import ogr, gdal
from osgeo.gdalconst import *
import numpy
import pandas as pd
from datetime import datetime, timedelta,date
import joblib
import requests
from dateutil.relativedelta import relativedelta
import argparse
import shutil
#from chirps import *
import joblib
#from getnasap import nasa, nasachirps

import queue
import threading
import time
import logging
# register all of the GDAL drivers
gdal.AllRegister()


def get_correc_nc(dt_s, dt_e, out_cor_nc):
    # Create a session for HTTP requests with custom parameters
    s = requests.Session()
    s.mount("https://data.chc.ucsb.edu", requests.adapters.HTTPAdapter(max_retries=10))

    # Calculate the difference in months between the start date and end date
    diff_month = (dt_e.year - dt_s.year) * 12 + (dt_e.month - dt_s.month)

    # Iterate over each month between the start date and end date, inclusive
    for n in range(diff_month+1):
        # Calculate the current month by adding n months to the start date
        yymm = dt_s + relativedelta(months=+n)
        yy = yymm.strftime("%Y")  # Extract the year in "YYYY" format
        mm = yymm.strftime("%m")  # Extract the month in "MM" format
        
        
      
        file_name = 'corr_chirps_' + yy + mm + '.nc'  # File name to download
        file_path = os.path.join(out_cor_nc, file_name)  # Full file path
        

        try:
            # Download the file from the specified URL
            response = s.get('https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf/p05/by_month/chirps-v2.0.'
                             + yy + '.' + mm + '.days_p05.nc', timeout=80)
            response.raise_for_status()  # Check if the request was successful (HTTP status code 200)

            if response.ok:
                # Create the destination folder if it doesn't already exist
                if not os.path.exists(out_cor_nc):
                    os.makedirs(out_cor_nc)

                # Save the downloaded file to the destination folder
                with open(file_path, 'wb') as file:
                    file.write(response.content)
                print(file_name + ' downloaded successfully.')

        except requests.exceptions.HTTPError as err:
            pass
            # raise SystemExit(err)

        response.close()  # Close the connection with the server




def get_prelim_nc(dt_s, dt_e, out_pre_nc):
    # Create a session for HTTP requests with custom parameters
    s = requests.Session()
    s.mount("https://data.chc.ucsb.edu", requests.adapters.HTTPAdapter(max_retries=10))

    # Iterate over each year between the start year and end year, inclusive
    for y in range(dt_s.year, dt_e.year + 1):
        single_y = str(y)  # Convert the current year to a string

        file_name = 'prelim_nc_' + single_y + '.nc'  # File name to download
        file_path = os.path.join(out_pre_nc, file_name)  # Full file path

        # Check if the file already exists in the destination folder
        #if os.path.exists(file_path):
        #    print(file_name + ' already exists. Skipping download.')
        #    continue  # Move to the next file without downloading

        try:
            # Download the file from the specified URL
            response = s.get('https://data.chc.ucsb.edu/products/CHIRPS-2.0/prelim/global_daily/fixed/netcdf/chirps-v2.0.'
                             + single_y + '.days_p05.nc', timeout=80)
            response.raise_for_status()  # Check if the request was successful (HTTP status code 200)

            if response.ok:
                # Create the destination folder if it doesn't already exist
                if not os.path.exists(out_pre_nc):
                    os.makedirs(out_pre_nc)

                # Save the downloaded file to the destination folder
                with open(file_path, 'wb') as file:
                    file.write(response.content)
                print(file_name + ' downloaded successfully.')

        except requests.exceptions.HTTPError as err:
            pass
            # raise SystemExit(err)

        response.close()  # Close the connection with the server


# Intended for long time series but few points (<10000)
def chirps1(in_file, in_nc_dir, outprec_file):
    nc_lst = os.listdir(in_nc_dir)  # Retrieve the list of all .sol files in the input folder.
    nc_lst.sort()  # Sort the files in chronological order.
    df_chirps = pd.DataFrame()  # Create an empty DataFrame to store precipitation data.

    with open(in_file, "r") as f1:
        
        in_pt = [line for line in f1.readlines() if line.strip()]  # Read lines from the input file, excluding empty lines.

        for row in in_pt[1:]:
            id = int(row.split(',')[0])  # Get the ID from the first column.
            lat = float(row.split(',')[1])  # Get the latitude from the second column.
            lon = float(row.split(',')[2])  # Get the longitude from the third column.
            time_lst = []  # List to store time values.
            precval = []  # List to store precipitation values.

            for nc_file in nc_lst:
                if nc_file.endswith(".nc"):  # Check if the file is in NetCDF format.
                    dsi = gdal.Open(in_nc_dir + "/" + nc_file, GA_ReadOnly)  # Open the NetCDF file in read-only mode.
                    if dsi is None:
                        print('Could not open NetCDF file')
                        sys.exit(1)  # Print an error message and exit if the file cannot be opened.
                    meta_nc = dsi.GetMetadata()  # Get the metadata of the NetCDF file.
                    date_start = meta_nc['time#units'][-14:]  # Get the start date from the metadata.
                    datetime_st = datetime.strptime(date_start, '%Y-%m-%d %H:%M:%S')  # Convert the start date to a datetime object.
                    bands_time = meta_nc['NETCDF_DIM_time_VALUES'][1:-1].split(',')  # Get the time dimension values.
                    bands_time = list(map(int, bands_time))  # Convert the values to integers.

                    gt = dsi.GetGeoTransform()  # Get the geotransformation of the NetCDF file.
                    px = int((lon - gt[0]) / gt[1])  # Calculate the X coordinate of the pixel corresponding to the longitude.
                    py = int((lat - gt[3]) / gt[5])  # Calculate the Y coordinate of the pixel corresponding to the latitude.
                    bands = dsi.RasterCount  # Get the number of bands in the NetCDF file.

                    for i in range(1, bands + 1):
                        d = dsi.GetRasterBand(i).ReadAsArray(px, py, 1, 1)  # Read the precipitation value for the given pixel.
                        dt_st = datetime_st + timedelta(days=bands_time[i - 1])  # Calculate the date corresponding to the band.
                        band_t = dt_st.strftime('%Y%j')  # Convert the date to YJJJ format (year and day of the year).
                        time_lst.append(band_t)  # Add the time value to the list.
                        if d is None:
                            d = numpy.float32([[-9999.0]])  # Replace missing values with a specific code.
                        precval.append(d[0][0])  # Add the precipitation value to the list.

            df_chirps[id] = precval  # Associate the list of precipitation values with the ID in the DataFrame.

    df_chirps = df_chirps.T  # Transpose the DataFrame to have dates as columns.
    df_chirps.index.name = 'ID'  # Assign a name to the DataFrame's index.
    df_chirps.columns = time_lst  # Associate the time values with column names.

    if not os.path.exists(os.path.dirname(outprec_file)):
        os.mkdir(os.path.dirname(outprec_file))  # Create the destination directory of the output file if it doesn't exist.
    joblib.dump(df_chirps, outprec_file)  # Save the DataFrame to a file in pickle format.



def precpkl(outdir_prec):
    # Load the "prec_corr.pkl" file into a dataframe df1
    df1 = joblib.load(outdir_prec + '/prec_corr.pkl')
    
    # Load the "prec_prelim.pkl" file into a dataframe df2
    df2 = joblib.load(outdir_prec + '/prec_prelim.pkl')
    
    # Extract the last date from the df1 dataframe (last column)
    last_dt_corr = df1.columns[-1]
    
    # Calculate the next date based on the last corrected date
    dt_s_pre = datetime.strptime(last_dt_corr, '%Y%j') + timedelta(days=1)
    
    # Convert the next date to a string format
    dt_st_p = dt_s_pre.strftime('%Y%j')
    
    # Select the columns from the df2 dataframe starting from the next date
    df3 = df2.loc[:, dt_st_p:]
    
    # Concatenate df1 and df3 dataframes along the column axis
    result = pd.concat([df1, df3], axis=1)
    
    # Save the result dataframe to a "prec.pkl" pkl file
    joblib.dump(result, outdir_prec + '/prec.pkl')


