#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 18:39:59 2023

@author: yacoub
"""
from datetime import datetime, timedelta,date

from Update_WTH import *
from CHIRPS import *
from GetNassaP import *
from Dssat_WTH import *
#s1 = datetime.now()


# Définissez la classe main
class main:
    def download(self, in_file, startDate, endDate, out_dir):
        dssat_wth(in_file, startDate, endDate, out_dir)
    
    def Update(self, in_file,in_dir, out_dir):
        update_wth(in_file, in_dir, out_dir)

# Create a class instance
m = main()

in_dir=out_dir = "C:/Old__CCAFSToolkit/CCAFSToolkit/weather_By_Python/Data/"


if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    print(f"Folder created: {out_dir}")
else:
    print(f"Folder already exists: {out_dir}")

#out_dir2 = "/Users/yacoub/Desktop/UF/spyder/WTH_Files/temp/NasaChirps_DataUpdate/"
#os.makedirs(in_dir, exist_ok=True)

in_file="C:/Old__CCAFSToolkit/CCAFSToolkit/weather_By_Python/Input.csv"

startDate = 20220601
endDate = 20220731


# Call the download() method with the specified values
m.download(in_file, startDate, endDate, out_dir)
#m.Update(in_file,out_dir, out_dir)


