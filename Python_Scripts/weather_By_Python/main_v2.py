# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 16:06:59 2023

@author: win10admin
"""


from datetime import datetime, timedelta,date

from Update_WTH import *
from CHIRPS import *
from GetNassaP import *
from Dssat_WTH import dssat_wth
#s1 = datetime.now()

def main():
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(dest='command')
    getwth = subparser.add_parser('get')
    updatewth = subparser.add_parser('update')

    getwth.add_argument('in_file', type=str, help='CSV file with the points required. It must contain ID(CellID), Latitude, Longitude')
    getwth.add_argument('startDate', type=int, help='Start date with format YYYYMMDD (e.g. 19841224)')
    getwth.add_argument('endDate', type=int, help='End date with format YYYYMMDD (e.g. 19841231)')
    getwth.add_argument('out_dir', type=str, help='Path of output directory for the new WTH files.')

    updatewth.add_argument('in_file', type=str, help='CSV file with the points required. It must contain ID(CellID), Latitude, Longitude.')
    updatewth.add_argument('in_dir', type=str, help='Path directory of current WTH files to update.')
    updatewth.add_argument('out_dir', type=str, help='Path of output directory for the new WTH files.')

    args = parser.parse_args()

    if args.command == 'get':
        dssat_wth(args.in_file, args.startDate, args.endDate, args.out_dir)
    elif args.command == 'update':
        update_wth(args.in_file, args.in_dir, args.out_dir)

if __name__ == "__main__":
    sys.exit(main())
    
#To download the weather data run this in the terminal:
# >python main_v2.py get "C:/Old__CCAFSToolkit/CCAFSToolkit/weather_By_Python/Input.csv" 20220101 20220331  "C:/Old__CCAFSToolkit/CCAFSToolkit/weather_By_Python/Data/"
#to update, run the folow line:   
#python main_v2.py update "C:/Old__CCAFSToolkit/CCAFSToolkit/weather_By_Python/Input.csv"  "C:/Old__CCAFSToolkit/CCAFSToolkit/weather_By_Python/Data/"  "C:/Old__CCAFSToolkit/CCAFSToolkit/weather_By_Python/Data/"