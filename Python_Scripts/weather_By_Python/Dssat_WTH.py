#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 18:36:07 2023

@author: yacoub
"""
from GetNassaP import *
get_correc_nc
def dssat_wth(in_file, startDate, endDate, out_dir):
    s1 = datetime.now()

    os.chdir(os.path.dirname(in_file))
    tempdir = os.path.dirname(in_file) + '/temp'
    nasa_outdir= tempdir+ '/nasap'
    os.makedirs(nasa_outdir, exist_ok=True)
    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)
        os.mkdir(tempdir)
    else:
        os.mkdir(tempdir)

    print('Getting NASA POWER data...')
    #nasa_outdir = os.path.join(tempdir, 'nasap')
    os.makedirs(nasa_outdir, exist_ok=True)
    
    #Getting NASA POWER data for the update period
    nasa(in_file, str(startDate), str(endDate), nasa_outdir)

    #Getting corrected data
    print('Getting corrected data from CHIRPS server...')
    out_cor_nc = tempdir + '/in_nc_cor'
    dt_s = datetime.strptime(str(startDate), '%Y%m%d')
    dt_e = datetime.strptime(str(endDate), '%Y%m%d')
    get_correc_nc(dt_s, dt_e, out_cor_nc)

    #Run chirps for corrected data
    outdir_prec = tempdir + '/prec_pkl'
    print('Processing CHIRPS data...')
    chirps1(in_file, out_cor_nc, outdir_prec + '/prec_corr.pkl')

    #Getting the latest day available in prec corrected data.
    df = joblib.load(outdir_prec + '/prec_corr.pkl')
    lastday_corr = df.columns[-1]
    dt_s_p = datetime.strptime(lastday_corr, '%Y%j') + timedelta(days=1)

    if dt_s_p < dt_e:
        #Getting preliminary data
        print('Getting preliminary data from CHIRPS server...')
        out_pre_nc = tempdir + '/in_nc_pre'
        get_prelim_nc(dt_s_p, dt_e, out_pre_nc)
        print('CHIRPS netCDF files in disk.')

        #Run chirps for preliminary data
        chirps1(in_file, out_pre_nc, outdir_prec + '/prec_prelim.pkl')

        #Merging corrected and preliminary precipitation data.
        precpkl(outdir_prec)
        print('CHIRPS processing data are complete.')

    else:
        os.rename(outdir_prec + '/prec_corr.pkl', outdir_prec + '/prec.pkl')

    #Fusing NASA POWER and CHIRPS with QC on SRAD.
    print('Building the WTH files...')
    nasachirps(in_file,startDate, endDate, nasa_outdir, outdir_prec + '/prec.pkl', out_dir)

    e1 = datetime.now()
    print("Time for execution is: ", str(e1-s1))
    
