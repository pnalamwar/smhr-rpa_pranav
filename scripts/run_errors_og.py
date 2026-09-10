#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import glob, os, sys, time
import numpy as np
import pandas as pd
from smh import Session

'''
Written: 4-11-2025
Last Updated: 4-14-2025

This code will go througha list of .smh files and set the systematic and statistical uncertainties in the stellar parameters. Then it propagates the errors to the abundances of every spectral model available
This code is designed to be run before make_abun_tables.py that will actually do some last minute corrections and produce an ascii table of abundances with errors!
'''

if __name__=="__main__":
    smh_fnames = glob.glob("/home/pnalamwa/SMHR_py38-mpl313/smhr/smh/Analysis/all_of_M5_KOA_with_propogated_errors_and_updated_R9R21_4-10-25/*/*.smh")
    error_file = pd.read_csv('/home/pnalamwa/CRC_codes_data/M5_Star_Data/final_stars_with_errors_with_R9R21.csv')

    
    #fname = sys.argv[1]
   # newfname = sys.argv[2]
   # assert not os.path.exists(newfname), f"{newfname} already exists! Stopping..."
    #eTeff = float(sys.argv[2])
    #elogg = float(sys.argv[3])
    #evt = float(sys.argv[4])
    #eMH = float(sys.argv[5])
    #assert os.path.exists(fname)
    #smh_fnames = [fname]

    for fname in smh_fnames:
        print("Processing {}".format(fname))
        start = time.time()
        #newfname = fname.replace("syntheses","errors")
        #if os.path.exists(fname):
            #print("Already done: {}".format(fname))
        session = Session.load(fname)

        #session.stellar_parameter_uncertainty_analysis() #systematic_errors=[83,0.041,0.020,0.0])
        #session.stellar_parameter_uncertainty_analysis(systematic_errors=[50, 0.15, 0.10, 0.05])
        #eTeff, elogg, evt, eMH = session.stellar_parameters_staterr


        #find the stellar parameter errors by associating it with the name of the star!
        star_name = fname.split('/')[-2]
        print('Star name is: ', star_name)

        index_star = error_file[error_file['KOA Name'] == star_name].index[0]
        print('star name is:', star_name, 'and star index is : ', index_star)  

        dteff = error_file.iloc[index_star]['T_eff_error']
        dlogg = error_file.iloc[index_star]['Log(g)_error']
        dvt = error_file.iloc[index_star]['xi_error']
        dfeh = 0.0 

        print('\n Stellar parameter errors are: ', dteff, dlogg, dvt, dfeh, '\n') 

        #Finally, set the stellar parameter errors and run the uncertainties on the abundances again!

        session.set_stellar_parameters_errors("sys", dteff, dlogg, dvt, dfeh) #teff, logg, vt, mh 
        session.set_stellar_parameters_errors("stat",0,0,0,0) 
        print("We set the systematic stellar parameter errors to: {}".format(session.stellar_parameters_err))
        session.compute_all_abundance_uncertainties() 
        
        try:
            session.save(fname, overwrite=True)
            print("Finished {} in {:.1f}s".format(fname, time.time()-start)) 
        except OSError:
            print('Error in saving the file: ', fname)
            pass

        
