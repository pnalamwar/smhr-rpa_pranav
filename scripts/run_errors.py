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
Last Updated: 7-27-2026 (updated for specifically actinide boost star)

This code will go througha list of .smh files and set the systematic and statistical uncertainties in the stellar parameters. Then it propagates the errors to the abundances of every spectral model available
This code is designed to be run before make_abun_tables.py that will actually do some last minute corrections and produce an ascii table of abundances with errors!
'''

if __name__=="__main__":
    smh_fnames = glob.glob("/home/pnalamwa/actinide_boost_star_studies/*actinide_boost_star_file_synth_with_errors.smh")
    
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

        dteff = 45
        dlogg = 0.145
        dvt = 0.052
        dfeh = 0.031

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

        
