#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import glob, os, sys, time
from smh import Session, utils
import astropy.table
from astropy.io import ascii
import numpy as np

'''
Written: 4-11-2025
Last Updated: 4-14-25

This file is used to correct the errors in SMHR and output the line data and abundance data into a usable format! 
Note that we must run run_errors.py prior to running this script!

run_errors.py will update the provided .smh file for a given star and ensure proper formatting and corrected systematic and statistical errors are accounted for.
make_abundance_tables.py (this file) adds in last second correctons to the line and abundance data before formatting the data into a usable ascii table for reporting purposes! 

Note that compared to just the output from SMHR, this new output will actually change the values for [X/H] and [X/Fe] for the lines and abundances. This is because for each line (and thus the final abundance 
of each element, the uncertainties are propagated via weights and a new weighted log(eps) is calculated for each line! The base log(eps) reported in the abundance pattern will match the original SMHR output, though!)
'''

hfs_species = [21.1, 25.0, 23.0, 23.1, 27.0, 56.1, 57.1, 60.1, 63.1]
#hfs_species = [56.1, 60.1, 63.1] 

def fix_lines(linetab,star):
    """
    Fix metadata of syntheses that were not imported correctly
    """
    t1 = ascii.read("/home/pnalamwa/SMHR_py38-mpl313/smhr/smh/data/linelists/Ji20_linelist_with_header.txt")
    t2 = t1 # = ascii.read("/home/kaitlin/streamstars/synthesis_list/s5_sorted_master_list_extra.txt")
    for i,line in enumerate(linetab):
        if np.isnan(line["expot"]) or np.isnan(line["loggf"]):
            species = line["species"]
            # Ignore molecules
            if species > 100: continue
            wave = line["wavelength"]
            # I don't know why this breaks but let's just hardcode it
            #if species==25.0 and wave==4034.0:
            #    linetab["wavelength"][i] = 4034.5
            #    linetab["expot"][i] = 0.00
            #    linetab["loggf"][i] = -0.81
            #    continue
            # try to find in t1
            dx = 10
            while True: 
                matches = (np.abs(t1["wavelength"] - wave) < dx) & (t1["species"]==species)
                if matches.sum() > 1:
                    dx -= .5
                elif matches.sum()<=1:
                    break
            if matches.sum()==1:
                linetab["wavelength"][i] = t1[matches]["wavelength"][0]
                linetab["expot"][i] = t1[matches]["expot"][0]
                linetab["loggf"][i] = t1[matches]["loggf"][0]
            elif matches.sum()>1:
                print("TOO MANY {} t1 on star {} index {} species {:.1f} wave {:.1f}".format(matches.sum(), star,line["index"],species,wave))
            else:
                # otherwise look in t2
                matches = (np.abs(t2["wavelength"] - wave) < 5) & (t2["species"]==species)
                if matches.sum()==1:
                    linetab["wavelength"][i] = t2[matches]["wavelength"][0]
                    linetab["expot"][i] = t2[matches]["expot"][0]
                    linetab["loggf"][i] = t2[matches]["loggf"][0]
                elif matches.sum()>1:
                    print("TOO MANY {} t2 on star {} index {} species {:.1f} wave {:.1f}".format(matches.sum(), star,line["index"],species,wave))
                else:
                    print("FAILED on star {} index {} species {:.1f} wave {:.1f}".format(star,line["index"],species,wave))
    return linetab

def fix_ch_lines(name, linetab):
    """
    Use the CH lines with [O/Fe] = +0.4 everywhere
    Add 0.1 dex to the per-line uncertainty to account for [O/Fe]
    """
    ii = linetab["species"]==106.0
    ix = np.where(ii)[0]
    assert len(ix) == 2, ix
    
    """
    if ii.sum() == 0:
        print("{} does not have CH measurements".format(name))
        return linetab
    print("Updating {} CH measurements".format(name))
    corrtab = ascii.read("ch_ofe_differences.org", format="fixed_width")
    ix = np.where(corrtab["star"]==name)[0]
    assert len(ix)==1, ix
    ix = ix[0]
    ab13 = corrtab[ix]["OFe04_4313"]
    ab23 = corrtab[ix]["OFe04_4323"]

    assert int(linetab[ix[0]]["wavelength"]) == 4310
    assert int(linetab[ix[1]]["wavelength"]) == 4323
    linetab["logeps"][ix[0]] = ab13
    linetab["logeps"][ix[1]] = ab23
    """
    #linetab["e_tot"][ix[0]] = np.sqrt(linetab["e_tot"][ix[0]]**2 + 0.10**2)
    #linetab["e_tot"][ix[1]] = np.sqrt(linetab["e_tot"][ix[1]]**2 + 0.10**2)
    #linetab["weight"][ix[0]] = 1/linetab["e_tot"][ix[0]]**2
    #linetab["weight"][ix[1]] = 1/linetab["e_tot"][ix[1]]**2
    linetab["e_stat"][ix[0]] = np.sqrt(linetab["e_stat"][ix[0]]**2 + 0.10**2)
    linetab["e_stat"][ix[1]] = np.sqrt(linetab["e_stat"][ix[1]]**2 + 0.10**2)
    linetab["e_tot"][ix[0]] = np.sqrt(linetab["e_stat"][ix[0]]**2 + linetab["e_sys"][ix[0]]**2)
    linetab["e_tot"][ix[1]] = np.sqrt(linetab["e_stat"][ix[1]]**2 + linetab["e_sys"][ix[1]]**2)
    ## Do the weights and fix e_tot later
    return linetab

def add_ba_errors(name, linetab):
    """
    Add 0.20 dex uncertainty to the 4554 and 4934 lines
    """
    ii = linetab["species"]==56.1
    ix = np.where(ii)[0]
    waves = linetab[ii]["wavelength"]
    for i, wave in zip(ix,waves):
        if int(wave)==4554 or int(wave)==4933:
            #linetab["e_tot"][i] = np.sqrt(linetab["e_tot"][i]**2 + 0.20**2)
            #linetab["weight"][i] = 1/linetab["e_tot"][i]**2
            linetab["e_stat"][i] = np.sqrt(linetab["e_stat"][i]**2 + 0.20**2)
            linetab["e_tot"][i] = np.sqrt(linetab["e_stat"][i]**2 + linetab["e_sys"][i]**2)
    return linetab

def add_cn_errors(name, linetab):
    """
    Add 0.30 dex uncertainty to the CN lines if detected
    """
    ii = linetab["species"]==607.0
    ix = np.where(ii)[0]
    for i in ix:
        #linetab["e_tot"][i] = np.sqrt(linetab["e_tot"][i]**2 + 0.30**2)
        #linetab["weight"][i] = 1/linetab["e_tot"][i]**2
        linetab["e_stat"][i] = np.sqrt(linetab["e_stat"][i]**2 + 0.30**2)
        linetab["e_tot"][i] = np.sqrt(linetab["e_stat"][i]**2 + linetab["e_sys"][i]**2)
    return linetab

def fix_o_loggf(name, linetab):
    """
    Update the loggfs and abundances for the forbidden lines.
    """
    # Storey + Zeippen 2000 for atomic data, Caffau+08 loggfs
    loggf6300 = -9.717
    loggf6363 = -10.185
    
    ii = linetab["species"]==8.0
    ix = np.where(ii)[0]
    for i in ix:
        if int(linetab["wavelength"][i]) == 6300:
            old_loggf = linetab["loggf"][i]
            new_loggf = loggf6300
        elif int(linetab["wavelength"][i]) == 6363:
            old_loggf = linetab["loggf"][i]
            new_loggf = loggf6363
        else:
            raise
        dabund = old_loggf - new_loggf
        linetab["loggf"][i] = new_loggf
        linetab["logeps"][i] = linetab["logeps"][i] + dabund
    return linetab

def add_al_errors(name, linetab):
    """
    Add 0.30 dex systematic uncertainty to the Al
    """
    ii = linetab["species"]==13.0
    ix = np.where(ii)[0]
    for i in ix:
        linetab["e_sys"][i] = np.sqrt(linetab["e_sys"][i]**2 + 0.30**2)
        linetab["e_tot"][i] = np.sqrt(linetab["e_stat"][i]**2 + linetab["e_sys"][i]**2)
    return linetab

def add_hfs_errors(name, linetab):
    for species in hfs_species:
        ii = linetab["species"]==species
        ix = np.where(ii)[0]
        for i in ix:
            linetab["e_sys"][i] = np.sqrt(linetab["e_sys"][i]**2 + 0.10**2)
            linetab["e_tot"][i] = np.sqrt(linetab["e_stat"][i]**2 + linetab["e_sys"][i]**2)
    return linetab
        

if __name__=="__main__":
    rho_Tg = 0.96
    rho_Tv = -0.82
    rho_gv = -0.87
    rho_TM = -0.37
    rho_gM = -0.21
    rho_vM = 0.01
    rhomat = utils._make_rhomat(rho_Tg=rho_Tg, rho_Tv=rho_Tv, rho_gv=rho_gv,
                                rho_TM=rho_TM, rho_gM=rho_gM, rho_vM=rho_vM)

    version_out = "new_errors" # name to add to new files to distinguish any old ones
    fnames = glob.glob("/home/pnalamwa/SMHR_py38-mpl313/smhr/smh/Analysis/all_of_M5_KOA_with_propogated_errors_and_updated_R9R21_4-10-25/*/*.smh")

    for fname in fnames:

        #set up the arrays and sessions for each star
        all_lines = []
        all_abunds = []

        name = os.path.basename(fname)[:-4]
        print(name)
        session = Session.load(fname)
        
        ## Extract line info
        linetab = utils.process_session_uncertainties_lines(session, rhomat)

        ## Manual updates to line info
        try:
            linetab = fix_ch_lines(name, linetab)
            linetab = add_cn_errors(name, linetab)

        except Exception as e:
            print("CH",e)
            print("CN",e)

        ## Add in the errors and changes

        linetab = add_ba_errors(name, linetab)
        linetab = fix_o_loggf(name, linetab)
        linetab = add_al_errors(name, linetab)
        linetab = add_hfs_errors(name, linetab)
        
        for species in [8.0, 106.0, 56.1, 607.0,
                        13.0] + hfs_species:
            ix = np.where((linetab["species"]==species) & np.isfinite(linetab["e_stat"]))[0]
            if len(ix) == 0: continue
            t = linetab[ix]
            delta = utils.struct2array(t["e_Teff","e_logg","e_vt","e_MH"].as_array())
            sigma_tilde = np.diag(t["e_tot"]**2) + (delta.dot(rhomat.dot(delta.T)))
            sigma_tilde_inv = np.linalg.inv(sigma_tilde)
            w = np.sum(sigma_tilde_inv, axis=1)
            linetab["weight"][ix] = w 
        
        ## Abundance summaries
        summarytab = utils.process_session_uncertainties_abundancesummary(linetab, rhomat)

        ## Fix up Fe I and II errors
        summarytab["e_XFe"][summarytab["species"]==26.0] = 0.
        summarytab["e_XFe"][summarytab["species"]==26.1] = 0.

        ## The covariances are nice and should be saved so they can be calculated as needed
        var_X, cov_XY = utils.process_session_uncertainties_covariance(summarytab, rhomat)
        
        ## Finish by adding limits
        linetab, summarytab = utils.process_session_uncertainties_limits(session, linetab, summarytab, rhomat)
        linetab = fix_lines(linetab,name)
        linetab = fix_o_loggf(name, linetab) # do it again for the upper limits; defined in a way that's fine
        
        #Get the data ready for saving purposes
        linetab["star"] = name
        summarytab["star"] = name
        all_lines.append(linetab)
        all_abunds.append(summarytab)
        all_lines = astropy.table.vstack(all_lines)
        all_abunds = astropy.table.vstack(all_abunds)

        #Now we save the data files

        file_path ='/'.join(fname.split("/")[:-1])
        print(file_path)

        all_lines.write(os.path.join(file_path, 'line_data_{}.txt'.format(version_out)), format="ascii.fixed_width_two_line", overwrite=True)
        all_abunds.write(os.path.join(file_path, 'abund_data_{}.txt'.format(version_out)), format="ascii.fixed_width_two_line", overwrite=True)
