"""
X-ray Fitting Tool (XFT)
Author: Nis C. Gellert
Technichal University of Denmark, DTU Space
Last updated: 19/09/2023
"""

import numpy as np
from datetime import date

def save(fname,angleScan,sample_structure, theta_m, reflectance_bf, data,rho_bf, sigma_bf, z_bf, minFOM, a_res, energy_i, logFitting, weightsFitting, ang_off,eng_off):

    if angleScan:
        file1 = open(fname,"w+") 
        today = date.today()
        file1.write("# Saved on %s\n\n" % today.strftime("%d/%m/%Y"))
        file1.write("# Sample structure:\n")
    
        for i in range(len(sample_structure)-1):
            file1.write("# %s rho = %s (g/cm^3), sigma = %s (nm), z = %s (nm) \n" % (sample_structure[i][0:2],  str(round(rho_bf[i],2)),  str(round(sigma_bf[i],2)),  str(round(z_bf[i],1))))   
        file1.write("# %s rho = %s (g/cm^3), sigma = %s (nm) \n\n" % (sample_structure[-1][0:2], str(round(rho_bf[-1],2)), str(round(sigma_bf[-1],2))))
        file1.write("# FOM = %.2E \n" % minFOM)
        file1.write("# Angular Offset on model (degree): %s \n" % ang_off)
        file1.write("# Incident engery (ev):  %s \n" % energy_i)
        file1.write("# Angular resolution (degree): %s \n" % a_res )
        file1.write("# Logaritmic fitting: %d \n" % (logFitting))
        file1.write("# Weighting: %s \n\n"% weightsFitting)    
        file1.write("# Theta (deg), Reflectance (model), Residual, Reflectance (data)\n" )
    
        for i in range(len(theta_m)):
             file1.write("%.4E %.9E %.9E %.9E\n" % (theta_m[i], reflectance_bf[i], data[i]-reflectance_bf[i], data[i]))
        file1.close() 
        
        
    else:
        file1 = open(fname,"w+") 
        today = date.today()
        file1.write("# Saved on %s\n\n" % today.strftime("%d/%m/%Y"))
        file1.write("# Sample structure:\n")
        
        for i in range(len(sample_structure)-1):
            file1.write("# %s rho = %s (g/cm^3), sigma = %s (nm), z = %s (nm) \n" % (sample_structure[i][0:2],  str(round(rho_bf[i],2)),  str(round(sigma_bf[i],2)),  str(round(z_bf[i],2))))   
        file1.write("# %s rho = %s (g/cm^3), sigma = %s (nm) \n\n" % (sample_structure[-1][0:2], str(round(rho_bf[-1],2)), str(round(sigma_bf[-1],2))))
        file1.write("# FOM = %.2E \n" % minFOM)
        file1.write("# Energy Offset on model (eV): %s \n" % eng_off)
        file1.write("# Incident angle (degree):  %s \n" % theta_m)
        file1.write("# Angular resolution (degree): %s \n" % a_res )
        file1.write("# Logaritmic fitting: %d \n" % (logFitting))
        file1.write("# Weighting: %s \n\n"% weightsFitting)    
        file1.write("# Energy (eV), Reflectance (model), Residual, Reflectance (data)\n" )
    
        for i in range(len(energy_i)):
             file1.write("%.4E %.9E %.9E %.9E\n" % (energy_i[i], reflectance_bf[i], data[i]-reflectance_bf[i], data[i]))
        file1.close() 
        
        
        
        
        
        
        
        
        