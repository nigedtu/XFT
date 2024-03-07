"""
X-ray Fitting Tool (XFT)
Author: Nis C. Gellert
Technichal University of Denmark, DTU Space
Last updated: 19/09/2023


Input units: 
    lambda_i: incident wavelength, nm
    density: density of material, g/cm3
    number_of_atoms and elements define the composition of the compound 
    default source for calculating optical constants: CXRO

Output: 
    Returns the refractive indices of a given material at the specified wavelengths

"""

import numpy as np
import matplotlib.pyplot as plt
import re
from unit_conversion import wavelength2energy, energy2wavelength
import math

def oc(lambda_i, density, number_of_atoms, elements, oc_source = "oc_source//CXRO//"):

    
    f1_interp = np.zeros(len(elements))
    f2_interp = np.zeros(len(elements))
    AW = np.zeros(len(elements))
    
    
    if np.isscalar(lambda_i): 
        # angle scan, one wavelength
        lambda_i= np.array([lambda_i]) # array of length 1
    
    n = np.zeros(len(lambda_i), dtype=np.complex)
    density = density*1.e+6*np.sum(number_of_atoms) #  [g/m3]
    lambda_i = lambda_i* 1e-9 # converting from nm to m
    
    for j in range(len(lambda_i)): # looping over all wavelengths given in lambda_i
        lambda_temp = lambda_i[j]
        
        for a in range(len(elements)):
            fname = (oc_source + elements[a] + ".txt")
            data = np.loadtxt(fname)
            # Header data: Atomic weight, density, relativistic corrections
            with open(fname, 'r') as the_file:
                all_data = [line.strip() for line in the_file.readlines()]
                AW[a] = re.findall('\d*\.?\d+',all_data[2])[1]
                f_rel =  re.findall('\d*\.?\d+',all_data[4])[2]
                f_NT =  re.findall('\d*\.?\d+',all_data[5])[1]
           
            lambda_source = data[:,3]*1.e-9 # wavelength [m]
            f1 = data[:,1] # [Unitless, e atom-1] Atomic form factor
            f2 = data[:,2] # [Unitless, e atom-1] Atomic form factor
            
            
            # Get index for scalar lambda value (angle scan)
            f1_interp[a] = np.interp(lambda_temp, np.flip(lambda_source), np.flip(f1))
            f2_interp[a] = np.interp(lambda_temp, np.flip(lambda_source), np.flip(f2))
            
            r0 =  2.8179403227e-15 # [m] Classical electron radius
            N0 = 6.0221409e+23 # [mol^-1] Avogadros number
        
            f1_interp[a] = f1_interp[a] # + f_rel# + f_NT # [Unitless, e atom-1] Atomic form factor
            
            f1_interp[a] = (number_of_atoms[a]/np.sum(number_of_atoms))*f1_interp[a]
            f2_interp[a] = (number_of_atoms[a]/np.sum(number_of_atoms))*f2_interp[a]
            AW[a] = AW[a]*number_of_atoms[a] #  g/mol
            
        
        delta = (N0*density*r0)/(np.sum(AW)*2*np.pi)*lambda_temp**2.*np.sum(f1_interp) #  [unit less]
        beta = (N0*density*r0)/(np.sum(AW)*2*np.pi)*lambda_temp**2*np.sum(f2_interp) # [unit less]
        n[j] = (1 - delta) + complex(0,1)*beta # [unit less]
        
        
        
    return n 
    

def oc_MO(energy_i, elements_m, oc_source = "oc_source//LLNL_MO//"):
    # Energy eV, elements "Ir", OC_source='oc_source/LLNL_MO/'
    energy_m = energy_i/1000 # keV
    energy_m = np.around(energy_m,1) # 
    #n = np.zeros(len([energy_m]), dtype=np.complex)
    #idx = np.zeros(len([energy_m]))
    #fname = (oc_source + elements_m + ".txt")
    #data = np.loadtxt(fname)
    #energy_source = data[:,0] # [keV]
    #energy_s=energy_source.tolist()
    #n = data[:,1] # [Unitless, e atom-1] Atomic form factor
    #k = data[:,2] # [Unitless, e atom-1] Atomic form factor
    if np.isscalar(energy_m): # Angle scan
        n = np.zeros(len([energy_m]), dtype=np.complex)
        idx = np.zeros(len([energy_m]))
        fname = (oc_source + elements_m + ".txt")
        data = np.loadtxt(fname)
        energy_source = data[:,0] # [keV]
        energy_s=energy_source.tolist()
        idx = [i for i, e in enumerate(energy_s) if e == energy_m]
        n = data[idx,1]
        k = data[idx,2]
        n = n + complex(0,1)*k
    else: # Energy scan
        n = np.zeros(len(energy_m), dtype=np.complex)
        idx = np.zeros(len(energy_m))
        fname = (oc_source + elements_m + ".txt")
        data = np.loadtxt(fname)
        energy_source = data[:,0] # [keV]
        energy_s=energy_source.tolist()
        for j in range(len(n)):
            idx[j] = [i for i, e in enumerate(energy_s) if e == energy_m[j]][0]
        idx_list = [int(i) for i in idx.tolist()]
        n = data[idx_list,1]
        k = data[idx_list,2]
        n = n + complex(0,1)*k
            
    #idx = [ int(x) for x in idx ]
    #[i for i, e in enumerate(energy) if e == enumerate(energy_i)]
    #min_i = energy_i[0]
    #max_i = energy_i[-1]
    #idx_min = energy_source.index(min_i)
    #idx_max = energy_source.index(max_i)
    
    #n = data[idx,1]
    #k = data[idx,2]
    #n = n + complex(0,1)*k
    return n 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    