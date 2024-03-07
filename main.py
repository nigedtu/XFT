"""
X-ray Fitting Tool (XFT)
Author: Nis C. Gellert
Technichal University of Denmark, DTU Space
Last updated: 19/09/2023
"""
import os

os.chdir(r"C:\Users\Your\Directory\for\XFT")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from differential_evolution import de, visualize_fit
from fresnel import fresnel
from unit_conversion import wavelength2energy, energy2wavelength
from optical_constants import oc
from figure_of_merit import fom
from save_output import save
from datetime import date
import time

os.chdir(r"C:\Your\Directory\for\usert_input")
from user_input_xft_evaluation import * # Choose desired user input 

theta_m = theta_i + ang_off # Theta model includes angle off set
energy_m = energy_i + eng_off # Energy model includes energy off set

if np.isscalar(energy_m) and xrange[1]>xrange[0] : # Energy in eV, angles in deg
    idx_range = np.where( (theta_m>= xrange[0]) &  (theta_m <= xrange[1]))
    theta_m = theta_m[idx_range]
    data = data[idx_range]
    angleScan = True
elif np.isscalar(theta_m) and xrange[1]>xrange[0] : # Energy in eV, angles in deg
    energy_m = energy_m#*1000
    idx_range = np.where( (energy_m>= xrange[0]) &  (energy_m <= xrange[1]))
    energy_m = energy_m[idx_range]
    data = data[idx_range] 
    angleScan = False
    
lambda_i = energy2wavelength(energy_m)

# Create model_space input for d7e()
number_of_layers = len(sample_structure) # including substrate layer
number_of_atoms =[] 
elements = []
densities = [] 
thickness = [] 
roughness = [] 

for layer in sample_structure:
        # layer =  [material, composition, z, sigma, rho]
        elements.append(layer[0])
        number_of_atoms.append(layer[1])
        thickness.append(layer[2])
        roughness.append(layer[3])
        densities.append(layer[4])
        
model_space = densities + roughness + thickness[:-1] # excluding substrate "thickness"

# plt.close('all')
print("Logaritmic fitting: %s" %(logFitting))
print("Weighting: ", weightsFitting)
print('Iterations: ',iterations)
print('Angular off-set: ',ang_off)
print("Optical constants from "+ optical_constants_dir)
start=time.time()
fitting = list(de(model_space, theta_m, lambda_i, data, elements, number_of_atoms, angleScan, a_res, 
                  mut=mutation_factor, crossp=crossover_probability,
                  pop_size=population_size, its=iterations, fitFlag=plotFit, log_FOM = logFitting,
                  weightFunc = weightsFitting, oc_source = optical_constants_dir, mut_sch = mutation_scheme))
end=time.time()
print('Time: %.2f s'%(end-start))
pop_best, best_fit_param, minFOM = fitting[-1]
rho_bf = best_fit_param[0:number_of_layers]
sigma_bf = best_fit_param[number_of_layers:2*number_of_layers]
z_bf = best_fit_param[2*number_of_layers:] 
# np.std(pop_best[:,0])
print("Densities:", rho_bf)
print("Roughness:", sigma_bf)
print("Thickness:", z_bf)

print("Std, rho pt, rho si, sigma pt, sigma sub, z pt:", [np.std(pop_best[:,0]),np.std(pop_best[:,1]),np.std(pop_best[:,2]),np.std(pop_best[:,3]),np.std(pop_best[:,4])])
print('Fitted data from filename: ',fname)

if save_as_txt[0]:   
    if angleScan:  
        n_bf = np.zeros(len(elements)+1, dtype=np.complex)+1 
        # Calculate refractive indices of the layer materials
        for j in range(len(elements)): 
            n_bf[j+1] = oc(lambda_i,rho_bf[j],number_of_atoms[j],elements[j], optical_constants_dir)                         
        reflectance_bf = fresnel(theta_m,lambda_i ,n_bf,z_bf,sigma_bf, ares = a_res)
        save(save_as_txt[1],angleScan,sample_structure, theta_m, reflectance_bf, data,rho_bf, sigma_bf, z_bf, minFOM, a_res, energy_m, logFitting, weightsFitting, ang_off, eng_off)

    else:
        n_bf = np.zeros((len(lambda_i),len(elements)+1), dtype=np.complex)+1
        for j in range(len(elements)):
            # Get optical constants
            n_bf[:,j+1] = oc(lambda_i,rho_bf[j],number_of_atoms[j],elements[j], optical_constants_dir)    
        reflectance_bf = fresnel(theta_m,lambda_i ,n_bf,z_bf,sigma_bf, ares = a_res)
        save(save_as_txt[1],angleScan,sample_structure, theta_m, reflectance_bf, data,rho_bf, sigma_bf, z_bf, minFOM, a_res, energy_m, logFitting, weightsFitting, ang_off,eng_off)
 
    
    


#%% Plot the best fit and residual

if angleScan:
    fig, (ax3, ax4) = plt.subplots(nrows=1, ncols=2,figsize=(11,4.5))
    fig.tight_layout(pad=4.0)

    #fig3, ax3 = plt.subplots()
    ax3.plot(theta_m, data, 'k.',label="Data",markersize = 3 )
    ax3.plot(theta_m, reflectance_bf ,'r-',label ='Best fit')
    #ax3.plot(theta_m, data2 ,'b-',label ='tester.txt')
    ax3.set_ylabel("Reflectance")
    ax3.set_yscale("log")
    ax3.set_xlabel("Grazing angle (deg)")
else:
    fig, (ax3, ax4) = plt.subplots(nrows=1, ncols=2,figsize=(11,4.5))
    fig.tight_layout(pad=4.0)
    ax3.plot(energy_m, data, 'k.',label="Data",markersize = 3 )
    ax3.plot(energy_m, reflectance_bf ,'r-',label ='Best fit')
    ax3.set_xlabel("Energy (keV)")       
ax3.legend()

ax4.plot(theta_m, reflectance_bf-data2,'k.-',linewidth = 0.4, label="Data - Best fit") # 
ax4.set_ylabel("Delta Reflectance")
if angleScan:
    ax4.set_xlabel("Grazing angle (deg)")
else:
    ax4.set_xlabel("Energy (keV)")       
ax4.legend()

best_FOM = fom(data, reflectance_bf, logarithmic = logFitting, weighting = weightsFitting)
print(best_FOM)



