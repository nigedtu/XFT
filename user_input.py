"""
X-ray Fitting Tool (XFT)
Author: Nis C. Gellert
Technichal University of Denmark, DTU Space
Last updated: 19/09/2023

===============================================================
The following is an example of how to fit the 
8.048 keV XRR data of a 30 nm Pt single layer created in IMD. 

"""

import numpy as np
import os
os.chdir(r"C:\Users\nige\DTU Space\PhD\XRR_Fitting_tool")

fname = "Pt_mf2_IMD.txt" # Name of data file you want to fit
fname_dir = "C:/Users/..Directory../XFT/Data/" # Directory of that file
save_as_txt = [False, fname_dir+fname[:-4]+"_XFT.txt"] # Save fitted data as txt 
data = abs(np.loadtxt( fname_dir+fname)[:,1]) # Normalized measured reflectance

'''# Define Instrument parameters: '''
theta_i = np.loadtxt(fname_dir+fname)[:,0] # [deg] Incident angle 
xrange = [0.25, 3.0] # min_x and max_x
energy_i = 8048# [eV]. Reflectometer energy. Array or variable. 
a_res = .00#7 # [degree] Instrument Angular resolution Note: Not sure if correct
e_res = 0.0 # [m]. Instrumental spectral resolution. Note: Not implemented
polarization_incident = 0.0 # Incident Polarization [-1:1].s
polarization_analyzer = 1.0 # Polarization analyzer sensitivity
ang_off = 0#(0.015) # [degree] Angular offset of data
eng_off = 0.0 # [eV] Energy offset 


'''# Define Differential evolution parameters: '''
mutation_factor = 0.7 # [0.1:2.0]
crossover_probability = 0.7 # [0:1]
population_size = 30 # [10:50]
iterations = 161 # [25 - 550] 
mutation_scheme = "Rand/1" # "Rand/1" or "Best/1"
plotFit = True # plots the best fit in every iteration
logFitting = False
weightsFitting = 'equal' # 'equal' or 'statistical'


'''# Define Model parameters: '''
optical_constants_dir = "oc_source/CXRO/" # [CXRO, NIST, LLNL]

#           [material,      composition,    z,      sigma,      rho]
Si_layer = [["Si"],    [1],  [8.0, 11.5],     [0.1, 1.4],[2.32, 2.32]] # Si = 2.33
NiV_layer = [["Ni","V"],      [93,7],  [8.0, 11.5], [0.5, 2.0], [7.9, 8.79]] # NiV = 8.79 - 8.9, Ni=8.9
Pt_layer = [["Pt"],    [1],  [3,40],  [0.1,2.0], [18.0,24.0]] # 21.45
SiO2_layer = [["Si","O"],    [1,2], [0.2,0.2],     [0.3, 0.3],[2.65, 2.65]] 

#substrate = [["Si"],    [1],  np.NaN,     sio2_layer[3],[2.33, 2.33]] 
substrate = [["Si"],    [1],  np.NaN,     [0.1, 0.9],[2.33, 2.33]] 

sample_structure = [Pt_layer,
                    substrate]



