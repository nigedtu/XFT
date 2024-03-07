"""
X-ray Fitting Tool (XFT)
Author: Nis C. Gellert
Technichal University of Denmark, DTU Space
Last updated: 19/09/2023
"""
import numpy as np

def wavelength2energy(wavelength):
    # wavelength is in nm - conversion to energy in eV 
    eV = 1.60217662*1e-19 # [J] per eV
    h = 6.62607004*1e-34 # [m^2*kg/s]
    c = 299792458 # [m/s]
    
    E = h*c/(wavelength*1e-9)
    energy = E/eV 
    return energy
    
def energy2wavelength(energy):
    # energy is in eV - conversion to wavelength in nm 
    eV = 1.60217662*1e-19 # [J] per eV
    h = 6.62607004*1e-34 # [m^2*kg/s]
    c = 299792458 # [m/s]
    
    E = energy*eV # converting from [eV] to [J]
    wavelength = h*c/E*1e9 # [nm] Instrument wavelength / energy
    return wavelength



def thickness2energy(thickness,theta):
    # thickness (nm), energy (keV), wavelength (nm), theta (degree).   
    energy = wavelength2energy(2*thickness*np.sin(np.deg2rad(theta)))/1000
    return energy

def energy2thickness(energy,theta):
    # thickness (nm), energy (keV), wavelength (nm), theta (degree).   
    thickness = wavelength2energy(energy*1000)/(2*np.sin(np.deg2rad(theta)))
    return thickness


def thickness2energy_cor(thickness,theta,delta,wavelength,n):
    # thickness (nm), energy (keV), wavelength (nm), theta (degree).   
    energy = wavelength2energy(2*thickness*np.sin(np.deg2rad(theta))*(1-((4*thickness**2*delta)/(n**2*wavelength**2))))/1000
    return energy