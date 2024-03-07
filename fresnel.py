"""
X-ray Fitting Tool (XFT)
Author: Nis C. Gellert
Technichal University of Denmark, DTU Space
Last updated: 19/09/2023

Input units: 
    theta: incident angle, degree (either a single value or an array)
    lambda_i: incident wavelength, nm (either a single value or an array)
    refractive index n (real and complex part) for each layer material
    thickness z, nm of each layer material
    roughness sigma, nm  of each layer material (including substrate roughness)

Ouput:
    reflectance 

"""


import numpy as np
import matplotlib.pyplot as plt
from optical_constants import oc

def fresnel(theta, lambda_i, n, z, sigma, det_sensitivity = 1, ares = 0, polarization = 0.00, mfc = 1):
    lambda_i = lambda_i*1e-9 # convert to m
    
    if np.isscalar(theta):
        theta = np.array([theta])
        reflectance = np.zeros((len(lambda_i))) 
        energyScan = True
    else:
        lambda_i = np.array([lambda_i])
        reflectance = np.zeros((len(theta))) 
        energyScan = False

        
    theta = theta*np.pi/180. # degrees to radians
    theta_i = np.zeros((len(theta), len(sigma)+1), dtype=np.complex)
    theta_i[:, 0] = np.pi/2.- np.array(theta) 
    theta_t = np.zeros((len(theta), len(sigma)+1), dtype=np.complex)
        
    # interface profile
    s = np.zeros((len(theta), len(sigma)), dtype=np.complex)
    w = np.zeros((len(theta), len(sigma)), dtype=np.complex) 
    
    # reflection and transmission
    r_s = np.zeros((len(theta), len(sigma)), dtype=np.complex)
    r_p = np.zeros((len(theta), len(sigma)), dtype=np.complex) 
    
    t_s = np.zeros((len(theta), len(sigma)), dtype=np.complex)
    t_p = np.zeros((len(theta), len(sigma)), dtype=np.complex) 
    
    mfc_r_s = np.zeros((len(theta), len(sigma)), dtype=np.complex)
    mfc_r_p = np.zeros((len(theta), len(sigma)), dtype=np.complex) 
    
    
    z = np.append(z , [0]) # substrate "thickness" added as 0 
    z = np.array(z)*1.0e-9 # convert from nm to meter
    sigma=np.array(sigma)*1.0e-9; # convert from nm to meter
 
    net_r_s = np.zeros((len(theta_i), len(sigma)+1),dtype=np.complex)
    net_r_p = np.zeros((len(theta_i), len(sigma)+1),dtype=np.complex)
    beta = np.zeros((len(theta_i), len(sigma)),dtype=np.complex)
    

    for k in range(len(lambda_i)): # looping over all wavelengths
        
        if energyScan :
            n_temp = n[k,:]
            lambda_temp = lambda_i[k]
        else:
           n_temp = n
           lambda_temp = lambda_i
        
            
        for i in range(len(sigma)):
            j = i+1
            #Calculate theta_i and theta_t for every interface 
            theta_t[:,i] = np.arcsin(n_temp[i] * np.sin(theta_i[:,i])/n_temp[j])
            theta_i[:,j] = theta_t[:,i]
           
            # interface profiles 
            s[:,i] = 4.* np.pi * np.sqrt(np.cos(theta_i[:,i]) * np.cos(theta_t[:,i]))/lambda_temp
            w[:,i] = np.exp((-s[:,i]**2. * sigma[i]**2.)/2.)
            
            # s polarization
            r_s[:,i] = (n_temp[i] * np.cos(theta_i[:,i]) - n_temp[j] * np.cos(theta_t[:,i])) / (n_temp[i] * np.cos(theta_i[:,i]) + n_temp[j] * np.cos(theta_t[:,i]))  
            t_s[:,i] = 2. * n_temp[i] * np.cos(theta_i[:,i]) / (n_temp[i] * np.cos(theta_i[:,i]) + n_temp[j] * np.cos(theta_t[:,i]))  
            
            # p polarization
            r_p[:,i] = (n_temp[i] * np.cos(theta_t[:,i]) - n_temp[j] * np.cos(theta_i[:,i])) / (n_temp[i] * np.cos(theta_t[:,i]) + n_temp[j] * np.cos(theta_i[:,i]))  
            t_p[:,i] = 2. * n_temp[i] * np.cos(theta_i[:,i]) / (n_temp[i] * np.cos(theta_t[:,i]) + n_temp[j] * np.cos(theta_i[:,i]))  
    
            # Modified Fresnel coefficients
            mfc_r_s[:,i] = r_s[:,i] * w[:,i]
            mfc_r_p[:,i] = r_p[:,i] * w[:,i]
            
    
        
        for i in range(len(sigma)-1, -1, -1):
            j = i+1
            beta[:, i] = 2.*np.pi*z[i]*n_temp[j]*np.cos(theta_i[:, j])/lambda_temp
            
            net_r_s[:, i] = (mfc_r_s[:, i] + net_r_s[:, j]*np.exp(np.complex(0, 2)*beta[:, i]))/(1 + mfc_r_s[:, i]*net_r_s[:, j]*np.exp(complex(0, 2)*beta[:, i]))
            
            net_r_p[:, i] = (mfc_r_p[:, i] + net_r_p[:, j]*np.exp(np.complex(0, 2)*beta[:, i]))/(1 + mfc_r_p[:, i]*net_r_p[:, j]*np.exp(np.complex(0, 2)*beta[:, i]))
            
        R_s = np.abs(net_r_s)**2.
        R_p = np.abs(net_r_p)**2.
        R = (R_s*det_sensitivity*(1+polarization) + R_p*(1-polarization))/(polarization*(det_sensitivity-1) + (det_sensitivity+1))
        
        
        if energyScan:
            reflectance[k] = R[0, 0]
        else:
           reflectance = R[:, 0]
        
        
        if ares != 0: # https://github.com/genx-dev/genx/blob/master/genx/models/lib/instrument.py
            ares = ares*np.pi/180. # degrees to radians
            #theta = theta * 180/np.pi
            
            q_step = theta[1]-theta[0]
            range_x = 5
            I = reflectance
            #resvector = np.arange(-range_x*ares,range_x*ares+q_step,q_step) # from linc, though this form is not gaussian
            resvector = np.arange(-range_x*ares,range_x*ares+q_step/2,q_step) # 
            weight = 1/np.sqrt(2*np.pi)/ares*np.exp(-(resvector)**2/(ares)**2/2)
            reflectance = np.convolve(np.r_[np.ones(resvector.shape)*I[0], I, np.ones(resvector.shape)\
                                            *I[-1]], weight/weight.sum(), mode=1)[resvector.shape[0]:-resvector.shape[0]]
           
            
    return reflectance


