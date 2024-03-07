"""
X-ray Fitting Tool (XFT)
Author: Nis C. Gellert
Technichal University of Denmark, DTU Space
Last updated: 19/09/2023
"""

import numpy as np
import matplotlib.pyplot as plt
from optical_constants import oc 
from fresnel import fresnel
from matplotlib.widgets import Button
from figure_of_merit import fom
from matplotlib.ticker import MaxNLocator
from unit_conversion import wavelength2energy, energy2wavelength
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
import time


stopIterations = False
def de(model_space, theta, wavelength, data, elements, number_of_atoms, angleScan,a_res,
       mut=0.8, crossp=0.7, pop_size=30, its=150, fitFlag = False, log_FOM = False, weightFunc = 'equal', oc_source = "oc_source//CXRO//", mut_sch = "Rand/1"):
    # Initializes the individuals of the population
    number_of_layers = len(elements)
    dimensions = len(model_space) # no. of model parameters
    pop = np.random.rand(pop_size, dimensions) # Normalized population between 0 and 1, pop_size x dimensions
    min_b, max_b = np.asarray(model_space).T
    diff = np.fabs(min_b - max_b)
    pop_denorm = min_b + pop * diff # denormalize by converting each component from [0, 1] to [min, max]
    
    # Check if there are coupled parameters, and the couple them
    for i in range(len(model_space)):
        for j in range(i+1,len(model_space)):
            if (model_space [j] is model_space[i]) == True:
                coupled_parameters = [i,j]
                pop[:,j] = pop[:,i]
                pop_denorm[:,j] = pop_denorm[:,i]

    if angleScan:
        n = np.zeros(len(elements)+1, dtype=np.complex)+1  # + 1 due to top vacuum
        reflectance = np.zeros([len(theta),pop_size])
        print("Fitting angle scan")
        print(" ")
    else:
        # For energy scans, n is computed for each wavelength and each layer material
        n = np.zeros((len(wavelength),len(elements)+1), dtype=np.complex)+1  # + 1 due to top vacuum
        reflectance = np.zeros([len(wavelength),pop_size])
        print("Fitting energy scan")
        print(" ")
    FOM =  np.zeros(pop_size)
    
   
    # Evaluation of the initial population
    # fitness = np.asarray([fobj(ind) for ind in pop_denorm])
    for i in range(pop_size):
        # We only use the denormalized population when evaluating reflectance and corresponding FOM
        rho = pop_denorm[i,range(0,number_of_layers)] # position 0-3 (3 inclusive)
        sigma = pop_denorm[i,range(number_of_layers, 2*number_of_layers)] # position 4-7
        z = pop_denorm[i,range(2*number_of_layers,len(model_space) )] # position 8-10
        
        if angleScan:  
            # Calculate refractive indices of the layer materials
            for j in range(len(elements)): 
                n[j+1] = oc(wavelength,rho[j],number_of_atoms[j],elements[j], oc_source)
            reflectance[:,i] = fresnel(theta,wavelength,n,z,sigma, ares = a_res)        
        else:     
            for j in range(len(elements)):
                # Get optical constants
                n[:,j+1] = oc(wavelength,rho[j],number_of_atoms[j],elements[j], oc_source)    
            reflectance[:,i] = fresnel(theta,wavelength,n,z,sigma, ares = a_res)
        
        FOM[i] = fom(data,reflectance[:,i], log_FOM, weightFunc) # figure of merit of each individual in the population, 1D array
    
    # Locate the best FOM from the initial population
    best_idx = np.argmin(FOM)
    best_denorm = pop_denorm[best_idx]
    reflectance_best = reflectance[:,best_idx] 
    fom_evolution = [] 
    iteration = [] 
    delta_FOM = []
    pop_evolution = np.zeros(dimensions) 

    for i in range(its):
        global stopIterations # Used for the stop&save button 
        if stopIterations == True: 
            break
        
        for j in range(pop_size):
            
           
            # List of indices idxs of the population excluding the current j
            idxs = [idx for idx in range(pop_size) if idx != j] 
            if (mut_sch == "Rand/1"):
                # choose 3 other vectors/individuals from the normalized population, excluding the current pop[j], without replacement
                a, b, c = pop[np.random.choice(idxs, 3, replace = False)] # "Rand/1"
            
                # create a mutation vector from the 3 vectors
                mutant = np.clip(a + mut * (b - c), 0, 1) # values outside the interval [0, 1] are clipped to the interval edges
            
            if (mut_sch == "Best/1"):
                a = pop[best_idx]
                b, c = pop[np.random.choice(idxs, 2, replace = False)] # "Best/1"
                mutant = np.clip(a + mut * (b - c), 0, 1)
                
                
            # Recombination 
            cross_points = np.random.rand(dimensions) < crossp 
            # boolean array with True elements when the random values are lower than the crossp probability    
            if not np.any(cross_points):
                # if there aren't any True elements in cross_point, select a random index and set True 
                cross_points[np.random.randint(0, dimensions)] = True
            
            # Create trial: where cross_points is True, trial takes on the value from mutant - otherwise use pop[j]
            trial = np.where(cross_points, mutant, pop[j])
            trial_denorm = min_b + trial * diff
            
            # Couple parameters in trial-vector
            # This method is not the best, since you should not need to find the couple parameters again. It is done above
            for ii in range(len(model_space)):
                for jj in range(ii+1,len(model_space)):
                    if (model_space [jj] is model_space[ii]) == True:
                        trial[jj] = trial[ii]
                        trial_denorm[jj] = trial_denorm[ii]
                        
            # Evaluate the reflectance and FOM with the trial vector
            rho_old = rho # new here
            rho = trial_denorm[range(0,number_of_layers)] # position 0-3 (3 inclusive)
            sigma = trial_denorm[range(number_of_layers, 2*number_of_layers)] # position 4-7
            z = trial_denorm[range(2*number_of_layers,len(model_space) )] # position 8-10
            
            # NOTE 20210827: YOOO! Is there not at error here??! The trial reflectance should not be save before comparing??? 
            if angleScan:  
            # Calculate refractive indices of the layer materials
                for k in range(len(elements)):
                # Get optical constants
                    if rho_old[k] != rho[k]: # Reducing computation time if new density is equal to previous
                        n[k+1] = oc(wavelength,rho[k],number_of_atoms[k],elements[k], oc_source)                              
                reflectance_trial = fresnel(theta,wavelength,n,z,sigma, ares = a_res)
            else:     
                for k in range(len(elements)):
                    # Get optical constants
                    if rho_old[k] != rho[k]: # Reducing computation time if new density is equal to previous
                        n[:,k+1] = oc(wavelength,rho[k],number_of_atoms[k],elements[k], oc_source)       
                reflectance_trial  = fresnel(theta,wavelength,n,z,sigma, ares = a_res)
            
            f_temp = fom(data,reflectance_trial,log_FOM,weightFunc)
                        
            # Compare with the current FOM of the j'th individual
            if f_temp < FOM[j]: 
                FOM[j] = f_temp
                pop[j] = trial
                reflectance[:,j] = reflectance_trial
                if f_temp < FOM[best_idx]: 
                    # compare with the best FOM among all individuals in the population
                    best_idx = j
                    best_denorm = trial_denorm
                    #print(min_b + trial * diff)
                    #reflectance_best = reflectance[:,best_idx]  # the best reflectance of the 
                    reflectance_best = reflectance_trial # the best reflectance  
                    
                    
        fom_evolution.append(FOM[best_idx])
        iteration.append(i)
        avg10_FOM = np.average(fom_evolution[-10:]) # take the average fom for the last 10 iterations
        delta_FOM.append(avg10_FOM -  FOM[best_idx])
        pop_evolution=np.vstack([pop_evolution,np.asmatrix(np.mean(pop,axis=0))]) # mean or std# HERE
        if fitFlag == True : 
            if i == 0:
                plt.ion()
                fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3,figsize=(18,6))
                fig.tight_layout(pad=4.0)
                
                
            # Plot the current best-fit
            if angleScan:
                #print("no")
                visualize_fit(theta, data, reflectance_best, iteration, fom_evolution,delta_FOM,ax1,ax2,ax3,angleScan, pop_evolution,pop)
                fig.canvas.draw()
            else:
                #print("yay")
                x = wavelength2energy(wavelength)/1000
                #x = wavelength
                visualize_fit(x, data, reflectance_best, iteration, fom_evolution,delta_FOM,ax1,ax2,ax3,angleScan, pop_evolution,pop)
                fig.canvas.draw()
       
        yield min_b + pop * diff, best_denorm, FOM[best_idx] 
        # for each iteration - yields the entire denormalized population (all model parameters), 
        # the best individual of the population, and the correspoding FOM 


# Stop plot and save 
def stop(self): 
    global stopIterations
    stopIterations = True
    print("Script stopped:")
    
def sliding_std_window(elements):
    std = np.zeros(len(elements))
    for i in range(len(elements)):
        std[i] = np.std(elements[0:i+1])
    return(std)   

def visualize_fit(x, data, reflectance, iter_no, FOM, delta_FOM ,ax1,ax2,ax3,angleScan, pars_norm,pop):
    fts = 15
    ax1.clear()
    ax1.plot(x, data,'k.',label="Data",markersize = 5, linewidth = 1)  
    ax1.set_ylabel("Reflectance", fontsize=fts)
    ax1.plot(x, reflectance,'r-',linewidth = 1.5,label = 'Best model') # plot best reflectance for each iteration
    ax1.set_title('Iteration %i: FOM = %.6e' %(iter_no[-1],FOM[-1]))
    ax1.legend(handlelength=3, fontsize=fts)
    ax1.grid()
    if angleScan:
        ax1.set_yscale("log")
        ax1.set_xlabel("Grazing incidence angle (deg)", fontsize=fts)
    else:
       #ax1.set_yscale("log")
       ax1.set_xlabel("Energy (keV)") 
    
    ax2.plot(iter_no, FOM,'b.-',linewidth = 0.4) # label = 'Iteration %i: FOM = %.5e' %(i,fit_fom),
    ax2.set_xlabel("Iteration no.", fontsize=fts)
    ax2.set_ylabel("FOM", fontsize=fts)
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.legend(["FOM"],handlelength=3, fontsize=fts)
    # ax2.legend()
    ax2.grid()
    ax2.set_yscale("log")
    
    
    
    rho_pt = pars_norm[:,0]
    rho_si = pars_norm[:,1]
    sig_pt = pars_norm[:,2]
    sig_sub = pars_norm[:,3]
    z_pt = pars_norm[:,4]
    #eb_z_pt = np.std(pop[:,5])
    
    eb_rho_pt = + np.std(pop[:,0])
    eb_sig_pt = + np.std(pop[:,2])
    eb_sig_sub = + np.std(pop[:,3])
    eb_z_pt = + np.std(pop[:,4])
    #print(pop[:,4])
    #print(pars_norm)
    iter_no = np.arange(0,len(z_pt),1)
    #print(z_pt)
    
    ax3.plot(sig_sub,marker='s',linestyle= '-',label ='Sub roughness',color='r',linewidth = 1.0) 
    ax3.errorbar(iter_no, sig_sub, eb_sig_sub,  marker = '^', ms = 0.0, color = 'red', linewidth = 0., elinewidth = 0.5)
    
    ax3.plot(sig_pt,marker='s',linestyle= '-',label ='Pt roughness',color='b',linewidth = 1.0) 
    ax3.errorbar(iter_no, sig_pt, eb_sig_pt,  marker = 's', ms = 0.0, color = 'blue', linewidth = 0., elinewidth = 0.5)
    
    ax3.plot(z_pt,marker='s',linestyle= '-',label ='Pt thickness',color='k',linewidth = 1.0) 
    ax3.errorbar(iter_no, z_pt, eb_z_pt,  marker = '*', ms = 0.0, color = 'black', linewidth = 0., elinewidth = 0.8)
    
    ax3.plot(rho_pt,marker='s',linestyle= '-',label ='Pt density',color='g',linewidth = 1.0) 
    ax3.errorbar(iter_no, rho_pt, eb_rho_pt,  marker = '*', ms = 0.0, color = 'green', linewidth = 0., elinewidth = 0.8)
    
    #sliding_std_window(z_pt)
    ax3.grid()
    ax3.set_ylim(0.0,1)
    #ax3.set_yscale("log")

   

    #ax3.plot(iter_no, np.array(delta_FOM)/np.array(FOM),'g.-',linewidth = 0.4) # label = 'Iteration %i: FOM = %.5e' %(i,fit_fom),
    ax3.set_xlabel("Iteration no.", fontsize=fts)
    #ax3.legend(handlelength=5)
    ax3.legend(("Sub roughness","Pt roughness","Pt thickness","Pt density" ),handlelength=3, fontsize=fts)
    #ax3.set_xlabel("Iteration no.", fontsize=fts)
    ax3.set_ylabel("Parameter evolution", fontsize=fts)
    #ax3.set_ylabel("delta_FOM/FOM")
    #ax3.xaxis.set_major_locator(MaxNLocator(integer=True))
   # if np.any(np.array(delta_FOM)/np.array(FOM) > 0 ):
   #     ax3.set_yscale("log")
        
    #ax_button = fig.add_axes([0.45, 0.05, 0.08, 0.05])
    
    plt.pause(0.01)
    #ax_Button = plt.axes([0.45, 0.05, 0.08, 0.05])  #posx, posy, width, height
    #stop_Button = Button(ax_button, 'Stop & Save', color='white', hovercolor='red')
    
    #stop_Button.on_clicked(stop) # When clicked, call stop function
    plt.show() # show on figure 
    
    


