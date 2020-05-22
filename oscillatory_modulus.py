
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from scipy import optimize
import gc
import os


class oscilatory_modulus():

    """
    This class method compute the properties of the output results by the LAMMPS
    
    Args:
    It requirs the path to the stored files

    """
    def __int__(self,path_to_LAMMPS,gamma = 0.001, period = 400):
        """ Constructor for the call to assign memebrs of the class"""
        self.path_to_LAMMPS = path_to_LAMMPS
        self.period = period
        self.gamma = gamma


    def shear_modulus(self,freq):
        """ 
        A function to compute shear modulus from LAMMPS softwar

        Args:
        gamma: Oscillatory shear amplitude default 0.001
        period: number of data points in each period defualt 400

        output:
        returns frequency(w), storage modulus (G1) and loss modulus G2     
        """

        # Check if the file is not empty
        if os.path.exists(self.path_to_LAMMPS + '/fix.1.tot'): 
            
            # Loading the dataset using numpy loadtxt
            data = np.loadtxt(self.path_to_LAMMPS  + '/fix.1.tot')

            # Defining the period and gamma
            gamma = self.gamma
            period = self.period      

            #number of time steps in each ferequency period
            n_period = int((len(data)+1)/period) # Number of periods
            start = int(n_period / 2);   # Start of the period for analysis
            end = n_period; # End of the period for analysis
    
            # Initializing the average stress
            stress_average=np.linspace(0.0,0.0,period)
            for i in range(period):
                for j in range(start,end):
                    stress_average[i]+=data[j*period+i,7]
                    
                stress_average[i]=stress_average[i]/(end-start)
                
            stress_average=stress_average-stress_average.mean()
            a0=stress_average.max()
            t=np.linspace(0,2.0*math.pi,period) 
            # Define a function for parameter estimation
            def curve_fit_func(x,a,phi):
                return a*np.sin(x+phi)
                            
            params, params_covariance=optimize.curve_fit(curve_fit_func, t, 
                                                         stress_average,p0=[a0,0])
                        
            sigma_max, phi=params # Storing the resuls for sigma = sigma_max*sin(wt+pji)

            G1 = sigma_max*math.cos(phi)/gamma # computing G'
            G2 = sigma_max*math.sin(phi)/gamma # computing G"
        
            w = 1/float(freq)
        
        return w, G1, G2

