#!/usr/bin/env python

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
            frequency(w)
            storage modulus(G1) 
            loss modulus G2     
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
            sigma_average=np.linspace(0.0,0.0,period)
            for i in range(period):
                for j in range(start,end):
                    sigma_average[i]+=data[j*period+i,7]
                    
                sigma_average[i]=sigma_average[i]/(end-start)
                
            sigma_average=sigma_average-sigma_average.mean()
            a0=sigma_average.max()
            t=np.linspace(0,2.0*math.pi,period) 

            # Define a function for parameter estimation
            def curve_fit_func(x,a,phi):
                return a*np.sin(x+phi)

            # Fitting the stress average data to y = sigma_max*sin(wt+phi)           
            params, params_covariance=optimize.curve_fit(curve_fit_func, t, 
                                                         sigma_average,p0=[a0,0])
                        
            sigma_max, phi=params # Storing the resuls for sigma = sigma_max*sin(wt+pji)

            G1 = sigma_max*math.cos(phi)/gamma # computing G'
            G2 = sigma_max*math.sin(phi)/gamma # computing G"
        
            w = 1/float(freq)
        
        return w, G1, G2, sigma_average, t


        def particle_movement(self):
            
            i = 0
            ii = 0

            time = np.zeros(file_count)

            disp24 = np.zeros([file_count,5])

            for i in range(file_count):

                cnt = i*1000
                file = pathes + "dump." + str(cnt)

                try:
                    fp = open(file)

                    array = fp.readlines()
                    ar = np.zeros([1,5])
                    time[ii] = [float(j) for j in array[1].split()][0]
                    for line in array:
                        if len(line.split()) == 5:
                            ar[0,:] = [float(j) for j in line.split()]
                            if ar[0,0] == P[I]:
                                disp24[ii,:]=ar[0,:]

                    ii+=1

                finally:
                    fp.close()  

            df_disp = pd.DataFrame(disp24, columns=["id","r","x","y","z"])
            df_disp['time'] = time
            df_disp = df_disp[["time","id","r","x","y","z"]]

            df_disp.to_csv(PI +".csv")