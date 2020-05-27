#!/usr/bin/env python

import numpy as np
import pandas as pd
import math
from scipy import optimize
import os


class OscilatoryShear():
    """
    This class method compute the properties of the output results by the LAMMPS
    
    Args:
    Path to the stored Lammps files

    """
    def __init__(self,path_to_LAMMPS):
        """ Constructor for the call to assign memebrs of the class"""
        self.path_to_LAMMPS = path_to_LAMMPS


    def shear_modulus(self,freq,gamma = 0.001, period = 400):
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
            
            # non-dimentional time 
            t=np.linspace(0,2.0*math.pi,period) 

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
            w = 1/float(freq)   # storing the frequency omega(w)
        
        return w, G1, G2, sigma_average, t


    def particle_movement(self,particle):
        
        """ 
        A function to save particle location vs time in a csv file
        
        Args:
            particle: particle id eg 1., 2.,
            
        output:
            csv file saved in the courent directory   
        """
        
        # saving all particles coordinates
        particle_coord = os.listdir(self.path_to_LAMMPS +'/coord')
        
        # number of files in the coord folder
        file_count = len(particle_coord)
        
        # initializing the simulation time and particle coordinate
        simulation_time = np.zeros(file_count)
        coordiates = np.zeros([file_count,5])
        
        ii = 0
        for coord in particle_coord:
        
            # Read each file in the coord folder
            file = self.path_to_LAMMPS +'/coord/' + coord
            try:
                fp = open(file)
                array = fp.readlines()
                ar = np.zeros([1,5])
                simulation_time[ii] = [float(j) for j in array[1].split()][0]
                for line in array:
                    if len(line.split()) == 5:
                        ar[0,:] = [float(j) for j in line.split()]
                        if ar[0,0] == particle:
                            coordiates[ii,:]=ar[0,:]
        
                ii+=1
        
            finally:
                fp.close()  
        
        df_coordiates = pd.DataFrame(coordiates, columns=["id","r","x","y","z"])
        df_coordiates['time'] = simulation_time
        df_coordiates = df_coordiates[["time","id","r","x","y","z"]]
        
        return df_coordiates


    def pairwise_forces(self,particles):
        
        """ 
        A function to save pair wise forces bewteen two particles
    
        Args:
            particle: pair particles id eg. 12, 13,
            
        output:
            csv file saved in the courent directory   
        """
        
        # saving all particles forces
        force_files = os.listdir(self.path_to_LAMMPS +'force' + str(particles) +'/')
        
        # number of files in the coord folder
        file_count = len(force_files)
        
        pathes = self.path_to_LAMMPS + 'force' + str(particles) +'/'
        
        #initializaing the simualtion time and forces in each frc.* file    
        simulation_time = np.zeros(file_count)
        force1 = np.zeros([file_count,4])
        force2 = np.zeros([file_count,4])
    
    
        ii = 0
        
        for force in force_files:
            
            file  = self.path_to_LAMMPS + 'force' + str(particles) +'/' + force
            # Check if the file is not empt
            if os.path.exists(file): 
                try:
                    fp = open(file)
            
                    array = fp.readlines()
                    
                    simulation_time[ii] = [float(j) for j in array[1].split()][0]
                    force1[ii,:] = [float(j) for j in array[9].split()]
                    force2[ii,:] = [float(j) for j in array[10].split()]
                    
                    ii +=1
                
                finally:
                    fp.close()
            else:
                print(file, "cant be open")
    
            
        # Construicting two dataframes for forces in frc.* file 
        df1 = pd.DataFrame(force1, columns=["id","fx","fy","fz"])
        df1['time'] = simulation_time
        df_force1 = df1[["time","id","fx","fy","fz"]]
        
        df2 = pd.DataFrame(force2, columns=["id","fx","fy","fz"])
        df2['time'] = simulation_time
        df_force2 = df2[["time","id","fx","fy","fz"]]        
        
        
        return df_force1, df_force2
    
    def plot_coord(self,particle):
        
        
        if not os.path.exists(self.path + 'Figures'):
            os.makedirs(self.path + 'Figures')
            
        df = self.df_coords(particle)
        df_col = df.columns
        i = 0
        for col in df_col[3:]:
            
            f,ax = plt.subplots(figsize=(12,8))
            ax.plot(df['time'],df[col],label = col ,color = 'C' + str(i))
            ax.set_xlabel('$\~t$',fontsize = 14)
            ax.set_ylabel('Movement in '+col+' direction',fontsize = 14)
            if col == 'z':
                plt.ylim([-1,1])
            plt.savefig( self.path + 'Figures/'+col+'_dir_P'+str(particle)+".tif")
            i += 1
    
        
        
    def plot_force(self,particles):
        if not os.path.exists(self.path + 'Figures'):
            os.makedirs(self.path + 'Figures')
            
        df = self.df_force(particles)
        df_col = df.columns
        i = 0
        for col in df_col[1:]:
            
            f,ax = plt.subplots(figsize=(12,8))
            ax.plot(df['time'],df[col],label = col,color = 'C' + str(i))
            ax.set_xlabel('$\~t$',fontsize = 14)
            ax.set_ylabel('Force '+col, fontsize = 14)
            
            plt.savefig( self.path + 'Figures/'+col+'_P'+str(int(df['id'][0]))+".tif")
            i += 1   
        
        
        
    
        
        
        
        
        
        
        