# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 13:52:48 2020

@author: omarg

Model calibration solver

The present script contain a class with different methods to calibrate 
SPAC parameters potr,sf,Kmax,B,C and Lp in SPAC-Tree. 

Calibration can be performed using spotpy package or Ordinary least Square 
(recommended) method

Allowed calibration parameters are in spac.csv
"""
import spotpy
#from spotpy.parameter import Uniform
import numpy as np
import pandas as pd
import main_energy
from scipy.optimize import least_squares


class spotpy_setup(object):
    
    """
    Class to perform calibration of SPAC-Tree parameters potr,sf,Kmax,B,C and Lp
    included in spac.csv file.
    
    Attributes
    calib_folder: folder that contain the neccessary files for the calibration
    outputFile: format csv
                 output file containing the target variable to compare simulations
    target_output: format str
                 Evaluation variable to perform comparisons
    """
    
    
    # potr=Uniform('potr',low=-2000,high=-800)
    # sf=Uniform('sf',low=0.002,high=0.005)
    # Kmax=Uniform('Kmax',low=0.0001,high=0.001)
    # B=Uniform('B',low=100,high=1000)
    # C=Uniform('C',low=1,high=5)
    # Lp=Uniform('Lp',low=1e-11,high=9.9e-11)
    
    
    def __init__(self):
        
        self.file=open('experiment_folder_name.txt','r')
        self.calib_folder=str(self.file.read())
        self.file.close() 
        self.outputFile='out_day.csv'
        self.target_output='T (mm day-1)'
        self.params=[spotpy.parameter.Uniform('potr',low=-2000,high=-800,optguess=-1500),
                        spotpy.parameter.Uniform('sf',low=0.001,high=0.007,optguess=0.003),
                        spotpy.parameter.Uniform('alfa',low=1,high=10,optguess=2),
                        spotpy.parameter.Uniform('beta',low=10,high=100,optguess=30),
                        spotpy.parameter.Uniform('sigma',low=0.1,high=0.4,optguess=0.25),
                        spotpy.parameter.Uniform('Lp',low=6e-9,high=9e-8,optguess=6e-8)
                        ]
        #self.params=[spotpy.parameter.Uniform('Lp',low=1e-9,high=9.9e-11,optguess=2e-11)]
        
        #self.mrun=model_run()
        
    # def read_file(self):
    #     file=open('experiment_folder_name.txt','r')
    #     calib_folder=str(file.read())
    #     file.close()
    
    #     return calib_folder
        
    
    def param_generator(self):

          return spotpy.parameter.generate(self.params)
      
      
    def update_spac_spotpy(self,param_generated):
          
        legth_param_gen=np.shape(param_generated)[0]
        spac=pd.read_csv(self.calib_folder+'input/tree/spac.csv',header=None,index_col=0)
        
        if spac.shape[1]==2:
            spac.drop(columns=2,inplace=True)
        
        for i in range(0,legth_param_gen):
            spac.loc[param_generated[i][1]]=param_generated[i][0]
            
        spac.to_csv(self.calib_folder+'input/tree/spac.csv',header=False)    
     
    def update_spac(self,parameters):
        
        
        
        parameters={'potr':parameters[0],
                    'sf':parameters[1],
                    'alfa':parameters[2],
                    'beta':parameters[3],
                    'sigma':parameters[4],
                    'Lp':parameters[5]}
        #parameters={'Lp':parameters[0]}
        spac=pd.read_csv(self.calib_folder+'input/tree/spac.csv',header=None,index_col=0)
        
        if spac.shape[1]==2:
            spac.drop(columns=2,inplace=True)
        
        for key in parameters:
            spac.loc[key]=parameters[key]
            
        spac.to_csv(self.calib_folder+'input/tree/spac.csv',header=False)
        
              
    def simulation(self,vector):
        
        #calib_folder=self.read_file()
        
        
        param_generated=vector
        self.update_spac(param_generated)
        main_energy.main()
        
        
        output=pd.read_csv(self.calib_folder+'output/'+self.outputFile)
        simulations=list(np.array(output[self.target_output]))
        
        return simulations
    
    def evaluation (self):
        obs_sim=pd.read_csv(self.calib_folder+'output/observed.csv')
        observations=list(obs_sim['Obs'])          
        return observations

    def objectivefunction (self,simulation,evaluation):
        
        objectivefunction= spotpy.objectivefunctions.rmse(evaluation,simulation)      
        return objectivefunction
    
    def fun_residuals(self,parameters):
        
        print(parameters)
        
        simulations=np.array(self.simulation(parameters))
        observations=np.array(self.evaluation())
        
        return simulations-observations
    
        
            
         
            



# class model_run():
    
#     def _init_(self,outputFile='out_day.csv',target_output='T (mm day-1)'):
        
        
#         self._outputFile=outputFile
#         self._target_output=target_output
    
#     """
#     run the model and change the input parameters for spac file according 
#     to the calib_param dictionary
    
#     Inputs
#     ------
#         calib_params: dictionary
        
#         modelName: string
#         Name of the model for calibration
        
#         outputFile: csv
#         File to obtain model simulation results
        
#         target_output: string
#         variable to compare with observed data
        
#     Output
#     ------
#         simulation: array
#         Data to be evaluated 
#     """
#     def read_file(self):
#         file=open('experiment_folder_name.txt','r')
#         calib_folder=str(file.read())
#         file.close()
        
#         return calib_folder
    
#     def update_spac (self,param_generated,calib_folder):
        
#         spac=pd.read_csv(calib_folder+'input/tree/spac.csv',header=None,index_col=0)
        
#         if spac.shape[1]==2:
#             spac.drop(columns=2,inplace=True)
        
#         for i in range(0,param_generated.shape[0]):
#             spac.loc[param_generated[i][1]]=param_generated[i][0]
            
#         spac.to_csv(calib_folder+'input/tree/spac.csv',header=False)    
    
#     def _run (self):
#         main_RintOptions.main()
        
#     def read_data (self,calib_folder,outputFile='out_day.csv',target_output='T (mm day-1)'):
        
#         output=pd.read_csv(calib_folder+'output/'+outputFile)
#         simulation=list(np.array(output[target_output]))
    
#         return simulation



if __name__ == '__main__': 
    
    mcalib=spotpy_setup()
    parameters=[-1500,0.003,2,30,0.25,6e-8]
    bounds=([-2000,0.001,1,10,0.1,6e-9],[-800,0.007,10,100,0.4,9e-8])
    result=least_squares(mcalib.fun_residuals,parameters,bounds=bounds)
    
    
    # obj_func_results=[]
    # sim_results=[]
    # vector_results=[]
    # spac_calib=[]
    # niterations=1000
    # calib=spotpy_setup()
    # error=0.01
    # n=0
    # obj_func=100
    # while obj_func>error:# range(0,niterations):
      
    #     vector=calib.param_generator()
    #     sim=calib.simulation(vector)
    #     obs=calib.evaluation()
    #     obj_func=calib.objectivefunction(sim, obs)
        
    #     #result={'obj_func':obj_func,'vector':vector[0],'simulation':sim}
        
    #     obj_func_results.append(obj_func)
    #     sim_results.append(sim)
    #     vector_results.append(vector)
    #     n=n+1
    #     if n>niterations:
    #         break
    #     print(n)
        
    # min_obj=np.min(np.array(obj_func_results))
    # row=np.where(obj_func_results==min_obj)
    # row=row[0][0]
    # final_vector=vector_results[row]
    
    # for i in range(0,final_vector.shape[0]): 
    #     spac_calib.append(final_vector[i][0])
       
       # results=[]
       # spot_setup=spotpy_setup()
       # rep=3
       # timeout=60
       # parallel = "seq"
       # dbformat = "csv"
      
      # sampler=spotpy.algorithms.mc(spot_setup,parallel=parallel, dbname='SPAC-MC', dbformat=dbformat)
      # sampler.sample(rep)
      # results=results.append(sampler.getdata())
    
    
    

