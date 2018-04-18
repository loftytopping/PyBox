##########################################################################################
#                                                                                        #
#    Example gas phase model. This takes an equation file, given in the KPP format,      #
#    and then sets up an ODE solver where initial concentrations of any specie can       #
#    be set. It also relies on some pre-defined rate coefficients and photolysis rates   #
#    taken from the MCM. These are explicitly written into the relevant Fortran modules  #
#    in the file Parse_eqn_file.write_rate_file_fortran() which is provided with a 
#    Fortran syntax version of pre-defined rates in  
#                                                                                        #
#    Mixed Python - Fortran version. This version uses the f2py module to re-write       #
#    the RHS calculations to exploit multi-core shared/dsitributed memory machine        #
#                                                                                        #
#                                                                                        #
#    Copyright (C) 2017  David Topping : david.topping@manchester.ac.uk                  #
#                                      : davetopp80@gmail.com                            #
#    Personal website: davetoppingsci.com                                                #
#                                                                                        #
#    This program does not yet have a license, meaning the deault copyright law applies. #
#    Only users who have access to the private repository that holds this file may       #
#    use it or develop it, but may not distribute it without explicit permission.        #
#                                                                                        #
#                                                                                        #
##########################################################################################

# VERSION 0.9
# - Not fit for public release.

# Developed using the Anaconda Python distribution and with the Assimulo ODE solver suite: http://www.jmodelica.org/assimulo
# Assimulo gives a Python front-end acces to stiff solvers so this can be merged with aerosol developments in a seperate branch.
# The Assimulo package does not allow extra argument passing, and thus defines the structure of the code. 
# In the import statements, all files developed specifically for this project as marked [•]

import numpy 
import Parse_eqn_file # [•] Needed to parse the .eqn file, name given in this file
import rate_coeff_conversion # [•] Converts standard text rate coefficients into Python/Fortran
import MCM_constants # [•] holds more pre-defined rate coefficients and photolysis rates not provided in .eqn file. 
import collections
import pdb
from datetime import datetime
import time
from ODE_solver import run_simulation # [•] Contains routines to run ODE solver
import os
import pickle
            
# Start of the main body of code
if __name__=='__main__':
   
    #-------------------------------------------------------------------------------------

    #1)Define starting ambient conditions
    temp=288.15
    RH=0.5 #RH/100%
    #Define a start time 
    hour_of_day=12.0
    start_time=hour_of_day*60*60 # seconds, used as t0 in solver
    #2)Generate constants used in rate of reaction calculations
    #Convert RH to concentration of water vapour molecules [this will change when in Parcel model mode]
    temp_celsius=temp-273.15
    # Saturation VP of water vapour, to get concentration of H20
    Psat=610.78*numpy.exp((temp_celsius/(temp_celsius+238.3))*17.2694)
    Pw=RH*Psat
    Wconc=0.002166*(Pw/(temp_celsius+273.16)) #kg/m3
    #Convert from m3 to cm3
    Wconc=Wconc*1.0e-6
    #Convert from kg to molecules/cc
    H2O=Wconc*(1.0/(18.0e-3))*6.0221409e+23
    
    #-------------------------------------------------------------------------------------
    #2) Parse equation file 
    # Do files already exist? If so, you can bypass this stage and proceed with simulation
    # Note, this is for you to manage. If unsure which files are available, re-parse and re-compile
    # This is important since the species-2-dict array maps extracted species to array numbers. This
    # can change with each parse

    filename='MCM_test'    

    files_exist = False
    
    if files_exist is False:

        # Parse equation file and store relevant dictionaries for later retrieval
        print_options=dict()
        print_options['Full_eqn']=0 #Set to 1 to print details of all equations and rate coefficients parsed [usefu for checking]

        # Define the .eqn file to be used in the following
        outputdict=Parse_eqn_file.extract_mechanism('MCM_test.eqn.txt',print_options)

        # Collect the dictionaries generated
        reaction_dict=outputdict['reaction_dict']
        rate_dict=outputdict['rate_dict']
        rate_dict_reactants=outputdict['rate_dict_reactants']
        rate_def=outputdict['rate_def']
        loss_dict=outputdict['loss_dict']
        gain_dict=outputdict['gain_dict']
        stoich_dict=outputdict['stoich_dict']
        species_dict=outputdict['species_dict']
        species_dict2array=outputdict['species_dict2array']
        equations=outputdict['max_equations']
        num_species=len(species_dict.keys())
        #pdb.set_trace()
        print("Saving the mechanism dictionary as a pickled object for later retrieval")
        # save the dictionary to a file for later retrieval - have to do each seperately.
        with open(filename+'_species_dict2array.pickle', 'wb') as handle:
            pickle.dump(species_dict2array, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(filename+'_species_dict.pickle', 'wb') as handle:
            pickle.dump(species_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        #pdb.set_trace()
    
        # Convert the rate coefficient expressions into Fortran commands
        print("Converting rate coefficient operations into Python file")
        #rate_dict=rate_coeff_conversion.convert_rate_mcm(rate_dict)
        # Convert rate definitions in original *.eqn.txt file into a form to be used in Fortran
        # In the production model we use Numba for speed. You can select an equivalent 
        rate_dict=rate_coeff_conversion.convert_rate_mcm_numba(rate_dict)
        Parse_eqn_file.write_rate_file_numba(filename,rate_dict)    
        
        # Create python modules for product multiplications and dydt function
        # Also create sparse matrices for both operations if not using Numba
        print("Creating Python functions and sparse matrices for product multiplications and dydt function")
        Parse_eqn_file.write_reactants_indices(filename,equations,num_species,species_dict2array,rate_dict_reactants,loss_dict)
        Parse_eqn_file.write_loss_gain_matrix(filename,equations,num_species,loss_dict,gain_dict,species_dict2array)

        # Create .npy file with indices for all RO2 species
        print("Creating file that holds RO2 species indices")
        Parse_eqn_file.write_RO2_indices(filename,species_dict2array)
        
    else:
        
        # You have already parsed the .eqn file and stored relevant information. 
        # Load the dictioanties here to pass into the ODE solver
        with open(filename+'_species_dict2array.pickle', 'rb') as f:
            species_dict2array = pickle.load(f) 
        with open(filename+'_species_dict.pickle', 'rb') as f:
            species_dict = pickle.load(f) 
        
    print("Ready to run simulation") 
    pdb.set_trace()
            
    # Now load the numpy arrays generated in the parsing script for use in the simulations. Again, these are named to
    # match this particular code. However the ordering is not important.
    RO2_indices=numpy.load(filename+'_RO2_indices.npy')    
    #-------------------------------------------------------------------------------------
    # Define initial concentrations, in pbb, of species using names from KPP file
    species_initial_conc=dict()
    species_initial_conc['OH']=20.0
    species_initial_conc['TOLUENE']=200.0

    # Save this information to a dictionary to pass to ODE solver
    input_dict=dict()
    input_dict['species_dict']=species_dict
    input_dict['species_dict2array']=species_dict2array
    input_dict['species_initial_conc']=species_initial_conc
    input_dict['equations']=equations
    #-------------------------------------------------------------------------------------
    #3) Run the simulation
    run_simulation(filename, start_time, temp, RH, RO2_indices, H2O, input_dict)
    #-------------------------------------------------------------------------------------
    
