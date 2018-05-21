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
#    Copyright (C) 2018  David Topping : david.topping@manchester.ac.uk                  #
#                                      : davetopp80@gmail.com                            #
#    Personal website: davetoppingsci.com                                                #
#                                                                                        #
#    All Rights Reserved.                                                                #
#    This file is part of PyBox.                                                         #
#                                                                                        #
#    PyBox is free software: you can redistribute it and/or modify it under              #
#    the terms of the GNU General Public License as published by the Free Software       #
#    Foundation, either version 3 of the License, or (at your option) any later          #
#    version.                                                                            #
#                                                                                        #
#    PyBox is distributed in the hope that it will be useful, but WITHOUT                #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              #
#    details.                                                                            #
#                                                                                        #
#    You should have received a copy of the GNU General Public License along with        #
#    PyBox.  If not, see <http://www.gnu.org/licenses/>.                                 #
#                                                                                        #
##########################################################################################
# Developed using the Anaconda Python 3 distribution and with the Assimulo ODE solver    # 
# suite: http://www.jmodelica.org/assimulo                                               #
# In the import statements, all files developed specifically for this project            #
# as marked [•]                                                                          #
##########################################################################################

# VERSION 1.0

# Developed using the Anaconda Python distribution and with the Assimulo ODE solver suite: http://www.jmodelica.org/assimulo
# Assimulo gives a Python front-end acces to stiff solvers so this can be merged with aerosol developments in a seperate branch.
# The Assimulo package does not allow extra argument passing, and thus defines the structure of the code. 
# In the import statements, all files developed specifically for this project as marked [•]

import numpy 
import sys
import os
sys.path.append(os.path.abspath('../..'))
from PyBox import Parse_eqn_file # [•] Needed to parse the .eqn file, name given in this file
from PyBox import rate_coeff_conversion # [•] Converts standard text rate coefficients into Python/Fortran
from PyBox import MCM_constants # [•] holds more pre-defined rate coefficients and photolysis rates not provided in .eqn file. 
import collections
import pdb
from datetime import datetime
import time
from ODE_solver import run_simulation # [•] Contains routines to run ODE solver
import pickle
from shutil import copy2
   
# Start of the main body of code
if __name__=='__main__':
   
    #-------------------------------------------------------------------------------------

    #1)Define starting ambient conditions
    temp=288.15
    RH=0.5 #RH/100%
    #Define a start time 
    hour_of_day=12.0
    start_time=hour_of_day*60*60 # seconds, used as t0 in solver
    simulation_time= 7200.0 # seconds
    batch_step=100.0 # seconds
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

    filename='MCM_APINENE'    

    files_exist = False
    
    if files_exist is False:

        # Delete any 'straggling' f90 or cython files
        for fname in os.listdir('.'):
            #if "f2py.cpython" in fname:
            #    os.remove(fname)
            if ".f90" in fname:
                os.remove(fname)
            if ".eqn.txt" in fname:
                os.remove(fname)

        # Copy mechanism file into working directory
        copy2('../mechanism_files/'+filename+'.eqn.txt','.')

        # Parse equation file and store relevant dictionaries for later retrieval
        print_options=dict()
        print_options['Full_eqn']=0 #Set to 1 to print details of all equations and rate coefficients parsed [usefu for checking]

        # Define the .eqn file to be used in the following
        outputdict=Parse_eqn_file.extract_mechanism(filename+'.eqn.txt',print_options)

        # Collect the dictionaries generated
        reaction_dict=outputdict['reaction_dict']
        rate_dict=outputdict['rate_dict']
        rate_dict_fortran=rate_coeff_conversion.convert_rate_mcm_fortran(rate_dict)
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
        with open(filename+'_equations.pickle', 'wb') as handle:
            pickle.dump(equations, handle, protocol=pickle.HIGHEST_PROTOCOL)
            
        # Generate some static Fortran libraries for use in ODE solver
        # Here you have the option to use OpenMP for spreading caclualtions across available cores
        # This will automatically map to the number of cores on your machine
        openMP=True
        
        #pdb.set_trace()
    
        # Convert the rate coefficient expressions into Fortran commands
        print("Converting rate coefficient operation into Fortran file")
        #rate_dict=rate_coeff_conversion.convert_rate_mcm(rate_dict)
        # Convert rate definitions in original *.eqn.txt file into a form to be used in Fortran
        rate_dict_fortran=rate_coeff_conversion.convert_rate_mcm_fortran(rate_dict)
        Parse_eqn_file.write_rate_file_fortran(filename,rate_dict_fortran,openMP)    
        print("Compiling rate coefficient file using f2py")
        #Parse_eqn_file.write_rate_file(filename,rate_dict,mcm_constants_dict)
        os.system("python f2py_rate_coefficient.py build_ext --inplace --fcompiler=gfortran")
        
        # Create Fortran file for calculating prodcts all of reactants for all reactions
        print("Creating Fortran file to calculate reactant contribution to equation")
        Parse_eqn_file.write_reactants_indices_fortran(filename,equations,species_dict2array,rate_dict_reactants,loss_dict,openMP)
        print("Compiling reactant product file using f2py")
        os.system("python f2py_reactant_conc.py build_ext --inplace --fcompiler=gfortran")
        
        # Create Fortran file for calculating dy_dt
        print("Creating Fortran file to calculate dy_dt for each reaction")
        Parse_eqn_file.write_loss_gain_fortran(filename,equations,num_species,loss_dict,gain_dict,species_dict2array,openMP)
        print("Compiling dydt file using f2py")
        os.system("python f2py_loss_gain.py build_ext --inplace --fcompiler=gfortran")

        # Create .npy file with indices for all RO2 species
        print("Creating file that holds RO2 species indices")
        Parse_eqn_file.write_RO2_indices(filename,species_dict2array)
        
        # Create jacobian 
        Parse_eqn_file.write_gas_jacobian_fortran(filename,equations,num_species,loss_dict,gain_dict,species_dict2array,rate_dict_reactants,openMP)
        print("Compiling jacobian function using f2py")      
        os.system("python f2py_jacobian.py build_ext --inplace --fcompiler=gfortran")

    else:
        
        # You have already parsed the .eqn file and stored relevant information. 
        # Load the dictioanties here to pass into the ODE solver
        with open(filename+'_species_dict2array.pickle', 'rb') as f:
            species_dict2array = pickle.load(f) 
        with open(filename+'_species_dict.pickle', 'rb') as f:
            species_dict = pickle.load(f) 
        with open(filename+'_equations.pickle', 'rb') as f:
            equations = pickle.load(f) 
        
    print("Ready to run simulation") 
    #pdb.set_trace()
            
    # Now load the numpy arrays generated in the parsing script for use in the simulations. Again, these are named to
    # match this particular code. However the ordering is not important.
    RO2_indices=numpy.load(filename+'_RO2_indices.npy')    
    #-------------------------------------------------------------------------------------
    # Define initial concentrations, in pbb, of species using names from KPP file
    species_initial_conc=dict()
    species_initial_conc['O3']=18.0
    species_initial_conc['APINENE']=30.0
    
    # Save this information to a dictionary to pass to ODE solver
    input_dict=dict()
    input_dict['species_dict']=species_dict
    input_dict['species_dict2array']=species_dict2array
    input_dict['species_initial_conc']=species_initial_conc
    input_dict['equations']=equations

    #Do you want to save the output from the simulation as a .npy file?
    save_output=True
    #-------------------------------------------------------------------------------------
    #3) Run the simulation
    run_simulation(filename, start_time, save_output, temp, RH, RO2_indices, H2O, input_dict, simulation_time, batch_step)
    #-------------------------------------------------------------------------------------
    

