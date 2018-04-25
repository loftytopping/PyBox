##########################################################################################
#                                                                                        #
#    Contains definitions of functions use in RHS of ODE                                 #
#                                                                                        #
#    Mixed Python - Fortran version. This version uses the f2py module to re-write       #
#    the RHS calculations to exploit multi-core shared/distributed memory machine        #
#                                                                                        #
#                                                                                        #
#    Copyright (C) 2018  David Topping : david.topping@manchester.ac.uk                  #
#                                      : davetopp80@gmail.com                            #
#    Personal website: davetoppingsci.com                                                #
#                                                                                        #
#    This program does not yet have a license, meaning the deault copyright law applies. #
#    I will add an appropriate open-source icense once made public with paper            #
#    Only users who have access to the private repository that holds this file may       #
#    use it, but may not distribute it without explicit permission.                      #
#                                                                                        #
#                                                                                        #
##########################################################################################

import numpy 
import pylab as P
import pdb

def run_simulation(start_time, temp, RH, RO2_indices, H2O, input_dict):

    from assimulo.solvers import RodasODE, CVode #Choose solver accoring to your need. 
    from assimulo.problem import Explicit_Problem
    
    # In this function, we import functions that have been pre-compiled for use in the ODE solver
    # The function that calculates the RHS of the ODE is also defined within this function, such
    # that it can be used by the Assimulo solvers 
    
    # The variables passed to this function are defined as follows:
    
    #-------------------------------------------------------------------------------------
    # define the ODE function to be called
    def dydt_func(t,y):

        #pdb.set_trace()
        # Calculate time of day
        time_of_day_seconds=start_time+t
        
        # make sure the y array is not a list. Assimulo uses lists
        y_asnumpy=numpy.array(y)
        
        #Calculate the concentration of RO2 species, using an index file created during parsing
        RO2=numpy.sum(y[RO2_indices])

        #Calculate reaction rate for each equation.
        # Note that H2O will change in parcel mode
        # The time_of_day_seconds is used for photolysis rates - need to change this if want constant values
        rates=evaluate_rates_fortran(RO2,H2O,temp,time_of_day_seconds)
        #pdb.set_trace()
        # Calculate product of all reactants and stochiometry for each reaction [A^a*B^b etc]        
        reactants=reactants_fortran(y_asnumpy)
        #pdb.set_trace()
        #Multiply product of reactants with rate coefficient to get reaction rate            
        reactants = numpy.multiply(reactants,rates)
        #pdb.set_trace()
        # Now use reaction rates with the loss_gain matri to calculate the final dydt for each compound
        # With the assimulo solvers we need to output numpy arrays
        dydt=loss_gain_fortran(reactants)
        #pdb.set_trace()
        
        return dydt
    #-------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------
    # define jacobian function to be called
    def jacobian(t,y):    
        # Different solvers might call jacobian at different stages, so we have to redo some calculations here 
        # Calculate time of day
        time_of_day_seconds=start_time+t
        
        # make sure the y array is not a list. Assimulo uses lists
        y_asnumpy=numpy.array(y)
        
        #Calculate the concentration of RO2 species, using an index file created during parsing
        RO2=numpy.sum(y[RO2_indices])

        #Calculate reaction rate for each equation.
        # Note that H2O will change in parcel mode
        rates=evaluate_rates_fortran(RO2,H2O,temp,time_of_day_seconds)
        #pdb.set_trace()
        # Now use reaction rates with the loss_gain matrix to calculate the final dydt for each compound
        # With the assimulo solvers we need to output numpy arrays
        dydt_dydt=jacobian_fortran(rates,y_asnumpy)
        #pdb.set_trace()
        return dydt_dydt
    #-------------------------------------------------------------------------------------

    #import static compilation of Fortran functions for use in ODE solver
    print("Importing pre-compiled Fortran modules")
    from rate_coeff_f2py import evaluate_rates as evaluate_rates_fortran
    from reactants_conc_f2py import reactants as reactants_fortran
    from loss_gain_f2py import loss_gain as loss_gain_fortran  
    from jacobian_f2py import jacobian as jacobian_fortran
    
    # 'Unpack' variables from input_dict
    species_dict=input_dict['species_dict']
    species_dict2array=input_dict['species_dict2array']
    species_initial_conc=input_dict['species_initial_conc']
    equations=input_dict['equations']
    
    #Specify some starting concentrations [ppt]
    Cfactor= 2.55e+10 #ppb-to-molecules/cc
    
    # Create variables required to initialise ODE
    num_species=len(species_dict.keys())
    y0 = [0]*num_species #Initial concentrations, set to 0
    t0 = 0.0 #T0
        
    # Define species concentrations in ppb
    # You have already set this in the front end script, and now we populate the y array with those concentrations
    for specie in species_initial_conc.keys():
        y0[species_dict2array[specie]]=species_initial_conc[specie]*Cfactor #convert from pbb to molcules/cc
        
    #Set the total_time of the simulation to 0 [havent done anything yet]
    total_time=0.0
    
    # Now run through the simulation in batches. 
    # I do this to enable testing of coupling processes. Some initial investigations with non-ideality in
    # the condensed phase indicated that even defining a maximum step was not enough for ODE solvers to 
    # overshoot a stable region. It also helps with in-simulation debugging. Its up to you if you want to keep this.
    # To not run in batches, just define one batch as your total simulation time. This will reduce any overhead with
    # initialising the solvers
    # Set total simulation time and batch steps in seconds
    
    # Note also that the current module outputs solver information after each batch step. This can be turned off and the
    # the batch step change for increased speed
    simulation_time= 3600.0
    batch_step=100.0
    t_array=[]
    time_step=0
    number_steps=int(simulation_time/batch_step) # Just cycling through 3 steps to get to a solution
    
    # Define a matrix that stores values as outputs from the end of each batch step. Again, you can remove
    # the need to run in batches. You can tell the Assimulo solvers the frequency of outputs.
    y_matrix=numpy.zeros((int(number_steps),len(y0)))
    
    print("Starting simulation")

    # In the following, we can 
    while total_time < simulation_time:
        
        if total_time == 0.0:
            #Define an Assimulo problem
            #Define an explicit solver
            exp_mod = Explicit_Problem(dydt_func,y0,t0, name = 'MCM simulation')
            
        else:
            y0 = y_output[-1,:] # Take the output from the last batch as the start of this
            exp_mod = Explicit_Problem(dydt_func,y0,t0, name = 'MCM simulation')
            
        # Define ODE parameters. 
        # Initial steps might be slower than mid-simulation. It varies.
        #exp_mod.jac = dydt_jac
        # Define which ODE solver you want to use
        exp_sim = CVode(exp_mod) 
        tol_list=[1.0e-3]*num_species
        exp_sim.atol = tol_list #Default 1e-6
        exp_sim.rtol = 0.03 #Default 1e-6
        exp_sim.inith = 1.0e-6 #Initial step-size
        #exp_sim.discr = 'Adams'
        exp_sim.maxh = 100.0
        # Use of a jacobian makes a big differece in simulation time. This is relatively 
        # easy to define for a gas phase - not sure for an aerosol phase with composition
        # dependent processes. 
        exp_sim.usejac = True # To be provided as an option in future update. 
        #exp_sim.fac1 = 0.05
        #exp_sim.fac2 = 50.0
        exp_sim.report_continuously = True
        exp_sim.maxncf = 1000
        #Sets the parameters        
        t_output, y_output = exp_sim.simulate(batch_step) #Simulate 'batch' seconds
        total_time+=batch_step
        t_array.append(total_time) # Save the output from the end step, of the current batch, to a matrix
        y_matrix[time_step,:]=y_output[-1,:]
                
        #now save this information into a matrix for later plotting.
        time_step+=1
        
    with_plots=True
    
    #pdb.set_trace()
    #Plot the change in concentration over time for a given specie. For the user to change / remove
    #In a future release I will add this as a seperate module
    if with_plots:
        P.plot(t_array,y_matrix[:,species_dict2array['APINENE']], marker='o')
        P.title(exp_mod.name)
        P.ylabel("State: $y_1$")
        P.xlabel("Time")
        P.show()
        
    