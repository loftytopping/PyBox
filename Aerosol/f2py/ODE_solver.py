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

def run_simulation(filename, save_output, start_time, temp, RH, RO2_indices, H2O, PInit, y_cond, input_dict):

    from assimulo.solvers import RodasODE, CVode #Choose solver accoring to your need. 
    from assimulo.problem import Explicit_Problem
    
    # In this function, we import functions that have been pre-compiled for use in the ODE solver
    # The function that calculates the RHS of the ODE is also defined within this function, such
    # that it can be used by the Assimulo solvers 
    
    # The variables passed to this function are defined as follows:
    
    #-------------------------------------------------------------------------------------    
    #-------------------------------------------------------------------------------------
    # define the ODE function to be called
    def dydt_func(t,y):

        dy_dt=numpy.zeros((total_length_y,1),)
        
        #pdb.set_trace()
        # Calculate time of day
        time_of_day_seconds=start_time+t
        
        #pdb.set_trace()
        # make sure the y array is not a list. Assimulo uses lists
        y_asnumpy=numpy.array(y)
        Model_temp = temp
        #pdb.set_trace()
        #Calculate the concentration of RO2 species, using an index file created during parsing
        RO2=numpy.sum(y[RO2_indices])

        #Calculate reaction rate for each equation.
        # Note that H2O will change in parcel mode
        # The time_of_day_seconds is used for photolysis rates - need to change this if want constant values
        rates=evaluate_rates_fortran(RO2,H2O,Model_temp,time_of_day_seconds)
        #pdb.set_trace()
        # Calculate product of all reactants and stochiometry for each reaction [A^a*B^b etc]        
        reactants=reactants_fortran(y_asnumpy[0:num_species-1])
        #pdb.set_trace()
        #Multiply product of reactants with rate coefficient to get reaction rate            
        reactants = numpy.multiply(reactants,rates)
        #pdb.set_trace()
        # Now use reaction rates with the loss_gain matri to calculate the final dydt for each compound
        # With the assimulo solvers we need to output numpy arrays
        dydt_gas=loss_gain_fortran(reactants)
        pdb.set_trace()
        
        dy_dt[0:num_species-1,0]=dydt_gas
        
        # Change the saturation vapour pressure of water
        # Need to re-think the change of organic vapour pressures with temperature.
        # At the moment this is kept constant as re-calulation using UManSysProp very slow
        sat_vap_water = numpy.exp((-0.58002206E4 / Model_temp) + 0.13914993E1 - \
        (0.48640239E-1 * Model_temp) + (0.41764768E-4 * (Model_temp**2.0E0))- \
        (0.14452093E-7 * (Model_temp**3.0E0)) + (0.65459673E1 * numpy.log(Model_temp)))
        sat_vp[-1]=(numpy.log10(sat_vap_water*9.86923E-6))
        Psat=numpy.power(10.0,sat_vp)    
        
        # Convert the concentration of each component in the gas phase into a partial pressure using the ideal gas law
        # Units are Pascals
        Pressure_gas=(y_asnumpy[0:num_species,]/NA)*8.314E+6*Model_temp #[using R]
    
        core_mass_array=numpy.multiply(ycore_asnumpy/NA,core_molw_asnumpy)
            
        ####### Calculate the thermal conductivity of gases according to the new temperature ########
        K_water_vapour = (5.69+0.017*(Model_temp-273.15))*1e-3*4.187 #[W/mK []has to be in W/m.K]
        # Use this value for all organics, for now. If you start using a non-zero enthalpy of
        # vapourisation, this needs to change. 
        therm_cond_air = K_water_vapour

        #----------------------------------------------------------------------------
        #F2c) Extract the current gas phase concentrations to be used in pressure difference calculations
        C_g_i_t = y_asnumpy[0:num_species,]
        #Set the values for oxidants etc to 0 as will force no mass transfer
        C_g_i_t[ignore_index]=0.0
    
        pdb.set_trace()
        
        test=dydt_partition_fortran(total_length_y,num_bins,num_species,y_asnumpy,ycore_asnumpy,core_dissociation,core_mass_array,y_density_array_asnumpy,core_density_array_asnumpy,ignore_index_fortran,y_mw,Psat,DStar_org_asnumpy,alpha_d_org_asnumpy,C_g_i_t,N_perbin,gamma_gas_asnumpy,Latent_heat_asnumpy,GRAV,Updraft,sigma,NA,kb,Rv,R_gas,Model_temp,cp,Ra,Lv_water_vapour)
        
        total_SOA_mass,aw_array,size_array,dy_dt_calc = dydt_partition_fortran(len(y),num_bins,num_species,\
        y_asnumpy,ycore_asnumpy,core_dissociation, \
        core_mass_array,y_density_array_asnumpy,core_density_array_asnumpy,ignore_index_fortran,y_mw,Psat, \
        DStar_org_asnumpy,alpha_d_org_asnumpy,C_g_i_t,N_perbin,gamma_gas_asnumpy,Latent_heat_asnumpy,GRAV, \
        Updraft,sigma,NA,kb,Rv,R_gas,Model_temp,cp,Ra,Lv_water_vapour)
        
        pdb.set_trace()
        
        # Add the calculated gains/losses to the complete dy_dt array
        dy_dt[0:num_species+(num_species*num_bins),0]+=dy_dt_calc[:]
    
        #----------------------------------------------------------------------------
        #F4) Now calculate the change in water vapour mixing ratio. 
        #To do this we need to know what the index key for the very last element is     
        #pdb.set_trace()  
        pdb.set_trace()        
        #print "elapsed time=", elapsedTime   
        dydt_func.total_SOA_mass=total_SOA_mass
        dydt_func.size_array=size_array
        dydt_func.temp=Model_temp
        dydt_func.RH=Pressure_gas[-1]/(Psat[-1]*101325.0)
        dydt_func.water_activity=aw_array
                
        #----------------------------------------------------------------------------        
        return dy_dt
    #-------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------

    #import static compilation of Fortran functions for use in ODE solver
    print("Importing pre-compiled Fortran modules")
    from rate_coeff_f2py import evaluate_rates as evaluate_rates_fortran
    from reactants_conc_f2py import reactants as reactants_fortran
    from loss_gain_f2py import loss_gain as loss_gain_fortran  
    from partition_f2py import dydt_partition as dydt_partition_fortran
    
    # 'Unpack' variables from input_dict
    species_dict=input_dict['species_dict']
    species_dict2array=input_dict['species_dict2array']
    species_initial_conc=input_dict['species_initial_conc']
    equations=input_dict['equations']
    num_species=input_dict['num_species']
    y_density_array_asnumpy=input_dict['y_density_array_asnumpy']
    y_mw=input_dict['y_mw']
    sat_vp=input_dict['sat_vp']
    Delta_H=input_dict['Delta_H']
    Latent_heat_asnumpy=input_dict['Latent_heat_asnumpy']
    DStar_org_asnumpy=input_dict['DStar_org_asnumpy']
    alpha_d_org_asnumpy=input_dict['alpha_d_org_asnumpy']
    gamma_gas_asnumpy=input_dict['gamma_gas_asnumpy']
    Updraft=input_dict['Updraft']
    GRAV=input_dict['GRAV']
    Rv=input_dict['Rv']
    Ra=input_dict['Ra']
    R_gas=input_dict['R_gas']
    R_gas_other=input_dict['R_gas_other']
    cp=input_dict['cp']
    sigma=input_dict['sigma']
    NA=input_dict['NA']
    kb=input_dict['kb']
    Lv_water_vapour=input_dict['Lv_water_vapour']
    ignore_index=input_dict['ignore_index']
    ignore_index_fortran=input_dict['ignore_index_fortran']
    ycore_asnumpy=input_dict['ycore_asnumpy']
    core_density_array_asnumpy=input_dict['core_density_array_asnumpy']
    y_cond=input_dict['y_cond_initial']
    num_bins=input_dict['num_bins']
    core_molw_asnumpy=input_dict['core_molw_asnumpy']
    core_dissociation=input_dict['core_dissociation']
    N_perbin=input_dict['N_perbin']
    
    pdb.set_trace()
    
    #Specify some starting concentrations [ppt]
    Cfactor= 2.55e+10 #ppb-to-molecules/cc
    
    # Create variables required to initialise ODE
    y0 = [0]*(num_species+num_species*num_bins) #Initial concentrations, set to 0
    t0 = 0.0 #T0
        
    # Define species concentrations in ppb fr the gas phase
    # You have already set this in the front end script, and now we populate the y array with those concentrations
    for specie in species_initial_conc.keys():
        if specie is not 'H2O':
            y0[species_dict2array[specie]]=species_initial_conc[specie]*Cfactor #convert from pbb to molcules/cc
        elif specie is 'H2O':
            y0[species_dict2array[specie]]=species_initial_conc[specie]

    # Now add the initial condensed phase [including water]
    #pdb.set_trace()
    y0[num_species:num_species+((num_bins)*num_species)]=y_cond[:]
    #pdb.set_trace()
    
    #Set the total_time of the simulation to 0 [havent done anything yet]
    total_time=0.0
    
    # Define a 'key' that represents the end of the composition variables to track
    total_length_y=len(y0)
    key=num_species+((num_bins)*num_species)-1
    
    #pdb.set_trace()
    
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
        tol_list=[1.0e-3]*len(y0)
        exp_sim.atol = tol_list #Default 1e-6
        exp_sim.rtol = 1.0e-6 #Default 1e-6
        exp_sim.inith = 1.0e-6 #Initial step-size
        #exp_sim.discr = 'Adams'
        exp_sim.maxh = 100.0
        # Use of a jacobian makes a big differece in simulation time. This is relatively 
        # easy to define for a gas phase - not sure for an aerosol phase with composition
        # dependent processes. 
        exp_sim.usejac = False # To be provided as an option in future update. 
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
        
    