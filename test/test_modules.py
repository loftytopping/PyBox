##########################################################################################
#                                                                                        #
#    Perform unit tests using saved data held within data folder                         #
#    This tests the expected output from modules created by Parse_eqn_file as used in    #
#    the ODE solvers. This is not to test the output from the ODE solvers in Assimulo    #
#    , or other [Scipy], since these are stand alone modules with rigirous dev. Also,    #
#    each solver has multiple options and would render comparisons complicated           #
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
import sys
import os
sys.path.append(os.path.abspath('../'))
sys.path.append(os.path.abspath('../Aerosol/'))
import Parse_eqn_file # [•] Needed to parse the .eqn file, name given in this file
import rate_coeff_conversion # [•] Converts standard text rate coefficients into Python/Fortran
import Size_distributions # [•] Create size distribution according to number of bins and core-material
import collections
import pdb
from datetime import datetime
import time
import Property_calculation
import pickle
import shutil
import numpy.testing as npt
import unittest

def setup(filename):

    """
    Function that is used to setup the parsed routines to check expected output via unitttests
    """

    # Create modules through parsing routine and compile where needed

    # Define standard variables used in testing
    temp=288.15
    RH=0.5 #RH/100%
    PInit=98000 #Pascals - Starting pressure of parcel expansion [if run in Parcel model mode]
    #Define a start time 
    hour_of_day=12.0
    start_time=hour_of_day*60*60 # seconds, used as t0 in solver
    time_of_day_seconds=start_time
    #Convert RH to concentration of water vapour molecules [this will change when in Parcel model mode]
    temp_celsius=temp-273.15
    # Saturation VP of water vapour, to get concentration of H20
    Psat=610.78*numpy.exp((temp_celsius/(temp_celsius+238.3))*17.2694)
    Pw=RH*Psat
    Updraft=0.0
    Wconc=0.002166*(Pw/(temp_celsius+273.16)) #kg/m3
    #Convert from m3 to cm3
    Wconc=Wconc*1.0e-6
    #Convert from kg to molecules/cc
    H2O=Wconc*(1.0/(18.0e-3))*6.0221409e+23
    Model_temp=temp
    
    Lv_water_vapour=2.5e3 # Latent heat of vapourisation of water [J/g] 
    Rv=461.0 #Individual gas constant of water vapour [J/Kg.K]
    Ra=287.0 #Gas constant for dry air [J/Kg.K]
    R_gas=8.3144598 #Ideal gas constant [kg m2 s-2 K-1 mol-1]
    R_gas_other=8.2057e-5 #Ideal gas constant [m3 bar K-1 mol-1]
    GRAV=9.8; #Gravitational acceleration [(m/2)2]
    cp=1005; #Specific heat capacity of air [J/Kg.K]
    sigma=72.0e-3 # Assume surface tension of water (mN/m)
    NA=6.0221409e+23 #Avogadros number
    kb=1.380648E-23 #Boltzmanns constant

    # check if __pycache__ exists
    if os.path.isdir("__pycache__"):
        my_dir = '__pycache__' # enter the dir name
        for fname in os.listdir(my_dir):
            if "_numba_" in fname:
                os.remove(os.path.join(my_dir, fname))
    # Delete any 'straggling' f90 or cython files
    for fname in os.listdir('.'):
        #if "f2py.cpython" in fname:
        #    os.remove(fname)
        if ".f90" in fname:
            os.remove(fname)
    
    for fname in os.listdir('.'):
        if ".npy" in fname:
            os.remove(fname)
        if ".pickle" in fname:
            os.remove(fname)
        if ".npz" in fname:
            os.remove(fname)
    
    # Define the .eqn file to be used in the following
    print_options=dict()
    print_options['Full_eqn']=0
    outputdict=Parse_eqn_file.extract_mechanism(filename+'.eqn.txt',print_options)
        
    # Now map these species onto SMILES according to the relevant .xml file that comes with the MCM. If this file changes
    # you will need to change the reference here
    outputdict=Parse_eqn_file.extract_smiles_species(outputdict,'../Aerosol/MCM.xml')

    #pdb.set_trace()
    # Collect the dictionaries generated
    #reaction_dict=outputdict['reaction_dict']
    rate_dict=outputdict['rate_dict']
    rate_dict_fortran=rate_coeff_conversion.convert_rate_mcm_fortran(rate_dict)
    rate_dict_reactants=outputdict['rate_dict_reactants']
    #rate_def=outputdict['rate_def']
    loss_dict=outputdict['loss_dict']
    gain_dict=outputdict['gain_dict']
    stoich_dict=outputdict['stoich_dict']
    species_dict=outputdict['species_dict']
    species_dict2array=outputdict['species_dict2array']
    equations=outputdict['max_equations']
    Pybel_object_dict=outputdict['Pybel_object_dict'] # Might be different size to mechanism extracted dicts [due to ignoring oxidants etc]
    SMILES_dict=outputdict['SMILES_dict']
    num_species=len(species_dict.keys())

    vp_method='nannoolal' # Saturation vapour pressure ['nannoolal': 'myrdal_and_yalkowsky': 'evaporation']
    bp_method='joback_and_reid' # Boiling point ['joback_and_reid': 'stein_and_brown': 'nannoolal']
    critical_method='nannoolal' # Critical properties for density ['nannoolal':'joback_and_reid']
    density_method='girolami' # Pure component liquid density ['girolami': 'schroeder':'le_bas':'tyn_and_calus']
    ignore_vp=False
    vp_cutoff=-6.0
        
    # Now calculate all properties that dictate gas-to-particle partitioning
    print("Calculating properties that dictate gas-to-particle partitioning")
    property_dict1=Property_calculation.Pure_component1(num_species,species_dict,
            species_dict2array,Pybel_object_dict,SMILES_dict,temp,vp_method,bp_method,
            critical_method,density_method,ignore_vp,vp_cutoff)


    y_density_array=property_dict1['y_density_array']
    y_density_array_asnumpy=numpy.array(y_density_array)
    y_mw=property_dict1['y_mw']
    sat_vp=property_dict1['sat_vp']
    Psat=numpy.power(10.0,sat_vp) 
    Delta_H=property_dict1['Delta_H']
    Latent_heat_gas=property_dict1['Latent_heat_gas']  
    Latent_heat_asnumpy=numpy.array(Latent_heat_gas)	
    ignore_index=property_dict1['ignore_index'] 
    ignore_index_fortran=property_dict1['ignore_index_fortran'] 
    include_index=sorted(property_dict1['include_index']) # This is an array that captures compounds to be included in 
    num_species_condensed=len(y_density_array)

    property_dict2=Property_calculation.Pure_component2(num_species_condensed,y_mw,R_gas,temp)
    alpha_d_org=property_dict2['alpha_d_org']
    alpha_d_org_asnumpy=numpy.array(alpha_d_org)
    DStar_org=property_dict2['DStar_org']
    DStar_org_asnumpy=numpy.array(DStar_org)
    mean_them_vel=property_dict2['mean_them_vel']
    gamma_gas=property_dict2['gamma_gas']
    gamma_gas_asnumpy=numpy.array(gamma_gas)
    

    num_bins=16
    y_cond=[0.0]*num_species_condensed*num_bins #array that contains all species in all bins
                                      #water is the final component
    total_conc=100 #Total particles per cc
    std=2.2 #Standard Deviation
    lowersize=0.01 #microns
    uppersize=1.0 #microns
    meansize=0.2 #microns
    #Create a number concentration for a lognormal distribution
    N_perbin,x=Size_distributions.lognormal(num_bins,total_conc,meansize,std,lowersize,uppersize)

    openMP=False

    # Convert the rate coefficient expressions into Fortran commands
    print("Converting rate coefficient operation into Fortran file")
    #rate_dict=rate_coeff_conversion.convert_rate_mcm(rate_dict)
    # Convert rate definitions in original *.eqn.txt file into a form to be used in Fortran
    rate_dict_fortran=rate_coeff_conversion.convert_rate_mcm_fortran(rate_dict)
    Parse_eqn_file.write_rate_file_fortran(filename,rate_dict_fortran,openMP)    
    print("Compiling rate coefficient file using f2py")
    #Parse_eqn_file.write_rate_file(filename,rate_dict,mcm_constants_dict)
    #os.system("python f2py_rate_coefficient.py build_ext --inplace")
    os.system('f2py -c -m rate_coeff_f2py Rate_coefficients.f90 --f90flags="-O3 -ffast-math -fopenmp" -lgomp')
        
    # Create Fortran file for calculating prodcts all of reactants for all reactions
    print("Creating Fortran file to calculate reactant contribution to equation")
    Parse_eqn_file.write_reactants_indices_fortran(filename,equations,species_dict2array,rate_dict_reactants,loss_dict,openMP)
    print("Compiling reactant product file using f2py")
    #os.system("python f2py_reactant_conc.py build_ext --inplace")
    os.system('f2py -c -m reactants_conc_f2py Reactants_conc.f90 --f90flags="-O3 -ffast-math -fopenmp" -lgomp')
        
    # Create Fortran file for calculating dy_dt
    print("Creating Fortran file to calculate dy_dt for each reaction")
    Parse_eqn_file.write_loss_gain_fortran(filename,equations,num_species,loss_dict,gain_dict,species_dict2array,openMP)
    print("Compiling dydt file using f2py")
    #os.system("python f2py_loss_gain.py build_ext --inplace")
    os.system('f2py -c -m loss_gain_f2py Loss_Gain.f90 --f90flags="-O3 -ffast-math -fopenmp" -lgomp')

    Parse_eqn_file.write_gas_jacobian_fortran(filename,equations,num_species,loss_dict,gain_dict,species_dict2array,rate_dict_reactants,openMP)
    print("Compiling jacobian function using f2py")      
    #os.system("python f2py_jacobian.py build_ext --inplace")
    os.system('f2py -c -m jacobian_f2py Jacobian.f90 --f90flags="-O3 -ffast-math -fopenmp" -lgomp')

    # Create .npy file with indices for all RO2 species
    print("Creating file that holds RO2 species indices")
    Parse_eqn_file.write_RO2_indices(filename,species_dict2array)
    RO2_indices=numpy.load(filename+'_RO2_indices.npy')  

    # Convert the rate coefficient expressions into Numba commands
    print("Converting rate coefficient operations into Python-Numba file")
    #rate_dict=rate_coeff_conversion.convert_rate_mcm(rate_dict)
    # Convert rate definitions in original *.eqn.txt file into a form to be used in Fortran
    # In the production model we use Numba for speed. You can select an equivalent 
    rate_dict=rate_coeff_conversion.convert_rate_mcm_numba(rate_dict)
    Parse_eqn_file.write_rate_file_numba(filename,rate_dict)    
    
    # Create python modules for product multiplications and dydt function
    # Also create sparse matrices for both operations if not using Numba
    print("Creating Python-Numba functions and sparse matrices for product multiplications and dydt function")
    Parse_eqn_file.write_reactants_indices(filename,equations,num_species,species_dict2array,rate_dict_reactants,loss_dict)
    Parse_eqn_file.write_loss_gain_matrix(filename,equations,num_species,loss_dict,gain_dict,species_dict2array)

    Parse_eqn_file.write_gas_jacobian_numba(filename,equations,num_species,loss_dict,gain_dict,species_dict2array,rate_dict_reactants)
    print("Creating jacobian function in Numba")      

    # Create a Fortran file for calculating gas-to-particle partitioning drivers
    print("Creating Fortran file to calculate gas-to-particle partitining for each compound")
    Parse_eqn_file.write_partitioning_section_fortran_ignore(num_species+num_species_condensed*num_bins,num_bins,num_species,num_species_condensed,include_index)
    print("Compiling gas-to-particle partitioning file using f2py")
    #os.system("python f2py_partition.py build_ext --inplace")        
    os.system('f2py -c -m partition_f2py Partitioning.f90 --f90flags="-O3 -ffast-math -fopenmp" -lgomp')


    # Load the compiled functions
    print("Importing Fortran compiled modules")
    from rate_coeff_f2py import evaluate_rates as evaluate_rates_fortran
    from reactants_conc_f2py import reactants as reactants_fortran
    from loss_gain_f2py import loss_gain as loss_gain_fortran  
    from jacobian_f2py import jacobian as jacobian_fortran
    print("Importing Numba modules [compiling]")
    from Rate_coefficients_numba import evaluate_rates as evaluate_rates_numba
    from Reactants_conc_numba import reactants as reactant_numba
    from Loss_Gain_numba import dydt as loss_gain_numba
    from Jacobian_numba import jacobian as jacobian_numba
    from partition_f2py import dydt_partition as dydt_partition_fortran    

    # Run the function tests and save data to files. These files wll be compared with the test output
    # in individual functions.
    y_asnumpy=numpy.array([1.0e12]*num_species)
    RO2=numpy.sum(y_asnumpy[RO2_indices])
    
    #numpy.save(filename+'_RO2_base', RO2)
    #shutil.move(filename+'_RO2_base.npy','./data')
    numpy.save(filename+'_RO2', RO2)

    # Fortran modules
    rates_fortran=evaluate_rates_fortran(RO2,H2O,temp,time_of_day_seconds)
    #numpy.save(filename+'rates_fortran_base', rates_fortran)
    #shutil.move(filename+'rates_fortran_base.npy','./data')
    numpy.save(filename+'rates_fortran', rates_fortran)
    
    # Calculate product of all reactants and stochiometry for each reaction [A^a*B^b etc]        
    reactants_fortran=reactants_fortran(y_asnumpy)
    #Multiply product of reactants with rate coefficient to get reaction rate            
    reactants_fortran = numpy.multiply(reactants_fortran,rates_fortran)
    #numpy.save(filename+'reactants_fortran_base', reactants_fortran)
    #shutil.move(filename+'reactants_fortran_base.npy','./data')
    numpy.save(filename+'reactants_fortran', reactants_fortran)
    # Now use reaction rates with the loss_gain matri to calculate the final dydt for each compound
    # With the assimulo solvers we need to output numpy arrays
    dydt_fortran=loss_gain_fortran(reactants_fortran)
    #numpy.save(filename+'dydt_fortran_base', dydt_fortran)
    #shutil.move(filename+'dydt_fortran_base.npy','./data')
    numpy.save(filename+'dydt_fortran', dydt_fortran)
    
    dydt_dydt_fortran=jacobian_fortran(rates_fortran,y_asnumpy)
    #numpy.save(filename+'dydt_dydt_fortran_base', dydt_dydt_fortran)
    #shutil.move(filename+'dydt_dydt_fortran_base.npy','./data')
    numpy.save(filename+'dydt_dydt_fortran', dydt_dydt_fortran)
    
    # Numba modules
    rates_numba=evaluate_rates_numba(time_of_day_seconds,RO2,H2O,temp,numpy.zeros((equations)),numpy.zeros((63)))
    #numpy.save(filename+'rates_numba_base', rates_numba)
    #shutil.move(filename+'rates_numba_base.npy','./data')
    numpy.save(filename+'rates_numba', rates_numba)
    
    # Calculate product of all reactants and stochiometry for each reaction [A^a*B^b etc]        
    reactants_numba=reactant_numba(y_asnumpy,equations,numpy.zeros((equations)))
    #Multiply product of reactants with rate coefficient to get reaction rate            
    reactants_numba = numpy.multiply(reactants_numba,rates_numba)
    #numpy.save(filename+'reactants_numba_base', reactants_numba)
    #shutil.move(filename+'reactants_numba_base.npy','./data')
    numpy.save(filename+'reactants_numba', reactants_numba)
    
    # Now use reaction rates with the loss_gain matri to calculate the final dydt for each compound
    # With the assimulo solvers we need to output numpy arrays
    dydt_numba=loss_gain_numba(numpy.zeros((len(y_asnumpy))),reactants_numba)
    #numpy.save(filename+'dydt_numba_base', dydt_numba)
    #shutil.move(filename+'dydt_numba_base.npy','./data')
    numpy.save(filename+'dydt_numba', dydt_numba)
    
    dydt_dydt_numba=jacobian_numba(rates_numba,y_asnumpy,numpy.zeros((len(y_asnumpy),len(y_asnumpy))))
    #numpy.save(filename+'dydt_dydt_numba_base', dydt_dydt_numba)
    #shutil.move(filename+'dydt_dydt_numba_base.npy','./data')
    numpy.save(filename+'dydt_dydt_numba', dydt_dydt_numba)
    
    y_asnumpy_full=numpy.zeros((num_species+num_species_condensed*num_bins,1),)
    y_asnumpy_full[:]=1.0e12
    Pressure_gas=(y_asnumpy[0:num_species,]/NA)*8.314E+6*Model_temp #[using R]
    #numpy.save(filename+'Pressure_gas_base', Pressure_gas)
    #shutil.move(filename+'Pressure_gas_base.npy','./data')
    numpy.save(filename+'Pressure_gas', Pressure_gas)
    
    y_core=[1.0e-3]*num_bins #Will hold concentration of core material, only initialise here [molecules/cc] 
    core_density_array=[1770.0]*num_bins #[kg/m3] - need to make sure this matches core definition above
    core_density_array_asnumpy=numpy.array(core_density_array)
    core_mw=[132.14]*num_bins #[g/mol]
    core_dissociation=3.0 #Define this according to choice of core type. Please note this value might change
    y_core=(4.0/3.0)*numpy.pi*numpy.power(numpy.array(x*1.0e-6),3.0) #4/3*pi*radius^3
    y_core=y_core*numpy.array(core_density_array) #mass per particle [kg]
    y_core=y_core/(numpy.array(core_mw)*1.0e-3) #moles per particle, changing mw from g/mol to kg/mol
    y_core=y_core*NA #molecules per particle
    y_core=y_core*numpy.array(N_perbin) #molecules/cc representing each size range
    ycore_asnumpy=numpy.array(y_core)
    core_molw_asnumpy=numpy.array(core_mw)

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
    #C_g_i_t[ignore_index]=0.0
    #C_g_i_t=C_g_i_t[include_index]

    #pdb.set_trace()         
    total_SOA_mass,aw_array,size_array,dy_dt_calc = dydt_partition_fortran(y_asnumpy_full,ycore_asnumpy,core_dissociation, \
        core_mass_array,y_density_array_asnumpy,core_density_array_asnumpy,ignore_index_fortran,y_mw,Psat, \
        DStar_org_asnumpy,alpha_d_org_asnumpy,C_g_i_t,N_perbin,gamma_gas_asnumpy,Latent_heat_asnumpy,GRAV, \
        Updraft,sigma,NA,kb,Rv,R_gas,Model_temp,cp,Ra,Lv_water_vapour)
    
    #numpy.save(filename+'total_SOA_mass_base', total_SOA_mass)
    #shutil.move(filename+'total_SOA_mass_base.npy','./data')
    numpy.save(filename+'total_SOA_mass', total_SOA_mass)
    
    #numpy.save(filename+'aw_array_base', aw_array)
    #shutil.move(filename+'aw_array_base.npy','./data')
    numpy.save(filename+'aw_array', aw_array)
    
    #numpy.save(filename+'size_array_base', size_array)
    #shutil.move(filename+'size_array_base.npy','./data')
    numpy.save(filename+'size_array', size_array)
    
    #numpy.save(filename+'dy_dt_calc_base', dy_dt_calc)
    #shutil.move(filename+'dy_dt_calc_base.npy','./data')
    numpy.save(filename+'dy_dt_calc', dy_dt_calc)
    
    #pdb.set_trace()
    
    pass

class TestParsing(unittest.TestCase):

    """
    Test the output from the derived functions used within the ODE solvers
    """
    
    def test_RO2(self):
        RO2_base=numpy.load('./data/'+filename+'_RO2_base.npy')
        RO2=numpy.load(filename+'_RO2.npy')
        numpy.allclose(RO2_base,RO2, rtol=1e-04, atol=1e-08)
    
    def test_rates_fortran(self):
        rates_fortran_base=numpy.load('./data/'+filename+'rates_fortran_base.npy')
        rates_fortran=numpy.load(filename+'rates_fortran.npy')
        numpy.allclose(rates_fortran_base, rates_fortran, rtol=1e-04, atol=1e-08)
    
    def test_reactants_fortran(self):
        reactants_fortran_base=numpy.load('./data/'+filename+'reactants_fortran_base.npy')
        reactants_fortran=numpy.load(filename+'reactants_fortran.npy')
        numpy.allclose(reactants_fortran_base, reactants_fortran, rtol=1e-04, atol=1e-08)
    
    def test_dydt_fortran(self):
        dydt_fortran_base=numpy.load('./data/'+filename+'dydt_fortran_base.npy')
        dydt_fortran=numpy.load(filename+'dydt_fortran.npy')
        numpy.allclose(dydt_fortran_base, dydt_fortran, rtol=1e-04, atol=1e-08)
    
    def test_dydt_dydt_fortran(self):
        dydt_dydt_fortran_base=numpy.load('./data/'+filename+'dydt_dydt_fortran_base.npy')    
        dydt_dydt_fortran=numpy.load(filename+'dydt_dydt_fortran.npy') 
        numpy.allclose(dydt_dydt_fortran_base, dydt_dydt_fortran, rtol=1e-04, atol=1e-08)   
    
    def test_rates_numba(self):
        rates_numba_base=numpy.load('./data/'+filename+'rates_numba_base.npy')
        rates_numba=numpy.load(filename+'rates_numba.npy')
        numpy.allclose(rates_numba_base, rates_numba, rtol=1e-04, atol=1e-08)

    #def test_rates_fortran_numba(self):
    #    rates_fortran_base=numpy.load('./data/'+filename+'rates_fortran_base.npy')
    #    rates_numba_base=numpy.load('./data/'+filename+'rates_numba_base.npy')
    #    npt.assert_almost_equal(rates_fortran_base, rates_numba_base, decimal=4)
    
    def test_reactants_numba(self):
        reactants_numba_base=numpy.load('./data/'+filename+'reactants_numba_base.npy')
        reactants_numba=numpy.load(filename+'reactants_numba.npy')
        numpy.allclose(reactants_numba_base, reactants_numba, rtol=1e-04, atol=1e-08)
    
    def test_dydt_numba(self):
        dydt_numba_base=numpy.load('./data/'+filename+'dydt_numba_base.npy')
        dydt_numba=numpy.load(filename+'dydt_numba.npy')
        numpy.allclose(dydt_numba_base, dydt_numba, rtol=1e-04, atol=1e-08)
    
    def test_dydt_dydt_numba(self):
        dydt_dydt_numba_base=numpy.load('./data/'+filename+'dydt_dydt_numba_base.npy')
        dydt_dydt_numba=numpy.load(filename+'dydt_dydt_numba.npy')
        numpy.allclose(dydt_dydt_numba_base, dydt_dydt_numba, rtol=1e-04, atol=1e-08)
    
    def test_pressure_gas(self):
        Pressure_gas_base=numpy.load('./data/'+filename+'Pressure_gas_base.npy')
        Pressure_gas=numpy.load(filename+'Pressure_gas.npy')
        numpy.allclose(Pressure_gas_base, Pressure_gas, rtol=1e-04, atol=1e-08)
    
    def test_SOA_mass(self):
        total_SOA_mass_base=numpy.load('./data/'+filename+'total_SOA_mass_base.npy')
        total_SOA_mass=numpy.load(filename+'total_SOA_mass.npy')
        numpy.allclose(total_SOA_mass_base, total_SOA_mass, rtol=1e-04, atol=1e-08)
    
    def test_aw_array(self):
        aw_array_base=numpy.load('./data/'+filename+'aw_array_base.npy')
        aw_array=numpy.load(filename+'aw_array.npy')
        numpy.allclose(aw_array_base, aw_array, rtol=1e-04, atol=1e-08)
    
    def test_size_array(self):
        size_array_base=numpy.load('./data/'+filename+'size_array_base.npy')
        size_array=numpy.load(filename+'size_array.npy')
        numpy.allclose(size_array_base, size_array, rtol=1e-04, atol=1e-08)
    
    def test_dy_dt_calc(self):
        dy_dt_calc_base=numpy.load('./data/'+filename+'dy_dt_calc_base.npy')
        dy_dt_calc=numpy.load(filename+'dy_dt_calc.npy')    
        numpy.allclose(dy_dt_calc_base, dy_dt_calc, rtol=1e-04, atol=1e-08)
    

# Start of the main body of code
if __name__=='__main__':
    
    filename='MCM_APINENE'    
    
    setup(filename)
    
    # Now run the testing suite
    unittest.main()
    
