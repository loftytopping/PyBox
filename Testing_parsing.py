import numpy 
import pylab as P
import nose
import Parse_eqn_file
import rate_coeff_conversion
import MCM_constants
import collections
import pdb
from scipy.sparse import csr_matrix
import glob
import pickle
import os

def load_sparse_csr(filename):
    loader = numpy.load('loss_gain_'+filename+'.npz')
    return csr_matrix((  loader['data'], loader['indices'], loader['indptr']),shape = loader['shape'])
    

#Define variables that are accessible to all functions defined within the run_example space
temp=288.15
RH=0.3 #RH/100%
hour_of_day=12.0
time_of_day_seconds=hour_of_day*60*60
TEMP=temp
#3)Generate constants used in rate of reaction calculations
#Convert RH to concentration of water vapour molecules [this will change when in Parcel model mode]
temp_celsius=temp-273.15
# Saturation VP of water vapour
Psat=610.78*numpy.exp((temp_celsius/(temp_celsius+238.3))*17.2694)
Pw=RH*Psat
Wconc=0.002166*(Pw/(temp_celsius+273.16)) #kg/m3
#Covert from m3 to cm3
Wconc=Wconc*1.0e-6
#Convert from kg to molecules/cc
H2O=Wconc*(1.0/(18.0e-3))*6.0221409e+23

#)Now extract generic MCM rate coefficient functions
mcm_constants_dict=MCM_constants.mcm_constants(time_of_day_seconds, temp, H2O)

filename='MCM_test'

print_options=dict()
print_options['Full_eqn']=0 #Set to 1 to print details of all equations and rate coefficients during 
outputdict=Parse_eqn_file.extract_mechanism(filename+'.eqn.txt',print_options)
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
#2)Now convert the rate coefficient expressions into Python readable commands
print("Converting rate coefficient operation into Python format")
rate_dict=rate_coeff_conversion.convert_rate_mcm(rate_dict)
print("Generating indices files for use in the simulation")
Parse_eqn_file.write_rate_file(rate_dict,mcm_constants_dict)
Parse_eqn_file.write_RO2_indices(filename,species_dict2array)
Parse_eqn_file.write_reactants_indices(filename,equations,species_dict2array,rate_dict_reactants,loss_dict)
Parse_eqn_file.write_loss_gain_matrix(filename,equations,num_species,loss_dict,gain_dict,species_dict2array)
loss_gain = load_sparse_csr(filename)
reactants_indices=numpy.load(filename+'_reactant_indices.npy')
reactants_stoich_indices=numpy.load(filename+'_reactants_stoich_indices.npy')
RO2_indices=numpy.load(filename+'_RO2_indices.npy')    

    
pdb.set_trace()
