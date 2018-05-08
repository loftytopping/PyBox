import sys
import numpy 
sys.path.append('/Users/davidtopping/Code/Git_repos/UManSysProp_public/')
from umansysprop import boiling_points
from umansysprop import vapour_pressures
from umansysprop import critical_properties
from umansysprop import liquid_densities
from umansysprop import partition_models
from umansysprop import activity_coefficient_models_dev as aiomfac
from umansysprop.forms import CoreAbundanceField

def Pure_component1(num_species,species_dict,species_dict2array,Pybel_object_dict,SMILES_dict,temp,vp_method,bp_method,critical_method,density_method,ignore_vp,vp_cutoff):
    
    y_density_array=[1000.0]*num_species
    y_mw=[200.0]*num_species
    sat_vp=[100.0]*num_species
    sat_vp_org=dict() #Recorded seperately and will not include the extension for water. This is used for any checks with equilibrium partitioning predictions
    y_gas=[0.0]*num_species
    species_step=0
    ignore_index=[] #Append to this to identify any compounds that do not have automated calculation of properties.
    include_index=[]
    include_dict=dict()
    ignore_index_fortran=numpy.zeros((num_species+1),)
    print("Calculating component properties using UManSysProp")
    
    # Which boiling point method has been chosen
    boiling_point = {
        'joback_and_reid': boiling_points.joback_and_reid,
        'stein_and_brown': boiling_points.stein_and_brown,
        'nannoolal':       boiling_points.nannoolal,
        }[bp_method]
    vapour_pressure = {
        'nannoolal':            vapour_pressures.nannoolal,
        'myrdal_and_yalkowsky': vapour_pressures.myrdal_and_yalkowsky,
        # Evaporation doesn't use boiling point
        'evaporation': lambda c, t, b: vapour_pressures.evaporation(c, t),
        }[vp_method]
    # Which density method has been chosen
    critical_property = {
        'nannoolal':          critical_properties.nannoolal,
        'joback_and_reid':    critical_properties.joback_and_reid,
        }[critical_method]
    liquid_density = {
        'girolami':      lambda c, t, p: liquid_densities.girolami(c),
        'schroeder':     liquid_densities.schroeder,
        'le_bas':        liquid_densities.le_bas,
        'tyn_and_calus': liquid_densities.tyn_and_calus,
        }[density_method]
    
    for compound in species_dict.values():
        if compound in SMILES_dict.keys():
            
            # Calculate a boiling point with Nanoolal for density methods
            b1 = boiling_points.nannoolal(Pybel_object_dict[SMILES_dict[compound]])
            y_density_array[species_dict2array[compound]]=liquid_density(Pybel_object_dict[SMILES_dict[compound]], temp, critical_property(Pybel_object_dict[SMILES_dict[compound]], b1))
            
            #y_density_array[species_dict2array[compound]]=(liquid_densities.girolami(Pybel_object_dict[SMILES_dict[compound]])*1.0E3) #Convert from g/cc to kg/m3
            #y_density_array.append(1400.0)
            y_mw[species_dict2array[compound]]=(Pybel_object_dict[SMILES_dict[compound]].molwt)
            #In the following you will need to select which vapour pressure method you like.
            
            # Calculate boiling point
            b = boiling_point(Pybel_object_dict[SMILES_dict[compound]])
            
            sat_vp[species_dict2array[compound]]=vapour_pressure(Pybel_object_dict[SMILES_dict[compound]], temp,b)
            
            if ignore_vp is True:
                if sat_vp[species_dict2array[compound]] > vp_cutoff:
                    ignore_index.append(species_dict2array[compound])
                    ignore_index_fortran[species_dict2array[compound]]=1.0
                else:
                    include_index.append(species_dict2array[compound])
                    
            #sat_vp[species_dict2array[compound]]=(vapour_pressures.nannoolal(Pybel_object_dict[SMILES_dict[compound]], temp, boiling_points.nannoolal(Pybel_object_dict[SMILES_dict[compound]])))
            #sat_vp_org[Pybel_object_dict[SMILES_dict[compound]]]=vapour_pressures.nannoolal(Pybel_object_dict[SMILES_dict[compound]], temp, boiling_points.nannoolal(Pybel_object_dict[SMILES_dict[compound]]))
            #y_gas.append(concentration_array[species_step])
        else:
            ignore_index.append(species_dict2array[compound])
            ignore_index_fortran[species_dict2array[compound]]=1.0
        species_step+=1
    Delta_H=[0.0]*num_species #Dont prescribe an enthalpy of vaporisation. For non-VBS model simulations with varying temperature, we recalculate using GCM
    Latent_heat_gas=[0.0]*num_species # - Still need to account for any latent heat release [Future work]

    return_dict=dict()
    return_dict['y_density_array']=y_density_array
    return_dict['y_mw']=y_mw
    return_dict['sat_vp']=sat_vp
    return_dict['Delta_H']=Delta_H
    return_dict['Latent_heat_gas']=Latent_heat_gas   
    return_dict['ignore_index']=ignore_index 
    return_dict['ignore_index_fortran']=ignore_index_fortran
    return_dict['include_index']=include_index

    return return_dict
    
def Pure_component2(num_species,y_mw,R_gas,temp):
    #----------------------------------------------------------------------------
    #2e) Kinetic properties of all species ######
    #sat_vp = [-12.0]*num_species # np.random.random_sample((num_species,))*-12.0 #log10 atm
    alpha_d_org = [0.1]*num_species # Accomodation coefficient [including water]
    # - Molecular diffusion coeffficient of each molecula in air. [cm2/s] - 
    #We usean approximation to initialise this, but these values will also depend
    #on the accomodation coefficient and size of particles condensing to. Following the
    #approach in Topping et al (2013), diffusivity is calculated as 
    DStar_org = 1.9E0*numpy.power(numpy.array(y_mw),-2.0E0/3.0E0)
    # - Mean thermal velocity of each molecule [m/s] - 
    mean_them_vel=numpy.power((8.0E0*R_gas*temp)/(numpy.pi*(numpy.array(y_mw)*1.0E-3)),0.5E0)
    # - Mean free path for each molecule [m] - 
    gamma_gas = ((3.0E0*DStar_org)/(mean_them_vel*1.0E2))*1.0E-2
    ##############################################
    
    return_dict=dict()
    return_dict['alpha_d_org']=alpha_d_org
    return_dict['DStar_org']=DStar_org
    return_dict['mean_them_vel']=mean_them_vel
    return_dict['gamma_gas']=gamma_gas
    
    return return_dict