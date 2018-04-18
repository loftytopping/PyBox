##########################################################################################
#                                                                                        #
#    Scripts to build model from chemical mechanism file and/or extract species          #
#                                                                                        #
#    Copyright (C) 2017  David Topping : david.topping@manchester.ac.uk                  #
#                                      : davetopp80@gmail.com                            #
#    Personal website: davetoppingsci.com                                                #
#                                                                                        #
#    This program does not have a license, meaning the deault copyright law applies.     #
#    Only users who have access to the private repository that holds this file may       #
#    use it or develop it, but may not distribute it.                                    #
#                                                                                        #
#                                                                                        #
##########################################################################################

# VERSION 0.9
# - Not fit for public release.

# Last modification 30/11/17
# Please note - All mechanism work here is based on the MCM. This includes definitions of contribution
# to RO2 concetrations etc.

import collections
import pdb
import re
import pybel # Needed for calculating olecular properties in aerosol model [Parse_eqn_file used for both gas and aerosol]
import numpy as np
import datetime
from scipy.sparse import lil_matrix
import os
from xml.dom import minidom
import multiprocessing

def extract_mechanism(filename,print_options): #filename='saprc99.eqn' #Change to the filename of interest

    # The following steps are used:
    # 1) Read the equation file

    # --------------------------------------------------------------------------------------------
    #[1] - Read the species equation file. Count how many equations there are
    # Find all instances of curly brackets {} which enclose an equation number

    print("Opening file ", str(filename), "for parsing")
  
    text=open(filename,'rU')
    with open (filename, "r") as myfile:
        data=myfile.read().replace('\n', '')

    char1='{'
    char2='}'

    max_equations=int(re.findall(r"\{(.*?)\}",data)[-1].strip('.')) 
    print("Calculating total number of equations = ",max_equations)

    # --------------------------------------------------------------------------------------------


    # --------------------------------------------------------------------------------------------
    #[2] - Cycle through all equations and extract:
    #        - Reactants and stochiometric coefficients
    #        - Products and stochiometric coefficients
    #        - Rate coefficients:
    #             = Coefficients and typical forms used in MCM models
    #
    # This information is stored in dictionaries 

    # The following loop looks into the .eqn file and finds the text sat between each equation. 
    # For example, the first step pulls out the equation information for equation 1, as defined 
    # by the KPP standard. We cannot just parse line by line since equations can span multiple 
    # lines according to KPP files. Each equation 'line' can be defined to start after a '}' and
    # end with a ';' For eaxmple:
    #{1}  NO2 + hv = NO + O :	0.533 ;
    #{2}  O + O2 = O3       :	2.183E-5 ;
    #{3}  NO + 2.0O3 = NO2 + O2 :	26.59 ;
    #{4}  RH + OH = RO2 + H2O :	3.775E+3 ;
    #
    #Dictionaries used to store information about reactants/products/equations
    reaction_dict=collections.defaultdict(lambda: collections.defaultdict())
    rate_dict=collections.defaultdict(lambda: collections.defaultdict())
    rate_def=collections.defaultdict()
    loss_dict=collections.defaultdict(lambda: collections.defaultdict())
    gain_dict=collections.defaultdict(lambda: collections.defaultdict())
    stoich_dict=collections.defaultdict(lambda: collections.defaultdict())
    rate_dict_reactants=collections.defaultdict(lambda: collections.defaultdict())
    species_dict=collections.defaultdict()
    species_dict2array=collections.defaultdict()
    species_hess_data=collections.defaultdict()
    species_hess_loss_data=collections.defaultdict()
    species_hess_gain_data=collections.defaultdict()

    # Extract all lines with equations on as a list
    eqn_list=re.findall(r"\}(.*?)\;",data)

    #Create an integer that stores number of unque species
    species_step=0

    #pdb.set_trace()

    print("Parsing each equation")
	
    for equation_step in range(max_equations):

        equation_full=eqn_list[equation_step] #pulls out text in between {n}..{n+1}
        equation_full=re.sub('\t','',equation_full)
        equation=equation_full.split(':',1)[0].split('=',1) #split the line into reactants and products
        reactants=equation[0].split('+') # extract content to the left of the previous split [reactants]
        reactants= [x.strip(' ') for x in reactants] #strip away all whitespace
        products=equation[1].split('+') # extract content to the right of the previous split [products]
        products = [x.strip(' ') for x in products] #strip away all whitespace
        #pdb.set_trace()
        #At the moment, we have not seperated the reactant/product from its stochiometric value
        ##rate_full=equation_full.split(':',1)[1].rsplit('\t',1)[1].split(';',1)[0] #pulls out the rate 
        rate_full=equation_full.split(':',1)[1]
        #but as a text string
        rate_full=rate_full.strip() #[x.strip(' ') for x in rate_full]
        rate_dict[equation_step]=rate_full
        #This assumes everyline, as in KPP, finishes with a ';' character

        #pdb.set_trace()

        #print "Extracting equation :"
        #print equation_step
        #if print_options['Full_eqn']==1:
        #    print "Full equation extracted"
        #    print equation_full
        #if print_options['Full_eqn']==1:
        #    print "Rate extracted :"
        #    print rate_full

        # Now cycle through all reactants and products, splitting the unique specie from its stochiometry
        # This information is then stored in dictionaries for use in the ODE solver

        # At the moment the reactants, and products, include joint stoichiometric information. 
        # This means we need to identify these numbers and then split the string again to ensure
        # the specie always remains unique. Thus, we may have saved something like
        # '2.0NO2' or '5ISOPOOH'. The use of integer versus float can vary so have to assume no care taken
        # in being consistent 

        reactant_step=0 #used to identify reactants by number, for any given reaction
        product_step=0

        for reactant in reactants: #This are the 'reactants' in this equation
            reactant=reactant.split()[0] #remove all tables, newlines, whitespace
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # - Extract stochiometry and unique reactant identifier
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            try: #Now try to extract a stochiometric coefficient if given
                temp=re.findall(r"[-+]?\d*\.\d+|\d+|\d+",reactant) #This extracts all numbers either side 
                                                                   #of some text. 
                stoich=temp[0] #This selects the first number extracted, if any.
                # Now need to work out if this value is before the variable
                # If after, we ignore this. EG. If '2' it could come from '2NO2' or just 'NO2'
                # If len(temp)==1 then we only have one number and can proceed with the following
                if len(temp)==1:
                    if reactant.index(stoich) == 0 :
                        reactant=reactant.split(stoich,1)[1] 
                        stoich=float(stoich)
                    else:
                        stoich=1.0
                elif len(temp)>1:
                    #If this is the case, we need to ensure the reactant extraction is unique. For example
                    #If the string is '2NO2' the above procedure extracts 'NO' as the unique reactant. 
                    #We therefore need to ensure that the reactant is 'NO2'. To do this we cut the value
                    #in temp[0] away from the original string. Lets assume that we can attach the first
                    #part of the text with the second number in temp. Thus
                    if reactant.index(stoich) == 0 :
                        reactant=reactant.split(stoich,1)[1]+temp[1]
                        stoich=float(stoich)
                    else:
                        stoich=1.0
            except:
                stoich=1.0
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # - Store stoichiometry, species flags and hessian info in dictionaries
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #Now store the stoichiometry and reactant in dictionaries
            if reactant not in ['hv']:
                stoich_dict[equation_step][reactant_step]=stoich
                rate_dict_reactants[equation_step][reactant_step]=reactant 
                # -- Update species dictionaries --
                if reactant not in species_dict.values(): #check to see if entry already exists
                    species_dict[species_step]=reactant # useful for checking all parsed species
                    species_dict2array[reactant]=species_step #useful for converting a dict to array
                    species_step+=1
                # -- Update hessian dictionaries --
                if reactant not in species_hess_loss_data:
                    species_hess_loss_data[reactant]=[]
                species_hess_loss_data[reactant].append(equation_step) #so from this we can work out a dx/dy
                # -- Update loss dictionaries -- 
                if equation_step in loss_dict[reactant]:
                    loss_dict[reactant][equation_step]+=stoich
                else:
                    loss_dict[reactant][equation_step]=stoich

            reactant_step+=1

        for product in products: #This are the 'reactants' in this equation
            product=product.split()[0] #remove all tables, newlines, whitespace
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # - Extract stochiometry and unique product identifier
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            try: #Now try to extract a stochiometric coefficient if given
                temp=re.findall(r"[-+]?\d*\.\d+|\d+|\d+",product) #This extracts all numbers either side 
                                                                   #of some text. 
                stoich=temp[0] #This selects the first number extracted, if any.
                # Now need to work out if this value is before the variable
                # If after, we ignore this. EG. If '2' it could come from '2NO2' or just 'NO2'
                # If len(temp)==1 then we only have one number and can proceed with the following
                if len(temp)==1:
                    if product.index(stoich) == 0 :
                        product=product.split(stoich,1)[1] 
                        stoich=float(stoich)
                    else:
                        stoich=1.0
                elif len(temp)>1:
                    #If this is the case, we need to ensure the reactant extraction is unique. For example
                    #If the string is '2NO2' the above procedure extracts 'NO' as the unique reactant. 
                    #We therefore need to ensure that the reactant is 'NO2'. To do this we cut the value
                    #in temp[0] away from the original string. Lets assume that we can attach the first
                    #part of the text with the second number in temp. Thus
                    if product.index(stoich) == 0 :
                        product=product.split(stoich,1)[1]+temp[1]
                        stoich=float(stoich)
                    else:
                        stoich=1.0
            except:
                stoich=1.0
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # - Store stoichiometry, species flags and hessian info in dictionaries
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #Now store the stoichiometry and reactant in dictionaries
            if product not in ['hv']:
                # -- Update species dictionaries --
                if product not in species_dict.values(): #check to see if entry already exists
                    species_dict[species_step]=product # useful for checking all parsed species
                    species_dict2array[product]=species_step #useful for converting a dict to array
                    species_step+=1
                # -- Update hessian dictionaries --
                if product not in species_hess_gain_data:
                    species_hess_gain_data[product]=[]
                species_hess_gain_data[product].append(equation_step) #so from this we can work out a dx/dy
                # -- Update loss dictionaries -- 
                if equation_step in gain_dict[reactant]:
                    gain_dict[product][equation_step]+=stoich
                else:
                    gain_dict[product][equation_step]=stoich

            product_step+=1

    #pdb.set_trace()
    # --------------------------------------------------------------------------------------------
    print("Total number of species = ", len(species_dict.keys()))
    print("Saving all equation information to dictionaries")
    output_dict=dict()
    output_dict['reaction_dict']=reaction_dict
    output_dict['rate_dict']=rate_dict
    output_dict['rate_dict_reactants']=rate_dict_reactants
    output_dict['rate_def']=rate_def
    output_dict['loss_dict']=loss_dict
    output_dict['gain_dict']=gain_dict
    output_dict['stoich_dict']=stoich_dict
    output_dict['species_dict']=species_dict
    output_dict['species_dict2array']=species_dict2array
    output_dict['species_hess_data']=species_hess_data
    output_dict['species_hess_loss_data']=species_hess_loss_data
    output_dict['species_hess_gain_data']=species_hess_gain_data
    output_dict['max_equations']=max_equations
    
    #Now you need to create a new function file for calculating rate coefficients. This will store the 
    #equation rates as 'hardcoded' representations. The benefit is not having to loop through dictionary
    #entries that can become very slow
    #call another funtion for this
    #rate_function(rate_dict)
    
    #Also create the function for the gas phase 

    return output_dict

def extract_smiles_species(output_dict, SMILES_filename):
    #Here we map the SMILES to a species name. For the MCM, this information is stored in an XML file.
    #This is the name of the filename you need to provide.
    #This information would be seperate from the information provided in extract species, and you would run this function instead
    #For development purposes the default XML file has been MCM331.xml
    
    #output_dict is generated through a moulde that reads an equation file
    
    print('Mapping species names to SMILES and Pybel objects')
    
    SMILES_dict=collections.defaultdict()
    Pybel_object_dict=collections.defaultdict()
    
    species_dict=output_dict['species_dict']
    
    # parse an xml file by name - use this to store species name versus SMILES
    dom = minidom.parse(SMILES_filename)
    for node in dom.getElementsByTagName('species'):
        species_name=node.attributes['species_name'].value
        if species_name not in ['OH','HO2','O3','hv']: 
            #Now check to see if the species is in your dictionary extracted from the mechanism
            if species_name in species_dict.values():
                #pdb.set_trace()
                #species_SMILES=node.attributes['smiles'].value
                #Right, the following is an overkill way t extract the SMILES. If the XML was written correctly, and the SMILES were in quotation marks, the node. method would be all that is needed
                try:
                    species_SMILES=re.search(r'<smiles>(.*?)</smiles>',node.getElementsByTagName('smiles').item(0).toxml()).group(1)
                except:
                    print('No SMILES entry for species ',species_name)
                SMILES_dict[species_name]=species_SMILES
                Pybel_object=pybel.readstring('smi',species_SMILES)
                Pybel_object_dict[species_SMILES]=Pybel_object
    
    output_dict['SMILES_dict']=SMILES_dict
    output_dict['Pybel_object_dict']=Pybel_object_dict
    
    #pdb.set_trace()
    
    #return species mapped to SMILES - bare in mind this might not be complete.
    
    return output_dict
        

def extract_species(filename):
    
    #This function is used to test box-models by extracting from a pre-defined chemical mechanism snapshot
    smiles_array=[]
    concentration_array=[]
    SMILES_dict=collections.defaultdict()
    Pybel_object_dict=collections.defaultdict()
    Pybel_object_activity=collections.defaultdict()
    
    species_dict=collections.defaultdict()
    species_dict2array=collections.defaultdict()

    #text=open(filename,'rU')
    
    #Create an integer that stores number of unique species
    species_step=0
   
    for line in open(filename, "r"):
        
        input = line.split()
        # Keep a list of the information
        smiles=input[0]
        #pdb.set_trace()
        if smiles not in species_dict.values():
            smiles_array.append(smiles)
            concentration_array.append(float(input[1]))
            # Now create Pybel objects which are used in all property predictive techniquexs
            Pybel_object=pybel.readstring('smi',input[0])
            #Pybel_object=pybel.readstring(b'smi',str(input[0],'utf-8'))
            Pybel_object_dict[smiles]=Pybel_object
            Pybel_object_activity[Pybel_object]=float(input[1])
            species_dict[species_step]=smiles # useful for checking all parsed species
            species_dict2array[Pybel_object]=species_step #useful for converting a dict to array
            species_step+=1
            
    print("Total number of species = ", species_step)
    print("Saving all equation information to dictionaries")
    output_dict=dict()
    output_dict['species_dict']=species_dict
    output_dict['species_dict2array']=species_dict2array
    output_dict['Pybel_object_dict']=Pybel_object_dict
    output_dict['smiles_array']=smiles_array
    output_dict['concentration_array']=concentration_array
    output_dict['species_number']=species_step
    output_dict['Pybel_object_activity']=Pybel_object_activity

    return output_dict
    
def write_rate_file(filename,rate_dict):
    
    #Put all of the rate coefficient functional forms into a new python file.
    f = open('Rate_coefficients.py','w')
    f.write('##################################################################################################### \n') # python will convert \n to os.linesep
    f.write('# Python function to hold expressions for calculating rate coefficients for a given equation number # \n') # python will convert \n to os.linesep
    f.write('#    Copyright (C) 2017  David Topping : david.topping@manchester.ac.uk                             # \n')
    f.write('#                                      : davetopp80@gmail.com                                       # \n')
    f.write('#    Personal website: davetoppingsci.com                                                           # \n')
    f.write('#                                                                                                   # \n')
    f.write('#    This program does not have a license, meaning the deault copyright law applies.                # \n')
    f.write('#    Only users who have access to the private repository that holds this file may                  # \n')
    f.write('#    use it or develop it, but may not distribute it.                                               # \n')
    f.write('#                                                                                                   # \n')
    f.write('#                                                                                                   # \n')
    f.write('##################################################################################################### \n')    
    f.write('# File created at %s \n' % datetime.datetime.now()) # python will convert \n to os.linesep
    f.write('\n') 
    f.write('import numpy \n') 
    #f.write('import numba')
    f.write('\n') 
    #f.write('@numba.jit(nopython=True)\n')
    f.write('def evaluate_rates(ttime,RO2,H2O,temp):\n') 
    # Now cycle through all of the mcm_constants_dict values.
    # Please note this includes photolysis rates that will change with time of day or condition in the chamber. You will
    # need to modify this potentially.
    f.write('    # Creating reference to constant values used in rate expressions\n') 


    f.write('    # constants used in calculation of reaction rates\n') 
    f.write('    M  = 2.55E+19  #Check this against pressure assumed in model\n') 
    f.write('    N2 = 0.79*M\n') 
    f.write('    O2 = 0.2095*M\n') 

    f.write('    # kro2no : ro2      + no      = ro      + no2\n') 
    f.write('    #        : ro2      + no      = rono2\n') 
    f.write('    # iupac 1992\n') 
    f.write('    KRONO2    = 2.40E-12*numpy.exp(360.0/temp)\n') 
	
    f.write('    # kro2ho2: ro2      + ho2     = rooh    + o2\n') 
    f.write('    # mcm protocol [1997]\n') 
    f.write('    KRO2HO2   = 2.91E-13*numpy.exp(1300.0/temp)\n') 
	
    f.write('    # kapho2 : rcoo2    + ho2     = products\n') 
    f.write('    # mcm protocol [1997]\n') 
    f.write('    KAPHO2    = 4.30E-13*numpy.exp(1040.0/temp)\n') 
	
    f.write('    # kapno  : rcoo2    + no      = products\n') 
    f.write('    # mej [1998]\n') 
    f.write('    KAPNO = 8.10E-12*numpy.exp(270.0/temp)\n') 
	
    f.write('    # kro2no3: ro2      + no3     = products\n') 
    f.write('    # mcm protocol [1997]\n') 
    f.write('    KRO2NO3   = 2.50E-12\n') 
	
    f.write('    # kno3al : no3      + rcho    = rcoo2   + hno3\n') 
    f.write('    # mcm protocol [1997]\n') 
    f.write('    KNO3AL    = 1.44E-12*numpy.exp(-1862.0/temp)\n') 

    f.write('    # kdec   : ro                 = products\n') 
    f.write('    # mcm protocol [1997]\n') 
    f.write('    KDEC      = 1.00E+06\n') 

    f.write('    KROSEC = 1.50E-14*numpy.exp(-200.0/temp)\n') 

    f.write('    KALKOXY=3.70E-14*numpy.exp(-460.0/temp)*O2\n') 

    f.write('    KALKPXY=1.80E-14*numpy.exp(-260.0/temp)*O2\n') 

    f.write('    # -------------------------------------------------------------------\n') 
    f.write('    # complex reactions\n') 
    f.write('    # -------------------------------------------------------------------\n') 

    f.write('    # kfpan kbpan\n') 
    f.write('    # formation and decomposition of pan\n') 
    f.write('    # iupac 2001 (mcmv3.2)\n') 
    f.write('    kc0     = 2.70E-28*M*(temp/300.0)**(-7.1)\n') 
    f.write('    kci     = 1.21E-11*(temp/300.0)**(-0.9)\n') 
    f.write('    krc     = kc0/kci\n') 
    f.write('    fcc     = 0.30\n') 
    f.write('    nc      = 0.75-(1.27*numpy.log10(fcc))\n') 
    f.write('    fc      = 10**(numpy.log10(fcc)/(1.0+((numpy.log10(krc))/nc)**2.0))\n') 
    f.write('    KFPAN   = (kc0*kci)*fc/(kc0+kci)\n') 

    f.write('    kd0     = 4.90E-03*M*numpy.exp(-12100.0/temp)\n') 
    f.write('    kdi     = 5.40E+16*numpy.exp(-13830.0/temp)\n') 
    f.write('    krd     = kd0/kdi\n') 
    f.write('    fcd     = 0.30\n') 
    f.write('    ncd     = 0.75-(1.27*numpy.log10(fcd))\n') 
    f.write('    fd      = 10.0**(numpy.log10(fcd)/(1.0+((numpy.log10(krd))/ncd)**2.0))\n') 
    f.write('    KBPAN   = (kd0*kdi)*fd/(kd0+kdi)\n') 

    f.write('    # kmt01  : o        + no      = no2\n') 
    f.write('    # iupac 2001 (mcmv3.2)\n') 
    f.write('    k10     = 1.00E-31*M*(temp/300.0)**(-1.6)\n') 

    f.write('    k1i     = 3.00E-11*(temp/300.0)**(0.3)\n') 
    f.write('    kr1     = k10/k1i\n') 
    f.write('    fc1     = 0.85\n') 
    f.write('    nc1     = 0.75-(1.27*numpy.log10(fc1))\n') 
    f.write('    f1      = 10.0**(numpy.log10(fc1)/(1.0+((numpy.log10(kr1)/nc1))**2.0))\n') 
    f.write('    KMT01   = (k10*k1i)*f1/(k10+k1i)\n') 

    f.write('    # kmt02  : o        + no2     = no3\n') 
    f.write('    # iupac 2001 (mcmv3.2)\n') 
    f.write('    k20     = 1.30E-31*M*(temp/300.0)**(-1.5)\n') 
    f.write('    k2i     = 2.30E-11*(temp/300.0)**(0.24)\n') 
    f.write('    kr2     = k20/k2i\n') 
    f.write('    fc2     = 0.6\n') 
    f.write('    nc2     = 0.75-(1.27*numpy.log10(fc2))\n') 
    f.write('    f2      = 10.0**(numpy.log10(fc2)/(1.0+((numpy.log10(kr2)/nc2))**2.0))\n') 
    f.write('    KMT02   = (k20*k2i)*f2/(k20+k2i)\n') 

    f.write('    # kmt03  : no2      + no3     = n2o5\n') 
    f.write('    # iupac 2006, mcmv3.2\n') 
    f.write('    k30     = 3.60E-30*M*(temp/300.0)**(-4.1)\n') 
    f.write('    k3i     = 1.90E-12*(temp/300.0)**(0.2)\n') 
    f.write('    kr3     = k30/k3i\n') 
    f.write('    fc3     = 0.35\n') 
    f.write('    nc3     = 0.75-(1.27*numpy.log10(fc3))\n') 
    f.write('    f3      = 10.0**(numpy.log10(fc3)/(1.0+((numpy.log10(kr3)/nc3))**2.0))\n') 
    f.write('    KMT03   = (k30*k3i)*f3/(k30+k3i)\n') 

    f.write('    # kmt04  : n2o5               = no2     + no3\n') 
    f.write('    # iupac 2006, mcmv3.2\n') 
    f.write('    k40     = 1.30E-03*M*(temp/300.0)**(-3.5)*numpy.exp(-11000.0/temp)\n') 
    f.write('    k4i     = 9.70E+14*(temp/300.0)**(0.1)*numpy.exp(-11080.0/temp)\n') 
    f.write('    kr4     = k40/k4i\n') 
    f.write('    fc4     = 0.35\n') 
    f.write('    nc4     = 0.75-(1.27*numpy.log10(fc4))\n') 
    f.write('    f4      = 10.0**(numpy.log10(fc4)/(1+((numpy.log10(kr4)/nc4))**2.0))\n') 
    f.write('    KMT04   = (k40*k4i)*f4/(k40+k4i)\n') 

    f.write('    # kmt05  : oh       + co(+o2) = ho2     + co2\n') 
    f.write('    # iupac 2006\n') 
    f.write('    KMT05  = (1.0 + (M/4.2E19))\n') 

    f.write('    # kmt06  : ho2      + ho2     = h2o2    + o2\n') 
    f.write('    # water enhancement factor\n') 
    f.write('    # iupac 1992\n') 
    f.write('    KMT06  = 1.0 + (1.40E-21*numpy.exp(2200.0/temp)*H2O)\n') 

    f.write('    # kmt06  = 1.0 + (2.00E-25*numpy.exp(4670.0/temp)*h2o)\n') 
    f.write('    # S+R 2005 values\n') 

    f.write('    # kmt07  : oh       + no      = hono\n') 

    f.write('    # iupac 2006, mcmv3.2\n') 
    f.write('    k70     = 7.40E-31*M*(temp/300.0)**(-2.4)\n') 
    f.write('    k7i     = 3.30E-11*(temp/300.0)**(-0.3)\n') 
    f.write('    kr7     = k70/k7i\n') 
    f.write('    fc7     = numpy.exp(-temp/1420.0)\n') 
    f.write('    nc7     = 0.75-(1.27*numpy.log10(fc7))\n') 
    f.write('    f7      = 10.0**(numpy.log10(fc7)/(1+((numpy.log10(kr7)/nc7))**2.0))\n') 
    f.write('    KMT07   = (k70*k7i)*f7/(k70+k7i)\n') 

    f.write('    # kmt08  : oh       + no2     = hno3\n') 

    f.write('    # iupac 2006, mcmv3.2\n') 
    f.write('    k80     = 3.30E-30*M*(temp/300.0)**(-3.0)\n') 
    f.write('    k8i     = 4.10E-11\n') 
    f.write('    kr8     = k80/k8i\n') 
    f.write('    fc8     = 0.4\n') 
    f.write('    nc8     = 0.75-(1.27*numpy.log10(fc8))\n') 
    f.write('    f8      = 10.0**(numpy.log10(fc8)/(1.0+((numpy.log10(kr8)/nc8))**2.0))\n') 
    f.write('    KMT08   = (k80*k8i)*f8/(k80+k8i)\n') 

    f.write('    # kmt09  : ho2      + no2     = ho2no2\n') 
    f.write('    # iupac 1997, mcmv3.2\n') 
    
    f.write('    k90     = 1.80E-31*M*(temp/300.0)**(-3.2)\n') 
    f.write('    k9i     = 4.70E-12\n') 
    f.write('    kr9     = k90/k9i\n') 
    f.write('    fc9     = 0.6\n') 
    f.write('    nc9     = 0.75-(1.27*numpy.log10(fc9))\n') 
    f.write('    f9      = 10.0**(numpy.log10(fc9)/(1.0+((numpy.log10(kr9)/nc9))**2.0))\n') 
    f.write('    KMT09   = (k90*k9i)*f9/(k90+k9i)\n') 

    f.write('    # kmt10  : ho2no2             = ho2     + no2\n') 
    f.write('    # iupac 1997, mcmv3.2\n') 

    f.write('    k100     = 4.10E-05*M*numpy.exp(-10650.0/temp)\n') 
    f.write('    k10i     = 4.80E+15*numpy.exp(-11170.0/temp)\n') 
    f.write('    kr10     = k100/k10i\n') 
    f.write('    fc10     = 0.6\n') 
    f.write('    nc10     = 0.75-(1.27*numpy.log10(fc10))\n') 
    f.write('    f10      = 10.0**(numpy.log10(fc10)/(1.0+((numpy.log10(kr10)/nc10))**2.0))\n') 
    f.write('    KMT10    = (k100*k10i)*f10/(k100+k10i)\n') 

    f.write('    # kmt11  : oh       + hno3    = h2o     + no3\n') 
    f.write('    # iupac 2006, mcmv3.2\n') 

    f.write('    k1     = 2.40E-14*numpy.exp(460.0/temp)\n') 
    f.write('    k3     = 6.50E-34*numpy.exp(1335.0/temp)\n') 
    f.write('    k4     = 2.70E-17*numpy.exp(2199.0/temp)\n') 
    f.write('    k2     = (k3*M)/(1.0+(k3*M/k4))\n') 
    f.write('    KMT11  = k1 + k2\n') 

    f.write('    # kmt12 iupac 2006, mcmv3.2\n') 

    f.write('    k120 = 4.50E-31*((temp/300.0)**(-3.9))*M\n') 
    f.write('    k12i = 1.30E-12*((temp/300.0)**(-0.7))\n') 
    f.write('    kr12 = k120/k12i\n') 
    f.write('    fc12 = 0.525\n') 
    f.write('    nc12 = 0.75-(1.27*numpy.log10(fc12))\n') 
    f.write('    f12  = 10.0**(numpy.log10(fc12)/(1.0+((numpy.log10(kr12)/nc12))**2.0))\n') 
    f.write('    KMT12    = (k120*k12i)*f12/(k120+k12i)\n') 

    f.write('    # kmt13  : ch3o2    + no2     = ch3o2no2\n') 
    f.write('    # iupac 2006\n') 

    f.write('    k130     = 2.50E-30*((temp/300.0)**(-5.5))*M\n') 
    f.write('    k13i     = 1.80E-11\n') 
    f.write('    kr13     = k130/k13i\n') 
    f.write('    fc13     = 0.36\n') 
    f.write('    nc13     = 0.75-(1.27*numpy.log10(fc13))\n') 
    f.write('    f13      = 10.0**(numpy.log10(fc13)/(1.0+((numpy.log10(kr13)/nc13))**2.0))\n') 
    f.write('    KMT13    = (k130*k13i)*f13/(k130+k13i)\n') 

    f.write('    # kmt14  : ch3o2no2           = ch3o2   + no2\n') 
    f.write('    # iupac 2006, mcmv3.2\n') 

    f.write('    k140     = 9.00E-05*numpy.exp(-9690.0/temp)*M\n') 
    f.write('    k14i     = 1.10E+16*numpy.exp(-10560.0/temp)\n') 
    f.write('    kr14     = k140/k14i\n') 
    f.write('    fc14     = 0.4\n') 
    f.write('    nc14     = 0.75-(1.27*numpy.log10(fc14))\n') 
    f.write('    f14      = 10.0**(numpy.log10(fc14)/(1.0+((numpy.log10(kr14)/nc14))**2.0))\n') 
    f.write('    KMT14    = (k140*k14i)*f14/(k140+k14i)\n') 

    f.write('    # kmt15 iupac 2006, mcmv3.2\n') 

    f.write('    k150 = 8.60E-29*((temp/300.0)**(-3.1))*M\n') 
    f.write('    k15i = 9.00E-12*((temp/300.0)**(-0.85))\n') 
    f.write('    kr15 = k150/k15i\n') 
    f.write('    fc15 = 0.48\n') 
    f.write('    nc15 = 0.75-(1.27*numpy.log10(fc15))\n') 
    f.write('    f15  = 10.0**(numpy.log10(fc15)/(1.0+((numpy.log10(kr15)/nc15))**2.0))\n') 
    f.write('    KMT15 = (k150*k15i)*f15/(k150+k15i)\n') 

    f.write('    # kmt16  :  oh  +  c3h6\n') 
    f.write('    # iupac 2006\n') 

    f.write('    k160     = 8.00E-27*((temp/300.0)**(-3.5))*M\n') 
    f.write('    k16i     = 3.00E-11*((temp/300.0)**(-1.0))\n') 
    f.write('    kr16     = k160/k16i\n') 
    f.write('    fc16     = 0.5\n') 
    f.write('    nc16     = 0.75-(1.27*numpy.log10(fc16))\n') 
    f.write('    f16      = 10.0**(numpy.log10(fc16)/(1.0+((numpy.log10(kr16)/nc16))**2.0))\n') 
    f.write('    KMT16    = (k160*k16i)*f16/(k160+k16i)\n') 

    f.write('    # kmt17 iupac 2006\n') 

    f.write('    k170 = 5.00E-30*((temp/300.0)**(-1.5))*M\n') 
    f.write('    k17i = 1.00E-12\n') 
    f.write('    kr17 = k170/k17i\n') 
    f.write('    fc17 = (0.17*numpy.exp(-51./temp))+numpy.exp(-1.0*temp/204.)\n') 
    f.write('    nc17 = 0.75-(1.27*numpy.log10(fc17))\n') 
    f.write('    f17  = 10.0**(numpy.log10(fc17)/(1.0+((numpy.log10(kr17)/nc17))**2.0))\n') 
    f.write('    KMT17 = (k170*k17i)*f17/(k170+k17i)\n') 

    f.write('    #Taken from yet another MCM version, will need converting to lower case\n') 
    f.write('    KRO2NO = 2.7E-12*numpy.exp(360/temp)\n') 

    f.write('    KRO2HO2 = 2.91E-13*numpy.exp(1300/temp)\n') 

    f.write('    KAPHO2 = 5.2E-13*numpy.exp(980/temp)\n') 

    f.write('    KAPNO = 7.5E-12*numpy.exp(290/temp)\n') 

    f.write('    KRO2NO3 = 2.3E-12\n') 

    f.write('    KNO3AL = 1.4E-12*numpy.exp(-1860/temp)\n') 

    f.write('    KDEC = 1.00E+06\n') 

    f.write('    KROPRIM = 2.50E-14*numpy.exp(-300/temp)\n') 

    f.write('    KROSEC = 2.50E-14*numpy.exp(-300/temp)\n') 

    f.write('    KCH3O2 = 1.03E-13*numpy.exp(365/temp)\n') 

    f.write('    K298CH3O2 = 3.5E-13\n') 

    f.write('    KMT18 = 9.5E-39*O2*numpy.exp(5270/temp)/(1+7.5E-29*O2*numpy.exp(5610/temp))\n') 

    f.write('    KROPRIM  = 3.70E-14*numpy.exp(-460.0/temp)\n') 

    f.write('    KROSEC   = 1.80E-14*numpy.exp(-260.0/temp)\n') 

    f.write('    # ************************************************************************\n') 
    f.write('    # define photolysis reaction rates using derwent method from mcm2box.fac\n') 
    f.write('    # ************************************************************************\n') 

    f.write('    # solar declination angle \n') 
    f.write('    dec = 23.79\n') 
    f.write('    # latitude\n') 
    f.write('    lat = 50.0\n') 
    f.write('    pi = 4.0*numpy.arctan(1.0)\n') 
    f.write('    # local hour angle - representing time of day\n') 
    f.write('    lha = (1.0+ttime/4.32E+4)*pi\n') 
    f.write('    radian = 180.0/pi\n') 
    f.write('    lat = lat/radian\n') 
    f.write('    dec = dec/radian\n') 
    f.write('    theta = numpy.arccos(numpy.cos(lha)*numpy.cos(dec)*numpy.cos(lat)+numpy.sin(dec)*numpy.sin(lat))\n') 
    f.write('    sinld = numpy.sin(lat)*numpy.sin(dec)\n') 
    f.write('    cosld = numpy.cos(lat)*numpy.cos(dec)\n') 
    f.write('    cosx = (numpy.cos(lha)*cosld)+sinld\n') 
    f.write('    cosx = numpy.cos(theta)\n') 
    f.write('    secx = 1.0E+0/(cosx+1.0E-30)\n') 

    f.write('    # Data taken from photolysis.txt. Calculations done in the form of:\n') 
    f.write('    # j(k) = l(k)*cosx**( mm(k))*numpy.exp(-nn(k)*secx)\n') 
    f.write('    J=[None]*62\n') 
    f.write('    #J          L           M          N\n') 
    f.write('    J[1]=6.073E-05*cosx**(1.743)*numpy.exp(-1.0*0.474*secx)\n') 
    f.write('    J[2]=4.775E-04*cosx**(0.298)*numpy.exp(-1.0*0.080*secx)\n') 
    f.write('    J[3]=1.041E-05*cosx**(0.723)*numpy.exp(-1.0*0.279*secx)\n') 
    f.write('    J[4]=1.165E-02*cosx**(0.244)*numpy.exp(-1.0*0.267*secx)\n') 
    f.write('    J[5]=2.485E-02*cosx**(0.168)*numpy.exp(-1.0*0.108*secx)\n') 
    f.write('    J[6]=1.747E-01*cosx**(0.155)*numpy.exp(-1.0*0.125*secx)\n') 
    f.write('    J[7]=2.644E-03*cosx**(0.261)*numpy.exp(-1.0*0.288*secx)\n') 
    f.write('    J[8]=9.312E-07*cosx**(1.230)*numpy.exp(-1.0*0.307*secx)\n') 
    f.write('    J[11]=4.642E-05*cosx**(0.762)*numpy.exp(-1.0*0.353*secx)\n') 
    f.write('    J[12]=6.853E-05*cosx**(0.477)*numpy.exp(-1.0*0.323*secx)\n') 
    f.write('    J[13]=7.344E-06*cosx**(1.202)*numpy.exp(-1.0*0.417*secx)\n') 
    f.write('    J[14]=2.879E-05*cosx**(1.067)*numpy.exp(-1.0*0.358*secx)\n') 
    f.write('    J[15]=2.792E-05*cosx**(0.805)*numpy.exp(-1.0*0.338*secx)\n') 
    f.write('    J[16]=1.675E-05*cosx**(0.805)*numpy.exp(-1.0*0.338*secx)\n') 
    f.write('    J[17]=7.914E-05*cosx**(0.764)*numpy.exp(-1.0*0.364*secx)\n') 
    f.write('    J[18]=1.140E-05*cosx**(0.396)*numpy.exp(-1.0*0.298*secx)\n') 
    f.write('    J[19]=1.140E-05*cosx**(0.396)*numpy.exp(-1.0*0.298*secx)\n') 
    f.write('    J[21]=7.992E-07*cosx**(1.578)*numpy.exp(-1.0*0.271*secx)\n') 
    f.write('    J[22]=5.804E-06*cosx**(1.092)*numpy.exp(-1.0*0.377*secx)\n') 
    f.write('    J[23]=1.836E-05*cosx**(0.395)*numpy.exp(-1.0*0.296*secx)\n') 
    f.write('    J[24]=1.836E-05*cosx**(0.395)*numpy.exp(-1.0*0.296*secx)\n') 
    f.write('    J[31]=6.845E-05*cosx**(0.130)*numpy.exp(-1.0*0.201*secx)\n') 
    f.write('    J[32]=1.032E-05*cosx**(0.130)*numpy.exp(-1.0*0.201*secx)\n') 
    f.write('    J[33]=3.802E-05*cosx**(0.644)*numpy.exp(-1.0*0.312*secx)\n') 
    f.write('    J[34]=1.537E-04*cosx**(0.170)*numpy.exp(-1.0*0.208*secx)\n') 
    f.write('    J[35]=3.326E-04*cosx**(0.148)*numpy.exp(-1.0*0.215*secx)\n') 
    f.write('    J[41]=7.649E-06*cosx**(0.682)*numpy.exp(-1.0*0.279*secx)\n') 
    f.write('    J[51]=1.588E-06*cosx**(1.154)*numpy.exp(-1.0*0.318*secx)\n') 
    f.write('    J[52]=1.907E-06*cosx**(1.244)*numpy.exp(-1.0*0.335*secx)\n') 
    f.write('    J[53]=2.485E-06*cosx**(1.196)*numpy.exp(-1.0*0.328*secx)\n') 
    f.write('    J[54]=4.095E-06*cosx**(1.111)*numpy.exp(-1.0*0.316*secx)\n') 
    f.write('    J[55]=1.135E-05*cosx**(0.974)*numpy.exp(-1.0*0.309*secx)\n') 
    f.write('    J[56]=7.549E-06*cosx**(1.015)*numpy.exp(-1.0*0.324*secx)\n') 
    f.write('    J[57]=3.363E-06*cosx**(1.296)*numpy.exp(-1.0*0.322*secx)\n') 
    f.write('    J[61]=7.537E-04*cosx**(0.499)*numpy.exp(-1.0*0.266*secx)\n') 
    f.write('    TEMP=temp\n') 
    #for key in mcm_constants_dict.keys():
    #    f.write('    %s =mcm_constants_dict[%s] \n' %(key,repr(str(key)))) 
    #f.write('\n')     
    #    # constants used in calculation of reaction rates
    #f.write('    M  = 2.55E+19\n')     #Check this against pressure assumed in model!
    #f.write('    N2 = 0.79*M\n')   
    #f.write('    O2 = 0.2095*M\n')   
    f.write('\n') 
    f.write('    # Creating numpy array to hold results of expressions taken from .eqn file \n') 
    f.write('    rate_values=numpy.zeros(%s)\n' %(len(rate_dict.keys()))) 
    # Now cycle through the rate_dict dictionary and print to file
    for key in rate_dict.keys():
        f.write('    rate_values[%s] =%s \n' %(key,rate_dict[key]))        
    f.write('\n') 
    f.write('    return rate_values \n')   
    f.close()  
   
def write_rate_file_numba(filename,rate_dict):
    
    #Put all of the rate coefficient functional forms into a new python file.
    f = open('Rate_coefficients_numba.py','w')
    f.write('##################################################################################################### \n') # python will convert \n to os.linesep
    f.write('# Python function to hold expressions for calculating rate coefficients for a given equation number # \n') # python will convert \n to os.linesep
    f.write('#    Copyright (C) 2017  David Topping : david.topping@manchester.ac.uk                             # \n')
    f.write('#                                      : davetopp80@gmail.com                                       # \n')
    f.write('#    Personal website: davetoppingsci.com                                                           # \n')
    f.write('#                                                                                                   # \n')
    f.write('#    This program does not have a license, meaning the deault copyright law applies.                # \n')
    f.write('#    Only users who have access to the private repository that holds this file may                  # \n')
    f.write('#    use it or develop it, but may not distribute it.                                               # \n')
    f.write('#                                                                                                   # \n')
    f.write('#                                                                                                   # \n')
    f.write('##################################################################################################### \n')    
    f.write('# File created at %s \n' % datetime.datetime.now()) # python will convert \n to os.linesep
    f.write('\n') 
    f.write('import numpy as np \n') 
    f.write('import numba as nb')
    f.write('\n') 
    f.write('@nb.jit(nb.f8(nb.f8), nopython=True) \n') 
    f.write('def numba_exp(x): \n') 
    f.write('    return np.exp(x) \n') 
    f.write('\n') 
    f.write('@nb.jit(nb.f8(nb.f8), nopython=True) \n') 
    f.write('def numba_log10(x): \n') 
    f.write('    return np.log10(x) \n') 
    f.write('\n') 
    f.write('@nb.jit(nb.f8(nb.f8), nopython=True) \n') 
    f.write('def numba_log(x): \n') 
    f.write('    return np.log(x) \n') 
    f.write('\n') 
    f.write('@nb.jit(nb.f8(nb.f8), nopython=True) \n') 
    f.write('def numba_abs(x): \n') 
    f.write('    return np.abs(x) \n') 
    f.write('\n') 
    f.write('@nb.jit(nb.f8(nb.f8), nopython=True) \n') 
    f.write('def numba_sqrt(x): \n') 
    f.write('    return np.sqrt(x) \n') 
    f.write('\n') 
    f.write('@nb.jit(nb.f8(nb.f8), nopython=True) \n') 
    f.write('def numba_arccos(x): \n') 
    f.write('    return np.arccos(x) \n') 
    f.write('\n') 
    f.write('@nb.jit(nb.f8(nb.f8), nopython=True) \n') 
    f.write('def numba_cos(x): \n') 
    f.write('    return np.cos(x) \n') 
    f.write('\n') 
    f.write('@nb.jit(nb.f8(nb.f8), nopython=True) \n') 
    f.write('def numba_sin(x): \n') 
    f.write('    return np.sin(x) \n') 
    f.write('\n') 
    f.write('@nb.jit(nb.f8(nb.f8), nopython=True) \n') 
    f.write('def numba_arctan(x): \n') 
    f.write('    return np.arctan(x) \n') 
    f.write('\n') 
    f.write('@nb.jit(\'float64[:](float64,float64,float64,float64,float64[:],float64[:])\', nopython=True)\n')
    f.write('def evaluate_rates(ttime,RO2,H2O,temp,rate_values,J):\n') 
    # Now cycle through all of the mcm_constants_dict values.
    # Please note this includes photolysis rates that will change with time of day or condition in the chamber. You will
    # need to modify this potentially.
    f.write('    # Creating reference to constant values used in rate expressions\n') 

    f.write('    # constants used in calculation of reaction rates\n') 
    f.write('    M  = 2.55E+19  #Check this against pressure assumed in model\n') 
    f.write('    N2 = 0.79*M\n') 
    f.write('    O2 = 0.2095*M\n') 

    f.write('    # kro2no : ro2      + no      = ro      + no2\n') 
    f.write('    #        : ro2      + no      = rono2\n') 
    f.write('    # iupac 1992\n') 
    f.write('    KRONO2    = 2.40E-12*numba_exp(360.0/temp)\n') 

    f.write('    # kro2ho2: ro2      + ho2     = rooh    + o2\n') 
    f.write('    # mcm protocol [1997]\n') 
    f.write('    KRO2HO2   = 2.91E-13*numba_exp(1300.0/temp)\n') 

    f.write('    # kapho2 : rcoo2    + ho2     = products\n') 
    f.write('    # mcm protocol [1997]\n') 
    f.write('    KAPHO2    = 4.30E-13*numba_exp(1040.0/temp)\n') 

    f.write('    # kapno  : rcoo2    + no      = products\n') 
    f.write('    # mej [1998]\n') 
    f.write('    KAPNO = 8.10E-12*numba_exp(270.0/temp)\n') 

    f.write('    # kro2no3: ro2      + no3     = products\n') 
    f.write('    # mcm protocol [1997]\n') 
    f.write('    KRO2NO3   = 2.50E-12\n') 

    f.write('    # kno3al : no3      + rcho    = rcoo2   + hno3\n') 
    f.write('    # mcm protocol [1997]\n') 
    f.write('    KNO3AL    = 1.44E-12*numba_exp(-1862.0/temp)\n') 

    f.write('    # kdec   : ro                 = products\n') 
    f.write('    # mcm protocol [1997]\n') 
    f.write('    KDEC      = 1.00E+06\n') 

    f.write('    KROSEC = 1.50E-14*numba_exp(-200.0/temp)\n') 

    f.write('    KALKOXY=3.70E-14*numba_exp(-460.0/temp)*O2\n') 

    f.write('    KALKPXY=1.80E-14*numba_exp(-260.0/temp)*O2\n') 

    f.write('    # -------------------------------------------------------------------\n') 
    f.write('    # complex reactions\n') 
    f.write('    # -------------------------------------------------------------------\n') 

    f.write('    # kfpan kbpan\n') 
    f.write('    # formation and decomposition of pan\n') 
    f.write('    # iupac 2001 (mcmv3.2)\n') 
    f.write('    kc0     = 2.70E-28*M*(temp/300.0)**(-7.1)\n') 
    f.write('    kci     = 1.21E-11*(temp/300.0)**(-0.9)\n') 
    f.write('    krc     = kc0/kci\n') 
    f.write('    fcc     = 0.30\n') 
    f.write('    nc      = 0.75-(1.27*numba_log10(fcc))\n') 
    f.write('    fc      = 10**(numba_log10(fcc)/(1.0+((numba_log10(krc))/nc)**2.0))\n') 
    f.write('    KFPAN   = (kc0*kci)*fc/(kc0+kci)\n') 

    f.write('    kd0     = 4.90E-03*M*numba_exp(-12100.0/temp)\n') 
    f.write('    kdi     = 5.40E+16*numba_exp(-13830.0/temp)\n') 
    f.write('    krd     = kd0/kdi\n') 
    f.write('    fcd     = 0.30\n') 
    f.write('    ncd     = 0.75-(1.27*numba_log10(fcd))\n') 
    f.write('    fd      = 10.0**(numba_log10(fcd)/(1.0+((numba_log10(krd))/ncd)**2.0))\n') 
    f.write('    KBPAN   = (kd0*kdi)*fd/(kd0+kdi)\n') 

    f.write('    # kmt01  : o        + no      = no2\n') 
    f.write('    # iupac 2001 (mcmv3.2)\n') 
    f.write('    k10     = 1.00E-31*M*(temp/300.0)**(-1.6)\n') 

    f.write('    k1i     = 3.00E-11*(temp/300.0)**(0.3)\n') 
    f.write('    kr1     = k10/k1i\n') 
    f.write('    fc1     = 0.85\n') 
    f.write('    nc1     = 0.75-(1.27*numba_log10(fc1))\n') 
    f.write('    f1      = 10.0**(numba_log10(fc1)/(1.0+((numba_log10(kr1)/nc1))**2.0))\n') 
    f.write('    KMT01   = (k10*k1i)*f1/(k10+k1i)\n') 

    f.write('    # kmt02  : o        + no2     = no3\n') 
    f.write('    # iupac 2001 (mcmv3.2)\n') 
    f.write('    k20     = 1.30E-31*M*(temp/300.0)**(-1.5)\n') 
    f.write('    k2i     = 2.30E-11*(temp/300.0)**(0.24)\n') 
    f.write('    kr2     = k20/k2i\n') 
    f.write('    fc2     = 0.6\n') 
    f.write('    nc2     = 0.75-(1.27*numba_log10(fc2))\n') 
    f.write('    f2      = 10.0**(numba_log10(fc2)/(1.0+((numba_log10(kr2)/nc2))**2.0))\n') 
    f.write('    KMT02   = (k20*k2i)*f2/(k20+k2i)\n') 

    f.write('    # kmt03  : no2      + no3     = n2o5\n') 
    f.write('    # iupac 2006, mcmv3.2\n') 
    f.write('    k30     = 3.60E-30*M*(temp/300.0)**(-4.1)\n') 
    f.write('    k3i     = 1.90E-12*(temp/300.0)**(0.2)\n') 
    f.write('    kr3     = k30/k3i\n') 
    f.write('    fc3     = 0.35\n') 
    f.write('    nc3     = 0.75-(1.27*numba_log10(fc3))\n') 
    f.write('    f3      = 10.0**(numba_log10(fc3)/(1.0+((numba_log10(kr3)/nc3))**2.0))\n') 
    f.write('    KMT03   = (k30*k3i)*f3/(k30+k3i)\n') 

    f.write('    # kmt04  : n2o5               = no2     + no3\n') 
    f.write('    # iupac 2006, mcmv3.2\n') 
    f.write('    k40     = 1.30E-03*M*(temp/300.0)**(-3.5)*numba_exp(-11000.0/temp)\n') 
    f.write('    k4i     = 9.70E+14*(temp/300.0)**(0.1)*numba_exp(-11080.0/temp)\n') 
    f.write('    kr4     = k40/k4i\n') 
    f.write('    fc4     = 0.35\n') 
    f.write('    nc4     = 0.75-(1.27*numba_log10(fc4))\n') 
    f.write('    f4      = 10.0**(numba_log10(fc4)/(1+((numba_log10(kr4)/nc4))**2.0))\n') 
    f.write('    KMT04   = (k40*k4i)*f4/(k40+k4i)\n') 

    f.write('    # kmt05  : oh       + co(+o2) = ho2     + co2\n') 
    f.write('    # iupac 2006\n') 
    f.write('    KMT05  = (1.0 + (M/4.2E19))\n') 

    f.write('    # kmt06  : ho2      + ho2     = h2o2    + o2\n') 
    f.write('    # water enhancement factor\n') 
    f.write('    # iupac 1992\n') 
    f.write('    KMT06  = 1.0 + (1.40E-21*numba_exp(2200.0/temp)*H2O)\n') 

    f.write('    # kmt06  = 1.0 + (2.00E-25*numba_exp(4670.0/temp)*h2o)\n') 
    f.write('    # S+R 2005 values\n') 

    f.write('    # kmt07  : oh       + no      = hono\n') 

    f.write('    # iupac 2006, mcmv3.2\n') 
    f.write('    k70     = 7.40E-31*M*(temp/300.0)**(-2.4)\n') 
    f.write('    k7i     = 3.30E-11*(temp/300.0)**(-0.3)\n') 
    f.write('    kr7     = k70/k7i\n') 
    f.write('    fc7     = numba_exp(-temp/1420.0)\n') 
    f.write('    nc7     = 0.75-(1.27*numba_log10(fc7))\n') 
    f.write('    f7      = 10.0**(numba_log10(fc7)/(1+((numba_log10(kr7)/nc7))**2.0))\n') 
    f.write('    KMT07   = (k70*k7i)*f7/(k70+k7i)\n') 

    f.write('    # kmt08  : oh       + no2     = hno3\n') 

    f.write('    # iupac 2006, mcmv3.2\n') 
    f.write('    k80     = 3.30E-30*M*(temp/300.0)**(-3.0)\n') 
    f.write('    k8i     = 4.10E-11\n') 
    f.write('    kr8     = k80/k8i\n') 
    f.write('    fc8     = 0.4\n') 
    f.write('    nc8     = 0.75-(1.27*numba_log10(fc8))\n') 
    f.write('    f8      = 10.0**(numba_log10(fc8)/(1.0+((numba_log10(kr8)/nc8))**2.0))\n') 
    f.write('    KMT08   = (k80*k8i)*f8/(k80+k8i)\n') 

    f.write('    # kmt09  : ho2      + no2     = ho2no2\n') 
    f.write('    # iupac 1997, mcmv3.2\n') 
    
    f.write('    k90     = 1.80E-31*M*(temp/300.0)**(-3.2)\n') 
    f.write('    k9i     = 4.70E-12\n') 
    f.write('    kr9     = k90/k9i\n') 
    f.write('    fc9     = 0.6\n') 
    f.write('    nc9     = 0.75-(1.27*numba_log10(fc9))\n') 
    f.write('    f9      = 10.0**(numba_log10(fc9)/(1.0+((numba_log10(kr9)/nc9))**2.0))\n') 
    f.write('    KMT09   = (k90*k9i)*f9/(k90+k9i)\n') 

    f.write('    # kmt10  : ho2no2             = ho2     + no2\n') 
    f.write('    # iupac 1997, mcmv3.2\n') 

    f.write('    k100     = 4.10E-05*M*numba_exp(-10650.0/temp)\n') 
    f.write('    k10i     = 4.80E+15*numba_exp(-11170.0/temp)\n') 
    f.write('    kr10     = k100/k10i\n') 
    f.write('    fc10     = 0.6\n') 
    f.write('    nc10     = 0.75-(1.27*numba_log10(fc10))\n') 
    f.write('    f10      = 10.0**(numba_log10(fc10)/(1.0+((numba_log10(kr10)/nc10))**2.0))\n') 
    f.write('    KMT10    = (k100*k10i)*f10/(k100+k10i)\n') 

    f.write('    # kmt11  : oh       + hno3    = h2o     + no3\n') 
    f.write('    # iupac 2006, mcmv3.2\n') 

    f.write('    k1     = 2.40E-14*numba_exp(460.0/temp)\n') 
    f.write('    k3     = 6.50E-34*numba_exp(1335.0/temp)\n') 
    f.write('    k4     = 2.70E-17*numba_exp(2199.0/temp)\n') 
    f.write('    k2     = (k3*M)/(1.0+(k3*M/k4))\n') 
    f.write('    KMT11  = k1 + k2\n') 

    f.write('    # kmt12 iupac 2006, mcmv3.2\n') 

    f.write('    k120 = 4.50E-31*((temp/300.0)**(-3.9))*M\n') 
    f.write('    k12i = 1.30E-12*((temp/300.0)**(-0.7))\n') 
    f.write('    kr12 = k120/k12i\n') 
    f.write('    fc12 = 0.525\n') 
    f.write('    nc12 = 0.75-(1.27*numba_log10(fc12))\n') 
    f.write('    f12  = 10.0**(numba_log10(fc12)/(1.0+((numba_log10(kr12)/nc12))**2.0))\n') 
    f.write('    KMT12    = (k120*k12i)*f12/(k120+k12i)\n') 

    f.write('    # kmt13  : ch3o2    + no2     = ch3o2no2\n') 
    f.write('    # iupac 2006\n') 

    f.write('    k130     = 2.50E-30*((temp/300.0)**(-5.5))*M\n') 
    f.write('    k13i     = 1.80E-11\n') 
    f.write('    kr13     = k130/k13i\n') 
    f.write('    fc13     = 0.36\n') 
    f.write('    nc13     = 0.75-(1.27*numba_log10(fc13))\n') 
    f.write('    f13      = 10.0**(numba_log10(fc13)/(1.0+((numba_log10(kr13)/nc13))**2.0))\n') 
    f.write('    KMT13    = (k130*k13i)*f13/(k130+k13i)\n') 

    f.write('    # kmt14  : ch3o2no2           = ch3o2   + no2\n') 
    f.write('    # iupac 2006, mcmv3.2\n') 

    f.write('    k140     = 9.00E-05*numba_exp(-9690.0/temp)*M\n') 
    f.write('    k14i     = 1.10E+16*numba_exp(-10560.0/temp)\n') 
    f.write('    kr14     = k140/k14i\n') 
    f.write('    fc14     = 0.4\n') 
    f.write('    nc14     = 0.75-(1.27*numba_log10(fc14))\n') 
    f.write('    f14      = 10.0**(numba_log10(fc14)/(1.0+((numba_log10(kr14)/nc14))**2.0))\n') 
    f.write('    KMT14    = (k140*k14i)*f14/(k140+k14i)\n') 

    f.write('    # kmt15 iupac 2006, mcmv3.2\n') 

    f.write('    k150 = 8.60E-29*((temp/300.0)**(-3.1))*M\n') 
    f.write('    k15i = 9.00E-12*((temp/300.0)**(-0.85))\n') 
    f.write('    kr15 = k150/k15i\n') 
    f.write('    fc15 = 0.48\n') 
    f.write('    nc15 = 0.75-(1.27*numba_log10(fc15))\n') 
    f.write('    f15  = 10.0**(numba_log10(fc15)/(1.0+((numba_log10(kr15)/nc15))**2.0))\n') 
    f.write('    KMT15 = (k150*k15i)*f15/(k150+k15i)\n') 

    f.write('    # kmt16  :  oh  +  c3h6\n') 
    f.write('    # iupac 2006\n') 

    f.write('    k160     = 8.00E-27*((temp/300.0)**(-3.5))*M\n') 
    f.write('    k16i     = 3.00E-11*((temp/300.0)**(-1.0))\n') 
    f.write('    kr16     = k160/k16i\n') 
    f.write('    fc16     = 0.5\n') 
    f.write('    nc16     = 0.75-(1.27*numba_log10(fc16))\n') 
    f.write('    f16      = 10.0**(numba_log10(fc16)/(1.0+((numba_log10(kr16)/nc16))**2.0))\n') 
    f.write('    KMT16    = (k160*k16i)*f16/(k160+k16i)\n') 

    f.write('    # kmt17 iupac 2006\n') 

    f.write('    k170 = 5.00E-30*((temp/300.0)**(-1.5))*M\n') 
    f.write('    k17i = 1.00E-12\n') 
    f.write('    kr17 = k170/k17i\n') 
    f.write('    fc17 = (0.17*numba_exp(-51./temp))+numba_exp(-1.0*temp/204.)\n') 
    f.write('    nc17 = 0.75-(1.27*numba_log10(fc17))\n') 
    f.write('    f17  = 10.0**(numba_log10(fc17)/(1.0+((numba_log10(kr17)/nc17))**2.0))\n') 
    f.write('    KMT17 = (k170*k17i)*f17/(k170+k17i)\n') 

    f.write('    #Taken from yet another MCM version, will need converting to lower case\n') 
    f.write('    KRO2NO = 2.7E-12*numba_exp(360/temp)\n') 

    f.write('    KRO2HO2 = 2.91E-13*numba_exp(1300/temp)\n') 

    f.write('    KAPHO2 = 5.2E-13*numba_exp(980/temp)\n') 

    f.write('    KAPNO = 7.5E-12*numba_exp(290/temp)\n') 

    f.write('    KRO2NO3 = 2.3E-12\n') 

    f.write('    KNO3AL = 1.4E-12*numba_exp(-1860/temp)\n') 

    f.write('    KDEC = 1.00E+06\n') 

    f.write('    KROPRIM = 2.50E-14*numba_exp(-300/temp)\n') 

    f.write('    KROSEC = 2.50E-14*numba_exp(-300/temp)\n') 

    f.write('    KCH3O2 = 1.03E-13*numba_exp(365/temp)\n') 

    f.write('    K298CH3O2 = 3.5E-13\n') 

    f.write('    KMT18 = 9.5E-39*O2*numba_exp(5270/temp)/(1+7.5E-29*O2*numba_exp(5610/temp))\n') 

    f.write('    KROPRIM  = 3.70E-14*numba_exp(-460.0/temp)\n') 

    f.write('    KROSEC   = 1.80E-14*numba_exp(-260.0/temp)\n') 

    f.write('    # ************************************************************************\n') 
    f.write('    # define photolysis reaction rates using derwent method from mcm2box.fac\n') 
    f.write('    # ************************************************************************\n') 

    f.write('    # solar declination angle \n') 
    f.write('    dec = 23.79\n') 
    f.write('    # latitude\n') 
    f.write('    lat = 50.0\n') 
    f.write('    pi = 4.0*numba_arctan(1.0)\n') 
    f.write('    # local hour angle - representing time of day\n') 
    f.write('    lha = (1.0+ttime/4.32E+4)*pi\n') 
    f.write('    radian = 180.0/pi\n') 
    f.write('    lat = lat/radian\n') 
    f.write('    dec = dec/radian\n') 
    f.write('    theta = numba_arccos(numba_cos(lha)*numba_cos(dec)*numba_cos(lat)+numba_sin(dec)*numba_sin(lat))\n') 
    f.write('    sinld = numba_sin(lat)*numba_sin(dec)\n') 
    f.write('    cosld = numba_cos(lat)*numba_cos(dec)\n') 
    f.write('    cosx = (numba_cos(lha)*cosld)+sinld\n') 
    f.write('    cosx = numba_cos(theta)\n') 
    f.write('    secx = 1.0E+0/(cosx+1.0E-30)\n') 

    f.write('    # Data taken from photolysis.txt. Calculations done in the form of:\n') 
    f.write('    # j(k) = l(k)*cosx**( mm(k))*numba_exp(-nn(k)*secx)\n') 
    #f.write('    J=[None]*62\n') 
    f.write('    #J          L           M          N\n') 
    f.write('    J[1]=6.073E-05*cosx**(1.743)*numba_exp(-1.0*0.474*secx)\n') 
    f.write('    J[2]=4.775E-04*cosx**(0.298)*numba_exp(-1.0*0.080*secx)\n') 
    f.write('    J[3]=1.041E-05*cosx**(0.723)*numba_exp(-1.0*0.279*secx)\n') 
    f.write('    J[4]=1.165E-02*cosx**(0.244)*numba_exp(-1.0*0.267*secx)\n') 
    f.write('    J[5]=2.485E-02*cosx**(0.168)*numba_exp(-1.0*0.108*secx)\n') 
    f.write('    J[6]=1.747E-01*cosx**(0.155)*numba_exp(-1.0*0.125*secx)\n') 
    f.write('    J[7]=2.644E-03*cosx**(0.261)*numba_exp(-1.0*0.288*secx)\n') 
    f.write('    J[8]=9.312E-07*cosx**(1.230)*numba_exp(-1.0*0.307*secx)\n') 
    f.write('    J[11]=4.642E-05*cosx**(0.762)*numba_exp(-1.0*0.353*secx)\n') 
    f.write('    J[12]=6.853E-05*cosx**(0.477)*numba_exp(-1.0*0.323*secx)\n') 
    f.write('    J[13]=7.344E-06*cosx**(1.202)*numba_exp(-1.0*0.417*secx)\n') 
    f.write('    J[14]=2.879E-05*cosx**(1.067)*numba_exp(-1.0*0.358*secx)\n') 
    f.write('    J[15]=2.792E-05*cosx**(0.805)*numba_exp(-1.0*0.338*secx)\n') 
    f.write('    J[16]=1.675E-05*cosx**(0.805)*numba_exp(-1.0*0.338*secx)\n') 
    f.write('    J[17]=7.914E-05*cosx**(0.764)*numba_exp(-1.0*0.364*secx)\n') 
    f.write('    J[18]=1.140E-05*cosx**(0.396)*numba_exp(-1.0*0.298*secx)\n') 
    f.write('    J[19]=1.140E-05*cosx**(0.396)*numba_exp(-1.0*0.298*secx)\n') 
    f.write('    J[21]=7.992E-07*cosx**(1.578)*numba_exp(-1.0*0.271*secx)\n') 
    f.write('    J[22]=5.804E-06*cosx**(1.092)*numba_exp(-1.0*0.377*secx)\n') 
    f.write('    J[23]=1.836E-05*cosx**(0.395)*numba_exp(-1.0*0.296*secx)\n') 
    f.write('    J[24]=1.836E-05*cosx**(0.395)*numba_exp(-1.0*0.296*secx)\n') 
    f.write('    J[31]=6.845E-05*cosx**(0.130)*numba_exp(-1.0*0.201*secx)\n') 
    f.write('    J[32]=1.032E-05*cosx**(0.130)*numba_exp(-1.0*0.201*secx)\n') 
    f.write('    J[33]=3.802E-05*cosx**(0.644)*numba_exp(-1.0*0.312*secx)\n') 
    f.write('    J[34]=1.537E-04*cosx**(0.170)*numba_exp(-1.0*0.208*secx)\n') 
    f.write('    J[35]=3.326E-04*cosx**(0.148)*numba_exp(-1.0*0.215*secx)\n') 
    f.write('    J[41]=7.649E-06*cosx**(0.682)*numba_exp(-1.0*0.279*secx)\n') 
    f.write('    J[51]=1.588E-06*cosx**(1.154)*numba_exp(-1.0*0.318*secx)\n') 
    f.write('    J[52]=1.907E-06*cosx**(1.244)*numba_exp(-1.0*0.335*secx)\n') 
    f.write('    J[53]=2.485E-06*cosx**(1.196)*numba_exp(-1.0*0.328*secx)\n') 
    f.write('    J[54]=4.095E-06*cosx**(1.111)*numba_exp(-1.0*0.316*secx)\n') 
    f.write('    J[55]=1.135E-05*cosx**(0.974)*numba_exp(-1.0*0.309*secx)\n') 
    f.write('    J[56]=7.549E-06*cosx**(1.015)*numba_exp(-1.0*0.324*secx)\n') 
    f.write('    J[57]=3.363E-06*cosx**(1.296)*numba_exp(-1.0*0.322*secx)\n') 
    f.write('    J[61]=7.537E-04*cosx**(0.499)*numba_exp(-1.0*0.266*secx)\n') 
    f.write('    TEMP=temp\n') 
    #for key in mcm_constants_dict.keys():
    #    f.write('    %s =mcm_constants_dict[%s] \n' %(key,repr(str(key)))) 
    #f.write('\n')     
    #    # constants used in calculation of reaction rates
    #f.write('    M  = 2.55E+19\n')     #Check this against pressure assumed in model!
    #f.write('    N2 = 0.79*M\n')   
    #f.write('    O2 = 0.2095*M\n')   
    f.write('\n') 
    f.write('    # Creating numpy array to hold results of expressions taken from .eqn file \n') 
    #f.write('    rate_values=numpy.zeros(%s)\n' %(len(rate_dict.keys()))) 
    # Now cycle through the rate_dict dictionary and print to file
    for key in rate_dict.keys():
        f.write('    rate_values[%s] =%s \n' %(key,rate_dict[key]))        
    f.write('\n') 
    f.write('    return rate_values \n')   
    f.close()  

def write_rate_file_fortran(filename,rate_dict_fortran,openMP):
    
    #Put all of the rate coefficient functional forms into a new Fortran file.
    f = open('Rate_coefficients.f90','w')
    #f.write('¡ Fortran function to hold expressions for calculating rate coefficients for a given equation number  \n') # python will convert \n to os.linesep
    #f.write('!    Copyright (C) 2018  David Topping : david.topping@manchester.ac.uk                              \n')
    #f.write('!                                      : davetopp80@gmail.com                                        \n')
    #f.write('!    Personal website: davetoppingsci.com                                                            \n')
    #f.write('!                                                                                                    \n')
    #f.write('!    This program does not have a license, meaning the deault copyright law applies.                 \n')
    #f.write('!    Only users who have access to the private repository that holds this file may                   \n')
    #f.write('!    use it or develop it, but may not distribute it.                                                \n')
    #f.write('!                                                                                                    \n')
    #f.write('!                                                                                                    \n')
    #f.write('!  \n')    
    #f.write('! File created at %s \n' % datetime.datetime.now()) # python will convert \n to os.linesep
    #f.write('\n') 
    #f.write('import numpy \n') 
    #f.write('\n') 
    f.write('subroutine evaluate_rates(RO2,H2O,TEMP,ttime,rate_values)\n') 
    f.write('    implicit none \n') 
    f.write('    REAL(8), intent(in) :: RO2,H2O,TEMP \n') 
    f.write('    !f2py intent(in) :: RO2,H2O,TEMP \n') 
    f.write('    REAL(8), intent(out), dimension(%s) :: rate_values \n' %(len(rate_dict_fortran.keys()))) 
    f.write('    !f2py intent(out) :: rate_values \n') 
    # Now cycle through all of the mcm_constants_dict values.
    # Please note this includes photolysis rates that will change with time of day or condition in the chamber. You will
    # need to modify this potentially.
    f.write('    ! Creating reference to constant values used in rate expressions\n') 
    
    # Previously looped through fixed values, but need to calculate for varying temp
    # for key in mcm_constants_dict.keys():
    #    #pdb.set_trace()
    #    if key != 'J':
    #        f.write('    REAL(8), PARAMETER :: %s = %s \n' %(key,mcm_constants_dict[key]))    

    f.write('    REAL(8) KRONO2, KRO2HO2, KAPHO2, KAPNO, KRO2NO3  \n ')
    f.write('    REAL(8) KNO3AL, KDEC, KROSEC, KALKOXY, KALKPXY  \n ')
    f.write('    REAL(8) kc0, kci, krc, nc, fc, fcc, KFPAN  \n ')
    f.write('    REAL(8) kd0, kdi, krd, fcd, ncd, fd, KBPAN  \n ')
    f.write('    REAL(8) k10, k1i, kr1, fc1, nc1, f1, KMT01  \n ')
    f.write('    REAL(8) k20, k2i, kr2, fc2, nc2, f2, KMT02  \n ')
    f.write('    REAL(8) k30, k3i, kr3, fc3, nc3, f3, KMT03  \n ')
    f.write('    REAL(8) k40, k4i, kr4, fc4, nc4, f4, KMT04  \n ')
    f.write('    REAL(8) KMT05, KMT06  \n ')
    f.write('    REAL(8) k70, k7i, kr7, fc7, nc7, f7, KMT07  \n ')
    f.write('    REAL(8) k80, k8i, kr8, fc8, nc8, f8, KMT08  \n ')
    f.write('    REAL(8) k90, k9i, kr9, fc9, nc9, f9, KMT09  \n ')
    f.write('    REAL(8) k100, k10i, kr10, fc10, nc10, f10, KMT10  \n ')
    f.write('    REAL(8) k1, k3, k4, k2, KMT11  \n ')
    f.write('    REAL(8) k120, k12i, kr12, fc12, nc12, f12, KMT12  \n ')
    f.write('    REAL(8) k130, k13i, kr13, fc13, nc13, f13, KMT13  \n ')
    f.write('    REAL(8) k140, k14i, kr14, fc14, nc14, f14, KMT14  \n ')
    f.write('    REAL(8) k150, k15i, kr15, fc15, nc15, f15, KMT15  \n ')
    f.write('    REAL(8) k160, k16i, kr16, fc16, nc16, f16, KMT16  \n ')
    f.write('    REAL(8) k170, k17i, kr17, fc17, nc17, f17, KMT17  \n ')
    f.write('    REAL(8) KRO2NO  \n ')
    f.write('    REAL(8) KROPRIM, KCH3O2, K298CH3O2, KMT18  \n ')
    f.write('\n') 
        
    # Now add an account for the variables used in changing photolysis rates
    
    f.write('    REAL(8), intent(IN)  :: ttime \n') 
    f.write('    !f2py intent(in) :: ttime \n') 
    
    f.write('    REAL(8) theta, secx, cosx \n') 
    f.write('    REAL(8) lat, pi, radian, dec, lha, sinld, cosld \n') 
     
    # Now create the photolysis varible array 'J'    
    f.write('    REAL(8), dimension(61) :: J \n') 
    
        # constants used in calculation of reaction rates
    f.write('    REAL(8), PARAMETER :: M  = 2.55E+19\n')     #Check this against pressure assumed in model!
    f.write('    REAL(8), PARAMETER :: N2 = %E \n' %(0.79*2.55E+19))   
    f.write('    REAL(8), PARAMETER :: O2 = %E \n' %(0.2095*2.55E+19))   
    f.write('\n') 
    f.write('    ! Creating numpy array to hold results of expressions taken from .eqn file \n') 
    
    # Populate values in Photolysis rates
    
    f.write('    ! solar declination angle from july 1st - harwell traj model \n') 
    f.write('    dec = 23.79\n') 
    f.write('    ! latitude\n') 
    f.write('    lat = 50.0\n') 
    f.write('    pi = 4.0*ATAN(1.0)\n') 
    f.write('    ! local hour angle - representing time of day\n') 
    f.write('    lha = (1.0+ttime/4.32d+4)*pi\n') 
    f.write('    radian = 180.0/pi\n') 
    f.write('    lat = lat/radian\n') 
    f.write('    dec = dec/radian\n') 
    f.write('    theta = ACOS(COS(lha)*COS(dec)*COS(lat)+SIN(dec)*SIN(lat))\n') 
    f.write('    sinld = SIN(lat)*SIN(dec)\n') 
    f.write('    cosld = COS(lat)*COS(dec)\n') 
    f.write('    cosx = (COS(lha)*cosld)+sinld\n') 
    f.write('    cosx = COS(theta)\n') 
    f.write('    secx = 1.0d+0/(cosx+1.0d-30)\n') 
    f.write('    J(1)=6.073E-05*cosx**(1.743)*EXP(-1.0*0.474*secx)\n') 
    f.write('    J(2)=4.775E-04*cosx**(0.298)*EXP(-1.0*0.080*secx)\n') 
    f.write('    J(3)=1.041E-05*cosx**(0.723)*EXP(-1.0*0.279*secx)\n') 
    f.write('    J(4)=1.165E-02*cosx**(0.244)*EXP(-1.0*0.267*secx)\n') 
    f.write('    J(5)=2.485E-02*cosx**(0.168)*EXP(-1.0*0.108*secx)\n') 
    f.write('    J(6)=1.747E-01*cosx**(0.155)*EXP(-1.0*0.125*secx)\n') 
    f.write('    J(7)=2.644E-03*cosx**(0.261)*EXP(-1.0*0.288*secx)\n') 
    f.write('    J(8)=9.312E-07*cosx**(1.230)*EXP(-1.0*0.307*secx)\n') 
    f.write('    J(11)=4.642E-05*cosx**(0.762)*EXP(-1.0*0.353*secx)\n') 
    f.write('    J(12)=6.853E-05*cosx**(0.477)*EXP(-1.0*0.323*secx)\n') 
    f.write('    J(13)=7.344E-06*cosx**(1.202)*EXP(-1.0*0.417*secx)\n') 
    f.write('    J(14)=2.879E-05*cosx**(1.067)*EXP(-1.0*0.358*secx)\n') 
    f.write('    J(15)=2.792E-05*cosx**(0.805)*EXP(-1.0*0.338*secx)\n') 
    f.write('    J(16)=1.675E-05*cosx**(0.805)*EXP(-1.0*0.338*secx)\n') 
    f.write('    J(17)=7.914E-05*cosx**(0.764)*EXP(-1.0*0.364*secx)\n') 
    f.write('    J(18)=1.140E-05*cosx**(0.396)*EXP(-1.0*0.298*secx)\n') 
    f.write('    J(19)=1.140E-05*cosx**(0.396)*EXP(-1.0*0.298*secx)\n') 
    f.write('    J(21)=7.992E-07*cosx**(1.578)*EXP(-1.0*0.271*secx)\n') 
    f.write('    J(22)=5.804E-06*cosx**(1.092)*EXP(-1.0*0.377*secx)\n') 
    f.write('    J(23)=1.836E-05*cosx**(0.395)*EXP(-1.0*0.296*secx)\n') 
    f.write('    J(24)=1.836E-05*cosx**(0.395)*EXP(-1.0*0.296*secx)\n') 
    f.write('    J(31)=6.845E-05*cosx**(0.130)*EXP(-1.0*0.201*secx)\n') 
    f.write('    J(32)=1.032E-05*cosx**(0.130)*EXP(-1.0*0.201*secx)\n') 
    f.write('    J(33)=3.802E-05*cosx**(0.644)*EXP(-1.0*0.312*secx)\n') 
    f.write('    J(34)=1.537E-04*cosx**(0.170)*EXP(-1.0*0.208*secx)\n') 
    f.write('    J(35)=3.326E-04*cosx**(0.148)*EXP(-1.0*0.215*secx)\n') 
    f.write('    J(41)=7.649E-06*cosx**(0.682)*EXP(-1.0*0.279*secx)\n') 
    f.write('    J(51)=1.588E-06*cosx**(1.154)*EXP(-1.0*0.318*secx)\n') 
    f.write('    J(52)=1.907E-06*cosx**(1.244)*EXP(-1.0*0.335*secx)\n') 
    f.write('    J(53)=2.485E-06*cosx**(1.196)*EXP(-1.0*0.328*secx)\n') 
    f.write('    J(54)=4.095E-06*cosx**(1.111)*EXP(-1.0*0.316*secx)\n') 
    f.write('    J(55)=1.135E-05*cosx**(0.974)*EXP(-1.0*0.309*secx)\n') 
    f.write('    J(56)=7.549E-06*cosx**(1.015)*EXP(-1.0*0.324*secx)\n') 
    f.write('    J(57)=3.363E-06*cosx**(1.296)*EXP(-1.0*0.322*secx)\n') 
    f.write('    J(61)=7.537E-04*cosx**(0.499)*EXP(-1.0*0.266*secx)\n') 

    f.write('    ! Calculating standard rate coefficients \n') 
    f.write('    ! iupac 1992 \n')
    f.write('    KRONO2    = 2.40E-12*EXP(360.0/TEMP) \n')
    #mcm_constants_dict['KRONO2']=kro2no
    f.write('    ! kro2ho2: ro2      + ho2     = rooh    + o2 \n')
    f.write('    ! mcm protocol [1997] \n')
    f.write('    KRO2HO2   = 2.91E-13*EXP(1300.0/TEMP) \n')
    #mcm_constants_dict['KRO2HO2']=kro2ho2
    f.write('    ! kapho2 : rcoo2    + ho2     = products \n')
    f.write('    ! mcm protocol [1997] \n')
    f.write('    KAPHO2    = 4.30E-13*EXP(1040.0/TEMP) \n')
    #mcm_constants_dict['KAPHO2']=kapho2 
    f.write('    ! kapno  : rcoo2    + no      = products \n')
    f.write('    ! mej [1998] \n')
    f.write('    KAPNO = 8.10E-12*EXP(270.0/TEMP) \n')
    #mcm_constants_dict['KAPNO']=kapno 
    f.write('    ! kro2no3: ro2      + no3     = products \n')
    f.write('    ! mcm protocol [1997] \n')
    f.write('    KRO2NO3   = 2.50E-12 \n')
    #mcm_constants_dict['KRO2NO3']=kro2no3
    f.write('    ! kno3al : no3      + rcho    = rcoo2   + hno3 \n')
    f.write('    ! mcm protocol [1997] \n')
    f.write('    KNO3AL    = 1.44E-12*EXP(-1862.0/TEMP) \n')
    #mcm_constants_dict['KNO3AL']=kno3al 
    f.write('    ! kdec   : ro                 = products \n')
    f.write('    ! mcm protocol [1997] \n')
    f.write('    KDEC      = 1.00E06 \n')
    #mcm_constants_dict['KDEC']=kdec
    f.write('    KROSEC = 1.50E-14*EXP(-200.0/TEMP) \n')
    #mcm_constants_dict['KROSEC']=krosec
    f.write('    KALKOXY=3.70E-14*EXP(-460.0/TEMP)*O2 \n')
    #mcm_constants_dict['KALKOXY']=kalkoxy
    f.write('    KALKPXY=1.80E-14*EXP(-260.0/TEMP)*O2 \n')
    #mcm_constants_dict['KALKPXY']=kalkpxy
    f.write('    ! ------------------------------------------------------------------- \n')
    f.write('    ! complex reactions \n')
    f.write('    ! ------------------------------------------------------------------- \n')
    f.write('    ! kfpan kbpan \n')
    f.write('    ! formation and decomposition of pan \n')
    f.write('    ! iupac 2001 (mcmv3.2) \n')
    f.write('    kc0     = 2.70E-28*M*(TEMP/300.0)**(-7.1) \n') 
    f.write('    kci     = 1.21E-11*(TEMP/300.0)**(-0.9) \n') 
    f.write('    krc     = kc0/kci \n') 
    f.write('    fcc     = 0.30 \n') 
    f.write('    nc      = 0.75-(1.27*LOG10(fcc)) \n') 
    f.write('    fc      = 10**(LOG10(fcc)/(1.0+((LOG10(krc))/nc)**2.0)) \n') 
    f.write('    KFPAN   = (kc0*kci)*fc/(kc0+kci) \n') 
    #mcm_constants_dict['KFPAN']=kfpan    

    f.write('    kd0     = 4.90E-03*M*EXP(-12100.0/TEMP) \n') 
    f.write('    kdi     = 5.40E16*EXP(-13830.0/TEMP) \n') 
    f.write('    krd     = kd0/kdi \n') 
    f.write('    fcd     = 0.30 \n') 
    f.write('    ncd     = 0.75-(1.27*LOG10(fcd)) \n') 
    f.write('    fd      = 10.0**(LOG10(fcd)/(1.0+((LOG10(krd))/ncd)**2.0)) \n') 
    f.write('    KBPAN   = (kd0*kdi)*fd/(kd0+kdi) \n') 
    #mcm_constants_dict['KBPAN']=kbpan
    
    f.write('    ! kmt01  : o        + no      = no2 \n')
    f.write('    ! iupac 2001 (mcmv3.2) \n')
    f.write('    k10     = 1.00E-31*M*(TEMP/300.0)**(-1.6) \n') 

    f.write('    k1i     = 3.00E-11*(TEMP/300.0)**(0.3) \n') 
    f.write('    kr1     = k10/k1i \n') 
    f.write('    fc1     = 0.85 \n') 
    f.write('    nc1     = 0.75-(1.27*LOG10(fc1)) \n') 
    f.write('    f1      = 10.0**(LOG10(fc1)/(1.0+((LOG10(kr1)/nc1))**2.0)) \n') 
    f.write('    KMT01   = (k10*k1i)*f1/(k10+k1i) \n') 
    #mcm_constants_dict['KMT01']=kmt01    

    f.write('    ! kmt02  : o        + no2     = no3 \n')
    f.write('    ! iupac 2001 (mcmv3.2) \n')
    f.write('    k20     = 1.30E-31*M*(TEMP/300.0)**(-1.5) \n') 
    f.write('    k2i     = 2.30E-11*(TEMP/300.0)**(0.24) \n') 
    f.write('    kr2     = k20/k2i \n') 
    f.write('    fc2     = 0.6 \n') 
    f.write('    nc2     = 0.75-(1.27*LOG10(fc2)) \n') 
    f.write('    f2      = 10.0**(LOG10(fc2)/(1.0+((LOG10(kr2)/nc2))**2.0)) \n') 
    f.write('    KMT02   = (k20*k2i)*f2/(k20+k2i) \n') 
    #mcm_constants_dict['KMT02']=kmt02
    

    f.write('    ! kmt03  : no2      + no3     = n2o5 \n')
    f.write('    ! iupac 2006, mcmv3.2 \n')
    f.write('    k30     = 3.60E-30*M*(TEMP/300.0)**(-4.1)  \n') 
    f.write('    k3i     = 1.90E-12*(TEMP/300.0)**(0.2) \n') 
    f.write('    kr3     = k30/k3i \n') 
    f.write('    fc3     = 0.35 \n') 
    f.write('    nc3     = 0.75-(1.27*LOG10(fc3)) \n') 
    f.write('    f3      = 10.0**(LOG10(fc3)/(1.0+((LOG10(kr3)/nc3))**2.0)) \n') 
    f.write('    KMT03   = (k30*k3i)*f3/(k30+k3i) \n') 
    #mcm_constants_dict['KMT03']=kmt03

    f.write('    ! kmt04  : n2o5               = no2     + no3 \n')
    f.write('    ! iupac 2006, mcmv3.2 \n')
    f.write('    k40     = 1.30E-03*M*(TEMP/300.0)**(-3.5)*EXP(-11000.0/TEMP) \n') 
    f.write('    k4i     = 9.70E14*(TEMP/300.0)**(0.1)*EXP(-11080.0/TEMP) \n') 
    f.write('    kr4     = k40/k4i \n') 
    f.write('    fc4     = 0.35 \n') 
    f.write('    nc4     = 0.75-(1.27*LOG10(fc4)) \n') 
    f.write('    f4      = 10.0**(LOG10(fc4)/(1+((LOG10(kr4)/nc4))**2.0)) \n') 
    f.write('    KMT04   = (k40*k4i)*f4/(k40+k4i) \n') 
    #mcm_constants_dict['KMT04']=kmt04

    f.write('    ! kmt05  : oh       + co(+o2) = ho2     + co2 \n')
    f.write('    ! iupac 2006 \n')
    f.write('    KMT05  = (1.0 + (M/4.2E19)) \n') 
    #mcm_constants_dict['KMT05']=kmt05

    f.write('    ! kmt06  : ho2      + ho2     = h2o2    + o2 \n')
    f.write('    ! water enhancement factor \n')
    f.write('    ! iupac 1992 \n')
    f.write('    KMT06  = 1.0 + (1.40E-21*EXP(2200.0/TEMP)*H2O) \n') 
    #mcm_constants_dict['KMT06']=kmt06

    f.write('    ! kmt06  = 1.0 + (2.00E-25*numpy.exp(4670.0/temp)*h2o) \n')
    f.write('    ! S+R 2005 values \n')

    f.write('    ! kmt07  : oh       + no      = hono \n')

    f.write('    ! iupac 2006, mcmv3.2 \n')
    f.write('    k70     = 7.40E-31*M*(TEMP/300.0)**(-2.4) \n') 
    f.write('    k7i     = 3.30E-11*(TEMP/300.0)**(-0.3) \n') 
    f.write('    kr7     = k70/k7i \n') 
    f.write('    fc7     = EXP(-TEMP/1420.0) \n') 
    f.write('    nc7     = 0.75-(1.27*LOG10(fc7)) \n') 
    f.write('    f7      = 10.0**(LOG10(fc7)/(1+((LOG10(kr7)/nc7))**2.0)) \n') 
    f.write('    KMT07   = (k70*k7i)*f7/(k70+k7i) \n') 
    #mcm_constants_dict['KMT07']=kmt07

    f.write('    ! kmt08  : oh       + no2     = hno3 \n')
    
    f.write('    ! iupac 2006, mcmv3.2 \n')
    f.write('    k80     = 3.30E-30*M*(TEMP/300.0)**(-3.0) \n') 
    f.write('    k8i     = 4.10E-11 \n') 
    f.write('    kr8     = k80/k8i \n') 
    f.write('    fc8     = 0.4 \n') 
    f.write('    nc8     = 0.75-(1.27*LOG10(fc8)) \n') 
    f.write('    f8      = 10.0**(LOG10(fc8)/(1.0+((LOG10(kr8)/nc8))**2.0)) \n') 
    f.write('    KMT08   = (k80*k8i)*f8/(k80+k8i) \n') 
    #mcm_constants_dict['KMT08']=kmt08

    f.write('    ! kmt09  : ho2      + no2     = ho2no2 \n')
    f.write('    ! iupac 1997, mcmv3.2 \n')

    f.write('    k90     = 1.80E-31*M*(TEMP/300.0)**(-3.2) \n') 
    f.write('    k9i     = 4.70E-12 \n') 
    f.write('    kr9     = k90/k9i \n') 
    f.write('    fc9     = 0.6 \n') 
    f.write('    nc9     = 0.75-(1.27*LOG10(fc9)) \n') 
    f.write('    f9      = 10.0**(LOG10(fc9)/(1.0+((LOG10(kr9)/nc9))**2.0)) \n') 
    f.write('    KMT09   = (k90*k9i)*f9/(k90+k9i) \n') 
    #mcm_constants_dict['KMT09']=kmt09

    f.write('    ! kmt10  : ho2no2             = ho2     + no2 \n')
    f.write('    ! iupac 1997, mcmv3.2 \n')

    f.write('    k100     = 4.10E-05*M*EXP(-10650.0/TEMP) \n') 
    f.write('    k10i     = 4.80E15*EXP(-11170.0/TEMP) \n') 
    f.write('    kr10     = k100/k10i \n') 
    f.write('    fc10     = 0.6 \n') 
    f.write('    nc10     = 0.75-(1.27*LOG10(fc10)) \n') 
    f.write('    f10      = 10.0**(LOG10(fc10)/(1.0+((LOG10(kr10)/nc10))**2.0)) \n') 
    f.write('    KMT10    = (k100*k10i)*f10/(k100+k10i) \n') 
    #mcm_constants_dict['KMT10']=kmt10

    f.write('    ! kmt11  : oh       + hno3    = h2o     + no3 \n')
    f.write('    ! iupac 2006, mcmv3.2 \n')
    
    f.write('    k1     = 2.40E-14*EXP(460.0/TEMP) \n') 
    f.write('    k3     = 6.50E-34*EXP(1335.0/TEMP) \n') 
    f.write('    k4     = 2.70E-17*EXP(2199.0/TEMP) \n') 
    f.write('    k2     = (k3*M)/(1.0+(k3*M/k4)) \n') 
    f.write('    KMT11  = k1 + k2 \n') 
    #mcm_constants_dict['KMT11']=kmt11

    f.write('    ! kmt12 iupac 2006, mcmv3.2 \n')

    f.write('    k120 = 4.50E-31*((TEMP/300.0)**(-3.9))*M \n') 
    f.write('    k12i = 1.30E-12*((TEMP/300.0)**(-0.7)) \n') 
    f.write('    kr12 = k120/k12i \n') 
    f.write('    fc12 = 0.525 \n') 
    f.write('    nc12 = 0.75-(1.27*LOG10(fc12)) \n') 
    f.write('    f12  = 10.0**(LOG10(fc12)/(1.0+((LOG10(kr12)/nc12))**2.0)) \n') 
    f.write('    KMT12    = (k120*k12i)*f12/(k120+k12i) \n') 
    #mcm_constants_dict['KMT12']=kmt12

    f.write('    ! kmt13  : ch3o2    + no2     = ch3o2no2 \n')
    f.write('    ! iupac 2006 \n')
    
    f.write('    k130     = 2.50E-30*((TEMP/300.0)**(-5.5))*M \n') 
    f.write('    k13i     = 1.80E-11 \n') 
    f.write('    kr13     = k130/k13i \n') 
    f.write('    fc13     = 0.36 \n') 
    f.write('    nc13     = 0.75-(1.27*LOG10(fc13)) \n') 
    f.write('    f13      = 10.0**(LOG10(fc13)/(1.0+((LOG10(kr13)/nc13))**2.0)) \n') 
    f.write('    KMT13    = (k130*k13i)*f13/(k130+k13i) \n') 
    #mcm_constants_dict['KMT13']=kmt13

    f.write('    ! kmt14  : ch3o2no2           = ch3o2   + no2 \n')
    f.write('    ! iupac 2006, mcmv3.2     \n')

    f.write('    k140     = 9.00E-05*EXP(-9690.0/TEMP)*M \n') 
    f.write('    k14i     = 1.10E16*EXP(-10560.0/TEMP) \n') 
    f.write('    kr14     = k140/k14i \n') 
    f.write('    fc14     = 0.4 \n') 
    f.write('    nc14     = 0.75-(1.27*LOG10(fc14)) \n') 
    f.write('    f14      = 10.0**(LOG10(fc14)/(1.0+((LOG10(kr14)/nc14))**2.0)) \n') 
    f.write('    KMT14    = (k140*k14i)*f14/(k140+k14i) \n') 
    #mcm_constants_dict['KMT14']=kmt14

    f.write('    ! kmt15 iupac 2006, mcmv3.2 \n')

    f.write('    k150 = 8.60E-29*((TEMP/300.0)**(-3.1))*M \n') 
    f.write('    k15i = 9.00E-12*((TEMP/300.0)**(-0.85)) \n') 
    f.write('    kr15 = k150/k15i \n') 
    f.write('    fc15 = 0.48 \n') 
    f.write('    nc15 = 0.75-(1.27*LOG10(fc15)) \n') 
    f.write('    f15  = 10.0**(LOG10(fc15)/(1.0+((LOG10(kr15)/nc15))**2.0)) \n') 
    f.write('    KMT15 = (k150*k15i)*f15/(k150+k15i) \n') 
    #mcm_constants_dict['KMT15']=kmt15

    f.write('    ! kmt16  :  oh  +  c3h6 \n')
    f.write('    ! iupac 2006 \n')

    f.write('    k160     = 8.00E-27*((TEMP/300.0)**(-3.5))*M \n') 
    f.write('    k16i     = 3.00E-11*((TEMP/300.0)**(-1.0)) \n') 
    f.write('    kr16     = k160/k16i \n') 
    f.write('    fc16     = 0.5 \n') 
    f.write('    nc16     = 0.75-(1.27*LOG10(fc16)) \n') 
    f.write('    f16      = 10.0**(LOG10(fc16)/(1.0+((LOG10(kr16)/nc16))**2.0)) \n') 
    f.write('    KMT16    = (k160*k16i)*f16/(k160+k16i) \n') 
    #mcm_constants_dict['KMT16']=kmt16

    f.write('    ! kmt17 iupac 2006 \n')

    f.write('    k170 = 5.00E-30*((TEMP/300.0)**(-1.5))*M \n') 
    f.write('    k17i = 1.00E-12 \n') 
    f.write('    kr17 = k170/k17i \n') 
    f.write('    fc17 = (0.17*EXP(-51.0/TEMP))+EXP(-1.0*TEMP/204.0) \n') 
    f.write('    nc17 = 0.75-(1.27*LOG10(fc17)) \n') 
    f.write('    f17  = 10.0**(LOG10(fc17)/(1.0+((LOG10(kr17)/nc17))**2.0)) \n') 
    f.write('    KMT17 = (k170*k17i)*f17/(k170+k17i) \n') 
    #mcm_constants_dict['KMT17']=kmt17


    f.write('    ! #Taken from yet another MCM version, will need converting to lower case \n') 
    f.write('    KRO2NO = 2.7E-12*EXP(360.0/TEMP) \n') 
    #mcm_constants_dict['KRO2NO']=KRO2NO

    f.write('    KRO2HO2 = 2.91E-13*EXP(1300.0/TEMP) \n') 
    #mcm_constants_dict['KRO2HO2']=KRO2HO2

    f.write('    KAPHO2 = 5.2E-13*EXP(980.0/TEMP) \n') 
    #mcm_constants_dict['KAPHO2']=KAPHO2

    f.write('    KAPNO = 7.5E-12*EXP(290.0/TEMP) \n') 
    #mcm_constants_dict['KAPNO']=KAPNO

    f.write('    KRO2NO3 = 2.3E-12 \n') 
    #mcm_constants_dict['KRO2NO3']=KRO2NO3

    f.write('    KNO3AL = 1.4E-12*EXP(-1860.0/TEMP) \n') 
    #mcm_constants_dict['KNO3AL']=KNO3AL

    f.write('    KDEC = 1.00E06 \n') 
    #mcm_constants_dict['KDEC']=KDEC
        

    f.write('    KROPRIM = 2.50E-14*EXP(-300.0/TEMP) \n') 
    #mcm_constants_dict['KROPRIM']=KROPRIM

    f.write('    KROSEC = 2.50E-14*EXP(-300.0/TEMP) \n') 
    #mcm_constants_dict['KROSEC']=KROSEC

    f.write('    KCH3O2 = 1.03E-13*EXP(365.0/TEMP) \n') 
    #mcm_constants_dict['KCH302']=KCH3O2

    f.write('    K298CH3O2 = 3.5E-13 \n') 
    #mcm_constants_dict['K298CH3O2']=K298CH3O2

    f.write('    KMT18 = 9.5E-39*O2*EXP(5270.0/TEMP)/(1+7.5E-29*O2*EXP(5610.0/TEMP)) \n') 
    #mcm_constants_dict['KMT18']=kmt18

    f.write('    KROPRIM  = 3.70E-14*EXP(-460.0/TEMP) \n') 
    #mcm_constants_dict['KROPRIM']=kroprim

    f.write('    KROSEC   = 1.80E-14*EXP(-260.0/TEMP) \n') 
    #mcm_constants_dict['KROSEC']=krosec
    
    #f.write('    rate_values=numpy.zeros(%s)\n' %(len(rate_dict.keys()))) 
    # Now cycle through the rate_dict dictionary and print to file 
    f.write('    ! Explicitly calculating rates for each reaction \n') 
    

    if openMP is True:
        # Workout the number of cores on this machine so you can use a speific number
        # of single clauses
        cores=multiprocessing.cpu_count()
        # Now work out the equations at which to define a new clause
        break_list=[int(x) for x in np.linspace(0,len(rate_dict_fortran.keys()),cores)][1:-1]

    check_step=0
    for key in rate_dict_fortran.keys():
        if openMP is True and check_step == 0:
            f.write('!$OMP PARALLEL DO &\n')
            f.write('!$OMP& SHARED(rate_values, J)  \n') 
            f.write('!$OMP SINGLE \n') 
        if openMP is True and check_step in break_list:
            f.write('!$OMP END SINGLE NOWAIT \n') 
            f.write('!$OMP SINGLE \n') 
        #pdb.set_trace()
        f.write('    rate_values(%s) =%s \n' %(key+1,rate_dict_fortran[key])) 
        check_step+=1
    if openMP is True:
        f.write('!$OMP END SINGLE NOWAIT \n') 
        f.write('!$OMP END PARALLEL DO \n') 
    f.write('\n') 
    f.write('end subroutine\n')   
    f.close()  
    

    
def write_RO2_indices(filename,species_dict2array):
    
    #The following is a list of species defined as contributing to total RO2 concetrations
    RO2_names = ['NBUTOLAO2','HO3C4O2','BU1ENO3O2','C43NO34O2','BZBIPERO2','CH3O2','C2H5O2','HOCH2CH2O2',
    'ETHENO3O2','C6H5C2H4O2','EBZBIPERO2','ISOPAO2','ISOPBO2','ISOPCO2','ISOPDO2',
    'NISOPO2','CH3CO3','C2H5CO3','NC3H7O2','IC3H7O2','HYPROPO2','IPROPOLO2','PRONO3BO2',
    'PRONO3AO2','C6H5CH2O2','TLBIPERO2','CH3COCH2O2','BUT2OLO2','C42NO33O2','IC4H9O2','TC4H9O2','IBUTOLBO2','TBUTOLO2',
    'MPRANO3O2','MPRBNO3O2','IPEAO2','IPEBO2','IPECO2','MXYBIPERO2','MXYLO2','MEKAO2',
    'MEKCO2','MEKBO2','NC4H9O2','SC4H9O2','HEXAO2','HEXBO2','HEXCO2','PEAO2','PEBO2','PECO2','OXYBIPERO2','OXYLO2',
    'PXYBIPERO2','PXYLO2','MPRKAO2','CO2C54O2','HO2C5O2','DIEKAO2','DIEKBO2','BZEMUCO2',
    'BZEMUCCO3','C5DIALO2','PHENO2','NPHENO2','EBZMUCO2','EBZMUCCO3','C715CO2O2','EBENZOLO2','NEBNZOLO2','HCOCO3','HMVKAO2',
    'HMVKBO2','MVKO2','MACO3','MACRO2','TLEMUCO2','TLEMUCCO3','C615CO2O2','CRESO2',
    'NCRESO2','MXYMUCO2','MXYMUCCO3','C726CO5O2','MXYOLO2','NMXYOLO2','OXYMUCO2',
    'OXYMUCCO3','MC6CO2O2','OXYOLO2','NOXYOLO2','PXYMUCO2','PXYMUCCO3','C6M5CO2O2','PXYOLO2','NPXYOLO2','HO3C3CO3',
    'CO3C4NO3O2','MALDIALCO3','EPXDLCO3','C3DIALO2','MALDIALO2','HOCH2CO3',
    'NO3CH2CO3','C6H5CH2CO3','C6DCARBBO2','C58O2','HC4ACO3','HC4CCO3','C57O2',
    'C59O2','NC4CO3','C510O2','HO1C3O2','CH3CHOHCO3','PRNO3CO3','C6H5CO3','C5CO14O2',
    'IBUTOLCO2','IBUTALBO2','IBUTALCO2','IPRCO3','IPRHOCO3','MPRBNO3CO3','M2BUOL2O2',
    'HM2C43O2','BUT2CO3','C52O2','ME2BUOLO2','H2M3C4O2','MIPKAO2','MIPKBO2','ME2BU2OLO2',
    'PROL11MO2','HO2M2C4O2','C3MCODBCO3','MXYLCO3','EPXMDLCO3','C3MDIALO2','HO1CO3C4O2','CO2C3CO3','BIACETO2',
    'NBUTOLBO2','BUTALO2','C3H7CO3','HO1C4O2','C5H11CO3','HO1C6O2','HEX2ONAO2','HEX2ONBO2','HEX2ONCO2','HO2C6O2','HO3C6O2',
    'HEX3ONDO2','HEX3ONCO2','HEX3ONBO2','HEX3ONAO2','C4CHOBO2','C4H9CO3','HO1C5O2','PE2ENEBO2','HO3C5O2','OXYLCO3','EPXM2DLCO3',
    'C4MCO2O2','PXYLCO3','CO23C54O2','HO2CO4C5O2','CO24C53O2','HO2C4CO3','HO2C4O2','HOCO3C54O2','CO3C4CO3','BZFUO2','NBZFUO2','HCOCOHCO3','EBFUO2','BUTALAO2',
    'NEBFUO2','C6DICARBO2','C7CO3OHO2','MVKOHBO2','MVKOHAO2','CO2H3CO3','ACO3','TLFUO2','NTLFUO2','C5DICARBO2','MC3CODBCO3','C4M2ALOHO2','C4MCODBCO3','C5MCO2OHO2',
    'MXYFUO2','C23O3MO2','NMXYFUO2','PXYFUO2','MCOCOMOXO2','NPXYFUO2','MC5CO2OHO2','MC4CODBCO3','OXYFUO2','C6OTKETO2',
    'NOXYFUO2','DMKOHO2','C4CO2O2','C51O2','HCOCH2O2','PBZQO2','NBZQO2','NCATECO2','NNCATECO2','CO3H4CO3','PEBQO2','NPEBQO2','ENCATECO2','ENNCATECO2','HOC2H4CO3',
    'MECOACETO2','PTLQO2','NPTLQO2','MNCATECO2','MNNCATECO2','HOIPRCO3','IBUDIALCO3','PROPALO2','HOIBUTCO3','HO2C43CO3','C56O2',
    'C53O2','C41CO3','PROL1MCO3','H2M2C3CO3','C54O2','CHOMOHCO3','MXYQO2','NMXYQO2','MXNCATECO2','MXNNCATCO2','HO2C3CO3','HOC3H6CO3','C63O2','CO2HOC61O2','CO24C6O2',
    'CO25C6O2','C61O2','CO23C65O2','HO3C5CO3','C6HO1CO3O2','C3COCCO3','PEN2ONE1O2','C6CO34O2','C6CO3OH5O2','HO3C4CO3','HO13C5O2','TMB1FUO2','NTMB1FUO2','OXYQO2',
    'NOXYQO2','OXNCATECO2','OXNNCATCO2','PXYQO2','NPXYQO2','PXNCATECO2','PXNNCATCO2',
    'CO2C4CO3','HO13C4O2','MALANHYO2','DNPHENO2','NDNPHENO2','DNEBNZLO2','NDNEBNZLO2','H13CO2CO3','DNCRESO2','NDNCRESO2',\
    'HOBUT2CO3','ACCOCOMEO2','MMALANHYO2','DNMXYOLO2','NDNMXYOLO2','CO3C5CO3','C6O4KETO2','DNOXYOLO2','NDNOXYOLO2',
    'TL4OHNO2O2','DNPXYOLO2','NDNPXYOLO2','CO23C4CO3','HCOCH2CO3','C5CO2OHCO3','ECO3CO3','C7OHCO2CO3','ACCOMECO3',
    'CH3COCO3','C6CO2OHCO3','C42CO3','H13C43CO3','C4COMOHCO3','C23O3MCO3','C23O3CCO3',
    'C7CO2OHCO3','C62O2','C5M2OHOCO3','C6MOHCOCO3','HO13C3CO3','C4CO2DBCO3','C7CO2DBCO3','C5CO2DBCO3','C5CO234O2',
    'C5CO34CO3','C4DBM2CO3','C5DBCO2CO3','C5CO23O2','EMALANHYO2']
    
    RO2_indices=[]
    for name in RO2_names:
        if name in species_dict2array:
            RO2_indices.append(species_dict2array[name])
            
    np.save(filename+'_RO2_indices', RO2_indices)
    
    # -----------------------------------------------------------------------------------------------------------------------------
    # Placeholder for future development - create specific RO2 file for speed improvement
    #Put all of the rate coefficient functional forms into a new python file.
    #f = open('RO2_conc.py','w')
    #f.write('##################################################################################################### \n') # python will convert \n to os.linesep
    #f.write('# Python function to hold indices for calculating RO2 concentration                                 # \n') # python will convert \n to os.linesep
    #f.write('#    Copyright (C) 2017  David Topping : david.topping@manchester.ac.uk                             # \n')
    #f.write('#                                      : davetopp80@gmail.com                                       # \n')
    #f.write('#    Personal website: davetoppingsci.com                                                           # \n')
    #f.write('#                                                                                                   # \n')
    #f.write('#    This program does not have a license, meaning the deault copyright law applies.                # \n')
    #f.write('#    Only users who have access to the private repository that holds this file may                  # \n')
    #f.write('#    use it or develop it, but may not distribute it.                                               # \n')
    #f.write('#                                                                                                   # \n')
    #f.write('#                                                                                                   # \n')
    #f.write('##################################################################################################### \n')    
    #f.write('# File created at %s \n' % datetime.datetime.now()) # python will convert \n to os.linesep
    #f.write('\n') 
    #f.write('import numpy as np\n') 
    #f.write('\n') 
    #f.write('def RO2_conc(y):\n') 
    ## Now cycle through all of the mcm_constants_dict values.
    ## Please note this includes photolysis rates that will change with time of day or condition in the chamber. You will
    ## need to modify this potentially.
    #f.write('    # Sum all relevant contributions from y array\n') 
    #step=0
    #for name in RO2_names:
    #    if name in species_dict2array:
    #        if step==0:
    #            f.write('    RO2_indices =[%s,\ \n' %(species_dict2array[name]))
    #            step+=1
    #        elif step >0: 
    #            f.write('    %s,\ \n' %(species_dict2array[name]))
    #            step+=1
    #f.write('    ]\ \n')       
    #f.write('\n')     
    #f.write('    return RO2_indices \n')   
    #f.close() 
    # -----------------------------------------------------------------------------------------------------------------------------
      
def write_reactants_indices(filename,equations,num_species,species_dict2array,rate_dict_reactants,loss_dict):

    #Save the indices to a file
    reactants_indices=[]
    reactants_stoich_indices=[]
    for equation_step in range(equations):
        reactants_indices.append([0]*len(rate_dict_reactants[equation_step].items()))
        reactants_stoich_indices.append([0]*len(rate_dict_reactants[equation_step].items()))
        step=0
        for reactant_step, reactant in rate_dict_reactants[equation_step].items():
            if reactant not in ['hv']:
                stoich=loss_dict[reactant][equation_step]
                reactants_indices[equation_step][step]=species_dict2array[reactant]
                reactants_stoich_indices[equation_step][step]=stoich
                step+=1
                
    np.save(filename+'_reactant_indices', reactants_indices)
    np.save(filename+'_reactants_stoich_indices', reactants_stoich_indices)
    
    #Now create a sparse matrix of the same information, storing the stoichiometry
    #Initialise a sparse matrix 
    reactants_indices_sparse=lil_matrix((equations, num_species), )
    for equation_step in range(equations):
        for reactant_step, reactant in rate_dict_reactants[equation_step].items():
            if reactant not in ['hv']:
                stoich=loss_dict[reactant][equation_step]
                reactants_indices_sparse[equation_step,species_dict2array[reactant]]=stoich
                
    #Now convert matrix type for use in numerical operations later
    reactants_indices_sparse = reactants_indices_sparse.tocsr()
    #Now save to disk for use later
    save_sparse_csr('reactants_indices_sparse_'+filename,reactants_indices_sparse)
    
    # -----------------------------------------------------------------------------------------------------------------------------
    # Placeholder for future development - create specific reactants file for speed improvement
    #The following commented out code puts all of the rate coefficient functional forms into a new python file.
    f = open('Reactants_conc_numba.py','w')
    f.write('##################################################################################################### \n') # python will convert \n to os.linesep
    f.write('# Python function to hold indices for calculating reactant contributions to each equation           # \n') # python will convert \n to os.linesep
    f.write('#    Copyright (C) 2017  David Topping : david.topping@manchester.ac.uk                             # \n')
    f.write('#                                      : davetopp80@gmail.com                                       # \n')
    f.write('#    Personal website: davetoppingsci.com                                                           # \n')
    f.write('#                                                                                                   # \n')
    f.write('#    This program does not have a license, meaning the deault copyright law applies.                # \n')
    f.write('#    Only users who have access to the private repository that holds this file may                  # \n')
    f.write('#    use it or develop it, but may not distribute it.                                               # \n')
    f.write('#                                                                                                   # \n')
    f.write('#                                                                                                   # \n')
    f.write('##################################################################################################### \n')    
    f.write('# File created at %s \n' % datetime.datetime.now()) # python will convert \n to os.linesep
    f.write('\n') 
    f.write('import numpy as np\n') 
    f.write('import numba as nb')
    f.write('\n')
    f.write('@nb.jit(\'float64[:](float64[:],int32,float64[:])\', nopython=True)\n')
    f.write('def reactants(y,equations,reactants_conc):\n')
    f.write('    # Calculate reactant and stochiometric contribution for each reaction\n')
    f.write('    # reactants_conc=numpy.zeros((len(equations)))\n')
    #Get the required concentration products 
    for equation_step in range(equations):
        step=0
        for reactant_step, reactant in rate_dict_reactants[equation_step].items():
            if reactant not in ['hv']:
                stoich=loss_dict[reactant][equation_step]
                if step==0:
                    if stoich == 1:
                        f.write('    reactants_conc[%s]=y[%s]'%(equation_step,species_dict2array[reactant]))
                        step=1
                    else:
                        f.write('    reactants_conc[%s]=y[%s]**%s'%(equation_step,species_dict2array[reactant],stoich))
                        step=1
                else:  
                    if stoich == 1:
                        f.write('*y[%s]'%(species_dict2array[reactant]))
                    else:
                        f.write('*y[%s]**%s'%(species_dict2array[reactant],stoich))
        f.write('\n')
    f.write('    return reactants_conc \n') 
    f.close()  
    #reactants_jac[equation_step][reactant]=y[arg8[reactant]]
    #Calculate the rate coefficient for this equation
    #pdb.set_trace()
    # -----------------------------------------------------------------------------------------------------------------------------
    
def write_reactants_indices_fortran(filename,equations,species_dict2array,rate_dict_reactants,loss_dict,openMP):

    #Save the indices to a file
    reactants_indices=[]
    reactants_stoich_indices=[]
    for equation_step in range(equations):
        reactants_indices.append([0]*len(rate_dict_reactants[equation_step].items()))
        reactants_stoich_indices.append([0]*len(rate_dict_reactants[equation_step].items()))
        step=0
        for reactant_step, reactant in rate_dict_reactants[equation_step].items():
            if reactant not in ['hv']:
                stoich=loss_dict[reactant][equation_step]
                reactants_indices[equation_step][step]=species_dict2array[reactant]
                reactants_stoich_indices[equation_step][step]=stoich
                step+=1
                
    np.save(filename+'_reactant_indices', reactants_indices)
    np.save(filename+'_reactants_stoich_indices', reactants_stoich_indices)

    if openMP is True:
        # Workout the number of cores on this machine so you can use a speific number
        # of single clauses
        cores=multiprocessing.cpu_count()
        # Now work out the equations at which to define a new clause
        break_list=[int(x) for x in np.linspace(0,equations,cores)][1:-1]
    
    # -----------------------------------------------------------------------------------------------------------------------------
    # Placeholder for future development - create specific reactants file for speed improvement
    #The following commented out code puts all of the rate coefficient functional forms into a new python file.
    f = open('Reactants_conc.f90','w')
    #f.write('!#################################################################################################### \n') # python will convert \n to os.linesep
    #f.write('! Fortran function to hold indices for calculating reactant contributions to each equation           # \n') # python will convert \n to os.linesep
    #f.write('!    Copyright (C) 2017  David Topping : david.topping@manchester.ac.uk                             # \n')
    #f.write('!                                      : davetopp80@gmail.com                                       # \n')
    #f.write('!    Personal website: davetoppingsci.com                                                           # \n')
    #f.write('!                                                                                                   # \n')
    #f.write('!    This program does not have a license, meaning the deault copyright law applies.                # \n')
    #f.write('!    Only users who have access to the private repository that holds this file may                  # \n')
    #f.write('!    use it or develop it, but may not distribute it.                                               # \n')
    #f.write('!                                                                                                   # \n')
    #f.write('!                                                                                                   # \n')
    #f.write('!#################################################################################################### \n')    
    #f.write('! File created at %s \n' % datetime.datetime.now()) # python will convert \n to os.linesep
    #f.write('\n') 
    #f.write('import numpy as np\n') 
    f.write('\n')
    f.write('subroutine reactants(y,reactants_conc)\n')
    f.write('    implicit none \n')     
    f.write('    ! Calculate reactant and stochiometric contribution for each reaction\n')
    f.write('    REAL(8), intent(in), dimension(%s) :: y \n' %(str(len(species_dict2array.keys()))) )     
    f.write('    !f2py intent(in) :: y \n') 
    f.write('    REAL(8), intent(out), dimension(%s) :: reactants_conc \n' %(str(equations)) ) 
    f.write('    !f2py intent(out) :: reactants_conc \n') 
    #Get the required concentration products 
    for equation_step in range(equations):
        step=0
        if openMP is True and equation_step == 0:
            f.write('!$OMP PARALLEL DO &\n')
            f.write('!$OMP& SHARED(y,reactants_conc)  \n') 
            f.write('!$OMP SINGLE \n') 
        if openMP is True and equation_step in break_list:
            f.write('!$OMP END SINGLE NOWAIT \n') 
            f.write('!$OMP SINGLE \n') 
        for reactant_step, reactant in rate_dict_reactants[equation_step].items():
            if reactant not in ['hv']:
                stoich=loss_dict[reactant][equation_step]
                if step==0:
                    if stoich == 1:
                        f.write('    reactants_conc(%s)=y(%s)'%(equation_step+1,species_dict2array[reactant]+1))
                        step=1
                    else:
                        f.write('    reactants_conc(%s)=y(%s)**%s'%(equation_step+1,species_dict2array[reactant]+1,stoich))
                        step=1
                else:  
                    if stoich == 1:
                        f.write('*y(%s)'%(species_dict2array[reactant]+1))
                    else:
                        f.write('*y(%s)**%s'%(species_dict2array[reactant]+1,stoich))
        f.write('\n')
    if openMP is True:
        f.write('!$OMP END SINGLE NOWAIT \n') 
        f.write('!$OMP END PARALLEL DO \n') 
    f.write('end subroutine\n') 
    f.close()  
    #reactants_jac[equation_step][reactant]=y[arg8[reactant]]
    #Calculate the rate coefficient for this equation
    #pdb.set_trace()
    # -----------------------------------------------------------------------------------------------------------------------------
        
def save_sparse_csr(filename,array):
    np.savez(filename,data = array.data ,indices=array.indices,indptr =array.indptr, shape=array.shape )
                
def write_loss_gain_matrix(filename,equations,num_species,loss_dict,gain_dict,species_dict2array):
    #Here we create and then store a sparse matrix for use in the loss_gain calculations
    #This is then loaded before the simulation to calculate loss and gain
            
    #Initialise a sparse matrix 
    loss_gain=lil_matrix((num_species, equations), )
    
    for species, species_step in species_dict2array.items():
        if species not in ['AIR','H2O','O2','hv']:
            #Now add values
            for equation_step, stoich in loss_dict[species].items():
                if stoich==1:
                    #f.write('-1.0*r[%s]'%(equation_step))
                    loss_gain[species_step,equation_step]=-1.0
                else:
                    #f.write('-1.0*r[%s]*%s'%(equation_step,stoich))
                    loss_gain[species_step,equation_step]=-1.0*stoich
            for equation_step, stoich in gain_dict[species].items():
                if stoich==1:
                    #f.write('+r[%s]'%(equation_step))
                    loss_gain[species_step,equation_step]=1.0
                else:
                    #f.write('+r[%s]*%s'%(equation_step,stoich))
                    loss_gain[species_step,equation_step]=1.0*stoich
                
    #Now convert matrix type for use in numerical operations later
    loss_gain = loss_gain.tocsr()
    #Now save to disk for use later
    save_sparse_csr('loss_gain_'+filename,loss_gain)
    
    # -----------------------------------------------------------------------------------------------------------------------------
    # Placeholder for future development - create specific reactants file for speed improvement
    f = open('Loss_Gain_numba.py','w')
    f.write('##################################################################################################### \n') # python will convert \n to os.linesep
    f.write('# Python function to hold indices for calculating loss and gain for each species                    # \n') # python will convert \n to os.linesep
    f.write('#    Copyright (C) 2017  David Topping : david.topping@manchester.ac.uk                             # \n')
    f.write('#                                      : davetopp80@gmail.com                                       # \n')
    f.write('#    Personal website: davetoppingsci.com                                                           # \n')
    f.write('#                                                                                                   # \n')
    f.write('#    This program does not have a license, meaning the deault copyright law applies.                # \n')
    f.write('#    Only users who have access to the private repository that holds this file may                  # \n')
    f.write('#    use it or develop it, but may not distribute it.                                               # \n')
    f.write('#                                                                                                   # \n')
    f.write('#                                                                                                   # \n')
    f.write('##################################################################################################### \n')    
    f.write('# File created at %s \n' % datetime.datetime.now()) # python will convert \n to os.linesep
    f.write('\n') 
    f.write('import numpy as np\n') 
    f.write('import numba as nb')
    f.write('\n')
    f.write('@nb.jit(\'float64[:](float64[:],float64[:])\', nopython=True)\n')
    f.write('def dydt(loss_gain_list,r):\n')
    f.write('    # Calculate loss and gain for each specie\n')
    f.write('    # loss_gain_list=[0]*num_species\n')
    for species, species_step in species_dict2array.items():
        
        if species not in ['AIR','H2O','O2','hv']:
            if species in loss_dict.keys() or species in gain_dict.keys():
                f.write('    loss_gain_list[%s]='%(species_step))
                if species in loss_dict.keys():
                    loss_step=0
                    for equation_step, stoich in loss_dict[species].items():
                        if loss_step==0:
                            if stoich==1:
                                f.write('-1.0*r[%s]'%(equation_step))
                            else:
                                f.write('-1.0*r[%s]*%s'%(equation_step,stoich))
                            loss_step+=1
                        else:
                            if stoich==1:
                                f.write('-r[%s]'%(equation_step))
                            else:
                                f.write('-r[%s]*%s'%(equation_step,stoich))
                            loss_step+=1
                        #loss[species_step]+=reactants[equation_step]*stoich
                if species in gain_dict.keys():
                    gain_step=0
                    for equation_step, stoich in gain_dict[species].items():
                        if gain_step==0:
                             if stoich==1:
                                f.write('+r[%s]'%(equation_step))
                             else:
                                f.write('+r[%s]*%s'%(equation_step,stoich))
                             gain_step+=1
                        else:
                             if stoich==1:
                                f.write('+r[%s]'%(equation_step))
                             else:
                                f.write('+r[%s]*%s'%(equation_step,stoich))
                             gain_step+=1
                        #gain[species_step]+=reactants[equation_step]*stoich
        f.write('\n')
         #dydt_list[species_step]+=(gain[species_step]-loss[species_step])
    f.write('    return loss_gain_list \n') 
    f.close()  
    #pdb.set_trace()
    # -----------------------------------------------------------------------------------------------------------------------------
    
def write_loss_gain_fortran(filename,equations,num_species,loss_dict,gain_dict,species_dict2array,openMP):
    #Here we create and then store a sparse matrix for use in the loss_gain calculations
    #This is then loaded before the simulation to calculate loss and gain
            
    f = open('Loss_Gain.f90','w')
    f.write('!##################################################################################################### \n') # python will convert \n to os.linesep
    f.write('! Fortran function to hold indices for calculating loss and gain for each species                    # \n') # python will convert \n to os.linesep
    f.write('!    Copyright (C) 2018  David Topping : david.topping@manchester.ac.uk                              # \n')
    f.write('!                                      : davetopp80@gmail.com                                        # \n')
    f.write('!    Personal website: davetoppingsci.com                                                            # \n')
    f.write('!                                                                                                    # \n')
    f.write('!    This program does not have a license, meaning the deault copyright law applies.                 # \n')
    f.write('!    Only users who have access to the private repository that holds this file may                   # \n')
    f.write('!    use it or develop it, but may not distribute it.                                                # \n')
    f.write('!                                                                                                    # \n')
    f.write('!                                                                                                    # \n')
    f.write('!##################################################################################################### \n')    
    f.write('! File created at %s \n' % datetime.datetime.now()) # python will convert \n to os.linesep
    f.write('\n') 
    #f.write('import numpy as np\n') 
    f.write('\n')
    f.write('subroutine loss_gain(r,loss_gain_array)\n')
    f.write('    ! Calculate loss and gain for each specie\n')
    #f.write('    loss_gain_list=[0]*num_species\n')
    #for species, species_step in species_dict2array.items():
    #    
    f.write('    implicit none \n')         
    f.write('    REAL(8), intent(in), dimension(%s) :: r \n' %(str(equations)) )     
    f.write('    !f2py intent(in) :: r \n') 
    f.write('    REAL(8), intent(out), dimension(%s) :: loss_gain_array \n' %(str(num_species)) ) 
    f.write('    !f2py intent(out) :: loss_gain_array \n') 

    if openMP is True:
        # Workout the number of cores on this machine so you can use a speific number
        # of single clauses
        cores=multiprocessing.cpu_count()
        # Now work out the points at which to define a new clause
        break_list=[int(x) for x in np.linspace(0,num_species,cores)][1:-1]

    check_step=0
    
    for species, species_step in species_dict2array.items():

        if openMP is True and check_step == 0:
            f.write('!$OMP PARALLEL DO &\n')
            f.write('!$OMP& SHARED(r,loss_gain_array)  \n') 
            f.write('!$OMP SINGLE \n') 
        if openMP is True and check_step in break_list:
            f.write('!$OMP END SINGLE NOWAIT \n') 
            f.write('!$OMP SINGLE \n') 
            
        if species not in ['AIR','H2O','O2','hv']:
            if species in loss_dict.keys() or species in gain_dict.keys():
                f.write('    loss_gain_array(%s)='%(species_step+1))
                flag=0
                if species in loss_dict.keys():
                    loss_step=0
                    for equation_step, stoich in loss_dict[species].items():
                        if flag > 4:
                            f.write(' & \n')
                            f.write('    ')
                            flag=0
                        if loss_step==0:
                            if stoich==1:
                                f.write('-1.0*r(%s)'%(equation_step+1))
                            else:
                                f.write('-1.0*r(%s)*%s'%(equation_step+1,stoich+1))
                            loss_step+=1
                        else:
                            if stoich==1:
                                f.write('-r(%s)'%(equation_step+1))
                            else:
                                f.write('-r(%s)*%s'%(equation_step+1,stoich+1))
                            loss_step+=1
                        flag+=1
                        
                        #loss[species_step]+=reactants[equation_step]*stoich
                if species in gain_dict.keys():
                    gain_step=0
                            
                    for equation_step, stoich in gain_dict[species].items():
                        if flag > 4:
                            f.write(' & \n')
                            f.write('    ')
                            flag=0
                        if gain_step==0:
                            if stoich==1:
                                f.write('+r(%s)'%(equation_step+1))
                            else:
                                f.write('+r(%s)*%s'%(equation_step+1,stoich+1))
                            gain_step+=1
                        else:
                            if stoich==1:
                                f.write('+r(%s)'%(equation_step+1))
                            else:
                                f.write('+r(%s)*%s'%(equation_step+1,stoich+1))
                            gain_step+=1
                        flag+=1
                        
                        #gain[species_step]+=reactants[equation_step]*stoich
        f.write(' \n')
        check_step+=1
        #dydt_list[species_step]+=(gain[species_step]-loss[species_step])
    if openMP is True:
        f.write('!$OMP END SINGLE NOWAIT \n') 
        f.write('!$OMP END PARALLEL DO \n') 
    f.write('end subroutine \n') 
    f.close()  
    #pdb.set_trace()
    # -----------------------------------------------------------------------------------------------------------------------------

def write_gas_jacobian_fortran(filename,equations,num_species,loss_dict,gain_dict,species_dict2array,rate_dict_reactants,OpenMP):
    
    # Function to calculate the jacobian for the gas phase only model. Designed to provide some
    # improvements in comp efficiency. Cant be done for the full aerosol model due to dependency on
    # properties/processes that will change
    
    f = open('Jacobian.f90','w')
    f.write('!##################################################################################################### \n') # python will convert \n to os.linesep
    f.write('! Fortran function to calculate jacobian function                                                    # \n') # python will convert \n to os.linesep
    f.write('!    Copyright (C) 2018  David Topping : david.topping@manchester.ac.uk                              # \n')
    f.write('!                                      : davetopp80@gmail.com                                        # \n')
    f.write('!    Personal website: davetoppingsci.com                                                            # \n')
    f.write('!                                                                                                    # \n')
    f.write('!    This program does not have a license, meaning the deault copyright law applies.                 # \n')
    f.write('!    Only users who have access to the private repository that holds this file may                   # \n')
    f.write('!    use it or develop it, but may not distribute it.                                                # \n')
    f.write('!                                                                                                    # \n')
    f.write('!§# \n')
    f.write('!##################################################################################################### \n')    
    f.write('! File created at %s \n' % datetime.datetime.now()) # python will convert \n to os.linesep
    f.write('\n') 
    #f.write('import numpy as np\n') 
    f.write('\n')
    f.write('subroutine jacobian(r,y,dy_dy)\n')
    f.write('    ! Calculate loss and gain for each specie\n')
    #f.write('    loss_gain_list=[0]*num_species\n')
    #for species, species_step in species_dict2array.items():
    #    
    f.write('    implicit none \n')         
    f.write('    REAL(8), intent(in), dimension(%s) :: r \n' %(str(equations)) )     
    f.write('    !f2py intent(in) :: r \n') 
    f.write('    REAL(8), intent(in), dimension(%s) :: y \n' %(str(num_species)) ) 
    f.write('    !f2py intent(in) :: loss_gain_array \n') 
    f.write('    REAL(8), intent(out), dimension(%s,%s) :: dy_dy \n' %(str(num_species),str(num_species)) ) 
    f.write('    !f2py intent(out) :: dy_dy \n') 
    f.write('    dy_dy(:,:)=0.0D0 \n') 
    if OpenMP is True:
        # Workout the number of cores on this machine so you can use a specific number
        # of single clauses
        cores=multiprocessing.cpu_count()
        # Now work out the points at which to define a new clause
        break_list=[int(x) for x in np.linspace(0,num_species,cores)][1:-1]

    check_step=0

    for species, species_step in species_dict2array.items():

        if OpenMP is True and check_step == 0:
            f.write('!$OMP PARALLEL DO &\n')
            f.write('!$OMP& SHARED(r,y,dy_dy)  \n') 
            f.write('!$OMP SINGLE \n') 
        if OpenMP is True and check_step in break_list:
            f.write('!$OMP END SINGLE NOWAIT \n') 
            f.write('!$OMP SINGLE \n') 

        #pdb.set_trace()

        if species not in ['AIR','H2O','O2','hv']:

            # Now we need to extract other species that affect this specie
            # To do this we collect all reactions that this is involved in.
            # Then, pull out each other specie in those reactins to calculate d-species_d-otherspecie
            
            if species in loss_dict.keys():
                equations_list_loss=[num for num, stoich in loss_dict[species].items()]
                # Now extract all reactants in these equations as a list
                reactants_list_loss=[]
                for equation in equations_list_loss:
                    reactants_list_loss=reactants_list_loss+[reactant for step, reactant in rate_dict_reactants[equation].items()]
                # Now turn this into a unique list
                reactants_list_loss=set(reactants_list_loss)
                reactants_list_loss=list(reactants_list_loss)
            
            if species in gain_dict.keys():
                equations_list_gain=[num for num, stoich in gain_dict[species].items()]
                # Now extract all reactants in these equations as a list
                reactants_list_gain=[]
                for equation in equations_list_gain:
                    reactants_list_gain=reactants_list_gain+[reactant for step, reactant in rate_dict_reactants[equation].items()]
                # Now turn this into a unique list
                reactants_list_gain=set(reactants_list_gain)
                reactants_list_gain=list(reactants_list_gain)
            
            # Now combine these two lists and again merge into a unique one
            final_list=reactants_list_loss+reactants_list_gain
            final_list=set(final_list)
            final_list=list(final_list)
            
            # Now cycle through each of these species to calculate dy_dt_dx
            
            for reactant in final_list:
                species_step2=species_dict2array[reactant]
                f.write('    dy_dy(%s,%s)='%(species_step+1,species_step2+1)) #+1 due to change in indexing from Python to Fortran
                equation_flag=0
                loss_step=0
                gain_step=0

                if len(equations_list_loss) > 0:
                    for equation in equations_list_loss:
                        reactants_list_loss_temp=[reactant for step, reactant in rate_dict_reactants[equation].items()]
                        # Check if 'otherspecie' is in this reaction
                        if reactant in reactants_list_loss_temp:
                            if equation_flag > 3: # this is just to ensure line length isnt a problem for Fortran compiler
                                f.write(' & \n')
                                f.write('    ')  
                                equation_flag=0
                            step=1
                            stoich_first=loss_dict[reactant][equation]
                            # we now check the stochiometry to see if we need a variable before and reduction in power
                            # if our specie only has a stoich of 1 we dont need to consider it in the proceeding calculations
                            if stoich_first==1:
                                for reactant_step, reactant2 in rate_dict_reactants[equation].items():
                                
                                    if reactant2 not in ['hv'] and reactant2 != reactant:
                                        stoich=loss_dict[reactant2][equation]
                                        if stoich == 1:
                                            if step == 1:
                                                f.write('-1.0*y(%s)'%(species_dict2array[reactant2]+1))
                                                step+=1
                                            else:
                                                f.write('*y(%s)'%(species_dict2array[reactant2]+1))
                                        else:
                                            if step == 1:
                                                f.write('-1.0*y(%s)**%s'%(species_dict2array[reactant2]+1,stoich))
                                                step+=1
                                            else:
                                                f.write('*y(%s)**%s'%(species_dict2array[reactant2]+1,stoich))
                            else: # If our specie has a stoich of > 1 we need to include it in the expression
                                #f.write('-%s*'%(stoich_first))
                                for reactant_step, reactant2 in rate_dict_reactants[equation].items():
                                    #step=1
                                    if reactant2 not in ['hv']:
                                        if reactant2 == reactant:
                                            if stoich_first == 2:
                                                if step == 1:
                                                    f.write('-%s*y(%s)'%(stoich_first,species_dict2array[reactant2]+1))
                                                    step+=1
                                                else:
                                                    f.write('%s*y(%s)'%(stoich_first,species_dict2array[reactant2]+1))
                                            elif stoich_first > 2:
                                                if step == 1:
                                                    f.write('-%s*y(%s)**%s'%(stoich_first,species_dict2array[reactant2]+1,stoich_first-1))
                                                    step+=1
                                                else:
                                                    f.write('%s*y(%s)**%s'%(stoich_first,species_dict2array[reactant2]+1,stoich_first-1))
                                        elif reactant2 != reactant:
                                            stoich=loss_dict[reactant2][equation]
                                            if stoich == 1:
                                                if step == 1:
                                                    f.write('-1.0*y(%s)'%(species_dict2array[reactant2]+1))
                                                    step+=1
                                                else:
                                                    f.write('*y(%s)'%(species_dict2array[reactant2]+1))
                                            else:
                                                if step == 1:
                                                    f.write('-1.0*y(%s)**%s'%(species_dict2array[reactant2]+1,stoich))
                                                    step+=1
                                                else:
                                                    f.write('*y(%s)**%s'%(species_dict2array[reactant2]+1,stoich)) 
                                        
                            # Now at the end of the reactant multiplication, introduce reaction rate
                            if step > 1:
                                f.write('*r(%s)'%(equation+1))
                            elif step == 1:
                                f.write('-1.0*r(%s)'%(equation+1))
                            equation_flag+=1
                            
                # Now repeat the above procedure but for reactions that lead to concentration gain    
                if len(equations_list_gain) > 0:
                    for equation in equations_list_gain:
                        reactants_list_gain_temp=[reactant for step, reactant in rate_dict_reactants[equation].items()]
                        if reactant in reactants_list_gain_temp:
                            if equation_flag > 3: # this is just to ensure line length isnt a problem for Fortran compiler
                                f.write(' & \n')
                                f.write('    ')  
                                equation_flag=0
                            stoich_first=loss_dict[reactant][equation]
                            #except:
                            #    pdb.set_trace()
                            #if loss_step=0:
                                #f.write('-')
                            #    loss_step+=
                            step=1
                            if stoich_first==1:
                                for reactant_step, reactant2 in rate_dict_reactants[equation].items():
                                    #step=1
                                    if reactant2 not in ['hv'] and reactant2 != reactant:
                                        stoich=loss_dict[reactant2][equation]
                                        if stoich == 1:
                                            if step == 1:
                                                f.write('+1.0*y(%s)'%(species_dict2array[reactant2]+1))
                                                step+=1
                                            else:
                                                f.write('*y(%s)'%(species_dict2array[reactant2]+1))
                                        else:
                                            if step == 1:
                                                f.write('+1.0*y(%s)**%s'%(species_dict2array[reactant2]+1,stoich))
                                                step+=1
                                            else:
                                                f.write('*y(%s)**%s'%(species_dict2array[reactant2]+1,stoich))
                            else:
                                #f.write('-%s*'%(stoich_first))
                                for reactant_step, reactant2 in rate_dict_reactants[equation].items():
                                    #step=1
                                    if reactant2 not in ['hv']:
                                        if reactant2 == reactant:
                                            if stoich_first == 2:
                                                if step == 1:
                                                    f.write('+%s*y(%s)'%(stoich_first,species_dict2array[reactant2]+1))
                                                    step+=1
                                                else:
                                                    f.write('%s*y(%s)'%(stoich_first,species_dict2array[reactant2]+1))
                                            elif stoich_first > 2:
                                                if step == 1:
                                                    f.write('+%s*y(%s)**%s'%(stoich_first,species_dict2array[reactant2]+1,stoich_first-1))
                                                    step+=1
                                                else:
                                                    f.write('%s*y(%s)**%s'%(stoich_first,species_dict2array[reactant2]+1,stoich_first-1))
                                        elif reactant2 != reactant:
                                            stoich=loss_dict[reactant2][equation]
                                            if stoich == 1:
                                                if step == 1:
                                                    f.write('+1.0*y(%s)'%(species_dict2array[reactant2]+1))
                                                    step+=1
                                                else:
                                                    f.write('*y(%s)'%(species_dict2array[reactant2]+1))
                                            else:
                                                if step == 1:
                                                    f.write('+1.0*y(%s)**%s'%(species_dict2array[reactant2]+1,stoich))
                                                    step+=1
                                                else:
                                                    f.write('*y(%s)**%s'%(species_dict2array[reactant2]+1,stoich)) 
                                        

                            if step > 1:
                                f.write('*r(%s)'%(equation+1))
                            elif step == 1:
                                f.write('+1.0*r(%s)'%(equation+1))
                            equation_flag+=1


                f.write(' \n')
                check_step+=1
    if OpenMP is True:
        f.write('!$OMP END SINGLE NOWAIT \n') 
        f.write('!$OMP END PARALLEL DO \n') 
    f.write('end subroutine \n') 
    f.close()  


def write_partitioning_section_fortran(total_length_y,num_bins,num_species):
    f = open('Partitioning.f90','w')
    f.write('!##################################################################################################### \n') # python will convert \n to os.linesep
    f.write('! Fortran function to calculating dy/dt according to gas-to-particle partitioning                    # \n') # python will convert \n to os.linesep
    f.write('!    Copyright (C) 2018  David Topping : david.topping@manchester.ac.uk                              # \n')
    f.write('!                                      : davetopp80@gmail.com                                        # \n')
    f.write('!    Personal website: davetoppingsci.com                                                            # \n')
    f.write('!                                                                                                    # \n')
    f.write('!    This program does not have a license, meaning the deault copyright law applies.                 # \n')
    f.write('!    Only users who have access to the private repository that holds this file may                   # \n')
    f.write('!    use it or develop it, but may not distribute it.                                                # \n')
    f.write('!                                                                                                    # \n')
    f.write('!                                                                                                    # \n')
    f.write('!##################################################################################################### \n')    
    f.write('! File created at %s \n' % datetime.datetime.now()) # python will convert \n to os.linesep
    f.write('\n') 
    #f.write('import numpy as np\n') 
    f.write('\n')
    f.write('subroutine dydt_partition(y,ycore,core_diss,core_mass_array,& \n') 
    f.write('                          density_input,core_density_array,& \n')
    f.write('                          ignore_index,mw_array,Psat,& \n')
    f.write('                          DStar_org,alpha_d_org,C_g_i_t,N_perbin,& \n')
    f.write('                          gamma_gas,Latent_heat,GRAV,Updraft,& \n')
    f.write('                          sigma,NA,kb,Rv,R_gas,Model_temp,& \n')
    f.write('                          cp,Ra,Lv_water_vapour,size_array,dy_dt) \n')
    f.write('\n') 
    f.write('    implicit none \n')     
    f.write('\n') 
    f.write('    !Define variables in order of appearance. Leave intent(out) until end \n')
    f.write('    REAL(8), intent(in), dimension(%s) :: y \n'%(total_length_y))
    f.write('    !f2py intent(in) :: y \n') 
    f.write('    REAL(8), intent(in), dimension(%s) :: ycore,core_mass_array,& \n'%(num_bins))
    f.write('                         core_density_array, N_perbin \n') 
    f.write('    REAL(8), intent(in), dimension(%s) :: density_input,DStar_org,& \n'%(num_species))
    f.write('                         alpha_d_org, gamma_gas, Latent_heat,C_g_i_t \n') 
    f.write('    INTEGER(8), intent(in), dimension(%s) :: ignore_index,mw_array,Psat \n'%(num_species))
    f.write('    !f2py intent(in) ::  ycore,core_mass_array \n')
    f.write('    !f2py intent(in) ::  core_density_array, N_perbin \n') 
    f.write('    !f2py intent(in) ::  density_input,DStar_org,& \n')
    f.write('    !f2py intent(in) ::  alpha_d_org, gamma_gas, Latent_heat, C_g_i_t \n') 
    f.write('    !f2py intent(in) ::  ignore_index,mw_array \n')
    f.write('    !Global variables - define in parent routine \n')     
    f.write('    REAL(8), intent(in) :: core_diss, GRAV,Updraft,cp,Ra,Lv_water_vapour \n')
    f.write('    REAL(8), intent(in) :: sigma,NA,kb,Rv,R_gas,Model_temp \n')
    f.write('    !f2py intent(in) :: core_diss,GRAV,Updraft,cp,Ra,Lv_water_vapour \n') 
    f.write('    !f2py intent(in) :: sigma,NA,kb,Rv,R_gas,Model_temp \n')     
    f.write('    REAL(8), intent(out), dimension(%s) :: size_array \n'%(num_bins))
    f.write('    REAL(8), intent(out), dimension(%s) :: dy_dt \n'%(total_length_y))
    f.write('    !f2py intent(out) :: size_array,dy_dt \n')     
    f.write('\n') 
    f.write('    !temporary variables used in each size bin calculation \n')     
    f.write('    REAL(8), dimension(%s) :: Inverse_Kn, Correction_part1, &\n'%(num_species)) 
    f.write('                              Correction_part2, Correction_part3,& \n')
    f.write('                              Correction, Kn, dm_dt,& \n')
    f.write('                              k_i_m_t_part1, k_i_m_t \n') 
    f.write('    REAL(8), dimension(%s,%s) :: dy_dt_gas_matrix \n'%(num_species,num_bins))
    f.write('    REAL(8) total_moles, density, diameter, radius, water_moles  \n')
    f.write('    REAL(8), dimension(%s) :: Pressure_eq, Cstar_i_m_t \n'%(num_species)) 
    f.write('    REAL(8), dimension(%s) :: kelvin_factor,y_mole_fractions  \n'%(num_species))
    f.write('    REAL(8), dimension(%s) :: temp_array  \n'%(num_species))
    f.write('    REAL(8), dimension(%s) :: mass_array, density_array, mass_fractions_array \n'%(num_species+1)) 
    f.write('    REAL(8) total_SOA_mass, total_water_cond_mass, total_mass \n')
    f.write('    INTEGER(8) size, index1, index2  \n')
    f.write('\n') 
    f.write('    REAL(8), PARAMETER :: Pi = 3.1415927 \n')   
    f.write('    INTEGER(8), PARAMETER :: num_species = %s \n'%(num_species))       
    f.write('\n')     
    f.write('!$OMP PARALLEL DO &\n') 
    f.write('!$OMP& PRIVATE(size,total_mass),& \n') 
    f.write('!$OMP& PRIVATE(total_moles,y_mole_fractions,temp_array), & \n') 
    f.write('!$OMP& PRIVATE(mass_array,mass_fractions_array,density_array), & \n') 
    f.write('!$OMP& PRIVATE(density,water_moles), & \n') 
    f.write('!$OMP& PRIVATE(Kn,Inverse_Kn,Correction_part1, Correction_part2), & \n') 
    f.write('!$OMP& PRIVATE(Correction_part3,Correction,kelvin_factor), & \n') 
    f.write('!$OMP& PRIVATE(Pressure_eq,Cstar_i_m_t,k_i_m_t_part1), & \n') 
    f.write('!$OMP& PRIVATE(k_i_m_t,dm_dt,water_core_mole_frac), & \n') 
    f.write('!$OMP& SHARED(total_SOA_mass,total_water_cond_mass,size_array), & \n') 
    f.write('!$OMP& SHARED(ycore,core_diss,mw_array,density_input), & \n') 
    f.write('!$OMP& SHARED(C_g_i_t,dy_dt,core_mass_array,R_gas,Psat,num_species), & \n') 
    f.write('!$OMP& SHARED(gamma_gas,NA,alpha_d_org,Model_temp,DStar_org), & \n') 
    f.write('!$OMP& SHARED(dy_dt_gas_matrix) \n') 
    f.write('    DO size=1,%s \n'%(num_bins)) 
    f.write('\n')     	
    f.write('        ! Select a slice of y that represents this size bin \n') 
    f.write('        temp_array(:)=y(num_species+1+((size-1)*num_species):num_species+& \n') 
    f.write('        ((size)*num_species))\n') 
    f.write('        total_moles=SUM(temp_array(:))+ycore(size)*core_diss\n') 
    f.write('        y_mole_fractions(:)=temp_array(:)/total_moles\n') 
    f.write('        ! Now mask the array to set the mole fracs of some species to 0\n') 
    f.write('        WHERE (ignore_index==1) y_mole_fractions = 0.0d0\n') 
    f.write('        mass_array(1:num_species)=(temp_array(:)/NA)*mw_array(:)\n') 
    f.write('        mass_array(num_species+1)=core_mass_array(size)\n') 
    f.write('        density_array(1:num_species)=density_input(1:num_species)\n') 
    f.write('        density_array(num_species+1)=core_density_array(size)\n') 
    f.write('\n')     	
    f.write('        total_SOA_mass=total_SOA_mass+SUM(mass_array(1:num_species-1))\n') 
    f.write('        total_water_cond_mass=total_water_cond_mass+mass_array(num_species)\n') 
    f.write('\n')     	
    f.write('        total_mass=SUM(mass_array(:))\n') 
    f.write('        mass_fractions_array(:)=mass_array(:)/total_mass\n') 
    f.write('\n')     		 
    f.write('        density=1.0D0/(SUM(mass_fractions_array(:)/density_array(:)))\n') 
    f.write('\n')     	
    f.write('        size_array(size)=((3.0D0*((total_mass*1.0D3)/ & \n') 
    f.write('        (N_perbin(size)*1.0D6)))/(4.0D0*Pi*density))**(1.0D0/3.0D0)\n') 
    f.write('        ! F3b) Knudsen number  \n') 
    f.write('        ! Calculate the Knudsen number for all condensing molecules based on this new size \n') 
    f.write('        ! This relies on mean free path for each species [cm] and particle radius [cm]\n') 
    f.write('        Kn(:)=gamma_gas(:)/size_array(size)\n') 
    f.write('\n')     	
    f.write('        ! F3c) Non-continuum regime correction  \n') 
    f.write('        ! Calculate a correction factor according to the continuum versus non-continuum regimes\n') 
    f.write('        ! Expression taken from Jacobson et al (2000), page 457. They reference:\n') 
    f.write('        ! Fuchs and Sutugin 1971\n') 
    f.write('        ! Pruppacher and Klett 1997\n') 
    f.write('        Inverse_Kn(:)=1.0D0/Kn(:)\n') 
    f.write('        Correction_part1(:)=(1.33D0+0.71D0*Inverse_Kn(:))/(1.0D0+Inverse_Kn(:))\n') 
    f.write('        Correction_part2(:)=(4.0D0*(1.0D0-alpha_d_org(:)))/(3.0D0*alpha_d_org(:))\n') 
    f.write('        Correction_part3(:)=1.0E0+(Correction_part2(:)+Correction_part2(:))*Kn(:)\n') 
    f.write('        Correction(:)=1.0D0/Correction_part3(:)\n') 
    f.write('\n')     	
    f.write('        ! F3d) Kelvin factor  \n') 
    f.write('        ! Now calculate a kelvin factor for every semi-volatile compound in this size bin\n') 
    f.write('        kelvin_factor(:)=EXP((4.0E0*mw_array(:)*1.0D-3*sigma)/ & \n')
    f.write('        (R_gas*Model_temp*size_array(size)*2.0D0*density)) \n') 
    f.write('\n')     		
    f.write('        ! F3f) Equilibrium pressure above droplets  \n') 
    f.write('        ! Now calculate an equilibrium RH for every compound in this size bin \n') 
    f.write('        ! This is currently in atmospheres.  \n') 
    f.write('        Pressure_eq(:)=kelvin_factor(:)*y_mole_fractions(:)*Psat(:)*101325.0D0 \n') 
    f.write('\n')     
    f.write('        ! F3g) Calculate the equilibrium concentration equivalent \n') 
    f.write('        Cstar_i_m_t(:)=Pressure_eq(:)*(NA/(8.3144598D6*Model_temp))  \n') 
    f.write('\n')    
    f.write('        ! F3h) The droplet growth equation  \n') 
    f.write('        ! The following calculates the change in mass, per particle of this size. \n') 
    f.write('        ! The equation relies on the following parameters \n')  
    f.write('        ! and units. Always check units \n') 
    f.write('        ! radius [m] \n') 
    f.write('        ! DStar_org [m2/s] - pass in [cm2/s] so needs to be converted by *1.0E-4 \n') 
    f.write('        ! Pressure_gas [Pascals] \n') 
    f.write('        ! Pressure_eq [Pascals] \n') 
    f.write('        ! R_gas [m3 Pascals /(K mol)] \n') 
    f.write('        ! molw [g/mol] \n') 
    f.write('        ! T [K] \n') 
    f.write('        ! The units of the equation should therefore be g/s \n') 
    f.write('        ! ASTEM version - should be the same as the Jacobson version once coverted to molecules/cc/s \n') 
    f.write('        k_i_m_t_part1(:) = DStar_org(:)*Correction(:) \n') 
    f.write('        k_i_m_t(:)=4.0D0*Pi*size_array(size)*1.0D2*N_perbin(size)*k_i_m_t_part1(:) \n')   
    f.write('        dm_dt(:)=k_i_m_t(:)*(C_g_i_t(:)-Cstar_i_m_t(:))   \n') 
    f.write('\n')    	
    f.write('        ! F3i) Now update the contribution to the ODEs being solved  \n') 
    f.write('        ! 1) Add contributory loss from the gas phase to particle phase [this includes water] \n') 
    f.write('        dy_dt_gas_matrix(1:num_species,size)=dy_dt(1:num_species)-dm_dt(1:num_species) \n') 
    f.write('        ! 2) Add a contributory gain to the particle phase from the gas phase \n') 
    f.write('        index1=num_species+1+((size-1)*num_species)  \n') 
    f.write('        index2=num_species+((size)*num_species)  \n') 
    f.write('        dy_dt(index1:index2)=dm_dt(1:num_species)  \n') 
    f.write('\n')    		
    f.write('    ENDDO \n') 
    f.write('!$OMP END PARALLEL DO \n') 
    f.write('   dy_dt(1:num_species)=SUM(dy_dt_gas_matrix, DIM=2)  \n')
    f.write('\n')    		    
    f.write('end subroutine \n') 
    f.close()  
	


