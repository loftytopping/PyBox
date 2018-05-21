##########################################################################################
#                                                                                        #
#    Scripts to convert a text string definition of a rate coefficient and convert to    #
#    a format that can be used by a Python and/or Fortran environment                    #
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
##########################################################################################

# This script takes a chemical equation file, following the standard KPP format, and generates 
# information used to create an ODE instance to solver for set of specific conditions. 
# The rules used to generate this instance have been based on standard KPP examples. 
# The left hand / right hand descriptions are relatively easy. 
# The definition of a rate coefficient however might change with time.
# For this reason, this file is likely to change until a generic variant can be created that
# covers all expected formats. This is dealt with in a seperate file

import pdb

                
def convert_rate_mcm(rate_dict):

    """ This function takes the defintions of rate coefficients and converts to Python command

    inputs:
    • rate_dict - parsed string representations of rate coefficients
    outputs:
    • rate_dict - converted defintions for use in Python
  
    """
	
    # KPP formats for a specific mechanism also tend to use pre-requisite constants. 
    # For example, the MCM often uses formats for rate coefficients as follows:
    # K234*J(0)*J(3)
    # J(4)+J(5)*EXP(234/TEMP)
    # where K234, J(2) etc is a given constant 
    # In this case the string is saved, but conversions for floats and 'exp' expressions
    # carried out and the expression then evaluated as a whole. For example	
    # J(4)*2.34D-4*EXP(234/TEMP)
    # is converted to
    # J(4)*2.34E-4*numpy.exp(234.0/TEMP)
    # which is then processed as an entire command where needed. Specific variables are saved
    # in a seperate file, updated as the MCM updates

    #Store conversions in lists
    math_list=[
        ('dabs','numpy.abs'),
        ('dsqrt','numpy.sqrt'),
        ('dlog','numpy.log'),
        ('LOG','numpy.log'),
        ('EXP','numpy.exp')
        ]
    map_list=[('D-','E-'),
        ('D+','E+')
        ]

    afterText = ['numpy.abs','numpy.sqrt','numpy.log','numpy.exp']
    new_rate_full=[]
    for equation_step, rate_full in rate_dict.items():
        
        for operation in math_list:
            temp_text=rate_full.replace(operation[0],operation[1])
        for syntax in map_list:
            temp_text=temp_text.replace(syntax[0],syntax[1])

        #Replace all () with [] in the first instance
        temp_text=temp_text.replace('(','[').replace(')',']')

        #Now deal with replacing generic () with [] unless it is a numpy command
        #Find the first occurance of [ and replace back to ) if numpy command
        searchText1 = '['
        searchText2 = ']'
        temp_text_list=list(temp_text)
        for substr in afterText:
            try:
                after_index = temp_text.index(substr)
                index1=temp_text.find(searchText1, after_index)
                index2=temp_text.find(searchText2, after_index)
                temp_text_list[index1]='('
                temp_text_list[index2]=')'
            except:
                pass
            temp_text=''.join(temp_text_list)
    
        #new_rate_full.append(temp_text)
        #if print_options==1:
        #    print temp_text

        rate_dict[equation_step]=temp_text

    return rate_dict

def convert_rate_mcm_numba(rate_dict):

    """ This function takes the defintions of rate coefficients and converts to Python [Numba] command

    inputs:
    • rate_dict - parsed string representations of rate coefficients
    outputs:
    • rate_dict - converted defintions for use in Python
  
    """

    # KPP formats for a specific mechanism also tend to use pre-requisite constants. 
    # For example, the MCM often uses formats for rate coefficients as follows:
    # K234*J(0)*J(3)
    # J(4)+J(5)*EXP(234/TEMP)
    # where K234, J(2) etc is a given constant 
    # In this case the string is saved, but conversions for floats and 'exp' expressions
    # carried out and the expression then evaluated as a whole. For example	
    # J(4)*2.34D-4*EXP(234/TEMP)
    # is converted to
    # J(4)*2.34E-4*numpy.exp(234.0/TEMP)
    # which is then processed as an entire command where needed. Specific variables are saved
    # in a seperate file, updated as the MCM updates

    #Store conversions in lists
    math_list=[
        ('dabs','numba_abs'),
        ('dsqrt','numba_sqrt'),
        ('dlog','numba_log'),
        ('LOG','numba_log'),
        ('EXP','numba_exp')
        ]
    map_list=[('D-','E-'),
        ('D+','E+'),
        ('D1','E+1'),
        ('D2','E+2'),
        ('D3','E+3'),
        ('D4','E+4'),
        ('D5','E+5'),
        ('D6','E+6'),
        ('D7','E+7'),
        ('D8','E+8'),
        ('D9','E+9'),
        ]

    afterText = ['numba_abs','numba_sqrt','numba_log','numba_exp']
    new_rate_full=[]
    for equation_step, rate_full in rate_dict.items():
        
        for operation in math_list:
            temp_text=rate_full.replace(operation[0],operation[1])
        for syntax in map_list:
            temp_text=temp_text.replace(syntax[0],syntax[1])

        temp_text_list=list(temp_text)
        indices= []
        index = -1
        while True:
            index = temp_text.find('J(', index + 1)
            if index == -1:  
                break  # All occurrences have been found
            indices.append(index)
        #Replace () with [] if associated with J for photolysis rates
        try:
            for num in indices:
                #after_index=temp_text.index('J(')
                #pdb.set_trace()
                index1=temp_text.find('(', num)
                index2=temp_text.find(')', num)
                temp_text_list[index1]='['
                temp_text_list[index2]=']'
        except:
            pass
        temp_text=''.join(temp_text_list)

        ##Replace all () with [] in the first instance
        #temp_text=temp_text.replace('(','[').replace(')',']')

        ##Now deal with replacing generic () with [] unless it is a numpy command
        ##Find the first occurance of [ and replace back to ) if numpy command
        #searchText1 = '['
        #searchText2 = ']'
        #temp_text_list=list(temp_text)
        #for substr in afterText:
        #    try:
        #        after_index = temp_text.index(substr)
        #        index1=temp_text.find(searchText1, after_index)
        #        index2=temp_text.find(searchText2, after_index)
        #        temp_text_list[index1]='('
        #        temp_text_list[index2]=')'
        #    except:
        #        pass
        #    temp_text=''.join(temp_text_list)
    
        #new_rate_full.append(temp_text)
        #if print_options==1:
        #    print temp_text

        rate_dict[equation_step]=temp_text

    return rate_dict


def convert_rate_mcm_fortran(rate_dict):

    """ This function takes the defintions of rate coefficients and converts to Fortran command

    inputs:
    • rate_dict - parsed string representations of rate coefficients
    outputs:
    • rate_dict - converted defintions for use in Fortran
  
    """


    # KPP formats for a specific mechanism also tend to use pre-requisite constants. 
    # For example, the MCM often uses formats for rate coefficients as follows:
    # K234*J(0)*J(3)
    # J(4)+J(5)*EXP(234/TEMP)
    # where K234, J(2) etc is a given constant 
    # In this case the string is saved, but conversions for floats and 'exp' expressions
    # carried out and the expression then evaluated as a whole. For example	
    # J(4)*2.34D-4*EXP(234/TEMP)
    # is converted to
    # J(4)*2.34E-4*numpy.exp(234.0/TEMP)
    # which is then processed as an entire command where needed. Specific variables are saved
    # in a seperate file, updated as the MCM updates

    #Store conversions in lists
    math_list=[
        ('dabs','ABS'),
        ('dsqrt','SQRT'),
        ('dlog','LOG')
        ]
    map_list=[('D-','E-'),
        ('D+','E+'),
        ('D1','E+1'),
        ('D2','E+2'),
        ('D3','E+3'),
        ('D4','E+4'),
        ('D5','E+5'),
        ('D6','E+6'),
        ('D7','E+7'),
        ('D8','E+8'),
        ('D9','E+9'),
        ]
        
    afterText = ['ABS','SQRT','LOG','EXP']

    new_rate_full=[]
    for equation_step, rate_full in rate_dict.items():
        
        for operation in math_list:
            temp_text=rate_full.replace(operation[0],operation[1])
        for syntax in map_list:
            temp_text=temp_text.replace(syntax[0],syntax[1])

        temp_text_list=list(temp_text)
        #indices= []
        #index = -1
        #while True:
        #    index = temp_text.find('J(', index + 1)
        #    if index == -1:  
        #        break  # All occurrences have been found
        #    indices.append(index)
        ##Replace () with [] if associated with J for photolysis rates
        #try:
        #    for num in indices:
        #        #after_index=temp_text.index('J(')
        #        #pdb.set_trace()
        #        index1=temp_text.find('(', num)
        #        index2=temp_text.find(')', num)
        #        temp_text_list[index1]='['
        #        temp_text_list[index2]=']'
        #except:
        #    pass
        temp_text=''.join(temp_text_list)
        
    #new_rate_full=[]
    #for equation_step, rate_full in rate_dict.items():
    #    
    #    for operation in math_list:
    #        temp_text=rate_full.replace(operation[0],operation[1])
    #    for syntax in map_list:
    #        temp_text=temp_text.replace(syntax[0],syntax[1])#

        #Replace all () with [] in the first instance
    #    temp_text=temp_text.replace('(','[').replace(')',']')

        #Now deal with replacing generic () with [] unless it is a numpy command
        #Find the first occurance of [ and replace back to ) if numpy command
    #    searchText1 = '['
    #    searchText2 = ']'
    #    temp_text_list=list(temp_text)
    #    for substr in afterText:
    #        try:
    #            after_index = temp_text.index(substr)
    #            index1=temp_text.find(searchText1, after_index)
    #            index2=temp_text.find(searchText2, after_index)
    #            temp_text_list[index1]='('
    #            temp_text_list[index2]=')'
    #        except:
    #            pass
    #        temp_text=''.join(temp_text_list)
    
        #new_rate_full.append(temp_text)
        #if print_options==1:
        #    print temp_text

        rate_dict[equation_step]=temp_text

    return rate_dict






