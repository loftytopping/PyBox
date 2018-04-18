# This script takes a chemical equation file, following the standard KPP format, and generates 
# information used to create an ODE instance to solver for set of specific conditions. 
# The rules used to generate this instance have been based on standard KPP examples. 
# The left hand / right hand descriptions are relatively easy. 
# The definition of a rate coefficient however might change with time.
# For this reason, this file is likely to change until a generic variant can be created that
# covers all expected formats. This is dealt with in a seperate file

def convert_rate_general(rate_dict,print_options): #rate_dict contains a list of rate constant text strings
                                           #as extracted by the Parse_eqn_file.py script

    # Create dictionaries used when calculating final rates in ODE solvers
    rate_dict_converted=collections.defaultdict(lambda: collections.defaultdict())

    for equation_step, rate_full in rate_dict.iteritems():
        print("Extracting equation :",equation_step)
        print("Given rate string of:",rate_full)

        #1) Check for standard rate coefficient definitions
        if 'ARR' in rate_full:
            rate_def[equation_step]='ARR'
            numbers=rate_full.split('(',1)[1].split(')',1)[0].split(',') #take out numbers in the brackets
            rate_dict[equation_step]['a']=float(numbers[0])
            rate_dict[equation_step]['b']=float(numbers[1])
            rate_dict[equation_step]['c']=float(numbers[2])
        elif 'FALL' in rate_full:
            rate_def[equation_step]='FALL'
            numbers=rate_full.split('(',1)[1].split(')',1)[0].split(',')
            rate_dict[equation_step]['a']=float(numbers[0])
            rate_dict[equation_step]['b']=float(numbers[1])
            rate_dict[equation_step]['c']=float(numbers[2])
            rate_dict[equation_step]['d']=float(numbers[3])
            rate_dict[equation_step]['e']=float(numbers[4])
            rate_dict[equation_step]['f']=float(numbers[5])
            rate_dict[equation_step]['g']=float(numbers[6])
        elif 'EP2' in rate_full:
            rate_def[equation_step]='EP2'
            numbers=rate_full.split('(',1)[1].split(')',1)[0].split(',')
            rate_dict[equation_step]['a']=float(numbers[0])
            rate_dict[equation_step]['b']=float(numbers[1])
            rate_dict[equation_step]['c']=float(numbers[2])
            rate_dict[equation_step]['d']=float(numbers[3])
            rate_dict[equation_step]['e']=float(numbers[4])
            rate_dict[equation_step]['f']=float(numbers[5])
        elif 'EP3' in rate_full:
            rate_def[equation_step]='EP3'
            numbers=rate_full.split('(',1)[1].split(')',1)[0].split(',')
            rate_dict[equation_step]['a']=float(numbers[0])
            rate_dict[equation_step]['b']=float(numbers[1])
            rate_dict[equation_step]['c']=float(numbers[2])
            rate_dict[equation_step]['d']=float(numbers[3])
        elif 'SUN' in rate_full:
            rate_def[equation_step]='SUN'
            rate_dict[equation_step]['a']=float(rate_full.split('*',1)[0])
            rate_dict[equation_step]['b']=float(rate_full.split('/',1)[1].split(')',1)[0])
        else:
            rate_def[equation_step]='constant'
            try:
                rate_dict[equation_step]['a']=float(rate_full.split('(',1)[1].split(')',1)[0])
            except: 
                rate_dict[equation_step]['a']=float(rate_full)
                
def convert_rate_mcm(rate_dict):

    #2) KPP formats for a specific mechanism also tend to use pre-requisite constants. 
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

    #2) KPP formats for a specific mechanism also tend to use pre-requisite constants. 
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
        ('D+','E+')
        ]

    afterText = ['numba_abs','numba_sqrt','numba_log','numba_exp']
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


def convert_rate_mcm_fortran(rate_dict):

    #2) KPP formats for a specific mechanism also tend to use pre-requisite constants. 
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
        ('D+','E+')
        ]

    afterText = ['ABS','SQRT','LOG','EXP']
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






