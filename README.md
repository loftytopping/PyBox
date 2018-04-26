# PyBox

This repository will hold model variants that ultimately couple gas/aerosol chemistry, gas-to-particle partitioning and varying phase state through the full humidity cycle in the atmosphere. The first phase of the project is to develop and profile a gas phase only model, using the [Master Chemical Mechanism (MCM)](http://mcm.leeds.ac.uk/MCM/) as the basis. The model will also relate component properties, using molecular structural information, through the [UManSysProp](http://umansysprop.seaes.manchester.ac.uk) informatics suite.  Any public release will occur according to agreement from any partner contributions and associated papers.

Please check the project wiki page for more inforamtion on updates and planned releases

The model works on the basis of reading a file that defines reactions in the gas phase. Within this repo, some examples are given as taken from the MCM. For example, take the mixed VOC file given by 'MCM_mixed_test.eqn.txt'. This contains the following:

##### {46.} 	 CH3OH + OH = HO2 + HCHO : 	2.85D-12*EXP(-345/TEMP) 	;
##### {47.} 	 C2H5OH + OH = C2H5O : 	3.0D-12*EXP(20/TEMP)*0.05 	;
##### {48.} 	 C2H5OH + OH = CH3CHO + HO2 : 	3.0D-12*EXP(20/TEMP)*0.9 	;
##### {49.} 	 C2H5OH + OH = HOCH2CH2O2 : 	3.0D-12*EXP(20/TEMP)*0.05 	;

Where the equation number is first defined, then the reactants/products along with a defined rate coefficient. Some reactions rely on coefficients defined elsewhere, according to the MCM version number. This file is first parsed in the file 'Parse_eqn_file.py', providing information that can be used to set up and solve the relevant ODEs. 

The model relies on the [Assimulo](http://www.jmodelica.org/assimulo) ODE solver package, and has been built/tested within the Python3 Anaconda environment. There are two versions provided.

## Python only
This is the default version in the main folder. Using the [Numba](https://numba.pydata.org) package, the set of functions that define the ODEs being solved are compiled before the first simulation. This allows an improvement in computational speed over pure python functions. I do have the option to use standard Numpy libraries and sparse matrices, but please check the wiki for news on this release for any educational/training purposes. You will therefore find the first simulation will take some time to compile the relevant libraries, but once compiled will provide a benefit. To run the model, simply run:
#### python Gas_simulation.py
from the command line. Please check the ambient conditions and species concentrations within both 'Gas_simulation.py' and 'ODE_solver.py'. In the former, you can define the starting conditions for species using names explicitly defined in the above equation file. For example, the default option is provided as:

    # Define initial concentrations, in pbb, of species using names from KPP file
    species_initial_conc=dict()
    species_initial_conc['O3']=18.0
    species_initial_conc['APINENE']=30.0
    species_initial_conc['BCARY']=20.0

## Python + Fortran


