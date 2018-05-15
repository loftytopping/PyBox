# PyBox

This repository holds model variants that couple gas/aerosol chemistry, gas-to-particle partitioning and varying phase state through the full humidity cycle in the atmosphere. The first phase of the project is to develop and profile a gas phase only model, using the [Master Chemical Mechanism (MCM)](http://mcm.leeds.ac.uk/MCM/) as the basis. The model will also relate component properties, using molecular structural information, through the [UManSysProp](http://umansysprop.seaes.manchester.ac.uk) informatics suite.  Any public release will occur according to agreement from any partner contributions and associated papers.

Please check the project wiki page for more inforamtion on updates and planned releases.

## Model overview
============================

The model works on the basis of reading a file that defines reactions in the gas phase. Within this repository, examples are taken from the MCM. For example, take the mixed VOC file given by 'MCM_mixed_test.eqn.txt'. This contains the following snippet:

##### {46.} 	 CH3OH + OH = HO2 + HCHO : 	2.85D-12*EXP(-345/TEMP) 	;
##### {47.} 	 C2H5OH + OH = C2H5O : 	3.0D-12*EXP(20/TEMP)*0.05 	;
##### {48.} 	 C2H5OH + OH = CH3CHO + HO2 : 	3.0D-12*EXP(20/TEMP)*0.9 	;
##### {49.} 	 C2H5OH + OH = HOCH2CH2O2 : 	3.0D-12*EXP(20/TEMP)*0.05 	;

Where the equation number is first defined, then the reactants/products along with a defined rate coefficient. Some reactions rely on coefficients defined elsewhere, according to the MCM version number. This file is first parsed using the file 'Parse_eqn_file.py', providing information that can be used to set up and solve the relevant ordinary differential equations (ODEs). 

## Folder Structure 
============================

### Directory layout

    .                           # Gas phase only model [using Numba]
    ├── f2py                    # Gas phase only model [using f2py Fortran to Python Interface Generator] 
    ├── Aerosol                 # Coupled gas and gas-to-particle partitioning routines
    |------f2py                 # Coupled gas and gas-to-particle partitioning routines [using f2py]
    ├── test                    # Automated unit tests
    |------data                 # Data used in the automated unit tests
    ├── LICENSE
    └── README.md
   
### Python only variants [using Numba]
This is the default version in the main folder. Using the [Numba](https://numba.pydata.org) package, the set of functions that define the ODEs being solved are compiled before the first simulation. This allows an improvement in computational speed over pure python functions. I do have the option to use standard Numpy libraries and sparse matrices, but please check the wiki for news on this release for any educational/training purposes. You will therefore find the first simulation will take some time to compile the relevant libraries, but once compiled will provide a benefit. To run the model, simply run:

> python Gas_simulation.py

from the command line. Please check the ambient conditions and species concentrations within both 'Gas_simulation.py' and 'ODE_solver.py'. In the former, you can define the starting conditions for species using names explicitly defined in the above equation file. For example, the default option is provided as:

    # Define initial concentrations, in pbb, of species using names from KPP file
    species_initial_conc=dict()
    species_initial_conc['O3']=18.0
    species_initial_conc['APINENE']=30.0
    species_initial_conc['BCARY']=20.0

### Python + Fortran [using f2py Fortran to Python Interface Generator] 
Whilst the above variant uses the Numba package, in the folder 'f2py' the same model is constructed using the [F2Py](https://docs.scipy.org/doc/numpy/f2py/)package, where functions that define the ODEs are converted into pre-compiled Fortran modules with the option to use [OpenMP](http://www.openmp.org) to exploit the number of cores available to you on any given platform. As before, please check the relevant files for defining initial conditions, species concetrations, and expect some compilation time during the first run.


## Dependencies
============================

-[Assimulo](http://www.jmodelica.org/assimulo). The model currently relies on the Assimulo ODE solver package.  This allows us to use multiple ODE solvers designed for stiff systems. See the project website for installation instruction.

-[UManSysProp](http://umansysprop.seaes.manchester.ac.uk). As described on the UManSysProp project page, this model was developed at the University of Manchester in order to automate predictions of pure component and mixture properties. This suite requires the Python interface to the [OpenBabel](https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html) package.

Other dependecies used in the [Anaconda Python 3.6 environment](https://www.anaconda.com/download/#macos), or now included in existing Python packages, are:

-[f2py](https://docs.scipy.org/doc/numpy-1.13.0/f2py/index.html), the Fortran to Python Interface Generator which is now included in the Scipy distribution.

-[Numba](https://numba.pydata.org) installed through the conda package manager to generate optimized machine code using the LLVM compiler infrastructure at import time, runtime, or statically.

-[GCC with support for OpenMP](https://gcc.gnu.org/wiki/openmp) if you would like to exploit multicore capabilities of your system in the Python+Fortran model variants described below.

## Automated unit tests 
============================

Within the folder tests, run the following command:

> python test_modules.py -v

This will use the unittest Python module to test the output of generated functions used within the ODE solvers against pre-generated .npy files provided in the data subfolder of tests.

