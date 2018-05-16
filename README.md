# PyBox

This repository holds model variants that couple gas-phase chemistry, gas-to-particle partitioning and varying phase state through the full humidity cycle in the atmosphere. The first phase of the project is to develop and profile a gas phase model, using the [Master Chemical Mechanism (MCM)](http://mcm.leeds.ac.uk/MCM/) as the basis, with an idealised sectional aerosol model. The model will also relate component properties, using molecular structural information, through the [UManSysProp](http://umansysprop.seaes.manchester.ac.uk) informatics suite.  Any public release beyond the base variant will occur according to new processes added, agreement from any partner contributions and/or associated peer-review papers.

Please check the project wiki page for more inforamtion on updates and planned releases.

## Model overview
============================

The model works on the basis of reading a file that defines reactions in the gas phase. For those familiar with using the [Kinetic PreProcessor (KPP) software](http://people.cs.vt.edu/~asandu/Software/Kpp/), this file defined the reactants and products for each reaction within a chemical mechanism along with an associated rate coefficient. For example, take the mixed VOC file given by 'MCM_mixed_test.eqn.txt' in the root directory of PyBox. This contains the following snippet of text:

##### {46.} 	 CH3OH + OH = HO2 + HCHO : 	2.85D-12*EXP(-345/TEMP) 	;
##### {47.} 	 C2H5OH + OH = C2H5O : 	3.0D-12*EXP(20/TEMP)*0.05 	;
##### {48.} 	 C2H5OH + OH = CH3CHO + HO2 : 	3.0D-12*EXP(20/TEMP)*0.9 	;
##### {49.} 	 C2H5OH + OH = HOCH2CH2O2 : 	3.0D-12*EXP(20/TEMP)*0.05 	;

Where the equation number is defined first, then the reactants/products along with a defined rate coefficient. Some reactions rely on coefficients defined elsewhere, according to the MCM version number. These are also included in PyBox. This equation file is first parsed using the file 'Parse_eqn_file.py', providing information that can be used to set up and solve the relevant ordinary differential equations (ODEs) to simulate the entire chemical mechanism.  Each component, or specie, in this chemical mechanism also has an associated record of chemical structure in the form of a [SMILES string](http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html). This information is carried in a .xml file, provided by the MCM, and stored in the root directory of PyBox. Why is this important? Well, this information is taken by the [UManSysProp](http://umansysprop.seaes.manchester.ac.uk) informatics suite and allows us to predict properties of each compound that helps us predict whether they are likely to remain in the gas phase or condense to an existing particulate phase through gas-to-particle partitioning. Before we take a look at the directory structure provided in this repository, lets deal with the dependencies.

## Dependencies
============================

PyBox has been, and is continually, built in the [Anaconda Python 3.6 environment](https://www.anaconda.com/download/#macos). This allows me to use the Numpy and Scipy modules contained within it amongst others listed below. However, even if you do use the Anaconda distribution, PyBox required specific packages:

- [Assimulo](http://www.jmodelica.org/assimulo). This is the numerical core of PyBox. The Assimulo ODE solver package allows us to use multiple ODE solvers designed for stiff systems, including the Rosenbrock method. As found on the project website, there are multiple [methods for installation](https://jmodelica.org/assimulo/installation.html) from both package managers to compiling from source.  From my own experience, it is better to build from source against the Anaconda Python environment. You will need to point to the location of the [Sundials solver suite](https://computation.llnl.gov/projects/sundials) and both BLAS and LAPACK. You can check if your Assimulo installation has worked by opening an interactive Python shell and typing:

> from assimulo.solvers import RodasODE, CVode

to test import of both the Rosenbrock and CVode ODE method.

- [UManSysProp](http://umansysprop.seaes.manchester.ac.uk). As described on the UManSysProp project page, this model was developed at the University of Manchester under a research grant in order to automate predictions of pure component and mixture properties. This suite requires the Python interface to the [OpenBabel](https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html) package and uses [Flask WTF](https://flask-wtf.readthedocs.io/en/stable/). You can clone the suite from the [project github page](https://github.com/loftytopping/UManSysProp_public). Once you have cloned the repository, you will need to add the location of the suite in the python script 'Property_calculation.py' within the 'Aerosol' directory of the PyBox repository. As with the Assimulo package, you can test this import by opening an interactive Python shell and typing:

> import sys

> sys.path.append('<-- add your path here -->/UManSysProp_public/')

> from umansysprop import boiling_points

> from umansysprop import vapour_pressures

> from umansysprop import critical_properties

> from umansysprop import liquid_densities

> from umansysprop import partition_models

> from umansysprop.forms import CoreAbundanceField

Other dependecies used in the [Anaconda Python 3.6 environment](https://www.anaconda.com/download/#macos), or now included in existing Python packages, include:

-[f2py](https://docs.scipy.org/doc/numpy-1.13.0/f2py/index.html), the Fortran to Python Interface Generator which is now included in the Scipy distribution.

-[Numba](https://numba.pydata.org) installed through the conda package manager to generate optimized machine code using the LLVM compiler infrastructure at import time, runtime, or statically.

-[GCC with support for OpenMP](https://gcc.gnu.org/wiki/openmp) if you would like to exploit multicore capabilities of your system in the Python+Fortran model variants described below.

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
    
As noted above, everything is driven by a gas-phase chemical mechanism. It therefore makes sense to have this model contruction within the root directory before we consider gas-to-particle partitioning. Currently there are two variants provided to solve the gas phase model:
   
#### 1) Python [using Numba]
This is the default version in the root directory. Recall the parsing of the equation file? After this, using the [Numba](https://numba.pydata.org) package, the set of functions that define the ODEs being solved are written as new Python files and then compiled before the first simulation. Numba does this as the modules are imported. This allows an improvement in computational speed over pure python functions. You will therefore find the first simulation will take some time to compile the relevant libraries, but once compiled will provide a benefit. Indeed, if you retain the specifica chemical mechanism, Numba will not need to re-compile even when you start with new initial conditions. The current version provides you with an example. Specifically, it is based on the MCM representation of the degredation of [Alpha-Pinene](https://en.wikipedia.org/wiki/Alpha-Pinene). To run the model, once you are happy all dependecies are installed, type the following from the root directory:

> python Gas_simulation.py

You can modify the ambient conditions and species concentrations within both 'Gas_simulation.py' and 'ODE_solver.py'. In the former, you can define the starting conditions for species using names explicitly defined in the above equation file. For example, the default option is provided as:

    # Define initial concentrations, in pbb, of species using names from KPP file
    species_initial_conc=dict()
    species_initial_conc['O3']=18.0
    species_initial_conc['APINENE']=30.0

#### 2) Python + Fortran [using f2py Fortran to Python Interface Generator] 
Whilst the above variant uses the Numba package, in the folder 'f2py' the same model is constructed using the [F2Py](https://docs.scipy.org/doc/numpy/f2py/)package, where functions that define the ODEs are converted into pre-compiled Fortran modules with the option to use [OpenMP](http://www.openmp.org) to exploit the number of cores available to you on any given platform. As before, please check the relevant files for defining initial conditions, species concetrations, and expect some compilation time during the first run. To run this simulation, type the following from the f2py directory:

> python Gas_simulation_f2py.py


## Automated unit tests 
============================

Within the folder tests, run the following command:

> python test_modules.py -v

This will use the unittest Python module to test the output of generated functions used within the ODE solvers against pre-generated .npy files provided in the data subfolder of tests.

