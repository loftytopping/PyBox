# PyBox  [![DOI](http://joss.theoj.org/papers/10.21105/joss.00755/status.svg)](https://doi.org/10.21105/joss.00755)

PyBox is a Python based box-model generator and simulator designed for atmospheric chemistry and aerosol studies. The first phase of the PyBox project is to develop a gas phase model, using the reaction information within the [Master Chemical Mechanism (MCM)](http://mcm.leeds.ac.uk/MCM/) as the basis, coupled with an idealised sectional aerosol model. PyBox also relates component properties, using molecular structural information, through the [UManSysProp](http://umansysprop.seaes.manchester.ac.uk) informatics suite.  Any public release will occur according to new processes added, agreement from any partner contributions and/or associated peer-review papers.

Please check the project wiki page for more information on updates and planned releases. You can also check the current PyBox documentation at https://pybox.readthedocs.io/en/latest/

This project is licensed under the terms of the GNU General Public License v3.0, as provided with this repository. 

# Table of contents
1. [Model overview](#Model-overview)
2. [Dependencies and installation](#Dependencies)
     * 2a. [Using your own machine](#own)
     * 2b. [Binder](#Binder)
     * 2c. [Docker](#Docker)
3. [Folder structure and running the model](#Folder-Structure)
4. [Unit tests](#Automated-unit-tests)
5. [Contributing](#Contributing)
6. [Code of Conduct](#Code-of-Conduct)
7. [Citation](#Citation)

## Model overview<a name="Model-overview"></a>

PyBox works on the basis of reading a file that defines reactions between compounds in the gas phase and the associated reaction coefficient. For example, take the MCM [Alpha-Pinene](https://en.wikipedia.org/wiki/Alpha-Pinene) chemical mechanism file  'MCM_APINENE.eqn.txt' stored in the 'mechanism_files' directory of PyBox. This contains the following snippet of text:

     {125.} 	 C96OOH + OH = C96O2 : 	1.90D-12*EXP(190/TEMP) 	;
     {126.} 	 C96OOH + OH = NORPINAL + OH : 	1.30D-11 	;
     {127.} 	 C96OOH = C96O + OH : 	J(41)+J(22) 	;
     {128.} 	 C96NO3 + OH = NORPINAL + NO2 : 	2.88D-12 	;
     {129.} 	 C96NO3 = C96O + NO2 : 	J(53)+J(22) 	;
     {130.} 	 C96O = C97O2 : 	4.20D+10*EXP(-3523/TEMP) 	;
     {131.} 	 C96OH + OH = NORPINAL + HO2 : 	7.67D-12 	;
     {132.} 	 C96OH = C96O + HO2 : 	J(22) 	;

Where the equation number is defined first, then the reactants/products along with a defined rate coefficient. This equation file is parsed by functions in 'Parse_eqn_file.py', providing information that can be used to set up and solve the relevant ordinary differential equations (ODEs) to simulate the evolution of the chemical mechanism.  Each component in this chemical mechanism also has an associated record of chemical structure in the form of a [SMILES string](http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html). This information is carried in a .xml file, provided by the MCM, and stored in the root directory of PyBox. Why is this important? Well, this information is taken by the [UManSysProp](http://umansysprop.seaes.manchester.ac.uk) informatics suite and allows us to predict properties of each compound that helps us predict whether they are likely to remain in the gas phase or condense to an existing particulate phase through gas-to-particle partitioning. Before we take a look at the directory structure provided in this repository, lets deal with the dependencies.

## Dependencies and installation <a name="Dependencies"></a>

PyBox has been built in the [Anaconda environment](https://www.anaconda.com/download/#macos). [Assimulo](http://www.jmodelica.org/assimulo) is *currently* the numerical core of PyBox. The Assimulo Ordinary Differential Equation (ODE) solver package allows us to use solvers designed for stiff systems. [UManSysProp](http://umansysprop.seaes.manchester.ac.uk) is used to automate predictions of pure component and mixture properties to allow gas-to-particle partitioning simulations *if you need to run the partitioning option*. 
 
### 2(a). Using your own machine<a name="own"></a>

Having a Python distribution on your own machine is attractive for a number of reasons, not least gaining familiarity with building projects in your own time. If you havent already, I would reccomend installing the Anaconda distribution. You can download a copy using [this link](https://www.anaconda.com/products/individual). That page will give you the option to download a version for Windows, Mac or Linux. Download the graphical installer and, typically, accept all options. Once you have installed this, now open a terminal. On Windows, go to the menu of options and find 'Anaconda Prompt' under the Anaconda folder. On a Mac, go to Finder -> Utilities -> Terminal. If on a Mac, when in this terminal when you type:

> Python

Do you see the reference to Anaconda? For example, you may see something *similar* to:

> Python 3.7.6 (default, Jan  8 2020, 20:23:39) [MSC v.1916 64 bit (AMD64)] :: Anaconda, Inc. on win32

Now we are going to create a virtual environment to run our notebooks in. Virtal environments are a great way of maintaining a 'work space' that is seperate to your default installation. For example, if you are going to start installing lots of bespoke modules, you may sometimes come across a clash of version numbers which then becomes tricky to maintain. In the worst case scenario, this would require a re-installation of Python. So lets create a virtual environment for our project. You can switch-on and switch-ff these virtual environments from the command line/terminal whenever you need them.

If you are on Windows, go back to the Anaconda prompt. If you are on a Mac or Linux, go back to ther terminal. First we need to clone this repository. We should use Git for this, becuase with Git you can keep pulling updates from this repository. If you do not already have Git installed on your machine, you can get it from the [download page](https://git-scm.com/downloads). Once you have installed this, at the prompt/terminal type:

> git clone https://github.com/loftytopping/PyBox.git

This will download the project to the location you are in already. You can change this location before running the above command, or move the folder later. Github also gives you the option to download a ZIP file of the entire project if you cannot or do not want to use Git. Once you have the project downloaded, open a command propmt/terminal and navigate to the project folder. We are now going to use the file 'environment.yml' to create a new virtual environment. Run the following command:

> conda env create -f environment.yml

You will see a number of packages being downloaded, eventually, by the conda package manager which is part of the Anaconda distribution. Accept any requests and, when finished, you will see a message that resembles the following:

    To activate this environment, use
    
        $ conda activate PyBox
    
    To deactivate an active environment, use
    
        $ conda deactivate
        
These are the commands for swithing on/off this new virtual environment. Let's switch it on. Type the following in the command prompt/terminal:

> conda activate PyBox

In the command prompt, you will see the name (PyBox) appear from (base). Now we can run the default gas phase example. Still within the project folder, type the following:

> python Gas_simulation.py 

### 2(b). Binder <a name="Binder"></a>

If you do not, or cannot, run Python from your own machine we have provided the ability for you to interact with these files using Binder. The Binder project offers an easy place to share computing environments to everyone. It allows users to specify custom environments and share them with a [single link](https://jupyter.org/binder). Indeed, if you click the link below this will spin-up an individual session for you. Please bare in mind it can take a while to start, and if idle for a short period these sessions will stop. However you can download your notebook file during the session. Everytime you start a Binder link, it will start from scratch.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/loftytopping/Env_modelling/master)


### 2(c). Using a Docker container <a name="Docker"></a>

As a method for a 'fully automated setup' I have provided the option to setup and run PyBox within a Docker container. I have provided a Dockerfile that will automatically build all dependencies within a new container based on the Ubuntu:16.04 image. To build the new image, assuming you have Docker installed, run the following command in the directory of the supplied Dockerfile:

> docker build -t pybox .

After this has completed [which may take some time], type the following to see your new image listed:

> docker images

To create and run a new container based on this image, with a name 'project_pybox', type:

> docker run --name=project_pybox -it pybox

This will take you in to the container. So, lets run the gas phase model in PyBox whilst you are there. Change directory to where PyBox is located:

> cd /Code/Git_repos/PyBox/

Lets run the default simulation:

> python Gas_simulation.py 

Dont worry about the error message regarding the Matplotlib plots. This is a result of working in a Docker container. For those not familiar with standard Docker commands, please check the brief instructions provided in the Docker_README.txt file where I give some additional examples on how to stop, restart and delete the PyBox container. 

## Folder structure and running the model <a name="Folder-Structure"></a>

Before we run PyBox, lets make sure you have the correct link to the UManSysProp suite. Once you have cloned the repository, you will need to add the location of it in the python script 'Property_calculation.py' within the 'Aerosol' directory of PyBox. As with the Assimulo package, you can test this import by opening an interactive Python shell and typing:

> import sys

> sys.path.append('<-- add your path here -->/UManSysProp_public/')

> from umansysprop import boiling_points

> from umansysprop import vapour_pressures

> from umansysprop import critical_properties

> from umansysprop import liquid_densities

> from umansysprop import partition_models

> from umansysprop.forms import CoreAbundanceField

If you are happy all dependencies are installed and working, to run PyBox 'out of the box', type the following in the root directory:

> python Gas_simulation.py

If you are not running within a Docker container, you will see a plot displaying the concentration of two compounds over time. To understand what this simulation has actually done, let us now discuss the repository structure.

### Directory layout

    .                           # Gas phase only model [using Numba]
    ├── f2py                    # Gas phase only model [using f2py Fortran to Python Interface Generator] 
    ├── Aerosol                 # Coupled gas and gas-to-particle partitioning routines
    |------f2py                 # Coupled gas and gas-to-particle partitioning routines [using f2py]
    ├── test                    # Automated unit tests
    |------data                 # Data used in the automated unit tests
    ├── mechanism_files         # Copies of chemical mechanisms
    ├── LICENSE
    └── README.md
    
Currently there are two versions of the gas phase only model; one held within the root directory and the other in folder 'f2py':
   
#### 1) Python [using Numba]
This is the default version in the root directory. Recall the parsing of the equation file? After new python scripts are created for use within the ODE solvers, the [Numba](https://numba.pydata.org) package then compiles these before the first simulation. Numba does this as the modules are imported. You will therefore find the initial pre-simulation stages of the first simulation will take some time, but not in subsequent model simulations if you wish to study a fixed chemical mechanism. In this case Numba will not need to re-compile even when you start with new initial conditions. Once you have conducted your first simulation, you may change the following within 'Gas_simulation.py':

    files_exist = False

to

    files_exist = True

The current version of PyBox provides you with an out-of-the-box example. It is based on the MCM representation of the degredation of [Alpha-Pinene](https://en.wikipedia.org/wiki/Alpha-Pinene). The Alpha-Pinene mechanism file is stored within the 'mechanism_files' folder and referenced in the 'Gas_simulation.py' file through: 

    filename='MCM_APINENE'

As already noted, to run the model, once you are happy all dependecies are installed, type the following from the root directory:

> python Gas_simulation.py

You can modify the ambient conditions and species concentrations in 'Gas_simulation.py'. First you can define ambient conditions and simulation time, the default given as :

    temp=288.15 # Kelvin
    RH=0.5 # RH/100% [range 0-0.99]
    #Define a start time 
    hour_of_day=12.0 # 24 hr format
    simulation_time= 7200.0 # seconds
    batch_step=100.0 # seconds
    
The 'batch_step' variable allows us to define when to stop/start/record outputs from our simulation for later use. The ODE methods can provide output at every internal time-step, so it is up to the user to use this information, or not, within the output of the 'ODE_solver.py' script.  Following this, the default option for species concentrations is provided as:

    # Define initial concentrations, in pbb, of species using names from KPP file
    species_initial_conc=dict()
    species_initial_conc['O3']=18.0
    species_initial_conc['APINENE']=30.0

If you are not running within a Docker container, if you run 'Gas_simulation.py' as provided you will see a simple plot of Alpha-Pinene concentration decay over 2 hours.

#### 2) Python + Fortran [using f2py Fortran to Python Interface Generator] 
Whilst the above variant uses the Numba package, in the folder 'f2py' the same model is constructed using the [F2Py](https://docs.scipy.org/doc/numpy/f2py/)package, where functions that define the ODEs are converted into pre-compiled Fortran modules with the option to use [OpenMP](http://www.openmp.org) to exploit the number of cores available to you on any given platform. As before, please check the relevant files for defining initial conditions, species concetrations, and expect some compilation time during the first run. To run this simulation, type the following from the f2py directory:

> python Gas_simulation_f2py.py

<img src="images/Example_deafult_gas_simulation.png" width="600">
<em>Example output from the default gas phase simulation of Alpha-Pinene</em>

### Aerosol

In the Aerosol folder you can find gas-to-particle partitioning frameworks. There are two examples provided. The first, within the f2py folder, simulates the partitioning of compounds to 16 size bins again from the Alpha-Pinene chemical mechanism as this evolves over time. This uses properties calculated from the UManSysProp suite. To run the model, once you are happy all dependencies are installed, type the following from the Aerosol/f2py directory:

> python Aerosol_simulation_f2py.py

In addition to the species concentrations and ambient conditions, you can change the size distribution and number of size bins in the following:

    num_bins=16 #Number of size bins
    total_conc=100 #Total particles per cc
    std=2.2 #Standard Deviation
    lowersize=0.01 #microns
    uppersize=1.0 #microns
    meansize=0.2 #microns
    
Please note this does require some knowledge of typical aerosol size distributions and reasonable number concentrations. Within the folder 'Fixed_yield' is a partitioning only model, using fixed total concentrations of compounds in the gas phase. It is important to note that much more work is planned on the aerosol model since there are multiple properties and processes that affect gas-to-particle partitioning. The current version is beyond the most basic used in atmospheric research. Nonetheless, PyBox is designed with the community in mind and my goal is to include all relevant processes. This ethos is captured in the proceeding note on contributing to the project.


<img src="images/Example_SOA_simulation.png" width="600">
<em>Example total organic aerosol loading from the default aerosol simulation of Alpha-Pinene</em>

<img src="images/Example_SOA_contributions.png" width="600">
<em>Example normalised contributions per size bin for the default fixed yield simulations</em>


## Unit tests<a name="Automated-unit-tests"></a> 

Within the folder tests, run the following command:

> python test_modules.py -v

This will use the unittest Python module to test the output of generated functions used within the ODE solvers against pre-generated .npy files provided in the data subfolder of tests.

## Contributing<a name="Contributing"></a>

Contributions to PyBox are more than welcome. Box-models of aerosol systems can rely on many different process representations. It is thus difficult to define a 'standard' full complexity model. There are many developments planned for PyBox, which you can follow from a scientific perspective in the project wiki. I am therefore very happy to discuss ideas for improvement and how to add/remove features.  There are two key rules to follow:

 - Any addition must include appropriate unit tests
 - Any addition from a scientific process perspective must include a link to a peer-reviewed paper before it is accepted into the public branch

Please use the issue tracker at https://github.com/loftytopping/PyBox/issues if you want to notify me of an issue or need support. If you want to contribute, please either create an issue or make a pull request. Alternatively, come and see us in Manchester and/or lets meet for a coffee and a chat!

## Code of Conduct<a name="Code-of-Conduct"></a>

Please note that this project is released with a Contributor Code of Conduct. By participating in this project you agree to abide by its [terms](code-of-conduct.md). There needs to be greater support and recognition for software development and developers. PyBox can act as a vehicle for enabling better collaboration and therefore better science.

## Citation<a name="Citation"></a>

If you use PyBox in any study we ask you reference our paper in the Journal of Open Source Software [![DOI](http://joss.theoj.org/papers/10.21105/joss.00755/status.svg)](https://doi.org/10.21105/joss.00755). Citation: Topping et al., (2018). PyBox: An automated box-model generator for atmospheric chemistry and aerosol simulations. . Journal of Open Source Software, 3(28), 755, https://doi.org/10.21105/joss.00755
