# PyBox

This repository will hold model variants that ultimately couple gas/aerosol chemistry, gas-to-particle partitioning and varying phase state through the full humidity cycle in the atmosphere. The first phase of the project is to develop and profile a gas phase only model, using the [Master Chemical Mechanism (MCM)](http://mcm.leeds.ac.uk/MCM/) as the basis. The model will also relate component properties, using molecular structural information, through the [UManSysProp](http://umansysprop.seaes.manchester.ac.uk) informatics suite.  Any public release will occur according to agreement from any partner contributions and associated papers.

Please check the project wiki page for more inforamtion on updates and planned releases

The model works on the basis of reading a file that defines reactions in the gas phase. Within this repo, some examples are given as taken from the MCM. For example, take the mixed VOC file given by 'MCM_mixed_test.eqn.txt'. This contains the following:

##### {46.} 	 CH3OH + OH = HO2 + HCHO : 	2.85D-12*EXP(-345/TEMP) 	;
##### {47.} 	 C2H5OH + OH = C2H5O : 	3.0D-12*EXP(20/TEMP)*0.05 	;
{48.} 	 C2H5OH + OH = CH3CHO + HO2 : 	3.0D-12*EXP(20/TEMP)*0.9 	;
{49.} 	 C2H5OH + OH = HOCH2CH2O2 : 	3.0D-12*EXP(20/TEMP)*0.05 	;

Where the equation number os first defined, then the reactants/products along with a defined rate coefficient. Some reactions rely on coefficients defined elsewhere, according to the MCM version number. This file is first parsed in the file 'Parse_eqn_file.py', providing information that can be used to set up and solve the relevant ODEs. 


