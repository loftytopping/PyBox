# PyBox

PyBox is a 0-D box model, where all species are homogeneously distributed, built around a chemical mechanism file; a file that represents all the individual chemical reactions taking place
starting from a wide range of precursors. By parsing a chemical equation file, obtained from the MCM project, PyBox creates files that account for the gas phase chemistry as well as automatically
calculating properties that dictate gas-to-particle partitioning through connection with the UManSysProp informatics suite. Written in Python, PyBox uses Numba or the
Fortran-to-Python-Generator f2py to perform calculations currently within a library of ODE solvers provided by the Assimulo package.

PyBox is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.                                                                            

PyBox is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details. You should have received a copy of the GNU General Public License along with PyBox.  If not, see <http://www.gnu.org/licenses/>