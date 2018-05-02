##########################################################################################
#                                                                                        #
#    Script to compile Fortran module from Python using f2py                             #
#                                                                                        #
#                                                                                        #
#    Copyright (C) 2018  David Topping : david.topping@manchester.ac.uk                  #
#                                      : davetopp80@gmail.com                            #
#    Personal website: davetoppingsci.com                                                #
#                                                                                        #
#    This program does not yet have a license, meaning the deault copyright law applies. #
#    I will add an appropriate open-source icense once made public with paper            #
#    Only users who have access to the private repository that holds this file may       #
#    use it, but may not distribute it without explicit permission.                      #
#                                                                                        #
#                                                                                        #
##########################################################################################

from numpy.distutils.core import Extension
ext = Extension (name = "partition_f2py",
                 sources = ["Partitioning.f90"], 
                 extra_compile_args=["-O3",
                                     "-ffast-math",
                                     "-fopenmp"],
                 extra_link_args=["-lgomp",
                                  "-static",
                                  "-llapack",
                                  "-lblas"])
if __name__ == '__main__':
    from numpy.distutils.core import setup 
    setup (name = "partition_f2py" ,
           ext_modules = [ext])

    