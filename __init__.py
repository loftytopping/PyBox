##########################################################################################
#                                                                                        #
#    Example gas phase model. This takes an equation file, given in the KPP format,      #
#    and then sets up an ODE solver where initial concentrations of any specie can       #
#    be set. It also relies on some pre-defined rate coefficients and photolysis rates   #
#    taken from the MCM. These are explicitly written into the relevant Fortran modules  #
#    in the file Parse_eqn_file.write_rate_file_fortran() which is provided with a 
#    Fortran syntax version of pre-defined rates in  
#                                                                                        #
#                                                                                        #
#    Copyright (C) 2017  David Topping : david.topping@manchester.ac.uk                  #
#                                      : davetopp80@gmail.com                            #
#    Personal website: davetoppingsci.com                                                #
#                                                                                        #
#    This program does not yet have a license, meaning the deault copyright law applies. #
#    Only users who have access to the private repository that holds this file may       #
#    use it or develop it, but may not distribute it without explicit permission.        #
#                                                                                        #
#                                                                                        #
##########################################################################################