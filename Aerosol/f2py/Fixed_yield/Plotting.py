##########################################################################################
#                                                                                        #
#    Plotting scripts for parcel model, designed to look at the evolution of aerosol     #
#    particles via a  box model.                                                         #
#    This will act as a base for adding additional component                             #
#    properties as time progresses.                                                      #
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

from scipy.stats import lognorm
from scipy import stats # Import the scipy.stats module
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm 
from matplotlib.animation import FuncAnimation
import pdb
import time as time_func

def stacked_bar(time,y_matrix,num_species, num_bins,molw_asnumpy,NA):
    #This function plots the stacked contributions from each component [or top x%]
    #across each size bin at a given point in the simulation. We plot boht the normalised and absolute
    #contributions to mass
    
    #y_matrix is a matrix of outputs from the simulation. The specific time stamps have 
    #been chosen prior to this function call. Thus, here, we cycle through all entries.

	
    for time_step in range(len(time)):
        
        #Define a range of colours to represent every compound
        color_array=iter(cm.rainbow(np.linspace(0,1,num_species-1)))
        
        time_stamp=time[time_step]
        y_asnumpy=np.array(y_matrix[time_step,:])
        y_abs_matrix=np.zeros((num_species-1,num_bins),)
        y_norm_matrix=np.zeros((num_species-1,num_bins),)
        y_cummulative=np.zeros((1,num_bins),)
        
        for size in range(num_bins):
        
            size+=1
            #The following produces g/cc - thus *1.0e12 to get to micrograms/m3
            y_mass_array=np.multiply(y_asnumpy[num_species+((size-1)*num_species):num_species+((size)*num_species),]/NA,molw_asnumpy)*1.0e12
            y_mass_array[np.where(y_mass_array<0.0)]=0.0
            #We need to now take this information as a transpose and populate a new matrix of concentrations
            #Recall the last entry is water so we can ignore.
            #Within error, you might have 'some' negative values here. For expediency, justse setting these to postive
            #but need to create a more generic response here.
            y_abs_matrix[:,size-1]=np.transpose(y_mass_array[0:num_species-1])
            y_norm_matrix[:,size-1]=np.transpose(y_mass_array[0:num_species-1]/(np.sum(y_mass_array[0:num_species-1])))
            
        #Now plot the results
        fig, ax = plt.subplots()
        bar_locations = np.arange(num_bins)
        for org in range(num_species-1):
            c=next(color_array)
            if org == 0:
                ax.bar(bar_locations, y_norm_matrix[org,:],color=c)
            else:
                y_cummulative[0,:]=y_cummulative[0,:]+y_norm_matrix[org-1,:]
                ax.bar(bar_locations, y_norm_matrix[org,:],bottom=np.array(y_cummulative[0,:]),color=c)
    
        plt.title('SOA mass = {:.2f}'.format(np.sum(np.sum(y_abs_matrix))))
        plt.ylabel('Normalised SOA contribution')
        plt.xlabel('Size bin')
        plt.show()
        #time_func.sleep(3) 
        pdb.set_trace()
        plt.close(fig)
    
    