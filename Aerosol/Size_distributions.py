##########################################################################################
#                                                                                        #
#    Module to create log-normal size distribution. Code copied/modified from            #
#    http://all-geo.org/volcan01010/2013/09/how-to-use-lognormal-distributions-in-python/#
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

# Last modification 28/5/17



from scipy.stats import lognorm
from scipy import stats # Import the scipy.stats module
import matplotlib.pyplot as plt
import numpy as np

def lognormal(num_bins,total_conc,meansize,std,lowersize,uppersize):
    
    # Definitions of different parameters. (float() used to avoid problems with Python integer division)
    M = float(meansize) # Geometric mean == median
    s = float(std) # Geometric standard deviation
    mu = np.log(M) # Mean of log(X)
    sigma = np.log(s) # Standard deviation of log(X)
    shape = sigma # Scipy's shape parameter
    scale = np.exp(mu) # Scipy's scale parameter
    median = np.exp(mu)
    mode = np.exp(mu - sigma**2) # Note that mode depends on both M and s
    mean = np.exp(mu + (sigma**2/2)) # Note that mean depends on both M and s
    x = np.linspace(lowersize, uppersize, num=400) # values for x-axis
    pdf = stats.lognorm.pdf(x, shape, loc=0, scale=scale) # probability distribution
    
    x_output = np.exp(np.linspace(np.log(lowersize), np.log(uppersize), num=num_bins))
    pdf_output= stats.lognorm.pdf(x_output, shape, loc=0, scale=scale)

    #plt.figure(figsize=(12,4.5))
    # Figure on linear scale
    #plt.subplot(121)
    #plt.plot(x, pdf)
    #plt.fill_between(x, pdf, where=(x < M/s), alpha=0.15)
    #plt.fill_between(x, pdf, where=(x > M*s), alpha=0.15)
    #plt.fill_between(x, pdf, where=(x < M/s**2), alpha=0.15)
    #plt.fill_between(x, pdf, where=(x > M*s**2), alpha=0.15)
    #plt.vlines(mode, 0, pdf.max(), linestyle=':', label='Mode')
    #plt.vlines(mean, 0, stats.lognorm.pdf(mean, shape, loc=0, scale=scale), linestyle='--', color='green', label='Mean')
    #plt.vlines(median, 0, stats.lognorm.pdf(median, shape, loc=0, scale=scale), color='blue', label='Median')
    #plt.ylim(ymin=0)
    #plt.xlabel('Radius (microns)')
    #plt.title('Linear scale')
    #leg=plt.legend()

    # Figure on logarithmic scale
    #plt.subplot(122)
    #plt.semilogx(x, pdf)
    #plt.fill_between(x, pdf, where=(x < M/s), alpha=0.15)
    #plt.fill_between(x, pdf, where=(x > M*s), alpha=0.15)
    #plt.fill_between(x, pdf, where=(x < M/s**2), alpha=0.15)
    ##plt.fill_between(x, pdf, where=(x > M*s**2), alpha=0.15)
    #plt.vlines(mode, 0, pdf.max(), linestyle=':', label='Mode')
    #plt.vlines(mean, 0, stats.lognorm.pdf(mean, shape, loc=0, scale=scale), linestyle='--', color='green', label='Mean')
    #plt.vlines(median, 0, stats.lognorm.pdf(median, shape, loc=0, scale=scale), color='blue', label='Median')
    #plt.ylim(ymin=0)
    #plt.xlabel('Radius (microns)')
    #plt.title('Logarithmic scale')
    #leg=plt.legend()

    #plt.show()    
    
    return (pdf_output/sum(pdf_output))*total_conc,x_output
    