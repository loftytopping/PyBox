import numpy as np
import scipy as sp
import numba as nb
#from cython_resample import cython_resample
from timeit import default_timer as timer

@nb.jit(nb.f8[:](nb.f8[:]))
def numba_cumsum(x):
    return np.cumsum(x)

@nb.jit(nb.f8(nb.f8), nopython=True)
def numba_exp(x):
    return np.exp(x)

@nb.jit(nb.f8(nb.f8), nopython=True)
def numba_log10(x):
    return np.log10(x)

@nb.autojit
def numba_resample(qs, xs, rands):
    n = qs.shape[0]
    lookup = numba_cumsum(qs)
    results = np.empty(n)

    for j in range(n):
        for i in range(n):
            if rands[j] < lookup[i]:
                results[j] = xs[i]
                break
    return results

def python_resample(qs, xs, rands):
    n = qs.shape[0]
    lookup = np.cumsum(qs)
    results = np.empty(n)

    for j in range(n):
        for i in range(n):
            if rands[j] < lookup[i]:
                results[j] = xs[i]
                break
    return results

def numpy_resample(qs, xs, rands):
    results = np.empty_like(qs)
    lookup = sp.cumsum(qs)
    for j, key in enumerate(rands):
        i = sp.argmax(lookup>key)
        results[j] = xs[i]
    return results

def evaluate_rates_normal(temp):
    # Creating reference to constant values used in rate expressions
    # constants used in calculation of reaction rates
    M  = 2.55E+19  #Check this against pressure assumed in model
    N2 = 0.79*M
    O2 = 0.2095*M
    # kro2no : ro2      + no      = ro      + no2
    #        : ro2      + no      = rono2
    # iupac 1992
    KRONO2    = 2.40E-12*np.exp(360.0/temp)
    # kro2ho2: ro2      + ho2     = rooh    + o2
    # mcm protocol [1997]
    KRO2HO2   = 2.91E-13*np.exp(1300.0/temp)
    # kapho2 : rcoo2    + ho2     = products
    # mcm protocol [1997]
    KAPHO2    = 4.30E-13*np.exp(1040.0/temp)
    # kapno  : rcoo2    + no      = products
    # mej [1998]
    KAPNO = 8.10E-12*np.exp(270.0/temp)
    # kro2no3: ro2      + no3     = products
    # mcm protocol [1997]
    KRO2NO3   = 2.50E-12
    # kno3al : no3      + rcho    = rcoo2   + hno3
    # mcm protocol [1997]
    KNO3AL    = 1.44E-12*np.exp(-1862.0/temp)
    # kdec   : ro                 = products
    # mcm protocol [1997]
    KDEC      = 1.00E+06
    KROSEC = 1.50E-14*np.exp(-200.0/temp)
    KALKOXY=3.70E-14*np.exp(-460.0/temp)*O2
    KALKPXY=1.80E-14*np.exp(-260.0/temp)*O2
    # -------------------------------------------------------------------
    # complex reactions
    # -------------------------------------------------------------------
    # kfpan kbpan
    # formation and decomposition of pan
    # iupac 2001 (mcmv3.2)
    kc0     = 2.70E-28*M*(temp/300.0)**(-7.1)
    kci     = 1.21E-11*(temp/300.0)**(-0.9)
    krc     = kc0/kci
    fcc     = 0.30
    nc      = 0.75-(1.27*np.log10(fcc))
    fc      = 10**(np.log10(fcc)/(1.0+((np.log10(krc))/nc)**2.0))
    KFPAN   = (kc0*kci)*fc/(kc0+kci)
    
    return KFPAN

@nb.jit('float64(float64)', nopython=True)
def evaluate_rates_numba(temp):
    # Creating reference to constant values used in rate expressions
    # constants used in calculation of reaction rates
    M  = 2.55E+19  #Check this against pressure assumed in model
    N2 = 0.79*M
    O2 = 0.2095*M
    # kro2no : ro2      + no      = ro      + no2
    #        : ro2      + no      = rono2
    # iupac 1992
    KRONO2    = 2.40E-12*numba_exp(360.0/temp)
    # kro2ho2: ro2      + ho2     = rooh    + o2
    # mcm protocol [1997]
    KRO2HO2   = 2.91E-13*numba_exp(1300.0/temp)
    # kapho2 : rcoo2    + ho2     = products
    # mcm protocol [1997]
    KAPHO2    = 4.30E-13*numba_exp(1040.0/temp)
    # kapno  : rcoo2    + no      = products
    # mej [1998]
    KAPNO = 8.10E-12*numba_exp(270.0/temp)
    # kro2no3: ro2      + no3     = products
    # mcm protocol [1997]
    KRO2NO3   = 2.50E-12
    # kno3al : no3      + rcho    = rcoo2   + hno3
    # mcm protocol [1997]
    KNO3AL    = 1.44E-12*numba_exp(-1862.0/temp)
    # kdec   : ro                 = products
    # mcm protocol [1997]
    KDEC      = 1.00E+06
    KROSEC = 1.50E-14*numba_exp(-200.0/temp)
    KALKOXY=3.70E-14*numba_exp(-460.0/temp)*O2
    KALKPXY=1.80E-14*numba_exp(-260.0/temp)*O2
    # -------------------------------------------------------------------
    # complex reactions
    # -------------------------------------------------------------------
    # kfpan kbpan
    # formation and decomposition of pan
    # iupac 2001 (mcmv3.2)
    kc0     = 2.70E-28*M*(temp/300.0)**(-7.1)
    kci     = 1.21E-11*(temp/300.0)**(-0.9)
    krc     = kc0/kci
    fcc     = 0.30
    nc      = 0.75-(1.27*numba_log10(fcc))
    fc      = 10**(numba_log10(fcc)/(1.0+((numba_log10(krc))/nc)**2.0))
    KFPAN   = (kc0*kci)*fc/(kc0+kci)

    return KFPAN


if __name__ == '__main__':
    n = 100
    xs = np.arange(n, dtype=np.float64)
    qs = np.array([1.0/n,]*n)
    rands = np.random.rand(n)

    start = timer()
    for i in range(1000):
        numba_resample(qs, xs, rands)
    end = timer()
    print("Timing Numba Function:",end-start)
    start = timer()
    for i in range(1000):
        python_resample(qs, xs, rands)
    end = timer()    
    print("Timing Python Function:",end-start)
    #%timeit python_resample(qs, xs, rands)
    start = timer()
    for i in range(1000):
        numpy_resample(qs, xs, rands)
    end = timer()    
    print("Timing Numpy Function:",end-start)
    #%timeit numpy_resample(qs, xs, rands)
    #print("Timing Cython Function:")
    #%timeit cython_resample(qs, xs, rands)

    temp=298.15
    start = timer()
    for i in range(1000):
        evaluate_rates_numba(temp)
    end = timer()
    print("Timing Numba Function:",end-start)
    start = timer()
    for i in range(1000):
        evaluate_rates_normal(temp)
    end = timer()    
    print("Timing Python Function:",end-start)
