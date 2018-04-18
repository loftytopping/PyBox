#Script containing all relevant MCM modules. This file was created from a manual copy of the
#kpp_constants.f90 file found on mcm.leeds.ac.uk/MCM folder 23/3/17

import numpy

def mcm_constants(time, temp, H2O):
    # calculates rate constants from arrhenius informtion
    # Documentation with MCM does not indicate unites of time or temp!!

    # F90 archive - REAL(dp) time, temp, M, N2, O2, RO2, H2O

    # F90 archive - OPEN(13,file='photolysis.txt', status='old')

    # calculate zenith angle for time of day
    zenith_dict=zenith(time)
    theta=zenith_dict['theta']
    secx=zenith_dict['secx']
    cosx=zenith_dict['cosx']

    #return a dictionary containing all relevant variables
    mcm_constants_dict=dict()

    # ************************************************************************
    # define generic reaction rates.
    # ************************************************************************

    # constants used in calculation of reaction rates
    M  = 2.55E+19  #Check this against pressure assumed in model
    N2 = 0.79*M
    O2 = 0.2095*M

    # kro2no : ro2      + no      = ro      + no2
    #        : ro2      + no      = rono2
    # iupac 1992
    kro2no    = 2.40E-12*numpy.exp(360.0/temp)
    mcm_constants_dict['KRONO2']=kro2no

    # kro2ho2: ro2      + ho2     = rooh    + o2
    # mcm protocol [1997]
    kro2ho2   = 2.91E-13*numpy.exp(1300.0/temp)
    mcm_constants_dict['KRO2HO2']=kro2ho2

    # kapho2 : rcoo2    + ho2     = products
    # mcm protocol [1997]
    kapho2    = 4.30E-13*numpy.exp(1040.0/temp)
    mcm_constants_dict['KAPHO2']=kapho2 

    # kapno  : rcoo2    + no      = products
    # mej [1998]
    kapno = 8.10E-12*numpy.exp(270.0/temp)
    mcm_constants_dict['KAPNO']=kapno 

    # kro2no3: ro2      + no3     = products
    # mcm protocol [1997]
    kro2no3   = 2.50E-12
    mcm_constants_dict['KRO2NO3']=kro2no3

    # kno3al : no3      + rcho    = rcoo2   + hno3
    # mcm protocol [1997]
    kno3al    = 1.44E-12*numpy.exp(-1862.0/temp)
    mcm_constants_dict['KNO3AL']=kno3al 

    # kdec   : ro                 = products
    # mcm protocol [1997]
    kdec      = 1.00E+06
    mcm_constants_dict['KDEC']=kdec

    krosec = 1.50E-14*numpy.exp(-200.0/temp)
    mcm_constants_dict['KROSEC']=krosec

    kalkoxy=3.70E-14*numpy.exp(-460.0/temp)*O2
    mcm_constants_dict['KALKOXY']=kalkoxy

    kalkpxy=1.80E-14*numpy.exp(-260.0/temp)*O2
    mcm_constants_dict['KALKPXY']=kalkpxy

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
    nc      = 0.75-(1.27*numpy.log10(fcc))
    fc      = 10**(numpy.log10(fcc)/(1.0+((numpy.log10(krc))/nc)**2.0))
    kfpan   = (kc0*kci)*fc/(kc0+kci)
    mcm_constants_dict['KFPAN']=kfpan

    kd0     = 4.90E-03*M*numpy.exp(-12100.0/temp)
    kdi     = 5.40E+16*numpy.exp(-13830.0/temp)
    krd     = kd0/kdi
    fcd     = 0.30
    ncd     = 0.75-(1.27*numpy.log10(fcd))
    fd      = 10.0**(numpy.log10(fcd)/(1.0+((numpy.log10(krd))/ncd)**2.0))
    kbpan   = (kd0*kdi)*fd/(kd0+kdi)
    mcm_constants_dict['KBPAN']=kbpan

    # kmt01  : o        + no      = no2
    # iupac 2001 (mcmv3.2)
    k10     = 1.00E-31*M*(temp/300.0)**(-1.6)

    k1i     = 3.00E-11*(temp/300.0)**(0.3)
    kr1     = k10/k1i
    fc1     = 0.85
    nc1     = 0.75-(1.27*numpy.log10(fc1))
    f1      = 10.0**(numpy.log10(fc1)/(1.0+((numpy.log10(kr1)/nc1))**2.0))
    kmt01   = (k10*k1i)*f1/(k10+k1i)
    mcm_constants_dict['KMT01']=kmt01

    # kmt02  : o        + no2     = no3
    # iupac 2001 (mcmv3.2)
    k20     = 1.30E-31*M*(temp/300.0)**(-1.5)
    k2i     = 2.30E-11*(temp/300.0)**(0.24)
    kr2     = k20/k2i
    fc2     = 0.6
    nc2     = 0.75-(1.27*numpy.log10(fc2))
    f2      = 10.0**(numpy.log10(fc2)/(1.0+((numpy.log10(kr2)/nc2))**2.0))
    kmt02   = (k20*k2i)*f2/(k20+k2i)
    mcm_constants_dict['KMT02']=kmt02

    # kmt03  : no2      + no3     = n2o5
    # iupac 2006, mcmv3.2
    k30     = 3.60E-30*M*(temp/300.0)**(-4.1)
    k3i     = 1.90E-12*(temp/300.0)**(0.2)
    kr3     = k30/k3i
    fc3     = 0.35
    nc3     = 0.75-(1.27*numpy.log10(fc3))
    f3      = 10.0**(numpy.log10(fc3)/(1.0+((numpy.log10(kr3)/nc3))**2.0))
    kmt03   = (k30*k3i)*f3/(k30+k3i)
    mcm_constants_dict['KMT03']=kmt03

    # kmt04  : n2o5               = no2     + no3
    # iupac 2006, mcmv3.2
    k40     = 1.30E-03*M*(temp/300.0)**(-3.5)*numpy.exp(-11000.0/temp)
    k4i     = 9.70E+14*(temp/300.0)**(0.1)*numpy.exp(-11080.0/temp)
    kr4     = k40/k4i
    fc4     = 0.35
    nc4     = 0.75-(1.27*numpy.log10(fc4))
    f4      = 10.0**(numpy.log10(fc4)/(1+((numpy.log10(kr4)/nc4))**2.0))
    kmt04   = (k40*k4i)*f4/(k40+k4i)
    mcm_constants_dict['KMT04']=kmt04

    # kmt05  : oh       + co(+o2) = ho2     + co2
    # iupac 2006
    kmt05  = (1.0 + (M/4.2E19))
    mcm_constants_dict['KMT05']=kmt05

    # kmt06  : ho2      + ho2     = h2o2    + o2
    # water enhancement factor
    # iupac 1992
    kmt06  = 1.0 + (1.40E-21*numpy.exp(2200.0/temp)*H2O)
    mcm_constants_dict['KMT06']=kmt06

    # kmt06  = 1.0 + (2.00E-25*numpy.exp(4670.0/temp)*h2o)
    # S+R 2005 values

    # kmt07  : oh       + no      = hono

    # iupac 2006, mcmv3.2
    k70     = 7.40E-31*M*(temp/300.0)**(-2.4)
    k7i     = 3.30E-11*(temp/300.0)**(-0.3)
    kr7     = k70/k7i
    fc7     = numpy.exp(-temp/1420.0)
    nc7     = 0.75-(1.27*numpy.log10(fc7))
    f7      = 10.0**(numpy.log10(fc7)/(1+((numpy.log10(kr7)/nc7))**2.0))
    kmt07   = (k70*k7i)*f7/(k70+k7i)
    mcm_constants_dict['KMT07']=kmt07

    # kmt08  : oh       + no2     = hno3

    # iupac 2006, mcmv3.2
    k80     = 3.30E-30*M*(temp/300.0)**(-3.0)
    k8i     = 4.10E-11
    kr8     = k80/k8i
    fc8     = 0.4
    nc8     = 0.75-(1.27*numpy.log10(fc8))
    f8      = 10.0**(numpy.log10(fc8)/(1.0+((numpy.log10(kr8)/nc8))**2.0))
    kmt08   = (k80*k8i)*f8/(k80+k8i)
    mcm_constants_dict['KMT08']=kmt08

    # kmt09  : ho2      + no2     = ho2no2
    # iupac 1997, mcmv3.2

    k90     = 1.80E-31*M*(temp/300.0)**(-3.2)
    k9i     = 4.70E-12
    kr9     = k90/k9i
    fc9     = 0.6
    nc9     = 0.75-(1.27*numpy.log10(fc9))
    f9      = 10.0**(numpy.log10(fc9)/(1.0+((numpy.log10(kr9)/nc9))**2.0))
    kmt09   = (k90*k9i)*f9/(k90+k9i)
    mcm_constants_dict['KMT09']=kmt09

    # kmt10  : ho2no2             = ho2     + no2
    # iupac 1997, mcmv3.2

    k100     = 4.10E-05*M*numpy.exp(-10650.0/temp)
    k10i     = 4.80E+15*numpy.exp(-11170.0/temp)
    kr10     = k100/k10i
    fc10     = 0.6
    nc10     = 0.75-(1.27*numpy.log10(fc10))
    f10      = 10.0**(numpy.log10(fc10)/(1.0+((numpy.log10(kr10)/nc10))**2.0))
    kmt10    = (k100*k10i)*f10/(k100+k10i)
    mcm_constants_dict['KMT10']=kmt10

    # kmt11  : oh       + hno3    = h2o     + no3
    # iupac 2006, mcmv3.2

    k1     = 2.40E-14*numpy.exp(460.0/temp)
    k3     = 6.50E-34*numpy.exp(1335.0/temp)
    k4     = 2.70E-17*numpy.exp(2199.0/temp)
    k2     = (k3*M)/(1.0+(k3*M/k4))
    kmt11  = k1 + k2
    mcm_constants_dict['KMT11']=kmt11

    # kmt12 iupac 2006, mcmv3.2

    k120 = 4.50E-31*((temp/300.0)**(-3.9))*M
    k12i = 1.30E-12*((temp/300.0)**(-0.7))
    kr12 = k120/k12i
    fc12 = 0.525
    nc12 = 0.75-(1.27*numpy.log10(fc12))
    f12  = 10.0**(numpy.log10(fc12)/(1.0+((numpy.log10(kr12)/nc12))**2.0))
    kmt12    = (k120*k12i)*f12/(k120+k12i)
    mcm_constants_dict['KMT12']=kmt12

    # kmt13  : ch3o2    + no2     = ch3o2no2
    # iupac 2006

    k130     = 2.50E-30*((temp/300.0)**(-5.5))*M
    k13i     = 1.80E-11
    kr13     = k130/k13i
    fc13     = 0.36
    nc13     = 0.75-(1.27*numpy.log10(fc13))
    f13      = 10.0**(numpy.log10(fc13)/(1.0+((numpy.log10(kr13)/nc13))**2.0))
    kmt13    = (k130*k13i)*f13/(k130+k13i)
    mcm_constants_dict['KMT13']=kmt13

    # kmt14  : ch3o2no2           = ch3o2   + no2
    # iupac 2006, mcmv3.2

    k140     = 9.00E-05*numpy.exp(-9690.0/temp)*M
    k14i     = 1.10E+16*numpy.exp(-10560.0/temp)
    kr14     = k140/k14i
    fc14     = 0.4
    nc14     = 0.75-(1.27*numpy.log10(fc14))
    f14      = 10.0**(numpy.log10(fc14)/(1.0+((numpy.log10(kr14)/nc14))**2.0))
    kmt14    = (k140*k14i)*f14/(k140+k14i)
    mcm_constants_dict['KMT14']=kmt14

    # kmt15 iupac 2006, mcmv3.2

    k150 = 8.60E-29*((temp/300.0)**(-3.1))*M
    k15i = 9.00E-12*((temp/300.0)**(-0.85))
    kr15 = k150/k15i
    fc15 = 0.48
    nc15 = 0.75-(1.27*numpy.log10(fc15))
    f15  = 10.0**(numpy.log10(fc15)/(1.0+((numpy.log10(kr15)/nc15))**2.0))
    kmt15 = (k150*k15i)*f15/(k150+k15i)
    mcm_constants_dict['KMT15']=kmt15

    # kmt16  :  oh  +  c3h6
    # iupac 2006

    k160     = 8.00E-27*((temp/300.0)**(-3.5))*M
    k16i     = 3.00E-11*((temp/300.0)**(-1.0))
    kr16     = k160/k16i
    fc16     = 0.5
    nc16     = 0.75-(1.27*numpy.log10(fc16))
    f16      = 10.0**(numpy.log10(fc16)/(1.0+((numpy.log10(kr16)/nc16))**2.0))
    kmt16    = (k160*k16i)*f16/(k160+k16i)
    mcm_constants_dict['KMT16']=kmt16

    # kmt17 iupac 2006

    k170 = 5.00E-30*((temp/300.0)**(-1.5))*M
    k17i = 1.00E-12
    kr17 = k170/k17i
    fc17 = (0.17*numpy.exp(-51./temp))+numpy.exp(-1.0*temp/204.)
    nc17 = 0.75-(1.27*numpy.log10(fc17))
    f17  = 10.0**(numpy.log10(fc17)/(1.0+((numpy.log10(kr17)/nc17))**2.0))
    kmt17 = (k170*k17i)*f17/(k170+k17i)
    mcm_constants_dict['KMT17']=kmt17

    #Taken from yet another MCM version, will need converting to lower case
    KRO2NO = 2.7E-12*numpy.exp(360/temp)
    mcm_constants_dict['KRO2NO']=KRO2NO

    KRO2HO2 = 2.91E-13*numpy.exp(1300/temp)
    mcm_constants_dict['KRO2HO2']=KRO2HO2

    KAPHO2 = 5.2E-13*numpy.exp(980/temp)
    mcm_constants_dict['KAPHO2']=KAPHO2

    KAPNO = 7.5E-12*numpy.exp(290/temp)
    mcm_constants_dict['KAPNO']=KAPNO

    KRO2NO3 = 2.3E-12
    mcm_constants_dict['KRO2NO3']=KRO2NO3

    KNO3AL = 1.4E-12*numpy.exp(-1860/temp)
    mcm_constants_dict['KNO3AL']=KNO3AL

    KDEC = 1.00E+06
    mcm_constants_dict['KDEC']=KDEC

    KROPRIM = 2.50E-14*numpy.exp(-300/temp)
    mcm_constants_dict['KROPRIM']=KROPRIM

    KROSEC = 2.50E-14*numpy.exp(-300/temp)
    mcm_constants_dict['KROSEC']=KROSEC

    KCH3O2 = 1.03E-13*numpy.exp(365/temp)
    mcm_constants_dict['KCH302']=KCH3O2

    K298CH3O2 = 3.5E-13
    mcm_constants_dict['K298CH3O2']=K298CH3O2

    kmt18 = 9.5E-39*O2*numpy.exp(5270/temp)/(1+7.5E-29*O2*numpy.exp(5610/temp))
    mcm_constants_dict['KMT18']=kmt18

    kroprim  = 3.70E-14*numpy.exp(-460.0/temp)
    mcm_constants_dict['KROPRIM']=kroprim

    krosec   = 1.80E-14*numpy.exp(-260.0/temp)
    mcm_constants_dict['KROSEC']=krosec

    # ************************************************************************
    # define photolysis reaction rates using derwent method from mcm2box.fac
    # ************************************************************************

    # Data taken from photolysis.txt. Calculations done in the form of:
    # j(k) = l(k)*cosx**( mm(k))*numpy.exp(-nn(k)*secx)
    J=[None]*62

    #J          L           M          N
    J[1]=6.073E-05*cosx**(1.743)*numpy.exp(-1.0*0.474*secx)
    J[2]=4.775E-04*cosx**(0.298)*numpy.exp(-1.0*0.080*secx)
    J[3]=1.041E-05*cosx**(0.723)*numpy.exp(-1.0*0.279*secx)
    J[4]=1.165E-02*cosx**(0.244)*numpy.exp(-1.0*0.267*secx)
    J[5]=2.485E-02*cosx**(0.168)*numpy.exp(-1.0*0.108*secx)
    J[6]=1.747E-01*cosx**(0.155)*numpy.exp(-1.0*0.125*secx)
    J[7]=2.644E-03*cosx**(0.261)*numpy.exp(-1.0*0.288*secx)
    J[8]=9.312E-07*cosx**(1.230)*numpy.exp(-1.0*0.307*secx)
    J[11]=4.642E-05*cosx**(0.762)*numpy.exp(-1.0*0.353*secx)
    J[12]=6.853E-05*cosx**(0.477)*numpy.exp(-1.0*0.323*secx)
    J[13]=7.344E-06*cosx**(1.202)*numpy.exp(-1.0*0.417*secx)
    J[14]=2.879E-05*cosx**(1.067)*numpy.exp(-1.0*0.358*secx)
    J[15]=2.792E-05*cosx**(0.805)*numpy.exp(-1.0*0.338*secx)
    J[16]=1.675E-05*cosx**(0.805)*numpy.exp(-1.0*0.338*secx)
    J[17]=7.914E-05*cosx**(0.764)*numpy.exp(-1.0*0.364*secx)
    J[18]=1.140E-05*cosx**(0.396)*numpy.exp(-1.0*0.298*secx)
    J[19]=1.140E-05*cosx**(0.396)*numpy.exp(-1.0*0.298*secx)
    J[21]=7.992E-07*cosx**(1.578)*numpy.exp(-1.0*0.271*secx)
    J[22]=5.804E-06*cosx**(1.092)*numpy.exp(-1.0*0.377*secx)
    J[23]=1.836E-05*cosx**(0.395)*numpy.exp(-1.0*0.296*secx)
    J[24]=1.836E-05*cosx**(0.395)*numpy.exp(-1.0*0.296*secx)
    J[31]=6.845E-05*cosx**(0.130)*numpy.exp(-1.0*0.201*secx)
    J[32]=1.032E-05*cosx**(0.130)*numpy.exp(-1.0*0.201*secx)
    J[33]=3.802E-05*cosx**(0.644)*numpy.exp(-1.0*0.312*secx)
    J[34]=1.537E-04*cosx**(0.170)*numpy.exp(-1.0*0.208*secx)
    J[35]=3.326E-04*cosx**(0.148)*numpy.exp(-1.0*0.215*secx)
    J[41]=7.649E-06*cosx**(0.682)*numpy.exp(-1.0*0.279*secx)
    J[51]=1.588E-06*cosx**(1.154)*numpy.exp(-1.0*0.318*secx)
    J[52]=1.907E-06*cosx**(1.244)*numpy.exp(-1.0*0.335*secx)
    J[53]=2.485E-06*cosx**(1.196)*numpy.exp(-1.0*0.328*secx)
    J[54]=4.095E-06*cosx**(1.111)*numpy.exp(-1.0*0.316*secx)
    J[55]=1.135E-05*cosx**(0.974)*numpy.exp(-1.0*0.309*secx)
    J[56]=7.549E-06*cosx**(1.015)*numpy.exp(-1.0*0.324*secx)
    J[57]=3.363E-06*cosx**(1.296)*numpy.exp(-1.0*0.322*secx)
    J[61]=7.537E-04*cosx**(0.499)*numpy.exp(-1.0*0.266*secx)
    mcm_constants_dict['J']=J

    return mcm_constants_dict

  #***************************************************************************

def zenith(ttime):

    #F90 - archive REAL(dp), INTENT(IN)  :: ttime
    #F90 - archive REAL(dp), INTENT(OUT) :: theta, secx, cosx
    #F90 - archive REAL(dp) lat, pi, radian, dec, lha, sinld, cosld

    # solar declination angle 
    dec = 23.79
    # latitude
    lat = 50.0
    pi = 4.0*numpy.arctan(1.0)
    # local hour angle - representing time of day
    lha = (1.0+ttime/4.32E+4)*pi
    radian = 180.0/pi
    lat = lat/radian
    dec = dec/radian
    theta = numpy.arccos(numpy.cos(lha)*numpy.cos(dec)*numpy.cos(lat)+numpy.sin(dec)*numpy.sin(lat))
    sinld = numpy.sin(lat)*numpy.sin(dec)
    cosld = numpy.cos(lat)*numpy.cos(dec)
    cosx = (numpy.cos(lha)*cosld)+sinld
    cosx = numpy.cos(theta)
    secx = 1.0E+0/(cosx+1.0E-30)

    zenith_dict=dict()
    zenith_dict['theta']=theta
    zenith_dict['cosx']=cosx
    zenith_dict['secx']=secx

    return zenith_dict

  #***************************************************************************
  