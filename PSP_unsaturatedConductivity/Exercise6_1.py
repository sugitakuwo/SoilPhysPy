#PSP_unsaturatedConductivity.py
from __future__ import print_function
from PSP_readDataFile import readDataFile
import matplotlib.pyplot as plt
import numpy as np
import math

NODATA = -9999

def ComputeCampbell(b,psi,psie,Ks):
    if psi < psie :
        Kvalue=Ks*(psie/psi)**(2.0e0+3.0e0/b)
    elif psi >= psie :
        Kvalue=Ks

    return Kvalue

def main():

#-- Define the parameters for Eq. (6.35) ----------------#    
#-- For sand ------------
    bd=1.7e0
    psied=-0.7e0#[J/kg]
    Ksd=5.8e-3#[kg s m^-3]

#-- For silt loam -----------
    btl=4.7e0
    psietl=-2.1e0
    Kstl=0.19e-3#[kg s m^-3]

#-- For clay ------------
    by=7.6e0
    psiey=-3.7e0
    Ksy=0.017e-3#[kg s m^-3]

    pmin=0.01
    pmax=3.0e5
    waterPot=np.linspace(np.log10(pmax),np.log10(pmin),30)
    waterPot=-10.0**waterPot
    np1=len(waterPot)#Number of conditons on water potential.

    conductivity = np.zeros([np1,3],float)

    for ii in range(3):
        if ii==0 :
            b1=bd
            psie1=psied
            Ks1=Ksd
        elif ii==1 :
            b1=btl
            psie1=psietl
            Ks1=Kstl
        elif ii==2 :
            b1=by
            psie1=psiey
            Ks1=Ksy
                
        for i in range(np1):
            conductivity[i,ii] = ComputeCampbell(b1, waterPot[i], psie1, Ks1)
    
    plt.figure(figsize=(10,8))        
    plt.loglog (-waterPot, conductivity[:,0], 'ko-', ms=8, label="Sand")
    plt.loglog (-waterPot, conductivity[:,1], 'rD-', ms=8, label="Silt loam")
    plt.loglog (-waterPot, conductivity[:,2], 'b^-', ms=8, label="Clay")
    plt.xlabel('Water Potential [J kg$^{-1}$]',fontsize=20,labelpad=2)
    plt.ylabel('Hydraulic Conductivity [kg s m$^{-3}$]',fontsize=20,labelpad=2)
    plt.tick_params(axis='both', which='major', labelsize=20,pad=6)
    plt.tick_params(axis='both', which='minor', labelsize=20,pad=6)
    plt.legend(loc='best',fontsize=14)
    plt.show()   
        
main()
