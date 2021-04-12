#PSP_unsaturatedConductivity.py
from __future__ import print_function
from PSP_readDataFile import readDataFile
import matplotlib.pyplot as plt
import numpy as np
import math

NODATA = -9999

def ComputeKs(rhob,b,my,mt):
    Kvalue=4.0e-3*(1.3e3/rhob)**(1.3e0*b)*np.exp(-6.9e0*my-3.7e0*mt)
    return Kvalue

def main():

#-- Define the parameters for Eq. (6.35) ----------------#    
#-- For silt loam -----------
    btl=4.7e0
    mytl=0.15e0
    mttl=0.65e0

#-- For clay ------------
    by=7.6e0
    myy=0.20e0
    mty=0.60e0

    rhomin=900.0
    rhomax=1500.0
    rhob1=np.linspace(rhomin,rhomax,30)
    np1=len(rhob1)#Number of conditons on bulk density

    conductivity = np.zeros([np1,2],float)

    for ii in range(2):
        if ii==0 :
            b1=btl
            my1=mytl
            mt1=mttl
        elif ii==1 :
            b1=by
            my1=myy
            mt1=mty
                
        for i in range(np1):
            conductivity[i,ii] = ComputeKs(rhob1[i], b1, my1, mt1)
    
    plt.figure(figsize=(10,8))        
    plt.plot (rhob1, conductivity[:,0], 'rD-', ms=8, label="Silt loam")
    plt.plot (rhob1, conductivity[:,1], 'b^-', ms=8, label="Clay")
    plt.xlabel('Bulk density[kg m$^{-3}$]',fontsize=20,labelpad=2)
    plt.ylabel('Hydraulic Conductivity [kg s m$^{-3}$]',fontsize=20,labelpad=2)
    plt.tick_params(axis='both', which='major', labelsize=20,pad=6)
    plt.tick_params(axis='both', which='minor', labelsize=20,pad=6)
    plt.legend(loc='best',fontsize=14)
    plt.show()   
        
main()
