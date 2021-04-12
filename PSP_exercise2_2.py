# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 21:35:31 2018

@author: Takuya
"""

#PSP_basicProperties
from __future__ import print_function, division


def computePorosity(bulkDensity, gravWaterContent):
    waterDensity=1000
    particleDensity=2650
    porosity=1.0-(bulkDensity/particleDensity)#Fluid filled ratio
    voidRatio=porosity/(1.0-porosity)
    waterContent=gravWaterContent*(bulkDensity/waterDensity)#mass wetness
    gasPorosity=porosity-waterContent
    degreeSaturation=waterContent/porosity
    print("\nTotal porosity [m^3/m^3] =",format(porosity, ".3f"))
    print("Void ratio [m^3/m^3] =",format(voidRatio, ".3f"))
    print("Volumetric water content [m^3/m^3] =",format(waterContent, ".3f"))
    print("Gas-filled porosity [m^3/m^3] =",format(gasPorosity, ".3f"))
    print("Degree of saturation [-] =",format(degreeSaturation, ".3f"))
    return

def computeSaturationWetness(bulkDensity):
    waterDensity=1000
    particleDensity=2650
    porosity=1.0-(bulkDensity/particleDensity)#Fluid filled ratio
    return (porosity/(bulkDensity/waterDensity))

def main():
    SoilDensity=2650.0#[kg/m^3]
    WDensity=1000.0#[kg/m^3]
    hcyl=6.0#[cm]
    dcyl=8.0#[cm]
    Vcyl=hcyl*3.1415*dcyl*dcyl/4.0*1.0e-6# [m^3]
    mcyl=45.0#[g]
    mwet=445.0-mcyl#[g]
    mdry=375.0-mcyl#[g]

    mwater=mwet-mdry#water mass
    gravWaterContent=mwater/mdry
    BulkDensity=mdry*1.0e-3/Vcyl#[kg/m^3]
    volWaterContent=gravWaterContent*BulkDensity/WDensity
    porosity=1.0-BulkDensity/SoilDensity

    print("w=",gravWaterContent)
    print("Bulk density =",BulkDensity)
    print("theta =",volWaterContent)
    print("Porosity =",porosity)
    
#    satMassWetness=computeSaturationWetness(bulkDensity)
#    if (gravWaterContent >= 0) and (gravWaterContent < satMassWetness):
#        computePorosity(bulkDensity, gravWaterContent)
#    else:
#        print("Wrong water content! value at saturation =", satMassWetness)
        

main()