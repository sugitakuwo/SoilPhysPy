#PSP_TTwaterContent.py
from __future__ import division
from math import sqrt

c = 299792458                       
airPermittivity = 1.00058986        

def getLiquidPermittivity(temperature):
    deltaT = temperature - 25.
    return(78.54 * (1-4.579E-03 * deltaT))
    
def getBulkPermittivity(probleLenght, travelTime, Vp):
    return(((c * Vp * travelTime) / (2. * probleLenght))**2)

def getWaterContentTopp(bulkPermittivity):
    return(-5.3E-02 + 2.92E-02 * bulkPermittivity - 5.5E-04 * bulkPermittivity**2
           + 4.3E-06 * bulkPermittivity**3)
    
def getWaterContentMalicki(bulkPermittivity, bulkDensity):
    bulkDensity /= 1000.
    return((sqrt(bulkPermittivity) - 0.819 - 0.168*bulkDensity - 0.159*bulkDensity**2)
            / (7.17 + 1.18*bulkDensity))
    
def getWaterContentMixModel(bulkPermittivity, bulkDensity, 
                            solidPermittivity, liquidPermittivity, alpha):
    porosity = 1. - bulkDensity/2650.
    numerator = bulkPermittivity**alpha - ((1. - porosity) * solidPermittivity**alpha 
                                           + porosity * airPermittivity**alpha)
    denominator =  liquidPermittivity**alpha - airPermittivity**alpha 
    return(numerator/denominator)