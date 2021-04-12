#PSP_gasDiffusion
from __future__ import print_function, division
from math import exp
from PSP_ThomasAlgorithm import ThomasBoundaryCondition
import PSP_grid as grid
import matplotlib.pyplot as plt
import numpy as np

def gasSolver(boundaryLayerCond, boundaryOxygenConc, dg, respRate, totalDepth, n):
    a  = np.zeros(n+2, float)  
    b  = np.zeros(n+2, float)  
    c  = np.zeros(n+2, float)  
    d  = np.zeros(n+2, float) 
    g  = np.zeros(n+2, float) 
    u  = np.zeros(n+2, float)  
    co = np.zeros(n+2, float)  
    
    g[0] = boundaryLayerCond
    co[0] = boundaryOxygenConc
    # vector depth [m]
    z = grid.linear(n, totalDepth)
    
    # initialize matrix
    for i in range(1, n+1):
        u[i] = respRate * exp(-z[i] / 0.3) * (z[i + 1] - z[i - 1]) / 2.0
        if i < n:
            g[i] = dg / (z[i + 1] - z[i])
        else:
            g[i] = 0
        a[i + 1] = -g[i]
        b[i] = g[i - 1] + g[i]
        c[i] = -g[i]
        d[i] = u[i]

    d[1] = d[1] + g[0] * co[0]
    
    ThomasBoundaryCondition(a, b, c, d, co, 1, n)
    
    return(z, co)


def main():
    R = 8.3143                     
    n = 20                        
    totalDepth = 1.0               
    bulkDensity = 1300.0#[1100.0,1300.0,1500.0,1700.0]#1300.0            
    particleDensity = 2650.         
    waterContent = [0.1,0.2,0.25,0.3]#0.2                 
    respRate = -0.001               
    oxygenDiff = 1.77e-5
    co2Diff = 1.39e-5           
    temperature = 25.             
    atmPressure = 101.3           
    boundaryLayerCond = 0.01#[0.00001,0.0001,0.01,1.0,100.0]#0.01
    numpara=len(waterContent)
#    label1=["w=0.1","w=0.2","w=0.3","w=0.4"]     
    
    # O2 concentration in air [g/m^3]
    boundaryOxygenConc = (0.21 * atmPressure * 1000. * 32. / (R * (temperature + 273.15)))
    # CO2 concentration in air [g/m^3]
    boundaryCO2Conc = (0.0025 * atmPressure * 1000. * 44. / (R * (temperature + 273.15))) 
    # plot results
    fig = plt.figure(figsize=(10,8))
    for ii in range(numpara):
        porosity = 1. - bulkDensity / particleDensity
        gasPorosity = porosity - waterContent[ii]
        print ("gasPorosity =",gasPorosity)
    
    #  binary diffusion coefficient [m2/s]
        binaryDiffCoeffO2 = (oxygenDiff * (101.3 / atmPressure) * ((temperature + 273.15) / 273.15)**1.75)
        binaryDiffCoeffCO2 = (co2Diff * (101.3 / atmPressure) * ((temperature + 273.15) / 273.15)**1.75)
    
        bg = 0.9           
        mg = 2.3           
        dgo2 = binaryDiffCoeffO2 * bg * gasPorosity**mg
        dgco2 = binaryDiffCoeffCO2 * bg * gasPorosity**mg
    
        z, O2co = gasSolver(boundaryLayerCond, boundaryOxygenConc, dgo2, respRate, totalDepth, n)
        z, CO2co = gasSolver(boundaryLayerCond, boundaryCO2Conc, dgco2, -respRate, totalDepth, n)
        O2co=O2co*(R * (temperature + 273.15))/(atmPressure * 1000. * 32.)*100
        CO2co=CO2co*(R * (temperature + 273.15))/(atmPressure * 1000. * 44.)*100
      
#        print ("node   depth [m]   O2 conc. [%]")
#        for i in range(n + 2):
#            print ("%3d    %6.2f      %.2f" %(i, z[i], O2co[i]))

#        print ("node   depth [m]   CO2 conc. [%]")
#        for i in range(n + 2):
#            print ("%3d    %6.2f      %.2f" %(i, z[i], CO2co[i]))

#    for i in range(n+1):
        plt.plot(O2co[:n+1], -z[:n+1], 'ko')
        plt.plot(CO2co[:n+1], -z[:n+1], 'r^')
        
    plt.xlabel('Concentration [%]',fontsize=20,labelpad=8)
    plt.ylabel('Depth [m]',fontsize=20,labelpad=8)
    plt.legend(loc="best")
    plt.tick_params(axis='both', which='major', labelsize=20,pad=8)
    plt.tick_params(axis='both', which='minor', labelsize=20,pad=8)
    plt.show()

main()
