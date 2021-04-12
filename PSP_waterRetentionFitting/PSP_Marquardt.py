#PSP_Marquardt.py
from __future__ import print_function, division
import numpy as np
from math import sqrt
from PSP_waterRetention import *

EPSILON = 0.00001
MAX_ITERATIONS_NR = 100

def Marquardt(waterRetentionCurve, v0, vmin, vmax, x, y):
    n = len(v0)
    Lambda0 = 0.01
    vFactor = 2.
    l = np.array([Lambda0]*n)    # damping parameters
    v = np.zeros(n, float)
    for i in range(n): v[i] = v0[i]
    
    nrIter = 1
    maxDiff = 1.0
    sse = norm(waterRetentionCurve, v0, x, y)
    
    while (maxDiff > EPSILON) and (nrIter < MAX_ITERATIONS_NR):  
        diff = LeastSquares(waterRetentionCurve, l, v, vmin, vmax, x, y)
        maxDiff = max(abs(diff))
        v_new = computeNewParameters(v, vmin, vmax, diff, l, vFactor)
        sse_new = norm(waterRetentionCurve, v_new, x, y)
        
        if (sse_new < sse):
            sse = sse_new
            for i in range(n): v[i] = v_new[i]
            l /= vFactor
        else:
            l *= vFactor
        nrIter += 1
    print ("iterations nr:", nrIter)
    print ("sum of squared residuals:", sse)  
    return(v)


def estimate(waterRetentionCurve, v, Psi):
    waterContent = np.zeros(len(Psi))
    if (waterRetentionCurve == CAMPBELL):
        Campbell(v, Psi, waterContent)
    elif (waterRetentionCurve == VAN_GENUCHTEN):
        VanGenuchten(v, Psi, waterContent)
    elif (waterRetentionCurve == RESTRICTED_VG):
        VanGenuchtenRestricted(v, Psi, waterContent)
    elif (waterRetentionCurve == IPPISCH_VG):
         IppischVanGenuchten(v, Psi, waterContent)
    elif (waterRetentionCurve == CAMPBELL_IPPISCH_VG):
         CampbellIppischVanGenuchten(v, Psi, waterContent)

    return(waterContent)    


def computeNewParameters(v, vmin, vmax, diff, l, factor):
    n = len(v)
    v_new = np.zeros(n, float)
    for i in range(n):
        v_new[i] = v[i] + diff[i]
        if (v_new[i] > vmax[i]):
            v_new[i] = vmax[i]
            l[i] *= factor
        if (v_new[i] < vmin[i]):
            v_new[i] = vmin[i]
            l[i] *= factor  
    return(v_new)


def norm(waterRetentionCurve, v, x, y):
    yEst = estimate(waterRetentionCurve, v, x)
    norm = 0
    for i in range(len(x)):
        dy = y[i] - yEst[i]
        norm += (dy*dy)
    return(norm)


def LeastSquares(waterRetentionCurve, l, v, vmin, vmax, x, y):
    n = len(v)
    m = len(x)
    p = np.resize(np.zeros(n, float),(n,m))
    a = np.resize(np.zeros(n, float),(n,n))
    z = np.zeros(n, float)
    g = np.zeros(n, float)
    v1 = np.zeros(n, float)
    diff = np.zeros(n, float)
    
    for i in range(n): v1[i] = v[i]
    est = estimate(waterRetentionCurve, v, x)                           
       
    for i in range(n):
        change = (vmax[i] - vmin[i]) * 0.01
        v1[i] +=  change  
        # get a new set of estimates                              
        yEst = estimate(waterRetentionCurve, v1, x)                      
        v1[i] -= change                                 
        for j in range(m):
            # compute derivatives
            p[i][j] = (yEst[j] - est[j]) / change       
    
    for i in range(n):
        for j in range(i, n):
            a[i][j] = 0
            for k in range(m):
                a[i][j] = a[i][j] + p[i][k] * p[j][k]
        z[i] = sqrt(a[i][i]) + EPSILON

    for i in range(n):
        g[i] = 0
        for k in range(m):
            g[i] = g[i] + p[i][k] * (y[k] - est[k])
        g[i] = g[i] / z[i]
        for j in range(i, n):
            a[i][j] = a[i][j]  / (z[i] * z[j])
       
    for i in range(n):
        a[i][i] = a[i][i] + l[i]
        for j in range(i+1, n):
            a[j][i] = a[i][j]

    for j in range(n-1):
        pivot = a[j][j]
        for i in range(j+1, n):
            mult = a[i][j] / pivot
            for k in range(j+1, n): a[i][k] -= mult * a[j][k]
            g[i] -=  mult * g[j]

    diff[n-1] = g[n-1] / a[n-1][n-1]
     
    for i in range(n-2, -1, -1):
        top = g[i]
        for k in range(i+1, n):
            top -= a[i][k] * diff[k]
        diff[i] = top / a[i][i]
    
    for i in range(n): 
        diff[i] /=  z[i]

    return(diff)