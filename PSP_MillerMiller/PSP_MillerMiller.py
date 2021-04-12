#PSP_MillerMiller.py
import numpy as np
import matplotlib.pyplot as plt
from PSP_soil import *

GAUSSIAN = 1
EXPONENTIAL = 2

def statisticalFunction(x, y, funcType):
    if (funcType == GAUSSIAN): 
        return np.exp(-0.25*np.pi*(x*x + y*y))           
    elif (funcType == EXPONENTIAL): 
        return np.exp(-np.sqrt(x*x + y*y))             

def heterogeneousField(n, dx, dy, lx, ly, funcType):         
    f = np.array([[complex(0.,0.)]*int(n+1)]*int(n+1))
    phi = np.zeros((n+1)**2)
    for i in range((n+1)**2): 
        phi[i]=2*np.pi*np.random.random_sample()      
    for i in range (int(n/2)+1):                                              
        for j in range (int(n/2)+1):
            r = statisticalFunction(i*dx/lx, j*dy/ly, funcType)
            f[-i +int(n/2)][-j+int(n/2)] = complex(r,0.)
            f[+i +int(n/2)][-j+int(n/2)] = complex(r,0.)
            f[-i +int(n/2)][+j+int(n/2)] = complex(r,0.)
            f[+i +int(n/2)][+j+int(n/2)] = complex(r,0.)
    fft_f = np.fft.fft2(f);                                          
    for i in range(n+1):                                                       
        for j in range(n+1):
            r = pow(2.*np.real(fft_f[i][j])**2,0.25);
            fft_f[i][j] = complex(r*np.cos(phi[i*(n+1)+j]), r*np.sin(phi[i*(n+1)+j]))
    f  = np.real(np.fft.ifft2(fft_f))                                          
    mean_f = np.mean(f)                                                         
    std_f = np.std(f, dtype=np.float64)
    f  = (f - mean_f) / std_f
    return(f)

def plotHydraulicProperties(waterRetentionCurve, sigma):  
    soil = Csoil()
    soil_minus_sigma = Csoil()
    soil_plus_sigma = Csoil()
    soil_minus_sigma.Campbell_he /= np.exp(-sigma*sigma)
    soil_minus_sigma.VG_alpha *= np.exp(-sigma*sigma)
    soil_minus_sigma.Ks *= np.exp(-2*sigma*sigma)
    soil_plus_sigma.Campbell_he /= np.exp(sigma*sigma)
    soil_plus_sigma.VG_alpha *= np.exp(sigma*sigma)
    soil_plus_sigma.Ks *= np.exp(2*sigma*sigma)
    
    n = 100
    psi = -np.logspace(-2, 5, n)
    theta = np.zeros(n); theta_minus = np.zeros(n); theta_plus = np.zeros(n)
    K = np.zeros(n); K_minus = np.zeros(n); K_plus = np.zeros(n)
    for i in range(n):
        theta[i] = waterContent(waterRetentionCurve, soil, psi[i])
        K[i] = hydraulicConductivity(waterRetentionCurve, soil, psi[i])
        theta_minus[i] = waterContent(waterRetentionCurve, soil_minus_sigma, psi[i])
        K_minus[i] = hydraulicConductivity(waterRetentionCurve, soil_minus_sigma, psi[i])
        theta_plus[i] = waterContent(waterRetentionCurve, soil_plus_sigma, psi[i])
        K_plus[i] = hydraulicConductivity(waterRetentionCurve, soil_plus_sigma, psi[i])

    plt.subplot(1,2,1)
    plt.plot(-psi, theta_minus, '--k', label='-$\sigma$')
    plt.plot(-psi, theta, '-k', label='$\mu$')
    plt.plot(-psi, theta_plus, ':k', label='+$\sigma$')
    plt.xlabel('$\psi$ [J kg$^{-1}$]',labelpad=8,fontsize=20)
    plt.ylabel('$\\theta$ [m$^3$ m$^{-3}$]',labelpad=8,fontsize=20)
    plt.xscale('log'), plt.legend()
    plt.tick_params(axis='both', which='major', labelsize=14,pad=8)

    
    plt.subplot(1,2,2)
    plt.plot(-psi, K_minus*3600., '--k', label='-$\sigma$')
    plt.plot(-psi, K*3600., '-k',label='$\mu$')
    plt.plot(-psi, K_plus*3600., ':k', label='+$\sigma$')
    plt.xscale('log'), plt.yscale('log')
    plt.xlabel('$\psi$ [J kg$^{-1}$]',labelpad=8,fontsize=20)
    plt.ylabel('K [kg s$^{-1}$ m$^{-3}$]',labelpad=8,fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=14,pad=8)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plotThetaField(waterRetentionCurve, f, n):
    soil = Csoil() 
    currentSoil = Csoil()
    theta = np.resize(np.zeros(((n+1)**2)*3),(3,n+1,n+1))
    psi = [-0., -15., -100.]
    for i in range(len(psi)):
        for k in range(n+1):
            for l in range(n+1):
                currentSoil.Campbell_he = soil.Campbell_he / np.exp(f[k,l])
                currentSoil.VG_alpha = soil.VG_alpha * np.exp(f[k,l])
                theta[i,k,l] = waterContent(waterRetentionCurve, currentSoil, psi[i])
                
        plt.subplot(1, len(psi), i+1),
        plt.imshow(theta[i], cmap=plt.gray(), vmin = 0.0, vmax = 0.5)
        plt.axis('off'), plt.title(r'$\theta$ ($\psi$='+str(psi[i])+')')
        cbar = plt.colorbar(orientation='horizontal', ticks = [0.0,0.5],
                                    shrink=0.8, aspect = 7., pad = 0.05)   
        cbar.set_ticklabels(['%.1f' %(0.0), '%.1f' %(0.5)])
    plt.show()


def plotConductivityField(waterRetentionCurve, f, n):
    soil = Csoil()
    currentSoil = Csoil()
    K = np.resize(np.zeros(((n+1)**2)*3),(3,n+1,n+1))
    psi = [-0., -15., -100.]
    for i in range(len(psi)):
        for k in range(n+1):
            for l in range(n+1):
                currentSoil.Campbell_he = soil.Campbell_he / np.exp(f[k,l])
                currentSoil.VG_alpha = soil.VG_alpha * np.exp(f[k,l])
                currentSoil.Ks = soil.Ks * np.exp(2*f[k,l])
                K[i,k,l] = hydraulicConductivity(waterRetentionCurve, currentSoil, psi[i])
             
        plt.subplot(1, len(psi), i+1),
        plt.imshow(np.log(K[i])/np.log(10.), cmap=plt.gray(), vmin = -10, vmax = -1)
        plt.axis('off'), plt.title(r'$\log_{10}$K ($\psi$=' + str(psi[i])+')')
        cbar = plt.colorbar(orientation='horizontal', ticks = [-10,-8,-6,-4,-2],
                                        shrink=0.8, aspect = 7., pad = 0.05)
        cbar.set_ticklabels([-10,-8,-6,-4,-2])
    plt.show()