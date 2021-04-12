#PSP_MillerMiller
from __future__ import print_function, division
from PSP_MillerMiller import *
    
def main():
    choice = 0
    print (CAMPBELL,' Campbell')
    print (VAN_GENUCHTEN,' van Genuchten')
    while (choice < CAMPBELL) or (choice > VAN_GENUCHTEN):
        choice = float(input("Choose water retention curve: "))
        if (choice < CAMPBELL) or (choice > VAN_GENUCHTEN):
            print('wrong choice.')
    waterRetentionCurve = choice
    
    sigma = 1.0
    plotHydraulicProperties(waterRetentionCurve, sigma) 
    
    choice = 0
    print (GAUSSIAN,' Gaussian model ')
    print (EXPONENTIAL,' exponential model')
    while (choice < GAUSSIAN) or (choice > EXPONENTIAL):
        choice = float(input("Choose autocovariance model: "))
        if (choice < GAUSSIAN) or (choice > EXPONENTIAL):
            print('wrong choice.')
    funcType = choice
    
    n = 128                    
    dx = dy = 0.01              
    lx = ly = 0.1                                 
    f = heterogeneousField(n, dx, dy, lx, ly, funcType)
    f *= (sigma*sigma)
    plt.imshow(f,cmap=plt.gray()), plt.axis('off'), plt.title(''), plt.show()
    
    plotThetaField(waterRetentionCurve, f, n)
    plotConductivityField(waterRetentionCurve, f, n)
main()