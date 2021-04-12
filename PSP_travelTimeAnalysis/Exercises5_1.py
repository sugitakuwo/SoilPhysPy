#PSP Exercises 5.1

from sys import exit
import numpy as np
import matplotlib.pyplot as plt


#--- Import the data --------------#
data1="watercontent_bulkdensity_temp20degC_ers4.txt"
data2="watercontent_ers_temp20degC_bdensity1350.txt"
data3="watercontent_temp_ers4_bdensity1350.txt"

data1=np.loadtxt(data1,float)
data2=np.loadtxt(data2,float)
data3=np.loadtxt(data3,float)

#-------plot data in Figure------------

#Setting rcParams ---------------
plt.rcParams["figure.figsize"]=(6.5,5.0)
plt.rcParams["font.family"]="sans-serif"
plt.rcParams["xtick.direction"]="in"
plt.rcParams["ytick.direction"]="in"
plt.rcParams["xtick.major.width"]=1.0
plt.rcParams["ytick.major.width"]=1.0
plt.rcParams["font.size"]=14
plt.rcParams["axes.linewidth"]=1.1
plt.rcParams["errorbar.capsize"]=3.0
#-------------------------------------
plt.rcParams["legend.fancybox"]=False# rounded corner
plt.rcParams["legend.framealpha"]=0#Transparency
plt.rcParams["legend.edgecolor"]='white'#edge color
#plt.rcParams["legend.handlelength"]=1.0#line length
#plt.rcParams["legend.labelspacing"]=1.1#spacing between legend
#plt.rcParams["legend.handletextpad"]=1.0#space between line and text
#plt.rcParams["legend.markerscale"]=1.0#marked scale
#-------------------------------------

#--- Water content vs bulkdensity ---#
plt.plot(data1[:,0],data1[:,1],"ko--",label="Bulk density")
plt.xlabel("Bulk density [kg/m$^{3}$]",fontsize=14)
plt.ylabel("Water content [-]",fontsize=14)
#plt.xlim(0.04,200.0)
plt.ylim(0.28,0.33)
#plt.xscale('log')
#plt.yscale('log')
plt.legend(loc="best")

plt.show()

#--- Water content vs dielectric constant for solid phase ---#
plt.plot(data2[:,0],data2[:,1],"b^--",label="Solid dielectric constant")
plt.xlabel("Solid dielectric constant [-]",fontsize=14)
plt.ylabel("Water content [-]",fontsize=14)
#plt.xlim(0.04,200.0)
plt.ylim(0.26,0.37)
#plt.xscale('log')
#plt.yscale('log')
plt.legend(loc="best")

plt.show()

#--- Water content vs temperature ---#
plt.plot(data3[:,0],data3[:,1],"rD--",label="Temperature")
plt.xlabel("Temperature [degC]",fontsize=14)
plt.ylabel("Water content [-]",fontsize=14)
#plt.xlim(0.04,200.0)
plt.ylim(0.28,0.33)
#plt.xscale('log')
#plt.yscale('log')
plt.legend(loc="best")

plt.show()