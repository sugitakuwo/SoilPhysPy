#PSP_soil.py

CAMPBELL = 1
VAN_GENUCHTEN = 2

class Csoil:
    def __init__(self):
        self.thetaS = 0.46              # [m^3 m^-3] saturated water content
        self.thetaR = 0.02              # [m^3 m^-3] residual water content
        self.Ks = 1.e-3                 # [kg s^-1 m^-3] saturated hydraulic conductivity
        self.Mualem_L = 0.5             # [-] tortuosity parameter
                                        # Campbell model:
        self.Campbell_he = -4.2         # [J kg^-1] air-entry potential 
        self.Campbell_b = 3.6           # [-] slope parameter 
                                        # van Genuchten modified model:
        self.VG_he = 3.8                # [J kg^-1] air-entry potential 
        self.VG_alpha = 0.29            # [kg J^1] related to the inverse of the air entry (>0) 
        self.VG_n = 1.33                # [-] measure of the pore-size distribution (>1)
        self.VG_m = 1.-(1./self.VG_n)   # [-] m parameter
        
def hydraulicConductivity(waterRetentionCurve, soil, signPsi):
    if (waterRetentionCurve == CAMPBELL):
        if (signPsi >= soil.Campbell_he):
            K = soil.Ks
        else:
            K = soil.Ks * (soil.Campbell_he / signPsi)**(2.+ 3./soil.Campbell_b)
    elif (waterRetentionCurve == VAN_GENUCHTEN):
        if (signPsi >= 0.):
            K = soil.Ks
        else:
            Se = 1. / (1. + (soil.VG_alpha * abs(signPsi))**soil.VG_n)**soil.VG_m
            K = soil.Ks * Se**soil.Mualem_L * (1.-(1.-Se**(1./soil.VG_m))**soil.VG_m)**2 
    return(K)

def waterContent(waterRetentionCurve, soil, signPsi):
    if (waterRetentionCurve == CAMPBELL):
        if (signPsi >= soil.Campbell_he):
            return soil.thetaS
        else:
            Se = (soil.Campbell_he / signPsi)**(1. / soil.Campbell_b)
            return Se * soil.thetaS
    elif(waterRetentionCurve == VAN_GENUCHTEN):
        if (signPsi >= 0.):
            return soil.thetaS
        else:
            Se = 1. / (1. + (soil.VG_alpha * abs(signPsi))**soil.VG_n)**soil.VG_m
            return Se*(soil.thetaS - soil.thetaR) + soil.thetaR