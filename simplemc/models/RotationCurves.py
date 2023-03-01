
from simplemc.cosmo.paramDefs import Anfw_par, rs_par#, eta_par, phi0_par, phi1_par, phi2_par
from simplemc.barotropic_functions import *
import numpy as np
from scipy import interpolate


class RotationCurves():
    def __init__(self, varya= True, varyb= True):#,varyc=True):,varyd=True,varye=True):
        """
        Class to constrain rotational curves profiles,
            Here we assume a NFW profile
        Parameters
        ----------
        varya
        varyb

        Returns
        -------

        """
        ## Example used to fit a straigth line two parameters: a, b
        self.varya = varya
        self.varyb = varyb
        #self.varyc = varyc
        #self.varyd = varyd
        #self.varye = varye

        self.Anfw = Anfw_par.value
        self.rs = rs_par.value
        #self.eta = eta_par.value
        #self.phi1 = phi1_par.value
        #self.phi2 = phi2_par.value

        data_path = "simplemc/data/Blok_McGaugh_&_Rubin_(2001)/"
        data = np.loadtxt(data_path+'ESO3020120.dat')
        self.vecRp_data = np.array([row[1] for row in data])


    # my free params (parameters/priors see ParamDefs.py)
    def freeParameters(self):
        l = []
        if (self.varya): l.append(Anfw_par)
        if (self.varyb): l.append(rs_par)
        #if (self.varyc): l.append(eta_par)
       # if (self.varyd): l.append(phi1_par)
        #if (self.varye): l.append(phi2_par)
        return l

    def printFreeParameters(self):
        print("Free parameters and values currently accepted:")
        self.printParameters(self.freeParameters())


    def printParameters(self, params):
        for p in params:
            print(p.name, '=' , p.value , '+/-' , p.error)


    def updateParams(self, pars):
        x = self.vecRp_data
        for p in pars:
            if p.name == "Anfw":
                self.Anfw = p.value
            elif p.name == "rs":
                self.rs = p.value
            #elif p.name == "eta":
            #    self.eta = p.value
            #elif p.name == "phi1":
            #    self.phi1 = p.value
            #elif p.name == "phi2":
            #    self.phi2 = p.value
        A = self.Anfw
        rs = self.rs
        #eta = self.eta
        #phi1 = self.phi1
        #phi2 = self.phi2
        kappa = 10.**(A)
        rho0 = rs
        #eta = eta
        #phi1 = 10.**(phi1)
        #phi2 = 10.**(phi2)
        #Vc = Vc_(x,m_a,eps,eta)
        #self.f = interpolate.interp1d(Vc[0],Vc[1],fill_value='extrapolate')
        Vc = Vc_(kappa,rho0)
        self.f = interpolate.interp1d(Vc[0],Vc[1],axis=0, fill_value="extrapolate")
        return True


    #def rotation(self,x):
        #def nfw(r, rs =1, A=0.05):
    #    A  = self.Anfw
    #    rs = self.rs
    #    return (A**2*(rs**3)/x)*(np.log((rs+x)/rs) - x/(rs+x))


    def rotation(self,x):
        #A = self.Anfw
        #rs = self.rs
        #phi0 = self.phi0
        #phi1 = self.phi1
        #phi2 = self.phi2

        #m_a = 10.**(A)
        #eps0 = 10.**(rs)
        #phi0 = 10.**(phi0)
        #phi1 = 10.**(phi1)
        #phi2 = 10.**(phi2)
        #Vc2 = Vc_inter(x,m_a,eps0,phi0,phi1,phi2)
        Vc_new = self.f(x)
        return np.sqrt(Vc_new)

    def prior_loglike(self):
        return 0 
