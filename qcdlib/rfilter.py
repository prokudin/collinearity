# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 21:36:05 2017

@author: prokudin
"""
#!/usr/bin/env python

import numpy as np
#############################################


class Rfilter(object):

    def __init__(self, hadron='pi+', fudge=[0, 0]):
        self.M = 0.938
        self.M2 = self.M**2
        self.set_Mh(hadron)
        

        self.fudge1 = fudge[0]
        self.fudge2 = fudge[1]

        self.MiT2 = 0.5 + 0.3 * fudge[0]
        self.MfT2 = 0.5 + 0.3 * fudge[1]
        self.kT = 0.3

        self.MX = 1.3
        self.Ma = 1.5
        self.Mb = 0.3
        self.deltaM = 0.3
        self.MJ = 0.3
        
        self.MiT = self.MiT2**0.5
        self.MfT = self.MfT2**0.5
        self.kT2 = self.kT**2
        
        # from paper, we need kf2, deltakT2, ki2
        # we need to have a way to vary them!
        self.kf2      = 0.3**2
        self.ki2      = 0.3**2
        self.deltakT  = 0.3 
        self.deltakT2 = self.deltakT**2

    def set_Mh(self, hadron):
        if hadron == 'pi+':
            self.Mh = 0.135  #  
        if hadron == 'pi-':
            self.Mh = 0.135  #  
        if hadron == 'pi0':
            self.Mh = 0.135  #  
        if hadron == 'k+':
            self.Mh = 0.493
        if hadron == 'k-':
            self.Mh = 0.493
        if hadron == 'k0':
            self.Mh = 0.493
        self.Mh2 = self.Mh**2

    def get_W2(self, x, Q2):
        W2 = Q2 * (1. - x) / x + self.M2
        return W2

    def get_MBT(self, PhT):
        MBT = np.sqrt( self.Mh2 + PhT**2 )
        return MBT
        
    def get_MiT(self, x, Q2):
        xn = self.get_xn(x, Q2)
#    self.MiT = np.sqrt( (xn*self.kT**2 + xn*self.MX**2 - (1-xn)*xn*self.M2)/(1-xn) )
        self.MiT = np.sqrt((xn * self.kT**2 + xn * (self.Ma + self.Mb / np.sqrt(xn))
                            ** 2 - (1 - xn) * xn * self.M2 + self.deltaM**2 * self.fudge1) / (1 - xn))
        return self.MiT

    def get_MfT(self):
        self.MfT = np.sqrt(self.kT**2 + self.MJ**2 +
                           self.deltaM**2 * self.fudge2)
        return self.MfT

    def get_xn(self, x, Q2):
        return 2 * x / (1 + np.sqrt(1 + 4 * x**2 * self.M2 / Q2))

    def get_yh(self, x, z, Q2, PhT, hadron, sign=-1):
        xn = self.get_xn(x, Q2)
        self.set_Mh(hadron)
        expy = Q2**0.5 * z * (Q2 - xn**2 * self.M2)\
            / (2 * self.M2 * xn**2 * (self.Mh2 + PhT**2)**0.5)\
            + sign * Q2**0.5 / (xn * self.M) * (z**2 * (Q2 - xn**2 * self.M2)**2
                                        / (4 * self.M2 * xn**2 * (self.Mh2 + PhT**2)) - 1)**0.5
        return np.log(expy)


    # zn
    def get_zn(self, x, z, Q2, PhT, hadron):
        xn = self.get_xn(x,Q2)
        self.set_Mh(hadron)
        MBT = self.get_MBT(PhT)
        zn =  xn * z/ (2 * x) * ( 1. + np.sqrt( 1.- 4. * self.M2 * MBT**2 * x**2 /  ( Q2**2 * z**2) ) )
        return zn
        
# rapidity of the target

    def get_yp(self, x, Q2):
        xn = self.get_xn(x, Q2)
        return np.log(np.sqrt(Q2) / (xn * self.M))

    def get_yi(self, Q2):
        return 0.5 * np.log(Q2 / self.MiT**2)

    def get_yf(self, Q2):
        return -0.5 * np.log(Q2 / self.MfT**2)

    def get_MhT(self, PhT):
        return np.sqrt(self.Mh2 + PhT**2)

    def get_R(self, x, z, Q2, PhT, hadron):
        self.set_Mh(hadron)
        MfT = self.get_MfT()
        MiT = self.get_MiT(x, Q2)
        yi = self.get_yi(Q2)
        yf = self.get_yf(Q2)
        MhT = self.get_MhT(PhT)
        yh = self.get_yh(x, z, Q2, PhT, hadron)
        zn = self.get_zn(x, z, Q2, PhT, hadron)
        znhat = zn/z
        # from paper..
        Ph_kf = 0.5 * MhT * MfT * (np.exp(yf - yh) + np.exp(yh - yf)) - znhat/zn * PhT**2 - PhT * self.kT
        Ph_ki = 0.5 * MhT * MiT * (np.exp(yi - yh) - np.exp(yh - yi)) - PhT * self.kT
        return np.abs(Ph_kf / Ph_ki)

        
# R0, Eq. (4.14)...
    def get_R0(self, Q2):
        """
        Collinearity ratio defined in the paper Eq. (4.15)
        """
        return  np.maximum(np.maximum(self.ki2/Q2,self.kf2/Q2),self.kT2/Q2)  
        
        
# We call R1 in the new paper what was R in the previous...
    def get_R1(self, x, z, Q2, PhT, hadron):
        """
        Collinearity ratio defined in the paper Eq. (4.15)
        """
        return  self.get_R( x, z, Q2, PhT, hadron)  

# R2 from Eq. (4.17)        
    def get_R2(self, x, z, Q2, PhT, hadron):
        self.set_Mh(hadron)
        zn = self.get_zn(x, z, Q2, PhT, hadron)
        znhat = zn/z
        qT = -PhT/zn
        # from paper..
        return np.abs( -(1.-znhat) - znhat*qT**2/Q2 - (1.-znhat)*self.kf2/(znhat*Q2) -self.deltakT2/(znhat*Q2) + 2.*qT*self.deltakT/Q2 )
    
#    def get_R0(self, x, z, Q2, PhT, hadron):
#        """
#        General hardness ratio defined in the paper Eq. (4.14)
#        """
#        return  ....  
    

if __name__ == '__main__':

    kin = Rfilter(fudge=[0, 0])

    x = 0.01
    z = 0.5
    Q2 = 4.0
    PhT = 0.1
    R0 = kin.get_R0(Q2)
    R1 = kin.get_R1(x, z, Q2, PhT, 'pi+')
    R2 = kin.get_R2(x, z, Q2, PhT, 'pi+')
    print('R0=', R0, 'R1=', R1, 'R2=', R2)
