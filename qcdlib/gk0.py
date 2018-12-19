#!/usr/bin/env python
import sys
import os
import numpy as np
import time
from qcdlib.core import CORE
from qcdlib.interpolator import INTERPOLATOR
from tools.config import conf

class GK:
    """
    GK function. We will use a simple form based on large bT dominance where Ktilde is constant
    """

    def __init__(self):
        self.aux = conf['aux']
        self.set_default_params()
#        self.setup()
 
    def set_default_params(self):
        # free params
        self.gk  = 0.1
        self.Q0  = 1.
 
    def setup(self):
        self.lam  = 0.25
        self.Nf = 3
        self.beta0 = (33.-2.*self.Nf)/(12.*np.pi)
        self.CF = 4./3.

        pass
#        self.gk = self._gk
#        self.Q0 = self._Q0

    def get_gk(self,Q2):
        return (Q2/self.Q0**2)**self.gk

    def get_state(self):
        return (self.gk,self.Q0)

    def set_state(self, state):
        self.gk = state[0]
        self.Q0 = state[1]


# Alpha_S at the first loop, Lambda = 0.25
    def alpha_s(self,Q):
        return 1./(beta0*np.log(Q**2/self.lam**2))

# Perturbative Sudakov
    def get_pertsud(self,Q2):
        Q=np.sqrt(Q2)
        x = np.log(np.log(Q/self.lam)/np.log(self.Q0/self.lam))
        return -self.CF/(2.*np.pi*self.beta0)*(-4.*np.log(self.Q0)+(3.+4.*np.log(self.lam))*x
                                          +np.log(Q)*(4.-4.*x))

if __name__ == '__main__':

    from qcdlib.aux import AUX

    conf['aux']  = AUX()
    conf['gk']  = GK()

 
    Q2 = 2.4
    print conf['gk'].get_gk(Q2)
















