import sys
import os
import numpy as np
from mpmath import fp
from scipy.special import gamma
from scipy.integrate import quad


class AUX:

    def __init__(self):

        self.set_constants()
        self.set_masses()
        self.set_couplings()
        self.set_evolution()

    def set_constants(self):

        self.CA = 3.0
        self.CF = 4.0 / 3.0
        self.TR = 0.5
        self.TF = 0.5
        self.euler = fp.euler

        # soft factors
        self.soft = {}

    def set_masses(self):

        self.me = 0.000511
        self.mmu = 0.105658
        self.mtau = 1.77684
        self.mu = 0.055
        self.md = 0.055
        self.ms = 0.2
        self.mc = 1.51
        self.mb = 4.92
        self.mZ = 91.1876
        self.mW = 80.398
        self.M = 0.93891897
        self.Mpi = 0.135
        self.Mk = 0.497

        self.me2 = self.me**2
        self.mmu2 = self.mmu**2
        self.mtau2 = self.mtau**2
        self.mu2 = self.mu**2
        self.md2 = self.md**2
        self.ms2 = self.ms**2
        self.mc2 = self.mc**2
        self.mb2 = self.mb**2
        self.mZ2 = self.mZ**2
        self.mW2 = self.mW**2
        self.M2 = self.M**2
        self.Mpi2 = self.Mpi**2

    def set_couplings(self):

        self.c2w = self.mW2 / self.mZ2
        self.s2w = 1.0 - self.c2w
        self.s2wMZ = 0.23116
        self.alfa = 1 / 137.036
        self.alphaSMZ = 0.118

    def set_evolution(self):
        self.Q02=1.0
        self.lam2=0.4

    def p2n(self, p):
        return np.copy(p[[0,3,4,1,2,5,6,7,8,9,10]])

    def charge_conj(self,p):
        return np.copy(p[[0,2,1,4,3,6,5,8,7,10,9]])

    def pip2piz(self, pip):
        piz = np.copy(pip)
        u = pip[1]
        ub = pip[2]
        d = pip[3]
        db = pip[4]
        piz[1] = 0.5 * (u + ub)
        piz[2] = 0.5 * (u + ub)
        piz[3] = 0.5 * (d + db)
        piz[4] = 0.5 * (d + db)
        return piz




