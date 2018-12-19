import sys,os
import numpy as np
import pylab as py
from scipy.interpolate import RectBivariateSpline


class INTERPOLATOR:

    def __init__(self,fname):

        L=open('%s/qcdlib/tables/%s.dat'%(os.environ['JAM3D'],fname)).readlines()
        x=[float(val) for val in L[3].split()]
        Q=[float(val) for val in L[4].split()]
        iflav=[int(val) for val in L[5].split()]
        L=[l.split() for l in L[6:] if '---' not in l]
        data=[[float(val) for val in l] for l in L]
        data=np.transpose(data)

        nx=len(x)
        nQ=len(Q)
        nf=len(iflav)

        self.T={}

        for i in range(nf): 
            table=np.zeros((nx,nQ))
            cnt=-1
            for ix in range(nx):
                for iQ in range(nQ):
                    cnt+=1
                    table[ix,iQ]=data[i][cnt]
            self.T[iflav[i]]=RectBivariateSpline(x,Q,table,kx=3, ky=4)

        self.IFLAV=[21,2,-2,1,-1,3,-3,4,-4,5,-5]

    def _get_f(self,iflav,x,Q2):
        return self.T[iflav](x,Q2**0.5)[0,0]/x

    def get_f(self,x,Q2):
        return np.array([self._get_f(iflav,x,Q2) for iflav in self.IFLAV])


