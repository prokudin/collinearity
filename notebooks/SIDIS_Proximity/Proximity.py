
"""
Created on Sun Dec 29 18:29:02 2019

@author: scott
"""

from ratios import ratios
import numpy as np


def Proximity(data,domain1,domain2,domainPhi,p,N):
    """
    data    = any Pandas DataFrame should do as long as it contains x,z,Q2, and pT
    domain1 = 
    domain2 = 
    p       = radius of sphere representing the demarcation of the proximity
    hadron  = any hadron found in the data
    """
    
    for i in range(len(data)):
        
        # Creates vectors of length N made from uniform random numbers in defined domains.
        zeta = np.random.uniform(data.loc[i,'z']-.1,data.loc[i,'z']+.1,N)
        ki = np.random.uniform(domain1[0],domain1[1],N)
        kf = np.random.uniform(domain1[0],domain1[1],N)
        kT = np.random.uniform(domain1[0],domain1[1],N)
        kibT = np.random.uniform(domain1[0],domain1[1],N)
        deltakT = np.random.uniform(domain2[0],domain2[1],N)
        xi = np.linspace(domain2[0],domain2[1],N)
        phi = np.linspace(domainPhi[0],domainPhi[1],N)
        
        if data.loc[i,'hadron'] == 'k-' or data.loc[i,'hadron'] == 'k+':
            Mh = 0.493
        else:
            Mh = 0.135
        M = 0.94
        
        R0,R1,R2 =  ratios(data.loc[i,'x'], data.loc[i,'z'], data.loc[i,'Q2'], data.loc[i,'pT'],Mh,M,zeta,xi,ki,kf,kT,kibT,deltakT,phi)
        
        data.loc[i,'proximity'] = sum(1 for i in range(N) if (R0[i]**2+R1[i]**2+R2[i]**2)<p)/N
        
        
    return data