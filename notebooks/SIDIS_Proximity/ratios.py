#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 18:26:56 2019

@author: scott
"""

from numpy import *

def ratios(x,z,Q2,pT,Mh,M,zeta,xi,ki2,kf,kT,kiT,deltakT,phi):
    
    R0 = (maximum(maximum(-ki2/Q2,kf**2/Q2),kT**2/Q2))
    
    R1 = 2*Q2*x*xi*z*(Mh**2+zeta**2*(deltakT**2+kf**2))*(sqrt((-4*M**2*x**2*(Mh**2+pT**2)+Q2**2*z**2) \
         /(Q2**2*z**2))+1)/(zeta*(Q2**2*xi**2*z**2*(sqrt((-4*M**2*x**2*(Mh**2+pT**2)+Q2**2*z**2) \
         /(Q2**2*z**2))+1)**2 - 4*Q2*kiT*pT*x*xi*z*(sqrt((-4*M**2*x**2*(Mh**2+pT**2)+Q2**2*z**2) \
         /(Q2**2*z**2))+1)*cos(phi)+4*x**2*(Mh**2+pT**2)*(ki2+kiT**2)))
    
    R2 = abs((Q2*z**2*(sqrt((-4.*M**2*x**2*(Mh**2+pT**2)+Q2**2*z**2)/(Q2**2*z**2))+1)**2-Q2*z*zeta \
         *(sqrt((4.*M**2*x**2+Q2)/Q2)+1)*(sqrt((-4.*M**2*x**2*(Mh**2+pT**2)+Q2**2*z**2)/(Q2**2*z**2))+1) \
         - pT**2*(sqrt((4.*M**2*x**2+Q2)/Q2)+1)**2 - zeta*(sqrt((4*M**2*x**2+Q2)/Q2)+1) \
         *(deltakT**2*zeta *(sqrt((4*M**2*x**2+Q2)/Q2)+1)+2*deltakT*pT*(sqrt((4.*M**2*x**2+Q2)/Q2)+1)  \
         +kf**2*(-z*(sqrt((-4*M**2*x**2*(Mh**2+pT**2)+Q2**2*z**2)/(Q2**2*z**2))+1) \
         +zeta*(sqrt((4*M**2*x**2+Q2)/Q2)+1))))/(Q2*z*zeta*(sqrt((4.*M**2*x**2+Q2)/Q2)+1) \
         *(sqrt((-4*M**2*x**2*(Mh**2+pT**2)+Q2**2*z**2)/(Q2**2*z**2))+1)))
    
    return R0,R1,R2