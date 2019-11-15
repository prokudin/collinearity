#!/usr/bin/env python
from __future__ import division
import numpy as np
import ratlib as rat
from numpy import sqrt,log,exp
import affinity as af






if __name__=="__main__":


    M= .94
    M_h=.14
    x_bj = .2
    Q = 3
    z_h = .3
    P_t = 0.1

    #N= 1000
    xi = .3 #np.random.uniform(0.3,0.4,N)
    zeta = .3 #np.random.uniform(0.3,0.4,N)
    delta_k_t = .3 #np.random.uniform(0.3,0.4,N)
    k_i_t = .2 #np.random.uniform(0,0.1,N)
    M_ki = .1 #np.random.uniform(0,0.1,N)
    M_kf = .1 #np.random.uniform(0,0.1,N)
    #M= .938
    #M_h=.135
    #eta = rat.get_eta( M,M_h,x_bj,z_h,Q,q_t)

    #xi = np.random.uniform(0.3,0.4,N)
    #zeta = np.random.uniform(0.3,0.4,N)
    #delta_k_t = np.random.uniform(0.3,0.4,N)
    #k_i_t = np.random.uniform(0,0.1,N)
    #M_ki = np.random.uniform(0,0.1,N)
    #M_kf = np.random.uniform(0,0.1,N)
   
    
    q_t = -P_t/z_h
    eta = .1


    R1 = rat.get_R1( M,M_h,x_bj,z_h,eta,Q,q_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf)
    R2 = rat.get_R2( M,M_h,x_bj,z_h,eta,Q,q_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf)
    



    
    #print(len(R1))
    print(R1)
    print(R2)

    
