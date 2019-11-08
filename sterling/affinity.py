#!/usr/bin/env python
from __future__ import division
import numpy as np
import ratlib as rat




def get_affinity(M,M_h,x_bj,z_h,eta,Q,q_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf):
	R1 = rat.get_R1( M,M_h,x_bj,z_h,eta,Q,q_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf)
	R2 = rat.get_R2( M,M_h,x_bj,z_h,eta,Q,q_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf)
	#inside=np.sum([1 for r1 in R1 if r1 < 0.019])
	N = 1000
	inside = np.sum([1 for i in range(N) if R1[i]<0.3 if R2[i]<.3]) #.019
	affinity = inside/N
	#print affinity
	#print inside
	return affinity


if __name__=="__main__":
	#PS
	#x_bj = .2
	#Q = 3
	#z_h = .3
	#P_t = 0.1
	#q_t
	#q_t = P_t/z_h

	#eta = 0.1 #get_eta(x_bj,Q,z_h,P_t) #define eta in ratlib called get_eta

	M= .94
	M_h=.14
	eta = rat.get_eta( M,M_h,x_bj,z_h,Q,q_t)
	
	xi = np.random.uniform(0.3,0.4,N)
	zeta = np.random.uniform(0.3,0.4,N)
	delta_k_t = np.random.uniform(0.3,0.4,N)
	k_i_t = np.random.uniform(0,0.1,N)
	M_ki = np.random.uniform(0,0.1,N)
	M_kf = np.random.uniform(0,0.1,N)
	affinity = get_affinity(M,M_h,x_bj,z_h,eta,Q,q_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf)
	print(affinity)


























