#!/usr/bin/env python
from __future__ import division
import numpy as np
import ratlib as rat
import pandas as pd
import affinity as af
import matplotlib.pyplot as plt
tab = pd.read_excel("expdata/expdata/1008.xlsx")
npts = tab["x"].values.size
#x_bj = .2
#Q = 3
# = .3
#P_t = 0.1
#q_t


#eta = 0.1 #get_eta(x_bj,Q,z_h,P_t) #define eta in ratlib called get_eta

M= .94
M_h=.14

N = 1000
xi = np.random.uniform(0.3,0.4,N)
zeta = np.random.uniform(0.3,0.4,N)
delta_k_t = np.random.uniform(0.3,0.4,N)
k_i_t = np.random.uniform(0,0.1,N)
M_ki = np.random.uniform(0,0.1,N)
M_kf = np.random.uniform(0,0.1,N)

for i in range(npts):
	x=tab["x"].values[i]
	z=tab["z"].values[i]
	Q=tab["Q2"].values[i]**0.5
	P_t=tab["pT"].values[i]

	#x_bj = x
	#z_h = z
	q_t = P_t/z
	eta = .1
	#eta = rat.get_eta( M,M_h,x,z,Q,q_t)
	affinity =af.get_affinity(M,M_h,x,z,eta,Q,q_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf)
	#print x,z,Q,P_t,affinity
	
	print affinity

	#plt.scatter(x,Q)
	#plt.show()




















