from __future__ import division
from numpy import sqrt,log,exp
import numpy as np
sq2 = np.sqrt(2)
get_xN=lambda M,M_h,x_bj,z_h,eta,Q,T_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf: (2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)))
get_zN=lambda M,M_h,x_bj,z_h,eta,Q,T_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf: (Q*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))/((Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)))
get_R1=lambda M,M_h,x_bj,z_h,eta,Q,T_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf: (2*Q**2*x_bj*zeta*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)*(-2*Q**3*T_t*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))**3*(Q*T_t*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))/(zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)) - delta_k_t)/(zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))**3*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)**3) + Q**2*(M_kf**2 + (Q*T_t*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))/(zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)) - delta_k_t)**2)*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))**2/((Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)**2) + Q**2*(M_h**2 + Q**2*T_t**2*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))**2/((Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)**2))*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))**2/(zeta**2*(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)**2))/(xi*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))*(Q**6*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))**2/((Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)**2) + 4*Q**5*T_t*k_i_t*x_bj*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))**2/(xi*(Q + sqrt(4*M**2*x_bj**2 + Q**2))**3*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)**2) + 4*Q**2*x_bj**2*(M_h**2 + Q**2*T_t**2*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))**2/((Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)**2))*(-M_ki**2 + k_i_t**2)/(xi**2*(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2))))
get_R2=lambda M,M_h,x_bj,z_h,eta,Q,T_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf: np.abs((M_kf**2/Q**2 - M_kf**2*zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)/(Q**3*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))) + Q*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))/(zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)) - 1 - zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)*(Q*T_t*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))/(zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)) - delta_k_t)**2/(Q**3*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2)))))
get_R3=lambda M,M_h,x_bj,z_h,eta,Q,T_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf: (-1/2*xi*zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)*(2*Q**4*x_bj*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))*(Q*T_t*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))/(zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)) - delta_k_t + k_i_t)**2/(xi*zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)) - (Q**2*(Q*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))/(zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)) - 1) - 2*Q*x_bj*(-M_ki**2 + k_i_t**2)/(xi*(Q + sqrt(4*M**2*x_bj**2 + Q**2))))*(2*Q**4*x_bj*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))/(xi*zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)) - Q**3*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))/(zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)) + 2*Q*x_bj*(M_kf**2 + (Q*T_t*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))/(zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)) - delta_k_t)**2)/(xi*(Q + sqrt(4*M**2*x_bj**2 + Q**2)))))/(Q**6*x_bj*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))))
get_R4=lambda M,M_h,x_bj,z_h,eta,Q,T_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf: ((1/4)*(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)*(4*M**2*Q**2*x_bj**2*(M_h**2 + Q**2*T_t**2*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))**2/((Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)**2))/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**6*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))**2/((Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)**2))/(Q**6*x_bj*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))))
get_yp=lambda M,M_h,x_bj,z_h,eta,Q,T_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf: ((1/2)*log((1/4)*(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2/(M**2*x_bj**2)))
get_yh=lambda M,M_h,x_bj,z_h,eta,Q,T_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf: ((1/2)*log(2*((1/2)*M_h**2 + (1/2)*Q**2*T_t**2*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))**2/((Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)**2))*(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)**2/(Q**4*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))**2)))
get_yi=lambda M,M_h,x_bj,z_h,eta,Q,T_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf: (log((1/2)*xi*(Q + sqrt(4*M**2*x_bj**2 + Q**2))/(x_bj*sqrt(-M_ki**2 + k_i_t**2))))
get_yf=lambda M,M_h,x_bj,z_h,eta,Q,T_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf: (log(zeta*sqrt(M_kf**2 + (Q*T_t*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))/(zeta*(Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)) - delta_k_t)**2)*(Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)/(Q**2*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2)))))
get_W2=lambda M,M_h,x_bj,z_h,eta,Q,T_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf: (-1/2*(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)*(2*Q**6*T_t**2*x_bj*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))**3/((Q + sqrt(4*M**2*x_bj**2 + Q**2))**4*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)**3) + (2*M**2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) - Q**2*(Q*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))/((Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)) - 1))*(2*Q**4*x_bj*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))/((Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)) - Q**3*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))/((Q + sqrt(4*M**2*x_bj**2 + Q**2))*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)) + 2*Q*x_bj*(M_h**2 + Q**2*T_t**2*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))**2/((Q + sqrt(4*M**2*x_bj**2 + Q**2))**2*(4*M**2*Q**2*T_t**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + Q**4)**2))/(Q + sqrt(4*M**2*x_bj**2 + Q**2))))/(Q**4*x_bj*(Q**4*z_h + sqrt(-16*M**4*M_h**2*Q**2*T_t**2*x_bj**4/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*M_h**2*Q**4*x_bj**2 + Q**8*z_h**2))))
get_Htmax=lambda M,M_h,x_bj,z_h,eta,Q,T_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf: (sqrt(-(4*M**2*M_h*Q**2*x_bj**2*exp(2*eta)/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M**2*Q**3*x_bj**2*exp(eta)/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + 2*M_h**2*Q**2*x_bj*exp(eta)/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + 2*M_h*Q**3*x_bj*exp(2*eta)/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) - 2*M_h*Q**3*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + M_h*Q**2 - 2*Q**4*x_bj*exp(eta)/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + Q**3*exp(eta))*(4*M**2*M_h*Q**2*x_bj**2*exp(2*eta)/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + 4*M**2*Q**3*x_bj**2*exp(eta)/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 2*M_h**2*Q**2*x_bj*exp(eta)/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + 2*M_h*Q**3*x_bj*exp(2*eta)/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) - 2*M_h*Q**3*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + M_h*Q**2 + 2*Q**4*x_bj*exp(eta)/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) - Q**3*exp(eta)))/(4*M**2*Q**2*x_bj**2*exp(2*eta)/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + 2*Q**3*x_bj*exp(2*eta)/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) - 2*Q**3*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + Q**2))
get_etamin=lambda M,M_h,x_bj,z_h,eta,Q,T_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf: (log((1/4)*(Q + sqrt(4*M**2*x_bj**2 + Q**2))*(-4*M**2*Q**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + 2*M_h**2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) - 2*Q**3*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + Q**2 - sqrt((4*M**2*Q**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M*M_h*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + 2*M_h**2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + 2*Q**3*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) - Q**2)*(4*M**2*Q**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + 4*M*M_h*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + 2*M_h**2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + 2*Q**3*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) - Q**2)))/(M_h*x_bj*(2*M**2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + Q**2))))
get_etamax=lambda M,M_h,x_bj,z_h,eta,Q,T_t,xi,zeta,delta_k_t,k_i_t,M_ki,M_kf: (log((1/4)*(Q + sqrt(4*M**2*x_bj**2 + Q**2))*(-4*M**2*Q**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + 2*M_h**2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) - 2*Q**3*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + Q**2 + sqrt((4*M**2*Q**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 - 4*M*M_h*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + 2*M_h**2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + 2*Q**3*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) - Q**2)*(4*M**2*Q**2*x_bj**2/(Q + sqrt(4*M**2*x_bj**2 + Q**2))**2 + 4*M*M_h*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + 2*M_h**2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + 2*Q**3*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) - Q**2)))/(M_h*x_bj*(2*M**2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)) + Q**2))))
get_eta = lambda M,M_h,x_bj,z_h,Q,T_t: log(((Q*z_h*(Q**2 -(2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)))**2*M**2))/(2*(2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)))**2*sqrt(M_h**2 + T_t**2*z_h**2)))+(Q/((2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)))*M))*(sqrt((z_h**2*(Q**2-(2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)))**2*M**2)**2)/(4*(2*Q*x_bj/(Q + sqrt(4*M**2*x_bj**2 + Q**2)))**2*M**2*(M_h**2+z_h**2*T_t**2))-1)))


# Nachtmann x
def xnac(xb,Mpp,QQ):
    val = 2*xb/(1+np.sqrt(1+4*xb**2*Mpp**2/((QQ)**2*1.0))) 
    return val

# Nachtmann z
def zn(zh,xb,QQ,Mpp,Mhh,pt):
    term1 = xnac(xb,Mpp,QQ)*zh/xb/2.0
    term2 =np.sqrt(Mhh**2+pt**2)
    term3 = 1+np.sqrt(1-(4*Mpp**2)*(term2**2)*xb**2/(zh**2*QQ**4))
    return term1*term3
