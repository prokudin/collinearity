#!/usr/bin/env python
import sys,os
import subprocess
import numpy as np
import pylab as py
import matplotlib.cm as cm
import sympy as sp
from sympy.printing.mathml import print_mathml,mathml
from sympy.parsing.sympy_parser import parse_expr
from sympy.tensor.array import MutableDenseNDimArray as Array
import copy
from sympy.utilities.lambdify import lambdastr

class Vec: 

    def __init__(self,symbol,rep='LC',mod=[1,1,1,1],theta=None):
        if rep=='LC': self.vec=self.LC(symbol,mod=mod,theta=theta)
        if rep=='MC': self.vec=self.MC(symbol,mod=mod,theta=theta)
        self.rep=rep

    def __add__(self,other):
        dum=Vec('dum',self.rep)
        dum.vec=[self.vec[i]+other.vec[i] for i in range(len(other.vec))]
        return dum

    def __sub__(self,other):
        #print self.vec
        #print other.vec
        dum=Vec('dum',self.rep)
        dum.vec=[self.vec[i]-other.vec[i] for i in range(len(other.vec))]
        return dum

    def __eq__(self,other):
        dum=Vec('dum',self.rep)
        dum.vec=[self.vec[i] for i in range(len(other.vec))]
        return dum
 
    def __mul__(self,other):
        if isinstance(other, Vec): 
            return self.dot(other.vec,self.vec)
        else:
            dum=Vec('dum',self.rep)
            dum.vec=[other*self.vec[i] for i in range(len(self.vec))]
            return dum
        
    def __rmul__(self,other):
        if isinstance(other, Vec): 
            return self.dot(other.vec,self.vec)
        else:
            dum=Vec('dum',self.rep)
            dum.vec=[other*self.vec[i] for i in range(len(self.vec))]
            return dum
        
    def LC(self,symbol,mod=[1,1,1,1],theta=None):
        p=sp.S('%s_p'%symbol)
        m=sp.S('%s_m'%symbol)
        t=sp.symbols('%s_t'%symbol,positive = True, real = True)
        if theta==None: theta=sp.S('theta_%s'%symbol)
        vec=[p,m,t*sp.cos(theta),t*sp.sin(theta)]
        self.p=p
        self.m=m
        self.t=t
        return [vec[i]*mod[i] for i in range(4)]
    
    def MC(self,symbol,mod=[1,1,1,1],theta=None):
        E=sp.symbols('%s_E'%symbol,positive = True, real = True)
        z=sp.symbols('%s_z'%symbol,real = True)
        t=sp.symbols('%s_t'%symbol,positive = True, real = True)
        if theta==None: theta=sp.S('theta_%s'%symbol)
        vec=[E,t*sp.cos(theta),t*sp.sin(theta),z]
        self.E=E
        self.z=z
        self.t=t
        self.theta=theta
        return [vec[i]*mod[i] for i in range(4)]

    def dot(self,A,B):
        if self.rep=='LC':
            return A[0]*B[1]+A[1]*B[0]-A[2]*B[2]-A[3]*B[3]
        elif self.rep=='MC':
            return A[0]*B[0]-A[1]*B[1]-A[2]*B[2]-A[3]*B[3]
        
    def set_plus(self,plus):
        self.vec[0]=plus
        
    def set_minus(self,minus):
        self.vec[1]=minus
        
    def set_perp(self,perp):
        self.vec[2:]=perp
        
    def set_E(self,E):
        self.vec[0]=E
        
    def set_x(self,x):
        self.vec[1]=x
        
    def set_y(self,y):
        self.vec[2]=y
        
    def set_z(self,z):
        self.vec[3]=z

    def get_E(self):
        return self.vec[0]
        
    def get_x(self):
        return self.vec[1]
        
    def get_y(self):
        return self.vec[2]
        
    def get_z(self):
        return self.vec[3]
    
    def change_plus(self,M):
        if  self.rep=='LC':
            perp2=(self.vec[2]**2+self.vec[3]**2)#.simplify()
            self.vec[0]=(perp2+M**2)/2/self.vec[1]
        
    def change_minus(self,M):
        if  self.rep=='LC':
            perp2=(self.vec[2]**2+self.vec[3]**2)#.simplify()
            self.vec[1]=(perp2+M**2)/2/self.vec[0]

    def LC2MC(self):
        if  self.rep=='LC':
            vec0=((self.vec[0]+self.vec[1])/sp.sqrt(2))#.simplify()
            vec3=((self.vec[0]-self.vec[1])/sp.sqrt(2))#.simplify()
            self.vec=[vec0,self.vec[2],self.vec[3],vec3]
            self.rep='MC'

    def MC2LC(self):
        if  self.rep=='MC':
            vec_p=((self.vec[0]+self.vec[3])/sp.sqrt(2))#.simplify()
            vec_m=((self.vec[0]-self.vec[3])/sp.sqrt(2))#.simplify()
            self.vec=[vec_p,vec_m,self.vec[1],self.vec[2]]
            self.rep='LC'

    def norm(self):
        return sp.sqrt(sum([self.vec[i]**2 for i in range(1,4)]))#.simplify()
        
    def uvec(self):
        dum=Vec('dum','MC')
        dum.set_E(sp.S(0))
        norm=self.norm()
        dum.vec[1:]=[self.vec[i]/norm for i in range(1,len(self.vec))]
        return dum
    
    def get_mass(self):
        if self.rep=='MC':
            return sp.sqrt(self.vec[0]**2-self.norm()**2)
        elif self.rep=='LC':
            return sp.sqrt(2*self.vec[0]*self.vec[1]-self.vec[2]**2-self.vec[3]**2)

    def get_sqmass(self):
        if self.rep=='MC':
            return self.vec[0]**2-self.norm()**2
        elif self.rep=='LC':
            return 2*self.vec[0]*self.vec[1]-self.vec[2]**2-self.vec[3]**2

    @staticmethod
    def dot3(A,B):
        return sum([A.vec[i]*B.vec[i] for i in range(1,4)])
    
    @staticmethod
    def cross(A,B):
        dum=Vec('dum','MC')
        dum.set_E(sp.S(0))
        x= A.vec[2]*B.vec[3]-A.vec[3]*B.vec[2]
        y=-A.vec[1]*B.vec[3]+A.vec[3]*B.vec[1]
        z= A.vec[1]*B.vec[2]-A.vec[2]*B.vec[1]
        dum.vec[1:]=[x,y,z]
        return dum
    
    def simplify_vec(self):
        for i in range(4):
            self.vec[i]=self.vec[i].simplify()

    def subs_vec(self,A,B):
        for i in range(4):
            self.vec[i]=self.vec[i].subs(A,B)
    
    def rotate(self,K,cos,sin):
        #--rootate around K 
        dum=Vec('dum','MC')
        UK=K.uvec()
        dum.vec[1:]=(self*cos+self.cross(UK,self)*sin+UK*(self.dot3(UK,self)*(1-cos))).vec[1:]
        dum.vec[0]=self.vec[0]
        return dum
        
    def get_cos(self,K):
        #--get cos along K
        UK=K.uvec()
        UV=self.uvec()
        return self.dot3(UV,UK)
    
    def get_sin(self,K):
        #--get sin along K
        UK=K.uvec()
        UV=self.uvec()
        return self.cross(UV,UK).norm()
        
    @staticmethod
    def get_boost_CS(A,bA):
        #--get cosh and sinsh for A -> boosted A (bA)
        #--A and bA are 2D vectors 0: energy, 1: momentum along boost dir
        cosh=(bA[1]*A[1]-bA[0]*A[0])/(A[1]**2-A[0]**2)
        if A[1]==0: sinh=(bA[1]-cosh*A[1])/A[0]
        else:       sinh=(bA[0]-cosh*A[0])/A[1]
        return cosh,sinh
        
    def boost(self,V,C=None,S=None):
        
        if C==None and S==None:
            #--boost to rest frame of V
            E=V.vec[0]
            P=V.norm()
            M=V.get_mass()
            C,S=self.get_boost_CS([E,P],[M,0])
            
        UV=V.uvec()
        par=self.dot3(self,UV)*UV
        per=self-par
        E=self.vec[0]
        npar=self.dot3(self,UV)#par.norm()
        #npar=par.norm()
        bE  =(C*E + S*npar)#.simplify()
        bpar=(S*E + C*npar)#.simplify()
        new=per+bpar*UV
        new.vec[0]=bE
        return new

def gen_params():
    params={}

    params['M']   = 0.938
    params['M_h'] = 0.139

    params['x_bj']= 0.1
    params['z_h'] = 0.1
    params['Q']   = 10.0
    params['T_t'] = 0.1
    params['eta'] = 0

    params['xi']  = 0.3
    params['zeta']= 0.3
    params['delta_k_t']=0.1
    params['k_i_t']=0.01
    params['M_ki']=0.1
    params['M_kf']=0.1
    return params

def evaluate(func,params=None,verb=False):
    
    if params==None: params=gen_params()
        
    msg=func.__doc__.replace(' ','').split('func(')[1].split('Expression')[0]
    msg=msg.replace('\n','').rstrip(')').split(',')
    if verb: print(msg)
    for _ in msg: 
        if _ not in params: 
            print('%s is missing in params'%_)
            return None
    
    args=[params[_] for _ in msg]
    return func(*args)

def evaluate2(func,params=None,verb=False):
    if params==None: params=gen_params()
    keys=['M','M_h','x_bj','z_h','eta','Q','T_t','xi','zeta','delta_k_t','k_i_t','M_ki','M_kf']
    args=[params[_] for _ in keys]
    return func(*args)

def get_massless(exp,replace=[],n=3):
    if len(replace)>0:
        for _ in replace:
            A,B=_
            exp=exp.subs(A,B)
    exp=exp.subs(Q,1/eps)        
    return sp.series(exp,eps,0,n=n).subs(eps,1/Q).removeO()

def gen_expressions():

    xb=sp.symbols('x_bj',positive = True, real = True)
    xN=sp.symbols('x_N',positive = True, real = True)
    zh=sp.symbols('z_h',positive = True, real = True)
    zN=sp.symbols('z_N',positive = True, real = True)
    xh=sp.symbols('x_h',positive = True, real = True)
    
    Q=sp.symbols('Q',positive = True, real = True)
    M=sp.symbols('M',positive = True, real = True)
    Mh=sp.symbols('M_h',positive = True, real = True)
    Mki=sp.symbols('M_ki',positive = True, real = True)
    Mkf=sp.symbols('M_kf',positive = True, real = True)
    Mx=sp.symbols('M_x',positive = True, real = True)
    
    MhT=sp.symbols('M_hT',positive = True, real = True)
    qt=sp.S('q_t')
    eps=sp.S('epsilon')
    eta=sp.S('eta')
    xi=sp.symbols('xi',positive = True, real = True)
    zeta=sp.symbols('zeta',positive = True, real = True)
    hxN=sp.symbols('\hat{x}_N',positive = True, real = True)
    hzN=sp.symbols('\hat{z}_N',positive = True, real = True)
    
    
    ki=Vec('k_i','LC',theta=0)
    kf=Vec('k_f','LC',theta=0)
    kx=Vec('k_x','LC',theta=0)
    
    q=Vec('q','LC',mod=[1,1,0,0],theta=0)
    P=Vec('P','LC',mod=[1,1,0,0],theta=0)
    H=Vec('H','LC',theta=0)
    T=Vec('T',theta=0,mod=[0,0,1,1])
    dk=Vec('delta_k',theta=0,mod=[0,0,1,1])
    
    q.set_plus(-Q/sp.sqrt(2))
    q.set_minus(Q/sp.sqrt(2))
    P.set_plus(Q/xN/sp.sqrt(2))
    P.set_minus(xN*M**2/sp.sqrt(2)/Q)
    H.set_minus(zN*q.vec[1])
    H.set_perp((-zN*T).vec[2:])
    H.change_plus(Mh)
    
    ki.set_plus((xN/hxN)*P.vec[0])
    kf.set_minus(H.vec[1]/(zN/hzN))
    kf.set_perp((-hzN*T+dk).vec[2:])
    
    ki.change_minus(Mki)
    kf.change_plus(Mkf)

    
    def build(exp):
        exp=exp.subs(hzN,zN/zeta)
        exp=exp.subs(hxN,xN/xi)
        exp=exp.subs(zN,zN2zh)
        exp=exp.subs(xN,xN2xb)
        exp=exp.subs(Mki**2,-Mki**2)
        args=(M,Mh,xb,zh,eta,Q,T.t,xi,zeta,dk.t,ki.t,Mki,Mkf)
        return sp.lambdify(args, exp)

    def buildstr(exp):
        exp=exp.subs(hzN,zN/zeta)
        exp=exp.subs(hxN,xN/xi)
        exp=exp.subs(zN,zN2zh)
        exp=exp.subs(xN,xN2xb)
        exp=exp.subs(Mki**2,-Mki**2)
        args=(M,Mh,xb,zh,eta,Q,T.t,xi,zeta,dk.t,ki.t,Mki,Mkf)
        out=lambdastr(args, exp)
        out=out.replace('math.','')
        return out


    data={} #--store values for benchmark
    
    code=['from __future__ import division'] 
    code.append('from numpy import sqrt,log,exp')    

    #--build xN and zN 
    
    zN2zh=sp.solve(zh-2*xb*(P*H)/Q**2,zN)[1]
    zN2zh=zN2zh.simplify()
    
    xN2xb=sp.solve(xb-Q**2/(2*P*q),xN)[0]
    xN2xb=xN2xb*(Q+sp.sqrt(Q**2+4*M**2*xb**2))
    xN2xb=xN2xb.expand()
    xN2xb/=(Q+sp.sqrt(Q**2+4*M**2*xb**2))
    
    get_xN=build(xN2xb)
    get_zN=build(zN2zh)
    data['zN']=evaluate(get_zN)
    data['xN']=evaluate(get_xN)
    
    code.append('get_xN=%s'%buildstr(xN2xb))
    code.append('get_zN=%s'%buildstr(zN2zh))
    
    #--build kx, k  
    kx=ki+q-kf
    kx.simplify_vec()
    #kx.vec=[_.subs(Mkf,0).subs(Mki,0).subs(ki.t,0) for _ in kx.vec]
    Mx2_val=kx.get_sqmass()
    
    k=kf-q
    
    #--build R1
    R1=(H*kf)/(H*ki)
    R1=R1.simplify()
    get_R1=build(R1)
    data['R1']=evaluate(get_R1)
    code.append('get_R1=%s'%buildstr(R1))
    
    
    #--build R2
    R2=(k.get_sqmass()/Q**2)
    R2=R2.simplify()
    get_R2=build(R2)
    data['R2']=evaluate(get_R2)
    code.append('get_R2=%s'%buildstr(R2))
    
    #--build R3
    R3=(kx.get_sqmass()/Q**2)
    R3=R3.simplify()
    get_R3=build(R3)
    data['R3']=evaluate(get_R3)
    code.append('get_R3=%s'%buildstr(R3))

    #--build R4
    R4=H*P/Q**2
    R4=R4.simplify()
    get_R4=build(R4)
    data['R4']=evaluate(get_R4)
    code.append('get_R4=%s'%buildstr(R4))

    #--build proton rapidity
    yp=1/sp.S(2)*sp.log(P.vec[0]/P.vec[1]).simplify()
    get_yp=build(yp)
    data['yp']=evaluate(get_yp)
    code.append('get_yp=%s'%buildstr(yp))
    
    #--build hadron rapidity
    yh=1/sp.S(2)*sp.log(H.vec[0]/H.vec[1]).simplify()
    get_yh=build(yh)
    data['yh']=evaluate(get_yh)
    code.append('get_yh=%s'%buildstr(yh))
    
    #--build incomming parton rapidity
    yi=(1/sp.S(2)*sp.log(ki.vec[0]/ki.vec[1])).simplify()
    get_yi=build(yi)
    data['yi']=evaluate(get_yi)
    code.append('get_yi=%s'%buildstr(yi))
    
    #--build outgoing parton rapidity
    yf=(1/sp.S(2)*sp.log(kf.vec[0]/kf.vec[1])).simplify()
    get_yf=build(yf)
    data['yf']=evaluate(get_yf)
    code.append('get_yf=%s'%buildstr(yf))

    #--build 
    W=P+q-H
    W2=W.get_sqmass().simplify()
    get_W2=build(W2)
    data['W2']=evaluate(get_W2)
    code.append('get_W2=%s'%buildstr(W2))

    #--get Htmax

    W=P+q
    W.LC2MC()
    W0=W.vec[0].simplify()
    W3=W.vec[3].simplify()

    zero=W0-sp.sqrt(Mh**2+H.t**2)*sp.cosh(eta)\
           -sp.sqrt(M**2+H.t**2+(W3-sp.sqrt(Mh**2+H.t**2)*sp.sinh(eta))**2)

    sol=sp.solve(zero,H.t)
    Htmax=sol[1]

    get_Htmax=build(Htmax)
    data['Htmax']=evaluate(get_Htmax)
    code.append('get_Htmax=%s'%buildstr(Htmax))


    #--get eta min
    zero=W0-Mh*sp.cosh(eta)-sp.sqrt(M**2+(W3-Mh*sp.sinh(eta))**2)

    sol=sp.solve(zero,eta)
    etamin=sol[0]
    etamax=sol[1]

    get_etamin=build(etamin)
    data['etamin']=evaluate(get_etamin)
    code.append('get_etamin=%s'%buildstr(etamin))

    get_etamax=build(etamax)
    data['etamax']=evaluate(get_etamax)
    code.append('get_etamax=%s'%buildstr(etamax))

    #--write expressions to file
    code=[l+'\n' for l in code]
    F=open('ratlib.py','w')
    F.writelines(code)
    F.close()

    return data

def test(data):
    import ratlib as rat
    
    for _ in data:
      if   _=='xN': diff=data[_]-evaluate2(rat.get_xN)
      elif _=='zN': diff=data[_]-evaluate2(rat.get_zN)
      elif _=='R1': diff=data[_]-evaluate2(rat.get_R1)
      elif _=='R2': diff=data[_]-evaluate2(rat.get_R2)
      elif _=='R3': diff=data[_]-evaluate2(rat.get_R3)
      elif _=='R4': diff=data[_]-evaluate2(rat.get_R4)
      elif _=='yp': diff=data[_]-evaluate2(rat.get_yp)
      elif _=='yh': diff=data[_]-evaluate2(rat.get_yh)
      elif _=='yi': diff=data[_]-evaluate2(rat.get_yi)
      elif _=='yf': diff=data[_]-evaluate2(rat.get_yf)
      elif _=='W2': diff=data[_]-evaluate2(rat.get_W2)

      print(_,data[_],diff)

if __name__=="__main__":

    data=gen_expressions()
    test(data)

    

