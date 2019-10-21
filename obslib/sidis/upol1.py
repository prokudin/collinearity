#!/usr/bin/env python
import sys
import os
import numpy as np
from .tools.tools import load_config
from .qcdlib.auxiliary import AUX  # renamed from "aux"
from .tools.config import conf

eu2, ed2 = 4/9., 1/9.
e2 = []
e2.append(0)    # g
e2.append(eu2)  # u
e2.append(eu2)  # ub
e2.append(ed2)  # d
e2.append(ed2)  # db
e2.append(ed2)  # s
e2.append(ed2)  # sb
e2.append(0)    # c
e2.append(0)    # cb
e2.append(0)    # b
e2.append(0)    # bb
e2 = np.array(e2)

def _get_FUU(x,z,Q2,pT,tar,had,F,D,w_tar,w_had):
    """
    We will use a simple form of TMD evolution based on large bT dominance where Ktilde is constant
    """

    K = x

    # Standard choice:
    evolution=conf['gk'].get_gk(Q2)*np.exp(-conf['gk'].get_pertsud(Q2))

    # Ted's idea on NP behavior:
    # evolution=conf['gk'].get_gk(Q2,x,z)*np.exp(-conf['gk'].get_pertsud(Q2))

    if had.endswith('+'):

        wq = z**2 * np.abs(w_tar) + np.abs(w_had)
        gauss = np.exp(-pT**2 / wq) / (np.pi * wq)
        return np.sum(e2*K*F*D*evolution*gauss)

    elif had.endswith('-'):

        D=conf['aux'].charge_conj(D)
        w_had=conf['aux'].charge_conj(w_had)
        wq = z**2 * np.abs(w_tar) + np.abs(w_had)
        gauss = np.exp(-pT**2 / wq) / (np.pi * wq)
        return np.sum(e2*K*F*D*evolution*gauss)

    elif had.endswith('0'):

        Dp=D
        Dm=conf['aux'].charge_conj(D)
        w_hadp=w_had
        w_hadm=conf['aux'].charge_conj(w_had)

        wqp = z**2 * np.abs(w_tar) + np.abs(w_hadp)
        wqm = z**2 * np.abs(w_tar) + np.abs(w_hadm)
        gaussp = np.exp(-pT**2 / wqp) / (np.pi * wqp)
        gaussm = np.exp(-pT**2 / wqm) / (np.pi * wqm)

        FUUp=np.sum(e2*K*F*Dp*evolution*gaussp)
        FUUm=np.sum(e2*K*F*Dm*evolution*gaussm)
        return 0.5*(FUUp+FUUm)

def get_FUU(x,z,Q2,pT,tar,had):
    """
    We will use a simple form of TMD evolution based on large bT dominance where Ktilde is constant
    """
    Q0 = conf['gk'].Q0

    # get collinear parts (proton and positive hadrons)
    F = conf['pdf'].get_C(x, Q0**2)
    if   'pi' in had:  D = conf['ffpi'].get_C(z, Q0**2)
    elif  'k' in had:  D = conf['ffk'].get_C(z, Q0**2)
    elif  'h' in had:  D = conf['ffpi'].get_C(z, Q0**2)
    F[0],D[0]=0,0  # set glue to zero

    # get widths (proton and positive hadrons)
    w_tar=conf['pdf'].get_widths(Q2)
    if   'pi' in had: w_had=np.abs(conf['ffpi'].get_widths(Q2))
    elif 'k'  in had: w_had=np.abs(conf['ffk'].get_widths(Q2))
    elif 'h'  in had: w_had=np.abs(conf['ffpi'].get_widths(Q2))

    # build structure function

    if tar=='p':

        return _get_FUU(x,z,Q2,pT,tar,had,F,D,w_tar,w_had)

    elif tar=='n':

        F=conf['aux'].p2n(F)
        w_tar=conf['aux'].p2n(w_tar)
        return  _get_FUU(x,z,Q2,pT,tar,had,F,D,w_tar,w_had)

    elif tar=='d':

      return 0.5*(get_FUU(x,z,Q2,pT,'p',had)+get_FUU(x,z,Q2,pT,'n',had))


if __name__ == '__main__':

    from .qcdlib.pdf0 import PDF
    from .qcdlib.ff0 import FF
    from .qcdlib.gk0 import GK
    conf['aux']= AUX()
    conf['pdf']=PDF()
    conf['ffpi']=FF('pi')
    conf['ffk']=FF('k')
    conf['gk']=GK()

    x = 0.25
    z = 0.5
    Q2 = 2.4
    mu2 = Q2
    pT = 0.3
    tar = 'p'
    had = 'pi+'

    print(get_FUU(x,z,Q2,pT,'p','pi+'))
    print(get_FUU(x,z,Q2,pT,'p','pi-'))
    print(get_FUU(x,z,Q2,pT,'p','pi0'))

    print(get_FUU(x,z,Q2,pT,'n','pi+'))
    print(get_FUU(x,z,Q2,pT,'n','pi-'))
    print(get_FUU(x,z,Q2,pT,'n','pi0'))

    print(get_FUU(x,z,Q2,pT,'d','pi+'))
    print(get_FUU(x,z,Q2,pT,'d','pi-'))
    print(get_FUU(x,z,Q2,pT,'d','pi0'))





