#!/usr/bin/env python
import sys
import os
import numpy as np
from tools.tools import load_config
from qcdlib.auxiliary import AUX  # renamed from "aux"
from tools.config import conf

eu2, ed2 = 4 / 9., 1 / 9.
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

def _get_FUT(x, z, Q2, pT, tar, had, F, D, w_tar, w_had):

    if 'pi' in had: Mh = conf['aux'].Mpi
    elif 'k' in had: Mh = conf['aux'].Mk

    if had.endswith('+'):

        wq = z**2 * np.abs(w_tar) + np.abs(w_had)
        K = 2 * x * z * pT * Mh / wq
        gauss = np.exp(-pT**2 / wq) / (np.pi * wq)
        return np.sum(e2 * K * F * D * gauss)

    elif had.endswith('-'):
        D = conf['aux'].charge_conj(D)
        w_had = conf['aux'].charge_conj(w_had)
        wq = z**2 * np.abs(w_tar) + np.abs(w_had)
        K = 2 * x * z * pT * Mh / wq
        gauss = np.exp(-pT**2 / wq) / (np.pi * wq)
        return np.sum(e2 * K * F * D * gauss)

    elif had.endswith('0'):

        Dp = D
        Dm = conf['aux'].charge_conj(D)
        w_hadp = w_had
        w_hadm = conf['aux'].charge_conj(w_had)

        wqp = z**2 * np.abs(w_tar) + np.abs(w_hadp)
        wqm = z**2 * np.abs(w_tar) + np.abs(w_hadm)
        gaussp = np.exp(-pT**2 / wqp) / (np.pi * wqp)
        gaussm = np.exp(-pT**2 / wqm) / (np.pi * wqm)
        Kp = 2 * x * z * pT * Mh / wqp
        Km = 2 * x * z * pT * Mh / wqm

        FUTp = np.sum(e2 * Kp * F * Dp * gaussp)
        FUTm = np.sum(e2 * Km * F * Dm * gaussm)
        return 0.5 * (FUTp + FUTm)

def get_FUT(x, z, Q2, pT, tar, had):

    # get collinear parts (proton and positive hadrons)
    F = conf['transversity'].get_C(x, Q2)
    if 'pi' in had: D = conf['collinspi'].get_C(z, Q2)
    elif 'k' in had: D = conf['collinsk'].get_C(z, Q2)
    else:
        raise ValueError('had=%s not supported in sidis collins AUT' % had)
    F[0], D[0] = 0, 0  # set glue to zero

    # get widths (proton and positive hadrons)
    w_tar = conf['transversity'].get_widths(Q2)
    if 'pi' in had: w_had = np.abs(conf['collinspi'].get_widths(Q2))
    elif 'k' in had: w_had = np.abs(conf['collinsk'].get_widths(Q2))

    # build structure function
    K = x
    if tar == 'p':

        return _get_FUT(x, z, Q2, pT, tar, had, F, D, w_tar, w_had)

    elif tar == 'n':

        F = conf['aux'].p2n(F)
        w_tar = conf['aux'].p2n(w_tar)
        return _get_FUT(x, z, Q2, pT, tar, had, F, D, w_tar, w_had)

    elif tar == 'd':

      return 0.5 * (get_FUT(x, z, Q2, pT, 'p', had) + get_FUT(x, z, Q2, pT, 'n', had))


if __name__ == '__main__':

    from .qcdlib.pdf1 import PDF
    from .qcdlib.ff1 import FF
    conf['aux'] = AUX()
    conf['transversity'] = PDF()
    conf['collinspi'] = FF('pi')
    conf['collinsk'] = FF('k')

    x = 0.25
    z = 0.5
    Q2 = 2.4
    mu2 = Q2
    pT = 0.3
    tar = 'p'
    had = 'pi+'

    print(get_FUT(x, z, Q2, pT, 'p', 'pi+'))
    print(get_FUT(x, z, Q2, pT, 'p', 'pi-'))
    print(get_FUT(x, z, Q2, pT, 'p', 'pi0'))

    print(get_FUT(x, z, Q2, pT, 'n', 'pi+'))
    print(get_FUT(x, z, Q2, pT, 'n', 'pi-'))
    print(get_FUT(x, z, Q2, pT, 'n', 'pi0'))

    print(get_FUT(x, z, Q2, pT, 'd', 'pi+'))
    print(get_FUT(x, z, Q2, pT, 'd', 'pi-'))
    print(get_FUT(x, z, Q2, pT, 'd', 'pi0'))



