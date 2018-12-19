#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
from tools.reader import _READER
from tools.config import conf


class READER(_READER):

    def __init__(self):
        pass

    def get_W2(self, tab, k):
        cols = tab.columns.values
        tab['W2'] = pd.Series(conf['aux'].M**2 + tab['Q2'] /
                              tab['x'] - tab['Q2'], index=tab.index)
        return tab

    def get_rap(self, tab, k):
        Q2 = tab['Q2']
        pT = tab['pT']
        x = tab['x']
        zh = tab['z']
        Q = Q2**0.5
        hadron = tab['hadron'][0]
        M = conf['aux'].M
        M2 = conf['aux'].M**2
        if 'h' in hadron:
            Mh = conf['aux'].Mpi
        if 'pi' in hadron:
            Mh = conf['aux'].Mpi
        if 'k' in hadron:
            Mh = conf['aux'].Mk
        MhT = np.sqrt(Mh**2 + pT**2)
        xn = 2 * x / (1 + np.sqrt(1 + 4 * x**2 * M2 / Q2))
        yp = 0.5 * np.log(Q2 / xn**2 / M2)
        yh = np.log(Q * zh * (Q2 - xn**2 * M2) / (2 * xn**2 * M2 * MhT)
                    - Q / (xn * M) * np.sqrt(zh**2 * (Q2 - xn**2 * M2)**2 / (4 * xn**2 * M2 * MhT**2) - 1))
        dy = yp - yh

        MiT = 0.5
        MfT = 0.5
        
        yi = np.log(Q/MiT)
        yf = -np.log(Q/MfT)

        R = np.sqrt( (MfT/MiT * (np.exp(yf-yh) + np.exp(yh-yf))/(np.exp(yi-yh)-np.exp(yh-yi)))**2 )

        tab['yh'] = pd.Series(yh, index=tab.index)
        tab['yp'] = pd.Series(yp, index=tab.index)
        tab['yi'] = pd.Series(yi, index=tab.index)
        tab['yf'] = pd.Series(yf, index=tab.index)
        tab['dy'] = pd.Series(yp - yh, index=tab.index)
        tab['R'] = pd.Series(R,index=tab.index)
        tab['lnR'] = pd.Series(np.log(R),index=tab.index)
        return tab

    def modify_table(self, tab, k):
        tab = self.get_W2(tab, k)
        tab = self.get_rap(tab, k)
        tab = self.apply_cuts(tab, k)
        return tab


if __name__ == "__main__":

    conf['datasets'] = {}
    conf['datasets']['sidis'] = {}

    conf['datasets']['sidis']['xlsx'] = {}
    conf['datasets']['sidis']['xlsx'][5000] = '../../database/sidis/expdata/5000.xlsx'
    conf['datasets']['sidis']['xlsx'][5008] = '../../database/sidis/expdata/5008.xlsx'

    conf['datasets']['sidis']['filters'] = {}
    conf['datasets']['sidis']['filters'][1] = {}
    conf['datasets']['sidis']['filters'][1]['list'] = range(1000, 2000)
    conf['datasets']['sidis']['filters'][1]['cond'] = []
    conf['datasets']['sidis']['filters'][1]['cond'].append("z<0.6")
    conf['datasets']['sidis']['filters'][1]['cond'].append("Q2>1.69")
    conf['datasets']['sidis']['filters'][1]['cond'].append("pT>0.2 and pT<0.9")

    TAB = READER().load_data_sets('sidis')
    print TAB
