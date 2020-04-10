#!/usr/bin/env python
import sys
import os
import numpy as np
import pylab as py
import pandas as pd

F = os.listdir('./')
F = [f for f in F if f.endswith('.dat')]

idx = 3000

for f in F:
    columns = f.replace('.dat', '').split('_')
    if 'compass' in f:
        dum1, dum2, hadron, proj, col, target, year = columns
    elif 'lba' in f:
        dum1, dum2, hadron, proj, col, year = columns
        target = 'proton'
    elif 'jlab' in f:
        dum1, dum2, col, target, hadron, proj, year = columns

    if hadron == 'pip':
        hadron = 'pi+'
    elif hadron == 'pim':
        hadron = 'pi-'
    elif hadron == 'pi0':
        hadron = 'pi0'
    elif hadron == 'Kp':
        hadron = 'k+'
    elif hadron == 'Km':
        hadron = 'k-'
    elif hadron == 'K0':
        hadron = 'k0'
    elif hadron == 'hp':
        hadron = 'h+'
    elif hadron == 'hm':
        hadron = 'h-'
    elif hadron == 'h0':
        hadron = 'h0'
    else:
        print 'ERR'
        sys.exit()

    with open(f) as fh:
        L = [l.strip() for l in fh.readlines()]

    for l in L:
        if '<' and '>' not in l:
            continue
        if 'stat' in l:
            H = l
        if 'sigma' in l:
            H = l

    H = H.replace('#', '')
    H = H.replace('Q^2', 'Q2')
    H = H.replace('Phperp', 'pT')
    H = H.replace('<', '')
    H = H.replace('>', '')
    H = H.replace('.', '')
    H = H.replace('and', '')
    H = H.replace('sigma', 'error')
    H = H.replace('Collins', 'value')
    H = H.split()

    L = [l for l in L if l.startswith('#') == False]
    L = [[float(x) for x in l.split()] for l in L if l != '']

    H.extend(['target', 'hadron', 'Experiment', 'dependence'])
    for i in range(len(L)):
        L[i].extend([target, hadron, col, proj])

    data = pd.DataFrame(L, columns=H)

    #writer = pd.ExcelWriter('../../expdata/%d.xlsx'%idx)
    #data.to_excel(writer, index=False,header=True)
    # writer.save()

    print '|%d|[link][?]|SIDIS|%s|%s|AUTcollins|%s|%s|' % (
        idx, target, hadron, col, proj)

    idx += 1
