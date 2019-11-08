#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
import cPickle
import zlib
import argparse


def save(data, name):
    compressed = zlib.compress(cPickle.dumps(data))
    with open(name, "wb") as f:
        f.writelines(compressed)


def load(name):
    with open(name, "rb") as compressed:
        data = cPickle.loads(zlib.decompress(compressed.read()))
    return data


def retrieve(tab, entry):
    columns = tab.columns.values
    val = None
    for c in columns:
        if entry.upper() == c.upper():
            val = tab[c].values[0]
            break
    return val


def get_info_file():
    L = os.listdir('expdata')
    T = []
    for l in L:
        fname = 'expdata/%s' % l
        idx = l.replace('.xlsx', '')
        tab = pd.read_excel(fname)
        tab = tab[:1]
        dep = retrieve(tab, 'dependence')
        tar = retrieve(tab, 'target')
        had = retrieve(tab, 'hadron')
        col = retrieve(tab, 'col')
        obs = retrieve(tab, 'obs')
        if dep == None:
            dep = '-'
        T.append([idx, tar, had, obs, col, dep])
    T = np.transpose(T)
    H = ['idx', 'tar', 'had', 'obs', 'col', 'dep']
    D = {}
    for i in range(len(H)):
        D[H[i]] = T[i]
    D = pd.DataFrame(D)
    save(D, '.info')


def gen_input_file(table):
    l = 'conf["datasets"]["sidis"]["xlsx"][%s]="../database/sidis/expdata/%s.xlsx" %s'
    for i in table.index:
        col = table.col[i]
        obs = table.obs[i]
        tar = table.tar[i]
        had = table.had[i]
        dep = table.dep[i]
        idx = table.idx[i]
        info = ' # %8s %9s %4s %4s' % (col, tar, had, dep)
        print l % (idx, idx, info)


if __name__ == "__main__":

    ap = argparse.ArgumentParser()
    msg = ' 0: generate .info file'
    msg += ' 1: print datasets info'
    msg += ' 2: print inputfile'
    ap.add_argument('option', type=int, help=msg)
    args = ap.parse_args()

    if args.option == 0:
        get_info_file()
    elif args.option > 0:
        try:
            tab = load('.info')
        except:
            raise ValueError(
                '.info file not found. Use option 0 to generate .info file')
        tab = tab.query('obs=="AUTcollins"')
        compass = tab.query('col=="compass"').sort_values(
            by=['tar', 'had', 'dep'])
        hermes = tab.query('col=="HERMES"').sort_values(
            by=['tar', 'had', 'dep'])
        if args.option == 1:
            print compass
            print hermes
        elif args.option == 2:
            gen_input_file(compass)
            gen_input_file(hermes)
