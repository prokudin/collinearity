import os
import numpy as np
from tools.tools import load, save, checkdir
from tools.config import conf
from tools.randomstr import id_generator
import nest
from multiprocessing import Process, Queue, Pool, Pipe


class MCSAMP:

    def __init__(self):
        pass

    def get_residuals(self, par):
        res, rres, nres = conf['resman'].get_residuals(par)
        if len(rres) != 0:
            res = np.append(res, rres)
        if len(nres) != 0:
            res = np.append(res, nres)
        return res

    def loglike(self, par):
        return -0.5 * np.sum(self.get_residuals(par)**2)

    def nll(self, par):
        res = self.get_residuals(par)
        return 0.5 * (np.sum(res**2))

    def linmap(self, q, pmin, pmax):
        return (pmax - pmin) * q + pmin

    def get_par_lims(self):
        plims = []
        for i in range(len(conf['parman'].order)):
            ii, k, kk = conf['parman'].order[i]
            if ii == 1:
                pmin = conf['params'][k][kk]['min']
                pmax = conf['params'][k][kk]['max']
            elif ii == 2:
                pmin = conf['datasets'][k]['norm'][kk]['min']
                pmax = conf['datasets'][k]['norm'][kk]['max']
            plims.append([pmin, pmax])
        return plims

    def single_run(self, path, nestout, factor=2, kappa=1, tol=1e-10, itmax=10, sample_size=1000, method='cov', nll_shift=0):

        # set nest params if not specified in conf
        if 'factor' not in conf:
            conf['factor'] = factor
        if 'tol' not in conf:
            conf['tol'] = tol
        if 'itmax' not in conf:
            conf['itmax'] = itmax
        if 'kappa' not in conf:
            conf['kappa'] = kappa
        if 'method' not in conf:
            conf['method'] = method
        if 'sample size' not in conf:
            conf['sample size'] = sample_size
        if 'nll shift' not in conf:
            conf['nll shift'] = nll_shift

        conf['nestout'] = nestout

        npar = len(conf['parman'].par)
        conf['nll'] = self.nll
        conf['par lims'] = self.get_par_lims()
        conf['num points'] = int(npar * conf['factor'])

        par = conf['parman'].par
        res = self.get_residuals(par)
        conf['nll shift'] = 0  # 0.5*(len(res)+conf['nll shift'])

        if 'cmd' in conf:
            cmd = conf['cmd']
        else:
            cmd = None
        nest.NEST().run(path, cmd)

        # results=nest.NEST().run()
        # save(results,path)

    def run(self):
        nruns = conf['nruns']
        for i in range(nruns):
            fname = id_generator(size=6)
            self.single_run(fname, None)

