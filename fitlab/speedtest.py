import sys
import os
import numpy as np
from timeit import default_timer as timer
from tools.tools import lprint
from tools.config import conf


class SPEEDTEST:

    def __init__(self):
        pass

    def run(self):
        resman = conf['resman']
        parman = conf['parman']
        args = conf['args']

        a = timer()
        np.random.seed(12345)
        for i in range(10):
            msg = '[%d/10]' % i
            lprint(msg)
            R = np.random.randn(len(parman.par))
            par = parman.par + 1e-3 * R
            res, rres, nres = resman.get_residuals(par)
        print '\nt-elapsed (sec): ', (timer() - a) / 10
        print 'npts=%d' % res.size
        print 'chi2=%f' % np.sum(res**2)
        print 'chi2/npts=%f' % (np.sum(res**2) / res.size)
        print 'chi2(r)=%f' % np.sum(rres**2)
        print 'chi2(norm)=%f' % np.sum(nres**2)
