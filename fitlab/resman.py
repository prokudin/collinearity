import os
import sys
import time
import logging
import fitlab.parallel as parallel
import numpy as np
import qcdlib
from qcdlib import pdf0,ff0,pdf1,ff1,gk0
import qcdlib.aux
import qcdlib.alphaS
import qcdlib.interpolator
import obslib.sidis.residuals
import obslib.sidis.reader
import obslib.sia.stfuncs
import obslib.sia.residuals
import obslib.sia.reader
import obslib.moments.reader
import obslib.moments.residuals
import obslib.AN_pp.AN_theory
import obslib.AN_pp.residuals
import obslib.AN_pp.reader
from parman import PARMAN
from tools.config import load_config, conf

class RESMAN:

    def __init__(self, mode='solo', ip=None, nworkers=None):

        # initial setup for parallelization
        self.mode = mode
        self.master = None
        self.slave = None
        self.broker = None

        if self.mode == 'parallel':
            self.master = parallel.Server(ip=ip)
            self.slave = parallel.Worker(ip=ip)
            self.broker = parallel.Broker()
        elif self.mode == 'master':
            self.master = parallel.Server(ip=ip)
        elif self.mode == 'slave':
            self.slave = parallel.Worker(ip=ip)

        # theory setups
        conf['aux'] = qcdlib.aux.AUX()
        self.setup_tmds()
        conf['parman'] = PARMAN()

        if 'datasets' in conf:

            if 'sidis'   in conf['datasets']: self.setup_sidis()
            if 'sia'     in conf['datasets']: self.setup_sia()
            if 'moments' in conf['datasets']: self.setup_moments()
            if 'AN'      in conf['datasets']: self.setup_AN()

        # final setups for paralleization
        if self.mode == 'parallel':
            self.broker.run_subprocess()
            self.slave.run_subprocess()

    def setup_tmds(self):

        if 'pdf'          in conf['params']: conf['pdf']          = pdf0.PDF()
        if 'gk'           in conf['params']: conf['gk']           = gk0.GK()
        if 'transversity' in conf['params']: conf['transversity'] = pdf1.PDF()
        if 'sivers'       in conf['params']: conf['sivers']       = pdf1.PDF()
        if 'boermulders'  in conf['params']: conf['boermulders']  = pdf1.PDF()

        if 'ffpi' in conf['params']: conf['ffpi'] = ff0.FF('pi')
        if 'ffk'  in conf['params']: conf['ffk']  = ff0.FF('k')
        if 'collinspi' in conf['params']: conf['collinspi'] = ff1.FF('pi')
        if 'collinsk'  in conf['params']: conf['collinsk']  = ff1.FF('k')

    def setup_sidis(self):
        conf['sidis tabs']    = obslib.sidis.reader.READER().load_data_sets('sidis')
        self.sidisres = obslib.sidis.residuals.RESIDUALS()

        if (self.slave):
            self.slave.add_mproc('sidis', self.sidisres.mproc)
        if (self.master):
            self.sidisres.mproc = self.master.wrap_mproc(
                'sidis', self.sidisres.mproc)

    def setup_sia(self):
        conf['sia tabs']    = obslib.sia.reader.READER().load_data_sets('sia')
        conf['sia stfuncs'] = obslib.sia.stfuncs.STFUNCS()
        self.siares = obslib.sia.residuals.RESIDUALS()

        if (self.slave):
            self.slave.add_mproc('sia', self.siares.mproc)
        if (self.master):
            self.siares.mproc = self.master.wrap_mproc(
                'sia', self.siares.mproc)

    def setup_moments(self):
        conf['moments tabs'] = obslib.moments.reader.READER().load_data_sets('moments')
        self.momres = obslib.moments.residuals.RESIDUALS()

        if (self.slave):
            self.slave.add_mproc('moments', self.momres.mproc)
        if (self.master):
            self.momres.mproc = self.master.wrap_mproc(
                'moments', self.momres.mproc)

    def setup_AN(self):
        conf['AN tabs']   = obslib.AN_pp.reader.READER().load_data_sets('AN')
        conf['AN theory'] = obslib.AN_pp.AN_theory.ANTHEORY()
        self.ANres = obslib.AN_pp.residuals.RESIDUALS()

        if (self.slave):
            self.slave.add_mproc('AN', self.ANres.mproc)
        if (self.master):
            self.ANres.mproc = self.master.wrap_mproc('AN', self.ANres.mproc)

    def get_residuals(self, par, calc=True, simple=False):
        """
        Get the residuals that result from the given parameters.

        Args:
            par (vector): A vector (numpy array) of the parameters for the fit.

        Returns:
            A 3-tuple of the residuals, *'r-residuals'*, and normalized
            residuals. The *r-residuals* result from the correlational
            considerations.
        """
        conf['parman'].set_new_params(par)

        if (self.master):
            self.master.assign_work()

        res, rres, nres = [], [], []
        if 'sidis' in conf['datasets']:
            out = self.sidisres.get_residuals(calc=calc, simple=simple)
            res = np.append(res, out[0])
            rres = np.append(rres, out[1])
            nres = np.append(nres, out[2])
        if 'sia' in conf['datasets']:
            out = self.siares.get_residuals(calc=calc, simple=simple)
            res = np.append(res, out[0])
            rres = np.append(rres, out[1])
            nres = np.append(nres, out[2])
        if 'moments' in conf['datasets']:
            out = self.momres.get_residuals(calc=calc, simple=simple)
            res = np.append(res, out[0])
            rres = np.append(rres, out[1])
            nres = np.append(nres, out[2])
        if 'AN' in conf['datasets']:
            out = self.ANres.get_residuals(calc=calc, simple=simple)
            res = np.append(res, out[0])
            rres = np.append(rres, out[1])
            nres = np.append(nres, out[2])
        return res, rres, nres

    def gen_report(self, verb=0, level=0):
        """
        Get a report.

        Returns:
            A list of the lines of the report.
        """
        L = []
        if 'sidis' in conf['datasets']:
            L.extend(self.sidisres.gen_report(verb, level))
        if 'sia' in conf['datasets']:
            L.extend(self.siares.gen_report(verb, level))
        if 'moments' in conf['datasets']:
            L.extend(self.momres.gen_report(verb, level))
        if 'AN' in conf['datasets']:
            L.extend(self.ANres.gen_report(verb, level))
        return L

    def run_worker(self):
        self.slave.run()

    def shutdown(self):
        """Release the resources used by the RESMAN instance."""
        if (self.master):
            self.master.finis()
            time.sleep(3)
        if (self.slave):
            self.slave.stop()
        if (self.broker):
            self.broker.stop()
