import sys
import os
import numpy as np
from numpy.random import choice, randn, uniform
from tools.config import conf, load_config
import pandas as pd

class PARMAN:

    def __init__(self):
        self.get_ordered_free_params()
        self.set_new_params(self.par, initial=True)

    def get_ordered_free_params(self):
        self.par = []
        self.order = []

        for k in conf['params']:
            for kk in conf['params'][k]:
                if conf['params'][k][kk]['fixed'] == False:
                    # p=uniform(conf['params'][k][kk]['min'],conf['params'][k][kk]['max'],1)[0]
                    self.par.append(conf['params'][k][kk]['value'])
                    self.order.append([1, k, kk])

        if 'datasets' in conf:
            for k in conf['datasets']:
                for kk in conf['datasets'][k]['norm']:
                    if conf['datasets'][k]['norm'][kk]['fixed'] == False:
                        self.par.append(conf['datasets'][k]
                                        ['norm'][kk]['value'])
                        self.order.append([2, k, kk])

    def set_new_params(self, parnew, initial=False):
        self.shifts = 0
        semaphore = {}

        for i in range(len(self.order)):
            ii, k, kk = self.order[i]
            if ii == 1:
                if k not in semaphore:
                    semaphore[k] = 0
                if conf['params'][k][kk]['value'] != parnew[i]:
                    conf['params'][k][kk]['value'] = parnew[i]
                    semaphore[k] = 1
                    self.shifts += 1
            elif ii == 2:
                if conf['datasets'][k]['norm'][kk]['value'] != parnew[i]:
                    conf['datasets'][k]['norm'][kk]['value'] = parnew[i]
                    self.shifts += 1

        if initial:
            for k in conf['params']:
                semaphore[k] = 1

        self.propagate_params(semaphore)

    def gen_report(self):
        """
        Get a report.

        Returns:
            A list of the lines of the report.
        """
        L = []

        for k in conf['params']:
            for kk in sorted(conf['params'][k]):
                if conf['params'][k][kk]['fixed'] == False:
                    if conf['params'][k][kk]['value'] < 0:
                        L.append('%-10s  %-20s  %10.5e' %
                                 (k, kk, conf['params'][k][kk]['value']))
                    else:
                        L.append('%-10s  %-20s   %10.5e' %
                                 (k, kk, conf['params'][k][kk]['value']))

        for k in conf['datasets']:
            for kk in conf['datasets'][k]['norm']:
                if conf['datasets'][k]['norm'][kk]['fixed'] == False:
                    L.append('%10s %10s %10d  %10.5e' % (
                        'norm', k, kk, conf['datasets'][k]['norm'][kk]['value']))
        return L

    def propagate_params(self, semaphore):
        if 'pdf'          in semaphore and semaphore['pdf']          == 1: self.set_pdf_params()
        if 'gk'           in semaphore and semaphore['gk']           == 1: self.set_gk_params()
        if 'transversity' in semaphore and semaphore['transversity'] == 1: self.set_transversity_params()
        if 'sivers'       in semaphore and semaphore['sivers']       == 1: self.set_sivers_params()
        if 'boermulders'  in semaphore and semaphore['boermulders']  == 1: self.set_boermulders_params()
        if 'ffpi'         in semaphore and semaphore['ffpi']         == 1: self.set_ffpi_params()
        if 'ffk'          in semaphore and semaphore['ffk']          == 1: self.set_ffk_params()
        if 'collinspi'    in semaphore and semaphore['collinspi']    == 1: self.set_collinspi_params()
        if 'collinsk'     in semaphore and semaphore['collinsk']     == 1: self.set_collinsk_params()

    def set_constraits(self, parkind):

        for k in conf['params'][parkind]:
            if conf['params'][parkind][k]['fixed'] == True:  continue
            if conf['params'][parkind][k]['fixed'] == False: continue
            ref_par = conf['params'][parkind][k]['fixed']
            conf['params'][parkind][k]['value'] = conf['params'][parkind][ref_par]['value']

    def set_pdf_params(self):
        self.set_constraits('pdf')
        conf['pdf']._widths1_uv  = conf['params']['pdf']['widths1_uv' ]['value']
        conf['pdf']._widths1_dv  = conf['params']['pdf']['widths1_dv' ]['value']
        conf['pdf']._widths1_sea = conf['params']['pdf']['widths1_sea']['value']
        conf['pdf']._widths2_uv  = conf['params']['pdf']['widths2_uv' ]['value']
        conf['pdf']._widths2_dv  = conf['params']['pdf']['widths2_dv' ]['value']
        conf['pdf']._widths2_sea = conf['params']['pdf']['widths2_sea']['value']
        conf['pdf'].setup()

    def set_gk_params(self):
        self.set_constraits('gk')
        conf['gk'].gk  = conf['params']['gk']['gk0']['value']
        conf['gk'].Q0   = conf['params']['gk']['Q0']['value']
        conf['gk'].setup()

        
    def set_transversity_params(self):
        self.set_constraits('transversity')

        conf['transversity']._widths1_uv  = conf['params']['transversity']['widths1_uv']['value']
        conf['transversity']._widths1_dv  = conf['params']['transversity']['widths1_dv']['value']
        conf['transversity']._widths1_sea = conf['params']['transversity']['widths1_sea']['value']

        conf['transversity']._widths2_uv  = conf['params']['transversity']['widths2_uv']['value']
        conf['transversity']._widths2_dv  = conf['params']['transversity']['widths2_dv']['value']
        conf['transversity']._widths2_sea = conf['params']['transversity']['widths2_sea']['value']

        iflav=0
        for flav in ['u','ub','d','db','s','sb']:
            iflav+=1
            ipar=-1
            for par in ['N0','a0','b0','c0','d0','N1','a1','b1','c1','d1']:
                ipar+=1
                if '%s %s 1'%(flav,par) in conf['params']['transversity']:
                    conf['transversity'].shape1[iflav][ipar] = conf['params']['transversity']['%s %s 1'%(flav,par)]['value']
                if '%s %s 2'%(flav,par) in conf['params']['transversity']:
                    conf['transversity'].shape2[iflav][ipar] = conf['params']['transversity']['%s %s 2'%(flav,par)]['value']

        conf['transversity'].setup()

    def set_sivers_params(self):
        self.set_constraits('sivers')

        conf['sivers']._widths1_uv  = conf['params']['sivers']['widths1_uv']['value']
        conf['sivers']._widths1_dv  = conf['params']['sivers']['widths1_dv']['value']
        conf['sivers']._widths1_sea = conf['params']['sivers']['widths1_sea']['value']

        conf['sivers']._widths2_uv  = conf['params']['sivers']['widths2_uv']['value']
        conf['sivers']._widths2_dv  = conf['params']['sivers']['widths2_dv']['value']
        conf['sivers']._widths2_sea = conf['params']['sivers']['widths2_sea']['value']

        iflav=0
        for flav in ['u','ub','d','db','s','sb']:
            iflav+=1
            ipar=-1
            for par in ['N0','a0','b0','c0','d0','N1','a1','b1','c1','d1']:
                ipar+=1
                if '%s %s 1'%(flav,par) in conf['params']['sivers']:
                    conf['sivers'].shape1[iflav][ipar] = conf['params']['sivers']['%s %s 1'%(flav,par)]['value']
                if '%s %s 2'%(flav,par) in conf['params']['sivers']:
                    conf['sivers'].shape2[iflav][ipar] = conf['params']['sivers']['%s %s 2'%(flav,par)]['value']

        conf['sivers'].setup()

    def set_boermulders_params(self):
        self.set_constraits('boermulders')

        conf['boermulders']._widths1_uv  = conf['params']['boermulders']['widths1_uv']['value']
        conf['boermulders']._widths1_dv  = conf['params']['boermulders']['widths1_dv']['value']
        conf['boermulders']._widths1_sea = conf['params']['boermulders']['widths1_sea']['value']

        conf['boermulders']._widths2_uv  = conf['params']['boermulders']['widths2_uv']['value']
        conf['boermulders']._widths2_dv  = conf['params']['boermulders']['widths2_dv']['value']
        conf['boermulders']._widths2_sea = conf['params']['boermulders']['widths2_sea']['value']

        iflav=0
        for flav in ['u','ub','d','db','s','sb']:
            iflav+=1
            ipar=-1
            for par in ['N0','a0','b0','c0','d0','N1','a1','b1','c1','d1']:
                ipar+=1
                if '%s %s 1'%(flav,par) in conf['params']['boermulders']:
                    conf['boermulders'].shape1[iflav][ipar] = conf['params']['boermulders']['%s %s 1'%(flav,par)]['value']
                if '%s %s 2'%(flav,par) in conf['params']['boermulders']:
                    conf['boermulders'].shape2[iflav][ipar] = conf['params']['boermulders']['%s %s 2'%(flav,par)]['value']

        conf['boermulders'].setup()

    def set_ffpi_params(self):
        self.set_constraits('ffpi')
        conf['ffpi']._widths1_fav  = conf['params']['ffpi']['widths1_fav']['value']
        conf['ffpi']._widths1_ufav = conf['params']['ffpi']['widths1_ufav']['value']
        conf['ffpi']._widths2_fav  = conf['params']['ffpi']['widths2_fav']['value']
        conf['ffpi']._widths2_ufav = conf['params']['ffpi']['widths2_ufav']['value']
        conf['ffpi'].setup()

    def set_ffk_params(self):
        self.set_constraits('ffk')
        conf['ffk']._widths1_fav   = conf['params']['ffk']['widths1_fav']['value']
        conf['ffk']._widths1_ufav  = conf['params']['ffk']['widths1_ufav']['value']
        conf['ffk']._widths2_fav   = conf['params']['ffk']['widths2_fav']['value']
        conf['ffk']._widths2_ufav  = conf['params']['ffk']['widths2_ufav']['value']
        conf['ffk'].setup()

    def set_collinspi_params(self):
        self.set_constraits('collinspi')
        conf['collinspi']._widths1_fav  = conf['params']['collinspi']['widths1_fav']['value']
        conf['collinspi']._widths1_ufav = conf['params']['collinspi']['widths1_ufav']['value']
        conf['collinspi']._widths2_fav  = conf['params']['collinspi']['widths2_fav']['value']
        conf['collinspi']._widths2_ufav = conf['params']['collinspi']['widths2_ufav']['value']

        iflav=0
        for flav in ['u','ub','d','db','s','sb']:
            iflav+=1
            ipar=-1
            for par in ['N0','a0','b0','c0','d0','N1','a1','b1','c1','d1']:
                ipar+=1
                if '%s %s 1'%(flav,par) in conf['params']['collinspi']:
                    conf['collinspi'].shape1[iflav][ipar] = conf['params']['collinspi']['%s %s 1'%(flav,par)]['value']
                if '%s %s 2'%(flav,par) in conf['params']['collinspi']:
                    conf['collinspi'].shape2[iflav][ipar] = conf['params']['collinspi']['%s %s 2'%(flav,par)]['value']
        conf['collinspi'].setup()

    def set_collinsk_params(self):
        self.set_constraits('collinsk')
        conf['collinsk']._widths1_fav   = conf['params']['collinsk']['widths1_fav']['value']
        conf['collinsk']._widths1_ufav  = conf['params']['collinsk']['widths1_ufav']['value']
        conf['collinsk']._widths2_fav   = conf['params']['collinsk']['widths2_fav']['value']
        conf['collinsk']._widths2_ufav  = conf['params']['collinsk']['widths2_ufav']['value']
        iflav=0
        for flav in ['u','ub','d','db','s','sb']:
            iflav+=1
            ipar=-1
            for par in ['N0','a0','b0','c0','d0','N1','a1','b1','c1','d1']:
                ipar+=1
                if '%s %s 1'%(flav,par) in conf['params']['collinsk']:
                    conf['collinsk'].shape1[iflav][ipar] = conf['params']['collinsk']['%s %s 1'%(flav,par)]['value']
                if '%s %s 2'%(flav,par) in conf['params']['collinsk']:
                    conf['collinsk'].shape2[iflav][ipar] = conf['params']['collinsk']['%s %s 2'%(flav,par)]['value']
        conf['collinsk'].setup()


