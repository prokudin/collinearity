#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
from tools import isnumeric
from config import conf


class _READER:

    def __init__(self):
        pass

    def apply_cuts(self, tab,k):
        if 'filters' in conf['datasets'][self.reaction]:
            filters = conf['datasets'][self.reaction]['filters']
            for _ in filters: tab = tab.query(_)
        return tab

    def load_data_sets(self, reaction):
        self.reaction = reaction
        if reaction not in conf['datasets']:
            return None
        XLSX = conf['datasets'][reaction]['xlsx']
        TAB = {}
        for k in XLSX:
            sys.stdout.write('\r')
            sys.stdout.write('loading %s data sets %d' % (reaction, k))
            sys.stdout.flush()
            fname = conf['datasets'][reaction]['xlsx'][k]
            if fname.startswith('.'):
                tab = pd.read_excel(fname)
            else:
                tab = pd.read_excel('%s/database/%s' %
                                    (os.environ['JAM3D'], fname))
            tab = self.modify_table(tab,k)
            if len(tab.index)==0: continue
            TAB[k] = tab.to_dict(orient='list')
            for kk in TAB[k]:
                if len(TAB[k][kk]) == 0:
                    continue
                if isnumeric(TAB[k][kk][0]):
                    TAB[k][kk] = np.array(TAB[k][kk])

        return TAB
