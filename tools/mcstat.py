#!/usr/bin/env python
import os
import sys
import time
import numpy as np
import pylab as py
import pandas as pd
from .tools import load, save
from .config import load_config,conf
from .inputmod import INPUTMOD

def chi2hist(runs):

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*5,nrows*3))
    
    ax=py.subplot(nrows,ncols,1)
    for k in runs:
        if k=='all': continue
        ax.hist(2*runs[k]['nll'],bins=50,range=None,histtype='step',label=str(k))
    ax.legend()
    ax.set_xlabel('chi2-npts')
    py.tight_layout()
    py.savefig('chi2hist.pdf')

class parhist:

    def __init__(self,runs,inputmod=None):

        self.inputmod=inputmod 
        self.order=self.get_ordered_free_params()
        self.tabs,self.tabs_a=self.get_tabs(runs)

        self.kind1=[]
        self.kind2=[]

        for _ in conf['params']:   self.kind1.append(_)
        for _ in conf['datasets']: self.kind2.append(_)    

    def get_ordered_free_params(self):
        order=[]

        for k in conf['params']:
            for kk in conf['params'][k]:
                if conf['params'][k][kk]['fixed']==False:
                    order.append([1,k,kk])

        if 'datasets' in conf:
            for k in conf['datasets']:
                for kk in conf['datasets'][k]['norm']:
                    if conf['datasets'][k]['norm'][kk]['fixed']==False:
                      order.append([2,k,kk])

        return order

    def get_tabs(self,runs):
        """
        create pandas data frame for the samples
        """
        tabs={}
        tabs_a={}
        for k in runs:
            tab,tab_a={},{}
            tab['nll']=runs[k]['nll']
            tab['weights']=runs[k]['weights']
            samples=np.transpose(runs[k]['samples'])
            active_p=np.transpose(runs[k]['active p'])
            for i in range(len(self.order)):
                _,kind,par=self.order[i]
                tab['%s:%s'%(kind,str(par))]=samples[i]
                tab_a['%s:%s'%(kind,str(par))]=active_p[i]
            tabs[k]=pd.DataFrame(tab)
            tabs_a[k]=pd.DataFrame(tab_a)
        return tabs,tabs_a

    def plot(self,tabs,tabs_a,entries,kind1,kind2,iRange=0):

        for i in range(len(entries)):
            self.cnt+=1
            if entries[i]==None: continue
            ax=py.subplot(self.nrows,self.ncols,self.cnt)
            kind,par=entries[i].split(':')
            for _ in kind1: 
              if kind==_:
                vmin=conf['params'][_][par]['min']
                vmax=conf['params'][_][par]['max']
                R=(vmin,vmax)
                E0=conf['params'][_][par]['value']
            for _ in kind2: 
              if kind==_:
                vmin=conf['datasets'][_]['norm'][int(par)]['min']
                vmax=conf['datasets'][_]['norm'][int(par)]['max']
                R=(vmin,vmax)
                E0=conf['datasets'][_]['norm'][int(par)]['value']
            if iRange==0:pass
            else: R=None         

            for _ in tabs:
                if _=='all': continue
                tab=tabs[_]
                ax.hist(tab[entries[i]],range=R,bins=100,weights=tab['weights'],\
                        histtype='step',label=str(_))   

            ax.hist(tabs['all'][entries[i]],range=R,bins=100,\
                    edgecolor='k',hatch='...',\
                    weights=tabs['all']['weights'],histtype='step',label='all')

            ax.plot(tabs_a['all'][entries[i]],np.zeros(tabs_a['all'][entries[i]].size),'ro')
            #ax.axvline(E0)
            ax.set_title('%s:%s'%(kind,par))

            if self.inputmod!=None:
              # update input.py            
              E=np.einsum('i,i',tabs['all']['weights'],tabs['all'][entries[i]])
              for _ in kind1: 
                if kind==_: self.inputmod.mod_par(_,par,'value',E)
              for _ in kind2: 
                if kind==_: self.inputmod.mod_norm(_,par,'value',E)

    def hist_widths(self,isave=False):

        entries=[]
        for kind in self.kind1:
            for par in conf['params'][kind]:
                for _ in ['widths1','widths2']:
                    if _ in par and conf['params'][kind][par]['fixed']==False:
                        entries.append('%s:%s'%(kind,par))

        if len(entries)==0: return 
        self.ncols=3
        self.nrows=len(entries)/self.ncols
        if len(entries)%self.ncols!=0: self.nrows+=1
        fig = py.figure(figsize=(self.ncols*3,self.nrows*1.5))
        self.cnt=0
        self.plot(self.tabs,self.tabs_a,entries,self.kind1,self.kind2,iRange=0)
        py.tight_layout()
        if isave:
            py.savefig('hist_widths.pdf')
            py.close()

    def hist_shape(self,isave=False):

        entries=[]
        for kind in self.kind1:
            entry=[]
            for par in sorted(conf['params'][kind]):
                if 'width' in par: continue
                for _ in ['N','a','b','c','d']:
                    if _ in par and conf['params'][kind][par]['fixed']==False:
                        entry.append('%s:%s'%(kind,par))
                    else:
                      entry.append(None)
            entries.extend(entry)
        if len(entries)==0: return 
        self.ncols=3
        self.nrows=len(entries)/self.ncols
        if len(entries)%self.ncols!=0: self.nrows+=1
        fig = py.figure(figsize=(self.ncols*3,self.nrows*1.5))
        self.cnt=0
        self.plot(self.tabs,self.tabs_a,entries,self.kind1,self.kind2,iRange=0)
        py.tight_layout()
        if isave:
            py.savefig('hist_shape.pdf')
            py.close()

    def hist_norm(self,isave=False):

        entries=[]
        for kind in self.kind2:
            for idx in conf['datasets'][kind]['norm']:
                if conf['datasets'][kind]['norm'][idx]['fixed']==False:
                    entries.append('%s:%s'%(kind,idx))

        if len(entries)==0: return 
        self.ncols=3
        self.nrows=len(entries)/self.ncols
        if len(entries)%self.ncols!=0: self.nrows+=1
        fig = py.figure(figsize=(self.ncols*3,self.nrows*1.5))
        self.cnt=0
        self.plot(self.tabs,self.tabs_a,entries,self.kind1,self.kind2,iRange=0)
        py.tight_layout()
        if isave:
            py.savefig('hist_norm.pdf')
            py.close()


















