#!/usr/bin/env python
import os
import sys
import time
import argparse
from .tools import load, save
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def get_runs(path2mc,veto=[]):

  # list mc file names
  F=os.listdir(path2mc)
  F=[f for f in F if f.endswith('mc') if '-' in f]
  if len(F)==0:
    raise ValueError('*.mc files not found')

  # split mc files by runs
  names=list(set([f.split('-')[0] for f in F]))

  # combine blocks for each run
  runs={}
  cnt=0
  for name in names:
    run_files=sorted([f for f in F if name in f])
    run_data=[load('%s/%s'%(path2mc,f)) for f in run_files]
    nap=run_data[0]['num active points']
    samples=[] 
    nll=[]
    for _ in run_data: 
      for p in _['samples']: samples.append(p)
      nll.extend(_['nll'])
    active_p=run_data[-1]['active p']
    nll=np.array(nll)
    samples=np.array(samples)
    runs[cnt]={'samples':samples,'nll':nll,'num active points':nap,'name':name,'active p':active_p}
    cnt+=1

  # remove veto runs
  if veto!=[]:
    for _ in veto:
      del runs[_]

  # combine samples
  nap=0
  samples=[] 
  active_p=[] 
  nll=[]
  for k in runs:
    for p in runs[k]['samples']: samples.append(p)
    for p in runs[k]['active p']: active_p.append(p)
    nll.extend(runs[k]['nll'])
    nap+=runs[k]['num active points'] 
  nll=np.array(nll)
  samples=np.array(samples)
  runs['all']={'samples':samples,'nll':nll,'num active points':nap,'name':'all','active p':active_p}
  print 'original samples size =',len(samples)
  return runs

def get_ordered_samples(nap,nll,samples):
  nllc=np.copy(nll)
  I=np.argsort(nll)
  nll=nll[I]
  # regularize
  min_nll=np.amin(nll)
  samples=samples[I]
  likelihood=np.exp(-(nll-min_nll)) 
  x=np.array([((nap-1.)/nap)**i for i in range(likelihood.size+1)])[::-1]
  dx=(0.5*(x[1:]-x[:-1]))
  weights=likelihood*dx
  weights/=np.sum(weights)
  samples=np.array([samples[i] for i in range(weights.size) if weights[i]>0])
  nll=np.array([nll[i] for i in range(weights.size) if weights[i]>0])
  weights=np.array([weights[i] for i in range(weights.size) if weights[i]>0])
  weights/=np.sum(weights)
  return nll,weights,samples

def impose_cdf_cut(mcdata,cdfcut):
  weights=np.copy(mcdata['weights'])
  samples=np.copy(mcdata['samples'])
  I=np.argsort(weights)
  weights=weights[I]
  samples=samples[I]

  cdf=[]
  for i in range(weights.size):
    cdf.append(np.sum(weights[:i+1]))
  II=[i  for i in range(weights.size)  if cdf[i]>cdfcut]
  weights=weights[II]
  weights/=np.sum(weights)
  samples=samples[II]

  return weights,samples

