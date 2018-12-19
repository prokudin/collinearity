import sys,os
import numpy as np
from tools.tools import checkdir,save,load
from tools.config import load_config, conf
import time
from scipy.optimize import minimize,leastsq
import json
from tools.bar import BAR
import itertools as it
import pandas as pd
from IPython.display import clear_output
from tools.inputmod import INPUTMOD
from scipy.optimize import least_squares

class ML:

  def __init__(self):
    if 'args' in conf:
      self.inputfile=inputfile=conf['args'].config
    self.cum_shifts=0
    self.iteration=0
    self.CHI2T=[]
    self.CHI2V=[]
    self.PARAMS=[]
    self.TOT=[]

  def get_stats(self,res,rres,nres,delay):

    shifts=conf['parman'].shifts
    etime = (time.time()-self.t0)/60

    npts=res.size
    chi2=np.sum(res**2)
    rchi2=np.sum(rres**2)
    nchi2=np.sum(nres**2)
    chi2tot=chi2+rchi2+nchi2
    dchi2=chi2tot-self.chi2tot
    if shifts>2: 
      if chi2tot<self.chi2tot:
        self.dchi2=self.chi2tot-chi2tot
        self.chi2tot=chi2tot


    status=[]
    status.append('JAM FITTER')
    status.append('count        = %d'%self.cnt)
    status.append('elapsed time(mins)=%f'%etime)
    status.append('shifts       = %d'%shifts)
    status.append('npts         = %d'%npts)
    dof = npts-len(conf['parman'].par)
    status.append('d.o.f.       = %d'%dof)
    status.append('chi2         = %f'%chi2)
    if dof>0:
      status.append('chi2/d.o.f   = %f'%(chi2/dof))
    status.append('rchi2        = %f'%rchi2)
    status.append('nchi2        = %f'%nchi2)
    status.append('chi2tot      = %f'%(chi2tot))
    status.append('dchi2(iter)  = %f'%self.dchi2)
    status.append('dchi2(local) = %f'%dchi2)
    status.append('')
    status.extend(conf['resman'].gen_report())
    parstatus = conf['parman'].gen_report()

    nstatus=len(status)
    nparstatus=len(parstatus)

    os.system('clear') 
    for i in range(max([nstatus,nparstatus])):
      data=[]
      if i<nstatus: data.append(status[i])
      else: data.append('')
      if i<nparstatus: data.append(parstatus[i])
      else: data.append('')
      print '%-120s  | %s'%tuple(data)
      
    
    if delay==True: 
      if 'args' in conf: self.gen_output()
      self.gen_output()
      #sys.exit()

    #clear_output(wait=True)
    #if conf['screen mode']=='curses':

    #  self.status.append('')
    #  if    delay==False: self.status.append('')
    #  else: self.status.append('Its over!!!')
    #  self.status.append('')
    #  self.status.append('press q to exit')

    #  conf['screen'].clear()
    #  conf['screen'].border(0)
    #  if delay==False: conf['screen'].nodelay(1)
    #  else: conf['screen'].nodelay(0)

    #  for i in range(len(self.status)):
    #    conf['screen'].addstr(i+2,2,self.status[i])

    #  for i in range(len(parstatus)):
    #    conf['screen'].addstr(i+2,80,parstatus[i])


    #  conf['screen'].refresh()
    #  if conf['screen'].getch()==ord('q'):
    #    curses.endwin()
    #    self.gen_output()
    #    sys.exit()

    #elif conf['screen mode']=='plain':

  def get_residuals(self,par,delay=False,status=True):
    res,rres,nres=conf['resman'].get_residuals(par)
    if status: 
      self.cnt+=1
      self.get_stats(res,rres,nres,delay)
    if len(rres)!=0: res=np.append(res,rres)
    if len(nres)!=0: res=np.append(res,nres)
    return res

  def gen_output(self):

    inputmod=INPUTMOD(conf['args'].config)

    for kind in conf['params']: 
      for par in conf['params'][kind]:
        value=conf['params'][kind][par]['value']
        inputmod.mod_par(kind,par,'value',value) 

    for reaction in conf['datasets']:
      for idx in conf['datasets'][reaction]['norm']:
        value=conf['datasets'][reaction]['norm'][idx]['value']
        inputmod.mod_norm(reaction,idx,'value',value) 

    inputmod.gen_input()

  def get_chi2(self,par,status=True):
    return np.sum(self.get_residuals(par,status=status)**2)

  def run_minimize(self):

    guess=conf['parman'].par
    order=conf['parman'].order
    conf['screen mode']='plain'

    bounds=[]
    for entry in order:
      i,k,kk=entry
      if i==1:
        bounds.append([conf['params'][k][kk]['min'],conf['params'][k][kk]['max']])
      elif i==2:
        bounds.append([None,None])

    if conf['screen mode']=='curses':
      conf['screen']=curses.initscr()

    self.chi2tot=1e1000
    self.dchi2=0
    self.t0 = time.time()
    self.cnt=0

    #self.get_residuals(guess)
    #sys.exit()

    res = minimize(self.get_chi2,guess, bounds=bounds,method='TNC')
    res=self.get_residuals(res.x,delay=True)

  def run_leastsq(self):

    guess=conf['parman'].par
    order=conf['parman'].order
    conf['screen mode']='plain'


    bounds=[]
    for entry in order:
      i,k,kk=entry
      if i==1:
        bounds.append([conf['params'][k][kk]['min'],conf['params'][k][kk]['max']])
      elif i==2:
        bounds.append([None,None])

    if conf['screen mode']=='curses':
      conf['screen']=curses.initscr()

    self.chi2tot=1e1000
    self.dchi2=0
    self.t0 = time.time()
    self.cnt=0

    #self.get_residuals(guess)
    #self.gen_output()
    #sys.exit()

    #print conf['lmdiff tol']
    #sys.exit()

    if 'lmdiff tol' in conf:
      fit=leastsq(self.get_residuals,guess,full_output = 1,ftol=conf['lmdiff tol'])
    else:
      fit=leastsq(self.get_residuals,guess,full_output = 1)
    res=self.get_residuals(fit[0],delay=True)

  def run_leastsq2(self):

    if 'flat par' in conf and conf['flat par']==True:
      guess=conf['parman'].gen_flat()
    else:
      guess=conf['parman'].par

    order=conf['parman'].order
    #conf['screen mode']='plain'

    bounds_min=[]
    bounds_max=[]
    for entry in order:
      i,k,kk=entry
      if i==1:
        p=conf['params'][k][kk]['value']
        pmin=conf['params'][k][kk]['min']
        pmax=conf['params'][k][kk]['max']
        if p<pmin or p>pmax: 
            raise ValueError('%s/%s outsize the limits'%(k,kk))
        bounds_min.append(conf['params'][k][kk]['min'])
        bounds_max.append(conf['params'][k][kk]['max'])
      elif i==2:
        bounds_min.append(conf['datasets'][k]['norm'][kk]['min'])
        bounds_max.append(conf['datasets'][k]['norm'][kk]['max'])

    self.chi2tot=1e1000
    self.dchi2=0
    self.t0 = time.time()
    self.cnt=0
    fit = least_squares(self.get_residuals, guess,bounds=(bounds_min,bounds_max),method='trf',ftol=conf['ftol'])
    res=self.get_residuals(fit.x,delay=True)

    if 'bootstrap' in conf and conf['bootstrap']:
      fname='%s.mcr'%id_generator(size=12)
      save(fit.x,fname)
      if 'cmd' in conf:
        try:
          os.system(conf['cmd'].replace('<<fname>>',fname))
        except:
          print 'could not execute %s'%cmd

  def run_test(self):

    guess=conf['parman'].par
    order=conf['parman'].order
    conf['screen mode']='plain'

    bounds=[]
    for entry in order:
      i,k,kk=entry
      if i==1:
        bounds.append([conf['params'][k][kk]['min'],conf['params'][k][kk]['max']])
      elif i==2:
        bounds.append([None,None])

    if conf['screen mode']=='curses':
      conf['screen']=curses.initscr()

    self.chi2tot=1e1000
    self.dchi2=0
    self.t0 = time.time()
    self.cnt=0
    self.get_residuals(guess)

    for l in conf['resman'].gen_report(verb=0,level=1):
      print l

  def analysis(self):
    self.gen_report()

  def gen_report(self):
    inputfile=conf['args'].config
    outdir='outputs/'+inputfile.split('/')[-1].replace('.py','')
    checkdir(outdir)
    par=conf['parman'].par
    self.t0 = time.time()
    self.cnt=0
    self.get_residuals(par,status=False)
    print conf['pdf'].sr
    report=conf['resman'].gen_report(verb=1,level=1)
    save(report,'%s/report-ML'%outdir)

  def simulation(self):
    par=conf['parman'].par
    self.get_residuals(par,delay=False,status=False)
    for obs in conf['datasets']:
      for idx in conf['datasets'][obs]['xlsx']:
        if idx<0: 
          fname=conf['datasets'][obs]['xlsx'][idx]
          tab=pd.read_excel(fname)
          tab.value=conf['%s tabs'%obs][idx]['thy']
          tab['stat_u']=tab.value*tab['stat(%)']
          writer = pd.ExcelWriter(fname)
          tab.to_excel(writer,'Sheet1',index=False)
          writer.save()

  def simulation_json(self):
    resman=conf['resman']
    parman=conf['parman']
    inputfile=conf['args'].config
    
    par=parman.par
    parman.set_new_params(par)


    M2=conf['aux'].M2
    M=np.sqrt(M2)

    data={}
    data['author']="N. Sato"
    data['generator']="JAM"
    data['reaction']="DIS"
    data['target']="p"
    data['Elab']="10.6"
    data['lepton']="e-"
    data['variables']=["x,y,Q2,F2,FL,FL,dsig/dxdy"]

    Elab=float(data['Elab'])
    xmin=1/(2*M*Elab)
    ymin=1/(2*M*Elab)

    data['axis']=[]
    data['axis'].append({ "name": "a", "bins":  200, "min":  xmin, "max":   0.999, "scale":"arb" ,"description":"Bjorken x"})
    data['axis'].append({ "name": "b", "bins":  200, "min":  ymin, "max":   0.999, "scale":"arb", "description":"y"},)

    Tab=json.dumps(data,sort_keys=True,indent=4, separators=(',', ': '))
    Tab=[Tab]
    Tab.append('%10s %10s %10s %10s %10s %10s %10s %10s %10s'%('ix','iy','x','y','Q2','F2','FL','F3','dsig/dxdy'))

    #Tab=["#!"+Tab]

    bin_center = lambda xmin,xmax,n: np.array([xmin + (2*i+1)*(xmax-xmin)/(2*n) for i in range(n)])

    xs =bin_center(data['axis'][0]['min'],data['axis'][0]['max'],data['axis'][0]['bins'])
    ys =bin_center(data['axis'][1]['min'],data['axis'][1]['max'],data['axis'][1]['bins'])

    I=range(xs.size)
    print 
    bar=BAR('generating json file',xs.size**2)
    for item in it.product(I,I):
      ix,iy=item

      x=xs[ix]
      y=ys[iy]
      Q2=x*y*(2*M*Elab)
      if Q2<1: continue

      F2=conf['stfuncs'].get_FXN(x,Q2,'F2',nucleon='proton',twist=2,tmc=False,evolve=True)
      FL=conf['stfuncs'].get_FXN(x,Q2,'FL',nucleon='proton',twist=2,tmc=False,evolve=True)
      F3=conf['stfuncs'].get_FXN(x,Q2,'F3',nucleon='proton',twist=2,tmc=False,evolve=True)

      if np.isnan(F2): continue
      if np.isnan(FL): continue
      if np.isnan(F3): continue
      
      s=M2+2*Elab*M2**0.5
      y=(Q2/2/x)/((s-M2)/2)

      if y<0 or y>1: 
        bar.next()
        continue 
        
      YP=1+(1-y)**2
      YM=1-(1-y)**2
      if   data['lepton']=="e-": sign=1
      elif data['lepton']=="e+": sign=-1

      xsec=2*np.pi*conf['aux'].alfa**2/x/y/Q2
      xsec*=(YP+2*x**2*y**2*M2/Q2)*F2-y**2*FL+sign*YM*x*F3

      if xsec<0: continue

      row='%10d %10d %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e'
      row=row%(ix,iy,x,y,Q2,F2,FL,F3,xsec)
      Tab.append(row)
      bar.next()
    bar.finish()
    print '\nsaving...'
    Tab=[row+'\n' for row in Tab]
    f = open("dis-jam.json",'w')
    f.writelines(Tab)
    f.close()

  def simulation_json2(self):
    resman=conf['resman']
    parman=conf['parman']
    inputfile=conf['args'].config
    
    par=parman.par
    parman.set_new_params(par)


    M2=conf['aux'].M2
    M=np.sqrt(M2)

    data={}
    data['author']="N. Sato"
    data['generator']="JAM"
    data['reaction']="DIS"
    data['target']="p"
    data['Elab']="10.6"
    data['lepton']="e-"
    data['variables']=["Counts","Err.Counts","xav","yav","q2av"]

    Elab=float(data['Elab'])
    xmin=1/(2*M*Elab)
    ymin=1/(2*M*Elab)

    data['axis']=[]
    data['axis'].append({ "name": "a", "bins":  100, "min":  xmin, "max":   0.999, "scale":"arb" ,"description":"Bjorken x"})
    data['axis'].append({ "name": "b", "bins":  50, "min":  ymin, "max":   0.999, "scale":"arb", "description":"y"},)

    Tab=json.dumps(data,sort_keys=True,indent=4, separators=(',', ': '))
    Tab=[Tab]

    #Tab=["#!"+Tab]

    bin_center = lambda xmin,xmax,n: np.array([xmin + (2*i+1)*(xmax-xmin)/(2*n) for i in range(n)])

    xs =bin_center(data['axis'][0]['min'],data['axis'][0]['max'],data['axis'][0]['bins'])
    ys =bin_center(data['axis'][1]['min'],data['axis'][1]['max'],data['axis'][1]['bins'])

    IX=range(xs.size)
    IY=range(ys.size)
    print 
    bar=BAR('generating json file',xs.size**2)
    for item in it.product(IX,IY):
      ix,iy=item

      x=xs[ix]
      y=ys[iy]
      Q2=x*y*(2*M*Elab)
      if Q2<1: continue

      F2=conf['stfuncs'].get_FXN(x,Q2,'F2',nucleon='proton',twist=2,tmc=False,evolve=True)
      FL=conf['stfuncs'].get_FXN(x,Q2,'FL',nucleon='proton',twist=2,tmc=False,evolve=True)
      F3=conf['stfuncs'].get_FXN(x,Q2,'F3',nucleon='proton',twist=2,tmc=False,evolve=True)

      if np.isnan(F2): continue
      if np.isnan(FL): continue
      if np.isnan(F3): continue
      
      s=M2+2*Elab*M2**0.5
      y=(Q2/2/x)/((s-M2)/2)

      if y<0 or y>1: 
        bar.next()
        continue 
        
      YP=1+(1-y)**2
      YM=1-(1-y)**2
      if   data['lepton']=="e-": sign=1
      elif data['lepton']=="e+": sign=-1

      xsec=2*np.pi*conf['aux'].alfa**2/x/y/Q2
      xsec*=(YP+2*x**2*y**2*M2/Q2)*F2-y**2*FL+sign*YM*x*F3

      xsec*=1e8

      if xsec<0: continue

      row='%10d %10d %10.4e %10.4e %10.4e %10.4e %10.4e '
      row=row%(ix,iy,xsec,xsec*0.05,x,y,Q2)
      Tab.append(row)
      bar.next()
    bar.finish()
    print '\nsaving...'
    Tab=[row+'\n' for row in Tab]
    f = open("jam.dat",'w')
    f.writelines(Tab)
    f.close()

  def bootstrap(self):

    if 'args' in conf:
      outputdir='%s/mcdata'%conf['args'].config.split('/')[-1].replace('.py','')
    checkdir(outputdir)

    guess=conf['parman'].par
    order=conf['parman'].order
    conf['screen mode']=None
    npar=len(guess)
    nruns=conf['num replicas']
    bar=BAR('running bootstrap run',nruns)
    for i in range(nruns):
      conf['parman'].get_ordered_free_params()
      conf['resman'].resample()
      self.chi2tot=1e1000
      self.dchi2=0
      self.t0 = time.time()
      self.cnt=0
      if 'lmdiff tol' in conf:
        fit=leastsq(lambda p: self.get_residuals(p,False,True),guess,full_output = 1,ftol=conf['lmdiff tol'])
      else:
        fit=leastsq(lambda p: self.get_residuals(p,False,True),guess,full_output = 1)
      save(fit[0],'%s/rep%d.dat'%(outputdir,i))
      bar.next()
    bar.finish()

  def hessian(self):
    p0=conf['parman'].par
    npar=len(conf['parman'].par)
    #idel=[i for i in range(npar) if conf['parman'].order[i][0]==2]
    #H=np.delete(H,idel,axis=0)
    #H=np.delete(H,idel,axis=1)
    hess=HESS(p0,lambda p: self.get_chi2(p,status=False)).data
    if 'args' in conf:
      outputdir='%s/hess'%conf['args'].outdir
      checkdir(outputdir)
      fname=conf['args'].config.split('/')[-1].replace('.py','')
      save(hess,'%s/%s.par'%(outputdir,fname))
    else:
      return hess
     



