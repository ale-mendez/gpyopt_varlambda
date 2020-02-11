#!/usr/bin/env python3

# Alejandra Mendez - 28/04/2019
# v1.2
#
#  * this program implements hyperoptimization of NLAMVAR parameters in
#     autosctructure code (max=20)
#  * the autovarlambda.f90 subroutine is compiled with f2py and incorporated
#     to this code
#  * the autostructure code is executed using system
#  * the number of configurations defined by "cfgs" in .yaml is considered 
#     in the computation
#  * multiple excited energies can be considered via "nener"
#  * the number of elements in array "weight" must be equal
#     to the amount of energies "nener" considered
#
#  * input variables:
#           cfgs     -- number of configurations to be included in the calculation
#           initer   -- initial number of evaluations (first prior)
#           maxevals -- number of evaluations in minimization
#           mtype    -- model type
#                         ="GP" : Gaussian Process
#           aftype   -- acquisition function type
#                         ="EI"     : Expected improvement
#                         ="MPI"    : Maximum probability of improvement
#                         ="GP-UCB" : Upper confidence bound
#           afweight -- acquisition weight (float)
#           gridtype -- grid type
#                         ="continuous"
#           nlamvar  -- number of lambda parameters to be varied
#           maxlam   -- maximum lambda value
#           minlam   -- minimum lambda value
#           mfunc    -- minimization function
#                         ="Er"    : sum of weighted relative errros
#                         ="Er**2" : sum of weighted square relative errors
#           minst    -- states to be included in minimization
#                         ="gr+ex" : ground and excited states
#                         ="ex"    : only excited states
#           nener    -- number of states to be considered in minimization
#           wi       -- minimization weight
#                         ="eq"  : all elements are weighted equally
#                         ="gi"  : use statistical weight gi=sum 2j+1
#                         ="inp" : read input values (below)
#           weight   -- weight in relative errors (number of elements == nener)
#
#
#  * input files:
#           varlam_gpyopt.yml     -- input file for .py
#           das_XXCFG               -- autostructure input
#           exactvalues.dat         -- exact/experimental energy data
#           f2py_autovarlambda.f90  -- subroutines for variation of lambda
#

# Quietly compile autovarlamba.f90
import os
import sys
# os.system("f2py -c -m autovarlambda SR.read_write_as.for f2py_autovarlambda.f90 --quiet")
# print("\n===>>> autovarlambda compiled with f2py <<<===")

# Import other stuff
import numpy as np
import autovarlambda
import time
import pandas as pd
import yaml
import GPy
import GPyOpt
from numpy.random import seed
#seed(123456)


import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import gridspec

colors=['k','r','g','b','m','y','c','tab:blue',
        'tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray']
orb=['1s','2s','2p','3s','3p','3d','4s','4p','4d','4f','5s','5p','5d','5f','5g']
symbs=['o','s','v','^','>','<','o','o','o','o','o','o','o','o','o']

# Input parameters
global gridtype, cfgs, maxevals, mtype, aftype, afweight, gridtype, igrid, initer, ipol
global nlamvar, maxlam, minlam, mfunc, ifun, minst, imin, nener, wi, igi, iseed
global ii,npolvar

#########################################################################################################
# Define functions
#########################################################################################################

def data_input():
    global gridtype, cfgs, maxevals, mtype, aftype, afweight, gridtype, igrid, initer, ipol
    global nlamvar, maxlam, minlam, mfunc, ifun, minst, imin, nener, wi, igi, iseed
    print("\n===>>> Input parameters <<<===\n")
    with open("inpvar_gpyopt.yml", 'r') as stream:
        data_loaded=yaml.load(stream)
    cfgs=data_loaded.get("cfgs")
    initer=data_loaded.get("initer")
    maxevals=data_loaded.get("maxevals")
    mtype=data_loaded.get("mtype")
    aftype=data_loaded.get("aftype")
    afweight=data_loaded.get("afweight")
    gridtype=data_loaded.get("gridtype")
    nlamvar=data_loaded.get("nlamvar")
    maxlam=data_loaded.get("maxlam")
    minlam=data_loaded.get("minlam")
    mfunc=data_loaded.get("mfunc")
    minst=data_loaded.get("minst")
    nener=data_loaded.get("nener")
    wi=data_loaded.get("wi")
    vpol=data_loaded.get("vpol")
    iseed=data_loaded.get("iseed")

    if cfgs is None: raise ValueError("Error: cfgs not defined.")
    if initer is None: raise ValueError("Error: initer not defined.")
    if maxevals is None: raise ValueError("Error: maxevals not defined.")
    if mtype is None: 
        mtype="GP"
        print("Warning: model type not defined. Use default: GP")
#    if mtype!="GP": raise ValueError("Error: "+mtype+" not implemented.")
    if aftype is None: 
        aftype="EI"
        print("Warning: acquisition function type not defined. Use default: EI")
    if aftype!="EI" and aftype!="MPI" and aftype!="GP-UCB": raise ValueError("Error: "+aftype+" not implemented.")
    if afweight is None: 
        afweight=1.
        print("Warning: acquisition function weight not defined. Use default: 1.")
    if gridtype=="continuous":
        igrid=0
    elif gridtype is None: 
        igrid=0
        print("Warning: grid type not defined. Use default: continuous")
    else: raise ValueError("Error: "+gridtype+" not implemented.")
    if nlamvar is None: raise ValueError("Error: nlamvar not defined.")
    if maxlam is None or minlam is None: raise ValueError("Error: maxlam and/or minlam not defined.")
    if mfunc=="Er":
        ifun=0
    elif mfunc=="Er**2":
        ifun=1
    elif mfunc is None: 
        imin=0
        print("Warning: minimization function not defined. Use default: Er ")
    else: raise ValueError("Error: "+mfunc+" not implemented.")
    if minst=="gr+ex":
        imin=0
    elif minst=="ex":
        imin=1
    elif minst is None: 
        ifun=0
        print("Warning: states to be minimized not defined. Use default: ground + excited")
    else: raise ValueError("Error: "+minst+" not implemented.")
    if nener is None: raise ValueError("Error: nener not defined.")
    if wi=="eq" or wi=="gi":
        if wi=="eq": igi=0
        if wi=="gi": igi=1
        weight=np.full(nener,1.)
    elif wi=="inp":
        inpweight=data_loaded.get("weight")
        linpw=len(inpweight)
        if linpw<nener: # fill the missing wi with the last input value
            for i in range(linpw,nener):
                inpweight.append(inpweight[linpw-1])
        weight=np.array(listweight)
        if inpweight is None or linpw==0:
            raise ValueError("Error: weight not defined.")
    elif wi is None:
        igi=1
        print("Warning: wi not defined. Use default: eq")
    else: raise ValueError("Error: "+wi+"not implemented.")
    ipol=0
    if vpol==True: 
        ipol=1
        autovarlambda.polbck.ipol=ipol
    if iseed is None: 
        iseed=123456.
        raise ValueError("Warning: iseed = 123456 ")

def print_input():
    global gridtype, cfgs, maxevals, mtype, aftype, afweight, gridtype, igrid, initer, ipol
    global nlamvar, maxlam, minlam, mfunc, ifun, minst, imin, nener, wi, igi, iseed
    print(" * initial evaluations: "+str(initer))
    print(" * evaluations: "+str(maxevals))
    print(" * model type: "+mtype)
    print(" * acquisition function type: "+aftype)
    print(" * acquisition function weight: "+str(afweight))
    print(" * grid type: "+gridtype)
    if ifun==0: print(" * minimization function: sum of weighted relative errors")
    if ifun==1: print(" * minimization function: sum of weighted square relative errors")
    if imin==0: print(" * states included in minimization: ground + excited")
    if imin==1: print(" * states included in minimization: excited")
    print(" * weight type: "+wi)
    if ipol==1: print(" * include polarization ON.")
    if ipol==0: print(" * no polarization.")
    print(" * iseed = ",iseed)

def atom_name(nzion):
    nzion=abs(nzion)
    atomname=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg']
    nname=len(atomname)
    if nzion>nname: raise ValueError("atomname not defined for Z="+str(nzion))
    for i in range(nname):
        if nzion==i+1: chatom=atomname[i] 
    return chatom

def pot_type(nzion):
    if nzion < 0:
        chtype="STO"
    else:
        chtype="TFDA"
    return chtype

# Define function for writing iteration number in a human-readable-way
def human_format(num):
    magnitude=0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000
    return '%.f%s' % (num, ['', 'K', 'M', 'G', 'T', 'P'][magnitude])

# Determine file name string vector from fortran subroutine data
def def_filename():
    global ncfg, filename,initer
    ncfg=autovarlambda.salgebbck.mxconf
    nzion=autovarlambda.sminimbck.nzion
    chatom=atom_name(nzion)
    chtype=pot_type(nzion)
    cfgs=str(ncfg)+"CFG"+chtype
    idrun=str(initer)+mtype+human_format(maxevals)
    filename=chatom+cfgs+"_"+idrun

# Open .out file
def open_fout():
    global fener
    fener=open(filename+".out","w")

# Define and initialize variables
def init_var():
    enere=np.zeros(nener)
    enerc=np.zeros(nener)
    neex=int()
    return enere, enerc, neex


############  SPACES  ############
def mylamspace(minlam,maxlam):
    space=[]
    for i in range(nlamvar):
        alam=str(i+1)
        space.append({'name': 'lam'+alam, 'type': gridtype, 'domain': (minlam,maxlam)})
    return space

def mypolspace(alfd_min,alfd_max,rcut_min,rcut_max,npolvar):
    space=[]
    for i in range(npolvar):
        apol=str(i)
        space.append({'name': 'alfd'+apol,'type': 'continuous', 'domain': (alfd_min,alfd_max)})
        space.append({'name': 'rcut'+apol,'type': 'continuous', 'domain': (rcut_min,rcut_max)})
    return space

############  KERNELS  ############
def mykernel(rbf,stdperiodic):
    if rbf==1: kernel=GPy.kern.RBF(input_dim=nlamvar,variance=varf,lengthscale=lf,ARD=ARD)
    if rbf==1 and stdperiodic==1:
        kernel1=GPy.kern.RBF(input_dim=nlamvar,variance=varf,lengthscale=lf,ARD=ARD)
        kernel2=GPy.kern.StdPeriodic(input_dim=nlamvar,variance=varf,period=None,lengthscale=None,ARD1=ARD,ARD2=ARD)
        kernel=kernel1*kernel2
    return kernel

def mypolkernel(rbf,stdperiodic):
    if rbf==1: kernel=GPy.kern.RBF(input_dim=2*npolvar,variance=varf,lengthscale=lf,ARD=ARD)
    return kernel
###################################

def error_relat(valexact,valcomp,neex):
    if valexact == 0.0 or valcomp == 0.0:
        error=5./neex
    else:
        error=abs((valexact-valcomp)/valexact)*100
    return error

def loss_sumwE():
    global imin
    neex=autovarlambda.eei_ls.ne
    enere=autovarlambda.eicompare.enere
    enerc=autovarlambda.eicompare.enerc
    loss=0.0
    if igi==1: weight=autovarlambda.cei_ls.gic
    if igi==0: weight=autovarlambda.cei_ls.gic*0+1.
    for i in range(imin,neex): # consider ground energy and/or-only excited states
        erp=error_relat(enere[i],enerc[i],neex)
        if ifun==1: erp=erp*erp
        loss=loss + weight[i]*erp
#         print(' {} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}'.format(i,enere[i],enerc[i],erp,weight[i],loss))
    return loss

def loss_sumaki():
    global imin
    ntrtot=autovarlambda.akicompare.ntrtot
    vakie=autovarlambda.akicompare.vakie
    vakic=autovarlambda.akicompare.vakic
    loss=0.0
    for i in range(imin,ntrtot): # consider a certain number of transitions
        erp=error_relat(vakie[i],vakic[i],ntrtot)
        if ifun==1: erp=erp*erp
        loss=loss+erp
#        print(' {} {:.3e} {:.3e} {:.4f} {:.4f}'.format(i,vakie[i],vakic[i],erp,loss))
    return loss
    
def loss_total(ii,x):
    global icost
    lx=len(x)
    lei=loss_sumwE()
    xp='{:4d}: [ '
    for i in range(lx):
        xp=xp+'{vp['+str(i)+']:.{prec}f} '
    xp=xp+'] // {l1:.{prec}f}'
    if icost==0:
        loss=lei
#        print(xp.format(ii,vp=x,l1=lei,prec=4),file=fener)
        print(xp.format(ii,vp=x,l1=lei,prec=4))
    if icost==1:
        laki=loss_sumaki()
        loss=lei+laki
        xp=xp+' + {l2:.{prec}f} = {lt:.{prec}f}'
#        print(xp.format(ii,vp=x,l1=lei,l2=laki,lt=loss,prec=4),file=fener)
        print(xp.format(ii,vp=x,l1=lei,l2=laki,lt=loss,prec=4))
    return loss

def check_errlog():
    ierr=0
    errlog=os.path.exists("errlog")
    if errlog == True: 
        ierr=1
        os.system("rm errlog")
    autovarlambda.errlogbck.ierr=ierr
    return

def py_runAS_lam(lam,nlamvar):
    check_errlog()
    autovarlambda.run_varlam(lam,nlamvar)
    check_errlog()
    autovarlambda.inp_comp()
    autovarlambda.compare_ei()
    autovarlambda.compare_aki()
    return

def py_runAS_pol(pol,npolvar):
    check_errlog()
    autovarlambda.run_varpol(pol,npolvar)
    check_errlog()
    autovarlambda.inp_comp()
    autovarlambda.compare_ei()
    autovarlambda.compare_aki()
    return

def var_lambda(x):
    global ii
    ii=ii+1
    lam=x[0,:]
    py_runAS_lam(lam,nlamvar)
    loss=loss_total(ii,lam)
    return loss

def var_pol(x):
    global ii,npolvar
    pol=x[0,:]
    ii=ii+1
    py_runAS_pol(pol,npolvar)
    loss=loss_total(ii,pol)
    return loss

def run_mybo_pol(iseed,space,kernel,initer,maxevals,varf,lf,noise_var,exact_feval,optimize_restarts,ARD,xi):
    seed(iseed)
    enere, enerc, neex=init_var()
    cpdas="cp das_" + str(cfgs) + "CFG das"
    os.system(cpdas)
    autovarlambda.open_files()
    autovarlambda.inp_das()
    autovarlambda.inp_obs()
    dummyne=autovarlambda.eei_ls.ne.copy()
    if nener != dummyne:
        autovarlambda.eei_ls.ne=nener
    def_filename()
    if ncfg != cfgs:
        print(" >>>> configuration mismatch !!!! ")
        sys.exit()
    open_fout()
    t0=time.time()
    global ii
    ii=0
    model=GPyOpt.models.GPModel(kernel=kernel,noise_var=noise_var,exact_feval=exact_feval,
                                optimizer='lbfgs',max_iters=1500,optimize_restarts=optimize_restarts,
                                verbose=False,ARD=ARD)
    dspace=GPyOpt.Design_space(space=space)
    objective=GPyOpt.core.task.SingleObjective(var_pol)
    initial_design=GPyOpt.experiment_design.initial_design('random',dspace,initer)
#    initial_design=GPyOpt.experiment_design.initial_design('latin',dspace,initer)
    acquisition_optimizer=GPyOpt.optimization.AcquisitionOptimizer(dspace, optimizer='lbfgs')
    acquisition=GPyOpt.acquisitions.AcquisitionEI(model,dspace,acquisition_optimizer,jitter=xi)
    evaluator=GPyOpt.core.evaluators.Sequential(acquisition)
    myBopt=GPyOpt.methods.ModularBayesianOptimization(model,dspace,objective,acquisition,
                                                      evaluator,X_init=initial_design,normalize_Y=True,
                                                      model_update_interval=1)
    myBopt.run_optimization(maxevals,eps=1e-7,verbosity=False)
    myBopt.plot_acquisition("polparam_"+str(initer)+"GP"+str(maxevals)+".eps")
    myBopt.save_evaluations(filename+".dat")
    t1=time.time()
    total=(t1-t0)/60.
    print_optpol(total,myBopt)
    return 

def print_optpol(total,myBopt):
    global fmin,weight,fener,npolvar,initer,maxevals,filename,iseed
    pol_best=myBopt.x_opt
    fmin=myBopt.fx_opt
    py_runAS_pol(pol_best,npolvar)
    neex=autovarlambda.eicompare.neex
    enere=autovarlambda.eicompare.enere
    enerc=autovarlambda.eicompare.enerc
    loss_best=loss_sumwE()
    erelat=[]
    erroregr=error_relat(enere[0],enerc[0],neex)
    erelat.append(erroregr)
    xalpha='\n alpha = '
    xrcut='  rcut = '
    ii=0
    for i in range(npolvar):
        j=i+ii
        xalpha=xalpha+'{vp['+str(j)+']:.{prec}f}'
        xrcut=xrcut+'{vp['+str(j+1)+']:.{prec}f}'
        if j<npolvar: 
            xalpha=xalpha+' , '
            xrcut=xrcut+' , '
        ii=ii+1

    print("\n    Initial iterations = {:5d}".format(initer),file=fener)
    print(" Number of evaluations = {:5d}".format(maxevals),file=fener)
    print("           Random seed = {:8d}".format(iseed),file=fener)
    print("            Total time = {:.{prec}f} min".format(total,prec=4),file=fener)
    print("-"*80,file=fener)
    print(" Best results:",file=fener)
    print("-"*80,file=fener)
    print(xalpha.format(vp=pol_best,prec=4),file=fener)
    print(xrcut.format(vp=pol_best,prec=4),file=fener)
    for i in range(1,neex):
        errorex=error_relat(enere[i],enerc[i],neex)
        erelat.append(errorex)
    print("\n Ji:",file=fener)
    for i in erelat:
        print('{:12.4f}'.format(i,),file=fener)
    print("\n            Total loss = {:12.4f} %".format(loss_best),file=fener)
    autovarlambda.print_ener()
    os.system("mv relat_error.dat "+filename+".erp")
    os.system("mv tmp "+filename+".das")
    fener.close()
    return

#########################################################################################################
# main program starts
#########################################################################################################

# call input data
data_input()

# define new values for importan parameters
alfd_min=0.001
alfd_max=0.100
rcut_min=0.50
rcut_max=1.50

noise_var=None
exact_feval=True
ARD=True
varf=2.5
lf=0.1
xi=0.0001
optimize_restarts=5

initer=12
maxevals=48 
nener=10
npolvar=3
icost=0

print_input()
space=mypolspace(alfd_min,alfd_max,rcut_min,rcut_max,npolvar)
kernel=mykernel(rbf=1,stdperiodic=0)
run_mybo_pol(iseed,space,kernel,initer,maxevals,varf,lf,noise_var,exact_feval,optimize_restarts,ARD,xi)

