#!/usr/bin/env python3

# Alejandra Mendez - 28/04/2019
# v1.2
#
#  * this program implements hyperoptimization of NLAMVAR parameters in
#     autosctructure code (max=20)
#  * the autovarlambda.f90 subroutine is compiled with f2py and incorporated
#     to this code
#  * the autostructure code is executed using system
#  * multiple configurations can be computed:
#     "cfgs" defined in .yaml are considered in the computation
#  * multiple excited energies can be considered via "nener"
#  * the number of elements in array "weight" must be equal
#     to the amount of energies "nener" considered
#
#  * input variables:
#           cfgs     -- configurations to be included in the calculation
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
# os.system("f2py -c -m autovarlambda f2py_autovarlambda.f90 --quiet")
# print("\n===>>> autovarlambda compiled with f2py <<<===")

# Import other stuff
import numpy as np
import matplotlib.pyplot as plt
import autovarlambda
import time
import pandas as pd
import yaml
import GPy
import GPyOpt
from numpy.random import seed

# Input parameters
print("\n===>>> Input parameters <<<===\n")
with open("varlam_gpyopt.yml", 'r') as stream:
    data_loaded=yaml.full_load(stream)

cfgs=data_loaded.get("cfgs")
ntcfg=len(cfgs)
if cfgs is None: raise ValueError("Error: cfgs not defined.")

maxevals=data_loaded.get("maxevals")
if maxevals is None: raise ValueError("Error: maxevals not defined.")

mtype=data_loaded.get("mtype")
if mtype=="GP": print(" * model type: "+mtype)
elif mtype is None: 
    mtype="GP"
    print("Warning: model type not defined. Use default: GP")
else: raise ValueError("Error: "+mtype+" not implemented.")

aftype=data_loaded.get("aftype")
if aftype=="EI" or aftype=="MPI" or aftype=="GP-UCB": print(" * acquisition function type: "+aftype)
elif aftype is None: 
    aftype="EI"
    print("Warning: acquisition function type not defined. Use default: EI")
else: raise ValueError("Error: "+mtype+" not implemented.")

afweight=data_loaded.get("afweight")
if afweight: print(" * acquisition function weight: "+str(afweight))
elif afweight is None: 
    afweight=1.
    print("Warning: acquisition function weight not defined. Use default: 1.")

gridtype=data_loaded.get("gridtype")
if gridtype=="continuous":
    igrid=0
    print(" * grid type: "+gridtype)
elif gridtype is None: 
    igrid=0
    print("Warning: grid type not defined. Use default: continuous")
else: raise ValueError("Error: "+gridtype+" not implemented.")

nlamvar=data_loaded.get("nlamvar")
if nlamvar is None: raise ValueError("Error: nlamvar not defined.")

maxlam=data_loaded.get("maxlam")
minlam=data_loaded.get("minlam")
if maxlam is None or minlam is None: raise ValueError("Error: maxlam and/or minlam not defined.")

mfunc=data_loaded.get("mfunc")
if mfunc=="Er":
    ifun=0
    print(" * minimization function: sum of weigthed relative errors")
elif mfunc=="Er**2":
    ifun=1
    print(" * minimization function: sum of weigthed square relative errors")
elif mfunc is None: 
    imin=0
    print("Warning: minimization function not defined. Use default: Er ")
else: raise ValueError("Error: "+mfunc+" not implemented.")

minst=data_loaded.get("minst")
if minst=="gr+ex":
    imin=0
    print(" * states included in minimization: ground + excited")
elif minst=="ex":
    imin=1
    print(" * states included in minimization: excited")
elif minst is None: 
    ifun=0
    print("Warning: states to be minimized not defined. Use default: ground + excited")
else: raise ValueError("Error: "+minst+" not implemented.")

nener=data_loaded.get("nener")
if nener is None: raise ValueError("Error: nener not defined.")

wi=data_loaded.get("wi")
if wi=="eq" or wi=="gi":
    if wi=="eq": igi=0
    if wi=="gi": igi=1
    weight=np.full(nener,1.)
    print(" * weight type: "+wi)
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

################################################################################

# Define useful functions
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

# Define functions to be used
def error_relat(valexact,valcomp):
    if valexact == 0.0 or valcomp == 0.0:
        error=999
    error=abs((valexact-valcomp)/valexact)*100
    return error

# Loss function defined as weighted sum of energies relative errors
# def loss_wsumE(enere,enerc,neex,weight):
#     loss=0.0
#     for i in range(neex):
#         diff=error_relat(enere[i],enerc[i])
#         loss=loss + weight[i]*diff
#     return loss
def loss_wsumE(enere,enerc,neex,weight):
    loss=0.0
    if igi==1: weight=autovarlambda.compbck.gic
    for i in range(imin,neex): # consider ground energy and/or-only excited states
        diff=error_relat(enere[i],enerc[i])
        if ifun==1: diff=diff*diff
        loss=loss + weight[i]*diff
    return loss

# Define function to be minimized
def var_lambda(x):
    minlam=0.1
    lam=x[0,:]
    autovarlambda.run_as(lam,nlamvar)
    neex=autovarlambda.exactbck.ne
    enere=autovarlambda.eneroutbck.enere
    enerc=autovarlambda.eneroutbck.enerc
    loss=loss_wsumE(enere,enerc,neex,weight)
    return loss

# Determine file name string vector from fortran subroutine data
def def_filename():
    ncfg=autovarlambda.salgebbck.mxconf
    nzion=autovarlambda.salgebbck.nzion
    chatom=atom_name(nzion)
    chtype=pot_type(nzion)
    cfgs=str(ncfg)+"CFG"+chtype
    filename=chatom+cfgs+"_"+mtype+human_format(maxevals)
    return ncfg, filename

# Open .out file
def open_fout(filename,ncfg):
    fener=open(filename+".out","w")
    return fener

# Define and initialize variables
def init_var(nener):
    enere=np.zeros(nener)
    enerc=np.zeros(nener)
    neex=int()
    return enere, enerc, neex

# Print best results
def print_eresults(lam_best,nlamvar):
    autovarlambda.run_as(lam_best,nlamvar)
    neex=autovarlambda.eneroutbck.neex
    enere=autovarlambda.eneroutbck.enere
    enerc=autovarlambda.eneroutbck.enerc
    loss_best=loss_wsumE(enere,enerc,neex,weight)
    erroregr=error_relat(enere[0],enerc[0])
    print("\n Number of evaluations={:10d}".format(maxevals),file=fener)
    print("                  Time={:10.3f} minutes".format(total),file=fener)
    print("-"*80,file=fener)
    print(" Best results:",file=fener)
    print("-"*80,file=fener)
    print("\n  lambda={} ".format(lam_best[0]),file=fener)
    print("\n  Ground State Energy ={:12.6f}".format(enerc[0]),
          "\n                  NIST={:12.6f}".format(enere[0]),
          "\n                   Er%={:12.4f} %".format(erroregr),file=fener)
    for i in range(1,neex):
        errorex=error_relat(enere[i],enerc[i])
        print("\n       {} Excit. Energy={:12.6f}".format(i,enerc[i]),
              "\n                  NIST={:12.6f}".format(enere[i]),
              "\n                   Er%={:12.6f} %".format(errorex),file=fener)
    print("\n            Total loss={:12.4f} %".format(loss_best),file=fener)
    autovarlambda.print_ener()
    os.system("mv relat_error.dat "+filename+".erp")
    fener.close()

################################################################################

# Initialize variables
enere, enerc, neex=init_var(nener)

# Define the search space domain
space=[]
for i in range(nlamvar):
    alam=str(i+1)
    if igrid==0:
        space.append({'name': 'lam'+alam, 'type': gridtype, 'domain': (minlam,maxlam)})

# Begin loop in "configurations"
for i in range(ntcfg):
    icfg=cfgs[i]
    print("\n===>>> RUN {:3d} CFGs <<<===\n".format(icfg))
    cpdas="cp das_" + str(icfg) + "CFG das"
    os.system(cpdas)
# Run initialization f2py subroutines
    autovarlambda.open_files()
    autovarlambda.inp_das()
    autovarlambda.inp_exact()

# if "ne" of exactvalues.dat is != than nener of .yaml
#  => "nener" is considered and ne is dumped
    dummyne=autovarlambda.exactbck.ne.copy()
    if nener != dummyne:
        print(" Warning: forcing neex == nener \n")
        autovarlambda.exactbck.ne=nener

# Initialize files and necessary arrays
    ncfg, filename=def_filename()
    if ncfg != icfg:
        print(" >>>> configuration mismatch !!!! ")
        sys.exit()
    fener=open_fout(filename,ncfg)

# Create object with BO method
    print("===>>> RUN GPyOpt <<<===")
    t0=time.time()
    myBopt=GPyOpt.methods.BayesianOptimization(f=var_lambda,
                                               domain=space,
                                               model_type=mtype,
                                               acquisition_type=aftype,
                                               normalize_Y=True,
                                               acquisition_weight=afweight)
    myBopt.run_optimization(maxevals,verbosity=False)
    t1=time.time()
    total=(t1-t0)/60.


# PRINT RESULTS
    lam_best=myBopt.x_opt
    lam_best=lam_best.reshape(-1,nlamvar)
    print("\n===>>> Print results <<<===\n")
    myBopt.save_evaluations(filename+".dat")
    myBopt.save_models(filename+".mod")
    myBopt.save_report(filename+".rep")
    print_eresults(lam_best,nlamvar)

# cat .out file
    icat=0
    if icat==1:
        with open(filename+".out","r") as fp:
            line=fp.readline()
            cnt=1
            while line:
                print(line,end='')
                line=fp.readline()
        print('')

os.system("rm CONFIG.DAT TERMS LEVELS olg ols oic das ")
os.system("mv tmp "+filename+".das")
