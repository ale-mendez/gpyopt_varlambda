#!/usr/bin/env python3

# Alejandra Mendez - 22/04/2019
# v1.1
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
#           itype    -- algorithm type
#                          = 1 : GP 
#                          = 2 : ...
#           igrid    -- type of space grid
#                          = 0 : continuous
#                          = 1 : ...
#           imin     -- consider ground plus/or-only excited states
#                          = 0 : includ ground state
#                          = 1 : only excited states
#           ifun     -- type of loss funtional
#                          = 0 : use loss as sum of weighted relative errros
#                          = 1 : use loss as sum of weighted square relative errors
#           maxevals -- number of evaluations in minimization
#           eps_up   -- superior deviation in lambda values
#           eps_down -- inferior deviation in lambda values
#           nlamvar  -- number of lambda values to be varied
#           nener    -- number of energies to be considered (ground + excited)
#           igi      -- type of weight to be used in minimizing functional
#                          = 0 : use weight input below
#                          = 1 : use statistical weight gi = sum 2j+1
#           weight   -- weight in relative errors (number of elements == nener)
#           cfgs     -- configurations to be included in the calculation
#
#
#  * input files:
#           varlam_hyperopt.yml     -- input file for .py
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
from hyperopt import hp, rand, tpe, Trials, fmin
import pandas as pd
import yaml
import GPy
import GPyOpt
from numpy.random import seed

# Default parameters
itype = 1
igrid = 0
imin = 0
ifun = 0
igi = 1

# Input parameters
with open("varlam_gpyopt.yml", 'r') as stream:
    data_loaded = yaml.load(stream)

itype = data_loaded.get("itype")
if itype==1:
    print(" ==> algorithm type: TPE")
if itype==2:
    print(" ==> algorithm type: random")

igrid = data_loaded.get("igrid") 
if igrid==0:
    print(" ==> grid type: uniform") 
if igrid==1:
    print(" ==> grid type: normal")

imin = data_loaded.get("imin") 
if imin==0:
    print(" ==> states included in loss function: ground + excited")
if imin==1:
    print(" ==> states included in loss function: excited")

ifun = data_loaded.get("ifun")
if ifun==0:
    print(" ==> loss functional: sum of weigthed relative errors")
if ifun==1:
    print(" ==> loss functional: sum of weigthed square relative errors")

maxevals = data_loaded.get("maxevals")
eps_up = data_loaded.get("eps_up")
eps_down = data_loaded.get("eps_down")
eps=eps_down/2.0
if eps_up > eps_down:
    eps=eps_up/2.0
nlamvar = data_loaded.get("nlamvar")
nener = data_loaded.get("nener")
igi = data_loaded.get("igi") 
if igi==0:
    print("using input weight")
    listweight = data_loaded.get("weight")
    NEMX = len(listweight)
    weight = np.array(listweight)
if igi==1:
    print("using statistical weight gi")
    NEMX = nener
    weight = np.full(nener,1.)
cfgs = data_loaded.get("cfgs")

################################################################################

# Define useful functions
def atom_name(nzion):
    nzion=abs(nzion)
    atomname=['0','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg']
    for i in range(len(atomname)):
        if nzion==i:
            chatom = atomname[i]
    return chatom
def pot_type(nzion):
    if nzion < 0:
        chtype="STO"
    else:
        chtype="TFDA"
    return chtype

# Define function for writing iteration number in a human-readable-way
def human_format(num):
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000
    return '%.f%s' % (num, ['', 'K', 'M', 'G', 'T', 'P'][magnitude])

# Define functions to be used
def error_relat(valexact,valcomp):
    if valexact == 0.0 or valcomp == 0.0:
        error=999
    error = abs((valexact-valcomp)/valexact)
    return error*100

# Loss function defined as weighted sum of energies relative errors
# def loss_wsumE(enere,enerc,neex,weight):
#     loss = 0.0
#     for i in range(neex):
#         diff = error_relat(enere[i],enerc[i])
#         loss = loss + weight[i]*diff
#     return loss
def loss_wsumE(enere,enerc,neex,weight):
    loss = 0.0
    if igi==1: # use statistical weight 
        weight=autovarlambda.compbck.gic
    for i in range(imin,neex): # consider ground energy and/or-only excited states
        diff = error_relat(enere[i],enerc[i])
        if ifun==1: # new functional 
            diff = diff*diff
        loss = loss + weight[i]*diff
    return loss

# Define function to be minimized
def var_lambda(x):
    minlam=0.1
    lam0 = x[0,:]
    nlam0 = len(lam0)
    for i in range(nlam0):
        if lam0[i] < minlam:
            lam0[i] = minlam
    autovarlambda.run_as(lam0,nlam0)
    neex = autovarlambda.exactbck.ne
    enere = autovarlambda.eneroutbck.enere
    enerc = autovarlambda.eneroutbck.enerc
    loss = loss_wsumE(enere,enerc,neex,weight)
    print(lam0,loss)
    return loss

# Determine file name string vector from fortran subroutine data
def def_filename():
    ncfg=autovarlambda.salgebbck.mxconf
    nzion=autovarlambda.salgebbck.nzion
    chatom=atom_name(nzion)
    chtype=pot_type(nzion)
    cfgs=str(ncfg)+"CFG"+chtype
    filename=chatom+cfgs+"_GP"+human_format(maxevals)
    return ncfg, filename

# Open .out file
def open_fout(filename,ncfg):
    fener = open(filename+".out","w")
    return fener

# Define and initialize variables
def init_var(NVMX,NEMX):
    lam0 = np.array(NVMX)
    nlam0 = NVMX
    enere = np.zeros(NEMX)
    enerc = np.zeros(NEMX)
    neex = int()
    return lam0, nlam0, enere, enerc, neex

def print_eresults():
    lam_best = myBopt.x_opt
    nlam0 = len(lam_best)
    autovarlambda.run_as(lam_best.reshape(-1,nlam0),nlam0)
    neex = autovarlambda.exactbck.ne
    enere = autovarlambda.eneroutbck.enere
    enerc = autovarlambda.eneroutbck.enerc
    loss_best = loss_wsumE(enere,enerc,neex,weight)
#    print("x*= {} f(x*)= {}".format(lam_best,loss_best))
    erroregr = error_relat(enere[0],enerc[0])
    print("\n===>>> print results in .out <<<===")
    print("\n Number of evaluations = {:10d}".format(maxevals),file=fener)
    print("                  Time = {:10.3f} minutes".format(total),file=fener)
    print("-"*80,file=fener)
    print(" Best results:",file=fener)
    print("-"*80,file=fener)
    print("\n  lambda = {} ".format(lam_best),file=fener)
    print("\n  Ground State Energy  = {:12.6f}".format(enerc[0]),
          "\n                  NIST = {:12.6f}".format(enere[0]),
          "\n                   Er% = {:12.4f} %".format(erroregr),file=fener)
    for i in range(1,neex):
        errorex = error_relat(enere[i],enerc[i])
        print("\n       {} Excit. Energy = {:12.6f}".format(i,enerc[i]),
              "\n                  NIST = {:12.6f}".format(enere[i]),
              "\n                   Er% = {:12.6f} %".format(errorex),file=fener)
    print("\n            Total loss = {:12.4f} %".format(loss_best),file=fener)
    fener.close()

################################################################################

# Initialize variables
lv = np.full(nlamvar, 1.)
NVMX=len(lv)
ntcfg = len(cfgs)
lam0, nlam0, enere, enerc, neex = init_var(NVMX,NEMX)

# Define the search space domain
bounds=[]
for i in range(NVMX):
    alam=str(i+1)
    if igrid==0:
        bounds.append({'name': 'lam'+alam, 'type': 'continuous', 'domain': (lv[i]-eps_down,lv[i]+eps_up)})



# Begin loop in "configurations"
for i in range(ntcfg):
    icfg = cfgs[i]
    print("\n===>>> RUN {:3d} CFGs <<<===\n".format(icfg))
    cpdas = "cp das_" + str(icfg) + "CFG das"
    os.system(cpdas)
# Run initialization f2py subroutines
    autovarlambda.open_files()
    autovarlambda.inp_das()
    autovarlambda.inp_exact()

# if "ne" of exactvalues.dat is != than nener of .yaml
#  => "nener" is considered and ne is dumped
    dummyne = autovarlambda.exactbck.ne.copy()
    if nener != dummyne:
        print(" ************ forcing: neex == nener ************\n")
        autovarlambda.exactbck.ne = nener

# Initialize files and necessary arrays
    ncfg, filename = def_filename()
    if ncfg != icfg:
        print(" >>>> configuration mismatch !!!! ")
        sys.exit()
    fener = open_fout(filename,ncfg)

# Creates three identical objects that we will later use to compare the optimization strategies 
    print("===>>> RUN GPyOpt <<<===")
    t0 = time.time()
    myBopt = GPyOpt.methods.BayesianOptimization(f=var_lambda,
                                                 domain=bounds,
                                                 model_type = 'GP',
                                                 acquisition_type='EI',  
                                                 normalize_Y = True,
                                                 acquisition_weight = 2)
# runs the optimization for the three methods
    myBopt.run_optimization(maxevals,verbosity=False)
    myBopt.save_evaluations(filename+".dat")
    myBopt.save_models(filename+".mod")
    myBopt.save_report(filename+".rep")
    t1 = time.time()
    total = (t1-t0)/60.

# run autostructure for best lambda found
    print_eresults()

# cat .out file
    with open(filename+".out","r") as fp:
       line = fp.readline()
       cnt = 1
       while line:
           print(line,end='')
           line = fp.readline()
    print('')

os.system("rm CONFIG.DAT LEVELS TERMS olg oic ols tmp das")
