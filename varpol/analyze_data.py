#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as mticker

plt.rc('text', usetex=True)
plt.rc('font', family='sans')
global nbins,Nmin,Nmax,amin,amax,rmin,rmax,jmin,jmax,imin,imax,initer,totevals,yimin_maj,yimin_min,xerp_maj,xerp_min

initer=int(input('initer:  '))
maxevals=int(input('maxevals:  '))
#srand=str(input('random (y/n):  '))
#slat=str(input('latin (y/n):  '))
srand='n'
slat='y'

totevals=initer+maxevals
imin=0
imax=totevals*1.10
if maxevals==0: imax=(totevals+60)*1.10
if maxevals==12: imax=(totevals+48)*1.10
if maxevals==24: imax=(totevals+36)*1.10
if maxevals==36: imax=(totevals+24)*1.10
if maxevals==48: imax=(totevals+12)*1.10

nbins=39
jmin=0.15
jmax=4.05
Nmin=0
Nmax=55

yimin_maj=int(imax/5.)
yimin_min=yimin_maj/5.

xerp_maj=0.5
xerp_min=0.1

########################################################################

def rename_params(alpha,rho):
    alpha=alpha.rename(columns={2: "l=0"});
    alpha=alpha.rename(columns={4: "l=1"});
    alpha=alpha.rename(columns={6: "l=2"});
    rho=rho.rename(columns={2: "l=0"});
    rho=rho.rename(columns={4: "l=1"});
    rho=rho.rename(columns={6: "l=2"});
    return alpha,rho

def find_xtrms(alpha,rho):
    axtrm=[]
    rxtrm=[]
    for i in range(3):
        lstr="l="+str(i)
        axtrm.append(min(alpha[lstr]))
        axtrm.append(max(alpha[lstr]))
        rxtrm.append(min(rho[lstr]))
        rxtrm.append(max(rho[lstr]))
    amax=max(axtrm)
    rmax=max(rxtrm)
    amin=min(axtrm)
    rmin=min(rxtrm)
    amax=amax+amax*0.05
    amin=amin-amax*0.05
    rmax=rmax+rmax*0.05
    rmin=rmin-rmax*0.05
    return amin,amax,rmin,rmax

def plot_Jmin(folder,Jmin):
    global jmin,jmax,nbins,Nmin,Nmax
    #plt.hist(Jmin['erp'],bins=10)
    #plt.show()
    fig=plt.figure(figsize=(8,8))
    ax=plt.subplot(111)
    bins=np.linspace(jmin,jmax,nbins+1)
    plt.hist(Jmin['erp'],bins=bins,color='orange')
    ax.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
    ax.tick_params(direction='in',which='minor',length=6,bottom=True,top=True,left=True,right=True);
    plt.xlabel("Error relativo \%",fontsize=35,labelpad=15)
    plt.ylabel("N",fontsize=35,labelpad=5)
    plt.xlim(jmin+0.075,jmax-0.075)
    plt.ylim(Nmin,Nmax)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(20));
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(5));
    ax.xaxis.set_major_locator(mticker.MultipleLocator(xerp_maj));
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(xerp_min));
    plt.subplots_adjust(left=0.15,right=0.9,bottom=0.14,top=0.9)   # <--
    plt.savefig("Jmin_"+folder+".eps",transparent=True)
    #plt.show()
    return

def plot_minspace(folder,alpha,rho):
    global amin,amax,rmin,rmax
    fig=plt.figure(figsize=(8,8))
    ax=plt.subplot(111)
    plt.plot(alpha["l=0"],rho["l=0"],'ko',markersize=12,alpha=0.75,label='$l=0$')
    plt.plot(alpha["l=1"],rho["l=1"],'rs',markersize=12,alpha=0.75,label='$l=1$')
    plt.plot(alpha["l=2"],rho["l=2"],'g^',markersize=12,alpha=0.75,label='$l=2$')
    ax.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
    ax.tick_params(direction='in',which='minor',length=6,bottom=True,top=True,left=True,right=True);
    ax.legend(loc='upper left',bbox_to_anchor=(0.01, 1.16, 0, 0),fontsize=30,ncol=3,handlelength=0.2,frameon=False)
    plt.xlabel(r"$\alpha_l$",fontsize=35,labelpad=15)
    plt.ylabel(r"$\rho_l$",fontsize=35,labelpad=15)
    plt.xlim(amin,amax)
    plt.ylim(rmin,rmax)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(0.2));
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(0.1));
    ax.xaxis.set_major_locator(mticker.MultipleLocator(0.05));
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(0.01));
    plt.subplots_adjust(left=0.15,right=0.9,bottom=0.14,top=0.9)   # <--
    plt.savefig("minspace_"+folder+".eps",transparent=True)
    #plt.show()
    return

def plot_imin(folder,iJmin):
    global jmin,jmax,imin,imax,initer,totevals,yimin_maj,yimin_min,xerp_maj,xerp_min
    fig=plt.figure(figsize=(8,8))
    ax=plt.subplot(111)
    plt.plot(iJmin['Y'],iJmin['Iteration'],'bo',markersize=12,alpha=0.75)
    ax.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
    ax.tick_params(direction='in',which='minor',length=6,bottom=True,top=True,left=True,right=True);
    plt.xlabel("Error relativo \%",fontsize=35,labelpad=15)
    if imax>100: ipad=5
    else: ipad=15
    plt.ylabel("Iteración",fontsize=35,labelpad=ipad)
    plt.xlim(jmin+0.075,jmax-0.075)
    plt.ylim(imin,imax)
    plt.hlines(initer,jmin,jmax,colors='k',linestyles='dashed')
    plt.hlines(totevals,jmin,jmax,colors='k',linestyles='dashed')
    ax.yaxis.set_major_locator(mticker.MultipleLocator(yimin_maj));
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(yimin_min));
    ax.xaxis.set_major_locator(mticker.MultipleLocator(xerp_maj));
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(xerp_min));
    plt.subplots_adjust(left=0.15,right=0.9,bottom=0.14,top=0.9)   # <--
    plt.savefig("imin_"+folder+".eps",transparent=True)
    #plt.show()
    return

def print_resume(folder,Jmin,alpha,rho):
    print("\n *** "+folder+" ***")
    Jminmax=max(Jmin['erp'])
    Jminmin=min(Jmin['erp'])
    print("J max =",Jminmax)
    print("J min =",Jminmin)
    lJmin=Jmin['erp'].tolist()
    iminmin=lJmin.index(Jminmin)
    alpha_min=alpha.iloc[iminmin].tolist()
    rho_min=rho.iloc[iminmin].tolist()
    print("alpha =",alpha_min[1:],"\n  rho =",rho_min[1:])
    return

filename="Be6CFGSTO_"+str(initer)+"GP"+str(maxevals)+".dat"
evenlist=[]
oddlist=[]
for i in range(200): 
    if (i % 2 == 0): evenlist.append(i) 
    else: oddlist.append(i) 

########################################################################
#                             L A T I N
########################################################################

if slat=='y':
    folder="latin"
    Jmin_lat=pd.read_csv(folder+"/total_loss.dat",sep='\s+',header=None,usecols=[3])
    Jmin_lat=Jmin_lat.rename(columns={3: "erp"});
    alpha_lat=pd.read_csv(folder+"/params.dat",sep='\s+',header=None,skiprows=oddlist,usecols=[0,2,4,6])
    rho_lat=pd.read_csv(folder+"/params.dat",sep='\s+',header=None,skiprows=evenlist,usecols=[0,2,4,6])
    alpha_lat,rho_lat=rename_params(alpha_lat,rho_lat)
    amin,amax,rmin,rmax=find_xtrms(alpha_lat,rho_lat)
    df=pd.read_csv(folder+"/seed_1/"+filename,sep='\t')
    iJmin_lat=pd.DataFrame(columns=df.columns)
    for i in range(100):
        fseed="/seed_"+str(i+1)+"/"
        df=pd.read_csv(folder+fseed+filename,sep='\t')
        iJmin_lat=iJmin_lat.append(df.iloc[df['Y'].idxmin(), :], ignore_index=True)
    plot_Jmin(folder,Jmin_lat)
    plot_minspace(folder,alpha_lat,rho_lat)
    plot_imin(folder,iJmin_lat)
    print_resume(folder,Jmin_lat,alpha_lat,rho_lat)


########################################################################
#                             R A N D O M
########################################################################

if srand=='y':
    folder="random"
    Jmin_rand=pd.read_csv(folder+"/total_loss.dat",sep='\s+',header=None,usecols=[3])
    Jmin_rand=Jmin_rand.rename(columns={3: "erp"});
    alpha_rand=pd.read_csv(folder+"/params.dat",sep='\s+',header=None,skiprows=oddlist,usecols=[0,2,4,6])
    rho_rand=pd.read_csv(folder+"/params.dat",sep='\s+',header=None,skiprows=evenlist,usecols=[0,2,4,6])
    alpha_rand,rho_rand=rename_params(alpha_rand,rho_rand)
    amin,amax,rmin,rmax=find_xtrms(alpha_rand,rho_rand)
    df=pd.read_csv(folder+"/seed_1/"+filename,sep='\t')
    iJmin_rand=pd.DataFrame(columns=df.columns)
    for i in range(100):
        fseed="/seed_"+str(i+1)+"/"
        df=pd.read_csv(folder+fseed+filename,sep='\t')
        iJmin_rand=iJmin_rand.append(df.iloc[df['Y'].idxmin(), :], ignore_index=True)
    plot_Jmin(folder,Jmin_rand)
    plot_minspace(folder,alpha_rand,rho_rand)
    plot_imin(folder,iJmin_rand)
    print_resume(folder,Jmin_rand,alpha_rand,rho_rand)

