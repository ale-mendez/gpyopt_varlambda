#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as mticker
import os

plt.rc('text', usetex=True)
plt.rc('font', family='sans')
plt.rcParams["figure.figsize"] = [8,7]

colors=['k','r','g','b','m',
        'y','c','tab:blue','tab:orange','tab:green',
        'tab:red','tab:purple','tab:brown','tab:pink','tab:gray',
        'tab:olive','tab:cyan','darkgoldenrod','coral','forestgreen']
orbs=['1s','2s','2p','3s','3p',
      '3d','4s','4p','4d','4f',
      '5s','5p','5d','5f','5g',
      '6s','6p','6d','6f','6g']
symbs=['o','s','v','^','>',
       'o','s','v','^','>',
       'o','s','v','^','>',
       'o','s','v','^','>']

sp_left=0.16
sp_right=0.9
sp_bottom=0.15
sp_top=0.98
xpad=15
ypad=10

#################################################################################
#                              F U N C T I O N S
#################################################################################
def load_Jmindata(folder):
    global ncfgs,initer,maxevals,min_lam,max_lam,filename
    filename="O"+str(ncfgs)+"CFGTFDA_"+str(initer)+"GP"+str(maxevals)
    df=pd.read_csv(folder+"/seed_1/"+filename+".dat",sep='\t')
    Jmin=pd.DataFrame(columns=df.columns)
    for i in range(10):
        fseed="/seed_"+str(i+1)+"/"
        df=pd.read_csv(folder+fseed+filename+".dat",sep='\t')
        Jmin=Jmin.append(df.iloc[df['Y'].idxmin(), :], ignore_index=True)
    Jmin,nlam=rename_Jmin(Jmin)
    min_lam,max_lam=find_xtrms(Jmin,nlam)
    return Jmin,nlam

def load_erpdata(folder):
    global filename,ntran,nener
    nskip=list(range(1,nener+1))
    df=pd.read_csv(folder+"/seed_1/"+filename+".erp",sep='\s+')
    erp=[]
    for i in range(ntran):
        dumerp=pd.DataFrame(columns=df.columns)
        erp.append(dumerp)
    for i in range(10):
        fseed="/seed_"+str(i+1)+"/"
        df=pd.read_csv(folder+fseed+filename+".erp",skiprows=nskip,sep='\s+')
        for j in range(ntran-1):
            erp[j]=erp[j].append(df.iloc[j, :],ignore_index=True)
    return erp

def rename_Jmin(Jmin):
    Jmin=Jmin.rename(columns={"Y": "erp"});
    nlam=len(Jmin.columns)-2
    for i in range(nlam):
        Jmin=Jmin.rename(columns={"var_"+str(i+1): orbs[i]});
    return Jmin,nlam

def find_xtrms(Jmin,nlam):
    lamxtrm=[]
    for i in range(nlam):
        iorb=orbs[i]
        lamxtrm.append(min(Jmin[iorb]))
        lamxtrm.append(max(Jmin[iorb]))
    max_lam=max(lamxtrm)
    min_lam=min(lamxtrm)
    max_lam=max_lam+max_lam*0.05
    min_lam=min_lam-min_lam*0.05
    return min_lam,max_lam

def print_resume(folder,Jmin,nlam):
    global erp_min,erp_max,filename
    rname='resumen_latin.dat'
    fresume=open(rname,'w')
    Jminmax=max(Jmin['erp'])
    Jminmin=min(Jmin['erp'])
    erp_min=Jminmin
    erp_max=Jminmax
    lJmin=Jmin['erp'].tolist()
    iminmin=lJmin.index(Jminmin)
    lam_min=Jmin.iloc[iminmin].tolist()
    xlambda=' lambda = '
    for i in range(nlam):
        xlambda=xlambda+' {vp['+str(i)+']:.{prec}f} '
    print("\n *** "+folder+" ***",file=fresume)
    print(" max(Ji) =",Jminmax,", min(Ji) =",Jminmin,file=fresume)
    print("\n Best results:\n",13*'=',file=fresume)
    print("      J = ",Jminmin,file=fresume)
    print(xlambda.format(vp=lam_min[2:],prec=4),file=fresume)
    print("   seed = ",iminmin+1,"\n",file=fresume)
    fprint=folder+"/seed_"+str(iminmin+1)+"/"+filename
    f = open(fprint+".erp", "r")
    text = f.read()
    print(text,file=fresume)
    f.close()
    os.system("enscript "+fprint+".erp -o - | ps2pdf - erpmin_"+folder+".pdf ")
    fresume.close()
    os.system("cat "+rname)
    return

def plot_Jmin(folder,Jmin):
    global erp_min,erp_max,nbins,Nmin,Nmax
    fig=plt.figure()
    ax=plt.subplot(111)
    bins=np.linspace(erp_min,erp_max,nbins+1)
    plt.hist(Jmin['erp'],bins=bins,color='orange')
    ax.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
    ax.tick_params(direction='in',which='minor',length=4,bottom=True,top=True,left=True,right=True);
    plt.xlabel("Costo",fontsize=35,labelpad=xpad)
    plt.ylabel("N",fontsize=35,labelpad=ypad)
#    plt.xlim(erp_min*0.9,erp_max*1.1)
    plt.ylim(Nmin,Nmax)
    #ax.xaxis.set_major_locator(mticker.MultipleLocator(xerp_maj));
    #ax.yaxis.set_major_locator(mticker.MultipleLocator(20));
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator(n=5));
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator(n=5))
    plt.subplots_adjust(left=sp_left,right=sp_right,bottom=sp_bottom,top=sp_top)   # <--
    plt.savefig("Jmin_"+folder+".eps",transparent=True)
    #plt.show()
    return

def plot_imin(folder,Jmin):
    global erp_min,erp_max,imin,imax,initer,totevals,yimin_maj,yimin_min,xerp_maj,xerp_min
    fig=plt.figure()
    ax=plt.subplot(111)
    plt.plot(Jmin['erp'],Jmin['Iteration'],'bo',markersize=12,alpha=0.75)
    ax.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
    ax.tick_params(direction='in',which='minor',length=4,bottom=True,top=True,left=True,right=True);
    if imax>100: ipad=5
    else: ipad=ypad
    plt.xlabel("Costo",fontsize=35,labelpad=xpad)
    plt.ylabel("Iteración",fontsize=35,labelpad=ipad)
    plt.hlines(initer,erp_min*0.9,erp_max*1.1,colors='k',linestyles='dashed')
    plt.hlines(totevals,erp_min*0.9,erp_max*1.1,colors='k',linestyles='dashed')
#    plt.xlim(erp_min*0.9,erp_max*1.1)
    #plt.ylim(imin,imax)
    #ax.xaxis.set_major_locator(mticker.MultipleLocator(xerp_maj));
    #ax.yaxis.set_major_locator(mticker.MultipleLocator(yimin_maj));
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator(n=5));
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator(n=5));
    plt.subplots_adjust(left=sp_left,right=sp_right,bottom=sp_bottom,top=sp_top)   # <--
    plt.savefig("imin_"+folder+".eps",transparent=True)
    #plt.show()
    return

def plot_minspace(folder,Jmin,nlam):
    global min_lam,max_lam

    if nlam<=3:
        fig=plt.figure()
        ax=plt.subplot(111)
        for i in range(nlam):
            iorb=orbs[i]
            plt.plot(Jmin.index,Jmin[iorb],c=colors[i],marker=symbs[i],linestyle='None',markersize=12,alpha=0.75,label=iorb)
        ax.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
        ax.tick_params(direction='in',which='minor',length=4,bottom=True,top=True,left=True,right=True);
        ax.legend(bbox_to_anchor=(1,1),loc='upper left',fontsize=30,
                  ncol=1,handlelength=0.2,labelspacing=0.1,borderpad=0.2,
                  handletextpad=0.4,frameon=False)
        ax.set_ylabel(r"$\lambda_{nl}$",fontsize=35,labelpad=ypad)
        ax.set_xlabel('Iteración',fontsize=35,labelpad=xpad)
        plt.subplots_adjust(left=0.15,right=0.88,bottom=sp_bottom,top=sp_top)   # <--
        plt.savefig("minspace_"+folder+".eps",transparent=True)

    if nlam>3 and nlam<=6:
        nrow=1; ncol=2; nplot=[0,3,nlam]
        fig, axs = plt.subplots(nrows=nrow, ncols=ncol,sharey=True)
        ii=0
        scatlam=[]
        for ax in axs.reshape(-1): 
            for i in range(nplot[ii],nplot[ii+1]):
                iorb=orbs[i]
                lami,=ax.plot(Jmin.index,Jmin[iorb],c=colors[i],marker=symbs[i],linestyle='None',markersize=12,alpha=0.75,label=iorb)
                scatlam.append(lami)
            ax.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
            ax.tick_params(direction='in',which='minor',length=4,bottom=True,top=True,left=True,right=True);
#             ax.yaxis.set_major_locator(mticker.MultipleLocator(0.2));
#             ax.xaxis.set_major_locator(mticker.MultipleLocator(0.05));
            ax.yaxis.set_minor_locator(mticker.AutoMinorLocator(n=5));
            ax.xaxis.set_minor_locator(mticker.AutoMinorLocator(n=5));
            ii=ii+1
            ax.set_ylabel(r"$\lambda_{nl}$",fontsize=35,labelpad=ypad)
        ax.set_xlabel('Iteración',fontsize=35,labelpad=xpad,position=(0,1))
        # Hide x labels and tick labels for top plots and y ticks for right plots.
        for ax in axs.flat:
            ax.label_outer()
        plt.legend(handles=scatlam,bbox_to_anchor=(1,1),loc='upper left',fontsize=30,
                   ncol=1,handlelength=0.2,labelspacing=0.1,borderpad=0.2,
                   handletextpad=0.4,frameon=False)
        for ax in fig.get_axes():
            ax.label_outer()
        fig.subplots_adjust(hspace=0.05, wspace=0.05)
        plt.subplots_adjust(left=0.15,right=0.88,bottom=sp_bottom,top=sp_top)   # <--
        plt.savefig("minspace_"+folder+".eps",transparent=True)
    
    if nlam>6 and nlam<=20:
        nrow=1; ncol=2; nplot=[0,3,6]
        fig, axs = plt.subplots(nrows=nrow, ncols=ncol,sharey=True)
        ii=0
        scatlam=[]
        for ax in axs.reshape(-1): 
            for i in range(nplot[ii],nplot[ii+1]):
                iorb=orbs[i]
                lami,=ax.plot(Jmin.index,Jmin[iorb],c=colors[i],marker=symbs[i],linestyle='None',markersize=12,alpha=0.75,label=iorb)
                scatlam.append(lami)
            ax.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
            ax.tick_params(direction='in',which='minor',length=4,bottom=True,top=True,left=True,right=True);
#             ax.yaxis.set_major_locator(mticker.MultipleLocator(0.2));
#             ax.xaxis.set_major_locator(mticker.MultipleLocator(0.05));
            ax.yaxis.set_minor_locator(mticker.AutoMinorLocator(n=5));
            ax.xaxis.set_minor_locator(mticker.AutoMinorLocator(n=5));
            ii=ii+1
            ax.set_ylabel(r"$\lambda_{nl}$",fontsize=35,labelpad=ypad)
        ax.set_xlabel('Iteración',fontsize=35,labelpad=xpad,position=(0,1))
        # Hide x labels and tick labels for top plots and y ticks for right plots.
        for ax in axs.flat:
            ax.label_outer()
        plt.legend(handles=scatlam,bbox_to_anchor=(1,1),loc='upper left',fontsize=30,
                   ncol=1,handlelength=0.2,labelspacing=0.1,borderpad=0.2,
                   handletextpad=0.4,frameon=False)
        for ax in fig.get_axes():
            ax.label_outer()
        fig.subplots_adjust(hspace=0.05, wspace=0.05)

        plt.subplots_adjust(left=0.15,right=0.88,bottom=sp_bottom,top=sp_top)   # <--
        plt.savefig("minspace1_"+folder+".eps",transparent=True)
        
        nrow=1; ncol=3; nplot=[7,11,15,20]
        fig, axs = plt.subplots(nrows=nrow, ncols=ncol,sharey=True)
        ii=0
        scatlam=[]
        for ax in axs.reshape(-1): 
            for i in range(nplot[ii],nplot[ii+1]):
                iorb=orbs[i]
                lami,=ax.plot(Jmin.index,Jmin[iorb],c=colors[i],marker=symbs[i],linestyle='None',markersize=12,alpha=0.75,label=iorb)
                scatlam.append(lami)
            ax.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
            ax.tick_params(direction='in',which='minor',length=4,bottom=True,top=True,left=True,right=True);
#             ax.yaxis.set_major_locator(mticker.MultipleLocator(0.2));
#             ax.xaxis.set_major_locator(mticker.MultipleLocator(0.05));
            ax.yaxis.set_minor_locator(mticker.AutoMinorLocator(n=5));
            ax.xaxis.set_minor_locator(mticker.AutoMinorLocator(n=5));
            ii=ii+1
            ax.set_ylabel(r"$\lambda_{nl}$",fontsize=35,labelpad=ypad)
        ax.set_xlabel('Iteración',fontsize=35,labelpad=xpad,position=(0,1))
        # Hide x labels and tick labels for top plots and y ticks for right plots.
        for ax in axs.flat:
            ax.label_outer()
        plt.legend(handles=scatlam,bbox_to_anchor=(1,1),loc='upper left',fontsize=30,
                   ncol=1,handlelength=0.2,labelspacing=0.1,borderpad=0.2,
                   handletextpad=0.4,frameon=False)
        for ax in fig.get_axes():
            ax.label_outer()
        fig.subplots_adjust(hspace=0.05, wspace=0.05)

        plt.subplots_adjust(left=0.15,right=0.88,bottom=sp_bottom,top=sp_top)   # <--
        plt.savefig("minspace2_"+folder+".eps",transparent=True)
    #plt.show()
    return

def plot_akierp(erp):
    global folder,ypad,xpad,akierp_min1,akierp_max1,akierp_min2,akierp_max2
    fig=plt.figure(figsize=(12,7))
    ax1=plt.subplot(121)
    tran1,=plt.plot(erp[0]['er%'],erp[0].index,c=colors[0],marker=symbs[0],linestyle='None',markersize=12,label="1")
    tran2,=plt.plot(erp[1]['er%'],erp[1].index,c=colors[1],marker=symbs[1],linestyle='None',markersize=12,label="2")
    tran3,=plt.plot(erp[2]['er%'],erp[2].index,c=colors[2],marker=symbs[2],linestyle='None',markersize=12,label="3")
    plt.ylabel("N",fontsize=35,labelpad=ypad)
    plt.xlabel(r"Error relativo \%",fontsize=35,labelpad=xpad,position=(1,0))
#    plt.xlim(akierp_min1,akierp_max1)
    ax2=plt.subplot(122)
    tran4,=plt.plot(erp[3]['er%'],erp[3].index,c=colors[3],marker=symbs[3],linestyle='None',markersize=12,label="4")
    tran5,=plt.plot(erp[4]['er%'],erp[4].index,c=colors[4],marker=symbs[4],linestyle='None',markersize=12,label="5")
#    plt.xlim(akierp_min2,akierp_max2)
    ax1.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
    ax1.tick_params(direction='in',which='minor',length=4,bottom=True,top=True,left=True,right=True);
    ax2.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
    ax2.tick_params(direction='in',which='minor',length=4,bottom=True,top=True,left=True,right=True);
    ax1.get_shared_y_axes().join(ax1, ax2)
    for ax in fig.get_axes():
        ax.label_outer()
    fig.subplots_adjust(hspace=0.05, wspace=0.05)
    plt.legend(handles=[tran1,tran2,tran3,tran4,tran5],bbox_to_anchor=(1,1),loc='upper left',fontsize=25,ncol=1,handlelength=0.2,frameon=False)
    ax1.xaxis.set_minor_locator(mticker.AutoMinorLocator(n=5))
    ax2.xaxis.set_minor_locator(mticker.AutoMinorLocator(n=5))
    plt.subplots_adjust(left=0.11,right=0.75,bottom=sp_bottom,top=sp_top)   # <--
    plt.savefig("akierp_"+folder+".eps",transparent=True)
#    plt.show()
    return

#################################################################################
#                             M A I N   P R O G R A M
#################################################################################

global nbins,Nmin,Nmax,lam_min,lam_max,erp_min,erp_max,imin,imax,yimin_maj,yimin_min,xerp_maj,xerp_min,akierp_min1,akierp_max1,akierp_min2,akierp_max2
global initer,totevals,ncfgs,filename,ntran

ncfgs=15
initer=60
maxevals=24
ntran=5
nener=8
initer=int(input('initer:  '))
maxevals=int(input('maxevals:  '))
#ntran=int(input('ntran:  '))
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

nbins=18
erp_min=0
erp_max=190
Nmin=0
Nmax=100

akierp_min1=0
akierp_max1=450
akierp_min2=0
akierp_max2=6000

yimin_maj=int(imax/4.)
yimin_min=yimin_maj/5.

xerp_maj=0.5
xerp_min=0.1

#################################################################################
#                                  L A T I N
#################################################################################

if slat=='y':
    folder="latin"
    Jmin_lat,nlam_lat=load_Jmindata(folder)
    erp_lat=load_erpdata(folder)
    print_resume(folder,Jmin_lat,nlam_lat)
    plot_Jmin(folder,Jmin_lat)
    plot_imin(folder,Jmin_lat)
    plot_minspace(folder,Jmin_lat,nlam_lat)
    plot_akierp(erp_lat)
