#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as mticker

plt.rc('text', usetex=True)
plt.rc('font', family='sans')

initer=36  
maxevals=36  
ymajor_imin=20
yminor_imin=5

nbins=7
jmin=0.15
jmax=0.85
Nmin=0
Nmax=110
totevals=initer+maxevals
imin=0
imax=totevals*1.10

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

########################################################################
#                             L A T I N
########################################################################

folder="latin"
total_loss_lat=pd.read_csv(folder+"/total_loss.dat",sep='\s+',header=None,usecols=[3])
total_loss_lat=total_loss_lat.rename(columns={3: "erp"});

evenlist=[]
oddlist=[]
for i in range(200): 
    if (i % 2 == 0): evenlist.append(i) 
    else: oddlist.append(i) 
alpha_lat=pd.read_csv(folder+"/params.dat",sep='\s+',header=None,skiprows=oddlist,usecols=[0,2,4,6])
rho_lat=pd.read_csv(folder+"/params.dat",sep='\s+',header=None,skiprows=evenlist,usecols=[0,2,4,6])
alpha_lat,rho_lat=rename_params(alpha_lat,rho_lat)
amin,amax,rmin,rmax=find_xtrms(alpha_lat,rho_lat)

filename="Be6CFGSTO_"+str(initer)+"GP"+str(maxevals)+".dat"
df=pd.read_csv(folder+"/seed_1/"+filename,sep='\t')
df2=pd.DataFrame(columns=df.columns)
for i in range(100):
    fseed="/seed_"+str(i+1)+"/"
    df=pd.read_csv(folder+fseed+filename,sep='\t')
    df2=df2.append(df.iloc[df['Y'].idxmin(), :], ignore_index=True)

print(folder)
print("max=",max(total_loss_lat['erp']))
print("min=",min(total_loss_lat['erp']))

########################################################################
#plt.hist(total_loss_lat['erp'],bins=10)
#plt.show()
fig=plt.figure(figsize=(8,8))
ax=plt.subplot(111)
bins=np.linspace(jmin,jmax,nbins+1)
plt.hist(total_loss_lat['erp'],bins=bins,color='orange')
ax.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
ax.tick_params(direction='in',which='minor',length=6,bottom=True,top=True,left=True,right=True);
plt.xlabel("Error relativo \%",fontsize=35,labelpad=15)
plt.ylabel("N",fontsize=35,labelpad=5)
plt.xlim(jmin+0.075,jmax-0.075)
plt.ylim(Nmin,Nmax)
ax.yaxis.set_major_locator(mticker.MultipleLocator(20));
ax.yaxis.set_minor_locator(mticker.MultipleLocator(5));
ax.xaxis.set_major_locator(mticker.MultipleLocator(0.1));
ax.xaxis.set_minor_locator(mticker.MultipleLocator(0.05));
plt.subplots_adjust(left=0.15,right=0.9,bottom=0.14,top=0.9)   # <--
plt.savefig("Jmin_"+folder+".eps",transparent=True)
#plt.show()

fig=plt.figure(figsize=(8,8))
ax=plt.subplot(111)
plt.plot(alpha_lat["l=0"],rho_lat["l=0"],'ko',markersize=12,alpha=0.75,label='$l=0$')
plt.plot(alpha_lat["l=1"],rho_lat["l=1"],'rs',markersize=12,alpha=0.75,label='$l=1$')
plt.plot(alpha_lat["l=2"],rho_lat["l=2"],'g^',markersize=12,alpha=0.75,label='$l=2$')
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

fig=plt.figure(figsize=(8,8))
ax=plt.subplot(111)
plt.plot(df2['Y'],df2['Iteration'],'bo',markersize=12,alpha=0.75)
ax.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
ax.tick_params(direction='in',which='minor',length=6,bottom=True,top=True,left=True,right=True);
plt.xlabel("Error relativo \%",fontsize=35,labelpad=15)
plt.ylabel("Iteración",fontsize=35,labelpad=15)
plt.xlim(jmin+0.075,jmax-0.075)
plt.ylim(imin,imax)
plt.hlines(initer,jmin,jmax,colors='k',linestyles='dashed')
plt.hlines(totevals,jmin,jmax,colors='k',linestyles='dashed')
ax.yaxis.set_major_locator(mticker.MultipleLocator(ymajor_imin));
ax.yaxis.set_minor_locator(mticker.MultipleLocator(yminor_imin));
ax.xaxis.set_major_locator(mticker.MultipleLocator(0.1));
ax.xaxis.set_minor_locator(mticker.MultipleLocator(0.05));
plt.subplots_adjust(left=0.15,right=0.9,bottom=0.14,top=0.9)   # <--
plt.savefig("imin_"+folder+".eps",transparent=True)
#plt.show()

########################################################################
#                             R A N D O M
########################################################################

folder="random"
total_loss_rand=pd.read_csv(folder+"/total_loss.dat",sep='\s+',header=None,usecols=[3])
total_loss_rand=total_loss_rand.rename(columns={3: "erp"});

evenlist=[]
oddlist=[]
for i in range(200): 
    if (i % 2 == 0): evenlist.append(i) 
    else: oddlist.append(i) 
alpha_rand=pd.read_csv(folder+"/params.dat",sep='\s+',header=None,skiprows=oddlist,usecols=[0,2,4,6])
rho_rand=pd.read_csv(folder+"/params.dat",sep='\s+',header=None,skiprows=evenlist,usecols=[0,2,4,6])
alpha_rand,rho_rand=rename_params(alpha_rand,rho_rand)
amin,amax,rmin,rmax=find_xtrms(alpha_rand,rho_rand)

filename="Be6CFGSTO_"+str(initer)+"GP"+str(maxevals)+".dat"
df=pd.read_csv(folder+"/seed_1/"+filename,sep='\t')
df2=pd.DataFrame(columns=df.columns)
for i in range(100):
    fseed="/seed_"+str(i+1)+"/"
    df=pd.read_csv(folder+fseed+filename,sep='\t')
    df2=df2.append(df.iloc[df['Y'].idxmin(), :], ignore_index=True)

print(folder)
print("max=",max(total_loss_rand['erp']))
print("min=",min(total_loss_rand['erp']))

########################################################################
#plt.hist(total_loss_rand['erp'],bins=10)
#plt.show()
fig=plt.figure(figsize=(8,8))
ax=plt.subplot(111)
bins=np.linspace(jmin,jmax,nbins+1)
plt.hist(total_loss_rand['erp'],bins=bins,color='orange')
ax.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
ax.tick_params(direction='in',which='minor',length=6,bottom=True,top=True,left=True,right=True);
plt.xlabel("Error relativo \%",fontsize=35,labelpad=15)
plt.ylabel("N",fontsize=35,labelpad=5)
plt.xlim(jmin+0.075,jmax-0.075)
plt.ylim(Nmin,Nmax)
ax.yaxis.set_major_locator(mticker.MultipleLocator(20));
ax.yaxis.set_minor_locator(mticker.MultipleLocator(5));
ax.xaxis.set_major_locator(mticker.MultipleLocator(0.1));
ax.xaxis.set_minor_locator(mticker.MultipleLocator(0.05));
plt.subplots_adjust(left=0.15,right=0.9,bottom=0.14,top=0.9)   # <--
plt.savefig("Jmin_"+folder+".eps",transparent=True)
#plt.show()

fig=plt.figure(figsize=(8,8))
ax=plt.subplot(111)
plt.plot(alpha_rand["l=0"],rho_rand["l=0"],'ko',markersize=12,alpha=0.75,label='$l=0$')
plt.plot(alpha_rand["l=1"],rho_rand["l=1"],'rs',markersize=12,alpha=0.75,label='$l=1$')
plt.plot(alpha_rand["l=2"],rho_rand["l=2"],'g^',markersize=12,alpha=0.75,label='$l=2$')
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

fig=plt.figure(figsize=(8,8))
ax=plt.subplot(111)
plt.plot(df2['Y'],df2['Iteration'],'bo',markersize=12,alpha=0.75)
ax.tick_params(direction='in',which='major',labelsize=30,length=8,bottom=True,top=True,left=True,right=True);
ax.tick_params(direction='in',which='minor',length=6,bottom=True,top=True,left=True,right=True);
plt.xlabel("Error relativo \%",fontsize=35,labelpad=15)
plt.ylabel("Iteración",fontsize=35,labelpad=15)
plt.xlim(jmin+0.075,jmax-0.075)
plt.ylim(imin,imax)
plt.hlines(initer,jmin,jmax,colors='k',linestyles='dashed')
plt.hlines(totevals,jmin,jmax,colors='k',linestyles='dashed')
ax.yaxis.set_major_locator(mticker.MultipleLocator(ymajor_imin));
ax.yaxis.set_minor_locator(mticker.MultipleLocator(yminor_imin));
ax.xaxis.set_major_locator(mticker.MultipleLocator(0.1));
ax.xaxis.set_minor_locator(mticker.MultipleLocator(0.05));
plt.subplots_adjust(left=0.15,right=0.9,bottom=0.14,top=0.9)   # <--
plt.savefig("imin_"+folder+".eps",transparent=True)
#plt.show()


