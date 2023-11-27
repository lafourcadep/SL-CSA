from scipy.optimize import Bounds, differential_evolution
import warnings, sys
warnings.filterwarnings("ignore")
from sklearn.metrics import mean_squared_error
import numpy as np
from math import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import tkinter
mpl.use('TkAgg')

mpl.rcParams['ytick.labelsize']=12
mpl.rcParams['xtick.labelsize']=12
mpl.rcParams['axes.labelsize']=20
mpl.rcParams['legend.fontsize']=12
mpl.rcParams['axes.titlesize']=20

def regressor_pote(temp_traj, params):
    a=params[0]
    b=params[1]
    c=params[2]
    res=a*np.power(temp_traj,2)+b*temp_traj+c*np.ones(len(temp_traj))
    return res

def cost_function(params, temp_traj, pote_traj):
    err = 0.
    threshold=params[6]

    temp_traj_inf=temp_traj[temp_traj<=threshold]
    pote_traj_inf=pote_traj[temp_traj<=threshold]

    temp_traj_sup=temp_traj[temp_traj>threshold]
    pote_traj_sup=pote_traj[temp_traj>threshold]
    
    params_inf=[params[0], params[1], params[2]]
    new_temps=np.linspace(0,threshold,len(temp_traj_inf))
    pote_traj_inf_hat = regressor_pote(new_temps,params_inf)

    params_sup=[params[3], params[4], params[5]]
    new_temps=np.linspace(threshold,3000,len(temp_traj_sup))
    pote_traj_sup_hat = regressor_pote(new_temps,params_sup)

    err_inf = mean_squared_error(pote_traj_inf,pote_traj_inf_hat)
    err_sup = mean_squared_error(pote_traj_sup,pote_traj_sup_hat)    
    
    return err_sup+err_inf

tmelt_list=[]
natoms=816
fig, ax = plt.subplots(2,2,figsize=(16,16), sharex=True)

data = np.loadtxt('npt_Fe_bcc/thermodynamic_state_Fe_bcc.dat')
temp_traj=data[:,1]
pote_traj=(data[:,2]-np.mean(data[:100,2]))/natoms
bnds = [(-1,1), (-1,1), (0,0.01), (-1,1), (-1,1), (0,0.01), (1900,2100)]
res = differential_evolution(cost_function,bounds=bnds,args=(temp_traj,pote_traj),tol=1e-5,workers=16)

ainf=res.x[0]
binf=res.x[1]
cinf=res.x[2]
asup=res.x[3]
bsup=res.x[4]
csup=res.x[5]
tmelt=res.x[6]
tmelt_list.append(tmelt)

new_temp_inf=np.linspace(0,tmelt,1000)
new_pote_inf=regressor_pote(new_temp_inf,[ainf,binf,cinf])

new_temp_sup=np.linspace(tmelt,3000,1000)
new_pote_sup=regressor_pote(new_temp_sup,[asup,bsup,csup])

ax[0,0].axvspan(tmelt-100,tmelt+100,alpha=0.25,color='cyan')
ax[0,0].plot(temp_traj,pote_traj,'o', markersize=2,color='darkred',alpha=0.2)
ax[0,0].plot(new_temp_inf,new_pote_inf,'k--', linewidth=5.0)
ax[0,0].plot(new_temp_sup,new_pote_sup,'k--', linewidth=5.0)
ax[0,0].text(0.1,0.5,r'$T_{fus}=%d$ K' %(int(np.floor(tmelt))),fontsize=20,transform=ax[0,0].transAxes)

data = np.loadtxt('npt_Al_fcc/thermodynamic_state_Al_fcc.dat')
temp_traj=data[:,1]
pote_traj=(data[:,2]-np.mean(data[:100,2]))/natoms
bnds = [(-1,1), (-1,1), (0,0.01), (-1,1), (-1,1), (0,0.01), (800,1200)]
res = differential_evolution(cost_function,bounds=bnds,args=(temp_traj,pote_traj),tol=1e-5,workers=16)

ainf=res.x[0]
binf=res.x[1]
cinf=res.x[2]
asup=res.x[3]
bsup=res.x[4]
csup=res.x[5]
tmelt=res.x[6]
tmelt_list.append(tmelt)

new_temp_inf=np.linspace(0,tmelt,1000)
new_pote_inf=regressor_pote(new_temp_inf,[ainf,binf,cinf])

new_temp_sup=np.linspace(tmelt,3000,1000)
new_pote_sup=regressor_pote(new_temp_sup,[asup,bsup,csup])

ax[0,1].axvspan(tmelt-100,tmelt+100,alpha=0.25,color='cyan')
ax[0,1].plot(temp_traj,pote_traj,'o', markersize=2,color='darkorange',alpha=0.2)
ax[0,1].plot(new_temp_inf,new_pote_inf,'k--', linewidth=5.0)
ax[0,1].plot(new_temp_sup,new_pote_sup,'k--', linewidth=5.0)
ax[0,1].text(0.1,0.5,r'$T_{fus}=%d$ K' %(int(np.floor(tmelt))),fontsize=20,transform=ax[0,1].transAxes)

data = np.loadtxt('npt_Zr_hcp/thermodynamic_state_Zr_hcp.dat')
temp_traj=data[:,1]
pote_traj=(data[:,2]-np.mean(data[:100,2]))/natoms
bnds = [(-1,1), (-1,1), (0,0.01), (-1,1), (-1,1), (0,0.01), (2300,2600)]
res = differential_evolution(cost_function,bounds=bnds,args=(temp_traj,pote_traj),tol=1e-5,workers=16)

ainf=res.x[0]
binf=res.x[1]
cinf=res.x[2]
asup=res.x[3]
bsup=res.x[4]
csup=res.x[5]
tmelt=res.x[6]
tmelt_list.append(tmelt)

new_temp_inf=np.linspace(0,tmelt,1000)
new_pote_inf=regressor_pote(new_temp_inf,[ainf,binf,cinf])

new_temp_sup=np.linspace(tmelt,3000,1000)
new_pote_sup=regressor_pote(new_temp_sup,[asup,bsup,csup])

ax[1,0].axvspan(tmelt-100,tmelt+100,alpha=0.25,color='cyan')
ax[1,0].plot(temp_traj,pote_traj,'o', markersize=2,color='forestgreen',alpha=0.2)
ax[1,0].plot(new_temp_inf,new_pote_inf,'k--', linewidth=5.0)
ax[1,0].plot(new_temp_sup,new_pote_sup,'k--', linewidth=5.0)
ax[1,0].text(0.1,0.5,r'$T_{fus}=%d$ K' %(int(np.floor(tmelt))),fontsize=20,transform=ax[1,0].transAxes)

data = np.loadtxt('npt_Si_dia/thermodynamic_state_Si_dia.dat')
temp_traj=data[:,1]
pote_traj=(data[:,2]-np.mean(data[:100,2]))/natoms
bnds = [(-1,1), (-1,1), (0,0.01), (-1,1), (-1,1), (0,0.01), (2400,2600)]
res = differential_evolution(cost_function,bounds=bnds,args=(temp_traj,pote_traj),tol=1e-5,workers=16)

ainf=res.x[0]
binf=res.x[1]
cinf=res.x[2]
asup=res.x[3]
bsup=res.x[4]
csup=res.x[5]
tmelt=res.x[6]
tmelt_list.append(tmelt)

new_temp_inf=np.linspace(0,tmelt,1000)
new_pote_inf=regressor_pote(new_temp_inf,[ainf,binf,cinf])

new_temp_sup=np.linspace(tmelt,3000,1000)
new_pote_sup=regressor_pote(new_temp_sup,[asup,bsup,csup])

ax[1,1].axvspan(tmelt-100,tmelt+100,alpha=0.25,color='cyan')
ax[1,1].plot(temp_traj,pote_traj,'o', markersize=2,color='lightblue',alpha=0.2)
ax[1,1].plot(new_temp_inf,new_pote_inf,'k--', linewidth=5.0)
ax[1,1].plot(new_temp_sup,new_pote_sup,'k--', linewidth=5.0)
ax[1,1].text(0.1,0.5,r'$T_{fus}=%d$ K' %(int(np.floor(tmelt))),fontsize=20,transform=ax[1,1].transAxes)

ax[0,0].set_xlim([0,3000])
ax[0,1].set_xlim([0,3000])
ax[1,0].set_xlim([0,3000])
ax[1,1].set_xlim([0,3000])

ax[0,0].set_xlabel(r'Temperature (K)')
ax[1,0].set_xlabel(r'Temperature (K)')
ax[0,1].set_xlabel(r'Temperature (K)')
ax[1,1].set_xlabel(r'Temperature (K)')

ax[0,0].set_ylabel(r'$\Delta E$ (eV/atom)')
ax[0,1].set_ylabel(r'$\Delta E$ (eV/atom)')
ax[1,0].set_ylabel(r'$\Delta E$ (eV/atom)')
ax[1,1].set_ylabel(r'$\Delta E$ (eV/atom)')

ax[0,0].set_title(r'BCC Fe')
ax[0,1].set_title(r'FCC Al')
ax[1,0].set_title(r'HCP Zr')
ax[1,1].set_title(r'c-DIA Si')

plt.tight_layout()
plt.savefig('energy_vs_temperature.png')

print("%5.4f %5.4f %5.4f %5.4f" %(tmelt_list[0],tmelt_list[1],tmelt_list[2],tmelt_list[3]))
