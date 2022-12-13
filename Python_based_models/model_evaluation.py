## Date:22ndNov 2022, Author: Koustav Jana ##
""" Code Description """ 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ITO_analytical_model_ver1_Rs import ITO_analytical_model as model
from plotter import lin_plot, logx_plot, logy_plot
from numpy import pi, log, log10, exp

# constants
kB = 8.617e-5 #eV/K
q = 1.6e-19 #C
mo = 9.11e-31 #kg
hred = 1.054571817e-34
h = hred*2*pi
epso = 8.85e-12

#important parameters
L = 1
W = 100
mu_band = 33.5*1e-4 #band mobility in m^2/V.s
Ttraps = 400
Ntraps = 5.5e16
Nch = 1e17
T = 300
phiM = 5.2

params = dict(L=L, W=W, mu_band=mu_band, Ntraps=Ntraps, Ttraps=Ttraps, Nch=Nch, T=T, tITO=4.5e-9,\
    tDE=5.3e-9, me_factor=0.3, kDE=16, kITO=9, phiM=phiM, chiS=4.3, Rs=1, Rd=1, Ndeep=Ntraps, Tdeep=50*Ntraps)

Cox = epso*params['kDE']/params['tDE'];

# df = [pd.read_excel('IdVg L2W100 Vd=0.1V.xlsx',usecols=[0,2]),\
#     pd.read_excel('IdVg L2W100 Vd=1V.xlsx',usecols=[0,2])]
df = [pd.read_excel('IdVg L2W100 Vd=0.1V.xlsx',usecols=[0,2])]

ds = [np.array(X) for X in df]
# Vds = [0.1, 1]
Vds = [0.1]

rng = range(201);
Vgmin = -0.5
ax_Vg_exp = [X[rng,0][X[rng,0] > Vgmin] for X in ds]
ax_Id_exp = [X[rng,1][X[rng,0] > Vgmin] for X in ds]
ax_mu_exp = [1e4*(np.gradient(X[rng,1][X[rng,0] > Vgmin])/np.gradient(X[rng,0][X[rng,0] > Vgmin]))/(Cox*(W/L)*Vds[i])\
             for i, X in enumerate(ds)]
ax_SS_exp = [(np.gradient(X[rng,0][X[rng,0] > Vgmin]*1e3)/np.gradient(log10(X[rng,1][X[rng,0] > Vgmin])))\
             for X in ds]

ax_Vg_sim = ax_Vg_exp
ax_Id_sim = [ np.array([model(params=params,inputs=dict(Vgs=Vg,Vds=Vd))[1] for Vg in Vg_arr])\
             for Vd, Vg_arr in zip(Vds, ax_Vg_sim) ]
ax_mu_sim = [1e4*(np.gradient(Id_arr)/np.gradient(Vg_arr))/(Cox*(W/L)*Vd)\
             for Vd, Vg_arr, Id_arr in zip(Vds, ax_Vg_sim, ax_Id_sim)]
ax_SS_sim = [(np.gradient(Vg_arr)*1e3)/np.gradient(log10(Id_arr))\
             for Vg_arr, Id_arr in zip(ax_Vg_sim, ax_Id_sim)]

ax_Vg = ax_Vg_exp + ax_Vg_sim
ax_Id = ax_Id_exp + ax_Id_sim
ax_mu = ax_mu_exp + ax_mu_sim
ax_SS = ax_SS_exp + ax_SS_sim
    
lin_plot(x=ax_Vg, y=ax_Id, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
     ylim=(1e-12, ax_Id[1][-1]*1.25), labels=['Exp', 'Model'], fnames=['lin_e.csv','lin_m.csv'])    
logy_plot(x=ax_Vg, y=ax_Id, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
     ylim=(1e-12, ax_Id[1][-1]*1.25), labels=['Exp', 'Model'], fnames=['log_e.csv','log_m.csv'])    
lin_plot(x=ax_Vg, y=ax_mu, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
     labels=['Exp', 'Model'], fnames=['mu_e.csv','mu_m.csv'])    
logx_plot(x=ax_Id, y=ax_SS, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
     xlim=(1e-11, 1e-6), ylim=(50,150), labels=['Exp', 'Model'], fnames=['SS_e.csv','SS_m.csv'])


# mu_band = 33.5*1e-4; #band mobility in m^2/V.s
# Ttraps = [400, 800, 1200]
# # Ttraps = 400;
# # Ntraps = [1.5e16, 5.5e16, 9.5e16];
# Ntraps = 5.5e16;
# Nch = 1e17;
# Vds = 0.1;
# params = dict(L=L, W=W, mu_band=mu_band, Ntraps=Ntraps, Nch=Nch, T=T, tITO=4.5e-9,\
#     tDE=5.3e-9, me_factor=0.3, kDE=16, kITO=9, phiM=phiM, chiS=4.3)
# pars = [params|{'Ttraps':elem} for elem in Ttraps];

# ax_Vg, ax_Id, ax_mu, ax_SS, ax_phi = [], [], [], [], []
# ax_Vg = [ax_Vg_exp[0] for elem in Ttraps]
# ax_Id = [np.array([model(params=par,inputs=dict(Vgs=Vg,Vds=Vds))[1] for Vg in Vg_arr]) \
#          for Vg_arr, par in zip(ax_Vg, pars)]
# ax_mu = [1e4*(np.gradient(Id_arr)/np.gradient(Vg_arr))/(Cox*(W/L)*Vds)\
#          for Vg_arr, Id_arr in zip(ax_Vg, ax_Id)]
# ax_SS = [(np.gradient(Vg_arr)*1e3)/np.gradient(log10(Id_arr))\
#          for Vg_arr, Id_arr in zip(ax_Vg, ax_Id)]


# lin_plot(x=ax_Vg, y=ax_Id, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
#      ylim=(1e-12, ax_Id[1][-1]*1.25), labels=[str(elem) for elem in Ttraps],\
#            fnames=['lin_1.csv','lin_2.csv','lin_3.csv'])    
# logy_plot(x=ax_Vg, y=ax_Id, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
#      ylim=(1e-12, ax_Id[1][-1]*1.25), labels=[str(elem) for elem in Ttraps],\
#            fnames=['log_1.csv','log_2.csv','log_3.csv'])    
# lin_plot(x=ax_Vg, y=ax_mu, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
#      labels=[str(elem) for elem in Ttraps], fnames=['mu_1.csv','mu_2.csv','mu_3.csv'])    
# logx_plot(x=ax_Id, y=ax_SS, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
#      xlim=(1e-11, 1e-6), ylim=(50,150), labels=[str(elem) for elem in Ttraps],\
#            fnames=['SS_1.csv','SS_2.csv','SS_3.csv'])