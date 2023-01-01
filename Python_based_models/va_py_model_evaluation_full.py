## Date:22ndNov 2022, Author: Koustav Jana ##
""" Code Description """ 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ITO_analytical_model_ver1 import ITO_analytical_model as model
from numpy import pi, log, log10, exp
from cascade_csv_reader_full import cascade_csv_reader, lin_plot, logx_plot, logy_plot
import fileinput as fi
  
lines = list(fi.input(files = 'evaluate_12_31_v1.lis'))

line_num = 0
first_word = lines[line_num].split()
first_word = None if len(first_word) == 0 else first_word[0]
acq_start = 0
acq_list = []
while True:
     if first_word == 'y':
          acq_start = 0
     if acq_start == 1:
          acq_list[-1].append(lines[line_num].split())
     if first_word == 'volt':
          line_num += 2
          acq_start = 1
          acq_list.append([])
     else:
          line_num += 1
     if line_num == len(lines):
          break
     first_word = lines[line_num].split()
     first_word = None if len(first_word) == 0 else first_word[0]

suffix_dict = {'a':1e-18, 'f':1e-15, 'p':1e-12, 'n':1e-9, 'u':1e-6, 'm':1e-3}
num_devices = len(acq_list[0][0]) - 1
print(num_devices)
Vg_acq_lists = [[float(x[0][:-1])*suffix_dict[x[0][-1]] if x[0][-1] in suffix_dict else float(x[0]) for x in xx] \
               for y in range(1, num_devices+1) for xx in acq_list]
Id_acq_lists = [[float(x[y][:-1])*suffix_dict[x[y][-1]] if x[y][-1] in suffix_dict else float(x[y]) for x in xx] \
               for y in range(1, num_devices+1) for xx in acq_list]
Vd_list = [0.2, 0.4, 0.6, 0.8, 1]
Vd_acq_lists = [[Vd for x in xx] for y in range(1, num_devices+1) for xx, Vd in zip(acq_list, Vd_list)]
print(len(Vd_acq_lists), len(Vg_acq_lists), len(Id_acq_lists))
     
          


# raise Exception("STOP HERE")

# constants
kB = 8.617e-5 #eV/K
q = 1.6e-19 #C
mo = 9.11e-31 #kg
hred = 1.054571817e-34
h = hred*2*pi
epso = 8.85e-12

#important parameters
L = 1.5
W = 20
mu_band = 39.5*1e-4 #band mobility in m^2/V.s
Ttraps = 500
Tdeep = 40000
Ntraps = 6.75e16
Nch = 1.34e17
T = 273.15+27
phiM = 5.2

frac = 1

params = dict(L=L, W=W, mu_band=mu_band, Ntraps=Ntraps*frac, Ttraps=Ttraps, Nch=Nch, T=T, tITO=4.5e-9,\
     tDE=5.3e-9, me_factor=0.3, kDE=16, kITO=9, phiM=phiM, chiS=4.3, Rs=500/W, Rd=500/W,\
     Ndeep=Ntraps*(1-frac), Tdeep=Tdeep)

Cox = epso*params['kDE']/params['tDE']

ax_Vg_exp, ax_Id_exp, ax_mu_exp, ax_SS_exp, ax_Vd_exp = [], [], [], [], []
fnames = ['ITO_IdVG_param_Vd [(1) _RoomT measurement for 2R sample with Vd param_; 12_18_2022 1_45_33 PM].csv']
# fnames = ['ITO_IdVG [(26) _RoomT measurement for 2R sample_; 12_18_2022 1_33_51 PM].csv']
# fnames = ['ITO_IdVG [(10) _T = -20C 2R sample_; 12_18_2022 3_17_46 PM].csv']
# fnames = ['ITO_IdVG [(14) _T = -10C 2R sample_; 12_18_2022 3_29_42 PM].csv']
# fnames = ['ITO_IdVG [(17) _T=0C 2R sample_; 12_18_2022 3_35_50 PM].csv']
# fnames = ['ITO_IdVG [(19) _T=10C 2R sample_; 12_18_2022 3_41_18 PM].csv']
# fnames = ['ITO_IdVG [(21) _T=20C 2R sample_; 12_18_2022 3_46_38 PM].csv']
# fnames = ['ITO_IdVG [(24) _T=30C 2R sample_; 12_18_2022 3_55_33 PM].csv']
# fnames = ['ITO_IdVG [(26) _T=40C 2R sample_; 12_18_2022 4_01_19 PM].csv']
# fnames = ['ITO_IdVG [(31) _T=50C 2R sample_; 12_18_2022 4_12_41 PM].csv']
# fnames = ['ITO_IdVG [(40) _T=60C 2R sample_; 12_18_2022 4_27_05 PM].csv']
# fnames = ['ITO_IdVG [(47) _T=70C 2R sample_; 12_18_2022 4_46_32 PM].csv']
# fnames = ['ITO_IdVG [(53) _T=80C 2R sample_; 12_18_2022 4_57_27 PM].csv']
# fnames = ['ITO_IdVG [(59) _T=90C 2R sample_; 12_18_2022 5_10_46 PM].csv']
# fnames = ['ITO_IdVG [(64) _T=100C 2R sample_; 12_18_2022 5_21_57 PM].csv']
# fnames = ['ITO_IdVG [(69) _T=110C 2R sample_; 12_18_2022 5_33_53 PM].csv']
# fnames = ['ITO_IdVG [(75) _T=120C 2R sample_; 12_18_2022 5_52_21 PM].csv']
# fnames = ['ITO_IdVG [(81) _T=130C 2R sample_; 12_18_2022 6_07_31 PM].csv']


for fname in fnames:
     Vg_data, Vd_data, Id_data = cascade_csv_reader(fname)
     epso = 8.85e-12
     Vg_min = -0.5
     L=1.5; W =20; 
     kDE=16; tDE=5.3e-9; Cox = epso*kDE/tDE;

     ax_Vg_exp += [d_Vg[d_Vg >= Vg_min] for d_Vg in Vg_data[:]]
     ax_Id_exp += [d_Id[d_Vg >= Vg_min] for d_Vg, d_Id in zip(Vg_data[:], Id_data[:])]
     ax_Vd_exp += [d_Vd[d_Vg >= Vg_min] for d_Vg, d_Vd in zip(Vg_data[:], Vd_data[:])]
     ax_mu_exp += \
          [1e4*(np.gradient(d_Id[d_Vg >= Vg_min])/np.gradient(d_Vg[d_Vg >= Vg_min]))/(Cox*(W/L)*d_Vd[d_Vg >= Vg_min])\
          for d_Vg, d_Vd, d_Id in zip(Vg_data[:], Vd_data[:], Id_data[:])]
     ax_SS_exp += \
          [1e3*(np.gradient(d_Vg[d_Vg >= Vg_min])/np.gradient(log10(d_Id[d_Vg >= Vg_min])))\
          for d_Vg, d_Id in zip(Vg_data[:], Id_data[:])]

ax_Vg_sim, ax_Vd_sim, ax_Id_sim, ax_mu_sim, ax_SS_sim = [], [], [], [], []
# ax_Vg_sim = ax_Vg_exp
# ax_Vd_sim = ax_Vd_exp
# ax_Id_sim = [ np.array([model(params=params,inputs=dict(Vgs=Vg,Vds=Vd))[1] for Vg, Vd in zip(Vg_arr, Vd_arr)])\
#              for Vd_arr, Vg_arr in zip(ax_Vd_sim, ax_Vg_sim) ]
# ax_mu_sim = [1e4*(np.gradient(Id_arr)/np.gradient(Vg_arr))/(Cox*(W/L)*Vd_arr)\
#              for Vd_arr, Vg_arr, Id_arr in zip(ax_Vd_sim, ax_Vg_sim, ax_Id_sim)]
# ax_SS_sim = [(np.gradient(Vg_arr)*1e3)/np.gradient(log10(Id_arr))\
#              for Vg_arr, Id_arr in zip(ax_Vg_sim, ax_Id_sim)]

ax_Vg_hspice, ax_Vd_hspice, ax_Id_hspice, ax_mu_hspice, ax_SS_hspice = [], [], [], [], []
ax_Vg_hspice = [np.array(Vg_acq_list) for Vg_acq_list in Vg_acq_lists]
ax_Vd_hspice = [np.array(Vd_acq_list) for Vd_acq_list in Vd_acq_lists]
ax_Id_hspice = [np.array(Id_acq_list) for Id_acq_list in Id_acq_lists]
ax_mu_hspice = [1e4*(np.gradient(Id_arr)/np.gradient(Vg_arr))/(Cox*(W/L)*Vd_arr)\
             for Vd_arr, Vg_arr, Id_arr in zip(ax_Vd_hspice, ax_Vg_hspice, ax_Id_hspice)]
ax_SS_hspice = [(np.gradient(Vg_arr)*1e3)/np.gradient(log10(Id_arr))\
             for Vg_arr, Id_arr in zip(ax_Vg_hspice, ax_Id_hspice)] 

ax_Vg = ax_Vg_exp + ax_Vg_sim + ax_Vg_hspice
ax_Id = ax_Id_exp + ax_Id_sim + ax_Id_hspice
ax_mu = ax_mu_exp + ax_mu_sim + ax_mu_hspice
ax_SS = ax_SS_exp + ax_SS_sim + ax_SS_hspice

c = [yy for xx in ['solid', 'dotted', 'dashed'] for yy in ['r', 'b', 'g', 'k', 'm']] 
s = [xx for xx in ['solid', 'dotted', 'dashed'] for yy in ['r', 'b', 'g', 'k', 'm']]
mask = [yy for xx in range(3) for yy in [False, False, False, False, True]]
print(c, s)
lin_plot(x=ax_Vg, y=ax_Id, c=c, s=s, mask=mask)
logy_plot(x=ax_Vg, y=ax_Id, c=c, s=s, mask=mask)
lin_plot(x=ax_Vg, y=ax_mu, c=c, s=s, mask=mask)
logy_plot(x=ax_Vg, y=ax_mu, c=c, s=s, mask=mask)
logx_plot(x=ax_Id, y=ax_SS, c=c, s=s, mask=mask,\
    xlim=(W*1e-13, W*1e-7), ylim=(45, 200))

# print(ax_Id_sim[0])










## ignore below
# # df = [pd.read_excel('IdVg L2W100 Vd=0.1V.xlsx',usecols=[0,2]),\
# #     pd.read_excel('IdVg L2W100 Vd=1V.xlsx',usecols=[0,2])]
# df = [pd.read_excel('IdVg L2W100 Vd=0.1V.xlsx',usecols=[0,2])]

# ds = [np.array(X) for X in df]
# # Vds = [0.1, 1]
# Vds = [0.1]

# rng = range(201);
# Vgmin = -0.5
# ax_Vg_exp = [X[rng,0][X[rng,0] > Vgmin] for X in ds]
# ax_Id_exp = [X[rng,1][X[rng,0] > Vgmin] for X in ds]
# ax_mu_exp = [1e4*(np.gradient(X[rng,1][X[rng,0] > Vgmin])/np.gradient(X[rng,0][X[rng,0] > Vgmin]))/(Cox*(W/L)*Vds[i])\
#              for i, X in enumerate(ds)]
# ax_SS_exp = [(np.gradient(X[rng,0][X[rng,0] > Vgmin]*1e3)/np.gradient(log10(X[rng,1][X[rng,0] > Vgmin])))\
#              for X in ds]

# ax_Vg_sim = ax_Vg_exp
# ax_Id_sim = [ np.array([model(params=params,inputs=dict(Vgs=Vg,Vds=Vd))[1] for Vg in Vg_arr])\
#              for Vd, Vg_arr in zip(Vds, ax_Vg_sim) ]
# ax_mu_sim = [1e4*(np.gradient(Id_arr)/np.gradient(Vg_arr))/(Cox*(W/L)*Vd)\
#              for Vd, Vg_arr, Id_arr in zip(Vds, ax_Vg_sim, ax_Id_sim)]
# ax_SS_sim = [(np.gradient(Vg_arr)*1e3)/np.gradient(log10(Id_arr))\
#              for Vg_arr, Id_arr in zip(ax_Vg_sim, ax_Id_sim)]

# ax_Vg = ax_Vg_exp + ax_Vg_sim
# ax_Id = ax_Id_exp + ax_Id_sim
# ax_mu = ax_mu_exp + ax_mu_sim
# ax_SS = ax_SS_exp + ax_SS_sim
    
# lin_plot(x=ax_Vg, y=ax_Id, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
#      ylim=(1e-12, ax_Id[1][-1]*1.25), labels=['Exp', 'Model'], fnames=['lin_e.csv','lin_m.csv'])    
# logy_plot(x=ax_Vg, y=ax_Id, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
#      ylim=(1e-12, ax_Id[1][-1]*1.25), labels=['Exp', 'Model'], fnames=['log_e.csv','log_m.csv'])    
# lin_plot(x=ax_Vg, y=ax_mu, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
#      labels=['Exp', 'Model'], fnames=['mu_e.csv','mu_m.csv'])    
# logx_plot(x=ax_Id, y=ax_SS, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
#      xlim=(1e-11, 1e-6), ylim=(50,150), labels=['Exp', 'Model'], fnames=['SS_e.csv','SS_m.csv'])


# # mu_band = 33.5*1e-4; #band mobility in m^2/V.s
# # Ttraps = [400, 800, 1200]
# # # Ttraps = 400;
# # # Ntraps = [1.5e16, 5.5e16, 9.5e16];
# # Ntraps = 5.5e16;
# # Nch = 1e17;
# # Vds = 0.1;
# # params = dict(L=L, W=W, mu_band=mu_band, Ntraps=Ntraps, Nch=Nch, T=T, tITO=4.5e-9,\
# #     tDE=5.3e-9, me_factor=0.3, kDE=16, kITO=9, phiM=phiM, chiS=4.3)
# # pars = [params|{'Ttraps':elem} for elem in Ttraps];

# # ax_Vg, ax_Id, ax_mu, ax_SS, ax_phi = [], [], [], [], []
# # ax_Vg = [ax_Vg_exp[0] for elem in Ttraps]
# # ax_Id = [np.array([model(params=par,inputs=dict(Vgs=Vg,Vds=Vds))[1] for Vg in Vg_arr]) \
# #          for Vg_arr, par in zip(ax_Vg, pars)]
# # ax_mu = [1e4*(np.gradient(Id_arr)/np.gradient(Vg_arr))/(Cox*(W/L)*Vds)\
# #          for Vg_arr, Id_arr in zip(ax_Vg, ax_Id)]
# # ax_SS = [(np.gradient(Vg_arr)*1e3)/np.gradient(log10(Id_arr))\
# #          for Vg_arr, Id_arr in zip(ax_Vg, ax_Id)]


# # lin_plot(x=ax_Vg, y=ax_Id, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
# #      ylim=(1e-12, ax_Id[1][-1]*1.25), labels=[str(elem) for elem in Ttraps],\
# #            fnames=['lin_1.csv','lin_2.csv','lin_3.csv'])    
# # logy_plot(x=ax_Vg, y=ax_Id, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
# #      ylim=(1e-12, ax_Id[1][-1]*1.25), labels=[str(elem) for elem in Ttraps],\
# #            fnames=['log_1.csv','log_2.csv','log_3.csv'])    
# # lin_plot(x=ax_Vg, y=ax_mu, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
# #      labels=[str(elem) for elem in Ttraps], fnames=['mu_1.csv','mu_2.csv','mu_3.csv'])    
# # logx_plot(x=ax_Id, y=ax_SS, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
# #      xlim=(1e-11, 1e-6), ylim=(50,150), labels=[str(elem) for elem in Ttraps],\
# #            fnames=['SS_1.csv','SS_2.csv','SS_3.csv'])