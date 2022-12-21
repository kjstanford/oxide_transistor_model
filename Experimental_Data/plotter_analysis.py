import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy import pi, log, log10, exp, sin
import csv

def lin_plot(x, y, c, s, xlim=0, ylim=0, lw=2.0):
    fig, ax = plt.subplots()
    for i in range(len(x)):
        ax.plot(x[i], y[i], linewidth=lw, color=c[i%len(c)], linestyle=s[i%len(s)])
    if not xlim == 0:
        ax.set(xlim=xlim)
    if not ylim == 0:
        ax.set(ylim=ylim)
    plt.show()

def logx_plot(x, y, c, s, xlim=0, ylim=0, lw=2.0):
    fig, ax = plt.subplots()
    for i in range(len(x)):
        ax.semilogx(x[i], y[i], linewidth=lw, color=c[i%len(c)], linestyle=s[i%len(s)])
    if not xlim == 0:
        ax.set(xlim=xlim)
    if not ylim == 0:
        ax.set(ylim=ylim)
    plt.show()

def logy_plot(x, y, c, s, xlim=0, ylim=0, lw=2.0):
    fig, ax = plt.subplots()
    for i in range(len(x)):
        ax.semilogy(x[i], y[i], linewidth=lw, color=c[i%len(c)], linestyle=s[i%len(s)])
    if not xlim == 0:
        ax.set(xlim=xlim)
    if not ylim == 0:
        ax.set(ylim=ylim)
    plt.show()

full_data_set = []
fname = 'ITO_IdVG [(26) _RoomT measurement for 2R sample_; 12_18_2022 1_33_51 PM].csv'
# fname = 'ITO_IdVG [(10) _T = -20C 2R sample_; 12_18_2022 3_17_46 PM].csv'
# fname = 'ITO_IdVG [(14) _T = -10C 2R sample_; 12_18_2022 3_29_42 PM].csv'
# fname = 'ITO_IdVG [(17) _T=0C 2R sample_; 12_18_2022 3_35_50 PM].csv'
# fname = 'ITO_IdVG [(19) _T=10C 2R sample_; 12_18_2022 3_41_18 PM].csv'
# fname = 'ITO_IdVG [(21) _T=20C 2R sample_; 12_18_2022 3_46_38 PM].csv'
# fname = 'ITO_IdVG [(24) _T=30C 2R sample_; 12_18_2022 3_55_33 PM].csv'
# fname = 'ITO_IdVG [(26) _T=40C 2R sample_; 12_18_2022 4_01_19 PM].csv'
# fname = 'ITO_IdVG [(31) _T=50C 2R sample_; 12_18_2022 4_12_41 PM].csv'
# fname = 'ITO_IdVG [(40) _T=60C 2R sample_; 12_18_2022 4_27_05 PM].csv'
# fname = 'ITO_IdVG [(47) _T=70C 2R sample_; 12_18_2022 4_46_32 PM].csv'
# fname = 'ITO_IdVG [(53) _T=80C 2R sample_; 12_18_2022 4_57_27 PM].csv'
# fname = 'ITO_IdVG [(59) _T=90C 2R sample_; 12_18_2022 5_10_46 PM].csv'
# fname = 'ITO_IdVG [(64) _T=100C 2R sample_; 12_18_2022 5_21_57 PM].csv'
# fname = 'ITO_IdVG [(69) _T=110C 2R sample_; 12_18_2022 5_33_53 PM].csv'
# fname = 'ITO_IdVG [(75) _T=120C 2R sample_; 12_18_2022 5_52_21 PM].csv'
# fname = 'ITO_IdVG [(81) _T=130C 2R sample_; 12_18_2022 6_07_31 PM].csv'

with open (fname, mode='r') as csv_file:
    csv_reader = csv.reader(csv_file)
    line_count = 0
    for row in csv_reader:
        # print(line_count, row)
        if row[0] in ['Dimension1', 'Dimension2', 'DataName', 'DataValue']:
            full_data_set.append(row)
        line_count += 1
print(f'Processed {line_count} lines')

line_count = 0
data_set = []
N = 0;
while line_count < len(full_data_set):
    row = full_data_set[line_count]
    if row[0] == 'Dimension1':
        N = int(row[1])
        data_set.append([])
        temp_data = []
        line_count += 2
        temp_data.append([str(X) for X in full_data_set[line_count][1:]])
        line_count += 1
        for elem in range(N):
            temp_data.append([float(X) for X in full_data_set[line_count][1:]])
            line_count += 1
        data_set[-1] = temp_data
print(len(data_set))

pd_data_set = [pd.DataFrame(X[1:], columns=X[0]) for X in data_set]

sweep_type = 1 
df = pd_data_set[:]
Id_data = [np.array(ds.loc[:sweep_type*int(N/2)-1,' absId']) for ds in df]
Vg_data = [np.array(ds.loc[:sweep_type*int(N/2)-1,' Vg']) for ds in df]
Vd_data = [np.array(ds.loc[:sweep_type*int(N/2)-1,' Vd']) for ds in df]

epso = 8.85e-12
Vg_min = -2
L=1.5; W =20; 
kDE=16; tDE=5.3e-9; Cox = epso*kDE/tDE;

ax_Vg_exp = [d_Vg[d_Vg >= Vg_min] for d_Vg in Vg_data]
ax_Id_exp = [d_Id[d_Vg >= Vg_min] for d_Vg, d_Id in zip(Vg_data, Id_data)]
ax_Vd_exp = [d_Vd[d_Vg >= Vg_min] for d_Vg, d_Vd in zip(Vg_data, Vd_data)]
ax_mu_exp = \
        [1e4*(np.gradient(d_Id[d_Vg >= Vg_min])/np.gradient(d_Vg[d_Vg >= Vg_min]))/(Cox*(W/L)*d_Vd[d_Vg >= Vg_min][0])\
        for d_Vg, d_Vd, d_Id in zip(Vg_data, Vd_data, Id_data)]



lin_plot(x=ax_Vg_exp, y=ax_Id_exp, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'])
logy_plot(x=ax_Vg_exp, y=ax_Id_exp, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'])
lin_plot(x=ax_Vg_exp, y=ax_mu_exp, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'])

print(max(ax_mu_exp[0][:-2]))


# ax_SS_exp = [1e3*(np.gradient(d_Vg)/np.gradient(log10(d_Id)))\
#      for d_Vg, d_Id in zip(ax_Vg_exp, ax_Id_exp)]

# logx_plot(x=ax_Id_exp, y=ax_SS_exp, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
#     xlim=(W*1e-13, W*1e-7), ylim=(55, 200))

