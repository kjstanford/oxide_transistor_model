import numpy as np
from numpy import pi, log, log10, exp, sin
from cascade_csv_reader import cascade_csv_reader, lin_plot, logx_plot, logy_plot

fnames = ['ITO_IdVG [(26) _RoomT measurement for 2R sample_; 12_18_2022 1_33_51 PM].csv',
'ITO_IdVG [(10) _T = -20C 2R sample_; 12_18_2022 3_17_46 PM].csv']
# 'ITO_IdVG [(14) _T = -10C 2R sample_; 12_18_2022 3_29_42 PM].csv']

ax_Vg_exp, ax_Id_exp, ax_mu_exp, ax_SS_exp, ax_Vd_exp = [], [], [], [], []
for fname in fnames:
    Vg_data, Vd_data, Id_data = cascade_csv_reader(fname)
    epso = 8.85e-12
    Vg_min = -2
    L=1.5; W =20; 
    kDE=16; tDE=5.3e-9; Cox = epso*kDE/tDE;

    ax_Vg_exp += [d_Vg[d_Vg >= Vg_min] for d_Vg in Vg_data]
    ax_Id_exp += [d_Id[d_Vg >= Vg_min] for d_Vg, d_Id in zip(Vg_data, Id_data)]
    ax_Vd_exp += [d_Vd[d_Vg >= Vg_min] for d_Vg, d_Vd in zip(Vg_data, Vd_data)]
    ax_mu_exp += \
        [1e4*(np.gradient(d_Id[d_Vg >= Vg_min])/np.gradient(d_Vg[d_Vg >= Vg_min]))/(Cox*(W/L)*d_Vd[d_Vg >= Vg_min][0])\
        for d_Vg, d_Vd, d_Id in zip(Vg_data, Vd_data, Id_data)]
    

lin_plot(x=ax_Vg_exp, y=ax_Id_exp, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'])
logy_plot(x=ax_Vg_exp, y=ax_Id_exp, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'])
lin_plot(x=ax_Vg_exp, y=ax_mu_exp, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'])



