import numpy as np
from numpy import pi, log, log10, exp, sin
from cascade_csv_reader import cascade_csv_reader, lin_plot, logx_plot, logy_plot

fnames = [ 'ITO_IdVG [(10) _T = -20C 2R sample_; 12_18_2022 3_17_46 PM].csv' ,
        'ITO_IdVG [(14) _T = -10C 2R sample_; 12_18_2022 3_29_42 PM].csv' ,
        'ITO_IdVG [(17) _T=0C 2R sample_; 12_18_2022 3_35_50 PM].csv' ,
        'ITO_IdVG [(19) _T=10C 2R sample_; 12_18_2022 3_41_18 PM].csv' ,
        'ITO_IdVG [(21) _T=20C 2R sample_; 12_18_2022 3_46_38 PM].csv' ,
        'ITO_IdVG [(24) _T=30C 2R sample_; 12_18_2022 3_55_33 PM].csv' ,
        'ITO_IdVG [(26) _T=40C 2R sample_; 12_18_2022 4_01_19 PM].csv' ,
        'ITO_IdVG [(31) _T=50C 2R sample_; 12_18_2022 4_12_41 PM].csv' ,
        'ITO_IdVG [(40) _T=60C 2R sample_; 12_18_2022 4_27_05 PM].csv' ,
        'ITO_IdVG [(47) _T=70C 2R sample_; 12_18_2022 4_46_32 PM].csv' ,
        'ITO_IdVG [(53) _T=80C 2R sample_; 12_18_2022 4_57_27 PM].csv' ,
        'ITO_IdVG [(59) _T=90C 2R sample_; 12_18_2022 5_10_46 PM].csv' ,
        'ITO_IdVG [(64) _T=100C 2R sample_; 12_18_2022 5_21_57 PM].csv' ,
        'ITO_IdVG [(69) _T=110C 2R sample_; 12_18_2022 5_33_53 PM].csv' ,
        'ITO_IdVG [(75) _T=120C 2R sample_; 12_18_2022 5_52_21 PM].csv' ,
        'ITO_IdVG [(81) _T=130C 2R sample_; 12_18_2022 6_07_31 PM].csv' ]

def evaluate_SS(ax_Vg, ax_Id, Id_UL, Id_LL):
    SS = []
    for Vg_data, Id_data in zip(ax_Vg, ax_Id):
        Vg_data = Vg_data[Id_data <= Id_UL]
        Id_data = Id_data[Id_data <= Id_UL]
        Vg_data = Vg_data[Id_data >= Id_LL]
        Id_data = Id_data[Id_data >= Id_LL]
        ans = 1e3*( Vg_data[-1] - Vg_data[0] ) / ( log10(Id_data[-1]) - log10(Id_data[0]) )
        SS.append(ans)
    return SS

T = -20
SS_T = []
for fname in fnames:
    Vg_data, Vd_data, Id_data = cascade_csv_reader(fname)
    # float(input('Id_UL = '))
    # float(input('Id_LL = '))
    list_SS = evaluate_SS(Vg_data, Id_data, 5e-11, 9e-12)
    print(T, sum(list_SS)/len(list_SS), list_SS)
    SS_T.append(sum(list_SS)/len(list_SS))
    T += 10
print(SS_T)

raise Exception('Stop Here')


# fnames = ['ITO_IdVG [(26) _RoomT measurement for 2R sample_; 12_18_2022 1_33_51 PM].csv']
fnames = ['ITO_IdVG [(10) _T = -20C 2R sample_; 12_18_2022 3_17_46 PM].csv']
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

ax_Vg_exp, ax_Id_exp, ax_mu_exp, ax_SS_exp, ax_Vd_exp = [], [], [], [], []
for fname in fnames:
    Vg_data, Vd_data, Id_data = cascade_csv_reader(fname)
    epso = 8.85e-12
    Vg_min = -1.2
    L=1.5; W =20; 
    kDE=16; tDE=5.3e-9; Cox = epso*kDE/tDE;

    ax_Vg_exp += [d_Vg[d_Vg >= Vg_min] for d_Vg in Vg_data]
    ax_Id_exp += [d_Id[d_Vg >= Vg_min] for d_Vg, d_Id in zip(Vg_data, Id_data)]
    ax_Vd_exp += [d_Vd[d_Vg >= Vg_min] for d_Vg, d_Vd in zip(Vg_data, Vd_data)]
    ax_mu_exp += \
        [1e4*(np.gradient(d_Id[d_Vg >= Vg_min])/np.gradient(d_Vg[d_Vg >= Vg_min]))/(Cox*(W/L)*d_Vd[d_Vg >= Vg_min][0])\
        for d_Vg, d_Vd, d_Id in zip(Vg_data, Vd_data, Id_data)]
    ax_SS_exp += \
          [1e3*(np.gradient(d_Vg[d_Vg >= Vg_min])/np.gradient(log10(d_Id[d_Vg >= Vg_min])))\
          for d_Vg, d_Id in zip(Vg_data, Id_data)]
    

lin_plot(x=ax_Vg_exp, y=ax_Id_exp, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'])
logy_plot(x=ax_Vg_exp, y=ax_Id_exp, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'])
lin_plot(x=ax_Vg_exp, y=ax_mu_exp, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'])
logx_plot(x=ax_Id_exp, y=ax_SS_exp, c=['r','b','r','b'], s=['solid','solid','dashdot','dashdot'],\
    xlim=(W*1e-13, W*1e-7), ylim=(45, 200))

