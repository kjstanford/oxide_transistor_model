import csv
import numpy as np
from numpy import pi, log, log10, exp, sin
import pandas as pd
import matplotlib.pyplot as plt

def lin_plot(x, y, c, s, xlim=0, ylim=0, lw=2.0):
    fig, ax = plt.subplots()
    for i in range(len(x)):
        ax.plot(x[i], y[i], linewidth=lw, color=c[i%len(c)], linestyle=s[i%len(s)])
    if not xlim == 0:
        ax.set(xlim=xlim)
    if not ylim == 0:
        ax.set(ylim=ylim)
    plt.show(block=True)

def logx_plot(x, y, c, s, xlim=0, ylim=0, lw=2.0):
    fig, ax = plt.subplots()
    for i in range(len(x)):
        ax.semilogx(x[i], y[i], linewidth=lw, color=c[i%len(c)], linestyle=s[i%len(s)])
    if not xlim == 0:
        ax.set(xlim=xlim)
    if not ylim == 0:
        ax.set(ylim=ylim)
    plt.show(block=True)

def logy_plot(x, y, c, s, xlim=0, ylim=0, lw=2.0):
    fig, ax = plt.subplots()
    for i in range(len(x)):
        ax.semilogy(x[i], y[i], linewidth=lw, color=c[i%len(c)], linestyle=s[i%len(s)])
    if not xlim == 0:
        ax.set(xlim=xlim)
    if not ylim == 0:
        ax.set(ylim=ylim)
    plt.show(block=True)

def cascade_csv_reader(fname):
    full_data_set = []
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
    N1 = 0
    N2 = 0
    while line_count < len(full_data_set):
        row = full_data_set[line_count]
        if row[0] == 'Dimension1':
            print(line_count)
            N1 = int(row[1])
            line_count += 1
            row = full_data_set[line_count]
            N2 = int(row[1])
            line_count += 1
            temp_data = [[[str(X) for X in full_data_set[line_count][1:]]] for ii in range(N2)]
            line_count += 1
            for ii in range(N2):
                print(line_count)
                for elem in range(N1):
                    temp_data[ii].append([float(X) for X in full_data_set[line_count][1:]])
                    line_count += 1
            data_set+=temp_data
    print(len(data_set))

    pd_data_set = [pd.DataFrame(X[1:], columns=X[0]) for X in data_set]

    sweep_type = 1 
    df = pd_data_set[:]
    Id_data = [np.array(ds.loc[:sweep_type*int(N1/2)-1,' absId']) for ds in df]
    Vg_data = [np.array(ds.loc[:sweep_type*int(N1/2)-1,' Vg']) for ds in df]
    Vd_data = [np.array(ds.loc[:sweep_type*int(N1/2)-1,' Vd']) for ds in df]

    return Vg_data, Vd_data, Id_data

# print(cascade_csv_reader('ITO_IdVG_param_Vd [(1) _RoomT measurement for 2R sample with Vd param_; 12_18_2022 1_45_33 PM].csv'))