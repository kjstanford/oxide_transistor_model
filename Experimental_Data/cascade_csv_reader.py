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

    return Vg_data, Vd_data, Id_data