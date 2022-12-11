import numpy as np
import matplotlib.pyplot as plt

def lin_plot(x, y, c, s, labels, fnames, xlim=0, ylim=0, lw=2.0):
    fig, ax = plt.subplots()
    for i in range(len(x)):
        ax.plot(x[i], y[i], label=labels[i], linewidth=lw, color=c[i], linestyle=s[i])
        np.savetxt(fnames[i],np.array([x[i],y[i]]),delimiter=',')
    if not xlim == 0:
        ax.set(xlim=xlim)
    if not ylim == 0:
        ax.set(ylim=ylim)
    plt.legend()
    plt.show(block=True)

def logx_plot(x, y, c, s, labels, fnames, xlim=0, ylim=0, lw=2.0):
    fig, ax = plt.subplots()
    for i in range(len(x)):
        ax.semilogx(x[i], y[i], label=labels[i], linewidth=lw, color=c[i], linestyle=s[i])
        np.savetxt(fnames[i],np.array([x[i],y[i]]),delimiter=',')
    if not xlim == 0:
        ax.set(xlim=xlim)
    if not ylim == 0:
        ax.set(ylim=ylim)
    plt.legend()
    plt.show(block=True)

def logy_plot(x, y, c, s, labels, fnames, xlim=0, ylim=0, lw=2.0):
    fig, ax = plt.subplots()
    for i in range(len(x)):
        ax.semilogy(x[i], y[i], label=labels[i], linewidth=lw, color=c[i], linestyle=s[i])
        np.savetxt(fnames[i],np.array([x[i],y[i]]),delimiter=',')
    if not xlim == 0:
        ax.set(xlim=xlim)
    if not ylim == 0:
        ax.set(ylim=ylim)
    plt.legend()
    plt.show(block=True)
