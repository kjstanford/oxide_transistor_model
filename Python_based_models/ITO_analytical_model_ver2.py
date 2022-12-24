## Date:22ndNov 2022, Author: Koustav Jana ##
""" Code Description: Dual Exponetial Trap Distribution + Self-Written NR code """ 
import numpy as np
import scipy.linalg as la
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from numpy import pi, exp, log, sin, log10
from scipy.special import hyp2f1
from scipy.optimize import fsolve

def NR_fsolve(f, df, z0):
    # print('Solving....')
    z = z0
    tol = 1e-10
    alpha = 0.8
    zprev = float('inf')
    c = 0
    while abs(z-zprev)>tol and c<200:
        zprev = z
        z = zprev - alpha*f(zprev)/df(zprev)
        c += 1
        # print('Iterating...')
    if abs(z-zprev)>tol:
        raise Exception('Could not solve :(')
    return z

# constants
kB = 8.617e-5 #eV/K
q = 1.6e-19 #C
mo = 9.11e-31 #kg
hred = 1.054571817e-34
h = hred*2*pi
epso = 8.85e-12

def ITO_analytical_model(inputs, params):
    # params
    mu_band = params['mu_band'] 
    Ntraps = params['Ntraps']
    Ttraps = params['Ttraps']
    Nch = params['Nch']
    T = params['T']
    tITO = params['tITO']
    tDE = params['tDE']
    me = params['me_factor']*mo
    kDE = params['kDE']
    kITO = params['kITO']
    phiM = params['phiM']
    chiS = params['chiS']
    L = params['L']
    W = params['W']
    # Rs = params['Rs']
    # Rd = params['Rd']
    Ndeep = params['Ndeep']
    Tdeep = params['Tdeep']

    # inputs
    Vgs = inputs['Vgs']
    Vds = inputs['Vds']
    
    Cox = epso*kDE/tDE
    g2D = (me*q)/(pi*hred**2)
    Etraps = kB*Ttraps
    Edeep = kB*Tdeep
    En = lambda n: (hred*pi*n)**2/(2*q*me*tITO**2);
    gtraps = lambda E: (Ntraps/Etraps)*exp((E-En(1))/Etraps)*(E<En(1))
    gdeep = lambda E: (Ndeep/Edeep)*exp((E-En(1))/Edeep)*(E<En(1))
    gfree = lambda E: g2D*((E>=En(1))+(E>=En(2))+(E>=En(3)))

    n2D = lambda delE: g2D*kB*T*log(1+exp(-(delE)/(kB*T)))
    dn2D = lambda delE: g2D*exp(-(delE)/(kB*T))/(1+exp(-(delE)/(kB*T)))
    n2D_approx = lambda delE: g2D*kB*T*exp(-(delE)/(kB*T))
    nfree = lambda Ef,phi: n2D(En(1)-phi-Ef)+n2D(En(2)-phi-Ef)+n2D(En(3)-phi-Ef) #phi is neg of energy
    dnfree = lambda Ef,phi: dn2D(En(1)-phi-Ef)+dn2D(En(2)-phi-Ef)+dn2D(En(3)-phi-Ef) #phi is neg of energy
    nfree_approx = lambda Ef,phi: n2D_approx(En(1)-phi-Ef)+n2D_approx(En(2)-phi-Ef)+n2D_approx(En(3)-phi-Ef) #phi is neg of energy

    xnplus1 = lambda x,n: x*hyp2f1(1,1/n,1+1/n,-(x**n))
    alphac = lambda Ef,phi: exp((En(1)-phi-Ef)/(kB*Ttraps))
    prefac = lambda Ef,phi: (Ntraps)*((alphac(Ef,phi))**(-1))
    ntraps = lambda Ef,phi: prefac(Ef,phi)*(xnplus1(alphac(Ef,phi),Ttraps/T)-xnplus1(0,Ttraps/T))
    dntraps = lambda Ef,phi: (Ntraps/(kB*T))*exp((En(1)-phi-Ef)/(kB*T))*hyp2f1(2,1+T/Ttraps,2+T/Ttraps,\
        -exp((En(1)-phi-Ef)/(kB*T))/(kB*T*(Ttraps/T)*(1+T/Ttraps)))
    ntraps_approx = lambda Ef,phi: prefac(Ef,phi)*(pi*(T/Ttraps)/sin(pi*(T/Ttraps)))

    d_xnplus1 = lambda x,n: x*hyp2f1(1,1/n,1+1/n,-(x**n))
    d_alphac = lambda Ef,phi: exp((En(1)-phi-Ef)/(kB*Tdeep))
    d_prefac = lambda Ef,phi: (Ndeep)*((d_alphac(Ef,phi))**(-1))
    ndeep = lambda Ef,phi: d_prefac(Ef,phi)*(d_xnplus1(d_alphac(Ef,phi),Tdeep/T)-d_xnplus1(0,Tdeep/T))
    dndeep = lambda Ef,phi: (Ndeep/(kB*T))*exp((En(1)-phi-Ef)/(kB*T))*hyp2f1(2,1+T/Tdeep,2+T/Tdeep,\
        -exp((En(1)-phi-Ef)/(kB*T))/(kB*T*(Tdeep/T)*(1+T/Tdeep)))

    Efermi = fsolve(lambda Ef: (nfree(Ef,0)+ntraps(Ef,0)+ndeep(Ef,0)-Nch)/Nch, 0)
    phiMS = phiM - (chiS-Efermi)

    def IdVg(Vg, Vds):
        # phi_fn = lambda V: fsolve(lambda phi: Cox*(Vg-phiMS-phi)\
        #      - q*(nfree(Efermi,phi-V)+ntraps(Efermi,phi-V)+ndeep(Efermi,phi-V)-Nch), 0);
        f = lambda phi, V: Cox*(Vg-phiMS-phi)\
            - q*(nfree(Efermi,phi-V)+ntraps(Efermi,phi-V)+ndeep(Efermi,phi-V)-Nch)
        df = lambda phi, V: -Cox \
            - q*(dnfree(Efermi,phi-V)+dntraps(Efermi,phi-V)+dndeep(Efermi,phi-V))
        phi_fn = lambda V: NR_fsolve(lambda phi: f(phi,V), lambda phi: df(phi,V), 0)
        Id_return = integrate.quad(lambda V: (W/L)*q*mu_band*nfree(Efermi,phi_fn(V)-V),0,Vds)
        # print('error for ',Vg,' is ',Id_return[1])
        return phi_fn(0), Id_return[0]

    phis, Id = IdVg(Vgs,Vds)

    return  phis, Id


# params = dict(L=1, W=1, mu_band=33.5*1e-4, Ntraps=5.5e16, Ttraps=400, Nch=1e17, T=300, tITO=4.5e-9, tDE=5.3e-9,\
#     me_factor=0.3, kDE=16, kITO=9, phiM=5.2, chiS=4.3, Rs=10, Rd=10, Ndeep=5.5e16, Tdeep = 1000)
# inputs = dict(Vgs=0.5, Vds=0.1)
# print(ITO_analytical_model(inputs=inputs, params=params)) 