'''
Plots BB curve for silicon
'''

import ROOT
import math
import os,sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib as mpl
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"],
    "axes.titlesize":"x-large"})
## for Palatino and other serif fonts use:
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
    "axes.titlesize":"x-large"})
## for Palatino and other serif fonts use:
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
    "axes.titlesize":"large"
})

# CONSTANTS!
m_e = 9.109 * pow(10,-31) # kg  SI
c   = 3 * pow(10,8) # m/s
z   = 1 # charge of incoming particle
e   = 1.602 * pow(10,-19) # C # SI
NA  = 6.022 * pow(10,23) # mol
Z   = 14 # atomic number of Si
rho = 2330 # kg / m^3 
A   = 28 * pow(10,-3) # atomic mass of silicon (kg mol-1)
n_e = Z * rho * NA/ A #0.498 * 2.33 * NA / (1.66* pow(10,-24))* pow(10,6)
epsilon_0 = 8.854 * pow(10,-12.) # F m^-1
I = 2.77 * pow(10,-17) # Joules  mean excitation potential
re = 2.8 * pow(10,-15) # m
hbar = 6.626 * pow(10,-34)/2/math.pi
alpha = e**2 / (4 * math.pi * epsilon_0 * hbar * c)

hnup = 4.9764 * pow(10,-18) # Joules, h * nu_p
S0   = 0.201
S1   = 2.872
a    = 0.149
md   = 3.255
delta0 = 0.14
M = 200 * m_e
C = -2 * math.log(I/hnup) - 1

def dedx(beta, flag_density=True, flag_beta = True):
    
    prefactor = 2 * math.pi * n_e *z**2 * e**4
    prefactor = prefactor / (m_e * c**2) / beta**2 
    prefactor = prefactor /pow(4 * math.pi * epsilon_0,2)
    
    Wmax = 2 * m_e *c**2 * bg(beta)**2 / (1 + 2*gamma(beta)*m_e/M + (m_e/M)**2)
    num = 2 * m_e * beta**2 * c**2 * Wmax
    den = I**2 * (1-beta**2)
    dedx = prefactor * (math.log(num/den))
    if flag_density:
        dedx += prefactor * (- density_stern(beta))
    if flag_beta:
        dedx += prefactor * (-2*beta**2)
    dedx = dedx * .01 # convert to cm 
    dedx = dedx * 6.24 * pow(10,12) # convert to MeV
    return dedx 

def gamma(beta):
    return 1./math.sqrt(1-beta**2)

def bg(beta):
    return beta*gamma(beta)

def density(beta):
    n = NA*rho*Z/A
    pe = math.sqrt(4*math.pi*n*pow(re,3)) * m_e * c**2 / alpha
    return 2*(math.log(pe/I) + math.log(bg(beta))-0.5)

def density_stern(beta):
    '''Taken from Leroy, Radiative Interactions in Matter'''
    if bg(beta) < pow(10,S0):
        delta = delta0 * pow(bg(beta)/pow(10,S0),2)
    elif bg(beta) < pow(10,S1):
        delta = 2 * math.log(bg(beta)) + C + a * pow(1/math.log(10) * math.log(pow(10,S1)/bg(beta)),md)
    else:
        delta = 2 * math.log(bg(beta)) + C
    return delta
    


def main():

    bpts           = np.linspace(0.2,.99999999,100)
#    bpts           = np.linspace(0.5,.99999999,100)
    xpts           = [bg(bpt) for bpt in bpts]
    ypts           = [dedx(bpt) for bpt in bpts]
    ypts_nodensity = [dedx(bpt,True,False) for bpt in bpts]
    ypts_class = [dedx(bpt,False,False) for bpt in bpts]
    
    plt.plot(xpts,ypts_class,'m--',linewidth=2, label = "Classical")
    plt.plot(xpts,ypts_nodensity,'b--',linewidth=2, label = "w. 2$\\beta^2$ term")
    plt.plot(xpts,ypts,'turquoise',linewidth=2, label = "w. $\\delta$ and 2$\\beta^2$ terms")
#    plt.plot(xpts,ypts,linewidth=2, label = "w. $\\delta$ and 2$\\beta^2$ terms")
    ax = plt.gca()
    plt.xscale("log")
    plt.yscale("log")
    plt.grid(True,which="both",linewidth=0.5)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(ScalarFormatter())
    ax.set_xlabel('$\\beta \\gamma$',{'fontname':'Helvetica'})
    ax.set_ylabel('$<$dE/dx$>$ [MeV/cm]',{'fontname':'Helvetica'})
    ax.xaxis.label.set_size(18)
    ax.yaxis.label.set_size(18)
    ax.legend(loc="center right")
    plt.savefig("test.pdf")

if __name__=="__main__":
    main()
