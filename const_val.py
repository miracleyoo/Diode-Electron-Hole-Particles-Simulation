# coding: utf-8
# Author: Zhongyang Zhang

from smop.libsmop import *
from numpy import dot as ndot
from copy import copy
from numpy import cos, sin, pi


def dot(A, B):
    if type(A) == float or type(B) == float:
        return A * B
    else:
        return ndot(A, B)


def max(a, b=None):
    if b is None:
        return np.amax(a)
    else:
        return np.max([a, b])


def min(a, b=None):
    if b is None:
        return np.amin(a)
    else:
        return np.min([a, b])

# -------Simulation Input Settings------------------------------------------
tsteps = 20000
dt = 1e-14
# include hole
dope_type = 1
cdop = 1e+23
cdop2 = 2e+20
cdop3 = 0
hd = 3
max_particles = 400000
de = 0.002
T = 300
# --------Simulation Settings-----------------------------------------------
dx = 1e-08
dy = 1e-08
# --------Device Geomrety---------------------------------------------------
Ttot = 8e-08
Ltot = 4e-07
Lp = Ltot / 2
Ln = Ltot / 2
# ---------Constants--------------------------------------------------------
bk = 1.38066e-23
q = 1.60219e-19
h = 1.05459e-34
emR = 9.10953e-31
eps_o = 8.85419e-12
# ---------GaAs Specific Constants------------------------------------------
Eg = 1.424
Egg = 0
Egl = 0.29
eC = matlabarray([Egg, Egl])
emG = dot(0.067, emR)
emL = dot(0.35, emR)
emh = dot(0.62, emR)
eml = dot(0.074, emR)
eM = matlabarray([emG, emL, emh, eml])
alpha_G = dot((1 / Eg), (1 - emG / emR) ** 2)
alpha_L = dot((1 / (Eg + Egl)), (1 - emL / emR) ** 2)
alpha = matlabarray([alpha_G, alpha_L, 0, 0])
eps_stat = dot(12.9, eps_o)
eps_inf = dot(10.92, eps_o)
eps_p = 1 / ((1 / eps_inf) - (1 / eps_stat))
qD = sqrt(dot(dot(q, q), cdop) / (dot(dot(eps_stat, bk), T)))
ni = 1.8e+12
contact_potential = dot(dot(bk, T) / q, log(dot(cdop, cdop) / ni ** 2))
hw0 = 0.03536
hwij = 0.03
hwe = copy(hwij)
# inverse band mass parameters
A = - 7.65
B = - 4.82
C = 7.7
g100 = B / A
g111 = sqrt((B / A) ** 2 + (C / A) ** 2 / 3)