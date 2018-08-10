# coding: utf-8
# Author: Zhongyang Zhang

from const_val import *
# pn_init_v2.m


@function
def pn_init_v2(max_particles=None, dope_type=None, dx=None, dy=None, Ltot=None, nx1=None, ny1=None, cdop=None, ppc=None,
               bk=None, T=None, q=None, h=None, alpha=None, eM=None, Gm=None, Lp=None, A=None, B=None, C=None, emR=None,
               hd=None, *args, **kwargs):
    varargin = pn_init_v2.varargin
    nargin = pn_init_v2.nargin

    # dope_type:  1:Uniform 2:Only S/D 3:Uniform+Contacts 4:S/D+Contacts
    # ppsp:       particles per super particle

    # Create Matrix/Vector for Particles and Valley Tracking
    particles = matlabarray(np.zeros((max_particles, 6)))
    # particles[max_particles, 6] = 0
    valley = matlabarray(np.zeros((max_particles, 1)))
    # valley[max_particles, 1] = 0
    valley[:, 1] = 9

    N = 0
    p_x = Lp / dx + 1
    # Create Doping Profile
    bg_charge = matlabarray(np.zeros((ny1, nx1)))
    # bg_charge[ny1, nx1] = 0
    if dope_type == 1:
        for i in arange(1, ny1).reshape(-1):
            for j in arange(1, nx1).reshape(-1):
                bg_charge[i, j] = cdop
                if j < p_x:
                    bg_charge[i, j] = dot(- 1, cdop)

                    # Calculate cpsp (charge per super particle) based on requested ppc (particles
                    # per cell) for the highest doped cell
    cpsp = (dot(dot(cdop, dx), dy)) / ppc
    xmax = dot((nx1 - 1), dx)
    ymax = dot((ny1 - 1), dy)
    n = 0
    for i in arange(1, ny1).reshape(-1):
        for j in arange(1, nx1).reshape(-1):
            npij = floor(dot(dot(abs(bg_charge[i, j]), dx), dy) / cpsp + 0.5)
            # half the number of particles
            if i == 1 or i == ny1:
                npij = npij / 2
            if j == 1 or j == nx1:
                npij = dot(npij / 2, hd)
            for m in arange(1, npij).reshape(-1):
                n = n + 1
                if n > max_particles:
                    disp('Max Number of Particles Exceeded!')
                    break
                if j < p_x:
                    iv = 4
                else:
                    iv = 1
                ei = dot(dot(- (dot(bk, T) / q), log(rand())), 1.5)
                if iv == 4:
                    cos_t = 1 - dot(2, rand())
                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                    phi = dot(dot(2, pi), rand())
                    g = ((B / A) ** 2 + dot((C / A) ** 2, (
                        dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2), sin(phi) ** 2)))) ** 0.5
                    ki = sqrt(dot(dot(dot(2, emR), ei), q) / (dot(dot(dot(abs(A), h), h), (1 + g))))
                    kx = dot(dot(ki, sin_t), cos(phi))
                    ky = dot(dot(ki, sin_t), sin(phi))
                    kz = dot(ki, cos_t)
                else:
                    ki = dot(dot((sqrt(dot(2, eM[iv])) / h), sqrt(dot(ei, (1 + dot(alpha[iv], ei))))), sqrt(q))
                    cos_t = 1 - dot(2, rand())
                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                    phi = dot(dot(2, pi), rand())
                    kx = dot(dot(ki, sin_t), cos(phi))
                    ky = dot(dot(ki, sin_t), sin(phi))
                    kz = dot(ki, cos_t)
                particles[n, 1] = kx
                particles[n, 2] = ky
                particles[n, 3] = kz
                particles[n, 4] = - log(rand()) / Gm[iv]
                particles[n, 5] = dot(dx, (rand() + j - 1.5))
                particles[n, 6] = dot(dy, (rand() + i - 1.5))
                particles[n, 7] = 0
                if i == 1:
                    particles[n, 6] = dot(dot(dy, 0.5), rand())
                if j == 1:
                    particles[n, 5] = dot(dot(dx, 0.5), rand())
                if i == ny1:
                    particles[n, 6] = ymax - dot(dot(dy, 0.5), rand())
                if j == nx1:
                    particles[n, 5] = xmax - dot(dot(dx, 0.5), rand())
                valley[n, 1] = iv
                N = N + 1

    # Number of Particles under p contact
    temp1 = particles[:, 5] < dot(0.5, dx)
    temp3 = valley[:, 1] != 9
    psi = logical_and(temp1, temp3)
    p_total = sum(psi)
    psi = multiply(psi, (arange(1, max_particles)).T)
    psi = psi * (psi != 0)
    # --------------------Added 11/1------------------------------------
    # Count particles belonging to each grid point...0.5*dx to left or right of
    # each grid location (except first and last grid pts)
    p_icpg = zeros(ny1, 1)
    for i in arange(1, length(psi)).reshape(-1):
        index = psi(i)
        y = particles[index, 6]
        if y < dot(0.5, dy):
            p_icpg[1, 1] = p_icpg(1, 1) + 1
        else:
            if y >= Ltot - dot(0.5, dx):
                p_icpg[end(), 1] = p_icpg(end(), 1) + 1
            else:
                yi = round(particles[index, 6] / dy) + 1
                p_icpg[yi, 1] = p_icpg(yi, 1) + 1

                # ------------------------------------------------------------------ 
                # Number of Particles near negative contact
    temp1 = particles[:, 5] > (xmax - dot(0.5, dx))
    pdi = logical_and(temp1, temp3)
    n_total = sum(pdi)
    pdi = multiply(pdi, (arange(1, max_particles)).T)
    pdi = pdi * (pdi != 0)
    # --------------------Added 11/1------------------------------------
    n_icpg = zeros(ny1, 1)
    for i in arange(1, length(pdi)).reshape(-1):
        index = pdi(i)
        y = particles[index, 6]
        if y < dot(0.5, dy):
            n_icpg[1, 1] = n_icpg(1, 1) + 1
        else:
            if y >= Ltot - dot(0.5, dx):
                n_icpg[end(), 1] = n_icpg(end(), 1) + 1
            else:
                xi = round(particles[index, 6] / dy) + 1
                n_icpg[xi, 1] = n_icpg(xi, 1) + 1
