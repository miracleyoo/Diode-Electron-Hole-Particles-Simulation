# coding: utf-8
# Author: Zhongyang Zhang

from smop.libsmop import *
from numpy import cos, sin, pi
from numpy import dot as ndot


def dot(A, B):
    if type(A) == float or type(B) == float:
        return A * B
    else:
        return ndot(A, B)


# pn_renew_v6.m

# function [particles,valley,s_enter,s_exit,s_real,d_enter,d_exit,d_real]
# =gaas_mesfet_renew_v6(particles,valley,Ls,Ld,Ltot,dx,dy,nx1,ny1,max_particles,
# bk,T,q,h,alpha,eM,GmG,tdt,s_pts,d_pts,bgc_tr,ppsp,s_enter,s_exit,s_real,d_enter,
# d_exit,d_real,ti,source_temp,drain_temp,s_icpg)

@function
def pn_renew_v6(particles=None, valley=None, Ttot=None, dx=None, dy=None, nx1=None, ny1=None, max_particles=None,
                p_icpg=None, n_icpg=None, bk=None, T=None, q=None, h=None, alpha=None, eM=None, emR=None, Gm=None,
                tdt=None, left_pts=None, right_pts=None, Ltot=None, A=None, B=None, C=None, ti=None, number=None,
                hd=None, *args, **kwargs):
    varargin = pn_renew_v6.varargin
    nargin = pn_renew_v6.nargin

    xmax = dot((nx1 - 1), dx)
    ymax = dot((ny1 - 1), dy)
    # --------------how many particles are under simulation--------------------
    temp1 = valley[:, 1] == 1
    temp2 = valley[:, 1] == 2
    temp3 = valley[:, 1] == 3
    temp4 = valley[:, 1] == 4
    temp1 = multiply(temp1, (arange(1, max_particles)).T)
    temp2 = multiply(temp2, (arange(1, max_particles)).T)
    temp3 = multiply(temp3, (arange(1, max_particles)).T)
    temp4 = multiply(temp4, (arange(1, max_particles)).T)
    temp1 = temp1*(temp1 != 0)
    temp2 = temp2*(temp2 != 0)
    temp3 = temp3*(temp3 != 0)
    temp4 = temp4*(temp4 != 0)
    number[ti, 1] = length(temp1)
    number[ti, 2] = length(temp2)
    number[ti, 3] = length(temp3)
    number[ti, 4] = length(temp4)
    # --------------positive contact Calculations-----------------------------------------
    # Find Particles In Half Cells at positive contact
    temp1 = particles[:, 5] < dot(0.5, dx)
    temp2 = valley[:, 1] == 3
    temp3 = valley[:, 1] == 4

    temp4 = logical_and(temp1, temp2)
    temp5 = logical_and(temp1, temp3)
    temp6 = valley[:, 1] != 9
    ph = logical_or(temp4, temp5)

    ph = multiply(ph, (arange(1, max_particles)).T)
    ph = ph*(ph != 0)
    # temp1=particles(:,5)<0.5*dx;
    temp2 = valley[:, 1] == 1
    temp3 = valley[:, 1] == 2

    temp4 = logical_and(temp1, temp2)
    temp5 = logical_and(temp1, temp3)
    # temp6=valley(:,1)~=9;

    pei = logical_or(temp4, temp5)

    pei = multiply(pei, (arange(1, max_particles)).T)
    pei = pei*(pei != 0)
    # Count particles belonging to each grid point...0.5*dx to left or right of
    # each grid location (except first and last grid pts)
    p_cpg = matlabarray(np.zeros((length(left_pts), 1)))
    for i in arange(1, length(ph)).reshape(-1):
        index = ph(i)
        y = particles(index, 6)
        if y < dot(0.5, dy):
            p_cpg[1, 1] = p_cpg(1, 1) + 1
        else:
            if y >= Ttot - dot(0.5, dy):
                p_cpg[end(), 1] = p_cpg(end(), 1) + 1
            else:
                yi = round(particles(index, 6) / dy) + 1
                p_cpg[yi, 1] = p_cpg(yi, 1) + 1

    for i in arange(1, length(pei)).reshape(-1):
        index = pei(i)
        y = particles(index, 6)
        if y < dot(0.5, dy):
            p_cpg[1, 1] = p_cpg(1, 1) - 1
        else:
            if y >= Ttot - dot(0.5, dy):
                p_cpg[end(), 1] = p_cpg(end(), 1) - 1
            else:
                yi = round(particles(index, 6) / dy) + 1
                p_cpg[yi, 1] = p_cpg(yi, 1) - 1

    p_delta = p_cpg - p_icpg

    p_added = dot(- 1, sum(p_delta))

    # Indices of Particles to Be Deleted
    v9i = logical_not(temp6)
    v9i = multiply(v9i, (arange(1, max_particles)).T)
    v9i = v9i*(v9i != 0)
    v9_count = 1
    for i in arange(1, length(p_delta)).reshape(-1):
        p_dif = abs(p_delta(i))
        if p_delta(i) < 0:
            for j in arange(1, p_dif).reshape(-1):
                index = v9i(v9_count)
                iv = 4
                ei = dot(dot(- (dot(bk, T) / q), log(rand())), 1.5)
                cos_t = 1 - dot(2, rand())
                sin_t = sqrt(1 - dot(cos_t, cos_t))
                phi = dot(dot(2, pi), rand)
                g = ((B / A) ** 2 + dot((C / A) ** 2, (
                    dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2), sin(phi) ** 2)))) ** 0.5
                ki = sqrt(dot(dot(dot(2, emR), ei), q) / (dot(dot(dot(abs(A), h), h), (1 + g))))
                kx = dot(dot(ki, sin_t), cos(phi))
                ky = dot(dot(ki, sin_t), sin(phi))
                kz = dot(ki, cos_t)
                if kx < 0:
                    kx = - kx
                particles[index, 1] = kx
                particles[index, 2] = ky
                particles[index, 3] = kz
                particles[index, 4] = tdt - log(rand()) / Gm(iv)
                if i == 1:
                    particles[index, 6] = dot(dot(0.5, dy), rand())
                else:
                    if i == length(p_delta):
                        particles[index, 6] = Ttot - dot(dot(0.5, dy), rand())
                    else:
                        particles[index, 6] = dot((i - 1.5), dy) + dot(dy, rand())
                particles[index, 5] = dot(dot(0.5, dx), rand())
                valley[index, 1] = iv
                v9_count = v9_count + 1
                # Need to Delete Particles
        else:
            if p_delta(i) > 0:
                # Get indices of particles belonging to current grid pt
                # If first grid pt, just get particles within first half cell
                if i == 1:
                    psi_temp = particles(ph, 6) < dot(0.5, dy)
                    psi_temp = multiply(psi_temp, (arange(1, length(ph))).T)
                    psi_temp = psi_temp*(psi_temp != 0)
                else:
                    if i == length(p_delta):
                        psi_temp = particles(ph, 6) > (Ttot - dot(0.5, dy))
                        psi_temp = multiply(psi_temp, (arange(1, length(ph))).T)
                        psi_temp = psi_temp*(psi_temp != 0)
                    else:
                        psi_temp1 = particles(ph, 6) > (dot((i - 1.5), dx))
                        psi_temp2 = particles(ph, 6) < (dot((i - 0.5), dx))
                        psi_temp = logical_and(psi_temp1, psi_temp2)
                        psi_temp = multiply(psi_temp, (arange(1, length(ph))).T)
                        psi_temp = psi_temp*(psi_temp != 0)
                # Go through and randomly delete particles
                for j in arange(1, p_dif).reshape(-1):
                    temp_index = 1 + floor(dot(length(psi_temp), rand()))
                    index = psi_temp(temp_index)
                    index2 = ph(index)
                    valley[index2] = 9
                    psi_temp[temp_index] = []

                    # ------------negative  Calculations--------------------------------------------
                    # Find Particles In Half Cells Under negative Grid Pts
    temp1 = particles[:, 5] > Ltot - dot(0.5, dx)
    temp2 = valley[:, 1] == 1
    temp3 = valley[:, 1] == 2
    temp4 = logical_and(temp1, temp2)
    temp5 = logical_and(temp1, temp3)
    temp6 = valley[:, 1] != 9
    nei = logical_or(temp4, temp5)
    nei = multiply(nei, (arange(1, max_particles)).T)
    nei = nei*(nei != 0)
    # temp1=particles(:,5)>Ltot-0.5*dx;
    temp2 = valley[:, 1] == 3
    temp3 = valley[:, 1] == 4
    temp4 = logical_and(temp1, temp2)
    temp5 = logical_and(temp1, temp3)
    temp6 = valley[:, 1] != 9
    npi = logical_or(temp4, temp5)
    npi = multiply(npi, (arange(1, max_particles)).T)
    npi = npi*(npi != 0)
    # d_cpg: drain count per grid
    n_cpg = matlabarray(np.zeros((length(right_pts), 1)))
    # n_cpg[length(right_pts), 1] = 0
    for i in arange(1, length(nei)).reshape(-1):
        index = nei(i)
        y = particles(index, 6)
        if y < dot(0.5, dy):
            n_cpg[1, 1] = n_cpg(1, 1) + 1
        else:
            if y >= Ttot - dot(0.5, dx):
                n_cpg[end(), 1] = n_cpg(end(), 1) + 1
            else:
                yi = round(particles(index, 6) / dy) + 1
                n_cpg[yi, 1] = n_cpg(yi, 1) + 1

    for i in arange(1, length(npi)).reshape(-1):
        index = npi(i)
        y = particles(index, 6)
        if y < dot(0.5, dy):
            n_cpg[1, 1] = n_cpg(1, 1) - 1
        else:
            if y >= Ttot - dot(0.5, dx):
                n_cpg[end(), 1] = n_cpg(end(), 1) - 1
            else:
                yi = round(particles(index, 6) / dy) + 1
                n_cpg[yi, 1] = n_cpg(yi, 1) - 1

    n_delta = n_cpg - n_icpg
    n_added = dot(- 1, sum(n_delta))

    # Indices of Particles to Be Deleted
    v9i = logical_not(temp6)
    v9i = multiply(v9i, (arange(1, max_particles)).T)
    v9i = v9i*(v9i != 0)
    v9_count = 1
    for i in arange(1, length(n_delta)).reshape(-1):
        n_dif = abs(n_delta(i))
        if n_delta(i) < 0:
            for j in arange(1, n_dif).reshape(-1):
                index = v9i(v9_count)
                iv = 1
                ei = dot(dot(- (dot(bk, T) / q), log(rand())), 1.5)
                ki = dot(dot((sqrt(dot(2, eM(iv))) / h), sqrt(dot(ei, (1 + dot(alpha(iv), ei))))), sqrt(q))
                cos_t = 1 - dot(2, rand())
                sin_t = sqrt(1 - dot(cos_t, cos_t))
                phi = dot(dot(2, pi), rand)
                kx = dot(dot(ki, sin_t), cos(phi))
                ky = dot(dot(ki, sin_t), sin(phi))
                kz = dot(ki, cos_t)
                if kx > 0:
                    kx = - kx
                particles[index, 1] = kx
                particles[index, 2] = ky
                particles[index, 3] = kz
                particles[index, 4] = tdt - log(rand()) / Gm(iv)
                # -----------------------
                if i == 1:
                    particles[index, 6] = dot(dot(0.5, dy), rand())
                else:
                    if i == length(n_delta):
                        particles[index, 6] = Ttot - dot(dot(0.5, dy), rand())
                    else:
                        particles[index, 6] = dot((i - 1.5), dy) + dot(dy, rand())
                particles[index, 5] = Ltot - dot(dot(0.5, dx), rand())
                valley[index, 1] = iv
                v9_count = v9_count + 1
                # Need to Delete Particles
        else:
            if n_delta(i) > 0:
                # Get indices of particles belonging to current grid pt
                # If first grid pt, just get particles within first half cell
                if i == 1:
                    pdi_temp = particles(nei, 6) < dot(0.5, dy)
                    pdi_temp = multiply(pdi_temp, (arange(1, length(nei))).T)
                    pdi_temp = pdi_temp*(pdi_temp != 0)
                else:
                    if i == length(n_delta):
                        pdi_temp = particles(nei, 5) > (Ttot - dot(0.5, dy))
                        pdi_temp = multiply(pdi_temp, (arange(1, length(nei))).T)
                        pdi_temp = pdi_temp*(pdi_temp != 0)
                    else:
                        pdi_temp1 = particles(nei, 6) > (dot((i - 1.5), dx))
                        pdi_temp2 = particles(nei, 6) < (dot((i - 0.5), dx))
                        pdi_temp = logical_and(pdi_temp1, pdi_temp2)
                        pdi_temp = multiply(pdi_temp, (arange(1, length(nei))).T)
                        pdi_temp = pdi_temp*(pdi_temp != 0)
                # Go through and randomly delete particles
                for j in arange(1, n_dif).reshape(-1):
                    temp_index = ceil(dot(length(pdi_temp), rand()))
                    index = pdi_temp(temp_index)
                    index2 = nei(index)
                    valley[index2] = 9
                    pdi_temp[temp_index] = []
