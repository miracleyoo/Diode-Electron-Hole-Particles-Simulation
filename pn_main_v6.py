# coding: utf-8
# Author: Zhongyang Zhang

# Generated with SMOP  0.41-beta

import pickle
from make_GaAs_hole_scatTable_v2 import *
from make_GaAs_scatTable import *
from pn_charge_v2 import *
from pn_init_v2 import *
from pn_poisson_v5 import *
from pn_renew_v6 import *
from pn_scat_v2 import *
from const_val import *

# pn_main_v6.m

Vp_all = zeros(1)
Vn_all = zeros(1)
POISSON_ITER_MAX = 12000
ppc = 8
path_str = 'C:\\Users\\sunz1\\Desktop\\boeing_dropbox\\nick_matlab_code\\hole_scattering\\pn_junction_data\\1'

# ---------Create Scattering Table------------------------------------------
scatGaAs, GmG, GmL = make_GaAs_scatTable(T, 0, de, 2, cdop, nargout=3)

scatGaAs_hole, Gmh, Gml = make_GaAs_hole_scatTable_v2(T, de, 2, cdop, nargout=3)
Gm = matlabarray([GmG, GmL, Gmh, Gml])
Gm_max = max(Gmh, Gml)
# ------------related to configuration--------------------------------------
nx1 = int(round(Ltot / dx) + 1)
nx = nx1 - 1
ny1 = int(round(Ttot / dy) + 1)
ny = ny1 - 1
bottom_pts = arange(2, nx)
left_pts = arange(1, dot(nx1, ny1), nx1)
right_pts = arange(nx1, dot(nx1, ny1), nx1)
top_pts = arange(dot(nx1, ny1) - nx1 + 2, dot(nx1, ny1) - 1)
p_icpg = [0] * length(left_pts)
n_icpg = [0] * length(right_pts)


for a in arange(1, length(Vp_all)).reshape(-1):
    for b in arange(1, length(Vn_all)).reshape(-1):
        Vp = Vp_all[a] - contact_potential / 2
        Vn = Vn_all[b] + contact_potential / 2
        particles, valley, bg_charge, cpsp, N, p_icpg, n_icpg, xmax, ymax = pn_init_v2(max_particles, dope_type, dx, dy,
                                                                                       Ltot, nx1, ny1, cdop, ppc, bk, T,
                                                                                       q, h, alpha, eM, Gm, Lp, A, B, C,
                                                                                       emR, hd, nargout=9)
        # --------Initial Charge/Field Computations---------------------------------
        charge_p, charge_n = pn_charge_v2(particles, valley, nx1, ny1, dx, dy, max_particles, cpsp, nargout=2)
        phi = zeros(ny1, nx1)
        fx, fy, phi, k = pn_poisson_v5(dx, dy, nx1, ny1, eps_stat, q, charge_p, charge_n, bg_charge, phi, Vp, Vn,
                                       nargout=4)
        number = zeros(tsteps, 4)
        p_enter = zeros(tsteps, 1)
        p_exit = zeros(tsteps, 1)
        p_real = zeros(tsteps, 1)
        n_enter = zeros(tsteps, 1)
        n_exit = zeros(tsteps, 1)
        n_real = zeros(tsteps, 1)
        current_n = zeros(tsteps, 1)
        current_p = zeros(tsteps, 1)
        input_ = zeros(ny1, nx1, 2, tsteps)
        ii = 1
        particle_num = zeros(tsteps + 1, 4)
        for ti in arange(1, tsteps).reshape(-1):
            
            p_temp = 0
            n_temp = 0
            t = dot((ti - 1), dt)
            tdt = t + dt
            for n in arange(1, max_particles).reshape(-1):
                if valley(n, 1) != 9:
                    ts = particles[n, 4]
                    t1 = t
                    while ts < tdt:

                        tau = ts - t1
                        iv = valley(n, 1)
                        if iv != 9:
                            if iv == 1 or iv == 2:
                                kx = particles[n, 1]
                                ky = particles[n, 2]
                                kz = particles[n, 3]
                                x = particles[n, 5]
                                y = particles[n, 6]
                                i = min(floor(y / dy) + 1, ny1)
                                j = min(floor(x / dx) + 1, nx1)
                                i = max(i, 1)
                                j = max(j, 1)
                                dkx = dot(dot(- (q / h), fx[i, j]), tau)
                                dky = dot(dot(- (q / h), fy[i, j]), tau)
                                sk = dot(kx, kx) + dot(ky, ky) + dot(kz, kz)
                                gk = dot(dot((dot(h, h) / (dot(2, eM(iv)))), sk), (1 / q))
                                x = x + dot(
                                    dot((h / eM(iv)), (kx + dot(0.5, dkx))) / (sqrt(1 + dot(dot(4, alpha(iv)), gk))),
                                    tau)
                                y = y + dot(
                                    dot((h / eM(iv)), (ky + dot(0.5, dky))) / (sqrt(1 + dot(dot(4, alpha(iv)), gk))),
                                    tau)
                                particles[n, 1] = kx + dkx
                                particles[n, 2] = ky + dky
                            else:
                                kx = particles[n, 1]
                                ky = particles[n, 2]
                                kz = particles[n, 3]
                                x = particles[n, 5]
                                y = particles[n, 6]
                                i = min(floor(y / dy) + 1, ny1)
                                j = min(floor(x / dx) + 1, nx1)
                                i = max(i, 1)
                                j = max(j, 1)
                                dkx = dot(dot((q / h), fx[i, j]), tau)
                                dky = dot(dot((q / h), fy[i, j]), tau)
                                if iv == 3:
                                    kf = sqrt((kx + dkx) ** 2 + (ky + dky) ** 2 + dot(kz, kz))
                                    cos_theta = kz / kf
                                    sin_theta = sqrt(1 - cos_theta ** 2)
                                    sin_phi = ky / kf / sin_theta
                                    cos_phi = (kx + dkx) / kf / sin_theta
                                    g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                        dot(sin_theta ** 2, cos_theta ** 2) + dot(dot(sin_theta ** 4, cos_phi ** 2),
                                                                                  sin_phi ** 2)))) ** 0.5
                                    mh = emR / (dot(abs(A), (1 - g)))
                                    x = x + dot(dot((h / mh), (kx + dot(0.5, dkx))), tau)
                                    y = y + dot(dot((h / mh), ky), tau)
                                    particles[n, 1] = kx + dkx
                                    particles[n, 2] = ky + dky
                                else:
                                    if iv == 4:
                                        kf = sqrt((kx + dkx) ** 2 + (ky + dky) ** 2 + dot(kz, kz))
                                        cos_theta = kz / kf
                                        sin_theta = sqrt(1 - cos_theta ** 2)
                                        sin_phi = ky / kf / sin_theta
                                        cos_phi = (kx + dkx) / kf / sin_theta
                                        g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                            dot(sin_theta ** 2, cos_theta ** 2) + dot(dot(sin_theta ** 4, cos_phi ** 2),
                                                                                      sin_phi ** 2)))) ** 0.5
                                        ml = emR / (dot(abs(A), (1 + g)))
                                        ef = dot(dot(dot((dot(dot(h, h), abs(A)) / (dot(2, emR))), kf ** 2), (1 / q)),
                                                 (1 + g))
                                        x = x + dot(dot((h / ml), (kx + dot(0.5, dkx))), tau)
                                        y = y + dot(dot((h / ml), ky), tau)
                                        particles[n, 1] = kx + dkx
                                        particles[n, 2] = ky + dky
                            # Boundary Condition-----the former change is incorrect, only kx or ky one has to change 921----------------
                            if x < 0:
                                valley[n, 1] = 9
                                if iv == 1 or iv == 2:
                                    p_temp = p_temp - 1
                                else:
                                    p_temp = p_temp + 1
                            else:
                                if x > xmax:
                                    valley[n, 1] = 9
                                    if iv == 1 or iv == 2:
                                        n_temp = n_temp + 1
                                    else:
                                        n_temp = n_temp - 1
                            if y > ymax:
                                y = ymax - (y - ymax)
                                particles[n, 2] = - particles[n, 2]
                            else:
                                if y < 0:
                                    y = - y
                                    particles[n, 2] = - particles[n, 2]
                            particles[n, 5] = x
                            particles[n, 6] = y
                            # Scatter----------------------------
                            if valley[n, 1] != 9:
                                particle, valley[n, 1] = pn_scat_v2(particles[n, :], valley(n, 1), scatGaAs,
                                                                    scatGaAs_hole, de, q, h, eM, alpha, qD, hw0, A, B,
                                                                    C, emR, n, hwij, Egl, Egg, hwe, g100, g111,
                                                                    nargout=2)
                                particles[n, :] = particle[1, :]
                            t1 = ts
                            ts = t1 - log(rand()) / Gm(iv)
                        else:
                            ts = tdt

                    tau = tdt - t1
                    iv = valley(n, 1)
                    if iv != 9:
                        
                        if iv == 1 or iv == 2:
                            kx = particles[n, 1]
                            ky = particles[n, 2]
                            kz = particles[n, 3]
                            x = particles[n, 5]
                            y = particles[n, 6]
                            i = min(floor(y / dy) + 1, ny1)
                            j = min(floor(x / dx) + 1, nx1)
                            i = max(i, 1)
                            j = max(j, 1)
                            dkx = dot(dot(- (q / h), fx[i, j]), tau)
                            dky = dot(dot(- (q / h), fy[i, j]), tau)
                            sk = dot(kx, kx) + dot(ky, ky) + dot(kz, kz)
                            gk = dot(dot((dot(h, h) / (dot(2, eM(iv)))), sk), (1 / q))
                            x = x + dot(
                                dot((h / eM(iv)), (kx + dot(0.5, dkx))) / (sqrt(1 + dot(dot(4, alpha(iv)), gk))), tau)
                            y = y + dot(
                                dot((h / eM(iv)), (ky + dot(0.5, dky))) / (sqrt(1 + dot(dot(4, alpha(iv)), gk))), tau)
                            particles[n, 1] = kx + dkx
                            particles[n, 2] = ky + dky
                        else:
                            t0 = time.time()
                            kx = particles[n, 1]
                            ky = particles[n, 2]
                            kz = particles[n, 3]
                            x = particles[n, 5]
                            y = particles[n, 6]
                            i = min(floor(y / dy) + 1, ny1)
                            j = min(floor(x / dx) + 1, nx1)
                            i = max(i, 1)
                            j = max(j, 1)
                            dkx = dot(dot((q / h), fx[i, j]), tau)
                            dky = dot(dot((q / h), fy[i, j]), tau)
                            if iv == 3:
                                kf = sqrt((kx + dkx) ** 2 + (ky + dky) ** 2 + dot(kz, kz))
                                cos_theta = kz / kf
                                sin_theta = sqrt(1 - cos_theta ** 2)
                                sin_phi = ky / kf / sin_theta
                                cos_phi = (kx + dkx) / kf / sin_theta
                                g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                    dot(sin_theta ** 2, cos_theta ** 2) + dot(dot(sin_theta ** 4, cos_phi ** 2),
                                                                              sin_phi ** 2)))) ** 0.5
                                mh = emR / (dot(abs(A), (1 - g)))
                                x = x + dot(dot((h / mh), (kx + dot(0.5, dkx))), tau)
                                y = y + dot(dot((h / mh), ky), tau)
                                particles[n, 1] = kx + dkx
                                particles[n, 2] = ky + dky
                            else:
                                if iv == 4:
                                    kf = sqrt((kx + dkx) ** 2 + (ky + dky) ** 2 + dot(kz, kz))
                                    cos_theta = kz / kf
                                    sin_theta = sqrt(1 - cos_theta ** 2)
                                    sin_phi = ky / kf / sin_theta
                                    cos_phi = (kx + dkx) / kf / sin_theta
                                    g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                        dot(sin_theta ** 2, cos_theta ** 2) + dot(dot(sin_theta ** 4, cos_phi ** 2),
                                                                                  sin_phi ** 2)))) ** 0.5
                                    ml = emR / (dot(abs(A), (1 + g)))
                                    ef = dot(dot(dot((dot(dot(h, h), abs(A)) / (dot(2, emR))), kf ** 2), (1 / q)),
                                             (1 + g))
                                    x = x + dot(dot((h / ml), (kx + dot(0.5, dkx))), tau)
                                    y = y + dot(dot((h / ml), ky), tau)
                                    particles[n, 1] = kx + dkx
                                    particles[n, 2] = ky + dky
                        # Boundary Condition-----the former change is incorrect, only kx or ky one has to change
                        # 921----------------
                        if x < 0:
                            valley[n, 1] = 9
                            if iv == 1 or iv == 2:
                                p_temp = p_temp - 1
                            else:
                                p_temp = p_temp + 1
                        else:
                            if x > xmax:
                                valley[n, 1] = 9
                                if iv == 1 or iv == 2:
                                    n_temp = n_temp + 1
                                else:
                                    n_temp = n_temp - 1
                        if y > ymax:
                            y = ymax - (y - ymax)
                            particles[n, 2] = - particles[n, 2]
                        else:
                            if y < 0:
                                y = - y
                                particles[n, 2] = - particles[n, 2]
                        particles[n, 5] = x
                        particles[n, 6] = y
                        particles[n, 4] = ts
            # Renew--------------------
            particles, valley, p_added, n_added, number = pn_renew_v6(particles, valley, Ttot, dx, dy, nx1, ny1,
                                                                      max_particles, p_icpg, n_icpg, bk, T, q, h, alpha,
                                                                      eM, emR, Gm, tdt, left_pts, right_pts, Ltot, A, B,
                                                                      C, ti, number, hd, nargout=5)
            p_real[ti, 1] = p_added - p_temp
            n_real[ti, 1] = n_added - n_temp
            # p_real is how many positive particles are injected in ti
            # n_real is how many negative particles are injected in ti
            # Charge Computation-------
            charge_p, charge_n = pn_charge_v2(particles, valley, nx1, ny1, dx, dy, max_particles, cpsp, nargout=2)
            
            fx, fy, phi, k = pn_poisson_v5(dx, dy, nx1, ny1, eps_stat, q, charge_p, charge_n, bg_charge, phi, Vp, Vn,
                                           nargout=4)
            if ti == 1:
                current_n[ti, 1] = current_n[ti, 1]
            else:
                current_n[ti, 1] = current_n[ti - 1, 1] + n_real[ti, 1]
            # net anode current----
            if ti == 1:
                current_p[ti, 1] = current_p[ti, 1]
            else:
                current_p[ti, 1] = current_p[ti - 1, 1] + p_real[ti, 1]
            index_v1 = (valley[:, 1] == 1)
            index_v1 = multiply(index_v1[:, 1], (arange(1, max_particles)).T)
            index_v1 = index_v1 * (index_v1 != 0)
            index_v2 = (valley[:, 1] == 2)
            index_v2 = multiply(index_v2[:, 1], (arange(1, max_particles)).T)
            index_v2 = index_v2 * (index_v2 != 0)
            index_v3 = (valley[:, 1] == 3)
            index_v3 = multiply(index_v3[:, 1], (arange(1, max_particles)).T)
            index_v3 = index_v3 * (index_v3 != 0)
            index_v4 = (valley[:, 1] == 4)
            index_v4 = multiply(index_v4[:, 1], (arange(1, max_particles)).T)
            index_v4 = index_v4 * (index_v4 != 0)
            particle_num[ti, :] = matlabarray([length(index_v1), length(index_v2), length(index_v3), length(index_v4)])
            jj = 'progress'
            fprintf('%s %f\\n', jj, ti)
        aa = dot(Vp_all[a], 100)
        bb = dot(Vn_all[b], 100)
        mtemp2 = 'Vp=' + str(aa) + 'Vn=' + str(bb) + 'dx=' + str(dot(dx, 1000000000.0)) + 'nm'
        pickle.dump(mtemp2, open('results.pkl', 'wb+'))
