# coding: utf-8
# Author: Zhongyang Zhang

from const_val import *


# make_GaAs_hole_scatTable_v2.m


@function
def make_GaAs_hole_scatTable_v2(T=None, de=None, Vmax=None, cdop=None, *args, **kwargs):
    varargin = make_GaAs_hole_scatTable_v2.varargin
    nargin = make_GaAs_hole_scatTable_v2.nargin

    # --------Electron Energy Steps for Formulas/Graphs-------------------------
    delt_Ek = de
    Ek_pts = int(Vmax / de)
    eV_axis = dot((arange(1, Ek_pts)), delt_Ek)
    scat_h = matlabarray(np.zeros((Ek_pts, 12)))
    # scat_h[Ek_pts, 12] = 0
    scat_l = matlabarray(np.zeros((Ek_pts, 12)))
    # scat_l[Ek_pts, 12] = 0

    # ---------General Constants------------------------------------------------
    kb = 1.38066e-23
    hwo = 0.03536
    hwoq = dot(hwo, q)
    eps_s = 12.9
    Ni = cdop
    m = emR
    ml = dot(0.082, m)
    mh = dot(0.45, m)
    adp = dot(7, q)
    rho = 5360
    s = 3860

    # inverse band mass parameters
    A = - 6.98
    B = - 4.5
    C = 6.2
    # -------acoustic phonon scattering-----------------------------------------
    for i in arange(1, Ek_pts).reshape(-1):
        ei = dot(delt_Ek, i)
        ef = ei
        gamma = dot(dot(dot(adp, adp), kb), T) / (dot(dot(dot(dot((dot(2, pi)) ** 2, h), rho), s), s))
        a = dot(dot(h, h), abs(A)) / (dot(2, emR))
        g = sqrt((B / A) ** 2 + (C / A) ** 2 / 3)
        # P=gamma*ei^0.5*0.25*(1+3*cos(k_theta)^2)/(a*(1-g))^1.5;
        P = dot(dot(dot(gamma, ei ** 0.5), 0.5), 1) / (dot(a, (1 - g))) ** 1.5
        scat_h[i, 11] = P / 1000000000.0
        scat_l[i, 12] = P / 1000000000.0
        g = B / A
        P = dot(dot(dot(gamma, ei ** 0.5), 0.5), 1) / (dot(a, (1 + g))) ** 1.5
        scat_l[i, 11] = P / 1000000000.0
        scat_h[i, 12] = P / 1000000000.0

    # nonoplar optical phonon scattering
    DK2 = 1.58e+22

    for i in arange(1, Ek_pts).reshape(-1):
        # final is heavy,abs
        ei = dot(delt_Ek, i)
        ef = ei + hw0
        N0 = 1 / (exp((dot(q, hw0)) / (dot(kb, T))) - 1)
        B0 = dot(h ** 2, DK2) / (dot(dot(2, rho), hw0))
        gamma = dot(dot(dot(2, pi) / h, B0), N0)
        a = dot(dot(h, h), abs(A)) / (dot(2, emR))
        g = sqrt((B / A) ** 2 + (C / A) ** 2 / 3)
        P = dot(dot(gamma / (dot(2, pi)) ** 3, ef ** 0.5), 0.5) / (dot(a, (1 - g))) ** 1.5
        scat_h[i, 7] = P / 1e+27
        scat_l[i, 7] = P / 1e+27
        # final is light, absorbsion
        ei = dot(delt_Ek, i)
        ef = ei + hw0
        N0 = 1 / (exp((dot(q, hw0)) / (dot(kb, T))) - 1)
        B0 = dot(h ** 2, DK2) / (dot(dot(2, rho), hw0))
        gamma = dot(dot(dot(2, pi) / h, B0), N0)
        a = dot(dot(h, h), abs(A)) / (dot(2, emR))
        g = B / A
        P = dot(dot(gamma / (dot(2, pi)) ** 3, ef ** 0.5), 0.5) / (dot(a, (1 + g))) ** 1.5
        scat_h[i, 9] = P / 1e+27
        scat_l[i, 9] = P / 1e+27
        # final is heavy, emission
        ei = dot(delt_Ek, i)
        if (ei - hw0) <= de:
            P = 0
            scat_h[i, 8] = P
            scat_l[i, 8] = P
        else:
            ef = ei - hw0
            N0 = 1 / (exp((dot(q, hw0)) / (dot(kb, T))) - 1)
            B0 = dot(h ** 2, DK2) / (dot(dot(2, rho), hw0))
            gamma = dot(dot(dot(2, pi) / h, B0), (N0 + 1))
            a = dot(dot(h, h), abs(A)) / (dot(2, emR))
            g = sqrt((B / A) ** 2 + (C / A) ** 2 / 3)
            P = dot(dot(gamma / (dot(2, pi)) ** 3, ef ** 0.5), 0.5) / (dot(a, (1 - g))) ** 1.5
            scat_h[i, 8] = P / 1e+27
            scat_l[i, 8] = P / 1e+27
        # final is light, emission
        ei = dot(delt_Ek, i)
        if (ei - hw0) <= de:
            P = 0
            scat_h[i, 10] = P
            scat_l[i, 10] = P
        else:
            ef = ei - hw0
            N0 = 1 / (exp((dot(q, hw0)) / (dot(kb, T))) - 1)
            B0 = dot(h ** 2, DK2) / (dot(dot(2, rho), hw0))
            gamma = dot(dot(dot(2, pi) / h, B0), (N0 + 1))
            a = dot(dot(h, h), abs(A)) / (dot(2, emR))
            g = B / A
            P = dot(dot(gamma / (dot(2, pi)) ** 3, ef ** 0.5), 0.5) / (dot(a, (1 + g))) ** 1.5
            scat_h[i, 10] = P / 1e+27
            scat_l[i, 10] = P / 1e+27

    # --------Impurity Scattering-----------------------------------------------
    for i in arange(1, Ek_pts).reshape(-1):
        # heavy holes to light holes---------------------------------
        ki = sqrt(dot(dot(dot(dot(2, mh), delt_Ek), i), q) / (dot(h, h)))
        kf = dot(sqrt(ml / mh), ki)
        beta = sqrt(dot(dot(Ni, q), q) / (dot(dot(dot(kb, T), eps_o), eps_s)))
        F = dot((dot(beta, beta) + dot(ki, ki) + dot(kf, kf)) / (dot(ki, kf)),
                log((dot(beta, beta) + (ki + kf) ** 2) / (dot(beta, beta) + (ki - kf) ** 2))) - 4
        P = dot(dot(dot(dot(3, q ** 4), Ni), ml), F) / (
            dot(dot(dot(dot(dot(dot(32, pi), h ** 3), eps_o ** 2), eps_s ** 2), ki ** 2), kf))
        scat_h[i, 1] = P
        ki = sqrt(dot(dot(dot(dot(2, ml), delt_Ek), i), q) / (dot(h, h)))
        kf = dot(sqrt(mh / ml), ki)
        F = dot((dot(beta, beta) + dot(ki, ki) + dot(kf, kf)) / (dot(ki, kf)),
                log((dot(beta, beta) + (ki + kf) ** 2) / (dot(beta, beta) + (ki - kf) ** 2))) - 4
        P = dot(dot(dot(dot(3, q ** 4), Ni), mh), F) / (
            dot(dot(dot(dot(dot(dot(32, pi), h ** 3), eps_o ** 2), eps_s ** 2), ki ** 2), kf))
        scat_l[i, 1] = P
        ki = sqrt(dot(dot(dot(dot(2, mh), delt_Ek), i), q) / (dot(h, h)))
        kf = ki
        F = dot((dot(beta, beta) + dot(dot(2, ki), ki)) / (dot(ki, ki)),
                log((dot(beta, beta)) / (dot(beta, beta) + dot(dot(4, ki), ki)))) + dot(4 / 3,
                (dot(3, beta ** 4) + dot(dot(12, beta ** 2), ki ** 2) + dot(8, ki ** 4))) / \
                (dot(beta ** 2,(beta ** 2 + dot(4,ki ** 2))))
        P = dot(dot(dot(dot(3, q ** 4), Ni), mh), F) / (
            dot(dot(dot(dot(dot(dot(32, pi), h ** 3), eps_o ** 2), eps_s ** 2), ki ** 2), kf))
        scat_h[i, 2] = P
        ki = sqrt(dot(dot(dot(dot(2, ml), delt_Ek), i), q) / (dot(h, h)))
        kf = ki
        F = dot((dot(beta, beta) + dot(dot(2, ki), ki)) / (dot(ki, ki)),
                log((dot(beta, beta)) / (dot(beta, beta) + dot(dot(4, ki), ki)))) + dot(4 / 3, (
            dot(3, beta ** 4) + dot(dot(12, beta ** 2), ki ** 2) + dot(8, ki ** 4))) / (dot(beta ** 2,
            (beta ** 2 + dot(4,ki ** 2))))
        P = dot(dot(dot(dot(3, q ** 4), Ni), ml), F) / (
            dot(dot(dot(dot(dot(dot(32, pi), h ** 3), eps_o ** 2), eps_s ** 2), ki ** 2), kf))
        scat_l[i, 2] = P

    # optical phonon scattering
    for i in arange(1, Ek_pts).reshape(-1):
        # heavy holes to heavy holes,aborbtion---------------------------------
        ep = 1 / ((1 / (dot(eps_inf, eps_o))) - (1 / (dot(eps_s, eps_o))))
        N0 = 1 / (exp((dot(q, hw0)) / (dot(kb, T))) - 1)
        c_const = dot((q / h) ** 3, (dot(hw0, mh))) / (dot(dot(4, pi), ep))
        ab_const = dot(c_const, N0)
        ei = dot(delt_Ek, i)
        ef = ei + hw0
        ki = sqrt(dot(dot(dot(dot(2, mh), delt_Ek), i), q) / (dot(h, h)))
        kf = sqrt(dot(dot(dot(2, mh), ef), q) / (dot(h, h)))
        fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)))
        theta = (dot(ki, ki) + dot(kf, kf)) / (dot(dot(2, ki), kf))
        hi = (1 + dot(dot(3, theta), (theta - 1 / fai))) / 4
        P = dot(dot(ab_const, fai), hi) / ki
        scat_h[i, 3] = P
        # light holes to light holes,aborbtion---------------------------------
        ep = 1 / ((1 / (dot(eps_inf, eps_o))) - (1 / (dot(eps_s, eps_o))))
        N0 = 1 / (exp((dot(q, hw0)) / (dot(kb, T))) - 1)
        c_const = dot((q / h) ** 3, (dot(hw0, ml))) / (dot(dot(4, pi), ep))
        ab_const = dot(c_const, N0)
        ei = dot(delt_Ek, i)
        ef = ei + hw0
        ki = sqrt(dot(dot(dot(dot(2, ml), delt_Ek), i), q) / (dot(h, h)))
        kf = sqrt(dot(dot(dot(2, ml), ef), q) / (dot(h, h)))
        fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)))
        theta = (dot(ki, ki) + dot(kf, kf)) / (dot(dot(2, ki), kf))
        hi = (1 + dot(dot(3, theta), (theta - 1 / fai))) / 4
        P = dot(dot(ab_const, fai), hi) / ki
        scat_l[i, 3] = P
        # heavy holes to light holes,aborbtion---------------------------------
        ep = 1 / ((1 / (dot(eps_inf, eps_o))) - (1 / (dot(eps_s, eps_o))))
        N0 = 1 / (exp((dot(q, hw0)) / (dot(kb, T))) - 1)
        c_const = dot((q / h) ** 3, (dot(hw0, ml))) / (dot(dot(4, pi), ep))
        ab_const = dot(c_const, N0)
        ei = dot(delt_Ek, i)
        ef = ei + hw0
        ki = sqrt(dot(dot(dot(dot(2, mh), delt_Ek), i), q) / (dot(h, h)))
        kf = sqrt(dot(dot(dot(2, ml), ef), q) / (dot(h, h)))
        fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)))
        theta = (dot(ki, ki) + dot(kf, kf)) / (dot(dot(2, ki), kf))
        hi = dot(3, (1 - dot(theta, (theta - 1 / fai)))) / 4
        P = dot(dot(ab_const, fai), hi) / ki
        scat_h[i, 4] = P
        # light holes to heavy holes,aborbtion---------------------------------
        ep = 1 / ((1 / (dot(eps_inf, eps_o))) - (1 / (dot(eps_s, eps_o))))
        N0 = 1 / (exp((dot(q, hw0)) / (dot(kb, T))) - 1)
        c_const = dot((q / h) ** 3, (dot(hw0, mh))) / (dot(dot(4, pi), ep))
        ab_const = dot(c_const, N0)
        ei = dot(delt_Ek, i)
        ef = ei + hw0
        ki = sqrt(dot(dot(dot(dot(2, ml), delt_Ek), i), q) / (dot(h, h)))
        kf = sqrt(dot(dot(dot(2, mh), ef), q) / (dot(h, h)))
        fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)))
        theta = (dot(ki, ki) + dot(kf, kf)) / (dot(dot(2, ki), kf))
        hi = dot(3, (1 - dot(theta, (theta - 1 / fai)))) / 4
        P = dot(dot(ab_const, fai), hi) / ki
        scat_l[i, 4] = P
        # heavy holes to heavy holes,emission---------------------------------
        ep = 1 / ((1 / (dot(eps_inf, eps_o))) - (1 / (dot(eps_s, eps_o))))
        N0 = 1 / (exp((dot(q, hw0)) / (dot(kb, T))) - 1)
        c_const = dot((q / h) ** 3, (dot(hw0, mh))) / (dot(dot(4, pi), ep))
        em_const = dot(c_const, (N0 + 1))
        ei = dot(delt_Ek, i)
        ef = ei - hw0
        if ef > de:
            ki = sqrt(dot(dot(dot(dot(2, mh), delt_Ek), i), q) / (dot(h, h)))
            kf = sqrt(dot(dot(dot(2, mh), ef), q) / (dot(h, h)))
            fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)))
            theta = (dot(ki, ki) + dot(kf, kf)) / (dot(dot(2, ki), kf))
            hi = (1 + dot(dot(3, theta), (theta - 1 / fai))) / 4
            P = dot(dot(em_const, fai), hi) / ki
            scat_h[i, 5] = P
        else:
            scat_h[i, 5] = 0
        # light holes to light holes,emission---------------------------------
        ep = 1 / ((1 / (dot(eps_inf, eps_o))) - (1 / (dot(eps_s, eps_o))))
        N0 = 1 / (exp((dot(q, hw0)) / (dot(kb, T))) - 1)
        c_const = dot((q / h) ** 3, (dot(hw0, ml))) / (dot(dot(4, pi), ep))
        em_const = dot(c_const, (N0 + 1))
        ei = dot(delt_Ek, i)
        ef = ei - hw0
        if ef > de:
            ki = sqrt(dot(dot(dot(dot(2, ml), delt_Ek), i), q) / (dot(h, h)))
            kf = sqrt(dot(dot(dot(2, ml), ef), q) / (dot(h, h)))
            fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)))
            theta = (dot(ki, ki) + dot(kf, kf)) / (dot(dot(2, ki), kf))
            hi = (1 + dot(dot(3, theta), (theta - 1 / fai))) / 4
            P = dot(dot(em_const, fai), hi) / ki
            scat_l[i, 5] = P
        else:
            scat_l[i, 5] = 0
        # heavy holes to light holes,emission---------------------------------
        ep = 1 / ((1 / (dot(eps_inf, eps_o))) - (1 / (dot(eps_s, eps_o))))
        N0 = 1 / (exp((dot(q, hw0)) / (dot(kb, T))) - 1)
        c_const = dot((q / h) ** 3, (dot(hw0, ml))) / (dot(dot(4, pi), ep))
        em_const = dot(c_const, (N0 + 1))
        ei = dot(delt_Ek, i)
        ef = ei - hw0
        if ef > de:
            ki = sqrt(dot(dot(dot(dot(2, mh), delt_Ek), i), q) / (dot(h, h)))
            kf = sqrt(dot(dot(dot(2, ml), ef), q) / (dot(h, h)))
            fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)))
            theta = (dot(ki, ki) + dot(kf, kf)) / (dot(dot(2, ki), kf))
            hi = dot(3, (1 - dot(theta, (theta - 1 / fai)))) / 4
            P = dot(dot(em_const, fai), hi) / ki
            scat_h[i, 6] = P
        else:
            scat_h[i, 6] = 0
        # light holes to heavy holes,emission---------------------------------
        ep = 1 / ((1 / (dot(eps_inf, eps_o))) - (1 / (dot(eps_s, eps_o))))
        N0 = 1 / (exp((dot(q, hw0)) / (dot(kb, T))) - 1)
        c_const = dot((q / h) ** 3, (dot(hw0, mh))) / (dot(dot(4, pi), ep))
        em_const = dot(c_const, (N0 + 1))
        ei = dot(delt_Ek, i)
        ef = ei - hw0
        if ef > de:
            ki = sqrt(dot(dot(dot(dot(2, ml), delt_Ek), i), q) / (dot(h, h)))
            kf = sqrt(dot(dot(dot(2, mh), ef), q) / (dot(h, h)))
            fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)))
            theta = (dot(ki, ki) + dot(kf, kf)) / (dot(dot(2, ki), kf))
            hi = dot(3, (1 - dot(theta, (theta - 1 / fai)))) / 4
            P = dot(dot(em_const, fai), hi) / ki
            scat_l[i, 6] = P
        else:
            scat_l[i, 6] = 0

    totScath = scat_h[:, 1] + scat_h[:, 2] + scat_h[:, 3] + scat_h[:, 4] + scat_h[:, 5] + scat_h[:, 6] + \
               scat_h[:, 7] + scat_h[:, 8] + scat_h[:, 9] + scat_h[:, 10] + scat_h[:, 11] + scat_h[:, 12]
    totScatl = scat_l[:, 1] + scat_l[:, 2] + scat_l[:, 3] + scat_l[:, 4] + scat_l[:, 5] + scat_l[:, 6] + \
               scat_l[:, 7] + scat_l[:, 8] + scat_l[:, 9] + scat_l[:, 10] + scat_l[:, 11] + scat_l[:, 12]
    Gmh = max(totScath[:, 1])
    Gml = max(totScatl[:, 1])
    scatGaAs_hole = matlabarray(np.zeros((12, Ek_pts, 2)))
    scatGaAs_hole[np.ix_([1], np.arange(1, Ek_pts + 1), [1])] = scat_h[:, 1].T
    for i in range(2, 13):
        scatGaAs_hole[np.ix_([i], np.arange(1, Ek_pts+1), [1])] = \
            scatGaAs_hole[np.ix_([i-1], np.arange(1, Ek_pts+1), [1])] + scat_h[:, i].T
    scatGaAs_hole[np.ix_(np.arange(1, 13), np.arange(1, Ek_pts+1), [1])] /= Gmh

    scatGaAs_hole[np.ix_([1], np.arange(1, Ek_pts + 1), [2])] = scat_l[:, 2].T
    for i in range(2, 13):
        scatGaAs_hole[np.ix_([i], np.arange(1, Ek_pts+1), [2])] = \
            scatGaAs_hole[np.ix_([i-1], np.arange(1, Ek_pts+1), [2])] + scat_l[:, i].T
    scatGaAs_hole[np.ix_(np.arange(1, 13), np.arange(1, Ek_pts+1), [2])] /= Gml
    return scatGaAs_hole, Gmh, Gml
    # scatGaAs_hole[1, :, 1] = scat_h[:, 1].T
    # scatGaAs_hole[2, :, 1] = scatGaAs_hole[1, :, 1] + scat_h[:, 2].T
    # scatGaAs_hole[3, :, 1] = scatGaAs_hole[2, :, 1] + scat_h[:, 3].T
    # scatGaAs_hole[4, :, 1] = scatGaAs_hole[3, :, 1] + scat_h[:, 4].T
    # scatGaAs_hole[5, :, 1] = scatGaAs_hole[4, :, 1] + scat_h[:, 5].T
    # scatGaAs_hole[6, :, 1] = scatGaAs_hole[5, :, 1] + scat_h[:, 6].T
    # scatGaAs_hole[7, :, 1] = scatGaAs_hole[6, :, 1] + scat_h[:, 7].T
    # scatGaAs_hole[8, :, 1] = scatGaAs_hole[7, :, 1] + scat_h[:, 8].T
    # scatGaAs_hole[9, :, 1] = scatGaAs_hole[8, :, 1] + scat_h[:, 9].T
    # scatGaAs_hole[10, :, 1] = scatGaAs_hole[9, :, 1] + scat_h[:, 10].T
    # scatGaAs_hole[11, :, 1] = scatGaAs_hole[10, :, 1] + scat_h[:, 11].T
    # scatGaAs_hole[12, :, 1] = scatGaAs_hole[11, :, 1] + scat_h[:, 12].T
    # scatGaAs_hole[:, :, 1] = scatGaAs_hole[:, :, 1] / Gmh
    # scatGaAs_hole[1, :, 2] = scat_l[:, 2].T
    # scatGaAs_hole[2, :, 2] = scatGaAs_hole[1, :, 2] + scat_l[:, 2].T
    # scatGaAs_hole[3, :, 2] = scatGaAs_hole[2, :, 2] + scat_l[:, 3].T
    # scatGaAs_hole[4, :, 2] = scatGaAs_hole[3, :, 2] + scat_l[:, 4].T
    # scatGaAs_hole[5, :, 2] = scatGaAs_hole[4, :, 2] + scat_l[:, 5].T
    # scatGaAs_hole[6, :, 2] = scatGaAs_hole[5, :, 2] + scat_l[:, 6].T
    # scatGaAs_hole[7, :, 2] = scatGaAs_hole[6, :, 2] + scat_l[:, 7].T
    # scatGaAs_hole[8, :, 2] = scatGaAs_hole[7, :, 2] + scat_l[:, 8].T
    # scatGaAs_hole[9, :, 2] = scatGaAs_hole[8, :, 2] + scat_l[:, 9].T
    # scatGaAs_hole[10, :, 2] = scatGaAs_hole[9, :, 2] + scat_l[:, 10].T
    # scatGaAs_hole[11, :, 2] = scatGaAs_hole[10, :, 2] + scat_l[:, 11].T
    # scatGaAs_hole[12, :, 2] = scatGaAs_hole[11, :, 2] + scat_l[:, 12].T
    # scatGaAs_hole[:, :, 2] = scatGaAs_hole[:, :, 2] / Gml
