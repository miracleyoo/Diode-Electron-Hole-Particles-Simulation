# coding: utf-8
# Author: Zhongyang Zhang

from smop.libsmop import *
from numpy import dot as ndot
from numpy import pi


def dot(A, B):
    if type(A) == float or type(B) == float:
        return A * B
    else:
        return ndot(A, B)


# make_GaAs_scatTable.m


@function
def make_GaAs_scatTable(T=None, spherical_only=None, de=None, Vmax=None, cimp=None, *args, **kwargs):
    varargin = make_GaAs_scatTable.varargin
    nargin = make_GaAs_scatTable.nargin

    # --------Electron Energy Steps for Formulas/Graphs-------------------------
    delt_Ek = de

    Ek_pts = int(Vmax / de)

    print(Ek_pts)
    eV_axis = dot((arange(1, Ek_pts)), delt_Ek)

    scatG = matlabarray(np.zeros((Ek_pts, 6)))
    # scatG[Ek_pts, 6] = 0
    scatL = matlabarray(np.zeros((Ek_pts, 8)))
    # scatL[Ek_pts, 8] = 0

    # ---------General Constants------------------------------------------------
    q = 1.60219e-19

    bk = 1.38066e-23

    h = 1.05459e-34

    emR = 9.10953e-31

    eps_o = 8.85419e-12

    # ---------GaAs Specific Constants------------------------------------------
    Eg = 1.424

    Egg = 0

    Egl = 0.29

    emG = dot(0.067, emR)

    emL = dot(0.35, emR)

    kconst_G = sqrt(dot(2, emG)) / h

    kconst_L = sqrt(dot(2, emL)) / h

    alpha_G = dot((1 / Eg), (1 - emG / emR) ** 2)

    alpha_L = dot((1 / (Eg + Egl)), (1 - emL / emR) ** 2)

    eps_stat = dot(12.9, eps_o)

    eps_inf = dot(10.92, eps_o)

    eps_p = 1 / ((1 / eps_inf) - (1 / eps_stat))
    # ---------GaAs Constants for Acoustic Phonon Scattering--------------------
    rho = 5360

    sv = 5240

    cl = dot(dot(rho, sv), sv)

    adp = dot(7, q)

    acoustic_const = (dot(dot(dot(dot(dot(2, pi), adp), adp), bk), T)) / (dot(h, cl))

    # --------GaAs Constants for Polar Optical Phonon Scattering----------------
    hwo = 0.03536

    hwoq = dot(hwo, q)
    wo = hwoq / h

    no = 1 / (exp((hwoq) / (dot(bk, T))) - 1)

    pop_constG = dot(((dot(dot(q, q), wo)) / (dot(dot(8, pi), eps_p))), kconst_G)

    pop_constL = dot(((dot(dot(q, q), wo)) / (dot(dot(8, pi), eps_p))), kconst_L)

    # -------GaAs Constants for Intervalley Scattering - Non Polar Optical Phonon
    hwij = 0.03

    hwijq = dot(hwij, q)
    wij = hwijq / h
    dij = dot(1e+11, q)

    Zgl = 4

    Zlg = 1

    nij = 1 / (exp(hwijq / (dot(bk, T))) - 1)

    npop_constGL = (dot(dot(dot(pi, dij), dij), Zgl)) / (dot(rho, wij))
    npop_constLG = (dot(dot(dot(pi, dij), dij), Zlg)) / (dot(rho, wij))
    # ------GaAs Constants for Intravalley Scattering - Non Polar Optical Phonon
    hwll = hwij
    hwllq = dot(hwll, q)
    Zll = Zgl - 1
    dll = dij
    wll = hwllq / h
    nll = 1 / (exp(hwllq / (dot(bk, T))) - 1)
    npop_constLL = (dot(dot(dot(pi, dll), dll), Zll)) / (dot(rho, wll))
    # -------GaAs Constants for Impurity Scattering-----------------------------
    qD = sqrt(dot(dot(q, q), cimp) / (dot(dot(eps_stat, bk), T)))

    NI = cimp
    imp_const = (dot(dot(dot(dot(dot(dot(2, pi), NI), q), q), q), q)) / (dot(dot(h, eps_stat), eps_stat))
    # ------------Acoustic Phonon Scattering------------------------------------
    for i in arange(1, Ek_pts).reshape(-1):
        Ek = dot(delt_Ek, i)
        # For Spherical Parabolic
        if spherical_only == 1:
            N_Ek = dot(dot((((dot(2, emG)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), sqrt(Ek)),
                       sqrt(q))
        else:
            gamma_part = dot(dot(sqrt(dot(Ek, (1 + dot(alpha_G, Ek)))), (1 + dot(dot(2, alpha_G), Ek))), sqrt(q))
            N_Ek = dot((((dot(2, emG)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), gamma_part)
        scatG[i, 1] = dot(acoustic_const, N_Ek)
        if spherical_only == 1:
            N_Ek = dot(dot((((dot(2, emL)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), sqrt(Ek)),
                       sqrt(q))
        else:
            gamma_part = dot(dot(sqrt(dot(Ek, (1 + dot(alpha_L, Ek)))), (1 + dot(dot(2, alpha_L), Ek))), sqrt(q))
            N_Ek = dot((((dot(2, emL)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), gamma_part)
        scatL[i, 1] = dot(acoustic_const, N_Ek)

    # -----------Polar Optical Phonon Scattering--------------------------------
    for i in arange(1, Ek_pts).reshape(-1):
        Ek = dot(delt_Ek, i)
        # CaseI: 1+hwo/Ek (Absorption)
        qmax = sqrt(Ek) + sqrt(Ek + hwo)
        qmin = sqrt(Ek + hwo) - sqrt(Ek)
        scatG[i, 2] = dot(((dot(dot(pop_constG, no), sqrt(dot(Ek, q)))) / (dot(Ek, q))), log(qmax / qmin))
        if Ek - hwo > 0:
            qmax = sqrt(Ek) + sqrt(Ek - hwo)
            qmin = sqrt(Ek) - sqrt(Ek - hwo)
            scatG[i, 3] = dot((dot(dot(pop_constG, (no + 1)), sqrt(dot(Ek, q))) / (dot(Ek, q))), log(qmax / qmin))
        else:
            scatG[i, 3] = 0
            # L Band------------------------------------------
            # CaseI: 1+hwo/Ek (Absorption)
        qmax = sqrt(Ek) + sqrt(Ek + hwo)
        qmin = sqrt(Ek + hwo) - sqrt(Ek)
        scatL[i, 2] = dot(((dot(dot(pop_constL, no), sqrt(dot(Ek, q)))) / (dot(Ek, q))), log(qmax / qmin))
        if Ek - hwo > 0:
            qmax = sqrt(Ek) + sqrt(Ek - hwo)
            qmin = sqrt(Ek) - sqrt(Ek - hwo)
            scatL[i, 3] = dot((dot(dot(pop_constL, (no + 1)), sqrt(dot(Ek, q))) / (dot(Ek, q))), log(qmax / qmin))
        else:
            scatL[i, 3] = 0

    # ---------Intervalley Scattering - Non Polar Optical Phonons---------------
    for i in arange(1, Ek_pts).reshape(-1):
        Ek = dot(delt_Ek, i)
        # CaseI (Absorption):
        Ek2 = Ek + hwij - (Egl - Egg)
        if Ek2 > 0:
            # For Spherical Parabolic
            if spherical_only == 1:
                N_Ek2 = dot(
                    dot((((dot(2, emL)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), sqrt(Ek2)),
                    sqrt(q))
            else:
                gamma_part = dot(dot(sqrt(dot(Ek2, (1 + dot(alpha_L, Ek2)))), (1 + dot(dot(2, alpha_L), Ek2))), sqrt(q))
                N_Ek2 = dot((((dot(2, emL)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), gamma_part)
            scatG[i, 4] = dot(dot(npop_constGL, N_Ek2), nij)
        else:
            scatG[i, 4] = 0
        # CaseII (Emission):
        Ek2 = Ek - hwij - (Egl - Egg)
        if Ek2 > 0:
            # For Spherical Parabolic
            if spherical_only == 1:
                N_Ek2 = dot(
                    dot((((dot(2, emL)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), sqrt(Ek2)),
                    sqrt(q))
            else:
                gamma_part = dot(dot(sqrt(dot(Ek2, (1 + dot(alpha_L, Ek2)))), (1 + dot(dot(2, alpha_L), Ek2))), sqrt(q))
                N_Ek2 = dot((((dot(2, emL)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), gamma_part)
            scatG[i, 5] = dot(dot(npop_constGL, N_Ek2), (nij + 1))
        else:
            scatG[i, 5] = 0
            # L to Gamma------------------------------------------------------------
            # CaseI (Absorption):
        Ek2 = Ek + hwij - (Egg - Egl)
        if Ek2 > 0:
            # For Spherical Parabolic
            if spherical_only == 1:
                N_Ek2 = dot(
                    dot((((dot(2, emG)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), sqrt(Ek2)),
                    sqrt(q))
            else:
                gamma_part = dot(dot(sqrt(dot(Ek2, (1 + dot(alpha_G, Ek2)))), (1 + dot(dot(2, alpha_G), Ek2))), sqrt(q))
                N_Ek2 = dot((((dot(2, emG)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), gamma_part)
            scatL[i, 4] = dot(dot(npop_constLG, N_Ek2), nij)
        else:
            scatL[i, 4] = 0
        # CaseII (Emission):
        Ek2 = Ek - hwij - (Egg - Egl)
        if Ek2 > 0:
            # For Spherical Parabolic
            if spherical_only == 1:
                N_Ek2 = dot(
                    dot((((dot(2, emG)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), sqrt(Ek2)),
                    sqrt(q))
            else:
                gamma_part = dot(dot(sqrt(dot(Ek2, (1 + dot(alpha_G, Ek2)))), (1 + dot(dot(2, alpha_G), Ek2))), sqrt(q))
                N_Ek2 = dot((((dot(2, emG)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), gamma_part)
            scatL[i, 5] = dot(dot(npop_constLG, N_Ek2), (nij + 1))
        else:
            scatL[i, 5] = 0

            # ---------Intravalley Scattering - Non Polar Optical Phonon----------------
            # For L Valley Only
    for i in arange(1, Ek_pts).reshape(-1):
        Ek = dot(delt_Ek, i)
        Ek2 = Ek + hwll
        if spherical_only == 1:
            N_Ek2 = dot(dot((((dot(2, emL)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), sqrt(Ek2)),
                        sqrt(q))
        else:
            gamma_part = dot(dot(sqrt(dot(Ek2, (1 + dot(alpha_L, Ek2)))), (1 + dot(dot(2, alpha_L), Ek2))), sqrt(q))
            N_Ek2 = dot((((dot(2, emL)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), gamma_part)
        scatL[i, 6] = dot(dot(npop_constLL, N_Ek2), nll)
        Ek2 = Ek - hwll
        if Ek2 > 0:
            if spherical_only == 1:
                N_Ek2 = dot(
                    dot((((dot(2, emL)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), sqrt(Ek2)),
                    sqrt(q))
            else:
                gamma_part = dot(dot(sqrt(dot(Ek2, (1 + dot(alpha_L, Ek2)))), (1 + dot(dot(2, alpha_L), Ek2))), sqrt(q))
                N_Ek2 = dot((((dot(2, emL)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), gamma_part)
            scatL[i, 7] = dot(dot(npop_constLL, N_Ek2), (1 + nll))
        else:
            scatL[i, 7] = 0

    # --------Impurity Scattering-----------------------------------------------
    for i in arange(1, Ek_pts).reshape(-1):
        # Gamma Band---------------------------------
        Ek = dot(delt_Ek, i)
        if spherical_only == 1:
            kpart = sqrt(dot(dot(dot(2, emG), Ek), q)) / h
            denom = dot(dot(qD, qD), (dot(dot(4, kpart), kpart) + dot(qD, qD)))
            N_Ek = dot(dot((((dot(2, emG)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), sqrt(Ek)),
                       sqrt(q))
        else:
            kpart = dot(sqrt(dot(Ek, (1 + dot(alpha_G, Ek)))), sqrt(q))
            denom = dot(dot(qD, qD), (dot(dot(dot(dot(4, kconst_G), kconst_G), kpart), kpart) + dot(qD, qD)))
            gamma_part = dot(dot(sqrt(dot(Ek, (1 + dot(alpha_G, Ek)))), (1 + dot(dot(2, alpha_G), Ek))), sqrt(q))
            N_Ek = dot((((dot(2, emG)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), gamma_part)
        scatG[i, 6] = dot(dot(imp_const, (1 / denom)), N_Ek)
        Ek = dot(delt_Ek, i)
        if spherical_only == 1:
            kpart = sqrt(dot(dot(dot(2, emL), Ek), q)) / h
            denom = dot(dot(qD, qD), (dot(dot(4, kpart), kpart) + dot(qD, qD)))
            N_Ek = dot(dot((((dot(2, emL)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), sqrt(Ek)),
                       sqrt(q))
        else:
            kpart = dot(sqrt(dot(Ek, (1 + dot(alpha_L, Ek)))), sqrt(q))
            denom = dot(dot(qD, qD), (dot(dot(dot(dot(4, kconst_L), kconst_L), kpart), kpart) + dot(qD, qD)))
            gamma_part = dot(dot(sqrt(dot(Ek, (1 + dot(alpha_L, Ek)))), (1 + dot(dot(2, alpha_L), Ek))), sqrt(q))
            N_Ek = dot((((dot(2, emL)) ** (3 / 2)) / (dot(dot(dot(dot(dot(4, pi), pi), h), h), h))), gamma_part)
        scatL[i, 8] = dot(dot(imp_const, (1 / denom)), N_Ek)

    # Total Scattering Rate = Sum(All Scatter Mechanisms)
    totScatG = scatG[:, 1] + scatG[:, 2] + scatG[:, 3] + scatG[:, 4] + scatG[:, 5] + scatG[:, 6]
    totScatL = scatL[:, 1] + scatL[:, 2] + scatL[:, 3] + scatL[:, 4] + scatL[:, 5] + scatL[:, 6] + scatL[:, 7] + scatL[:, 8]
    GmG = max(totScatG[:, 1])
    GmL = max(totScatL[:, 1])
    # scatGaAs[8, Ek_pts, 2] = 0
    scatGaAs = matlabarray(np.zeros((8, Ek_pts, 2)))
    scatGaAs[np.ix_([1], np.arange(1, Ek_pts+1), [1])] = scatG[:, 1].T
    for i in range(2, 7):
        scatGaAs[np.ix_([i], np.arange(1, Ek_pts+1), [1])] = scatGaAs[np.ix_([i-1], np.arange(1, Ek_pts+1), [1])] + scatG[:, i].T
    scatGaAs[np.ix_(np.arange(1, 7), np.arange(1, Ek_pts+1), [1])] /= GmG

    # scatGaAs[2, :, 1] = scatGaAs[1, :, 1] + scatG[:, 2].T
    # scatGaAs[3, :, 1] = scatGaAs[2, :, 1] + scatG[:, 3].T
    # scatGaAs[4, :, 1] = scatGaAs[3, :, 1] + scatG[:, 4].T
    # scatGaAs[5, :, 1] = scatGaAs[4, :, 1] + scatG[:, 5].T
    # scatGaAs[6, :, 1] = scatGaAs[5, :, 1] + scatG[:, 6].T
    # scatGaAs[:, :, 1] = scatGaAs[:, :, 1] / GmG

    mtemp3 = {2: 4, 3: 5, 4: 2, 5: 3, 6: 6, 7: 7, 8: 8}
    scatGaAs[np.ix_([1], np.arange(1, Ek_pts+1), [2])] = scatL[:, 1].T
    for i in range(2, 7):
        scatGaAs[np.ix_([i], np.arange(1, Ek_pts+1), [1])] = scatGaAs[np.ix_([i - 1], np.arange(1, Ek_pts+1), [2])] + scatL[:, mtemp3[i]].T
    scatGaAs[np.ix_(np.arange(1, 7), np.arange(1, Ek_pts+1), [2])] /= GmG
    # scatGaAs[1, :, 2] = scatL[:, 1].T
    # scatGaAs[2, :, 2] = scatGaAs[1, :, 2] + scatL[:, 4].T
    # scatGaAs[3, :, 2] = scatGaAs[2, :, 2] + scatL[:, 5].T
    # scatGaAs[4, :, 2] = scatGaAs[3, :, 2] + scatL[:, 2].T
    # scatGaAs[5, :, 2] = scatGaAs[4, :, 2] + scatL[:, 3].T
    # scatGaAs[6, :, 2] = scatGaAs[5, :, 2] + scatL[:, 6].T
    # scatGaAs[7, :, 2] = scatGaAs[6, :, 2] + scatL[:, 7].T
    # scatGaAs[8, :, 2] = scatGaAs[7, :, 2] + scatL[:, 8].T
    # scatGaAs[:, :, 2] = scatGaAs[:, :, 2] / GmL
    return scatGaAs, GmG, GmL
