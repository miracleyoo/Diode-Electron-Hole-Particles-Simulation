# coding: utf-8
# Author: Zhongyang Zhang

# Generated with SMOP  0.41-beta
from const_val import *


# pn_scat_v2.m


@function
def pn_scat_v2(particle=None, valley=None, scatGaAs=None, scatGaAs_hole=None, de=None, q=None, h=None, eM=None,
               alpha=None, qD=None, hw0=None, A=None, B=None, C=None, m0=None, n=None, hwij=None, Egl=None, Egg=None,
               hwe=None, g100=None, g111=None, *args, **kwargs):
    varargin = pn_scat_v2.varargin
    nargin = pn_scat_v2.nargin

    # function [particle,valleyf]=bulk_scat(particle,valley,scatTable,eM,alpha,xs,xd,qD,hwo,hwij,hwe,Egl,Egg)
    # particle: vector for particle properties (one row of particles matrix)
    # valley: valley index for particle (1:Gamma, 2:L, 9:Delete)
    # scatTable: table with precalculated scattering rates
    # eM: effective electron mass vector
    # alpha: non-parabolicity factor (vector for Gamma and L)
    # xs: x-coord of end of cathode region
    # xd: x-coord of start of anode region
    # qD: inverse debye length

    kx = particle[1, 1]
    ky = particle[1, 2]
    kz = particle[1, 3]
    a = dot(dot(h, h), abs(A)) / (dot(2, m0))
    b = qD
    skx = dot(kx, kx)
    sky = dot(ky, ky)
    skz = dot(kz, kz)
    sk = abs(skx + sky + skz)
    ki = sqrt(sk)
    iemax = length(scatGaAs[1, :, 1])

    kxf = kx
    kyf = ky
    kzf = kz
    valleyf = valley
    g1 = 0
    ei = 0
    if sk != 0:
        gk = dot(dot((dot(h, h) / (dot(2, eM(valley)))), sk), (1 / q))
        if valley == 1 or valley == 2:
            ei = (sqrt(1 + dot(dot(4, alpha(valley)), gk)) - 1) / (dot(2, alpha(valley)))
            particle[1, 7] = ei
        else:
            if valley == 3:
                cos_theta = kz / ki
                sin_theta = sqrt(1 - cos_theta ** 2)
                cos_phi = kx / (dot(ki, sin_theta))
                sin_phi = ky / (dot(ki, sin_theta))
                g1 = ((B / A) ** 2 + dot((C / A) ** 2, (
                    dot(sin_theta ** 2, cos_theta ** 2) + dot(dot(sin_theta ** 4, cos_phi ** 2), sin_phi ** 2)))) ** 0.5
                ei = dot(dot(dot((dot(dot(h, h), abs(A)) / (dot(2, m0))), sk), (1 / q)), (1 - g1))
                ef = ei
                particle[1, 7] = ei
            else:
                if valley == 4:
                    cos_theta = kz / ki
                    sin_theta = sqrt(1 - cos_theta ** 2)
                    cos_phi = kx / (dot(ki, sin_theta))
                    sin_phi = ky / (dot(ki, sin_theta))
                    g1 = ((B / A) ** 2 + dot((C / A) ** 2, (
                        dot(sin_theta ** 2, cos_theta ** 2) + dot(dot(sin_theta ** 4, cos_phi ** 2),
                                                                  sin_phi ** 2)))) ** 0.5
                    ei = dot(dot(dot((dot(dot(h, h), abs(A)) / (dot(2, m0))), sk), (1 / q)), (1 + g1))
                    ef = ei
                    particle[1, 7] = ei
        ie = abs(floor(ei / de) + 1)
        if ie > iemax:
            ie = iemax
        r1 = rand()
        # Gamma Valley
        if valley == 1:
            # Acoustic Scattering
            if r1 < scatGaAs[1, ie, valley]:
                kf = ki
                phi = dot(dot(2, pi), rand())
                cos_t = 1 - dot(2, rand())
                sin_t = sqrt(1 - dot(cos_t, cos_t))
                kxf = dot(dot(kf, sin_t), cos(phi))
                kyf = dot(dot(kf, sin_t), sin(phi))
                kzf = dot(kf, cos_t)
                #             valleyf=1;
                # POP Absorption
            else:
                if r1 < scatGaAs[2, ie, valley]:
                    ef = ei + hw0
                    ff = dot(2, sqrt(dot(ei, ef))) / ((sqrt(ei) - sqrt(ef)) ** 2)
                    if ff > 0:
                        kf = dot(dot((sqrt(dot(2, eM(valley))) / h), sqrt(dot(ef, (1 + dot(alpha(valley), ef))))),
                                 sqrt(q))
                        phi = dot(dot(2, pi), rand())
                        cos_t = (1 + ff - (1 + dot(2, ff)) ** rand()) / ff
                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                        sin_a = sqrt(skx + sky) / ki
                        cos_a = kz / ki
                        sin_b = kx / sqrt(skx + sky)
                        cos_b = ky / sqrt(skx + sky)
                        kxr = dot(dot(kf, sin_t), cos(phi))
                        kyr = dot(dot(kf, sin_t), sin(phi))
                        kzr = dot(kf, cos_t)
                        kxf = dot(cos_b, kxr) + dot(dot(cos_a, sin_b), kyr) + dot(dot(sin_a, sin_b), kzr)
                        kyf = dot(- sin_b, kxr) + dot(dot(cos_a, cos_b), kyr) + dot(dot(sin_a, cos_b), kzr)
                        kzf = dot(- sin_a, kyr) + dot(cos_a, kzr)
                        #             valleyf=1;
                        # POP Emission
                else:
                    if r1 < scatGaAs[3, ie, valley]:
                        ef = ei - hw0
                        if ef > 0:
                            ff = dot(2, sqrt(dot(ei, ef))) / ((sqrt(ei) - sqrt(ef)) ** 2)
                            if ff > 0:
                                kf = dot(
                                    dot((sqrt(dot(2, eM(valley))) / h), sqrt(dot(ef, (1 + dot(alpha(valley), ef))))),
                                    sqrt(q))
                                phi = dot(dot(2, pi), rand())
                                cos_t = (1 + ff - (1 + dot(2, ff)) ** rand()) / ff
                                sin_t = sqrt(1 - dot(cos_t, cos_t))
                                sin_a = sqrt(skx + sky) / ki
                                cos_a = kz / ki
                                sin_b = kx / sqrt(skx + sky)
                                cos_b = ky / sqrt(skx + sky)
                                kxr = dot(dot(kf, sin_t), cos(phi))
                                kyr = dot(dot(kf, sin_t), sin(phi))
                                kzr = dot(kf, cos_t)
                                kxf = dot(cos_b, kxr) + dot(dot(cos_a, sin_b), kyr) + dot(dot(sin_a, sin_b), kzr)
                                kyf = dot(- sin_b, kxr) + dot(dot(cos_a, cos_b), kyr) + dot(dot(sin_a, cos_b), kzr)
                                kzf = dot(- sin_a, kyr) + dot(cos_a, kzr)
                                #             valleyf=1;
                                # NPOP Absorption
                    else:
                        if r1 < scatGaAs[4, ie, valley]:
                            ef = ei + hwij - (Egl - Egg)
                            valleyf = 2
                            kf = dot(dot((sqrt(dot(2, eM(valleyf))) / h), sqrt(dot(ef, (1 + dot(alpha(valleyf), ef))))),
                                     sqrt(q))
                            phi = dot(dot(2, pi), rand())
                            cos_t = 1 - dot(2, rand())
                            sin_t = sqrt(1 - dot(cos_t, cos_t))
                            kxf = dot(dot(kf, sin_t), cos(phi))
                            kyf = dot(dot(kf, sin_t), sin(phi))
                            kzf = dot(kf, cos_t)
                        else:
                            if r1 < scatGaAs[5, ie, valley]:
                                ef = ei - hwij - (Egl - Egg)
                                valleyf = valley
                                if ef > 0:
                                    valleyf = 2
                                    kf = dot(dot((sqrt(dot(2, eM(valleyf))) / h),
                                                 sqrt(dot(ef, (1 + dot(alpha(valleyf), ef))))), sqrt(q))
                                    phi = dot(dot(2, pi), rand())
                                    cos_t = 1 - dot(2, rand())
                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                    kxf = dot(dot(kf, sin_t), cos(phi))
                                    kyf = dot(dot(kf, sin_t), sin(phi))
                                    kzf = dot(kf, cos_t)
                                    # Impurity
                            else:
                                if r1 < scatGaAs[6, ie, valley]:
                                    kf = ki
                                    r = rand()
                                    phi = dot(dot(2, pi), rand())
                                    cos_t = 1 - (dot(2, r)) / (1 + dot((1 - r), (dot(2, ki) / qD) ** 2))
                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                    sin_a = sqrt(skx + sky) / ki
                                    cos_a = kz / ki
                                    sin_b = kx / sqrt(skx + sky)
                                    cos_b = ky / sqrt(skx + sky)
                                    kxr = dot(dot(kf, sin_t), cos(phi))
                                    kyr = dot(dot(kf, sin_t), sin(phi))
                                    kzr = dot(kf, cos_t)
                                    kxf = dot(cos_b, kxr) + dot(dot(cos_a, sin_b), kyr) + dot(dot(sin_a, sin_b), kzr)
                                    kyf = dot(- sin_b, kxr) + dot(dot(cos_a, cos_b), kyr) + dot(dot(sin_a, cos_b), kzr)
                                    kzf = dot(- sin_a, kyr) + dot(cos_a, kzr)
                                    valleyf = 1
                                    # L Valley--------------------------------------------------------------
        else:
            if valley == 2:
                # Acoustic
                if r1 < scatGaAs[1, ie, valley]:
                    kf = ki
                    phi = dot(dot(2, pi), rand())
                    cos_t = 1 - dot(2, rand())
                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                    kxf = dot(dot(kf, sin_t), cos(phi))
                    kyf = dot(dot(kf, sin_t), sin(phi))
                    kzf = dot(kf, cos_t)
                    valleyf = 2
                else:
                    if r1 < scatGaAs[2, ie, valley]:
                        ef = ei + hwij + (Egl - Egg)
                        valleyf = 1
                        kf = dot(dot((sqrt(dot(2, eM(valleyf))) / h), sqrt(dot(ef, (1 + dot(alpha(valleyf), ef))))),
                                 sqrt(q))
                        phi = dot(dot(2, pi), rand())
                        cos_t = 1 - dot(2, rand())
                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                        kxf = dot(dot(kf, sin_t), cos(phi))
                        kyf = dot(dot(kf, sin_t), sin(phi))
                        kzf = dot(kf, cos_t)
                    else:
                        if r1 < scatGaAs[3, ie, valley]:
                            ef = ei - hwij + (Egl - Egg)
                            valleyf = valley
                            if ef > 0:
                                valleyf = 1
                                kf = dot(
                                    dot((sqrt(dot(2, eM(valleyf))) / h), sqrt(dot(ef, (1 + dot(alpha(valleyf), ef))))),
                                    sqrt(q))
                                phi = dot(dot(2, pi), rand())
                                cos_t = 1 - dot(2, rand())
                                sin_t = sqrt(1 - dot(cos_t, cos_t))
                                kxf = dot(dot(kf, sin_t), cos(phi))
                                kyf = dot(dot(kf, sin_t), sin(phi))
                                kzf = dot(kf, cos_t)
                                # POP Absorption
                        else:
                            if r1 < scatGaAs[4, ie, valley]:
                                ef = ei + hw0
                                ff = dot(2, sqrt(dot(ei, ef))) / ((sqrt(ei) - sqrt(ef)) ** 2)
                                if ff > 0:
                                    kf = dot(dot((sqrt(dot(2, eM(valley))) / h),
                                                 sqrt(dot(ef, (1 + dot(alpha(valley), ef))))), sqrt(q))
                                    phi = dot(dot(2, pi), rand())
                                    cos_t = (1 + ff - (1 + dot(2, ff)) ** rand()) / ff
                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                    sin_a = sqrt(skx + sky) / ki
                                    cos_a = kz / ki
                                    sin_b = kx / sqrt(skx + sky)
                                    cos_b = ky / sqrt(skx + sky)
                                    kxr = dot(dot(kf, sin_t), cos(phi))
                                    kyr = dot(dot(kf, sin_t), sin(phi))
                                    kzr = dot(kf, cos_t)
                                    kxf = dot(cos_b, kxr) + dot(dot(cos_a, sin_b), kyr) + dot(dot(sin_a, sin_b), kzr)
                                    kyf = dot(- sin_b, kxr) + dot(dot(cos_a, cos_b), kyr) + dot(dot(sin_a, cos_b), kzr)
                                    kzf = dot(- sin_a, kyr) + dot(cos_a, kzr)
                                valleyf = 2
                            else:
                                if r1 < scatGaAs[5, ie, valley]:
                                    ef = ei - hw0
                                    if ef > 0:
                                        ff = dot(2, sqrt(dot(ei, ef))) / ((sqrt(ei) - sqrt(ef)) ** 2)
                                        if ff > 0:
                                            kf = dot(dot((sqrt(dot(2, eM(valley))) / h),
                                                         sqrt(dot(ef, (1 + dot(alpha(valley), ef))))), sqrt(q))
                                            phi = dot(dot(2, pi), rand())
                                            cos_t = (1 + ff - (1 + dot(2, ff)) ** rand()) / ff
                                            sin_t = sqrt(1 - dot(cos_t, cos_t))
                                            sin_a = sqrt(skx + sky) / ki
                                            cos_a = kz / ki
                                            sin_b = kx / sqrt(skx + sky)
                                            cos_b = ky / sqrt(skx + sky)
                                            kxr = dot(dot(kf, sin_t), cos(phi))
                                            kyr = dot(dot(kf, sin_t), sin(phi))
                                            kzr = dot(kf, cos_t)
                                            kxf = dot(cos_b, kxr) + dot(dot(cos_a, sin_b), kyr) + dot(dot(sin_a, sin_b),
                                                                                                      kzr)
                                            kyf = dot(- sin_b, kxr) + dot(dot(cos_a, cos_b), kyr) + dot(
                                                dot(sin_a, cos_b), kzr)
                                            kzf = dot(- sin_a, kyr) + dot(cos_a, kzr)
                                    valleyf = 2
                                else:
                                    if r1 < scatGaAs[6, ie, valley]:
                                        ef = ei + hwe
                                        valleyf = 2
                                        kf = dot(dot((sqrt(dot(2, eM(valleyf))) / h),
                                                     sqrt(dot(ef, (1 + dot(alpha(valleyf), ef))))), sqrt(q))
                                        phi = dot(dot(2, pi), rand())
                                        cos_t = 1 - dot(2, rand())
                                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                                        kxf = dot(dot(kf, sin_t), cos(phi))
                                        kyf = dot(dot(kf, sin_t), sin(phi))
                                        kzf = dot(kf, cos_t)
                                    else:
                                        if r1 < scatGaAs[7, ie, valley]:
                                            ef = ei - hwe
                                            valleyf = 2
                                            kf = dot(dot((sqrt(dot(2, eM(valleyf))) / h),
                                                         sqrt(dot(ef, (1 + dot(alpha(valleyf), ef))))), sqrt(q))
                                            phi = dot(dot(2, pi), rand())
                                            cos_t = 1 - dot(2, rand())
                                            sin_t = sqrt(1 - dot(cos_t, cos_t))
                                            kxf = dot(dot(kf, sin_t), cos(phi))
                                            kyf = dot(dot(kf, sin_t), sin(phi))
                                            kzf = dot(kf, cos_t)
                                        else:
                                            if r1 < scatGaAs[8, ie, valley]:
                                                #             if x < Ll || x > (Ll+Lg)
                                                kf = ki
                                                r = rand()
                                                phi = dot(dot(2, pi), rand())
                                                cos_t = 1 - (dot(2, r)) / (1 + dot((1 - r), (dot(2, ki) / qD) ** 2))
                                                sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                sin_a = sqrt(skx + sky) / ki
                                                cos_a = kz / ki
                                                sin_b = kx / sqrt(skx + sky)
                                                cos_b = ky / sqrt(skx + sky)
                                                kxr = dot(dot(kf, sin_t), cos(phi))
                                                kyr = dot(dot(kf, sin_t), sin(phi))
                                                kzr = dot(kf, cos_t)
                                                kxf = dot(cos_b, kxr) + dot(dot(cos_a, sin_b), kyr) + dot(
                                                    dot(sin_a, sin_b), kzr)
                                                kyf = dot(- sin_b, kxr) + dot(dot(cos_a, cos_b), kyr) + dot(
                                                    dot(sin_a, cos_b), kzr)
                                                kzf = dot(- sin_a, kyr) + dot(cos_a, kzr)
                                                valleyf = 2
            else:
                if valley == 3:
                    valley = valley - 2
                    if r1 < scatGaAs_hole(1, ie, valley):
                        ef = ei
                        phi = dot(dot(2, pi), rand())
                        cos_t = 1 - dot(2, rand())
                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                        g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                            dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2), sin(phi) ** 2)))) ** 0.5
                        kf = sqrt(dot(dot(dot(2, m0), q), ef) / h / h / abs(A) / (1 + g2))
                        kxf = dot(dot(kf, sin_t), cos(phi))
                        kyf = dot(dot(kf, sin_t), sin(phi))
                        kzf = dot(kf, cos_t)
                        valleyf = 4
                    else:
                        if r1 < scatGaAs_hole(2, ie, valley):
                            kf = ki
                            phi = dot(dot(2, pi), rand())
                            cos_t = 1 - dot(2, rand())
                            sin_t = sqrt(1 - dot(cos_t, cos_t))
                            kxf = dot(dot(kf, sin_t), cos(phi))
                            kyf = dot(dot(kf, sin_t), sin(phi))
                            kzf = dot(kf, cos_t)
                            valleyf = 3
                        else:
                            if r1 < scatGaAs_hole(3, ie, valley):
                                r = rand()
                                ef = ei + hw0
                                phi = dot(dot(2, pi), rand())
                                cos_t = 1 - dot(2, rand())
                                sin_t = sqrt(1 - dot(cos_t, cos_t))
                                g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                    dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                      sin(phi) ** 2)))) ** 0.5
                                kf = sqrt(dot(dot(dot(2, m0), q), ef) / h / h / abs(A) / (1 - g2))
                                kxf = dot(dot(kf, sin_t), cos(phi))
                                kyf = dot(dot(kf, sin_t), sin(phi))
                                kzf = dot(kf, cos_t)
                                cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                G = dot(0.25, (1 + dot(3, cos_alpha ** 2)))
                                cos_theta = kz / ki
                                fa = dot(((1 - g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta),
                                                                              sqrt(dot((1 - g2) / (1 - g1), ef) / ei))),
                                         G)
                                fb = dot(sqrt(1 - g2), ((1 - g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                    dot((1 - g2) / (1 - g1), ef) / ei)) + dot(dot(dot(a, b), b), (1 - g2))) ** 2)
                                f = fa / fb
                                fmax = (1 + ef / ei - dot(dot(2, sqrt(ef / ei)), cos_theta)) / \
                                       (1 - g111) ** 0.5 / (1 + ef / ei - dot(dot(2, sqrt(ef / ei)), cos_theta) +
                                                            dot(dot(dot(a, b), b), (1 - g111))) ** 2
                                valleyf = 3
                                i = 1
                                while r > f / fmax:

                                    phi = dot(dot(2, pi), rand())
                                    cos_t = 1 - dot(2, rand())
                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                    g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                        dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                          sin(phi) ** 2)))) ** 0.5
                                    kf = sqrt(dot(dot(dot(2, m0), q), ef) / h / h / abs(A) / (1 - g2))
                                    kxf = dot(dot(kf, sin_t), cos(phi))
                                    kyf = dot(dot(kf, sin_t), sin(phi))
                                    kzf = dot(kf, cos_t)
                                    cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                    G = dot(0.25, (1 + dot(3, cos_alpha ** 2)))
                                    fa = dot(((1 - g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                        dot((1 - g2) / (1 - g1), ef) / ei))), G)
                                    fb = dot(sqrt(1 - g2), ((1 - g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                        dot((1 - g2) / (1 - g1), ef) / ei)) + dot(dot(dot(a, b), b), (1 - g2))) ** 2)
                                    f = fa / fb
                                    i = i + 1
                                    if i > 2:
                                        kxf = kx
                                        kyf = ky
                                        kzf = kz
                                        ef = ei
                                        break

                            else:
                                if r1 < scatGaAs_hole(4, ie, valley):
                                    r = rand()
                                    ef = ei + hw0
                                    phi = dot(dot(2, pi), rand())
                                    cos_t = 1 - dot(2, rand())
                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                    g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                        dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                          sin(phi) ** 2)))) ** 0.5
                                    kf = sqrt(dot(dot(dot(2, m0), q), ef) / h / h / abs(A) / (1 + g2))
                                    kxf = dot(dot(kf, sin_t), cos(phi))
                                    kyf = dot(dot(kf, sin_t), sin(phi))
                                    kzf = dot(kf, cos_t)
                                    cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                    G = dot(0.75, (1 - cos_alpha ** 2))
                                    cos_theta = kz / ki
                                    fa = dot(((1 + g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                        dot((1 + g2) / (1 - g1), ef) / ei))), G)
                                    fb = dot(sqrt(1 + g2), ((1 + g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                        dot((1 + g2) / (1 - g1), ef) / ei)) + dot(dot(dot(a, b), b), (1 + g2))) ** 2)
                                    f = fa / fb
                                    fmaxa = dot(0.75, ((1 + g111) / (1 - g111) + (ef / ei) - dot(dot(2, cos_theta),
                                                                                                 sqrt(dot((1 + g111) / (
                                                                                                     1 - g111),
                                                                                                          ef) / ei))))
                                    fmaxb = dot(sqrt(1 + g100), (
                                        (1 + g100) / (1 - g100) + (ef / ei) - dot(dot(2, cos_theta), sqrt(
                                            dot((1 + g100) / (1 - g100), ef) / ei)) + dot(dot(dot(a, b), b),
                                                                                          (1 + g100))) ** 2)
                                    fmax = fmaxa / fmaxb
                                    valleyf = 4
                                    i = 1
                                    while r > f / fmax:

                                        phi = dot(dot(2, pi), rand())
                                        cos_t = 1 - dot(2, rand())
                                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                                        g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                            dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                              sin(phi) ** 2)))) ** 0.5
                                        kf = sqrt(dot(dot(dot(2, m0), q), ef) / h / h / abs(A) / (1 + g2))
                                        kxf = dot(dot(kf, sin_t), cos(phi))
                                        kyf = dot(dot(kf, sin_t), sin(phi))
                                        kzf = dot(kf, cos_t)
                                        cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                        G = dot(0.75, (1 - cos_alpha ** 2))
                                        fa = dot(((1 + g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                            dot((1 + g2) / (1 - g1), ef) / ei))), G)
                                        fb = dot(sqrt(1 + g2), ((1 + g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta),
                                                                                                    sqrt(dot(
                                                                                                        (1 + g2) / (
                                                                                                            1 - g1),
                                                                                                        ef) / ei)) + dot(
                                            dot(dot(a, b), b), (1 + g2))) ** 2)
                                        f = fa / fb
                                        i = i + 1
                                        if i > 2:
                                            kxf = kx
                                            kyf = ky
                                            kzf = kz
                                            ef = ei
                                            valleyf = 3
                                            break

                                else:
                                    if r1 < scatGaAs_hole(5, ie, valley):
                                        r = rand()
                                        ef = ei - hw0
                                        if ef < 0:
                                            valleyf = 3
                                            fprintf('5h is wrong,ef=%f\\n', ef)
                                        else:
                                            phi = dot(dot(2, pi), rand())
                                            cos_t = 1 - dot(2, rand())
                                            sin_t = sqrt(1 - dot(cos_t, cos_t))
                                            g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                                  sin(phi) ** 2)))) ** 0.5
                                            kf = sqrt(dot(dot(dot(2, m0), q), ef) / h / h / abs(A) / (1 - g2))
                                            kxf = dot(dot(kf, sin_t), cos(phi))
                                            kyf = dot(dot(kf, sin_t), sin(phi))
                                            kzf = dot(kf, cos_t)
                                            g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                                  sin(phi) ** 2)))) ** 0.5
                                            cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                            G = dot(0.25, (1 + dot(3, cos_alpha ** 2)))
                                            cos_theta = kz / ki
                                            fa = dot(((1 - g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                                dot((1 - g2) / (1 - g1), ef) / ei))), G)
                                            fb = dot(sqrt(1 - g2), (
                                                (1 - g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                                    dot((1 - g2) / (1 - g1), ef) / ei)) + dot(dot(dot(a, b), b),
                                                                                              (1 - g2))) ** 2)
                                            f = fa / fb
                                            valleyf = 3
                                            fmax = (1 + ef / ei - dot(dot(2, sqrt(ef / ei)), cos_theta)) / \
                                                   (1 - g111) ** 0.5 / (1 + ef / ei - dot(dot(2, sqrt(ef / ei)),
                                                                                          cos_theta) +
                                                                        dot(dot(dot(a, b), b), (1 - g111))) ** 2
                                            i = 1
                                            while r > f / fmax:

                                                phi = dot(dot(2, pi), rand())
                                                cos_t = 1 - dot(2, rand())
                                                sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                    dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                                      sin(phi) ** 2)))) ** 0.5
                                                kf = sqrt(dot(dot(dot(2, m0), q), ef) / h / h / abs(A) / (1 - g2))
                                                kxf = dot(dot(kf, sin_t), cos(phi))
                                                kyf = dot(dot(kf, sin_t), sin(phi))
                                                kzf = dot(kf, cos_t)
                                                g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                    dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                                      sin(phi) ** 2)))) ** 0.5
                                                cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                                G = dot(0.25, (1 + dot(3, cos_alpha ** 2)))
                                                fa = dot(((1 - g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                                    dot((1 - g2) / (1 - g1), ef) / ei))), G)
                                                fb = dot(sqrt(1 - g2), (
                                                    (1 - g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                                        dot((1 - g2) / (1 - g1), ef) / ei)) + dot(dot(dot(a, b), b),
                                                                                                  (1 - g2))) ** 2)
                                                f = fa / fb
                                                i = i + 1
                                                if i > 2:
                                                    kxf = kx
                                                    kyf = ky
                                                    kzf = kz
                                                    ef = ei
                                                    break

                                    else:
                                        if r1 < scatGaAs_hole(6, ie, valley):
                                            r = rand()
                                            ef = ei - hw0
                                            if ef < 0:
                                                fprintf('6h is wrong,ef=%f\\n', ef)
                                                ef = ef + de
                                                valleyf = 4
                                            else:
                                                phi = dot(dot(2, pi), rand())
                                                cos_t = 1 - dot(2, rand())
                                                sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                    dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                                      sin(phi) ** 2)))) ** 0.5
                                                kf = sqrt(dot(dot(dot(2, m0), q), ef) / h / h / abs(A) / (1 + g2))
                                                kxf = dot(dot(kf, sin_t), cos(phi))
                                                kyf = dot(dot(kf, sin_t), sin(phi))
                                                kzf = dot(kf, cos_t)
                                                g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                    dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                                      sin(phi) ** 2)))) ** 0.5
                                                cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                                G = dot(0.75, (1 - cos_alpha ** 2))
                                                cos_theta = kz / ki
                                                fa = dot(((1 + g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                                    dot((1 + g2) / (1 - g1), ef) / ei))), G)
                                                fb = dot(sqrt(1 + g2), (
                                                    (1 + g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                                        dot((1 + g2) / (1 - g1), ef) / ei)) + dot(dot(dot(a, b), b),
                                                                                                  (1 + g2))) ** 2)
                                                f = fa / fb
                                                g111 = sqrt((B / A) ** 2 + (C / A) ** 2 / 3)
                                                g100 = B / A
                                                fmaxa = dot(0.75, (
                                                    (1 + g111) / (1 - g111) + (ef / ei) - dot(dot(2, cos_theta), sqrt(
                                                        dot((1 + g111) / (1 - g111), ef) / ei))))
                                                fmaxb = dot(sqrt(1 + g100), (
                                                    (1 + g100) / (1 - g100) + (ef / ei) - dot(dot(2, cos_theta), sqrt(
                                                        dot((1 + g100) / (1 - g100), ef) / ei)) + dot(dot(dot(a, b), b),
                                                                                                      (1 + g100))) ** 2)
                                                fmax = fmaxa / fmaxb
                                                valleyf = 4
                                                i = 1
                                                while r > f / fmax:

                                                    phi = dot(dot(2, pi), rand())
                                                    cos_t = 1 - dot(2, rand())
                                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                    g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                        dot(sin_t ** 2, cos_t ** 2) + dot(
                                                            dot(sin_t ** 4, cos(phi) ** 2),
                                                            sin(phi) ** 2)))) ** 0.5
                                                    kf = sqrt(dot(dot(dot(2, m0), q), ef) / h / h / abs(A) / (1 + g2))
                                                    kxf = dot(dot(kf, sin_t), cos(phi))
                                                    kyf = dot(dot(kf, sin_t), sin(phi))
                                                    kzf = dot(kf, cos_t)
                                                    g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                        dot(sin_t ** 2, cos_t ** 2) + dot(
                                                            dot(sin_t ** 4, cos(phi) ** 2),
                                                            sin(phi) ** 2)))) ** 0.5
                                                    cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                                    G = dot(0.75, (1 - cos_alpha ** 2))
                                                    fa = dot(((1 + g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta),
                                                                                                  sqrt(dot((1 + g2) / (
                                                                                                      1 - g1),
                                                                                                           ef) / ei))),
                                                             G)
                                                    fb = dot(sqrt(1 + g2), (
                                                        (1 + g2) / (1 - g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                                            dot((1 + g2) / (1 - g1), ef) / ei)) + dot(dot(dot(a, b), b),
                                                                                                      (1 + g2))) ** 2)
                                                    f = fa / fb
                                                    i = i + 1
                                                    if i > 2:
                                                        kxf = kx
                                                        kyf = ky
                                                        kzf = kz
                                                        valleyf = 3
                                                        ef = ei
                                                        break

                                        else:
                                            if r1 < scatGaAs_hole(7, ie, valley):
                                                r = rand()
                                                flag = 0
                                                phi = dot(dot(2, pi), rand())
                                                cos_t = 1 - dot(2, rand())
                                                sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                    dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                        dot(sin_t ** 4, cos(phi) ** 2), sin(phi) ** 2)))) ** 0.5
                                                fmax = 1 / (dot(a, (1 - sqrt((B / A) ** 2 + (C / A) ** 2 / 3)))) ** 1.5
                                                f = 1 / (dot(a, (1 - g))) ** 1.5
                                                i = 1
                                                while r > f / fmax:

                                                    phi = dot(dot(2, pi), rand())
                                                    cos_t = 1 - dot(2, rand())
                                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                    g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                        dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                            dot(sin_t ** 4, cos(phi) ** 2), sin(phi) ** 2)))) ** 0.5
                                                    f = 1 / (dot(a, (1 - g))) ** 1.5
                                                    i = i + 1
                                                    if i > 2:
                                                        flag = 1
                                                        break

                                                if flag == 0:
                                                    ef = ei + hw0
                                                    kf = sqrt(dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 - g)) / h
                                                    kxf = dot(dot(kf, sin_t), cos(phi))
                                                    kyf = dot(dot(kf, sin_t), sin(phi))
                                                    kzf = dot(kf, cos_t)
                                                    valleyf = 3
                                                else:
                                                    valleyf = 3
                                            else:
                                                if r1 < scatGaAs_hole(8, ie, valley):
                                                    ef = ei - hw0
                                                    if ef < 0:
                                                        valleyf = 3
                                                    else:
                                                        r = rand()
                                                        flag = 0
                                                        phi = dot(dot(2, pi), rand())
                                                        cos_t = 1 - dot(2, rand())
                                                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                        g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                            dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                dot(sin_t ** 4, cos(phi) ** 2), sin(phi) ** 2)))) ** 0.5
                                                        fmax = 1 / (dot(a, (
                                                            1 - sqrt((B / A) ** 2 + (C / A) ** 2 / 3)))) ** 1.5
                                                        f = 1 / (dot(a, (1 - g))) ** 1.5
                                                        i = 1
                                                        while r > f / fmax:

                                                            phi = dot(dot(2, pi), rand())
                                                            cos_t = 1 - dot(2, rand())
                                                            sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                            g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                    dot(sin_t ** 4, cos(phi) ** 2),
                                                                    sin(phi) ** 2)))) ** 0.5
                                                            f = 1 / (dot(a, (1 - g))) ** 1.5
                                                            i = i + 1
                                                            if i > 2:
                                                                flag = 1
                                                                break

                                                        if flag == 0:
                                                            kf = sqrt(
                                                                dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 - g)) / h
                                                            kxf = dot(dot(kf, sin_t), cos(phi))
                                                            kyf = dot(dot(kf, sin_t), sin(phi))
                                                            kzf = dot(kf, cos_t)
                                                            valleyf = 3
                                                        else:
                                                            valleyf = 3
                                                else:
                                                    if r1 < scatGaAs_hole(9, ie, valley):
                                                        r = rand()
                                                        flag = 0
                                                        phi = dot(dot(2, pi), rand())
                                                        cos_t = 1 - dot(2, rand())
                                                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                        g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                            dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                dot(sin_t ** 4, cos(phi) ** 2), sin(phi) ** 2)))) ** 0.5
                                                        fmax = 1 / (dot(a, (1 + B / A))) ** 1.5
                                                        f = 1 / (dot(a, (1 + g))) ** 1.5
                                                        i = 1
                                                        while r > f / fmax:

                                                            phi = dot(dot(2, pi), rand())
                                                            cos_t = 1 - dot(2, rand())
                                                            sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                            g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                    dot(sin_t ** 4, cos(phi) ** 2),
                                                                    sin(phi) ** 2)))) ** 0.5
                                                            f = 1 / (dot(a, (1 - g))) ** 1.5
                                                            i = i + 1
                                                            if i > 2:
                                                                flag = 1
                                                                break

                                                        if flag == 0:
                                                            ef = ei + hw0
                                                            #                 g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                                                            kf = sqrt(
                                                                dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 + g)) / h
                                                            kxf = dot(dot(kf, sin_t), cos(phi))
                                                            kyf = dot(dot(kf, sin_t), sin(phi))
                                                            kzf = dot(kf, cos_t)
                                                            valleyf = 4
                                                        else:
                                                            valleyf = 3
                                                    else:
                                                        if r1 < scatGaAs_hole(10, ie, valley):
                                                            ef = ei - hw0
                                                            if ef < 0:
                                                                valleyf = 3
                                                            else:
                                                                r = rand()
                                                                flag = 0
                                                                phi = dot(dot(2, pi), rand())
                                                                cos_t = 1 - dot(2, rand())
                                                                sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                                g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                    dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                        dot(sin_t ** 4, cos(phi) ** 2),
                                                                        sin(phi) ** 2)))) ** 0.5
                                                                fmax = 1 / (dot(a, (1 + B / A))) ** 1.5
                                                                f = 1 / (dot(a, (1 + g))) ** 1.5
                                                                i = 1
                                                                while r > f / fmax:

                                                                    phi = dot(dot(2, pi), rand())
                                                                    cos_t = 1 - dot(2, rand())
                                                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                                    g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                        dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                            dot(sin_t ** 4, cos(phi) ** 2),
                                                                            sin(phi) ** 2)))) ** 0.5
                                                                    f = 1 / (dot(a, (1 - g))) ** 1.5
                                                                    i = i + 1
                                                                    if i > 2:
                                                                        flag = 1
                                                                        break

                                                                if flag == 0:
                                                                    #                     g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                                                                    kf = sqrt(dot(dot(dot(ef, 2), m0), q) / abs(A) / (
                                                                        1 + g)) / h
                                                                    kxf = dot(dot(kf, sin_t), cos(phi))
                                                                    kyf = dot(dot(kf, sin_t), sin(phi))
                                                                    kzf = dot(kf, cos_t)
                                                                    valleyf = 4
                                                                else:
                                                                    valleyf = 3
                                                        else:
                                                            if r1 < scatGaAs_hole(11, ie, valley):
                                                                r = rand()
                                                                flag = 0
                                                                phi = dot(dot(2, pi), rand())
                                                                cos_t = 1 - dot(2, rand())
                                                                sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                                g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                    dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                        dot(sin_t ** 4, cos(phi) ** 2),
                                                                        sin(phi) ** 2)))) ** 0.5
                                                                fmax = 1 / (dot(a, (
                                                                    1 - sqrt((B / A) ** 2 + (C / A) ** 2 / 3)))) ** 1.5
                                                                f = 1 / (dot(a, (1 - g))) ** 1.5
                                                                i = 1
                                                                while r > f / fmax:

                                                                    phi = dot(dot(2, pi), rand())
                                                                    cos_t = 1 - dot(2, rand())
                                                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                                    g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                        dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                            dot(sin_t ** 4, cos(phi) ** 2),
                                                                            sin(phi) ** 2)))) ** 0.5
                                                                    f = 1 / (dot(a, (1 - g))) ** 1.5
                                                                    i = i + 1
                                                                    if i > 2:
                                                                        flag = 1
                                                                        break

                                                                if flag == 0:
                                                                    ef = ei
                                                                    #                 g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                                                                    kf = sqrt(dot(dot(dot(ef, 2), m0), q) / abs(A) / (
                                                                        1 - g)) / h
                                                                    kxf = dot(dot(kf, sin_t), cos(phi))
                                                                    kyf = dot(dot(kf, sin_t), sin(phi))
                                                                    kzf = dot(kf, cos_t)
                                                                    valleyf = 3
                                                                else:
                                                                    valleyf = 3
                                                            else:
                                                                if r1 < scatGaAs_hole(12, ie, valley):
                                                                    r = rand()
                                                                    flag = 0
                                                                    phi = dot(dot(2, pi), rand())
                                                                    cos_t = 1 - dot(2, rand())
                                                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                                    g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                        dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                            dot(sin_t ** 4, cos(phi) ** 2),
                                                                            sin(phi) ** 2)))) ** 0.5
                                                                    fmax = 1 / (dot(a, (1 + B / A))) ** 1.5
                                                                    f = 1 / (dot(a, (1 + g))) ** 1.5
                                                                    i = 1
                                                                    while r > f / fmax:

                                                                        phi = dot(dot(2, pi), rand())
                                                                        cos_t = 1 - dot(2, rand())
                                                                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                                        g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                            dot(dot(dot(sin_t, sin_t), cos_t),
                                                                                cos_t) + dot(
                                                                                dot(sin_t ** 4, cos(phi) ** 2),
                                                                                sin(phi) ** 2)))) ** 0.5
                                                                        f = 1 / (dot(a, (1 - g))) ** 1.5
                                                                        i = i + 1
                                                                        if i > 2:
                                                                            flag = 1
                                                                            break

                                                                    if flag == 0:
                                                                        ef = ei
                                                                        #                 g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                                                                        kf = sqrt(
                                                                            dot(dot(dot(ef, 2), m0), q) / abs(A) / (
                                                                                1 + g)) / h
                                                                        kxf = dot(dot(kf, sin_t), cos(phi))
                                                                        kyf = dot(dot(kf, sin_t), sin(phi))
                                                                        kzf = dot(kf, cos_t)
                                                                        valleyf = 4
                                                                    else:
                                                                        valleyf = 3
                                                                        # -----------------------------------
                else:
                    if valley == 4:
                        valley = valley - 2
                        if r1 < scatGaAs_hole(1, ie, valley):
                            ef = ei
                            phi = dot(dot(2, pi), rand())
                            cos_t = 1 - dot(2, rand())
                            sin_t = sqrt(1 - dot(cos_t, cos_t))
                            g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                  sin(phi) ** 2)))) ** 0.5
                            kf = sqrt(dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 - g)) / h
                            kxf = dot(dot(kf, sin_t), cos(phi))
                            kyf = dot(dot(kf, sin_t), sin(phi))
                            kzf = dot(kf, cos_t)
                            valleyf = 3
                        else:
                            if r1 < scatGaAs_hole(2, ie, valley):
                                ef = ei
                                phi = dot(dot(2, pi), rand())
                                cos_t = 1 - dot(2, rand())
                                sin_t = sqrt(1 - dot(cos_t, cos_t))
                                g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                    dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                      sin(phi) ** 2)))) ** 0.5
                                kf = sqrt(dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 + g)) / h
                                kxf = dot(dot(kf, sin_t), cos(phi))
                                kyf = dot(dot(kf, sin_t), sin(phi))
                                kzf = dot(kf, cos_t)
                                valleyf = 4
                            else:
                                if r1 < scatGaAs_hole(3, ie, valley):
                                    r = rand()
                                    ef = ei + hw0
                                    phi = dot(dot(2, pi), rand())
                                    cos_t = 1 - dot(2, rand())
                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                    g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                        dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                          sin(phi) ** 2)))) ** 0.5
                                    kf = sqrt(dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 + g)) / h
                                    kxf = dot(dot(kf, sin_t), cos(phi))
                                    kyf = dot(dot(kf, sin_t), sin(phi))
                                    kzf = dot(kf, cos_t)
                                    g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                        dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                          sin(phi) ** 2)))) ** 0.5
                                    cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                    G = dot(0.25, (1 + dot(3, cos_alpha ** 2)))
                                    cos_theta = kz / ki
                                    fa = dot(((1 + g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                        dot((1 + g2) / (1 + g1), ef) / ei))), G)
                                    fb = dot(sqrt(1 + g2), ((1 + g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                        dot((1 + g2) / (1 + g1), ef) / ei)) + dot(dot(dot(a, b), b), (1 + g2))) ** 2)
                                    f = fa / fb
                                    fmax = (1 + ef / ei - dot(dot(2, sqrt(ef / ei)), cos_theta)) / (1 + g100) ** 0.5 / \
                                           (1 + ef / ei - dot(dot(2, sqrt(ef / ei)), cos_theta) + dot(dot(dot(a, b), b),
                                                                                                      (1 + g100))) ** 2
                                    valleyf = 4
                                    i = 1
                                    while r > f / fmax:

                                        phi = dot(dot(2, pi), rand())
                                        cos_t = 1 - dot(2, rand())
                                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                                        g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                            dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                              sin(phi) ** 2)))) ** 0.5
                                        kf = sqrt(dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 + g)) / h
                                        kxf = dot(dot(kf, sin_t), cos(phi))
                                        kyf = dot(dot(kf, sin_t), sin(phi))
                                        kzf = dot(kf, cos_t)
                                        g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                            dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                              sin(phi) ** 2)))) ** 0.5
                                        cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                        G = dot(0.25, (1 + dot(3, cos_alpha ** 2)))
                                        fa = dot(((1 + g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                            dot((1 + g2) / (1 + g1), ef) / ei))), G)
                                        fb = dot(sqrt(1 + g2), ((1 + g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta),
                                                                                                    sqrt(dot(
                                                                                                        (1 + g2) / (
                                                                                                        1 + g1),
                                                                                                        ef) / ei)) + dot(
                                            dot(dot(a, b), b),
                                            (1 + g2))) ** 2)
                                        f = fa / fb
                                        fmax = (1 + ef / ei - dot(dot(2, sqrt(ef / ei)), cos_theta)) / (1 + g100) \
                                                                                                       ** 0.5 / (
                                                                                                                1 + ef / ei - dot(
                                                                                                                    dot(
                                                                                                                        2,
                                                                                                                        sqrt(
                                                                                                                            ef / ei)),
                                                                                                                    cos_theta) +
                                                                                                                dot(dot(
                                                                                                                    dot(
                                                                                                                        a,
                                                                                                                        b),
                                                                                                                    b),
                                                                                                                    (
                                                                                                                    1 + g100))) ** 2
                                        valleyf = 4
                                        i = i + 1
                                        if i > 2:
                                            kxf = kx
                                            kyf = ky
                                            kzf = kz
                                            break

                                else:
                                    if r1 < scatGaAs_hole(4, ie, valley):
                                        r = rand()
                                        ef = ei + hw0
                                        phi = dot(dot(2, pi), rand())
                                        cos_t = 1 - dot(2, rand())
                                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                                        g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                            dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                              sin(phi) ** 2)))) ** 0.5
                                        kf = sqrt(dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 - g)) / h
                                        kxf = dot(dot(kf, sin_t), cos(phi))
                                        kyf = dot(dot(kf, sin_t), sin(phi))
                                        kzf = dot(kf, cos_t)
                                        g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                            dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                              sin(phi) ** 2)))) ** 0.5
                                        cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                        G = dot(0.75, (1 - cos_alpha ** 2))
                                        cos_theta = kz / ki
                                        fa = dot(((1 - g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                            dot((1 - g2) / (1 + g1), ef) / ei))), G)
                                        fb = dot(sqrt(1 - g2), ((1 - g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta),
                                                                sqrt(dot((1 - g2) / (1 + g1), ef) / ei)) +
                                                                dot(dot(dot(a, b), b), (1 - g2))) ** 2)
                                        f = fa / fb
                                        fmaxa = dot(0.75, ((1 - g100) / (1 + g100) + (ef / ei) - dot(dot(2, cos_theta),
                                                sqrt(dot((1 - g100) / (1 + g100), ef) / ei))))
                                        fmaxb = dot(sqrt(1 - g111), (
                                            (1 - g111) / (1 + g111) + (ef / ei) - dot(dot(2, cos_theta), sqrt(
                                                dot((1 - g111) / (1 + g111), ef) / ei)) + dot(dot(dot(a, b), b),
                                                                                              (1 - g111))) ** 2)
                                        fmax = fmaxa / fmaxb
                                        i = 1
                                        valleyf = 3
                                        while r > f / fmax:

                                            phi = dot(dot(2, pi), rand())
                                            cos_t = 1 - dot(2, rand())
                                            sin_t = sqrt(1 - dot(cos_t, cos_t))
                                            g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                                  sin(phi) ** 2)))) ** 0.5
                                            kf = sqrt(dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 - g)) / h
                                            kxf = dot(dot(kf, sin_t), cos(phi))
                                            kyf = dot(dot(kf, sin_t), sin(phi))
                                            kzf = dot(kf, cos_t)
                                            g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                                  sin(phi) ** 2)))) ** 0.5
                                            cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                            G = dot(0.75, (1 - cos_alpha ** 2))
                                            fa = dot(((1 - g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                                dot((1 - g2) / (1 + g1), ef) / ei))), G)
                                            fb = dot(sqrt(1 - g2), (
                                                (1 - g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                                    dot((1 - g2) / (1 + g1), ef) / ei)) + dot(dot(dot(a, b), b),
                                                                                              (1 - g2))) ** 2)
                                            f = fa / fb
                                            fmaxa = dot(0.75, (
                                                (1 - g100) / (1 + g100) + (ef / ei) - dot(dot(2, cos_theta), sqrt(
                                                    dot((1 - g100) / (1 + g100), ef) / ei))))
                                            fmaxb = dot(sqrt(1 - g111), (
                                                (1 - g111) / (1 + g111) + (ef / ei) - dot(dot(2, cos_theta), sqrt(
                                                    dot((1 - g111) / (1 + g111), ef) / ei)) + dot(dot(dot(a, b), b),
                                                                                                  (1 - g111))) ** 2)
                                            fmax = fmaxa / fmaxb
                                            valleyf = 3
                                            i = i + 1
                                            if i > 2:
                                                kxf = kx
                                                kyf = ky
                                                kzf = kz
                                                valleyf = 4
                                                break

                                    else:
                                        if r1 < scatGaAs_hole(5, ie, valley):
                                            r = rand()
                                            ef = ei - hw0
                                            if ef < 0:
                                                valleyf = 4
                                            else:
                                                phi = dot(dot(2, pi), rand())
                                                cos_t = 1 - dot(2, rand())
                                                sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                    dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                                      sin(phi) ** 2)))) ** 0.5
                                                kf = sqrt(dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 + g)) / h
                                                kxf = dot(dot(kf, sin_t), cos(phi))
                                                kyf = dot(dot(kf, sin_t), sin(phi))
                                                kzf = dot(kf, cos_t)
                                                g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                    dot(sin_t ** 2, cos_t ** 2) + dot(dot(sin_t ** 4, cos(phi) ** 2),
                                                                                      sin(phi) ** 2)))) ** 0.5
                                                cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                                G = dot(0.25, (1 + dot(3, cos_alpha ** 2)))
                                                cos_theta = kz / ki
                                                fa = dot(((1 + g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                                    dot((1 + g2) / (1 + g1), ef) / ei))), G)
                                                fb = dot(sqrt(1 + g2), (
                                                    (1 + g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                                        dot((1 + g2) / (1 + g1), ef) / ei)) + dot(dot(dot(a, b), b),
                                                                                                  (1 + g2))) ** 2)
                                                f = fa / fb
                                                fmax = (1 + ef / ei - dot(dot(2, sqrt(ef / ei)), cos_theta)) / \
                                                       (1 + g100) ** 0.5 / (1 + ef / ei - dot(dot(2, sqrt(ef / ei)),
                                                                                              cos_theta) + dot(
                                                    dot(dot(a, b), b), (1 + g100))) ** 2
                                                i = 1
                                                valleyf = 4
                                                while r > f / fmax:

                                                    phi = dot(dot(2, pi), rand())
                                                    cos_t = 1 - dot(2, rand())
                                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                    g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                        dot(sin_t ** 2, cos_t ** 2) + dot(
                                                            dot(sin_t ** 4, cos(phi) ** 2),
                                                            sin(phi) ** 2)))) ** 0.5
                                                    kf = sqrt(dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 + g)) / h
                                                    kxf = dot(dot(kf, sin_t), cos(phi))
                                                    kyf = dot(dot(kf, sin_t), sin(phi))
                                                    kzf = dot(kf, cos_t)
                                                    g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                        dot(sin_t ** 2, cos_t ** 2) + dot(
                                                            dot(sin_t ** 4, cos(phi) ** 2),
                                                            sin(phi) ** 2)))) ** 0.5
                                                    cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                                    G = dot(0.25, (1 + dot(3, cos_alpha ** 2)))
                                                    fa = dot(((1 + g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta),
                                                                                                  sqrt(dot((1 + g2) / (
                                                                                                      1 + g1),
                                                                                                           ef) / ei))),
                                                             G)
                                                    fb = dot(sqrt(1 + g2), (
                                                        (1 + g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                                            dot((1 + g2) / (1 + g1), ef) / ei)) + dot(dot(dot(a, b), b),
                                                                                                      (1 + g2))) ** 2)
                                                    f = fa / fb
                                                    fmax = (1 + ef / ei - dot(dot(2, sqrt(ef / ei)), cos_theta)) / \
                                                           (1 + g100) ** 0.5 / (1 + ef / ei - dot(dot(2, sqrt(ef / ei)),
                                                                                                  cos_theta) + dot(
                                                        dot(dot(a, b), b), (1 + g100))) ** 2
                                                    valleyf = 4
                                                    i = i + 1
                                                    if i > 2:
                                                        kxf = kx
                                                        kyf = ky
                                                        kzf = kz
                                                        valleyf = 4
                                                        break

                                        else:
                                            if r1 < scatGaAs_hole(6, ie, valley):
                                                r = rand()
                                                ef = ei - hw0
                                                if ef < 0:
                                                    valleyf = 4
                                                else:
                                                    phi = dot(dot(2, pi), rand())
                                                    cos_t = 1 - dot(2, rand())
                                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                    g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                        dot(sin_t ** 2, cos_t ** 2) + dot(
                                                            dot(sin_t ** 4, cos(phi) ** 2),
                                                            sin(phi) ** 2)))) ** 0.5
                                                    kf = sqrt(dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 - g)) / h
                                                    kxf = dot(dot(kf, sin_t), cos(phi))
                                                    kyf = dot(dot(kf, sin_t), sin(phi))
                                                    kzf = dot(kf, cos_t)
                                                    g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                        dot(sin_t ** 2, cos_t ** 2) + dot(
                                                            dot(sin_t ** 4, cos(phi) ** 2),
                                                            sin(phi) ** 2)))) ** 0.5
                                                    cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz, kzf)) / ki / kf
                                                    G = dot(0.75, (1 - cos_alpha ** 2))
                                                    cos_theta = kz / ki
                                                    fa = dot(((1 - g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta),
                                                                                                  sqrt(dot((1 - g2) / (
                                                                                                      1 + g1),
                                                                                                           ef) / ei))),
                                                             G)
                                                    fb = dot(sqrt(1 - g2), (
                                                        (1 - g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                                            dot((1 - g2) / (1 + g1), ef) / ei)) + dot(dot(dot(a, b), b),
                                                                                                      (1 - g2))) ** 2)
                                                    f = fa / fb
                                                    fmaxa = dot(0.75, (
                                                        (1 - g100) / (1 + g100) + (ef / ei) - dot(dot(2, cos_theta),
                                                                                                  sqrt(
                                                                                                      dot((1 - g100) / (
                                                                                                          1 + g100),
                                                                                                          ef) / ei))))
                                                    fmaxb = dot(sqrt(1 - g111), (
                                                        (1 - g111) / (1 + g111) + (ef / ei) - dot(dot(2, cos_theta),
                                                                                                  sqrt(
                                                                                                      dot((1 - g111) / (
                                                                                                          1 + g111),
                                                                                                          ef) / ei)) + dot(
                                                            dot(dot(a, b), b),
                                                            (1 - g111))) ** 2)
                                                    fmax = fmaxa / fmaxb
                                                    i = 1
                                                    valleyf = 3
                                                    while r > f / fmax:

                                                        phi = dot(dot(2, pi), rand())
                                                        cos_t = 1 - dot(2, rand())
                                                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                        g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                            dot(sin_t ** 2, cos_t ** 2) + dot(
                                                                dot(sin_t ** 4, cos(phi) ** 2), sin(phi) ** 2)))) ** 0.5
                                                        kf = sqrt(dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 - g)) / h
                                                        kxf = dot(dot(kf, sin_t), cos(phi))
                                                        kyf = dot(dot(kf, sin_t), sin(phi))
                                                        kzf = dot(kf, cos_t)
                                                        g2 = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                            dot(sin_t ** 2, cos_t ** 2) + dot(
                                                                dot(sin_t ** 4, cos(phi) ** 2), sin(phi) ** 2)))) ** 0.5
                                                        cos_alpha = (dot(kx, kxf) + dot(ky, kyf) + dot(kz,
                                                                                                       kzf)) / ki / kf
                                                        G = dot(0.75, (1 - cos_alpha ** 2))
                                                        fa = dot(((1 - g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta),
                                                                                                      sqrt(dot(
                                                                                                          (1 - g2) / (
                                                                                                              1 + g1),
                                                                                                          ef) / ei))),
                                                                 G)
                                                        fb = dot(sqrt(1 - g2), (
                                                            (1 - g2) / (1 + g1) + ef / ei - dot(dot(2, cos_theta), sqrt(
                                                                dot((1 - g2) / (1 + g1), ef) / ei)) + dot(
                                                                dot(dot(a, b), b),
                                                                (1 - g2))) ** 2)
                                                        f = fa / fb
                                                        fmaxa = dot(0.75, (
                                                            (1 - g100) / (1 + g100) + (ef / ei) - dot(dot(2, cos_theta),
                                                                                                      sqrt(dot(
                                                                                                          (1 - g100) / (
                                                                                                              1 + g100),
                                                                                                          ef) / ei))))
                                                        fmaxb = dot(sqrt(1 - g111), (
                                                            (1 - g111) / (1 + g111) + (ef / ei) - dot(dot(2, cos_theta),
                                                                                                      sqrt(dot(
                                                                                                          (1 - g111) / (
                                                                                                              1 + g111),
                                                                                                          ef) / ei)) + dot(
                                                                dot(dot(a, b), b), (1 - g111))) ** 2)
                                                        fmax = fmaxa / fmaxb
                                                        valleyf = 3
                                                        i = i + 1
                                                        if i > 2:
                                                            kxf = kx
                                                            kyf = ky
                                                            kzf = kz
                                                            valleyf = 4
                                                            break

                                            else:
                                                if r1 < scatGaAs_hole(7, ie, valley):
                                                    r = rand()
                                                    flag = 0
                                                    phi = dot(dot(2, pi), rand())
                                                    cos_t = 1 - dot(2, rand())
                                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                    g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                        dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                            dot(sin_t ** 4, cos(phi) ** 2), sin(phi) ** 2)))) ** 0.5
                                                    fmax = 1 / (dot(a, (1 + B / A))) ** 1.5
                                                    f = 1 / (dot(a, (1 + g))) ** 1.5
                                                    i = 1
                                                    while r > f / fmax:

                                                        phi = dot(dot(2, pi), rand())
                                                        cos_t = 1 - dot(2, rand())
                                                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                        g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                            dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                dot(sin_t ** 4, cos(phi) ** 2), sin(phi) ** 2)))) ** 0.5
                                                        f = 1 / (dot(a, (1 + g))) ** 1.5
                                                        i = i + 1
                                                        if i > 2:
                                                            flag = 1
                                                            break

                                                    if flag == 0:
                                                        ef = ei + hw0
                                                        kf = sqrt(dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 - g)) / h
                                                        kxf = dot(dot(kf, sin_t), cos(phi))
                                                        kyf = dot(dot(kf, sin_t), sin(phi))
                                                        kzf = dot(kf, cos_t)
                                                        valleyf = 3
                                                    else:
                                                        valleyf = 4
                                                else:
                                                    if r1 < scatGaAs_hole(8, ie, valley):
                                                        r = rand()
                                                        ef = ei - hw0
                                                        if ef < 0:
                                                            valleyf = 4
                                                        else:
                                                            flag = 0
                                                            phi = dot(dot(2, pi), rand())
                                                            cos_t = 1 - dot(2, rand())
                                                            sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                            g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                    dot(sin_t ** 4, cos(phi) ** 2),
                                                                    sin(phi) ** 2)))) ** 0.5
                                                            fmax = 1 / (dot(a, (1 + B / A))) ** 1.5
                                                            f = 1 / (dot(a, (1 + g))) ** 1.5
                                                            i = 1
                                                            while r > f / fmax:

                                                                phi = dot(dot(2, pi), rand())
                                                                cos_t = 1 - dot(2, rand())
                                                                sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                                g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                    dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                        dot(sin_t ** 4, cos(phi) ** 2),
                                                                        sin(phi) ** 2)))) ** 0.5
                                                                f = 1 / (dot(a, (1 + g))) ** 1.5
                                                                i = i + 1
                                                                if i > 2:
                                                                    flag = 1
                                                                    break

                                                            if flag == 0:
                                                                kf = sqrt(
                                                                    dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 - g)) / h
                                                                kxf = dot(dot(kf, sin_t), cos(phi))
                                                                kyf = dot(dot(kf, sin_t), sin(phi))
                                                                kzf = dot(kf, cos_t)
                                                                valleyf = 3
                                                            else:
                                                                valleyf = 4
                                                    else:
                                                        if r1 < scatGaAs_hole(9, ie, valley):
                                                            r = rand()
                                                            flag = 0
                                                            phi = dot(dot(2, pi), rand())
                                                            cos_t = 1 - dot(2, rand())
                                                            sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                            g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                    dot(sin_t ** 4, cos(phi) ** 2),
                                                                    sin(phi) ** 2)))) ** 0.5
                                                            fmax = 1 / (dot(a, (1 + B / A))) ** 1.5
                                                            f = 1 / (dot(a, (1 + g))) ** 1.5
                                                            i = 1
                                                            while r > f / fmax:

                                                                phi = dot(dot(2, pi), rand())
                                                                cos_t = 1 - dot(2, rand())
                                                                sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                                g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                    dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                        dot(sin_t ** 4, cos(phi) ** 2),
                                                                        sin(phi) ** 2)))) ** 0.5
                                                                f = 1 / (dot(a, (1 + g))) ** 1.5
                                                                i = i + 1
                                                                if i > 2:
                                                                    flag = 1
                                                                    break

                                                            if flag == 0:
                                                                ef = ei + hw0
                                                                kf = sqrt(
                                                                    dot(dot(dot(ef, 2), m0), q) / abs(A) / (1 + g)) / h
                                                                kxf = dot(dot(kf, sin_t), cos(phi))
                                                                kyf = dot(dot(kf, sin_t), sin(phi))
                                                                kzf = dot(kf, cos_t)
                                                                valleyf = 4
                                                            else:
                                                                valleyf = 4
                                                        else:
                                                            if r1 < scatGaAs_hole(10, ie, valley):
                                                                ef = ei - hw0
                                                                if ef < 0:
                                                                    valleyf = 4
                                                                else:
                                                                    r = rand()
                                                                    flag = 0
                                                                    phi = dot(dot(2, pi), rand())
                                                                    cos_t = 1 - dot(2, rand())
                                                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                                    g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                        dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                            dot(sin_t ** 4, cos(phi) ** 2),
                                                                            sin(phi) ** 2)))) ** 0.5
                                                                    #                 a=h*h*abs(A)/(2*m0);
                                                                    fmax = 1 / (dot(a, (1 + B / A))) ** 1.5
                                                                    f = 1 / (dot(a, (1 + g))) ** 1.5
                                                                    i = 1
                                                                    while r > f / fmax:

                                                                        phi = dot(dot(2, pi), rand())
                                                                        cos_t = 1 - dot(2, rand())
                                                                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                                        g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                            dot(dot(dot(sin_t, sin_t), cos_t),
                                                                                cos_t) + dot(
                                                                                dot(sin_t ** 4, cos(phi) ** 2),
                                                                                sin(phi) ** 2)))) ** 0.5
                                                                        f = 1 / (dot(a, (1 + g))) ** 1.5
                                                                        i = i + 1
                                                                        if i > 2:
                                                                            flag = 1
                                                                            break

                                                                    if flag == 0:
                                                                        kf = sqrt(
                                                                            dot(dot(dot(ef, 2), m0), q) / abs(A) / (
                                                                                1 + g)) / h
                                                                        kxf = dot(dot(kf, sin_t), cos(phi))
                                                                        kyf = dot(dot(kf, sin_t), sin(phi))
                                                                        kzf = dot(kf, cos_t)
                                                                        valleyf = 4
                                                                    else:
                                                                        valleyf = 4
                                                            else:
                                                                if r1 < scatGaAs_hole(11, ie, valley):
                                                                    r = rand()
                                                                    flag = 0
                                                                    phi = dot(dot(2, pi), rand())
                                                                    cos_t = 1 - dot(2, rand())
                                                                    sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                                    g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                        dot(dot(dot(sin_t, sin_t), cos_t), cos_t) + dot(
                                                                            dot(sin_t ** 4, cos(phi) ** 2),
                                                                            sin(phi) ** 2)))) ** 0.5
                                                                    #             a=h*h*abs(A)/(2*m0);
                                                                    fmax = 1 / (dot(a, (1 + B / A))) ** 1.5
                                                                    f = 1 / (dot(a, (1 + g))) ** 1.5
                                                                    i = 1
                                                                    while r > f / fmax:

                                                                        phi = dot(dot(2, pi), rand())
                                                                        cos_t = 1 - dot(2, rand())
                                                                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                                        g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                            dot(dot(dot(sin_t, sin_t), cos_t),
                                                                                cos_t) + dot(
                                                                                dot(sin_t ** 4, cos(phi) ** 2),
                                                                                sin(phi) ** 2)))) ** 0.5
                                                                        f = 1 / (dot(a, (1 + g))) ** 1.5
                                                                        i = i + 1
                                                                        if i > 2:
                                                                            flag = 1
                                                                            break

                                                                    if flag == 0:
                                                                        ef = ei
                                                                        kf = sqrt(
                                                                            dot(dot(dot(ef, 2), m0), q) / abs(A) / (
                                                                                1 + g)) / h
                                                                        kxf = dot(dot(kf, sin_t), cos(phi))
                                                                        kyf = dot(dot(kf, sin_t), sin(phi))
                                                                        kzf = dot(kf, cos_t)
                                                                        valleyf = 4
                                                                    else:
                                                                        valleyf = 4
                                                                else:
                                                                    if r1 < scatGaAs_hole(12, ie, valley):
                                                                        r = rand()
                                                                        flag = 0
                                                                        phi = dot(dot(2, pi), rand())
                                                                        cos_t = 1 - dot(2, rand())
                                                                        sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                                        g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                            dot(dot(dot(sin_t, sin_t), cos_t),
                                                                                cos_t) + dot(
                                                                                dot(sin_t ** 4, cos(phi) ** 2),
                                                                                sin(phi) ** 2)))) ** 0.5
                                                                        #             a=h*h*abs(A)/(2*m0);
                                                                        fmax = 1 / (dot(a, (1 + B / A))) ** 1.5
                                                                        f = 1 / (dot(a, (1 + g))) ** 1.5
                                                                        i = 1
                                                                        while r > f / fmax:

                                                                            phi = dot(dot(2, pi), rand())
                                                                            cos_t = 1 - dot(2, rand())
                                                                            sin_t = sqrt(1 - dot(cos_t, cos_t))
                                                                            g = ((B / A) ** 2 + dot((C / A) ** 2, (
                                                                                dot(dot(dot(sin_t, sin_t), cos_t),
                                                                                    cos_t) + dot(
                                                                                    dot(sin_t ** 4, cos(phi) ** 2),
                                                                                    sin(phi) ** 2)))) ** 0.5
                                                                            f = 1 / (dot(a, (1 + g))) ** 1.5
                                                                            i = i + 1
                                                                            if i > 2:
                                                                                flag = 1
                                                                                break

                                                                        if flag == 0:
                                                                            ef = ei
                                                                            kf = sqrt(
                                                                                dot(dot(dot(ef, 2), m0), q) / abs(A) / (
                                                                                    1 - g)) / h
                                                                            kxf = dot(dot(kf, sin_t), cos(phi))
                                                                            kyf = dot(dot(kf, sin_t), sin(phi))
                                                                            kzf = dot(kf, cos_t)
                                                                            valleyf = 3
                                                                        else:
                                                                            valleyf = 4

    particle[1, 1] = kxf
    particle[1, 2] = kyf
    particle[1, 3] = kzf
