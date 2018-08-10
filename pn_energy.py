# coding: utf-8
# Author: Zhongyang Zhang

# Generated with SMOP  0.41-beta
from pn_main_v6 import valley, particles
from const_val import *


# pn_energy.m

load('ig=50id=0.mat')
x_charge = matlabarray(np.zeros((max_particles, 1)))
# x_charge[max_particles, 1] = 0
result = matlabarray(np.zeros((max_particles, 1)))
# result[max_particles, 1] = 0
for n in arange(1, max_particles).reshape(-1):
    iv = valley(n, 1)
    if iv == 1 or iv == 2:
        kx = particles[n, 1]
        ky = particles[n, 2]
        kz = particles[n, 3]
        skx = dot(kx, kx)
        sky = dot(ky, ky)
        skz = dot(kz, kz)
        sk = abs(skx + sky + skz)
        ki = sqrt(sk)
        cos_theta = kz / ki
        sin_theta = sqrt(1 - cos_theta ** 2)
        cos_phi = kx / (dot(ki, sin_theta))
        sin_phi = ky / (dot(ki, sin_theta))
        g1 = ((B / A) ** 2 + dot((C / A) ** 2, (
            dot(sin_theta ** 2, cos_theta ** 2) + dot(dot(sin_theta ** 4, cos_phi ** 2), sin_phi ** 2)))) ** 0.5
        gk = dot(dot((dot(h, h) / (dot(2, eM(iv)))), sk), (1 / q))
        if iv == 1 or iv == 2:
            ei = (sqrt(1 + dot(dot(4, alpha(iv)), gk)) - 1) / (dot(2, alpha(iv)))
            particles[n, 7] = ei
        else:
            if iv == 3:
                ei = dot(dot(dot((dot(dot(h, h), abs(A)) / (dot(2, emR))), sk), (1 / q)), (1 - g1))
                particles[n, 7] = ei
            else:
                if iv == 4:
                    ei = dot(dot(dot((dot(dot(h, h), abs(A)) / (dot(2, emR))), sk), (1 / q)), (1 + g1))
                    particles[n, 7] = ei
        # -----connect position and energy together-----
        result[n, 1] = particles[n, 5]
        result[n, 2] = particles[n, 7]

i = 1
temp1 = result[:, 1] != 0
ar = []
for n in arange(1, max_particles).reshape(-1):
    if temp1[n, 1] == 1:
        ar.append(result[n, :])
        i = i + 1

# figure
# scatter(ar[:, 1], ar[:, 2])
