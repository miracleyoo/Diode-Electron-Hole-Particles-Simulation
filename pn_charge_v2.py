# coding: utf-8
# Author: Zhongyang Zhang

# Generated with SMOP  0.41-beta
from const_val import *


# pn_charge_v2.m


@function
def pn_charge_v2(particles=None, valley=None, nx1=None, ny1=None, dx=None, dy=None, max_particles=None, cpsp=None,
                 *args, **kwargs):
    varargin = pn_charge_v2.varargin
    nargin = pn_charge_v2.nargin

    charge_p = matlabarray(np.zeros((ny1, nx1)))
    charge_n = matlabarray(np.zeros((ny1, nx1)))

    # charge_p[ny1, nx1] = 0
    # charge_n[ny1, nx1] = 0
    temp_n = 0
    for n in arange(1, max_particles).reshape(-1):
        if valley(n) != 9:
            x = particles(n, 5) / dx
            y = particles(n, 6) / dy
            i = min(floor(y) + 1, ny1 - 1)
            j = min(floor(x) + 1, nx1 - 1)
            i = max(i, 1)
            j = max(j, 1)
            yb = i - 1
            xb = j - 1
            y1 = 1 - (y - yb)
            x1 = 1 - (x - xb)
            if valley(n) == 1 or valley(n) == 2:
                charge_n[i, j] = charge_n(i, j) + dot(x1, y1)
                charge_n[i + 1, j] = charge_n(i + 1, j) + dot(x1, (1 - y1))
                charge_n[i, j + 1] = charge_n(i, j + 1) + dot(y1, (1 - x1))
                charge_n[i + 1, j + 1] = charge_n(i + 1, j + 1) + dot((1 - x1), (1 - y1))
                temp_n = temp_n + 1
            else:
                if valley(n) == 3 or valley(n) == 4:
                    #             jj='n=';#914
                    #             fprintf('#s #f\n',jj,n);#914
                    charge_p[i, j] = charge_p(i, j) + dot(x1, y1)
                    charge_p[i + 1, j] = charge_p(i + 1, j) + dot(x1, (1 - y1))
                    charge_p[i, j + 1] = charge_p(i, j + 1) + dot(y1, (1 - x1))
                    charge_p[i + 1, j + 1] = charge_p(i + 1, j + 1) + dot((1 - x1), (1 - y1))
                    temp_n = temp_n + 1

    for i in arange(1, ny1).reshape(-1):
        for j in arange(1, nx1).reshape(-1):
            charge_p[i, j] = dot(charge_p(i, j), cpsp) / dx / dy
            charge_n[i, j] = dot(charge_n(i, j), cpsp) / dx / dy
            if i == 1 or i == ny1:
                charge_p[i, j] = dot(charge_p(i, j), 2)
                charge_n[i, j] = dot(charge_n(i, j), 2)
            if j == 1 or j == nx1:
                charge_p[i, j] = dot(charge_p(i, j), 2)
                charge_n[i, j] = dot(charge_n(i, j), 2)

    charge_p[ny1, nx1] = charge_p(ny1, nx1 - 1)

    charge_p[ny1, 1] = charge_p(ny1 - 1, 1)

    charge_n[ny1, nx1] = charge_n(ny1, nx1 - 1)

    charge_n[ny1, 1] = charge_n(ny1 - 1, 1)
