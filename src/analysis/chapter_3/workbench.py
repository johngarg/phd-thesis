#!/usr/bin/env python3

import numpy as np


def make_ckm():
    lam = 0.2257
    A = 0.814
    rho = 0.135
    eta = 0.349
    return np.matrix(
        [
            [1 - lam ** 2 / 2, lam, A * lam ** 3 * (rho - 1j * eta)],
            [-lam, 1 - lam ** 2 / 2, A * lam ** 2],
            [A * lam ** 3 * (1 - rho - 1j * eta), -A * lam ** 2, 1],
        ]
    )


def kd(i, j):
    if i == j:
        return 1
    return 0


CKM = make_ckm()
SW2 = 0.23142
VEV = 246
GF = 1.16638e-5
HBAR = 6.582119e-25  # GeV sec
