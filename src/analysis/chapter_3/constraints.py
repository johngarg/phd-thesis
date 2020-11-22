#!/usr/bin/env python3

from workbench import CKM, SW2, HBAR, GF, kd

import rundec
import flavio
import numpy as np
from wilson import Wilson

CRD = rundec.CRunDec()
ASMZ = 0.1181
MZ = 91.1876
MT = 173.21
MB = 4.7
MC = 1.6
MW = 80.4
MMU = 0.106
MTAU = 1.777

ASMT5 = CRD.AlphasExact(ASMZ, MZ, MT, 5, 3)
ASMT6 = CRD.AlphasExact(ASMZ, MZ, MT, 6, 3)
ASMW = CRD.AlphasExact(ASMZ, MZ, MW, 5, 3)
ASMB = CRD.AlphasExact(ASMZ, MZ, 4.2, 5, 3)

ASMB4 = CRD.AlphasExact(ASMZ, MZ, MB, 4, 3)
ASMB5 = CRD.AlphasExact(ASMZ, MZ, MB, 5, 3)
ASMC3 = CRD.AlphasExact(ASMZ, MZ, MC, 3, 3)
ASMC4 = CRD.AlphasExact(ASMZ, MZ, MC, 4, 2)
# AS2 = CRD.AlphasExact(ASMZ, MZ, 2.0, 4, 2)


def beta0(nf):
    return 11 - 2 * nf / 3


def rknunu(m, x):
    rr = 1.91 / (m / 1000) ** 2
    term_2 = (2 * rr / 3) * (1 / np.conj(CKM[2, 1])) * (x.H * x)[1, 2]
    term_3 = (
        (rr ** 2 / 3) * (1 / np.abs(CKM[2, 1]) ** 2) * (x.H * x)[2, 2] * (x.H * x)[1, 1]
    )
    return 1 - np.real(term_2) + np.real(term_3)


def cbs(m, x):
    aem = 1 / 127
    asmlq = CRD.AlphasExact(ASMZ, MZ, m, 3)
    running = (ASMT5 / ASMW) ** (6 / 23) * (asmlq / ASMT6) ** (2 / 7)
    g = np.sqrt(4 * np.pi * aem) / np.sqrt(SW2)
    s0 = 2.3
    # only scanning over real values of couplings, so we don't care about phi_s
    # and can take the magnitude
    return np.abs(
        1
        + (running / (g ** 4 * s0))
        * (MW / m) ** 2
        * ((x.H * x)[1, 2] / np.conj(CKM[2, 1])) ** 2
    )


def Dsmunu(m, x, y):
    fDs = 249e-3
    tauDs = 500e-15
    z = x * CKM.H
    mDs = 1.9683
    mc = 1.27
    ms = 0.095
    inv_gamma_total = tauDs / HBAR

    gamma = -8
    asmlq = CRD.AlphasExact(ASMZ, MZ, m, 3)
    running_sd = (
        (ASMB4 / ASMC4) ** (gamma / (2 * beta0(4)))
        * (ASMT5 / ASMB5) ** (gamma / (2 * beta0(5)))
        * (asmlq / ASMT6) ** (gamma / (2 * beta0(6)))
    )

    prefactor = GF ** 2 * mDs * np.abs(CKM[1, 1]) ** 2 / (8 * np.pi)

    op_pref = 1 / (4 * np.sqrt(2) * GF * CKM[1, 1])
    ops = 0
    for i in (0, 1, 2):
        ops += (
            abs(
                kd(1, i)
                + op_pref * np.conj(z[1, 1]) * x[i, 1] / m ** 2
                - (mDs ** 2 / (MMU * (mc + ms)))
                * op_pref
                * running_sd
                * y[1, 1]
                * x[i, 1]
                / m ** 2
            )
            ** 2
        )
    return (
        inv_gamma_total
        * prefactor
        * fDs ** 2
        * MMU ** 2
        * (1 - MMU ** 2 / mDs ** 2) ** 2
        * ops
    )


def Bmunu(m, x, y):
    fB = 186e-3
    tauB = 1.638e-12
    z = x * CKM.H
    mu = 0.0022
    mb = 4.18
    mB = 5.26
    inv_gamma_total = tauB / HBAR

    asmlq = CRD.AlphasExact(ASMZ, MZ, m, 3)
    running_s = (ASMT5 / ASMB) ** (-12 / 23) * (asmlq / ASMT6) ** (-12 / 21)
    prefactor = GF ** 2 * mB * np.abs(CKM[0, 2]) ** 2 / (8 * np.pi)

    op_pref = 1 / (4 * np.sqrt(2) * GF * CKM[0, 2])
    ops = 0
    for i in (0, 1, 2):
        ops += (
            abs(
                kd(2, i)
                + op_pref * np.conj(z[2, 0]) * x[i, 2] / m ** 2
                - (mB ** 2 / (MMU * (mu + mb)))
                * op_pref
                * running_s
                * y[2, 0]
                * x[i, 2]
                / m ** 2
            )
            ** 2
        )
    return (
        inv_gamma_total
        * prefactor
        * fB ** 2
        * MMU ** 2
        * (1 - MMU ** 2 / mB ** 2) ** 2
        * ops
    )


def all_rdstar_ratios(m, x, y):
    """Function to RD(*) ratios with flavio"""
    GFVcb = 1.16638e-5 * 41.1e-3
    prefactor = 1 / (2 * np.sqrt(2) * GFVcb * m ** 2)

    z = x * CKM.H

    z32 = z[2, 1]
    z22 = z[1, 1]
    y32 = y[2, 1]
    y22 = y[1, 1]
    x13, x23, x33 = [float(i) for i in x[:, 2]]
    op_dict = {
        # tau
        "CVL_bctaunue": z32 * x13,
        "CVL_bctaunumu": z32 * x23,
        "CVL_bctaunutau": z32 * x33,
        "CSL_bctaunue": y32 * x13,
        "CSL_bctaunumu": y32 * x23,
        "CSL_bctaunutau": y32 * x33,
        # muon
        "CVL_bcmunue": z22 * x13,
        "CVL_bcmunumu": z22 * x23,
        "CVL_bcmunutau": z22 * x33,
        "CSL_bcmunue": y22 * x13,
        "CSL_bcmunumu": y22 * x23,
        "CSL_bcmunutau": y22 * x33,
    }

    op_dict["CT_bctaunue"] = -0.25 * op_dict["CSL_bctaunue"]
    op_dict["CT_bctaunumu"] = -0.25 * op_dict["CSL_bctaunumu"]
    op_dict["CT_bctaunutau"] = -0.25 * op_dict["CSL_bctaunutau"]
    op_dict["CT_bcmunue"] = -0.25 * op_dict["CSL_bcmunue"]
    op_dict["CT_bcmunumu"] = -0.25 * op_dict["CSL_bcmunumu"]
    op_dict["CT_bcmunutau"] = -0.25 * op_dict["CSL_bcmunutau"]

    wwet = Wilson(
        {k: v * prefactor for k, v in op_dict.items()},
        scale=m,
        eft="WET",
        basis="flavio",
    )

    rdstar = flavio.np_prediction("Rtaul(B->D*lnu)", wwet)
    rd = flavio.np_prediction("Rtaul(B->Dlnu)", wwet)

    rdmue = flavio.np_prediction("Rmue(B->Dlnu)", wwet)
    rdstarmue = flavio.np_prediction("Rmue(B->D*lnu)", wwet)

    return rd, rdstar, rdmue, rdstarmue


def rdstar_ratios(m, x, y):
    """Function to RD(*) ratios with flavio"""
    GFVcb = 1.16638e-5 * 41.1e-3
    prefactor = 1 / (2 * np.sqrt(2) * GFVcb * m ** 2)

    z = x * CKM.H

    z32 = z[2, 1]
    z22 = z[1, 1]
    y32 = y[2, 1]
    y22 = y[1, 1]
    x13, x23, x33 = [float(i) for i in x[:, 2]]
    op_dict = {
        # tau
        "CVL_bctaunue": z32 * x13,
        "CVL_bctaunumu": z32 * x23,
        "CVL_bctaunutau": z32 * x33,
        "CSL_bctaunue": y32 * x13,
        "CSL_bctaunumu": y32 * x23,
        "CSL_bctaunutau": y32 * x33,
    }

    op_dict["CT_bctaunue"] = -0.25 * op_dict["CSL_bctaunue"]
    op_dict["CT_bctaunumu"] = -0.25 * op_dict["CSL_bctaunumu"]
    op_dict["CT_bctaunutau"] = -0.25 * op_dict["CSL_bctaunutau"]

    wwet = Wilson(
        {k: v * prefactor for k, v in op_dict.items()},
        scale=m,
        eft="WET",
        basis="flavio",
    )

    rdstar = flavio.np_prediction("Rtaul(B->D*lnu)", wwet)
    rd = flavio.np_prediction("Rtaul(B->Dlnu)", wwet)

    return rd, rdstar
