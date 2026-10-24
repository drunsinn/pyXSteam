#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IAPWS R12-08(2008)
Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance

"""

import sys
import math
import logging

from .tables import R12_08
from .IAPWS_R6 import R6_p_rhoT, eq_phi_r_delta, eq_phi_r_deltadelta

logger = logging.getLogger(__name__)


def _valid_rhoT(rho: float, T: float, *, allow_zero_density: bool = False) -> bool:
    density_is_valid = rho >= 0 if allow_zero_density else rho > 0
    return math.isfinite(rho) and math.isfinite(T) and density_is_valid and T > 0


def R12_my_dash_0(T_dash: float) -> float:
    """eq 11, viscosity in the dilute-gas limit"""
    if not math.isfinite(T_dash) or T_dash <= 0:
        return float("NaN")
    try:
        numerator = 100 * math.sqrt(T_dash)
        denominator = 0
        for i, H in enumerate(R12_08.Table1_H):
            denominator += H / T_dash**i
        return numerator / denominator
    except (OverflowError, ZeroDivisionError, ValueError):
        return float("NaN")


def R12_my_dash_1(rho_dash: float, T_dash: float) -> float:
    """eq 12, contribution to viscosity due to finite density"""
    if not math.isfinite(rho_dash) or rho_dash < 0 or not math.isfinite(T_dash) or T_dash <= 0:
        return float("NaN")
    try:
        coeff_t = (1 / T_dash) - 1
        sum_ij = 0.0
        H = R12_08.Table2_H()
        for j, _ in enumerate(H[0]):
            for i, _ in enumerate(H):
                if H[i][j] != 0.0:
                    sum_ij += math.pow(coeff_t, i) * H[i][j] * math.pow((rho_dash - 1), j)
        return math.exp(rho_dash * sum_ij)
    except (OverflowError, ZeroDivisionError, ValueError):
        return float("NaN")


def R12_helpereq_rhoT(rho: float, T: float):
    if not _valid_rhoT(rho, T, allow_zero_density=True):
        return float("NaN")
    try:
        delta = rho / R12_08.RHO_STAR
        tau = R12_08.T_STAR / T

        if delta == 1:
            logger.debug("change delta slightly to avoid div 0 at delta = 1.0")
            delta += sys.float_info.epsilon

        phi_r_delta = eq_phi_r_delta(tau, delta)
        phi_r_deltadelta = eq_phi_r_deltadelta(tau, delta)

        part_1 = delta * phi_r_delta + math.pow(delta, 2) * phi_r_deltadelta
        part_2 = 1 + delta * phi_r_delta

        return (1.0 + (part_1 / part_2)) * rho
    except (OverflowError, ZeroDivisionError, ValueError):
        return float("NaN")


def R12_xi(rho: float, T: float) -> float:
    """eq 20, correlation length ξ"""
    if not _valid_rhoT(rho, T):
        return float("NaN")

    try:
        return _R12_xi_valid(rho, T)
    except (OverflowError, ZeroDivisionError, ValueError):
        return float("NaN")


def _R12_xi_valid(rho: float, T: float) -> float:
    # T_b
    T_dash = T / R12_08.T_STAR
    # rho_b
    rho_dash = rho / R12_08.RHO_STAR

    # T_r
    T_R = R12_08.Table3_T_dash_R * R12_08.T_STAR

    # p_b, p_dash
    p_b = R6_p_rhoT(rho, T) / R12_08.P_STAR

    p_b_R = R6_p_rhoT(rho, T_R) / R12_08.P_STAR

    zeta_T = (rho_dash / p_b) * rho / R12_helpereq_rhoT(rho, T)
    zeta_T_R = (rho_dash / p_b_R) * rho / R12_helpereq_rhoT(rho, T_R)

    Delta_Chi_dash = rho_dash * (zeta_T - zeta_T_R * (R12_08.Table3_T_dash_R / T_dash))

    if Delta_Chi_dash >= 0:
        part_1 = Delta_Chi_dash / R12_08.Table3_Gamma_0
        part_2 = R12_08.Table3_ny / R12_08.Table3_gamma
        return R12_08.Table3_xi_0 * math.pow(part_1, part_2)
    return 0.0


def R12_Y(rho: float, T: float) -> float:
    if not _valid_rhoT(rho, T):
        return float("NaN")

    xi = R12_xi(rho, T)
    if not math.isfinite(xi):
        return float("NaN")

    q_C_xi = R12_08.Table3_q_C * xi
    q_D_xi = R12_08.Table3_q_D * xi

    if xi <= 0.3817016416:
        part_1 = (1 / 5) * q_C_xi * math.pow(q_D_xi, 5)
        part_2 = 1.0
        part_2 += -1 * q_C_xi
        part_2 += math.pow(q_C_xi, 2)
        part_2 += -1 * (765 / 504) * math.pow(q_D_xi, 2)
        res = part_1 * part_2
    else:
        Psi_D = math.acos(math.pow(1 + math.pow(R12_08.Table3_q_D, 2) * math.pow(xi, 2), -1 / 2))
        w = math.sqrt(math.fabs((q_C_xi - 1) / (q_C_xi + 1))) * math.tan(Psi_D / 2)
        if q_C_xi > 1:
            L = math.log((1 + w) / (1 - w))
        else:
            L = 2 * math.atan(math.fabs(w))

        part_1 = (1 / 12) * math.sin(3 * Psi_D)
        part_2 = -1 * (1 / (4 * q_C_xi)) * math.sin(2 * Psi_D)
        part_3 = (1 / (math.pow(q_C_xi, 2))) * (1 - (5 / 4) * math.pow(q_C_xi, 2)) * math.sin(Psi_D)
        part_4 = -math.pow(q_C_xi, -3) * (
            (1 - (3 / 2) * math.pow(q_C_xi, 2)) * Psi_D - math.pow(math.fabs(math.pow(q_C_xi, 2) - 1), 3 / 2) * L
        )
        res = part_1 + part_2 + part_3 + part_4
    return res


def R12_my_dash_2(rho: float, T: float) -> float:
    """eq 14, critical enhancement of the viscosity"""
    if not _valid_rhoT(rho, T, allow_zero_density=True):
        return float("NaN")
    if rho == 0:
        return 1.0
    if 645.91 < T < 650.77 or 245.8 < rho < 405.3:  # eq 13
        # outside these boundaries, the critical enhancement contributes less than 2% to the viscosity value
        Y = R12_Y(rho, T)
        return math.exp(R12_08.Table3_x_u * Y)
    return 1.0


def my_rhoT(rho: float, T: float, industrial_application: bool = True) -> float:
    """Return viscosity in Pa·s for finite T > 0 and density >= 0.

    Positive finite states outside the release's recommended validity range are
    evaluated by extrapolation, as described by IAPWS R12-08. Invalid numeric
    inputs return NaN.
    """
    if not _valid_rhoT(rho, T, allow_zero_density=True):
        return float("NaN")

    rho_dash = rho / R12_08.RHO_STAR
    T_dash = T / R12_08.T_STAR
    try:
        my_dash = R12_my_dash_0(T_dash) * R12_my_dash_1(rho_dash, T_dash)
        if not industrial_application:
            my_dash *= R12_my_dash_2(rho, T)
        viscosity = my_dash * R12_08.MU_STAR
    except (OverflowError, ZeroDivisionError, ValueError):
        return float("NaN")
    return viscosity if math.isfinite(viscosity) else float("NaN")
