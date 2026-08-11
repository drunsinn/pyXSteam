# -*- coding: utf-8 -*-

import math
from .tables import R15_11, R12_08
from .IAPWS_R6 import R6_cp_rhoT, R6_cv_rhoT, R6_p_rhoT
from .IAPWS_R12 import R12_helpereq_rhoT, R12_xi, my_rhoT as R12_my_rhoT


def R15_11_correlation(T: float, rho: float) -> float:
    """R15-11 first partial result of correlating equation Eq 15
    :param T: temperature in [K]
    :param rho: density in [m³ / kg]
    """
    t_dash = T / R15_11.T_star  # Eq 7
    rho_dash = rho / R15_11.rho_star  # Eq 9

    # thermal conductivity in the dilute-gas limit
    sum_1 = 0.0
    for i, L in enumerate(R15_11.Table1_L):
        sum_1 += L / t_dash**i
    lambda_dash_0 = math.sqrt(t_dash) / sum_1  # Eq 16

    # contribution to thermal conductivity due to finite density
    sum_2 = 0.0
    for i, Li in enumerate(R15_11.Table2_L_inv):
        sub_sum_1 = 0.0
        for j, Lji in enumerate(Li):
            sub_sum_1 += Lji * (rho_dash - 1) ** j
        sum_2 += (1 / t_dash - 1) ** i * sub_sum_1
    lambda_dash_1 = math.exp(rho_dash * sum_2)  # Ep 17

    return lambda_dash_0 * lambda_dash_1  # Eq 15 (part)


def R15_11_Z(y: float, rho_dash: float, kappa: float) -> float:
    if y < 1.2e-7:
        return 0.0

    kappa_inv = 1.0 / kappa if kappa != 0.0 else 1.0
    denominator = (1.0 / y) + (y * y) / (3.0 * rho_dash) if rho_dash != 0 else (1.0 / y)
    return 2.0 / (math.pi * y) * (
        ((1.0 - kappa_inv) * math.atan(y) + kappa_inv * y)
        - (1.0 - math.exp(-1.0 / (denominator * denominator)))
    )


def R15_11_critical_enhancement(T: float, rho: float) -> float:
    """R15-11 second partial result of correlating equation Eq 15
    :param T: temperature in [K]
    :param rho: density in [m³ / kg]
    """
    t_dash = T / R15_11.T_star
    rho_dash = rho / R15_11.rho_star

    T_R = R15_11.Table3_T_dash_R * R15_11.T_star
    p_b = R6_p_rhoT(rho, T) / R15_11.p_star
    p_b_R = R6_p_rhoT(rho, T_R) / R15_11.p_star

    zeta_T = (rho_dash / p_b) * rho / R12_helpereq_rhoT(rho, T)
    zeta_T_R = (rho_dash / p_b_R) * rho / R12_helpereq_rhoT(rho, T_R)

    delta_chi = rho_dash * (zeta_T - zeta_T_R * (R15_11.Table3_T_dash_R / t_dash))
    if delta_chi <= 0.0:
        return 0.0

    xi = R12_xi(rho, T)
    if xi <= 0.0:
        return 0.0
    y = xi / R15_11.Table3_inv_q_dash_D

    cp = R6_cp_rhoT(rho, T) / R15_11.R
    cv = R6_cv_rhoT(rho, T) / R15_11.R
    mu = R12_my_rhoT(rho, T, industrial_application=True) / R12_08.MU_STAR
    kappa = cp / cv if cv != 0.0 else 1.0

    return R15_11.Table3_GAMMA * rho_dash * cp * t_dash / mu * R15_11_Z(y, rho_dash, kappa)


def tc_ptrho(p: float, T: float, rho: float) -> float:
    """R15-11 calculate thermal conductivity as a function of pressure, temperature and density.

    :param p: pressure in [MPa]
    :param T: temperature in [K]
    :param rho: density in [m³ / kg]

    :return: thermal conductivity in [W / (m K)]
    """
    return (R15_11_correlation(T, rho) + R15_11_critical_enhancement(T, rho)) / 1000.0
