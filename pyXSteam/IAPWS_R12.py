#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IAPWS R12-08(2008)
Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance
http://www.iapws.org/relguide/visc.pdf
"""
import math
import logging

from . import Constants

logger = logging.getLogger(__name__)


def _my_dash_0(T_dash: float) -> float:
    # viscosity in the dilute-gas limit, eq 11, page 5
    H = [1.67752, 2.20462, 0.6366564, -0.241605]
    numerator = 100 * math.sqrt(T_dash)
    denominator = 0
    for i in range(0, 4):
        denominator += H[i] / T_dash**i
    return numerator / denominator


def _my_dash_1(T_dash: float, rho_dash: float) -> float:
    # contribution to viscosity due to finite density, eq 12, page 5
    # first index is i, second index is j
    H = [[0.0 for j in range(0, 7)] for i in range(0, 6)]
    H[0][0] = 5.20094e-1
    H[1][0] = 8.50895e-2
    H[2][0] = -1.08374
    H[3][0] = -2.89555e-1
    H[0][1] = 2.22531e-1
    H[1][1] = 9.99115e-1
    H[2][1] = 1.88797
    H[3][1] = 1.26613
    H[5][1] = 1.20573e-1
    H[0][2] = -2.81378e-1
    H[1][2] = -9.06851e-1
    H[2][2] = -7.72479e-1
    H[3][2] = -4.89837e-1
    H[4][2] = -2.57040e-1
    H[0][3] = 1.61913e-1
    H[1][3] = 2.57399e-1
    H[0][4] = -3.25372e-2
    H[3][4] = 6.98452e-2
    H[4][5] = 8.72102e-3
    H[3][6] = -4.35673e-3
    H[5][6] = -5.93264e-4

    sum_t = 0
    for i in range(0, 6):  # i
        sum_rho = 0
        for j in range(0, 7):  # j
            sum_rho += H[i][j] * ((rho_dash - 1) ** j)
        sum_t += (((1 / T_dash) - 1) ** i) * sum_rho
    return math.exp(rho_dash * sum_t)


def _l(w: float, q_C: float, xi: float) -> float:
    # eq 18
    if (q_C * xi) > 1:
        return math.log((1 + w) / (1 - w))
    return 2 * math.atan(abs(w))


def _y(xi: float) -> float:
    inv_q_C = 1.9  # table 3, critical-region constants
    inv_q_D = 1.1  # table 3, critical-region constants
    q_C = 1 / inv_q_C
    q_D = 1 / inv_q_D

    if 0 <= xi <= 0.3817016416:
        coefficient_0 = (1 / 5) * q_C * xi * (q_D * xi) ** 5
        summand_0 = q_C * xi
        summand_1 = (q_C * xi) ** 2
        summand_2 = (765 / 504) * ((q_D * xi) ** 2)
        coefficient_1 = 1 - summand_0 + summand_1 - summand_2
        return coefficient_0 * coefficient_1

    psi_D = math.acos((1 + (q_D**2) * (xi**2)) ** (-1 / 2))  # ep 17

    summand_0 = 1 / 12 * math.sin(3 * psi_D)  # ep16, part 1

    summand_1 = (1 / (4 * q_C * xi)) * math.sin(2 * psi_D)  # ep16, part 2

    summand_2 = (1 / (q_C * xi) ** 2) * (1 - (5 / 4) * (q_C * xi) ** 2) * math.sin(psi_D)  # ep16, part 3

    w = abs((q_C * xi - 1) / (q_C * xi + 1)) ** (1 / 2) * math.tan(psi_D / 2)  # ep 19
    sum3_sub0 = (1 - (3 / 2 * (q_C * xi) ** 2)) * psi_D
    sum3_sub1 = abs((q_C * xi) ** 2 - 1) ** 3 / 2 * _l(w, q_C, xi)
    summand_3 = (1 / (q_C * xi) ** 3) * (sum3_sub0 - sum3_sub1)  # ep16, part 4

    y_total = summand_0 - summand_1 + summand_2 - summand_3  # ep 16

    logger.debug("value for Y: %f", y_total)

    return y_total


def _sigma(rho_dash: float, T_dash: float) -> float:
    # eq 21a
    # TODO
    logger.debug("rho_dash %f T_dash %f", rho_dash, T_dash)
    raise NotImplementedError("use derivative of density on pressure at constant temperature")
    return -1


def _delta_chi_dash(rho_dash: float, T_dash: float) -> float:
    T_dash_R = 1.5  # table 3, critical-region constants

    sigma_0 = _sigma(rho_dash, T_dash)
    sigma_1 = _sigma(rho_dash, T_dash_R)

    return rho_dash * (sigma_0 - sigma_1 * (T_dash_R / T_dash))


def _my_dash_2(T: float, rho: float, T_dash: float, rho_dash: float) -> float:
    if 645.91 < T < 650.77 and 245.8 < rho < 405.3:  # eq 13, page 6
        logger.debug("values for T and rho are within critical region")
        # critical region
        x_my = 0.068  # critical exponent for viscosity

        ny = 0.630  # table 3, critical-region constants
        gamma = 1.239  # table 3, critical-region constants
        xi_0 = 1.239  # table 3, critical-region constants
        gamma_0 = 0.06  # table 3, critical-region constants

        delta_chi_dash = _delta_chi_dash(rho_dash, T_dash)

        if delta_chi_dash < 0:
            delta_chi_dash = 0

        xi = xi_0 * (delta_chi_dash / gamma_0) ** (ny / gamma)

        my_dash_2 = math.exp(x_my * _y(xi))  # eq 14
        logger.debug("value for my_dash_2:%f", my_dash_2)
        return my_dash_2

    logger.debug("values for T and rho are outside of critical region, use my_dash_2=1.0")
    return 1.0  # section 2.8, page 8


def eq10(T: float, rho: float, industrial: bool = True) -> float:
    """ """
    logger.debug("input values T=%fK rho=%fkm/m^3" % (T, rho))
    # p_t = Constants.__TRIPLE_POINT_PRESSURE__
    # T_t = Constants.__TRIPLE_POINT_TEMPERATURE__

    # T_m = -1  # TODO: pressure-dependent melting temperature

    T_star = 647.096
    rho_star = 322.0
    # p_star = 22.064
    my_star = 1.00e-6

    T_dash = T / T_star
    rho_dash = rho / rho_star
    # p_dash = p / p_star

    # if 0 < p < p_t and T_t <= T <= 1173.15:
    #     pass
    # elif p_t <= p <= 300 and T_m <= T <= 1173.15:
    #     pass
    # elif 300 < p <= 350 and T_m <= T <= 873.15:
    #     pass
    # elif 350 < p <= 500 and T_m <= T <= 433.15:
    #     pass
    # elif 500 <= p <= 1000 and T_m <= T <= 373.15:
    #     pass
    # else:
    #     logger.warning("eq10 not valid for p=%f and T=%f", p, T)

    my_dash_0 = _my_dash_0(T_dash)
    my_dash_1 = _my_dash_1(T_dash, rho_dash)

    if industrial:  # section 3
        my_dash_2 = 1.0
    else:
        my_dash_2 = _my_dash_2(T, rho, T_dash, rho_dash)

    my_dash = my_dash_0 * my_dash_1 * my_dash_2  # eq 10

    logger.debug("calculated values µ_0=%f µ_1=%f µ_2=%f µ_dash=%f" % (my_dash_0, my_dash_1, my_dash_2, my_dash))

    my = my_dash * my_star  # eq8, value is in µPa*s
    # my = my * 10E6 # in Pa*s
    return my
