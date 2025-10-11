# -*- coding: utf-8 -*-

import math
from .tables import R15_11


def R15_11_correlation(T: float, rho: float) -> float:
    """R15-11 first partial result of correlating equation Eq 15
    :param T: temperature in [K]
    :param rho: density in [m³ / kg]
    """
    t_dash = T / 647.26  # Eq 7
    rho_dash = rho / 322.0  # Eq 9

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


def R15_11_critical_enhancement(T: float, rho: float) -> float:
    """R15-11 second partial result of correlating equation Eq 15
    :param T: temperature in [K]
    :param rho: density in [m³ / kg]
    """
    # t_dash = T / 647.26  # Eq 7
    # rho_dash = rho / 322.0  # Eq 9
    # GAMMA = 177.8514
    # cp =
    # lambda_dash_2 = GAMMA * (rho_dash * )
    # return lambda_dash_2
    raise NotImplementedError()


def tc_ptrho(p: float, T: float, rho: float) -> float:
    """R15-11 calculate thermal conductivity as a function of preasure temperature and density

    Section 5.2 Thermal Conductivity (IAPWS formulation 1985)

    Revised release on the IAPWS formulation 1985 for the Thermal Conductivity of ordinary water IAPWS, September 1998

    :param p: preasure in [MPa]
    :param T: temperature in [K]
    :param rho: density in [m³ / kg]

    :return: surface tension in [mN/m]
    """
    raise NotImplementedError()

    # # ver2.6 Start corrected bug
    # if T < FREEZING_TEMPERATURE_H2O:
    #     logger.warning("Temperature out of range of validity")
    #     return float("NaN")
    # if T < 500 + FREEZING_TEMPERATURE_H2O:
    #     if p > 100:
    #         logger.warning("Preasure out of range of validity")
    #         return float("NaN")
    # if T <= 650 + FREEZING_TEMPERATURE_H2O:
    #     if p > 70:
    #         logger.warning("Preasure out of range of validity")
    #         return float("NaN")
    # else:  # T <= 800 + __FREEZING_POINT_H2O__:
    #     if p > 40:
    #         logger.warning("Preasure out of range of validity")
    #         return float("NaN")
    # # ver2.6 End corrected bug

    # t_dash = T / 647.26  # Eq 7
    # p_dash = p / 22.064  # Eq 8
    # rho_dash = rho / 322.0  # Eq 9
    # lambda_dash_1_2 = R15_11_correlation(T, rho)

    # tc0 = T**0.5 * (0.0102811 + 0.0299621 * T + 0.0156146 * (T**2) - 0.00422464 * (T**3))  # Page 9, Eq 9

    # tc1 = -0.397070 + 0.400302 * rho + 1.06 * math.exp(-0.171587 * ((rho + 2.392190) ** 2))  # Page 9, Eq 10

    # dT = abs(T - 1) + 0.00308976  # Page 9, Eq 12
    # Q = 2 + 0.0822994 / (dT ** (3 / 5))  # Page 10, Eq 13
    # if T >= 1:  # Page 10, Eq 14
    #     s = 1 / dT
    # else:
    #     s = 10.0932 / (dT ** (3 / 5))

    # tc2 = (
    #     (0.0701309 / (T**10) + 0.0118520) * (rho ** (9 / 5)) * math.exp(0.642857 * (1 - rho ** (14 / 5)))
    #     + 0.00169937 * s * (rho**Q) * math.exp((Q / (1 + Q)) * (1 - rho ** (1 + Q)))
    #     - 1.02 * math.exp(-4.11717 * (T ** (3 / 2)) - 6.17937 / (rho**5))
    # )  # Page 9, Eq 11
    # return tc0 + tc1 + tc2  # Page 9, Eq 8
