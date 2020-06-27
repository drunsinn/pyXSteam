#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import logging
from . import Constants

logger = logging.getLogger(__name__)
"""
IAPWS R4-84(2007)
Revised Release on Viscosity and Thermal Conductivity of Heavy Water Substance
http://www.iapws.org/relguide/TransD2O-2007.pdf
"""

def myHW_rhoT_R4(rho, T):
    """Viscosity as a function of density and temperature for heavy water
    substance

    Args:
        rho (float): density value for heavy water
        T (float): temperature value for heavy water

    Returns:
        my (float): viscosity or NaN if arguments are out of range
    """
    logger.debug("myHW_rhoT_R4 input: ρ {} kg/m^3, T {} K".format(rho, T))

    T_star = 643.847 # K
    rho_star = 358 # kg / m^3
    my_star = 55.2651 # µ Pa s

    T_dash = T / T_star
    rho_dash = rho / rho_star

    A = [1.0, 0.940695, 0.578377, -0.202044]
    sum_A_T = 0
    for i in range(0, 4):
        sum_A_T += A[i] / (T_dash ** i)
    my_0_dash = math.sqrt(T_dash) / sum_A_T

    B = list()
    B.append([ 0.4864192,  -0.2448372,  -0.8702035,  0.8716056,  -1.051126, 0.3458395]) # j= 0
    B.append([ 0.3509007,   1.315436,    1.297752,   1.353448,    0.0,      0.0]) # j= 1
    B.append([-0.2847572,  -1.037026,   -1.287846,   0.0,         0.0,     -0.02148229]) # j= 2
    B.append([ 0.07013759,  0.4660127,   0.2292075, -0.4857462,   0.0,      0.0]) # j= 3
    B.append([ 0.01641220, -0.02884911,  0.0,        0.1607171,   0.0,     -0.009603846]) # j= 4
    B.append([-0.01163815, -0.008239587, 0.0,        0.0,         0.0,      0.004559914]) # j= 5
    B.append([ 0.0,         0.0,         0.0,       -0.003886659, 0.0,      0.0]) # j= 6

    sum = 0
    for i in range(0, 6):
        part_T = ((1.0 / T_dash) - 1.0) ** i

        sum_B = 0
        for j in range(0, 7):
            sum_B += B[j][i] * (rho_dash - 1.0) ** j
        sum += part_T * sum_B

    my_1_dash = math.exp(rho_dash * sum)

    my_dash = my_0_dash * my_1_dash

    my = my_dash * my_star
    logger.debug("myHW_rhoT_R4 result: my {}".format(my))

    return my

def tcHW_rhoT_R4(rho, T):
    """Thermal conductivity as a function of density and temperature for heavy water
    substance

    Note: using tc instead of lambda to minimize the confusion with the python function

    Args:
        rho (float): density value for heavy water
        T (float): temperature value for heavy water

    Returns:
        λ (float): thermal conductivity or NaN if arguments are out of range
    """
    logger.debug("tcHW_rhoT_R4 input: ρ {} kg/m^3, T {} K".format(rho, T))

    T_star = 643.847 # K
    rho_star = 358 # kg / m^3
    tc_star = 0.742128 # mW/(m K)

    T_dash = T / T_star
    rho_dash = rho / rho_star

    A = [1.0, 37.3223, 22.5485, 13.0465, 0.0, -2.60735]
    B_e = -2.506
    B = [-167.31, 483.656, -191.039, 73.0358, -7.57467]
    C_1 = 35429.6
    C_2 = 5000.0E6
    C_T1 = 0.144847
    C_T2 = -5.64493
    C_R1 = -2.80000
    C_R2 = -0.080738543
    C_R3 = -17.9430
    rho_r1 = 0.125698
    D_1 = -741.112
    tau = T_dash / (math.fabs(T_dash - 1.1) + 1.1) # B15
    f_1 = math.exp(C_T1 * T_dash + C_T2 * (T_dash**2)) # B12
    f_2 = math.exp(C_R1 * (rho_dash-1.0)**2) + C_R2 * math.exp(C_R3*(rho_dash-rho_r1)**2) # B13
    f_3 = 1 + math.exp(60.0 * (tau-1.0) + 20.0) # B14
    f_4 = 1 + math.exp(100.0 * (tau-1.0) + 15.0) # B14
    part_C2 = (C_2 * f_1**4) / f_3
    part_f2 = (3.5 * f_2) / f_4

    # equation B8
    sum = 0
    for i in range(0, 6):
        sum += A[i] * (T_dash**i)
    tc_o = sum

    # equation B9
    sum = 0
    for i in range(1, 5):
        sum += B[i] * (rho_dash**i)
    delta_tc = B[0] * (1.0-math.exp(B_e * rho_dash)) + sum

    # equation B10
    delta_tc_c = C_1 * f_1 * f_2 * (1.0 + f_2**2 * (part_C2 + part_f2))

    # equation B11
    delta_tc_L = D_1 * f_1**1.2 * (1.0 - math.exp(-1.0 * (rho_dash/2.5)**10))

    # equation B7
    tc_dash = tc_o + delta_tc + delta_tc_c + delta_tc_L

    tc = tc_dash * tc_star
    logger.debug("tcHW_rhoT_R4 result: λ {}".format(tc))
    return tc
