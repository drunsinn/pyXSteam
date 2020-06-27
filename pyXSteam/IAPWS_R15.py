#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import logging
from . import RegionSelection
from .Regions import Region1, Region2, Region3, Region4, Region5
from . import Constants

logger = logging.getLogger(__name__)
"""
IAPWS R15-11
Release on the IAPWS Formulation 2011 for the Thermal Conductivity of Ordinary Water Substance
http://www.iapws.org/relguide/ThCond.pdf
"""

def tc_ptrhocp_R15(p, T, rho, my, cp, cv):
    """Thermal conductivity as a function of preasure, temperature, density,
    viscosity, specific isobaric heat capacity an specific isochoric heat capacity

    Source: IAPWS R15-11
    Release on the IAPWS Formulation 2011 for the Thermal Conductivity of Ordinary Water Substance
    http://www.iapws.org/relguide/ThCond.pdf

    and

    http://twt.mpei.ac.ru/mcs/worksheets/iapws/wspTCPT.xmcd

    Args:
        p (float): density value
        T (float): temperature value
        rho (float): density value
        my (float): viscosity value
        cp (float): specific isobaric heat capacity value
        cv (float): specific isochoric heat capacity value

    Returns:
        my (float): viscosity or NaN if arguments are out of range
    """
    logger.debug("tc_ptrhocp_R15 input: p {} MPa, T {} K, rho {}, my {}, cp {}, vc {}".format(p, T, rho, my, cp, cv))

    T_dash = T / Constants.__CRITICAL_TEMPERATURE_IAPWS_R15_11__ # Page 3 / EQ 1 in [K] and Page 3 / EQ 7
    rho_dash = rho / Constants.__CRITICAL_DENSITY_IAPWS_R15_11__ # Page 3 / EQ 3 in [kg/m^3] and Page 3 / EQ 9
    c_p_dash = cp / Constants.__SPECIFIC_GAS_CONSTANT_IAPWS_R15_11__ # Page 3 / EQ 6 in [kJ/kg] and Page 3 / EQ 12
    my_dash = my / 1E-6 # Page 3 / EQ 5 in [Pas] and Page 3 / EQ 11

    if c_p_dash < 0 or c_p_dash > 1.0E13:
        c_p_dash = 1.0E13 # Page 12 / Footnote 2
        c_p = c_p_dash * Constants.__SPECIFIC_GAS_CONSTANT_IAPWS_R15_11__

    kappa = cp / cv # Page 3 / EQ 13
    tc_star = 1E-3 # Page 3 / EQ 4 in [W/Km]

    # Page 8 / Table 3
    Lambda = 177.8514
    q_D_dash = 0.40E-9 ** -1 # nm
    ny = 0.630
    gamma = 1.239
    xi_0 = 0.13E-9
    Gamma_0 = 0.06
    T_R_dash = 1.5

    #print(" - cp {} cv {} μ {} ρ {}".format(cp, cv, my, rho))
    #print(" - T_rel {} ρ_rel {} μ_rel {} kappa {}".format(T_dash, rho_dash, my_dash, kappa))

    tc0_L = (2.443221E-3, 1.323095E-2, 6.770357E-3, -3.454586E-3, 4.096266E-4) # Page 6 / Table 1
    sum = 0.0
    for k in range(0, 5):
        sum += tc0_L[k] / (T_dash ** k)
    tc0 = math.sqrt(T_dash) / sum # Page 5 / EQ 16
    # print(" - - lambda_0 {}".format(tc0))

    tc1_L = list() # Page 6 / Table 2
    tc1_L.append(( 1.60397357,     2.33771842,     2.19650529,  -1.21051378,  -2.7203370)) # j = 0
    tc1_L.append((-0.646013523,   -2.78843778,    -4.54580785,   1.60812989,   4.57586331)) # j = 1
    tc1_L.append(( 0.111443906,    1.53616167,     3.55777244,  -0.621178141, -3.18369245)) # j = 2
    tc1_L.append(( 0.102997357,   -0.463045512,   -1.40944978,   0.0716373224, 1.1168348)) # j = 3
    tc1_L.append((-0.0504123634,   0.0832827019,   0.275418278,  0.0,         -0.19268305)) # j = 4
    tc1_L.append(( 0.00609859258, -0.00719201245, -0.0205938816, 0.0,          0.012913842)) # j = 5
    tc1_sum1 = 0.0
    for i in range(0, 5):
        part_1 = ((1 / T_dash) - 1.0) ** i
        tc1_sum2 = 0.0
        for j in range(0, 6):
            tc1_sum2 += tc1_L[j][i] * (rho_dash - 1.0) ** j
        tc1_sum1 += part_1 * tc1_sum2
    tc1 = math.exp(rho_dash * tc1_sum1) # Page 6 / EQ 17

    #print(" - - lambda_1 {}".format(tc1))

    A = list()
    if rho_dash <= 0.310559006: # Page 12 / EQ 26
        A = (6.53786807199516, -5.61149954923348,  3.39624167361325,  -2.27492629730878, 10.2631854662709,  1.97815050331519) # j = 0
    elif rho_dash > 0.310559006 and rho_dash <= 0.776397516:
        A = (6.52717759281799, -6.30816983387575,  8.08379285492595,  -9.82240510197603, 12.1358413791395, -5.54349664571295) # j = 1
    elif rho_dash > 0.776397516 and rho_dash <= 1.242236025:
        A = (5.35500529896124, -3.96415689925446,  8.91990208918795, -12.0338729505790,  9.19494865194302, -2.16866274479712) # j = 2
    elif rho_dash > 1.242236025 and rho_dash <= 1.863354037:
        A = (1.55225959906681,  0.464621290821181, 8.93237374861479, -11.0321960061126,  6.16780999933360, -0.965458722086812) # j = 3
    else: # rho_dash > 1.863354037
        A = (1.11999926419994,  0.595748562571649, 9.88952565078920, -10.3255051147040,  4.66861294457414, -0.503243546373828) # j = 4
    zeta_TR_roh_sum = 0
    for i in range(0, 6):
        zeta_TR_roh_sum += A[i] * (rho_dash ** i)
    zeta_TR_roh = 1.0 / zeta_TR_roh_sum
    if zeta_TR_roh < 0 or zeta_TR_roh > 1.0E13:
        zeta_TR_roh = 1.0E13 # Page 12 / Footnote 2
    print(" - ρ_dash: {} ζ(T_R_dash,ρ_dash):{} ".format(rho_dash, zeta_TR_roh))

    delta_v_delta_T = wsp(p,T)
    delta_roh_delta_T = - roh**2 * delta_v_delta_T
    zeta_T_roh = delta_roh_delta_T * (p_ref / roh_ref)

    delta_Xi_dash = rho_dash * (zeta_T_roh - zeta_TR_roh * (T_R_dash / T_dash))
    if delta_Xi_dash < 0:
        delta_Xi_dash = 0

    xi = xi_0 * (delta_Xi_dash / Gamma_0) ** (ny / gamma) # Page 7 / EQ 22
    y = q_D_dash * xi
    if y < 1.2E-7: # Page 7 / EQ 21
        Z = 0.0
    else:
        Zpart1 = 2.0 / (math.pi * y)
        Zpart2 = (1.0 - (1.0/kappa)) * math.atan(y) + (1.0/kappa) * y
        Zpart3 = 1.0 - math.exp(-1.0 / ((1.0 / y) + (y**2.0 / (3.0 * rho_dash**2.0))))
        Z = Zpart1 * (Zpart2 - Zpart3) # Page 7 / EQ 19
    tc2 = Lambda * ( (rho_dash * c_p_dash * T_dash) / (my_dash) ) * Z # Page 7 / EQ 19

    #print(" - ξ {} Z {} λ0_dash {} λ1_dash {} λ2_dash {} ".format(xi, Z, tc0, tc1, tc2))

    tc_dash = tc0 * tc1 + tc2 # Page 5 / EQ 15

    tc = tc_dash * tc_star
    logger.debug("tc_ptrhocp_R15 result: λ {}".format(tc))

    return tc
