#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Content of the Tables from related documents

Sources:

* IAPWS R15-11(2011) Release on the IAPWS Formulation 2011 for the Thermal Conductivity of
Ordinary Water Substance

"""


class R15_11:
    """class with table data from IAPWS R15-11"""

    # reference temperature
    T_star = 647.096  # K

    # reference pressure
    p_star = 22.064  # MPa

    # reference density
    rho_star = 322.0  # kg / m³

    # reference thermal conductivity
    lambda_star = 1e-3  # W / K m

    # reference viscosity
    my_star = 1e-6  # Pa s

    # specific gas constant:
    R = 0.46151805  # kJ / kg K

    # Table 1: Coefficients L_k in Eq. (16) for λ_0 (T)
    Table1_L = [+2.443221e-3, +1.323095e-2, +6.770357e-3, -3.454586e-3, +4.096266e-4]

    # Table 2: Coefficients L_ij in Eq. (17) for λ_1 (T,ρ)
    Table2_L = [
        [+1.60397357, +2.33771842, +2.19650529, -1.21051378, -2.7203370],  # j = 0
        [-0.646013523, -2.78843778, -4.54580785, +1.60812989, +4.57586331],  # j = 1
        [+0.111443906, +1.53616167, +3.55777244, -0.621178141, -3.18369245],  # j = 2
        [+0.102997357, -0.463045512, -1.40944978, +0.0716373224, +1.1168348],  # j = 3
        [-0.0504123634, +0.0832827019, +0.275418278, +0, -0.19268305],  # j = 4
        [+0.00609859258, -0.00719201245, -0.0205938816, +0, +0.012913842],  # j = 5
    ]
    Table2_L_inv = [
        [+1.60397357, -0.646013523, +0.111443906, +0.102997357, -0.0504123634, +0.00609859258],  # i = 0
        [+2.33771842, -2.78843778, +1.53616167, -0.463045512, +0.0832827019, -0.00719201245],  # i = 1
        [+2.19650529, -4.54580785, +3.55777244, -1.40944978, +0.275418278, -0.0205938816],  # i = 2
        [-1.21051378, +1.60812989, -0.621178141, +0.0716373224, +0, +0],  # i = 3
        [-2.7203370, +4.57586331, -3.18369245, +1.1168348, -0.19268305, +0.012913842],  # i = 4
    ]

    # Table 3: Critical-region constants
    Table3_GAMMA = 177.8514
    Table3_inv_q_dash_D = 0.40  # in nm
    Table3_nu = 0.630
    Table3_gamma = 1.239
    Table3_xi_0 = 0.13  # in nm
    Table3_TAU_0 = 0.06
    Table3_T_dash_R = 1.5

    # Table 4: Sample points for computer-program verification of the correlating equation, Eq. (15). At these points, λ2 = 0.
    Table4_T = [298.15, 298.15, 298.15, 873.15]  # in K
    Table4_rho = [0.0, 998.0, 1200.0, 0.0]  # in kg / m³
    Table4_lambda = [18.4341883, 607.712868, 799.038144, 79.1034659]  # in mW / m K

    # Table 5: Sample points for computer-program verification of the correlating equation, Eq. (15), including the critical-enhancement contribution λ_2 . For all points, λ_0(647.35 K) = 51.5764797.
    Table5_T = [647.35, 647.35, 647.35, 647.35, 647.35, 647.35, 647.35, 647.35]
    Table5_rho = [1.0, 122.0, 222.0, 272.0, 322.0, 372.0, 422.0, 750.0]
    Table5_lambda_dash_1 = [1.0068497, 2.1445173, 3.4840736, 4.2233708, 4.9681953, 5.6961250, 6.3973429, 11.5870532]
    Table5_lambda_dash_2 = [0.0001300, 20.3162320, 188.091206, 540.133176, 1187.51354, 356.53333, 118.931062, 3.3419303]
    Table5_lambda = [51.9298924, 130.922885, 367.787459, 757.959776, 1443.75556, 650.319402, 448.883487, 600.961346]

    # Table 6: Coefficients Aij in Eq. (25) for ζ (T R , ρ )
    # TODO
