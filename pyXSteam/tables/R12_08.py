#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Content of the Tables from related documents

Sources:

* IAPWS R12-08 Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance

"""


class R12_08:
    """class with table data from IAPWS R12-08"""

    T_STAR = 647.096  # reference temperature [K], eq 1
    RHO_STAR = 322.0  # reference density [kg / m^3], eq 2
    P_STAR = 22.064  # reference pressure [Pa], eq 3
    MU_STAR = 1.0e-6  # reference viscosity [Pa s], eq 4

    # table 1: Coefficients H_i for my_dash_0
    Table1_H = [
        1.67752,
        2.20462,
        0.6366564,
        -0.241605,
    ]

    @staticmethod
    def Table2_H():
        # table 2: Coefficients H_ij for my_dash_1
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
        return H

    # Table 3: Critical-Region Constants
    Table3_x_u = 0.068  # critical exponent for viscosity
    Table3_inv_q_C = 1.9  # nm
    Table3_inv_q_D = 1.1  # nm
    Table3_q_C = 1 / 1.9e-9
    Table3_q_D = 1 / 1.1e-9
    Table3_ny = 0.630
    Table3_gamma = 1.239
    Table3_xi_0 = 0.13  # nm
    Table3_Gamma_0 = 0.06
    Table3_T_dash_R = 1.5

    # Table 4: Sample points for computer-program verification of the correlating equation, Eq. (10), with µ_2 =1.
    Table4_T = [298.15, 298.15, 373.15, 433.15, 433.15, 873.15, 873.15, 873.15, 1173.15, 1173.15, 1173.15]
    Table4_rho = [998, 1200, 1000, 1, 1000, 1, 100, 600, 1, 100, 400]
    Table4_my = [
        889.735100e-6,
        1437.649467e-6,
        307.883622e-6,
        14.538324e-6,
        217.685358e-6,
        32.619287e-6,
        35.802262e-6,
        77.430195e-6,
        44.217245e-6,
        47.640433e-6,
        64.154608e-6,
    ]

    # Table 5: Sample points for computer-program verification of the correlating equation, Eq. (10), in the region near the critical point.
    Table5_T = [647.35, 647.35, 647.35, 647.35, 647.35, 647.35]  # in K
    Table5_rho = [122, 222, 272, 322, 372, 422]  # in kg / m^3
    Table5_xi = [
        0.309247,
        1.571405,
        5.266522,
        16.590209,
        5.603768,
        1.876244,
    ]  # in nm
    Table5_my2dash = [
        1.00000289,  # Correlation length ξ < 0.3817016416 nm so Y is evaluated with Eq. (15).
        1.00375120,
        1.03416789,
        1.09190440,
        1.03665871,
        1.00596332,
    ]
    Table5_my = [25.520677e-6, 31.337589e-6, 36.228143e-6, 42.961579e-6, 45.688204e-6, 49.436256e-6]  # Pa*s
