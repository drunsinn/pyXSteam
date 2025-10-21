#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Content of the Tables from related documents

Sources:

* IAPWS SR1-86(1992) Revised Supplementary Release on Satutation Properties of
Ordinary Water Substances

"""


class SR1_86:
    """class with table data from IAPWS SR1-86(1992)"""

    # Reference constants
    T_c = 647.096  # in K
    p_c = 22.064  # in MPa
    rho_c = 322  # in kg / m³
    alpha_0 = 1000  # in J / kg
    phi_0 = alpha_0 / T_c

    # coefficients for Eq 1: Vapour pressure
    a = [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502]

    # coefficients for Eq 2: density of saturated liquid
    b = [1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -6.74694450e5]

    # coefficients for Eq 3: density of saturated vapour
    c = [-2.03150240, -2.68302940, -5.38626492, -17.2991605, -44.7586581, -63.9201063]

    # coefficients for Eq 4 and 5: auxiliary equations
    d = [-5.65134998e-8, 2690.66631, 127.287297, -135.003439, 0.981825814]
    d_alpha = -1135.905627715
    d_phi = 2319.5246

    h_dash_t = 0.611786  # specific enthalpy of the liquid at the triple point in J / kg

    # Table 1: Thermodynamic property values calculated at three selected temperatures
    Table1_T = [273.16, 373.1243, 647.096]  # in K
    Table1_p = [611.657, 0.101325e6, 22.064e6]  # in Pa
    Table1_drho_dT = [44.436693, 3.616e3, 268e3]  # in Pa / K
    Table1_rho_dash = [999.789, 958.365, 322.0]  # in kg / m³
    Table1_rho_double_dash = [0.00485426, 0.597586, 322.0]  # in kg / m³
    Table1_alpha = [-11.529101, 417.65e3, 1548.0e3]  # in J / kg
    Table1_h_dash = [0.611786, 419.05e3, 2086.6e3]  # in J / kg
    Table1_h_double_dash = [2500.5e3, 2675.7e3, 2086.6e3]  # in J / kg
    Table1_phi = [-0.04, 1.303e3, 3.578e3]  # in J / kg K
    Table1_s_dash = [0.0, 1.307e3, 4.410e3]  # in J / kg K
    Table1_s_double_dash = [9.154e3, 7.355e3, 4.410e3]  # in J / kg K
