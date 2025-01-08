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
