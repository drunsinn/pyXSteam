#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Content of the Tables from related documents

Sources:

* IAPWS R14-08(2011) Revised Release on the Pressure along the Melting
and Sublimation Curves of Ordinary Water Substance

"""


class R14_08:
    class Ih:
        LOWER_BORDER = 273.16  # [K]
        UPPER_BORDER = 251.165  # [K]

        T_star = 273.16  # [K]
        p_star = 611.657e-6  # [MPa]
        coeff_a = [+0.119539337e7, +0.808183159e5, +0.333826860e4]
        coeff_b = [+0.300000e1, +0.257500e2, +0.103750e3]

    class III:
        LOWER_BORDER = 251.165  # [K]
        UPPER_BORDER = 256.164  # [K]

        T_star = 251.165  # [K]
        p_star = 208.566  # [MPa]
        coeff = 0.299948

    class V:
        LOWER_BORDER = 256.164  # [K]
        UPPER_BORDER = 273.31  # [K]

        T_star = 256.164  # [K]
        p_star = 350.1  # [MPa]
        coeff = 1.18721

    class VI:
        LOWER_BORDER = 273.31  # [K]
        UPPER_BORDER = 355.0  # [K]

        T_star = 273.31  # [K]
        p_star = 632.4  # [MPa]
        coeff = 1.07476

    class VII:
        LOWER_BORDER = 355.0  # [K]
        UPPER_BORDER = 715.0  # [K]

        T_star = 355.0  # [K]
        p_star = 2216.0  # [MPa]
        coeff_1 = 0.173683e1
        coeff_2 = 0.544606e-1
        coeff_3 = 0.806106e-7

    class Sublimation:
        LOWER_BORDER = 50.0  # [K]
        UPPER_BORDER = 273.16  # [K]

        T_star = 273.16  # [K]
        p_star = 611.657e-6  # [MPa]

        coeff_a = [-0.212144006e2, +0.273203819e2, -0.610598130e1]
        coeff_b = [+0.333333333e-2, +0.120666667e1, +0.170333333e1]
