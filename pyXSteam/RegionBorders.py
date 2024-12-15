#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Section 4: Region Borders
"""
import logging

from .IF97_Tables import R3

logger = logging.getLogger(__name__)


def pB23_T(T: float) -> float:
    """
    calculate preasure from temperature for boundary between region 2 and 3

    Section 4.1 Boundary between region 2 and 3.

    Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties
    of Water and Steam 1997 Section 4 Auxiliary Equation for the Boundary
    between Regions 2 and 3

    Eq 5, Page 5

    :param T: temperature in [K]

    :return: preasure in [MPa]
    """
    return 348.05185628969 - 1.1671859879975 * T + 1.0192970039326e-03 * (T**2)


def TB23_p(p: float) -> float:
    """
    calculate temperature from preasure for boundary between region 2 and 3

    Section 4.1 Boundary between region 2 and 3.

    Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties
    of Water and Steam 1997 Section 4 Auxiliary Equation for the Boundary
    between Regions 2 and 3

    Eq 6, Page 6

    :param p: preasure in [MPa]

    :return: temperature in [K]
    """
    return 572.54459862746 + ((p - 13.91883977887) / 1.0192970039326e-03) ** 0.5


def hB13_s(s: float) -> float:
    """
    calculate enthalpy from specific entropy for border between region 1 and 3

    Section 4.3 Region boundary 1 to 3  & 3to2 as a functions of s

    Supplementary Release on Backward Equations ( ) , p h s for Region 3,
    Chapter 4.5 page 23.

    :param s: specific entropy in [kJ / (kg K)]

    :return: enthalpy in [kJ / kg]
    """
    Ii = R3.Sub_psh3_Table23_I
    Ji = R3.Sub_psh3_Table23_J
    ni = R3.Sub_psh3_Table23_n
    Sigma = s / 3.8
    eta = 0
    for i in range(0, 6):
        eta = eta + ni[i] * (Sigma - 0.884) ** Ii[i] * (Sigma - 0.864) ** Ji[i]
    return eta * 1700


def TB23_hs(h: float, s: float) -> float:
    """
    calculate temperature from specific entropy and enthalpy for border
    between region 2 and 3

    Section 4.3 Region boundary 1 to 3  & 3 to 2 as a functions of s

    Supplementary Release on Backward Equations () , p h s for Region 3,
    Chapter 4.6 page 25.

    :param h: enthalpy in [kJ / kg]
    :param s: Specific entropy in [kJ / (kg K)]

    :return: temperature in [K]
    """
    Ii = R3.Sub_psh3_Table25_I
    Ji = R3.Sub_psh3_Table25_J
    ni = R3.Sub_psh3_Table25_n
    Sigma = s / 5.3
    eta = h / 3000
    teta = 0
    for i in range(0, 25):
        teta = teta + ni[i] * (eta - 0.727) ** Ii[i] * (Sigma - 0.864) ** Ji[i]
    return teta * 900
