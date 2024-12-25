#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Section 4: Region Borders
"""
import logging

from .Tables import IF97, Sub_psh3

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
    return IF97.Table1_n[0] + IF97.Table1_n[1] * T + IF97.Table1_n[2] * (T**2)


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
    return IF97.Table1_n[3] + ((p - IF97.Table1_n[4]) / IF97.Table1_n[2]) ** 0.5


def hB13_s(s: float) -> float:
    """
    calculate enthalpy from specific entropy for border between region 1 and 3

    Section 4.3 Region boundary 1 to 3  & 3 to 2 as a functions of s

    Supplementary Release on Backward Equations ( ) , p h s for Region 3,
    Chapter 4.5 page 23.

    :param s: specific entropy in [kJ / (kg K)]

    :return: enthalpy in [kJ / kg]
    """
    Sigma = s / 3.8
    eta = 0
    for I, J, n in zip(Sub_psh3.Table23_I, Sub_psh3.Table23_J, Sub_psh3.Table23_n):
        eta = eta + n * (Sigma - 0.884) ** I * (Sigma - 0.864) ** J
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
    Sigma = s / 5.3
    eta = h / 3000
    teta = 0
    for I, J, n in zip(Sub_psh3.Table23_I, Sub_psh3.Table23_J, Sub_psh3.Table23_n):
        teta = teta + n * (eta - 0.727) ** I * (Sigma - 0.864) ** J
    return teta * 900
