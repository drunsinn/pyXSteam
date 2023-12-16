#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Section 4: Region Borders
"""
import logging

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
    Ii = [0, 1, 1, 3, 5, 6]
    Ji = [0, -2, 2, -12, -4, -3]
    ni = [
        0.913965547600543,
        -4.30944856041991e-05,
        60.3235694765419,
        1.17518273082168e-18,
        0.220000904781292,
        -69.0815545851641,
    ]
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
    Ii = [
        -12,
        -10,
        -8,
        -4,
        -3,
        -2,
        -2,
        -2,
        -2,
        0,
        1,
        1,
        1,
        3,
        3,
        5,
        6,
        6,
        8,
        8,
        8,
        12,
        12,
        14,
        14,
    ]
    Ji = [
        10,
        8,
        3,
        4,
        3,
        -6,
        2,
        3,
        4,
        0,
        -3,
        -2,
        10,
        -2,
        -1,
        -5,
        -6,
        -3,
        -8,
        -2,
        -1,
        -12,
        -1,
        -12,
        1,
    ]
    ni = [
        6.2909626082981e-04,
        -8.23453502583165e-04,
        5.15446951519474e-08,
        -1.17565945784945,
        3.48519684726192,
        -5.07837382408313e-12,
        -2.84637670005479,
        -2.36092263939673,
        6.01492324973779,
        1.48039650824546,
        3.60075182221907e-04,
        -1.26700045009952e-02,
        -1221843.32521413,
        0.149276502463272,
        0.698733471798484,
        -2.52207040114321e-02,
        1.47151930985213e-02,
        -1.08618917681849,
        -9.36875039816322e-04,
        81.9877897570217,
        -182.041861521835,
        2.61907376402688e-06,
        -29162.6417025961,
        1.40660774926165e-05,
        7832370.62349385,
    ]
    Sigma = s / 5.3
    eta = h / 3000
    teta = 0
    for i in range(0, 25):
        teta = teta + ni[i] * (eta - 0.727) ** Ii[i] * (Sigma - 0.864) ** Ji[i]
    return teta * 900
