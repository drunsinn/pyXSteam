#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IAPWS R14-08(2011)
Pressure along the Melting and Sublimation Curves of Ordinary Water Substance
http://www.iapws.org/relguide/MeltSub2011.pdf
"""
import math
import logging

logger = logging.getLogger(__name__)

__TYPE_ICE_Ih__ = 1
__TYPE_ICE_III__ = 3
__TYPE_ICE_V__ = 5
__TYPE_ICE_VI__ = 6
__TYPE_ICE_VII__ = 7


def pmelt_T_iceIh(T):
    """
    EQ 1 / Melting pressure of ice Ih
    """
    T_star = 273.16
    p_star = 611.657e-6
    theta = T / T_star
    a = (0.119539337e7, 0.808183159e5, 0.333826860e4)
    b = (0.300000e1, 0.257500e2, 0.103750e3)
    temp_sum = 0
    for i in range(0, 3):
        temp_sum += a[i] * (1 - theta ** b[i])
    pi_melt = 1 + temp_sum
    return pi_melt * p_star


def pmelt_T_iceIII(T):
    """
    EQ 2 / Melting pressure of ice III
    """
    T_star = 251.165
    p_star = 208.566
    theta = T / T_star
    pi_melt = 1 - 0.299948 * (1.0 - theta**60)
    return pi_melt * p_star


def pmelt_T_iceV(T):
    """
    EQ 3 / Melting pressure of ice V
    """
    T_star = 256.164
    p_star = 350.1
    theta = T / T_star
    pi_melt = 1 - 1.18721 * (1.0 - theta**8)
    return pi_melt * p_star


def pmelt_T_iceVI(T):
    """
    EQ 4 / Melting pressure of ice VI
    """
    T_star = 273.31
    p_star = 632.4
    theta = T / T_star
    pi_melt = 1 - 1.07476 * (1.0 - theta**4.6)
    return pi_melt * p_star


def pmelt_T_iceVII(T):
    """
    EQ 5 / Melting pressure of ice VII
    """
    T_star = 355.0
    p_star = 2216.0
    theta = T / T_star
    p1 = 0.173683e1 * (1 - (theta**-1))
    p2 = 0.544606e-1 * (1 - (theta**5))
    p3 = 0.806106e-7 * (1 - (theta**22))
    pi_melt = math.exp(p1 - p2 + p3)
    return pi_melt * p_star


def psubl_T(T):
    """
    EQ 6 / Sublimation Pressure
    """
    T_star = 273.16
    p_star = 611.657e-6
    a = (-0.212144006e2, 0.273203819e2, -0.610598130e1)
    b = (0.333333333e-2, 0.120666667e1, 0.170333333e1)
    theta = T / T_star
    temp_sum = 0
    for i in range(0, 3):
        temp_sum += a[i] * theta ** b[i]
    pi_subl = math.exp((theta**-1) * temp_sum)
    return pi_subl * p_star
