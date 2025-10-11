#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IAPWS R14-08(2011)
Revised Release on the Pressure along the Melting and Sublimation Curves of Ordinary Water Substance
http://www.iapws.org/relguide/MeltSub2011.pdf
"""
import math
import logging

from .tables import R14_08

logger = logging.getLogger(__name__)


def pmelt_T_iceIh(T: float) -> float:
    """R14-08(2011) calculate melting preassure for ice of type Ih, based on IAPWS R14-08(2011) EQ 1

    :param T: temperature in [K]

    :return: melting preasure in [MPa]
    """
    logger.debug("calculating 'melting preasure of ice type Ih' for T=%f", T)
    theta = T / R14_08.Ih.T_star
    pi_melt = 1
    for a, b in zip(R14_08.Ih.coeff_a, R14_08.Ih.coeff_b):
        pi_melt += a * (1 - theta**b)
    p_melt = pi_melt * R14_08.Ih.p_star
    logger.debug("result for 'melting preasure of ice type Ih': %f", p_melt)
    return p_melt


def pmelt_T_iceIII(T: float) -> float:
    """R14-08(2011) calculate melting preassure for ice of type III, based on IAPWS R14-08(2011) EQ 2

    :param T: temperature in [K]

    :return: melting preasure in [MPa]
    """
    logger.debug("calculating 'melting preasure of ice type III' for T=%f", T)
    theta = T / R14_08.III.T_star
    pi_melt = 1 - R14_08.III.coeff * (1.0 - theta**60)
    p_melt = pi_melt * R14_08.III.p_star
    logger.debug("result for 'melting preasure of ice type III': %f", p_melt)
    return p_melt


def Tmelt_p_iceIII(p: float) -> float:
    """inverse of `pmelt_T_iceIII`
    :param T: melting preasure in [MPa]

    :return: temperature in [K]
    """
    p_star = R14_08.III.p_star
    pi = p / p_star
    Tmelt = math.pow((pi - 1.0 + R14_08.III.coeff) / R14_08.III.coeff, 1 / 60) * R14_08.III.T_star
    return Tmelt


def pmelt_T_iceV(T: float) -> float:
    """R14-08(2011) calculate melting preassure for ice of type V, based on IAPWS R14-08(2011) EQ 3

    :param T: temperature in [K]

    :return: melting preasure in [MPa]
    """
    logger.debug("calculating 'melting preasure of ice type V' for T=%f", T)
    theta = T / R14_08.V.T_star
    pi_melt = 1 - R14_08.V.coeff * (1.0 - theta**8)
    p_melt = pi_melt * R14_08.V.p_star
    logger.debug("result for 'melting preasure of ice type V': %f", p_melt)
    return p_melt


def Tmelt_p_iceV(p: float) -> float:
    """inversee of `pmelt_T_iceV`
    :param T: melting preasure in [MPa]

    :return: temperature in [K]
    """
    pi = p / R14_08.V.p_star
    Tmelt = math.pow((pi - 1.0 + R14_08.V.coeff) / R14_08.V.coeff, 1 / 8) * R14_08.V.T_star
    return Tmelt


def pmelt_T_iceVI(T: float) -> float:
    """R14-08(2011) calculate melting preassure for ice of type VI, based on IAPWS R14-08(2011) EQ 4

    :param T: temperature in [K]

    :return: melting preasure in [MPa]
    """
    logger.debug("calculating 'melting preasure of ice type VI' for T=%f", T)
    theta = T / R14_08.VI.T_star
    pi_melt = 1 - R14_08.VI.coeff * (1.0 - theta**4.6)
    p_melt = pi_melt * R14_08.VI.p_star
    logger.debug("result for 'melting preasure of ice type VI': %f", p_melt)
    return p_melt


def Tmelt_p_iceVI(p: float) -> float:
    """inversee of `pmelt_T_iceVI`
    :param T: melting preasure in [MPa]

    :return: temperature in [K]
    """
    pi = p / R14_08.VI.p_star
    Tmelt = math.pow((pi - 1.0 + R14_08.VI.coeff) / R14_08.VI.coeff, 1 / 4.6) * R14_08.VI.T_star
    return Tmelt


def pmelt_T_iceVII(T: float) -> float:
    """R14-08(2011) calculate melting preassure for ice of type VII, based on IAPWS R14-08(2011) EQ 5

    :param T: temperature in [K]

    :return: melting preasure in [MPa]
    """
    logger.debug("calculating 'melting preasure of ice type VII' for T=%f", T)
    theta = T / R14_08.VII.T_star
    p1 = R14_08.VII.coeff_1 * (1 - (theta**-1))
    p2 = R14_08.VII.coeff_2 * (1 - (theta**5))
    p3 = R14_08.VII.coeff_3 * (1 - (theta**22))
    pi_melt = math.exp(p1 - p2 + p3)
    p_melt = pi_melt * R14_08.VII.p_star
    logger.debug("result for 'melting preasure of ice type VII': %f", p_melt)
    return p_melt


def psubl_T(T: float) -> float:
    """R14-08(2011) calculate sublimation preassure for ice of type, based on IAPWS R14-08(2011) EQ 6

    :param T: temperature in [K]

    :return: sublimation preasure in [MPa]
    """
    logger.debug("calculating 'sublimation preasure of ice' for T=%f", T)
    theta = T / R14_08.Sublimation.T_star
    temp_sum = 0
    for a, b in zip(R14_08.Sublimation.coeff_a, R14_08.Sublimation.coeff_b):
        temp_sum += a * theta**b
    pi_subl = math.exp((theta**-1) * temp_sum)
    p_subl = pi_subl * R14_08.Sublimation.p_star
    logger.debug("result for 'sublimation preasure of ice': %f", p_subl)
    return p_subl
