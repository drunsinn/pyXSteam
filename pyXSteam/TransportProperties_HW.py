#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""functions for transport properties of heavy water"""

import logging

from .Constants import TRIPLE_POINT_TEMPERATURE, CRITICAL_TEMPERATURE_D20_1992
from .Constants import __REFERENCE_TEMPERATURE_D20_R17_R18__ as T_STAR
from .Constants import __REFERENCE_PREASSURE_D20_R17_R18__ as P_STAR
from .Constants import __REFERENCE_DENSITY_D20_R17_R18__ as RHO_STAR
from .Constants import __REFERENCE_VISCOSITY_D20_R17_R18__ as MY_STAR
from .Constants import __REFERENCE_THERMAL_CONDUCTIVITY_D20_R18__ as LAMBDA_STAR
from .Constants import __SPECIFIC_GAS_CONSTANT_D20_R17_R18__ as R
from .Constants import __TRIPLE_POINT_TEMPERATURE_D2O_RESHW_2018__ as P_t
from .Constants import __TRIPLE_POINT_PRESSURE_D20_RESHW_2018__ as T_t

logger = logging.getLogger(__name__)


def my_AllRegions_pT() -> float:
    """R17-20 calculate viscosity of all liquide states for heavy water, replaces IAPWS R4-84(2007)"""
    # 0 < p ≤ pt pt ≤ p ≤ 100 MPa 100 MPa < p ≤ 200 MPa 200 MPa < p ≤ 960 MPa and and and and Tt ≤ T ≤ 775 K
    # Tm(p) ≤ T ≤ 775 K
    # Tm(p) ≤ T ≤ 473 K Tm(p) ≤ T ≤ 373 K
    pass


def tc_() -> float:
    """R18-21 calculate thermal conductivity of heavy water, replaces IAPWS R4-84(2007)"""
    pass


def surface_tension_T(T: float) -> float:
    """R5-85(1994) calculate surface tention as a function of temperature for heavy water

    IAPWS Release on Surface Tension of Heavy Water Substance
    http://www.iapws.org/relguide/surfd2o.pdf

    :param T: temperature in Kelvin

    :return: surface tension in mN/m
    """
    B = 238.0  # N/m
    bb = -0.639
    my = 1.25

    if TRIPLE_POINT_TEMPERATURE <= T <= CRITICAL_TEMPERATURE_D20_1992:
        tau = 1 - T / CRITICAL_TEMPERATURE_D20_1992
        sigma = B * tau**my * (1 + bb * tau)
        return sigma

    logger.warning("Temperature out of range of validity")
    return float("NaN")
