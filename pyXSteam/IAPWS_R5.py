#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IAPWS R5-85(1994)
Revised Release on Surface Tension of Heavy Water Substance
http://www.iapws.org/relguide/surfd2o.pdf
"""
import math
import logging

from .Constants import TRIPLE_POINT_TEMPERATURE, CRITICAL_TEMPERATURE_D20_1992

logger = logging.getLogger(__name__)

def surface_tension_T(T: float) -> float:
    """
    IAPWS Release on Surface Tension of Heavy Water Substance
    http://www.iapws.org/relguide/surfd2o.pdf

    :param T: temperature in Kelvin

    :return: surface tension in mN/m
    """
    B = 0.2358 # N/m
    bb = -0.639
    my = 1.25
    if TRIPLE_POINT_TEMPERATURE < T < CRITICAL_TEMPERATURE_D20_1992:
        tau = 1 - T / CRITICAL_TEMPERATURE_D20_1992
        return B * tau**my * (1 + bb * tau)
    else:
        logger.warning("Temperature out of range of validity")
        return float("NaN")


