#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Section 3: Region Selection
"""
import math
import logging
from .RegionBorders import pB23_T, TB23_p, hB13_s, TB23_hs
from .Regions import Region1, Region2, Region3, Region4, Region5
from .Constants import (
    CRITICAL_TEMPERATURE,
    TRIPLE_POINT_PRESSURE,
    FREEZING_TEMPERATURE_H2O,
    DiagramRegion,
)

logger = logging.getLogger(__name__)


def select_region_pT(p: float, T: float) -> DiagramRegion:
    """
    Section 3.1 Regions as a function of p and T
    Select diagram region based on values of preasure and temperature
    Returns `DiagramRegion.NILL` if arguments lead to no valid region

    :param p: preasure in [MPa]
    :param T: temperature in [K]

    :return: diagram region
    """
    if 0.000611 < p < 50.0 and 1073.15 < T < 2273.15:
        return DiagramRegion.R5
    if FREEZING_TEMPERATURE_H2O < T <= 1073.15 and 0.000611 < p <= 100:
        if T > 623.15:
            if p > pB23_T(T):
                region_pT_number = DiagramRegion.R3
                if T < CRITICAL_TEMPERATURE:
                    ps = Region4.p4_T(T)
                    if math.fabs(p - ps) < 0.00001:
                        region_pT_number = DiagramRegion.R4
            else:
                region_pT_number = DiagramRegion.R2
        else:
            ps = Region4.p4_T(T)
            if math.fabs(p - ps) < 0.00001:
                region_pT_number = DiagramRegion.R4
            elif p > ps:
                region_pT_number = DiagramRegion.R1
            else:
                region_pT_number = DiagramRegion.R2
    else:
        logger.warning("Temperature outside valid area")
        region_pT_number = DiagramRegion.NILL

    return region_pT_number


def select_region_ph(p: float, h: float) -> DiagramRegion:
    """
    Section 3.2 Regions as a function of p and h
    Select diagram region based on values of preasure and enthalpy
    Returns DiagramRegion.NILL if arguments lead to no valid region

    :param p: preasure in [MPa]
    :param h: enthalpy in [kJ / kg]

    :return: diagram region
    """
    # Check if outside pressure limits
    if (p < TRIPLE_POINT_PRESSURE) or (p > 100):
        logger.warning("Preasure outside valid area")
        return DiagramRegion.NILL

    # Check if outside low h.
    # Linear adaption to h1_pt()+2 to speed up calcualations.
    if h < (0.963 * p + 2.2):
        if h < Region1.h1_pT(p, FREEZING_TEMPERATURE_H2O):
            logger.warning("Enthalpy outside valid area")
            return DiagramRegion.NILL

    if p < 16.5292:  # Below region 3, Check region 1, 4, 2, 5
        # Check Region 1
        Ts = Region4.T4_p(p)
        hL = 109.6635 * math.log(p) + 40.3481 * p + 734.58  # Approximate function for hL_p
        if math.fabs(h - hL) < 100:  # if approximate is not god enough use real function
            hL = Region1.h1_pT(p, Ts)
        if h <= hL:
            return DiagramRegion.R1
        # Check Region 4
        hV = 45.1768 * math.log(p) - 20.158 * p + 2804.4  # Approximate function for hV_p
        if math.fabs(h - hV) < 50:  # if approximate is not god enough use real function
            hV = Region2.h2_pT(p, Ts)
        if h < hV:
            return DiagramRegion.R4
        # Check upper limit of region 2 Quick Test
        if h < 4000:
            return DiagramRegion.R2
        # Check region 2 (Real value)
        h_45 = Region2.h2_pT(p, 1073.15)
        if h <= h_45:
            return DiagramRegion.R2
        # Check region 5
        if p > 10:
            logger.warning("Preasure outside valid area")
            return DiagramRegion.NILL
        h_5u = Region5.h5_pT(p, 2273.15)
        if h < h_5u:
            return DiagramRegion.R5
        logger.warning("Enthalpy outside valid area")
        return DiagramRegion.NILL
    # for p > 16.5292
    # Check if in region1
    if h < Region1.h1_pT(p, 623.15):
        return DiagramRegion.R1
    # Check if in region 3 or 4 (Below Reg 2)
    if h < Region2.h2_pT(p, TB23_p(p)):
        # Region 3 or 4
        if p > Region3.psat3_h(h):
            return DiagramRegion.R3
        return DiagramRegion.R4
    # Check if region 2
    if h < Region2.h2_pT(p, 1073.15):
        return DiagramRegion.R2
    logger.warning("Preasure outside valid area")
    return DiagramRegion.NILL


def select_region_ps(p: float, s: float) -> DiagramRegion:
    """
    Section 3.3 Regions as a function of p and s
    Select diagram region based on values of preasure and specific entropy
    Returns `DiagramRegion.NILL` if arguments lead to no valid region

    :param p: preasure in [MPa]
    :param s: specific entropy in [kJ / (kg K)]

    :return: diagram region
    """
    if (p < TRIPLE_POINT_PRESSURE) or (p > 100) or (s < 0) or (s > Region5.s5_pT(p, 2273.15)):
        logger.warning("Preasure or Entropy outside valid area")
        return DiagramRegion.NILL
    # Check region 5
    if s > Region2.s2_pT(p, 1073.15):
        if p <= 10:
            return DiagramRegion.R5
        logger.warning("Preasure outside valid area")
        return DiagramRegion.NILL
    # Check region 2
    if p > 16.529:
        # Between 5.047 & 5.261. Use to speed up !
        ss = Region2.s2_pT(p, TB23_p(p))
    else:
        ss = Region2.s2_pT(p, Region4.T4_p(p))
    if s > ss:
        return DiagramRegion.R2
    # Check region 3
    ss = Region1.s1_pT(p, 623.15)
    if (p > 16.529) and (s > ss):
        if p > Region3.psat3_s(s):
            return DiagramRegion.R3
        return DiagramRegion.R4
    # Check region 4 (Not inside region 3)
    if (p < 16.529) and (s > Region1.s1_pT(p, Region4.T4_p(p))):
        return DiagramRegion.R4
    # Check region 1
    if (p > TRIPLE_POINT_PRESSURE) and (s > Region1.s1_pT(p, FREEZING_TEMPERATURE_H2O)):
        return DiagramRegion.R1
    # ToDo: Check if Defaulting to region 1 is correct here
    return DiagramRegion.R1


def select_region_hs(h: float, s: float) -> DiagramRegion:
    """
    Section 3.4 Regions as a function of h and s
    Select diagram region based on values of enthalpy and specific entropy
    Returns `DiagramRegion.NILL` if arguments lead to no valid region

    :param h: enthalpy in [kJ / kg]
    :param s: specific entropy in [kJ / (kg K)]

    :return: diagram region
    """
    if s < -0.0001545495919:
        logger.warning("Entropy outside valid area")
        return DiagramRegion.NILL
    # Check linear adaption to p=0.000611. if below region 4.
    hMin = ((-0.0415878 - 2500.89262) / (-0.00015455 - 9.155759)) * s
    if (s < 9.155759395) and (h < hMin):
        logger.warning("Entalpy or Entropy outside valid area")
        return DiagramRegion.NILL
    # Kolla 1 eller 4. (+liten bit ???ver B13)
    if -0.0001545495919 <= s <= 3.77828134:
        if h < Region4.h4_s(s):
            return DiagramRegion.R4
        if s < 3.397782955:  # 100MPa line is limiting
            TMax = Region1.T1_ps(100, s)
            hMax = Region1.h1_pT(100, TMax)
            if h < hMax:
                return DiagramRegion.R1
            logger.warning("Entalpy outside valid area")
            return DiagramRegion.NILL
        # The point is either in region 4,1,3. Check B23
        hB = hB13_s(s)
        if h < hB:
            return DiagramRegion.R1
        TMax = Region3.T3_ps(100, s)
        vmax = Region3.v3_ps(100, s)
        hMax = Region3.h3_rhoT(1 / vmax, TMax)
        if h < hMax:
            return DiagramRegion.R3
        logger.warning("Entalpy outside valid area")
        return DiagramRegion.NILL

    # Kolla region 2 eller 4. (???vre delen av omr???de b23-> max)
    if 5.260578707 <= s <= 11.9212156897728:
        if s > 9.155759395:  # Above region 4
            Tmin = Region2.T2_ps(0.000611, s)
            hMin = Region2.h2_pT(0.000611, Tmin)
            # function adapted to h(1073.15,s)
            hMax = -0.07554022 * s**4 + 3.341571 * s**3 - 55.42151 * s**2 + 408.515 * s + 3031.338
            if hMin < h < hMax:
                return DiagramRegion.R2
            logger.warning("Entalpy outside valid area")
            return DiagramRegion.NILL
        hV = Region4.h4_s(s)
        if h < hV:  # Region 4. Under region 3.
            return DiagramRegion.R4
        if s < 6.04048367171238:
            TMax = Region2.T2_ps(100, s)
            hMax = Region2.h2_pT(100, TMax)
        else:
            # function adapted to h(1073.15,s)
            hMax = -2.988734 * s**4 + 121.4015 * s**3 - 1805.15 * s**2 + 11720.16 * s - 23998.33
        if h < hMax:  # Region 2. ???ver region 4.
            return DiagramRegion.R2
        logger.warning("Entalpy outside valid area")
        return DiagramRegion.NILL
    # Check region 3 or 4 below the critical point.
    if 3.77828134 <= s <= 4.41202148223476:
        hL = Region4.h4_s(s)
        if h < hL:
            return DiagramRegion.R4
        TMax = Region3.T3_ps(100, s)
        vmax = Region3.v3_ps(100, s)
        hMax = Region3.h3_rhoT(1 / vmax, TMax)
        if h < hMax:
            return DiagramRegion.R3
        logger.warning("Entalpy outside valid area")
        return DiagramRegion.NILL
    # Check region 3 or 4 from the critical point to the upper part of B23
    if 4.41202148223476 <= s <= 5.260578707:
        hV = Region4.h4_s(s)
        if h < hV:
            return DiagramRegion.R4
        # Check if we are under the B23 validity area.
        if s <= 5.048096828:
            TMax = Region3.T3_ps(100, s)
            vmax = Region3.v3_ps(100, s)
            hMax = Region3.h3_rhoT(1 / vmax, TMax)
            if h < hMax:
                return DiagramRegion.R3
            logger.warning("Entalpy outside valid area")
            return DiagramRegion.NILL
        # In the area of B23 in s area.
        if h > 2812.942061:  # Above B23 in h_led
            if s > 5.09796573397125:
                TMax = Region2.T2_ps(100, s)
                hMax = Region2.h2_pT(100, TMax)
                if h < hMax:
                    return DiagramRegion.R2
                logger.warning("Entalpy outside valid area")
                return DiagramRegion.NILL
            logger.warning("Entropy outside valid area")
            return DiagramRegion.NILL
        if h < 2563.592004:  # Below B23 in h_led but we have already checked above for hV2c3b
            return DiagramRegion.R3
        # We are in the B23 field in both s and h joints.
        Tact = TB23_hs(h, s)
        pact = Region2.p2_hs(h, s)
        pBound = pB23_T(Tact)
        if pact > pBound:
            return DiagramRegion.R3
        return DiagramRegion.R2
    logger.warning("Entropy and Entalpy outside valid area")
    return DiagramRegion.NILL


def select_region_prho(p: float, rho: float) -> DiagramRegion:
    """
    Section 3.5 Regions as a function of p and rho
    Select diagram region based on values of preasure and density
    Returns `DiagramRegion.NILL` if arguments lead to no valid region

    :param p: preasure in [MPa]
    :param rho: density in [kg / m³]

    :return: diagram region
    """
    v = 1 / rho
    if (p < TRIPLE_POINT_PRESSURE) or (p > 100):
        logger.warning("Preasure outside valid area")
        return DiagramRegion.NILL
    if p < 16.5292:  # Below region 3, Check region 1,4,2
        if v < Region1.v1_pT(p, FREEZING_TEMPERATURE_H2O):
            # Observe that this is not actually min of v. Not valid for Water
            # if temperature is smaler than 4°C since it is ligther.
            logger.warning("Specific volume outside valid area")
            return DiagramRegion.NILL
        if v <= Region1.v1_pT(p, Region4.T4_p(p)):
            return DiagramRegion.R1
        if v < Region2.v2_pT(p, Region4.T4_p(p)):
            return DiagramRegion.R4
        if v <= Region2.v2_pT(p, 1073.15):
            return DiagramRegion.R2
        if p > 10:  # Above region 5
            logger.warning("Preasure outside valid area")
            return DiagramRegion.NILL
        if v <= Region5.v5_pT(p, 2073.15):
            return DiagramRegion.R5
    else:  # Check region 1,3,4,3,2 (Above the lowest point of region 3.)
        if v < Region1.v1_pT(p, FREEZING_TEMPERATURE_H2O):
            # Observe that this is not actually min of v. Not valid for Water
            # if temperature is smaler than 4°C since it is ligther.
            logger.warning("Specific volume outside valid area")
            return DiagramRegion.NILL
        if v < Region1.v1_pT(p, 623.15):
            return DiagramRegion.R1
        # Check if in region 3 or 4 (Below Reg 2)
        if v < Region2.v2_pT(p, TB23_p(p)):
            # Region 3 or 4
            if p > 22.064:  # Above region 4
                return DiagramRegion.R3
            if (v < Region3.v3_ph(p, Region4.h4L_p(p))) or (v > Region3.v3_ph(p, Region4.h4V_p(p))):  # Uses iteration!!
                return DiagramRegion.R3
            return DiagramRegion.R4
        # Check if region 2
        if v < Region2.v2_pT(p, 1073.15):
            return DiagramRegion.R2
    logger.warning("Preasure and Density outside valid area")
    return DiagramRegion.NILL
