# -*- coding: utf-8 -*-
"""
Section 3: Region Selection
"""
import math
import logging
from . import RegionBorders
from .Regions import Region1, Region2, Region3, Region4, Region5

logger = logging.getLogger(__name__)


def region_pT(p, T):
    """
    Section 3.1 Regions as a function of pT
    # function region_pT = region_pT(p, T)
    """
    if (T > 1073.15) and (p < 10) and (T < 2273.15) and (p > 0.000611):
        region_pT = 5
    elif (T <= 1073.15) and (T > 273.15) and (p <= 100) and (p > 0.000611):
        if T > 623.15:
            if p > RegionBorders.B23p_T(T):
                region_pT = 3
                if T < 647.096:
                    ps = Region4.p4_T(T)
                    if math.fabs(p - ps) < 0.00001:
                        region_pT = 4
            else:
                region_pT = 2
        else:
            ps = Region4.p4_T(T)
            if math.fabs(p - ps) < 0.00001:
                region_pT = 4
            elif p > ps:
                region_pT = 1
            else:
                region_pT = 2
    else:
        logger.warning('Temperature outside valid area')
        region_pT = 0  # %**Error, Outside valid area

    return region_pT


def region_ph(p, h):
    """Section 3.2 Regions as a function of ph
    # function region_ph = region_ph(p, h)
    """
    # %Check if outside pressure limits
    if (p < 0.000611657) or (p > 100):
        logger.warning('Preasure outside valid area')
        return 0
    # %Check if outside low h.
    if h < (0.963 * p + 2.2):  # %Linear adaption to h1_pt()+2 to speed up calcualations.
        if h < Region1.h1_pT(p, 273.15):
            logger.warning('Enthalpy outside valid area')
            return 0
    if p < 16.5292:  # % Bellow region 3, Check region 1, 4, 2, 5
        # % Check Region 1
        Ts = Region4.T4_p(p)
        hL = 109.6635 * math.log(p) + 40.3481 * p + 734.58  # % Approximate function for hL_p
        if math.fabs(h - hL) < 100:  # % if approximate is not god enough use real function
            hL = Region1.h1_pT(p, Ts)
        if h <= hL:
            return 1
        # % Check Region 4
        hV = 45.1768 * math.log(p) - 20.158 * p + 2804.4  # % Approximate function for hV_p
        if math.fabs(h - hV) < 50:  # % if approximate is not god enough use real function
            hV = Region2.h2_pT(p, Ts)
        if h < hV:
            return 4
        # % Check upper limit of region 2 Quick Test
        if h < 4000:
            return 2
        # % Check region 2 (Real value)
        h_45 = Region2.h2_pT(p, 1073.15)
        if h <= h_45:
            return 2
        # % Check region 5
        if p > 10:
            logger.warning('Preasure outside valid area')
            return 0
        h_5u = Region5.h5_pT(p, 2273.15)
        if h < h_5u:
            return 5
        logger.warning('Enthalpy outside valid area')
        return 0
    else:  # for p > 16.5292
        # % Check if in region1
        if h < Region1.h1_pT(p, 623.15):
            # region_ph = 1;
            return 1
        # % Check if in region 3 or 4 (Bellow Reg 2)
        if h < Region2.h2_pT(p, RegionBorders.B23T_p(p)):
            # % Region 3 or 4
            if p > RegionBorders.p3sat_h(h):
                return 3
            else:
                return 4
        # % Check if region 2
        if h < Region2.h2_pT(p, 1073.15):
            return 2
    logger.warning('Preasure outside valid area')
    return 0


def region_ps(p, s):
    """
    Section 3.3 Regions as a function of ps
    # function region_ps = region_ps(p, s)
    """
    if (p < 0.000611657) or (p > 100) or (s < 0) or (s > Region5.s5_pT(p, 2273.15)):
        logger.warning('Preasure or Entropy outside valid area')
        return 0
    # % Check region 5
    if s > Region2.s2_pT(p, 1073.15):
        if p <= 10:
            return 5
        else:
            logger.warning('Preasure outside valid area')
            return 0
    # % Check region 2
    if p > 16.529:
        ss = Region2.s2_pT(p, RegionBorders.B23T_p(p))  # % Between 5.047 & 5.261. Use to speed up !
    else:
        ss = Region2.s2_pT(p, Region4.T4_p(p))
    if s > ss:
        return 2
    # % Check region 3
    ss = Region1.s1_pT(p, 623.15)
    if (p > 16.529) and (s > ss):
        if p > RegionBorders.p3sat_s(s):
            return 3
        else:
            return 4
    # % Check region 4 (Not inside region 3)
    if (p < 16.529) and (s > Region1.s1_pT(p, Region4.T4_p(p))):
        return 4
    # % Check region 1
    if (p > 0.000611657) and (s > Region1.s1_pT(p, 273.15)):
        return 1
    return 1


def region_hs(h, s):
    """
    # Section 3.4 Regions as a function of hs
    # function region_hs = region_hs(h, s)
    """
    if s < -0.0001545495919:
        logger.warning('Entropy outside valid area')
        return 0
    # %Check linear adaption to p=0.000611. if bellow region 4.
    hMin = (((-0.0415878 - 2500.89262) / (-0.00015455 - 9.155759)) * s)
    if (s < 9.155759395) and (h < hMin):
        logger.warning('Entalpy or Entropy outside valid area')
        return 0
    # %******Kolla 1 eller 4. (+liten bit ???ver B13)
    if (s >= -0.0001545495919) and (s <= 3.77828134):
        if h < Region4.h4_s(s):
            return 4
        elif s < 3.397782955:  # %100MPa line is limiting
            TMax = Region1.T1_ps(100, s)
            hMax = Region1.h1_pT(100, TMax)
            if h < hMax:
                return 1
            else:
                logger.warning('Entalpy outside valid area')
                return 0
        else:  # %The point is either in region 4,1,3. Check B23
            hB = RegionBorders.hB13_s(s)
            if h < hB:
                return 1
            TMax = Region3.T3_ps(100, s)
            vmax = Region3.v3_ps(100, s)
            hMax = Region3.h3_rhoT(1 / vmax, TMax)
            if h < hMax:
                return 3
            else:
                logger.warning('Entalpy outside valid area')
                return 0

    # %******Kolla region 2 eller 4. (???vre delen av omr???de b23-> max)
    if (s >= 5.260578707) and (s <= 11.9212156897728):
        if s > 9.155759395:  # %Above region 4
            Tmin = Region2.T2_ps(0.000611, s)
            hMin = Region2.h2_pT(0.000611, Tmin)
            # %function adapted to h(1073.15,s)
            hMax = -0.07554022 * s ** 4 + 3.341571 * s ** 3 - 55.42151 * s ** 2 + 408.515 * s + 3031.338
            if (h > hMin) and (h < hMax):
                return 2
            else:
                logger.warning('Entalpy outside valid area')
                return 0
        hV = Region4.h4_s(s)
        if h < hV:  # %Region 4. Under region 3.
            return 4
        if s < 6.04048367171238:
            TMax = Region2.T2_ps(100, s)
            hMax = Region2.h2_pT(100, TMax)
        else:
            # %function adapted to h(1073.15,s)
            hMax = -2.988734 * s ** 4 + 121.4015 * s ** 3 - 1805.15 * s ** 2 + 11720.16 * s - 23998.33
        if h < hMax:  # %Region 2. ???ver region 4.
            return 2
        else:
            logger.warning('Entalpy outside valid area')
            return 0
    # %Kolla region 3 eller 4. Under kritiska punkten.
    if (s >= 3.77828134) and (s <= 4.41202148223476):
        hL = Region4.h4_s(s)
        if h < hL:
            return 4
        TMax = Region3.T3_ps(100, s)
        vmax = Region3.v3_ps(100, s)
        hMax = Region3.h3_rhoT(1 / vmax, TMax)
        if h < hMax:
            return 3
        else:
            logger.warning('Entalpy outside valid area')
            return 0
    # %Kolla region 3 eller 4 fr???n kritiska punkten till ???vre delen av b23
    if (s >= 4.41202148223476) and (s <= 5.260578707):
        hV = Region4.h4_s(s)
        if h < hV:
            return 4
        # %Kolla om vi ???r under b23 giltighetsomr???de.
        if s <= 5.048096828:
            TMax = Region3.T3_ps(100, s)
            vmax = Region3.v3_ps(100, s)
            hMax = Region3.h3_rhoT(1 / vmax, TMax)
            if h < hMax:
                return 3
            else:
                logger.warning('Entalpy outside valid area')
                return 0
        else:  # %Inom omr???det f???r B23 i s led.
            if h > 2812.942061:  # %Ovanf???r B23 i h_led
                if s > 5.09796573397125:
                    TMax = Region2.T2_ps(100, s)
                    hMax = Region2.h2_pT(100, TMax)
                    if h < hMax:
                        return 2
                    else:
                        logger.warning('Entalpy outside valid area')
                        return 0
                else:
                    logger.warning('Entropy outside valid area')
                    return 0
            if h < 2563.592004:  # %Nedanf???r B23 i h_led men vi har redan kollat ovanf???r hV2c3b
                return 3
            # %Vi ???r inom b23 omr???det i b???de s och h led.
            Tact = RegionBorders.TB23_hs(h, s)
            pact = Region2.p2_hs(h, s)
            pBound = RegionBorders.B23p_T(Tact)
            if pact > pBound:
                return 3
            else:
                return 2
    logger.warning('Entropy and Entalpy outside valid area')
    return 0


def region_prho(p, rho):
    """
    # Section 3.5 Regions as a function of p and rho
    # function region_prho = region_prho(p, rho)
    """
    v = 1 / rho
    if (p < 0.000611657) or (p > 100):
        logger.warning('Preasure outside valid area')
        return 0
    if p < 16.5292:  # %Bellow region 3, Check region 1,4,2
        if v < Region1.v1_pT(p, 273.15):  # %Observe that this is not actually min of v. Not valid Water of 4???C is ligther.
            logger.warning('Specific volume outside valid area')
            return 0
        if v <= Region1.v1_pT(p, Region4.T4_p(p)):
            return 1
        if v < Region2.v2_pT(p, Region4.T4_p(p)):
            return 4
        if v <= Region2.v2_pT(p, 1073.15):
            return 2
        if p > 10:  # %Above region 5
            logger.warning('Preasure outside valid area')
            return 0
        if v <= Region5.v5_pT(p, 2073.15):
            return 5
    else:  # %Check region 1,3,4,3,2 (Above the lowest point of region 3.)
        if v < Region1.v1_pT(p, 273.15):  # %Observe that this is not actually min of v. Not valid Water of 4???C is ligther.
            logger.warning('Specific volume outside valid area')
            return 0
        if v < Region1.v1_pT(p, 623.15):
            return 1
        # %Check if in region 3 or 4 (Bellow Reg 2)
        if v < Region2.v2_pT(p, RegionBorders.B23T_p(p)):
            # %Region 3 or 4
            if p > 22.064:  # %Above region 4
                return 3
            if (v < Region3.v3_ph(p, Region4.h4L_p(p))) or (v > Region3.v3_ph(p, Region4.h4V_p(p))):  # %Uses iteration!!
                return 3
            else:
                return 4
        # %Check if region 2
        if v < Region2.v2_pT(p, 1073.15):
            return 2
    logger.warning('Preasure and Density outside valid area')
    return 0
