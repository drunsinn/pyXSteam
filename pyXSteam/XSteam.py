#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import logging
from . import RegionSelection
from .Regions import Region1, Region2, Region3, Region4, Region5
from . import TransportProperties
from . import Constants
from .UnitConverter import UnitConverter
from . import IAPWS_R14


class XSteam(object):
    """Main pyXSteam object. Abstract of all other functions to allow auto selection of
    the correct region for each set of parameters.

    Args:
        unitSystem (int): set the unit system used for input and output values.
            Can be eather 0 (UNIT_SYSTEM_BARE), 1 (UNIT_SYSTEM_MKS) or 2 (UNIT_SYSTEM_FLS).
    """

    # Copy constant Values to expose them to the User
    UNIT_SYSTEM_BARE = UnitConverter.__UNIT_SYSTEM_BARE__
    UNIT_SYSTEM_MKS = UnitConverter.__UNIT_SYSTEM_MKS__
    UNIT_SYSTEM_FLS = UnitConverter.__UNIT_SYSTEM_FLS__

    TYPE_ICE_Ih = IAPWS_R14.__TYPE_ICE_Ih__
    TYPE_ICE_III = IAPWS_R14.__TYPE_ICE_III__
    TYPE_ICE_V = IAPWS_R14.__TYPE_ICE_V__
    TYPE_ICE_VI = IAPWS_R14.__TYPE_ICE_VI__
    TYPE_ICE_VII = IAPWS_R14.__TYPE_ICE_VII__

    def __init__(self, unitSystem=UnitConverter.__UNIT_SYSTEM_BARE__):
        self.logger = logging.getLogger(__name__)
        self.unitConverter = UnitConverter(unitSystem)
        self.logger.info('initialised pyXSteam with Unit System "{}"'.format(self.unitConverter))

    def specificGasConstant(self):
        """returns the specific Gas Constant R in kJ kg^-1 K^-1"""
        return Constants.__SPECIFIC_GAS_CONSTANT__

    def criticalTemperatur(self):
        """returns the specific temperature with conversion to the selected unit system"""
        return self.unitConverter.fromSIunit_T(Constants.__CRITICAL_TEMPERATURE__)

    def criticalPressure(self):
        """returns the specific pressure with conversion to the selected unit system"""
        return self.unitConverter.fromSIunit_p(Constants.__CRITICAL_PRESSURE__)

    def criticalDensity(self):
        """returns the specific density with conversion to the selected unit system"""
        return self.unitConverter.fromSIunit_p(Constants.__CRITICAL_DENSITY__)

    def triplePointTemperatur(self):
        """returns the temperature of the triple point with conversion to the selected unit system"""
        return self.unitConverter.fromSIunit_T(Constants.__TRIPLE_POINT_TEMPERATURE__)

    def triplePointPressure(self):
        """returns the Pressure of the triple poin with conversion to the selected unit systemt"""
        return self.unitConverter.fromSIunit_p(Constants.__TRIPLE_POINT_PRESSURE__)

    def zeroPointTemperature(self):
        """returns the absolute zero temperature with conversion to the selected unit system"""
        return self.unitConverter.fromSIunit_T(0.0)

# ***********************************************************************************************************
# Section 1.2 temperature
    def tsat_p(self, p):
        """Saturation-temperature as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            tsat (float): saturation temperature or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            return self.unitConverter.fromSIunit_T(Region4.T4_p(p))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

    def tsat_s(self, s):
        """Saturation-temperature as a function of entropy

        Args:
            s (float): entropy value

        Returns:
            tsat (float): saturation temperature or NaN if arguments are out of range
        """
        s = self.unitConverter.toSIunit_s(s)
        if (s > -0.0001545495919) and (s < 9.155759395):
            ps = Region4.p4_s(s)
            return self.unitConverter.fromSIunit_T(Region4.T4_p(ps))
        else:
            self.logger.warning("Entropy value of {:f} is out of range".format(s))
            return float("NaN")

    def t_ph(self, p, h):
        """temperature as a function of pressure and enthalpy

        Args:
            p (float): preasure value
            h (float): enthalpy value

        Returns:
            t (float): temperature or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        region = RegionSelection.region_ph(p, h)
        if region == 1:
            return self.unitConverter.fromSIunit_T(Region1.T1_ph(p, h))
        elif region == 2:
            return self.unitConverter.fromSIunit_T(Region2.T2_ph(p, h))
        elif region == 3:
            return self.unitConverter.fromSIunit_T(Region3.T3_ph(p, h))
        elif region == 4:
            return self.unitConverter.fromSIunit_T(Region4.T4_p(p))
        elif region == 5:
            return self.unitConverter.fromSIunit_T(Region5.T5_ph(p, h))
        else:
            self.logger.warning("Region switch t_ph returned unknown value: {:d}".format(region))
            return float("NaN")

    def t_ps(self, p, s):
        """temperature as a function of pressure and entropy

        Args:
            p (float): preasure value
            s (float): entropy value
            h (float): enthalpy value

        Returns:
            t (float): temperature or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        region = RegionSelection.region_ps(p, s)
        if region == 1:
            return self.unitConverter.fromSIunit_T(Region1.T1_ps(p, s))
        elif region == 2:
            return self.unitConverter.fromSIunit_T(Region2.T2_ps(p, s))
        elif region == 3:
            return self.unitConverter.fromSIunit_T(Region3.T3_ps(p, s))
        elif region == 4:
            return self.unitConverter.fromSIunit_T(Region4.T4_p(p))
        elif region == 5:
            return self.unitConverter.fromSIunit_T(Region5.T5_ps(p, s))
        else:
            self.logger.warning("Region switch t_ps returned unknown value: {:d}".format(region))
            return float("NaN")

    def t_hs(self, h, s):
        """temperature as a function of enthalpy and entropy

        Args:
            h (float): enthalpy value
            s (float): entropy value

        Returns:
            t (float): temperature or NaN if arguments are out of range
        """
        h = self.unitConverter.toSIunit_h(h)
        s = self.unitConverter.toSIunit_s(s)
        region = RegionSelection.region_hs(h, s)
        if region == 1:
            p1 = Region1.p1_hs(h, s)
            return self.unitConverter.fromSIunit_T(Region1.T1_ph(p1, h))
        elif region == 2:
            p2 = Region2.p2_hs(h, s)
            return self.unitConverter.fromSIunit_T(Region2.T2_ph(p2, h))
        elif region == 3:
            p3 = Region3.p3_hs(h, s)
            return self.unitConverter.fromSIunit_T(Region3.T3_ph(p3, h))
        elif region == 4:
            return self.unitConverter.fromSIunit_T(Region4.T4_hs(h, s))
        elif region == 5:
            self.logger.error("Functions of hs is not available in region 5")
            return float("NaN")
        else:
            self.logger.warning("Region switch t_hs returned unknown value: {:d}".format(region))
            return float("NaN")

# ***********************************************************************************************************
# Section 1.3 Pressure (p)
    def psat_s(self, s):
        """Saturation-Pressure as a function of entropy

        Args:
            s (float): entropy value

        Returns:
            psat (float): saturation pressure or NaN if arguments are out of range
        """
        s = self.unitConverter.toSIunit_s(s)
        if (s > -0.0001545495919) and (s < 9.155759395):
            return self.unitConverter.fromSIunit_p(Region4.p4_s(s))
        else:
            self.logger.warning("Entropy value of {:f} is out of range".format(s))
            return float("NaN")

    def psat_t(self, t):
        """Saturation-Pressure as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            psat (float): saturation pressure or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T < 647.096) and (T > 273.1):
            return self.unitConverter.fromSIunit_p(Region4.p4_T(T))
        else:
            self.logger.warning("Temperature value {:f} out of range".format(T))
            return float("NaN")

    def p_hs(self, h, s):
        """Pressure  as a function of enthalpy and entropy

        Args:
            h (float): enthalpy value
            s (float): entropy value

        Returns:
            p (float): pressure or NaN if arguments are out of range
        """
        h = self.unitConverter.toSIunit_h(h)
        s = self.unitConverter.toSIunit_s(s)
        region = RegionSelection.region_hs(h, s)
        if region == 1:
            return self.unitConverter.fromSIunit_p(Region1.p1_hs(h, s))
        elif region == 2:
            return self.unitConverter.fromSIunit_p(Region2.p2_hs(h, s))
        elif region == 3:
            return self.unitConverter.fromSIunit_p(Region3.p3_hs(h, s))
        elif region == 4:
            tSat = Region4.T4_hs(h, s)
            return self.unitConverter.fromSIunit_p(Region4.p4_T(tSat))
        elif region == 5:
            self.logger.warning('functions of hs is not available in region 5')
            return float("NaN")
        else:
            self.logger.warning("Region switch p_hs returned unknown value: {:d}".format(region))
            return float("NaN")

    def p_hrho(self, h, rho):
        """Pressure as a function of h and rho.

        Very unaccurate for solid water region since it's almost incompressible!

        Not valid for water or sumpercritical since water rho does not change very much with p. Uses iteration to find p.

        Args:
            h (float): enthalpy value
            rho (float): density value

        Returns:
            p (float): pressure or NaN if arguments are out of range
        """
        if rho <= 0.0:
            raise Exception('rho out of range')
            return float("NaN")
        h = self.unitConverter.toSIunit_h(h)
        High_Bound = self.unitConverter.fromSIunit_p(100)
        Low_Bound = self.unitConverter.fromSIunit_p(0.000611657)
        ps = self.unitConverter.fromSIunit_p(10)
        rhos = 1 / self.v_ph(ps, h)
        while math.fabs(rho - rhos) > 0.0000001:
            # rhos = 1 / XSteam('v_ph', ps, h)
            rhos = 1 / self.v_ph(ps, h)
            if rhos >= rho:
                High_Bound = ps
            else:
                Low_Bound = ps
            ps = (Low_Bound + High_Bound) / 2
        return ps

# ***********************************************************************************************************
# Section 1.4 Enthalpy (h)
    def hV_p(self, p):
        """Saturated vapour enthalpy as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            hV (float): saturated vapour enthalpy or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            return self.unitConverter.fromSIunit_h(Region4.h4V_p(p))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

    def hL_p(self, p):
        """Saturated liquid enthalpy as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            hL (float): saturated liquid enthalpy or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            return self.unitConverter.fromSIunit_h(Region4.h4L_p(p))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

    def hV_t(self, t):
        """Saturated vapour enthalpy as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            hV (float): saturated vapour enthalpy or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            p = Region4.p4_T(T)
            return self.unitConverter.fromSIunit_h(Region4.h4V_p(p))
        else:
            self.logger.warning('Temperature out of range')
            return float("NaN")

    def hL_t(self, t):
        """Saturated liquid enthalpy as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            hL (float): saturated liquid enthalpy or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            p = Region4.p4_T(T)
            return self.unitConverter.fromSIunit_h(Region4.h4L_p(p))
        else:
            self.logger.warning('Temperature out of range')
            return float("NaN")

    def h_pt(self, p, t):
        """Entalpy as a function of pressure and temperature

        Args:
            p (float): preasure value
            t (float): temperature value

        Returns:
            h (float): enthalpy or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        region = RegionSelection.region_pT(p, T)
        if region == 1:
            return self.unitConverter.fromSIunit_h(Region1.h1_pT(p, T))
        elif region == 2:
            return self.unitConverter.fromSIunit_h(Region2.h2_pT(p, T))
        elif region == 3:
            return self.unitConverter.fromSIunit_h(Region3.h3_pT(p, T))
        elif region == 4:
            self.logger.warning('function h_pt is not available in region 4')
            return float("NaN")
        elif region == 5:
            return self.unitConverter.fromSIunit_h(Region5.h5_pT(p, T))
        else:
            self.logger.warning("Region switch h_pt returned unknown value: {:d}".format(region))
            return float("NaN")

    def h_ps(self, p, s):
        """Entalpy as a function of pressure and entropy

        Args:
            p (float): preasure value
            s (float): entropy value

        Returns:
            h (float): enthalpy or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        region = RegionSelection.region_ps(p, s)
        if region == 1:
            return self.unitConverter.fromSIunit_h(Region1.h1_pT(p, Region1.T1_ps(p, s)))
        elif region == 2:
            return self.unitConverter.fromSIunit_h(Region2.h2_pT(p, Region2.T2_ps(p, s)))
        elif region == 3:
            return self.unitConverter.fromSIunit_h(Region3.h3_rhoT(1 / Region3.v3_ps(p, s), Region3.T3_ps(p, s)))
        elif region == 4:
            xs = Region4.x4_ps(p, s)
            return self.unitConverter.fromSIunit_h(xs * Region4.h4V_p(p) + (1 - xs) * Region4.h4L_p(p))
        elif region == 5:
            return self.unitConverter.fromSIunit_h(Region5.h5_pT(p, Region5.T5_ps(p, s)))
        else:
            self.logger.warning("Region switch h_ps returned unknown value: {:d}".format(region))
            return float("NaN")

    def h_px(self, p, x):
        """Entalpy as a function of pressure and vapour fraction

        Args:
            p (float): preasure value
            x (float): vapour fraction value

        Returns:
            h (float): enthalpy or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        x = self.unitConverter.toSIunit_x(x)
        if (x > 1) or (x < 0) or (p >= 22.064):
            self.logger.warning('Vapor fraction and/or pressure out of range')
            return float("NaN")
        hL = Region4.h4L_p(p)
        hV = Region4.h4V_p(p)
        return hL + x * (hV - hL)

    def h_prho(self, p, rho):
        """Entalpy as a function of pressure and density. Observe for low temperatures (liquid) this equation has 2 solutions

        Args:
            p (float): preasure value
            rho (float): density value

        Returns:
            h (float): enthalpy or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        rho = 1 / self.unitConverter.toSIunit_v(1 / float(rho))
        region = RegionSelection.region_prho(p, rho)
        if region == 1:
            return self.unitConverter.fromSIunit_h(Region1.h1_pT(p, Region1.T1_prho(p, rho)))
        elif region == 2:
            return self.unitConverter.fromSIunit_h(Region2.h2_pT(p, Region2.T2_prho(p, rho)))
        elif region == 3:
            return self.unitConverter.fromSIunit_h(Region3.h3_rhoT(rho, Region3.T3_prho(p, rho)))
        elif region == 4:
            if p < 16.529:
                vV = Region2.v2_pT(p, Region4.T4_p(p))
                vL = Region1.v1_pT(p, Region4.T4_p(p))
            else:
                vV = Region3.v3_ph(p, Region4.h4V_p(p))
                vL = Region3.v3_ph(p, Region4.h4L_p(p))
            hV = Region4.h4V_p(p)
            hL = Region4.h4L_p(p)
            x = (1 / rho - vL) / (vV - vL)
            return self.unitConverter.fromSIunit_h((1 - x) * hL + x * hV)
        elif region == 5:
            return self.unitConverter.fromSIunit_h(Region5.h5_pT(p, Region5.T5_prho(p, rho)))
        else:
            self.logger.warning("Region switch h_prho returned unknown value: {:d}".format(region))
            return float("NaN")

    def h_tx(self, t, x):
        """Entalpy as a function of temperature and vapour fraction

        Args:
            t (float): temperature value
            x (float): vapour fraction value

        Returns:
            h (float): enthalpy or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        x = self.unitConverter.toSIunit_x(x)
        if (x > 1) or (x < 0) or (T >= 647.096):
            self.logger.warning('Vapor fraction and/or pressure out of range')
            return float("NaN")
        p = Region4.p4_T(T)
        hL = Region4.h4L_p(p)
        hV = Region4.h4V_p(p)
        return hL + x * (hV - hL)

# ***********************************************************************************************************
# Section 1.5 Specific Volume (v)
    def vV_p(self, p):
        """Saturated vapour volume as a function of pressure

        Args:
            p (float): preasure value
            x (float): vapour fraction value

        Returns:
            vV (float): saturated vapour volume or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_v(Region2.v2_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_v(Region3.v3_ph(p, Region4.h4V_p(p)))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

    def vL_p(self, p):
        """Saturated liquid volume as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            vV (float): saturated liquid volume or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_v(Region1.v1_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_v(Region3.v3_ph(p, Region4.h4L_p(p)))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

    def vV_t(self, t):
        """Saturated vapour volume as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            vV (float): saturated vapour volume or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if T <= 623.15:
                return self.unitConverter.fromSIunit_v(Region2.v2_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_v(Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T))))
        else:
            self.logger.warning('temperature out of range')
            return float("NaN")

    def vL_t(self, t):
        """Saturated liquid volume as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            vV (float): saturated liquid volume or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if T <= 623.15:
                return self.unitConverter.fromSIunit_v(Region1.v1_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_v(Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T))))
        else:
            self.logger.warning('temperature out of range')
            return float("NaN")

    def v_pt(self, p, t):
        """Specific volume as a function of pressure and temperature

        Args:
            p (float): preasure value
            t (float): temperature value

        Returns:
            v (float): specific volume or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        region = RegionSelection.region_pT(p, T)
        if region == 1:
            return self.unitConverter.fromSIunit_v(Region1.v1_pT(p, T))
        elif region == 2:
            return self.unitConverter.fromSIunit_v(Region2.v2_pT(p, T))
        elif region == 3:
            return self.unitConverter.fromSIunit_v(Region3.v3_ph(p, Region3.h3_pT(p, T)))
        elif region == 4:
            self.logger.warning('function v_pt is not available in region 4')
            return float("NaN")
        elif region == 5:
            return self.unitConverter.fromSIunit_v(Region5.v5_pT(p, T))
        else:
            self.logger.warning("Region switch v_pt returned unknown value: {:d}".format(region))
            return float("NaN")

    def v_ph(self, p, h):
        """Specific volume as a function of pressure and enthalpy

        Args:
            p (float): preasure value
            h (float): enthalpy value

        Returns:
            v (float): specific volume or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        region = RegionSelection.region_ph(p, h)
        if region == 1:
            return self.unitConverter.fromSIunit_v(Region1.v1_pT(p, Region1.T1_ph(p, h)))
        elif region == 2:
            return self.unitConverter.fromSIunit_v(Region2.v2_pT(p, Region2.T2_ph(p, h)))
        elif region == 3:
            return self.unitConverter.fromSIunit_v(Region3.v3_ph(p, h))
        elif region == 4:
            xs = Region4.x4_ph(p, h)
            if p < 16.529:
                v4v = Region2.v2_pT(p, Region4.T4_p(p))
                v4L = Region1.v1_pT(p, Region4.T4_p(p))
            else:
                v4v = Region3.v3_ph(p, Region4.h4V_p(p))
                v4L = Region3.v3_ph(p, Region4.h4L_p(p))
            return self.unitConverter.fromSIunit_v((xs * v4v + (1 - xs) * v4L))
        elif region == 5:
            Ts = Region5.T5_ph(p, h)
            return self.unitConverter.fromSIunit_v(Region5.v5_pT(p, Ts))
        else:
            self.logger.warning("Region switch v_ph returned unknown value: {:d}".format(region))
            return float("NaN")

    def v_ps(self, p, s):
        """Specific volume as a function of pressure and entropy

        Args:
            p (float): preasure value
            s (float): entropy value

        Returns:
            v (float): specific volume or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        region = RegionSelection.region_ps(p, s)
        if region == 1:
            return self.unitConverter.fromSIunit_v(Region1.v1_pT(p, Region1.T1_ps(p, s)))
        elif region == 2:
            return self.unitConverter.fromSIunit_v(Region2.v2_pT(p, Region2.T2_ps(p, s)))
        elif region == 3:
            return self.unitConverter.fromSIunit_v(Region3.v3_ps(p, s))
        elif region == 4:
            xs = Region4.x4_ps(p, s)
            if p < 16.529:
                v4v = Region2.v2_pT(p, Region4.T4_p(p))
                v4L = Region1.v1_pT(p, Region4.T4_p(p))
            else:
                v4v = Region3.v3_ph(p, Region4.h4V_p(p))
                v4L = Region3.v3_ph(p, Region4.h4L_p(p))
            return self.unitConverter.fromSIunit_v((xs * v4v + (1 - xs) * v4L))
        elif region == 5:
            Ts = Region5.T5_ps(p, s)
            return self.unitConverter.fromSIunit_v(Region5.v5_pT(p, Ts))
        else:
            self.logger.warning("Region switch v_ps returned unknown value: {:d}".format(region))
            return float("NaN")

# ***********************************************************************************************************
# Section 1.6 Density (rho)
# Density is calculated as 1/v. See section 1.5 Volume
    def rhoV_p(self, p):
        """Saturated vapour density as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            rhoV (float): saturated vapour density
        """
        return 1 / self.vV_p(p)

    def rhoL_p(self, p):
        """Saturated liquid density as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            rhoV (float): saturated liquid density
        """
        return 1 / self.vL_p(p)

    def rhoV_t(self, t):
        """Saturated vapour density as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            rhoV (float): saturated vapour density
        """
        return 1 / self.vV_t(t)

    def rhoL_t(self, t):
        """Saturated liquid density as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            rhoV (float): saturated liquid density
        """
        return 1 / self.vL_t(t)

    def rho_pt(self, p, t):
        """Density as a function of pressure and temperature

        Args:
            p (float): preasure value
            t (float): temperature value

        Returns:
            rho (float): density
        """
        return 1 / self.v_pt(p, t)

    def rho_ph(self, p, h):
        """Density as a function of pressure and enthalpy

        Args:
            p (float): preasure value
            h (float): enthalpy value

        Returns:
            rho (float): density
        """
        return 1 / self.v_ph(p, h)

    def rho_ps(self, p, s):
        """Density as a function of pressure and entropy

        Args:
            p (float): preasure value
            s (float): entropy value

        Returns:
            rho (float): density
        """
        return 1 / self.v_ps(p, s)

# ***********************************************************************************************************
# Section 1.7 Specific entropy (s)
    def sV_p(self, p):
        """Saturated vapour entropy as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            sV (float): Saturated vapour entropy or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_s(Region2.s2_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_s(Region3.s3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p)))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

    def sL_p(self, p):
        """Saturated liquid entropy as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            sL (float): Saturated liquid entropy or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_s(Region1.s1_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_s(Region3.s3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p)))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

    def sV_t(self, t):
        """Saturated vapour entropy as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            sV (float): Saturated vapour entropy or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if T <= 623.15:
                return self.unitConverter.fromSIunit_s(Region2.s2_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_s(Region3.s3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T)))), T))
        else:
            self.logger.warning('temperature out of range')
            return float("NaN")

    def sL_t(self, t):
        """Saturated liquid entropy as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            sL (float): Saturated liquid entropy or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if T <= 623.15:
                return self.unitConverter.fromSIunit_s(Region1.s1_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_s(Region3.s3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T)))), T))
        else:
            self.logger.warning('temperature out of range')
            return float("NaN")

    def s_pt(self, p, t):
        """Specific entropy as a function of pressure and temperature (Returns saturated vapour entalpy if mixture)

        Args:
            p (float): preasure value
            t (float): temperature value

        Returns:
            s (float): entropy or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        region = RegionSelection.region_pT(p, T)
        if region == 1:
            return self.unitConverter.fromSIunit_s(Region1.s1_pT(p, T))
        elif region == 2:
            return self.unitConverter.fromSIunit_s(Region2.s2_pT(p, T))
        elif region == 3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self.unitConverter.fromSIunit_s(Region3.s3_rhoT(rhos, T))
        elif region == 4:
            self.logger.warning('function s_pt is not available in region 4')
            return float("NaN")
        elif region == 5:
            return self.unitConverter.fromSIunit_s(Region5.s5_pT(p, T))
        else:
            self.logger.warning("Region switch s_pt returned unknown value: {:d}".format(region))
            return float("NaN")

    def s_ph(self, p, h):
        """Specific entropy as a function of pressure and enthalpy

        Args:
            p (float): preasure value
            h (float): enthalpy value

        Returns:
            s (float): entropy or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        region = RegionSelection.region_ph(p, h)
        if region == 1:
            T = Region1.T1_ph(p, h)
            return self.unitConverter.fromSIunit_s(Region1.s1_pT(p, T))
        elif region == 2:
            T = Region2.T2_ph(p, h)
            return self.unitConverter.fromSIunit_s(Region2.s2_pT(p, T))
        elif region == 3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self.unitConverter.fromSIunit_s(Region3.s3_rhoT(rhos, Ts))
        elif region == 4:
            Ts = Region4.T4_p(p)
            xs = Region4.x4_ph(p, h)
            if p < 16.529:
                s4v = Region2.s2_pT(p, Ts)
                s4L = Region1.s1_pT(p, Ts)
            else:
                v4v = Region3.v3_ph(p, Region4.h4V_p(p))
                s4v = Region3.s3_rhoT(1 / v4v, Ts)
                v4L = Region3.v3_ph(p, Region4.h4L_p(p))
                s4L = Region3.s3_rhoT(1 / v4L, Ts)
            return self.unitConverter.fromSIunit_s((xs * s4v + (1 - xs) * s4L))
        elif region == 5:
            T = Region5.T5_ph(p, h)
            return self.unitConverter.fromSIunit_s(Region5.s5_pT(p, T))
        else:
            self.logger.warning("Region switch s_ph returned unknown value: {:d}".format(region))
            return float("NaN")

# ***********************************************************************************************************
# Section 1.8 Specific internal energy (u)
    def uV_p(self, p):
        """Saturated vapour internal energy as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            uV (float): saturated vapour internal energy or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_u(Region2.u2_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_u(Region3.u3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p)))
        else:
            return float("NaN")

    def uL_p(self, p):
        """Saturated liquid internal energy as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            uL (float): saturated liquid internal energy
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_u(Region1.u1_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_u(Region3.u3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p)))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

    def uV_t(self, t):
        """Saturated vapour internal energy as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            uV (float): saturated vapour internal energy or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if T <= 623.15:
                return self.unitConverter.fromSIunit_u(Region2.u2_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_u(Region3.u3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T)))), T))
        else:
            self.logger.warning('temperature out of range')
            return float("NaN")

    def uL_t(self, t):
        """Saturated liquid internal energy as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            uL (float): saturated liquid internal energy or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if T <= 623.15:
                return self.unitConverter.fromSIunit_u(Region1.u1_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_u(Region3.u3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T)))), T))
        else:
            self.logger.warning('temperature out of range')
            return float("NaN")

    def u_pt(self, p, t):
        """Specific internal energy as a function of pressure and temperature

        Args:
            p (float): preasure value
            t (float): temperature value

        Returns:
            u (float): specific internal energy or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        region = RegionSelection.region_pT(p, T)
        if region == 1:
            return self.unitConverter.fromSIunit_u(Region1.u1_pT(p, T))
        elif region == 2:
            return self.unitConverter.fromSIunit_u(Region2.u2_pT(p, T))
        elif region == 3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self.unitConverter.fromSIunit_u(Region3.u3_rhoT(rhos, T))
        elif region == 4:
            self.logger.warning('function u_pt is not available in region 4')
            return float("NaN")
        elif region == 5:
            return self.unitConverter.fromSIunit_u(Region5.u5_pT(p, T))
        else:
            self.logger.warning("Region switch u_pt returned unknown value: {:d}".format(region))
            return float("NaN")

    def u_ph(self, p, h):
        """Specific internal energy as a function of pressure and enthalpy

        Args:
            p (float): preasure value
            h (float): enthalpy value

        Returns:
            u (float): specific internal energy or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        region = RegionSelection.region_ph(p, h)
        if region == 1:
            Ts = Region1.T1_ph(p, h)
            return self.unitConverter.fromSIunit_u(Region1.u1_pT(p, Ts))
        elif region == 2:
            Ts = Region2.T2_ph(p, h)
            return self.unitConverter.fromSIunit_u(Region2.u2_pT(p, Ts))
        elif region == 3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self.unitConverter.fromSIunit_u(Region3.u3_rhoT(rhos, Ts))
        elif region == 4:
            Ts = Region4.T4_p(p)
            xs = Region4.x4_ph(p, h)
            if p < 16.529:
                u4v = Region2.u2_pT(p, Ts)
                u4L = Region1.u1_pT(p, Ts)
            else:
                v4v = Region3.v3_ph(p, Region4.h4V_p(p))
                u4v = Region3.u3_rhoT(1 / v4v, Ts)
                v4L = Region3.v3_ph(p, Region4.h4L_p(p))
                u4L = Region3.u3_rhoT(1 / v4L, Ts)
            return self.unitConverter.fromSIunit_u((xs * u4v + (1 - xs) * u4L))
        elif region == 5:
            Ts = Region5.T5_ph(p, h)
            return self.unitConverter.fromSIunit_u(Region5.u5_pT(p, Ts))
        else:
            self.logger.warning("Region switch u_ph returned unknown value: {:d}".format(region))
            return float("NaN")

    def u_ps(self, p, s):
        """Specific internal energy as a function of pressure and entropy

        Args:
            p (float): preasure value
            s (float): entropy value

        Returns:
            u (float): specific internal energy or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        region = RegionSelection.region_ps(p, s)
        if region == 1:
            Ts = Region1.T1_ps(p, s)
            return self.unitConverter.fromSIunit_u(Region1.u1_pT(p, Ts))
        elif region == 2:
            Ts = Region2.T2_ps(p, s)
            return self.unitConverter.fromSIunit_u(Region2.u2_pT(p, Ts))
        elif region == 3:
            rhos = 1 / Region3.v3_ps(p, s)
            Ts = Region3.T3_ps(p, s)
            return self.unitConverter.fromSIunit_u(Region3.u3_rhoT(rhos, Ts))
        elif region == 4:
            if p < 16.529:
                uLp = Region1.u1_pT(p, Region4.T4_p(p))
                uVp = Region2.u2_pT(p, Region4.T4_p(p))
            else:
                uLp = Region3.u3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p))
                uVp = Region3.u3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p))
            xs = Region4.x4_ps(p, s)
            return self.unitConverter.fromSIunit_u((xs * uVp + (1 - xs) * uLp))
        elif region == 5:
            Ts = Region5.T5_ps(p, s)
            return self.unitConverter.fromSIunit_u(Region5.u5_pT(p, Ts))
        else:
            self.logger.warning("Region switch u_ps returned unknown value: {:d}".format(region))
            return float("NaN")

# ***********************************************************************************************************
# Section 1.9 Specific isobaric heat capacity (Cp)
    def CpV_p(self, p):
        """Saturated vapour heat capacity as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            CpV (float): saturated vapour heat capacity or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_Cp(Region2.Cp2_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_Cp(Region3.Cp3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p)))
        else:
            self.logger.warning('preassure out of range')
            return float("NaN")

    def CpL_p(self, p):
        """Saturated liquid heat capacity as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            CpL (float): saturated liquid heat capacity or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_Cp(Region1.Cp1_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_Cp(Region3.Cp3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p)))
        else:
            self.logger.warning('preassure out of range')
            return float("NaN")

    def CpV_t(self, t):
        """Saturated vapour heat capacity as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            CpV (float): saturated vapour heat capacity or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if (T <= 623.15):
                return self.unitConverter.fromSIunit_Cp(Region2.Cp2_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_Cp(Region3.Cp3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T)))), T))
        else:
            self.logger.warning('temperature out of range')
            return float("NaN")

    def CpL_t(self, t):
        """Saturated liquid heat capacity as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            CpL (float): saturated liquid heat capacity or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if (T <= 623.15):
                return self.unitConverter.fromSIunit_Cp(Region1.Cp1_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_Cp(Region3.Cp3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T)))), T))
        else:
            self.logger.warning('temperature out of range')
            return float("NaN")

    def Cp_pt(self, p, t):
        """Specific isobaric heat capacity as a function of pressure and temperature

        Args:
            p (float): preasure value
            t (float): temperature value

        Returns:
            Cp (float): specific isobaric heat capacity or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        region = RegionSelection.region_pT(p, T)
        if region == 1:
            return self.unitConverter.fromSIunit_Cp(Region1.Cp1_pT(p, T))
        elif region == 2:
            return self.unitConverter.fromSIunit_Cp(Region2.Cp2_pT(p, T))
        elif region == 3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self.unitConverter.fromSIunit_Cp(Region3.Cp3_rhoT(rhos, T))
        elif region == 4:
            self.logger.warning('function Cp_pt is not available in region 4')
            return float("NaN")
        elif region == 5:
            return self.unitConverter.fromSIunit_Cp(Region5.Cp5_pT(p, T))
        else:
            self.logger.warning("Region switch Cp_pt returned unknown value: {:d}".format(region))
            return float("NaN")

    def Cp_ph(self, p, h):
        """Specific isobaric heat capacity as a function of pressure and enthalpy

        Args:
            p (float): preasure value
            h (float): enthalpy value

        Returns:
            Cp (float): specific isobaric heat capacity or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        region = RegionSelection.region_ph(p, h)
        if region == 1:
            Ts = Region1.T1_ph(p, h)
            return self.unitConverter.fromSIunit_Cp(Region1.Cp1_pT(p, Ts))
        elif region == 2:
            Ts = Region2.T2_ph(p, h)
            return self.unitConverter.fromSIunit_Cp(Region2.Cp2_pT(p, Ts))
        elif region == 3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self.unitConverter.fromSIunit_Cp(Region3.Cp3_rhoT(rhos, Ts))
        elif region == 4:
            self.logger.warning('function Cp_ph is not available in region 4')
            return float("NaN")
        elif region == 5:
            Ts = Region5.T5_ph(p, h)
            return self.unitConverter.fromSIunit_Cp(Region5.Cp5_pT(p, Ts))
        else:
            self.logger.warning("Region switch Cp_ph returned unknown value: {:d}".format(region))
            return float("NaN")

    def Cp_ps(self, p, s):
        """Specific isobaric heat capacity as a function of pressure and entropy

        Args:
            p (float): preasure value
            s (float): entropy value

        Returns:
            Cp (float): specific isobaric heat capacity or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        region = RegionSelection.region_ps(p, s)
        if region == 1:
            Ts = Region1.T1_ps(p, s)
            return self.unitConverter.fromSIunit_Cp(Region1.Cp1_pT(p, Ts))
        elif region == 2:
            Ts = Region2.T2_ps(p, s)
            return self.unitConverter.fromSIunit_Cp(Region2.Cp2_pT(p, Ts))
        elif region == 3:
            rhos = 1 / Region3.v3_ps(p, s)
            Ts = Region3.T3_ps(p, s)
            return self.unitConverter.fromSIunit_Cp(Region3.Cp3_rhoT(rhos, Ts))
        elif region == 4:
            self.logger.warning('function Cp_ps is not available in region 4')
            return float("NaN")
        elif region == 5:
            Ts = Region5.T5_ps(p, s)
            return self.unitConverter.fromSIunit_Cp(Region5.Cp5_pT(p, Ts))
        else:
            self.logger.warning("Region switch Cp_ps returned unknown value: {:d}".format(region))
            return float("NaN")

# ***********************************************************************************************************
# Section 1.10 Specific isochoric heat capacity (Cv)
    def CvV_p(self, p):
        """Saturated vapour isochoric heat capacity as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            CvV (float): saturated vapour isochoric heat capacity or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_Cv(Region2.Cv2_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_Cv(Region3.Cv3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p)))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

    def CvL_p(self, p):
        """Saturated liquid isochoric heat capacity as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            CvL (float): saturated liquid isochoric heat capacity or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_Cv(Region1.Cv1_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_Cv(Region3.Cv3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p)))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

    def CvV_t(self, t):
        """Saturated vapour isochoric heat capacity as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            CvV (float): saturated vapour isochoric heat capacity or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if (T <= 623.15):
                return self.unitConverter.fromSIunit_Cv(Region2.Cv2_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_Cv(Region3.Cv3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T)))), T))
        else:
            self.logger.warning('temperature out of range')
            return float("NaN")

    def CvL_t(self, t):
        """Saturated liquid isochoric heat capacity as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            CvL (float): saturated liquid isochoric heat capacity or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if (T <= 623.15):
                return self.unitConverter.fromSIunit_Cv(Region1.Cv1_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_Cv(Region3.Cv3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T)))), T))
        else:
            self.logger.warning('temperature out of range')
            return float("NaN")

    def Cv_pt(self, p, t):
        """Specific isochoric heat capacity as a function of pressure and temperature

        Args:
            p (float): preasure value
            t (float): temperature value

        Returns:
            Cv (float): specific isochoric heat capacity or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        region = RegionSelection.region_pT(p, T)
        if region == 1:
            return self.unitConverter.fromSIunit_Cv(Region1.Cv1_pT(p, T))
        elif region == 2:
            return self.unitConverter.fromSIunit_Cv(Region2.Cv2_pT(p, T))
        elif region == 3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self.unitConverter.fromSIunit_Cv(Region3.Cv3_rhoT(rhos, T))
        elif region == 4:
            self.logger.warning('function Cv_pt is not available in region 4')
            return float("NaN")
        elif region == 5:
            return self.unitConverter.fromSIunit_Cv(Region5.Cv5_pT(p, T))
        else:
            self.logger.warning("Region switch Cv_pt returned unknown value: {:d}".format(region))
            return float("NaN")

    def Cv_ph(self, p, h):
        """Specific isochoric heat capacity as a function of pressure and enthalpy

        Args:
            p (float): preasure value
            h (float): enthalpy value

        Returns:
            Cv (float): specific isochoric heat capacity or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        region = RegionSelection.region_ph(p, h)
        if region == 1:
            Ts = Region1.T1_ph(p, h)
            return self.unitConverter.fromSIunit_Cv(Region1.Cv1_pT(p, Ts))
        elif region == 2:
            Ts = Region2.T2_ph(p, h)
            return self.unitConverter.fromSIunit_Cv(Region2.Cv2_pT(p, Ts))
        elif region == 3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self.unitConverter.fromSIunit_Cv(Region3.Cv3_rhoT(rhos, Ts))
        elif region == 4:
            self.logger.warning('function Cv_ph is not available in region 4')
            return float("NaN")
        elif region == 5:
            Ts = Region5.T5_ph(p, h)
            return self.unitConverter.fromSIunit_Cv(Region5.Cv5_pT(p, Ts))
        else:
            self.logger.warning("Region switch Cv_ph returned unknown value: {:d}".format(region))
            return float("NaN")

    def Cv_ps(self, p, s):
        """Specific isochoric heat capacity as a function of pressure and entropy

        Args:
            p (float): preasure value
            s (float): entropy value

        Returns:
            Cv (float): specific isochoric heat capacity or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        region = RegionSelection.region_ps(p, s)
        if region == 1:
            Ts = Region1.T1_ps(p, s)
            return self.unitConverter.fromSIunit_Cv(Region1.Cv1_pT(p, Ts))
        elif region == 2:
            Ts = Region2.T2_ps(p, s)
            return self.unitConverter.fromSIunit_Cv(Region2.Cv2_pT(p, Ts))
        elif region == 3:
            rhos = 1 / Region3.v3_ps(p, s)
            Ts = Region3.T3_ps(p, s)
            return self.unitConverter.fromSIunit_Cv(Region3.Cv3_rhoT(rhos, Ts))
        elif region == 4:
            self.logger.warning('function Cv_ps is not available in region 4')
            return float("NaN")  # (xs * CvVp + (1 - xs) * CvLp) / Cv_scale - Cv_offset
        elif region == 5:
            Ts = Region5.T5_ps(p, s)
            return self.unitConverter.fromSIunit_Cv(Region5.Cv5_pT(p, Ts))
        else:
            self.logger.warning("Region switch Cv_ps returned unknown value: {:d}".format(region))
            return float("NaN")

# ***********************************************************************************************************
# Section 1.11 Speed of sound
    def wV_p(self, p):
        """Saturated vapour speed of sound as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            wV (float): speed of sound in saturated vapour or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_w(Region2.w2_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_w(Region3.w3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p)))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

    def wL_p(self, p):
        """Saturated liquid speed of sound as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            wL (float): speed of sound in saturated liquid or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_w(Region1.w1_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_w(Region3.w3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p)))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

    def wV_t(self, t):
        """Saturated vapour speed of sound as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            wV (float): speed of sound in saturated vapour or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if (T <= 623.15):
                return self.unitConverter.fromSIunit_w(Region2.w2_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_w(Region3.w3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T)))), T))
        else:
            self.logger.warning('temperature out of range')
            return float("NaN")

    def wL_t(self, t):
        """Saturated liquid speed of sound as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            wL (float): speed of sound in saturated liquid or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if (T <= 623.15):
                return self.unitConverter.fromSIunit_w(Region1.w1_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_w(Region3.w3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T)))), T))
        else:
            self.logger.warning('temperature out of range')
            return float("NaN")

    def w_pt(self, p, t):
        """Speed of sound as a function of pressure and temperature

        Args:
            p (float): preasure value
            t (float): temperature value

        Returns:
            w (float): speed of sound or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        region = RegionSelection.region_pT(p, T)
        if region == 1:
            return self.unitConverter.fromSIunit_w(Region1.w1_pT(p, T))
        elif region == 2:
            return self.unitConverter.fromSIunit_w(Region2.w2_pT(p, T))
        elif region == 3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self.unitConverter.fromSIunit_w(Region3.w3_rhoT(rhos, T))
        elif region == 4:
            self.logger.warning('function w_pt is not available in region 4')
            return float("NaN")
        elif region == 5:
            return self.unitConverter.fromSIunit_w(Region5.w5_pT(p, T))
        else:
            self.logger.warning("Region switch w_pt returned unknown value: {:d}".format(region))
            return float("NaN")

    def w_ph(self, p, h):
        """Speed of sound as a function of pressure and enthalpy

        Args:
            p (float): preasure value
            h (float): enthalpy value

        Returns:
            w (float): speed of sound or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        region = RegionSelection.region_ph(p, h)
        if region == 1:
            Ts = Region1.T1_ph(p, h)
            return self.unitConverter.fromSIunit_w(Region1.w1_pT(p, Ts))
        elif region == 2:
            Ts = Region2.T2_ph(p, h)
            return self.unitConverter.fromSIunit_w(Region2.w2_pT(p, Ts))
        elif region == 3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self.unitConverter.fromSIunit_w(Region3.w3_rhoT(rhos, Ts))
        elif region == 4:
            self.logger.warning('function w_ph is not available in region 4')
            return float("NaN")
        elif region == 5:
            Ts = Region5.T5_ph(p, h)
            return self.unitConverter.fromSIunit_w(Region5.w5_pT(p, Ts))
        else:
            self.logger.warning("Region switch w_ph returned unknown value: {:d}".format(region))
            return float("NaN")

    def w_ps(self, p, s):
        """Speed of sound as a function of pressure and entropy

        Args:
            p (float): preasure value
            s (float): entropy value

        Returns:
            w (float): speed of sound or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        region = RegionSelection.region_ps(p, s)
        if region == 1:
            Ts = Region1.T1_ps(p, s)
            return self.unitConverter.fromSIunit_w(Region1.w1_pT(p, Ts))
        elif region == 2:
            Ts = Region2.T2_ps(p, s)
            return self.unitConverter.fromSIunit_w(Region2.w2_pT(p, Ts))
        elif region == 3:
            rhos = 1 / Region3.v3_ps(p, s)
            Ts = Region3.T3_ps(p, s)
            return self.unitConverter.fromSIunit_w(Region3.w3_rhoT(rhos, Ts))
        elif region == 4:
            self.logger.warning('function w_ps is not available in region 4')
            return float("NaN")  # (xs * wVp + (1 - xs) * wLp) / w_scale - w_offset
        elif region == 5:
            Ts = Region5.T5_ps(p, s)
            return self.unitConverter.fromSIunit_w(Region5.w5_pT(p, Ts))
        else:
            self.logger.warning("Region switch w_ps returned unknown value: {:d}".format(region))
            return float("NaN")

# ***********************************************************************************************************
# Section  1.12 Viscosity
# Viscosity is not part of IAPWS Steam IF97. Equations from
# "Revised Release on the IAPWS Formulation 1985 for the Viscosity of Ordinary Water Substance", 2003 are used.
# Viscosity in the mixed region (4) is interpolated according to the density. This is not true since it will be two fases.
    def my_pt(self, p, t):
        """Viscosity as a function of pressure and temperature

        Args:
            p (float): preasure value
            t (float): temperature value

        Returns:
            my (float): viscosity or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        region = RegionSelection.region_pT(p, T)
        if region == 4:
            self.logger.warning('function my_pt is not available in region 4')
            return float("NaN")
        elif region in [1, 2, 3, 5]:
            return self.unitConverter.fromSIunit_my(TransportProperties.my_AllRegions_pT(p, T))
        else:
            self.logger.warning("Region switch my_pt returned unknown value: {:d}".format(region))
            return float("NaN")

    def my_ph(self, p, h):
        """Viscosity as a function of pressure and enthalpy

        Args:
            p (float): preasure value
            h (float): enthalpy value

        Returns:
            my (float): viscosity or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        region = RegionSelection.region_ph(p, h)
        if region in [1, 2, 3, 5]:
            return self.unitConverter.fromSIunit_my(TransportProperties.my_AllRegions_ph(p, h))
        elif region == 4:
            self.logger.warning('function my_pt is not available in region 4')
            return float("NaN")
        else:
            self.logger.warning("Region switch my_ph returned unknown value: {:d}".format(region))
            return float("NaN")

    def my_ps(self, p, s):
        """Viscosity as a function of pressure and entropy

        Args:
            p (float): preasure value
            s (float): entropy value

        Returns:
            my (float): viscosity
        """
        h = self.h_ps(p, s)
        return self.my_ph(p, h)

# ***********************************************************************************************************
# Section 1.13 Prandtl
    def pr_pt(self, p, t):
        """Prandtl number as a function of preasuere and temperature

        Args:
            p (float): preasure value
            t (float): temperature value

        Returns:
            pr (float): prandtl number
        """
        Cp = self.unitConverter.toSIunit_Cp(self.Cp_pt(p, t))
        my = self.unitConverter.toSIunit_my(self.my_pt(p, t))
        tc = self.unitConverter.toSIunit_tc(self.tc_pt(p, t))
        return Cp * 1000 * my / tc

    def pr_ph(self, p, h):
        """Prandtl number as a function of preasuere and enthalpy

        Args:
            p (float): preasure value
            h (float): enthalpy value

        Returns:
            pr (float): prandtl number
        """
        Cp = self.unitConverter.toSIunit_Cp(self.Cp_ph(p, h))
        my = self.unitConverter.toSIunit_my(self.my_ph(p, h))
        tc = self.unitConverter.toSIunit_tc(self.tc_ph(p, h))
        return Cp * 1000 * my / tc

# ***********************************************************************************************************
# Section 1.15 Surface tension
    def st_t(self, t):
        """Surface tension for two phase water/steam as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            st (float): surface tension
        """
        T = self.unitConverter.toSIunit_T(t)
        return self.unitConverter.fromSIunit_st(TransportProperties.Surface_Tension_T(T))

    def st_p(self, p):
        """Surface tension for two phase water/steam as a function of preasure

        Args:
            p (float): preasure value

        Returns:
            st (float): surface tension
        """
        T = self.tsat_p(p)
        T = self.unitConverter.toSIunit_T(T)
        return self.unitConverter.fromSIunit_st(TransportProperties.Surface_Tension_T(T))

# ***********************************************************************************************************
# Section 1.16 Thermal conductivity
# Revised release on the IAPS Formulation 1985 for the Thermal Conductivity of ordinary water substance (IAPWS 1998)
    def tcL_p(self, p):
        """Saturated liquid thermal conductivity as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            tcL (float): saturated liquid thermal conductivity
        """
        t = self.tsat_p(p)
        v = self.vL_p(p)
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        v = self.unitConverter.toSIunit_v(v)
        rho = 1.0 / v
        return self.unitConverter.fromSIunit_tc(TransportProperties.tc_ptrho(p, T, rho))

    def tcV_p(self, p):
        """Saturated vapour thermal conductivity as a function of pressure

        Args:
            p (float): preasure value

        Returns:
            tcV (float): saturated vapour thermal conductivity
        """
        ps = p
        T = self.tsat_p(p)
        v = self.vV_p(ps)
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(T)
        v = self.unitConverter.toSIunit_v(v)
        rho = 1.0 / v
        return self.unitConverter.fromSIunit_tc(TransportProperties.tc_ptrho(p, T, rho))

    def tcL_t(self, t):
        """Saturated vapour thermal conductivity as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            tcL (float): saturated liquid thermal conductivity
        """
        Ts = t
        p = self.psat_t(Ts)
        v = self.vL_t(Ts)
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(Ts)
        v = self.unitConverter.toSIunit_v(v)
        rho = 1 / v
        return self.unitConverter.fromSIunit_tc(TransportProperties.tc_ptrho(p, T, rho))

    def tcV_t(self, t):
        """Saturated liquid thermal conductivity as a function of temperature

        Args:
            t (float): temperature value

        Returns:
            tcV (float): saturated vapour thermal conductivity
        """
        Ts = t
        p = self.psat_t(Ts)
        v = self.vV_t(Ts)
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(Ts)
        v = self.unitConverter.toSIunit_v(v)
        rho = 1 / v
        return self.unitConverter.fromSIunit_tc(TransportProperties.tc_ptrho(p, T, rho))

    def tc_pt(self, p, t):
        """Thermal conductivity as a function of pressure and temperature

        Args:
            p (float): preasure value
            t (float): temperature value

        Returns:
            tc (float): thermal conductivity
        """

        Ts = t
        ps = p
        v = self.v_pt(ps, Ts)
        p = self.unitConverter.toSIunit_p(ps)
        T = self.unitConverter.toSIunit_T(Ts)
        v = self.unitConverter.toSIunit_v(v)
        rho = 1 / v

        return self.unitConverter.fromSIunit_tc(TransportProperties.tc_ptrho(p, T, rho))

    def tc_ph(self, p, h):
        """Thermal conductivity as a function of pressure and enthalpy

        Args:
            p (float): preasure value
            h (float): enthalpy value

        Returns:
            tc (float): thermal conductivity
        """
        hs = h
        ps = p
        v = self.v_ph(ps, hs)
        T = self.t_ph(ps, hs)
        p = self.unitConverter.toSIunit_p(ps)
        T = self.unitConverter.toSIunit_T(T)
        v = self.unitConverter.toSIunit_v(v)
        rho = 1 / v
        return self.unitConverter.fromSIunit_tc(TransportProperties.tc_ptrho(p, T, rho))

    def tc_hs(self, h, s):
        """Thermal conductivity as a function of enthalpy and entropy

        Args:
            h (float): enthalpy value
            s (float): entropy value

        Returns:
            tc (float): thermal conductivity
        """
        hs = h
        p = self.p_hs(hs, s)
        ps = p
        v = self.v_ph(ps, hs)
        T = self.t_ph(ps, hs)
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(T)
        v = self.unitConverter.toSIunit_v(v)
        rho = 1 / v
        return self.unitConverter.fromSIunit_tc(TransportProperties.tc_ptrho(p, T, rho))

# ***********************************************************************************************************
# Section 1.17 Vapour fraction
    def x_ph(self, p, h):
        """Vapour fraction as a function of pressure and enthalpy

        Args:
            p (float): preasure value
            h (float): enthalpy value

        Returns:
            x (float): vapour fraction or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        if (p > 0.000611657) and (p < 22.06395):
            return self.unitConverter.fromSIunit_x(Region4.x4_ph(p, h))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

    def x_ps(self, p, s):
        """Vapour fraction as a function of pressure and entropy

        Args:
            p (float): preasure value
            s (float): entropy value

        Returns:
            x (float): vapour fraction or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        if (p > 0.000611657) and (p < 22.06395):
            return self.unitConverter.fromSIunit_x(Region4.x4_ps(p, s))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

# ***********************************************************************************************************
# Section 1.18 Vapour Volume Fraction
    def vx_ph(self, p, h):
        """Vapour volume fraction as a function of pressure and enthalpy

        Args:
            p (float): preasure value
            h (float): enthalpy value

        Returns:
            vx (float): vapour volume fraction or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                vL = Region1.v1_pT(p, Region4.T4_p(p))
                vV = Region2.v2_pT(p, Region4.T4_p(p))
            else:
                vL = Region3.v3_ph(p, Region4.h4L_p(p))
                vV = Region3.v3_ph(p, Region4.h4V_p(p))
            xs = Region4.x4_ph(p, h)
            return self.unitConverter.fromSIunit_vx((xs * vV / (xs * vV + (1 - xs) * vL)))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

    def vx_ps(self, p, s):
        """Vapour volume fraction as a function of pressure and entropy

        Args:
            p (float): preasure value
            s (float): entropy value

        Returns:
            vx (float): vapour volume fraction or NaN if arguments are out of range
        """
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                vL = Region1.v1_pT(p, Region4.T4_p(p))
                vV = Region2.v2_pT(p, Region4.T4_p(p))
            else:
                vL = Region3.v3_ph(p, Region4.h4L_p(p))
                vV = Region3.v3_ph(p, Region4.h4V_p(p))
            xs = Region4.x4_ps(p, s)
            return self.unitConverter.fromSIunit_vx((xs * vV / (xs * vV + (1 - xs) * vL)))
        else:
            self.logger.warning('pressure out of range')
            return float("NaN")

# ***********************************************************************************************************
# Revised Release on the Pressure along the Melting and Sublimation Curves of Ordinary Water Substance
# Release IAPWS R14-08(2011)
# http://www.iapws.org/relguide/MeltSub2011.pdf
    def pmelt_t(self, t, hint=None):
        """Pressure along the melting curve as a function of temperature. Based on IAPWS R14-08(2011)
        http://www.iapws.org/relguide/MeltSub2011.pdf
        Because of the shape of the meltin curve it is not possible to automaticaly select the correct region
        automaticaly. Therfore the optional hint-parameter is used to tell the function which area you are interested in.
        The hint-values are namend after the ice types.
        XSteam.TYPE_ICE_Ih = 1
        XSteam.TYPE_ICE_III = 3
        XSteam.TYPE_ICE_V = 5
        XSteam.TYPE_ICE_VI = 6
        XSteam.TYPE_ICE_VII = 7
        If the hint is not one of the values above or None(Default), an Exception is raised

        Args:
            t (float): temperatur value
            hint (int): (optional) hint for the selection logic to decide which part of the melting curve to use

        Returns:
            p (float): preassure or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)

        if hint is None:
            if T >= 251.165 and T < 256.164:
                self.logger.error('cant select ice type based on temperatur, hint reqired')
                return float("NaN")
            elif T >= 256.164 and T < 273.31:
                self.logger.debug('chose ice type V based on temperature')
                return self.unitConverter.fromSIunit_p(IAPWS_R14.pmelt_T_iceV(T))
            elif T >= 273.31 and T < 355:
                self.logger.debug('chose ice type VI based on temperature')
                return self.unitConverter.fromSIunit_p(IAPWS_R14.pmelt_T_iceVI(T))
            elif T >= 355 and T < 751:
                self.logger.debug('chose ice type VII based on temperature')
                return self.unitConverter.fromSIunit_p(IAPWS_R14.pmelt_T_iceVII(T))
            else:
                self.logger.warning('temperature out of range')
                return float("NaN")

        elif hint is self.TYPE_ICE_Ih:
            if T >= 251.165 and T < 273.16:
                return self.unitConverter.fromSIunit_p(IAPWS_R14.pmelt_T_iceIh(T))
            else:
                self.logger.warning('pressure out of range')
                return float("NaN")

        elif hint is self.TYPE_ICE_III:
            if T >= 251.165 and T < 256.164:
                return self.unitConverter.fromSIunit_p(IAPWS_R14.pmelt_T_iceIII(T))
            else:
                self.logger.warning('pressure out of range')
                return float("NaN")

        elif hint is self.TYPE_ICE_V:
            if T >= 256.164 and T < 273.31:
                return self.unitConverter.fromSIunit_p(IAPWS_R14.pmelt_T_iceV(T))
            else:
                self.logger.warning('pressure out of range')
                return float("NaN")

        elif hint is self.TYPE_ICE_VI:
            if T >= 273.31 and T < 355:
                return self.unitConverter.fromSIunit_p(IAPWS_R14.pmelt_T_iceVI(T))
            else:
                self.logger.warning('pressure out of range')
                return float("NaN")

        elif hint is self.TYPE_ICE_VII:
            if T >= 355 and T < 751:
                return self.unitConverter.fromSIunit_p(IAPWS_R14.pmelt_T_iceVII(T))
            else:
                self.logger.warning('pressure out of range')
                return float("NaN")

        else:
            self.logger.error('unknown hint, cant select ice type')
            raise ValueError('unknown hint')

    def psubl_t(self, t):
        """Pressure along the sublimation curve as a function of temperature. Based on IAPWS R14-08(2011)
        http://www.iapws.org/relguide/MeltSub2011.pdf

        Args:
            t (float): temperatur value

        Returns:
            p (float): preassure or NaN if arguments are out of range
        """
        T = self.unitConverter.toSIunit_T(t)
        if T >= 50 and T < 273.16:
            return self.unitConverter.fromSIunit_p(IAPWS_R14.psubl_T(T))
        else:
            self.logger.warning('temperature out of range')
            return float("NaN")
