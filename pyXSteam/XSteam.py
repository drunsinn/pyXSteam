#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Main module for pyXSteam"""

import math
import logging

from .RegionSelection import (
    select_region_pT,
    select_region_ph,
    select_region_ps,
    select_region_hs,
    select_region_prho,
)
from .Regions import Region1, Region2, Region3, Region4, Region5
from .TransportProperties import (
    my_AllRegions_pT,
    my_AllRegions_ph,
    tc_ptrho,
    surface_tension_T,
)
from .Constants import (
    SPECIFIC_GAS_CONSTANT,
    CRITICAL_TEMPERATURE,
    CRITICAL_DENSITY,
    CRITICAL_PRESSURE,
    TRIPLE_POINT_TEMPERATURE,
    TRIPLE_POINT_PRESSURE,
    FREEZING_TEMPERATURE_H2O,
    UnitSystem,
    IceType,
    DiagramRegion,
)
from .UnitConverter import UnitConverter
from .IAPWS_R14 import (
    pmelt_T_iceIh,
    pmelt_T_iceIII,
    pmelt_T_iceV,
    pmelt_T_iceVI,
    pmelt_T_iceVII,
    psubl_T,
)
from .IAPWS_R12 import my_rhoT
from .IAPWS_R6 import (
    R6_cp_rhoT as _R6_cp_rhoT,
    R6_cv_rhoT as _R6_cv_rhoT,
    R6_delta_T_rhoT as _R6_delta_T_rhoT,
    R6_h_rhoT as _R6_h_rhoT,
    R6_isen_t_p_rhoT as _R6_isen_t_p_rhoT,
    R6_joulethomson_rhoT as _R6_joulethomson_rhoT,
    R6_p_rhoT as _R6_p_rhoT,
    R6_phase_equilibrium_condition as _R6_phase_equilibrium_condition,
    R6_s_rhoT as _R6_s_rhoT,
    R6_sec_virial_coefficients as _R6_sec_virial_coefficients,
    R6_third_virial_coefficients as _R6_third_virial_coefficients,
    R6_u_rhoT as _R6_u_rhoT,
    R6_w_rhoT as _R6_w_rhoT,
)


class XSteam:
    """Main pyXSteam object. Abstract of all other functions to allow auto selection of
    the correct region for each set of parameters.

    :param unitSystem: unit system used for input and output values. For supported values
        see the enum UnitSystem.
    """

    UNIT_SYSTEM_BARE = UnitSystem.BARE
    UNIT_SYSTEM_MKS = UnitSystem.MKS
    UNIT_SYSTEM_FLS = UnitSystem.FLS

    TYPE_ICE_Ih = IceType.Ih
    TYPE_ICE_III = IceType.III
    TYPE_ICE_V = IceType.V
    TYPE_ICE_VI = IceType.VI
    TYPE_ICE_VII = IceType.VII

    def __init__(self, unitSystem: UnitSystem = UnitSystem.BARE) -> None:
        """
        Create a property calculator using the requested unit convention.

        :param unitSystem: input/output unit system; MKS and BARE use bar and
            °C, while FLS uses psi and °F.
        """
        self.logger: logging.Logger = logging.getLogger(__name__)
        self._unit_converter = UnitConverter(unitSystem)
        self.logger.debug("initialised pyXSteam with Unit System %s", self._unit_converter)

    def specificGasConstant(self) -> float:
        """Return the specific gas constant R in kJ/(kg·K)."""
        return SPECIFIC_GAS_CONSTANT

    def criticalTemperatur(self) -> float:
        """Return the critical temperature in the selected temperature unit."""
        return self._unit_converter.fromSIunit_T(CRITICAL_TEMPERATURE)

    def criticalPressure(self) -> float:
        """Return the critical pressure in the selected pressure unit."""
        return self._unit_converter.fromSIunit_p(CRITICAL_PRESSURE)

    def criticalDensity(self) -> float:
        """Return the critical density (kg/m³ for MKS/BARE; lb/ft³ for FLS)."""
        return self._unit_converter.fromSIunit_rho(CRITICAL_DENSITY)

    def triplePointTemperatur(self) -> float:
        """Return triple-point temperature in the selected temperature unit."""
        return self._unit_converter.fromSIunit_T(TRIPLE_POINT_TEMPERATURE)

    def triplePointPressure(self) -> float:
        """Return the triple-point pressure in the selected pressure unit."""
        return self._unit_converter.fromSIunit_p(TRIPLE_POINT_PRESSURE)

    def zeroPointTemperature(self) -> float:
        """Return absolute-zero temperature in the selected temperature unit."""
        return self._unit_converter.fromSIunit_T(0.0)

    def tsat_p(self, p: float) -> float:
        """
        Section 1.2 temperature
        Saturation-temperature as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: saturation temperature or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            return self._unit_converter.fromSIunit_T(Region4.T4_p(p))
        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def tsat_s(self, s: float) -> float:
        """
        Section 1.2 temperature
        Saturation-temperature as a function of entropy

        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: saturation temperature or NaN if arguments are out of range
        """
        s = self._unit_converter.toSIunit_s(s)
        if -0.0001545495919 < s < 9.155759395:
            ps = Region4.p4_s(s)
            return self._unit_converter.fromSIunit_T(Region4.T4_p(ps))

        self.logger.debug("entropy value %f is out of range", s)
        return float("NaN")

    def t_ph(self, p: float, h: float) -> float:
        """
        Section 1.2 temperature
        temperature as a function of pressure and enthalpy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)

        :return: temperature or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        h = self._unit_converter.toSIunit_h(h)
        region = select_region_ph(p, h)
        if region == DiagramRegion.R1:
            return self._unit_converter.fromSIunit_T(Region1.T1_ph(p, h))
        if region == DiagramRegion.R2:
            return self._unit_converter.fromSIunit_T(Region2.T2_ph(p, h))
        if region == DiagramRegion.R3:
            return self._unit_converter.fromSIunit_T(Region3.T3_ph(p, h))
        if region == DiagramRegion.R4:
            return self._unit_converter.fromSIunit_T(Region4.T4_p(p))
        if region == DiagramRegion.R5:
            return self._unit_converter.fromSIunit_T(Region5.T5_ph(p, h))

        self.logger.debug(
            "Region switch t_ph returned unknown value %d for input p %f and h %f",
            region,
            p,
            h,
        )
        return float("NaN")

    def t_ps(self, p: float, s: float) -> float:
        """
        Section 1.2 temperature
        temperature as a function of pressure and entropy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: temperature or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        s = self._unit_converter.toSIunit_s(s)
        region = select_region_ps(p, s)
        if region == DiagramRegion.R1:
            return self._unit_converter.fromSIunit_T(Region1.T1_ps(p, s))
        if region == DiagramRegion.R2:
            return self._unit_converter.fromSIunit_T(Region2.T2_ps(p, s))
        if region == DiagramRegion.R3:
            return self._unit_converter.fromSIunit_T(Region3.T3_ps(p, s))
        if region == DiagramRegion.R4:
            return self._unit_converter.fromSIunit_T(Region4.T4_p(p))
        if region == DiagramRegion.R5:
            return self._unit_converter.fromSIunit_T(Region5.T5_ps(p, s))

        self.logger.debug(
            "Region switch t_ps returned unknown value %d for input p %f and s %f",
            region,
            p,
            s,
        )
        return float("NaN")

    def t_hs(self, h: float, s: float) -> float:
        """
        Section 1.2 temperature
        temperature as a function of enthalpy and entropy

        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)
        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: temperature or NaN if arguments are out of range
        """
        h = self._unit_converter.toSIunit_h(h)
        s = self._unit_converter.toSIunit_s(s)
        region = select_region_hs(h, s)
        if region == DiagramRegion.R1:
            p1 = Region1.p1_hs(h, s)
            return self._unit_converter.fromSIunit_T(Region1.T1_ph(p1, h))
        elif region == DiagramRegion.R2:
            p2 = Region2.p2_hs(h, s)
            return self._unit_converter.fromSIunit_T(Region2.T2_ph(p2, h))
        elif region == DiagramRegion.R3:
            p3 = Region3.p3_hs(h, s)
            return self._unit_converter.fromSIunit_T(Region3.T3_ph(p3, h))
        elif region == DiagramRegion.R4:
            return self._unit_converter.fromSIunit_T(Region4.T4_hs(h, s))
        elif region == DiagramRegion.R5:
            self.logger.debug(
                "functions t_hs is not available in region 5 for input h %f and s %f",
                h,
                s,
            )
            return float("NaN")
        else:
            self.logger.debug(
                "Region switch t_hs returned unknown value %d for input h %f and s %f",
                region,
                h,
                s,
            )
            return float("NaN")

    def psat_s(self, s: float) -> float:
        """
        Section 1.3 Pressure
        Saturation-Pressure as a function of entropy

        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: saturation pressure or NaN if arguments are out of range
        """
        s = self._unit_converter.toSIunit_s(s)
        if -0.0001545495919 < s < 9.155759395:
            return self._unit_converter.fromSIunit_p(Region4.p4_s(s))

        self.logger.debug("entropy value %f out of range", s)
        return float("NaN")

    def psat_t(self, t: float) -> float:
        """
        Section 1.3 Pressure
        Saturation-Pressure as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturation pressure or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if 273.1 < T < CRITICAL_TEMPERATURE:
            return self._unit_converter.fromSIunit_p(Region4.p4_T(T))
        else:
            self.logger.debug("temperature %f out of range", T)
            return float("NaN")

    def p_hs(self, h: float, s: float) -> float:
        """
        Section 1.3 Pressure
        Pressure  as a function of enthalpy and entropy

        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)
        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: pressure or NaN if arguments are out of range
        """
        h = self._unit_converter.toSIunit_h(h)
        s = self._unit_converter.toSIunit_s(s)
        region = select_region_hs(h, s)
        if region == DiagramRegion.R1:
            return self._unit_converter.fromSIunit_p(Region1.p1_hs(h, s))
        elif region == DiagramRegion.R2:
            return self._unit_converter.fromSIunit_p(Region2.p2_hs(h, s))
        elif region == DiagramRegion.R3:
            return self._unit_converter.fromSIunit_p(Region3.p3_hs(h, s))
        elif region == DiagramRegion.R4:
            tSat = Region4.T4_hs(h, s)
            return self._unit_converter.fromSIunit_p(Region4.p4_T(tSat))
        elif region == DiagramRegion.R5:
            self.logger.debug(
                "functions p_hs is not available in region 5 for input h %f and s %f",
                h,
                s,
            )
            return float("NaN")
        else:
            self.logger.debug(
                "Region switch p_hs returned unknown value %d for input h %f and s %f",
                region,
                h,
                s,
            )
            return float("NaN")

    def p_hrho(self, h: float, rho: float) -> float:
        """
        Section 1.3 Pressure
        Pressure as a function of h and rho.

        Very unaccurate for solid water region since it's almost incompressible!

        Not valid for water or sumpercritical since water rho does not change very much with p. Uses iteration to find p.

        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)
        :param rho: density in the selected system (kg/m³ for MKS/BARE; lb/ft³ for FLS)

        :raises ValueError: value of density is zero or negative

        :return: pressure or NaN if arguments are out of range
        """
        if rho <= 0.0:
            raise ValueError("rho out of range")
        h = self._unit_converter.toSIunit_h(h)
        High_Bound = self._unit_converter.fromSIunit_p(100)
        Low_Bound = self._unit_converter.fromSIunit_p(TRIPLE_POINT_PRESSURE)
        ps = self._unit_converter.fromSIunit_p(10)
        rhos = 1 / self.v_ph(ps, h)
        step_counter = 0
        while math.fabs(rho - rhos) > 0.0000001:
            step_counter += 1
            last_rhos = rhos

            rhos = 1 / self.v_ph(ps, h)

            if last_rhos == rho:
                self.logger.warning(
                    "p_hrho stopped iterating after %d steps because values did not converge for input values h %f and rho %f",
                    step_counter,
                    h,
                    rho,
                )
                break

            if rhos >= rho:
                High_Bound = ps
            else:
                Low_Bound = ps
            ps = (Low_Bound + High_Bound) / 2
        return ps

    def hV_p(self, p: float) -> float:
        """
        Section 1.4 enthalpy
        Saturated vapor enthalpy as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: saturated vapor enthalpy or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            return self._unit_converter.fromSIunit_h(Region4.h4V_p(p))
        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def hL_p(self, p: float) -> float:
        """
        Section 1.4 enthalpy
        Saturated liquid enthalpy as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: saturated liquid enthalpy or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            return self._unit_converter.fromSIunit_h(Region4.h4L_p(p))
        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def hV_t(self, t: float) -> float:
        """
        Section 1.4 enthalpy
        Saturated vapor enthalpy as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturated vapor enthalpy or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if FREEZING_TEMPERATURE_H2O < T < CRITICAL_TEMPERATURE:
            p = Region4.p4_T(T)
            return self._unit_converter.fromSIunit_h(Region4.h4V_p(p))
        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def hL_t(self, t: float) -> float:
        """
        Section 1.4 enthalpy
        Saturated liquid enthalpy as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturated liquid enthalpy or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if FREEZING_TEMPERATURE_H2O < T < CRITICAL_TEMPERATURE:
            p = Region4.p4_T(T)
            return self._unit_converter.fromSIunit_h(Region4.h4L_p(p))
        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def h_pt(self, p: float, t: float) -> float:
        """
        Section 1.4 enthalpy
        Enthalpy as a function of pressure and temperature

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: enthalpy or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        T = self._unit_converter.toSIunit_T(t)
        region = select_region_pT(p, T)
        if region == DiagramRegion.R1:
            return self._unit_converter.fromSIunit_h(Region1.h1_pT(p, T))
        if region == DiagramRegion.R2:
            return self._unit_converter.fromSIunit_h(Region2.h2_pT(p, T))
        if region == DiagramRegion.R3:
            return self._unit_converter.fromSIunit_h(Region3.h3_pT(p, T))
        if region == DiagramRegion.R4:
            self.logger.debug(
                "function h_pt is not available in region 4 for input p %f and T %f",
                p,
                T,
            )
            return float("NaN")
        if region == DiagramRegion.R5:
            return self._unit_converter.fromSIunit_h(Region5.h5_pT(p, T))

        self.logger.debug(
            "Region switch h_pt returned unknown value %d for input p %f and T %f",
            region,
            p,
            T,
        )
        return float("NaN")

    def h_ps(self, p: float, s: float) -> float:
        """
        Section 1.4 enthalpy
        Enthalpy as a function of pressure and entropy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: enthalpy or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        s = self._unit_converter.toSIunit_s(s)
        region = select_region_ps(p, s)
        if region == DiagramRegion.R1:
            return self._unit_converter.fromSIunit_h(Region1.h1_pT(p, Region1.T1_ps(p, s)))
        if region == DiagramRegion.R2:
            return self._unit_converter.fromSIunit_h(Region2.h2_pT(p, Region2.T2_ps(p, s)))
        if region == DiagramRegion.R3:
            return self._unit_converter.fromSIunit_h(Region3.h3_rhoT(1 / Region3.v3_ps(p, s), Region3.T3_ps(p, s)))
        if region == DiagramRegion.R4:
            xs = Region4.x4_ps(p, s)
            return self._unit_converter.fromSIunit_h(xs * Region4.h4V_p(p) + (1 - xs) * Region4.h4L_p(p))
        if region == DiagramRegion.R5:
            return self._unit_converter.fromSIunit_h(Region5.h5_pT(p, Region5.T5_ps(p, s)))

        self.logger.debug(
            "Region switch h_ps returned unknown value %d for input p %f and s %f ",
            region,
            p,
            s,
        )
        return float("NaN")

    def h_px(self, p: float, x: float) -> float:
        """
        Section 1.4 enthalpy
        Enthalpy as a function of pressure and vapor fraction

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param x: vapor quality, dimensionless mass fraction in [0, 1]

        :return: enthalpy or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        x = self._unit_converter.toSIunit_x(x)
        if (x > 1) or (x < 0) or (p >= 22.064):
            self.logger.debug("Vapor fraction %f and/or pressure %f out of range", x, p)
            return float("NaN")
        hL = Region4.h4L_p(p)
        hV = Region4.h4V_p(p)
        return self._unit_converter.fromSIunit_h(hL + x * (hV - hL))

    def h_prho(self, p: float, rho: float) -> float:
        """
        Section 1.4 enthalpy
        Enthalpy as a function of pressure and density. Observe for low temperatures (liquid) this equation has 2 solutions

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param rho: density in the selected system (kg/m³ for MKS/BARE; lb/ft³ for FLS)

        :return: enthalpy or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        rho = 1 / self._unit_converter.toSIunit_v(1 / float(rho))
        region = select_region_prho(p, rho)
        if region == DiagramRegion.R1:
            return self._unit_converter.fromSIunit_h(Region1.h1_pT(p, Region1.T1_prho(p, rho)))
        if region == DiagramRegion.R2:
            return self._unit_converter.fromSIunit_h(Region2.h2_pT(p, Region2.T2_prho(p, rho)))
        if region == DiagramRegion.R3:
            return self._unit_converter.fromSIunit_h(Region3.h3_rhoT(rho, Region3.T3_prho(p, rho)))
        if region == DiagramRegion.R4:
            if p < 16.529:
                vV = Region2.v2_pT(p, Region4.T4_p(p))
                vL = Region1.v1_pT(p, Region4.T4_p(p))
            else:
                vV = Region3.v3_ph(p, Region4.h4V_p(p))
                vL = Region3.v3_ph(p, Region4.h4L_p(p))
            hV = Region4.h4V_p(p)
            hL = Region4.h4L_p(p)
            x = (1 / rho - vL) / (vV - vL)
            return self._unit_converter.fromSIunit_h((1 - x) * hL + x * hV)
        if region == DiagramRegion.R5:
            return self._unit_converter.fromSIunit_h(Region5.h5_pT(p, Region5.T5_prho(p, rho)))

        self.logger.debug(
            "Region switch h_prho returned unknown value %d for input p %f and rho %f",
            region,
            p,
            rho,
        )
        return float("NaN")

    def h_tx(self, t: float, x: float) -> float:
        """
        Section 1.4 enthalpy
        Enthalpy as a function of temperature and vapor fraction

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :param x: vapor quality, dimensionless mass fraction in [0, 1]

        :return: enthalpy or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        x = self._unit_converter.toSIunit_x(x)
        if (x > 1) or (x < 0) or (T >= CRITICAL_TEMPERATURE):
            self.logger.debug("Vapor fraction %f and/or temperature %f out of range", x, T)
            return float("NaN")
        p = Region4.p4_T(T)
        hL = Region4.h4L_p(p)
        hV = Region4.h4V_p(p)
        return self._unit_converter.fromSIunit_h(hL + x * (hV - hL))

    def vV_p(self, p: float) -> float:
        """
        Section 1.5 Specific Volume
        Saturated vapor volume as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param x: vapor quality, dimensionless mass fraction in [0, 1]

        :return: saturated vapor volume or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            if p < 16.529:
                return self._unit_converter.fromSIunit_v(Region2.v2_pT(p, Region4.T4_p(p)))
            return self._unit_converter.fromSIunit_v(Region3.v3_ph(p, Region4.h4V_p(p)))
        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def vL_p(self, p: float) -> float:
        """
        Section 1.5 Specific Volume
        Saturated liquid volume as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: saturated liquid volume or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            if p < 16.529:
                return self._unit_converter.fromSIunit_v(Region1.v1_pT(p, Region4.T4_p(p)))
            return self._unit_converter.fromSIunit_v(Region3.v3_ph(p, Region4.h4L_p(p)))
        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def vV_t(self, t: float) -> float:
        """
        Section 1.5 Specific Volume
        Saturated vapor volume as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturated vapor volume or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if FREEZING_TEMPERATURE_H2O < T < CRITICAL_TEMPERATURE:
            if T <= 623.15:
                return self._unit_converter.fromSIunit_v(Region2.v2_pT(Region4.p4_T(T), T))
            return self._unit_converter.fromSIunit_v(Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T))))
        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def vL_t(self, t: float) -> float:
        """
        Section 1.5 Specific Volume
        Saturated liquid volume as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturated liquid volume or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if FREEZING_TEMPERATURE_H2O < T < CRITICAL_TEMPERATURE:
            if T <= 623.15:
                return self._unit_converter.fromSIunit_v(Region1.v1_pT(Region4.p4_T(T), T))
            return self._unit_converter.fromSIunit_v(Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T))))
        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def v_pt(self, p: float, t: float) -> float:
        """
        Section 1.5 Specific Volume
        Specific volume as a function of pressure and temperature

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: specific volume or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        T = self._unit_converter.toSIunit_T(t)
        region = select_region_pT(p, T)
        if region == DiagramRegion.R1:
            return self._unit_converter.fromSIunit_v(Region1.v1_pT(p, T))
        if region == DiagramRegion.R2:
            return self._unit_converter.fromSIunit_v(Region2.v2_pT(p, T))
        if region == DiagramRegion.R3:
            return self._unit_converter.fromSIunit_v(Region3.v3_ph(p, Region3.h3_pT(p, T)))
        if region == DiagramRegion.R4:
            self.logger.debug(
                "function v_pt is not available in region 4 for input p %f and T %f",
                p,
                T,
            )
            return float("NaN")
        if region == DiagramRegion.R5:
            return self._unit_converter.fromSIunit_v(Region5.v5_pT(p, T))

        self.logger.debug(
            "Region switch v_pt returned unknown value %d for input p %f and T %f",
            region,
            p,
            T,
        )
        return float("NaN")

    def v_ph(self, p: float, h: float) -> float:
        """
        Section 1.5 Specific Volume
        Specific volume as a function of pressure and enthalpy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)

        :return: specific volume or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        h = self._unit_converter.toSIunit_h(h)
        region = select_region_ph(p, h)
        if region == DiagramRegion.R1:
            return self._unit_converter.fromSIunit_v(Region1.v1_pT(p, Region1.T1_ph(p, h)))
        if region == DiagramRegion.R2:
            return self._unit_converter.fromSIunit_v(Region2.v2_pT(p, Region2.T2_ph(p, h)))
        if region == DiagramRegion.R3:
            return self._unit_converter.fromSIunit_v(Region3.v3_ph(p, h))
        if region == DiagramRegion.R4:
            xs = Region4.x4_ph(p, h)
            if p < 16.529:
                v4v = Region2.v2_pT(p, Region4.T4_p(p))
                v4L = Region1.v1_pT(p, Region4.T4_p(p))
            else:
                v4v = Region3.v3_ph(p, Region4.h4V_p(p))
                v4L = Region3.v3_ph(p, Region4.h4L_p(p))
            return self._unit_converter.fromSIunit_v((xs * v4v + (1 - xs) * v4L))
        if region == DiagramRegion.R5:
            Ts = Region5.T5_ph(p, h)
            return self._unit_converter.fromSIunit_v(Region5.v5_pT(p, Ts))

        self.logger.debug(
            "Region switch v_ph returned unknown value %d for input p %f and h %f",
            region,
            p,
            h,
        )
        return float("NaN")

    def v_ps(self, p: float, s: float) -> float:
        """
        Section 1.5 Specific Volume
        Specific volume as a function of pressure and entropy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: specific volume or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        s = self._unit_converter.toSIunit_s(s)
        region = select_region_ps(p, s)
        if region == DiagramRegion.R1:
            return self._unit_converter.fromSIunit_v(Region1.v1_pT(p, Region1.T1_ps(p, s)))
        if region == DiagramRegion.R2:
            return self._unit_converter.fromSIunit_v(Region2.v2_pT(p, Region2.T2_ps(p, s)))
        if region == DiagramRegion.R3:
            return self._unit_converter.fromSIunit_v(Region3.v3_ps(p, s))
        if region == DiagramRegion.R4:
            xs = Region4.x4_ps(p, s)
            if p < 16.529:
                v4v = Region2.v2_pT(p, Region4.T4_p(p))
                v4L = Region1.v1_pT(p, Region4.T4_p(p))
            else:
                v4v = Region3.v3_ph(p, Region4.h4V_p(p))
                v4L = Region3.v3_ph(p, Region4.h4L_p(p))
            return self._unit_converter.fromSIunit_v((xs * v4v + (1 - xs) * v4L))
        if region == DiagramRegion.R5:
            Ts = Region5.T5_ps(p, s)
            return self._unit_converter.fromSIunit_v(Region5.v5_pT(p, Ts))

        self.logger.debug(
            "Region switch v_ps returned unknown value %d for input p %f and s %f",
            region,
            p,
            s,
        )
        return float("NaN")

    def rhoV_p(self, p: float) -> float:
        """
        Section 1.6 Density (rho)
        Note: Density is calculated as 1/v. See section 1.5 Volume
        Saturated vapor density as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: saturated vapor density
        """
        return 1 / self.vV_p(p)

    def rhoL_p(self, p: float) -> float:
        """
        Section 1.6 Density (rho)
        Note: Density is calculated as 1/v. See section 1.5 Volume
        Saturated liquid density as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: saturated liquid density
        """
        return 1 / self.vL_p(p)

    def rhoV_t(self, t: float) -> float:
        """
        Section 1.6 Density (rho)
        Note: Density is calculated as 1/v. See section 1.5 Volume
        Saturated vapor density as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturated vapor density
        """
        return 1 / self.vV_t(t)

    def rhoL_t(self, t: float) -> float:
        """
        Section 1.6 Density (rho)
        Note: Density is calculated as 1/v. See section 1.5 Volume
        Saturated liquid density as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturated liquid density
        """
        return 1 / self.vL_t(t)

    def rho_pt(self, p: float, t: float) -> float:
        """
        Section 1.6 Density (rho)
        Note: Density is calculated as 1/v. See section 1.5 Volume
        Density as a function of pressure and temperature

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: density
        """
        return 1 / self.v_pt(p, t)

    def rho_ph(self, p: float, h: float) -> float:
        """
        Section 1.6 Density (rho)
        Note: Density is calculated as 1/v. See section 1.5 Volume
        Density as a function of pressure and enthalpy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)

        :return: density
        """
        return 1 / self.v_ph(p, h)

    def rho_ps(self, p: float, s: float) -> float:
        """
        Section 1.6 Density (rho)
        Note: Density is calculated as 1/v. See section 1.5 Volume
        Density as a function of pressure and entropy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: density
        """
        return 1 / self.v_ps(p, s)

    def sV_p(self, p: float) -> float:
        """
        Section 1.7 Specific entropy
        Saturated vapor entropy as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: Saturated vapor entropy or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            if p < 16.529:
                return self._unit_converter.fromSIunit_s(Region2.s2_pT(p, Region4.T4_p(p)))
            return self._unit_converter.fromSIunit_s(Region3.s3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p)))
        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def sL_p(self, p: float) -> float:
        """
        Section 1.7 Specific entropy
        Saturated liquid entropy as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: Saturated liquid entropy or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            if p < 16.529:
                return self._unit_converter.fromSIunit_s(Region1.s1_pT(p, Region4.T4_p(p)))
            return self._unit_converter.fromSIunit_s(Region3.s3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p)))
        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def sV_t(self, t: float) -> float:
        """
        Section 1.7 Specific entropy
        Saturated vapor entropy as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: Saturated vapor entropy or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if FREEZING_TEMPERATURE_H2O < T < CRITICAL_TEMPERATURE:
            if T <= 623.15:
                return self._unit_converter.fromSIunit_s(Region2.s2_pT(Region4.p4_T(T), T))
            return self._unit_converter.fromSIunit_s(
                Region3.s3_rhoT(
                    1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T)))),
                    T,
                )
            )
        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def sL_t(self, t: float) -> float:
        """
        Section 1.7 Specific entropy
        Saturated liquid entropy as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: Saturated liquid entropy or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if FREEZING_TEMPERATURE_H2O < T < CRITICAL_TEMPERATURE:
            if T <= 623.15:
                return self._unit_converter.fromSIunit_s(Region1.s1_pT(Region4.p4_T(T), T))
            return self._unit_converter.fromSIunit_s(
                Region3.s3_rhoT(
                    1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T)))),
                    T,
                )
            )
        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def s_pt(self, p: float, t: float) -> float:
        """
        Section 1.7 Specific entropy
        Specific entropy as a function of pressure and temperature (Returns saturated vapor entalpy if mixture)

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: entropy or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        T = self._unit_converter.toSIunit_T(t)
        region = select_region_pT(p, T)
        if region == DiagramRegion.R1:
            return self._unit_converter.fromSIunit_s(Region1.s1_pT(p, T))
        if region == DiagramRegion.R2:
            return self._unit_converter.fromSIunit_s(Region2.s2_pT(p, T))
        if region == DiagramRegion.R3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self._unit_converter.fromSIunit_s(Region3.s3_rhoT(rhos, T))
        if region == DiagramRegion.R4:
            self.logger.debug("function s_pt is not available in region 4 (p %f, T %f)", p, T)
            return float("NaN")
        if region == DiagramRegion.R5:
            return self._unit_converter.fromSIunit_s(Region5.s5_pT(p, T))

        self.logger.debug(
            "Region switch s_pt returned unknown value %d for input p %f and T %f",
            region,
            p,
            T,
        )
        return float("NaN")

    def s_ph(self, p: float, h: float) -> float:
        """
        Section 1.7 Specific entropy
        Specific entropy as a function of pressure and enthalpy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)

        :return: entropy or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        h = self._unit_converter.toSIunit_h(h)
        region = select_region_ph(p, h)
        if region == DiagramRegion.R1:
            T = Region1.T1_ph(p, h)
            return self._unit_converter.fromSIunit_s(Region1.s1_pT(p, T))
        if region == DiagramRegion.R2:
            T = Region2.T2_ph(p, h)
            return self._unit_converter.fromSIunit_s(Region2.s2_pT(p, T))
        if region == DiagramRegion.R3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self._unit_converter.fromSIunit_s(Region3.s3_rhoT(rhos, Ts))
        if region == DiagramRegion.R4:
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
            return self._unit_converter.fromSIunit_s((xs * s4v + (1 - xs) * s4L))
        if region == DiagramRegion.R5:
            T = Region5.T5_ph(p, h)
            return self._unit_converter.fromSIunit_s(Region5.s5_pT(p, T))
        self.logger.debug(
            "Region switch s_ph returned unknown value %d for input p %f and h %f",
            region,
            p,
            h,
        )
        return float("NaN")

    def uV_p(self, p: float) -> float:
        """
        Section 1.8 Specific internal energy
        Saturated vapor internal energy as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: saturated vapor internal energy or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            if p < 16.529:
                return self._unit_converter.fromSIunit_u(Region2.u2_pT(p, Region4.T4_p(p)))
            return self._unit_converter.fromSIunit_u(Region3.u3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p)))
        self.logger.debug("preasure %f out of range for uV_p(p)", p)
        return float("NaN")

    def uL_p(self, p: float) -> float:
        """
        Section 1.8 Specific internal energy
        Saturated liquid internal energy as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: saturated liquid internal energy
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            if p < 16.529:
                return self._unit_converter.fromSIunit_u(Region1.u1_pT(p, Region4.T4_p(p)))
            return self._unit_converter.fromSIunit_u(Region3.u3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p)))
        self.logger.debug("pressure %f out of range for uL_p(p)", p)
        return float("NaN")

    def uV_t(self, t: float) -> float:
        """
        Section 1.8 Specific internal energy
        Saturated vapor internal energy as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturated vapor internal energy or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if FREEZING_TEMPERATURE_H2O < T < CRITICAL_TEMPERATURE:
            if T <= 623.15:
                return self._unit_converter.fromSIunit_u(Region2.u2_pT(Region4.p4_T(T), T))
            return self._unit_converter.fromSIunit_u(
                Region3.u3_rhoT(
                    1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T)))),
                    T,
                )
            )
        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def uL_t(self, t: float) -> float:
        """
        Section 1.8 Specific internal energy
        Saturated liquid internal energy as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturated liquid internal energy or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if FREEZING_TEMPERATURE_H2O < T < CRITICAL_TEMPERATURE:
            if T <= 623.15:
                return self._unit_converter.fromSIunit_u(Region1.u1_pT(Region4.p4_T(T), T))

            return self._unit_converter.fromSIunit_u(
                Region3.u3_rhoT(
                    1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T)))),
                    T,
                )
            )

        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def u_pt(self, p: float, t: float) -> float:
        """
        Section 1.8 Specific internal energy
        Specific internal energy as a function of pressure and temperature

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: specific internal energy or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        T = self._unit_converter.toSIunit_T(t)
        region = select_region_pT(p, T)

        if region == DiagramRegion.R1:
            return self._unit_converter.fromSIunit_u(Region1.u1_pT(p, T))

        if region == DiagramRegion.R2:
            return self._unit_converter.fromSIunit_u(Region2.u2_pT(p, T))

        if region == DiagramRegion.R3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self._unit_converter.fromSIunit_u(Region3.u3_rhoT(rhos, T))

        if region == DiagramRegion.R4:
            self.logger.debug(
                "function u_pt is not available in region 4 for input p %f and T %f",
                p,
                T,
            )
            return float("NaN")

        if region == DiagramRegion.R5:
            return self._unit_converter.fromSIunit_u(Region5.u5_pT(p, T))

        self.logger.debug(
            "Region switch u_pt returned unknown value %d for input p %f and T %f",
            region,
            p,
            T,
        )
        return float("NaN")

    def u_ph(self, p: float, h: float) -> float:
        """
        Section 1.8 Specific internal energy
        Specific internal energy as a function of pressure and enthalpy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)

        :return: specific internal energy or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        h = self._unit_converter.toSIunit_h(h)
        region = select_region_ph(p, h)

        if region == DiagramRegion.R1:
            Ts = Region1.T1_ph(p, h)
            return self._unit_converter.fromSIunit_u(Region1.u1_pT(p, Ts))

        if region == DiagramRegion.R2:
            Ts = Region2.T2_ph(p, h)
            return self._unit_converter.fromSIunit_u(Region2.u2_pT(p, Ts))

        if region == DiagramRegion.R3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self._unit_converter.fromSIunit_u(Region3.u3_rhoT(rhos, Ts))

        if region == DiagramRegion.R4:
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
            return self._unit_converter.fromSIunit_u((xs * u4v + (1 - xs) * u4L))

        if region == DiagramRegion.R5:
            Ts = Region5.T5_ph(p, h)
            return self._unit_converter.fromSIunit_u(Region5.u5_pT(p, Ts))

        self.logger.debug(
            "Region switch u_ph returned unknown value %d for input p %f and h %f",
            region,
            p,
            h,
        )
        return float("NaN")

    def u_ps(self, p: float, s: float) -> float:
        """
        Section 1.8 Specific internal energy
        Specific internal energy as a function of pressure and entropy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: specific internal energy or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        s = self._unit_converter.toSIunit_s(s)
        region = select_region_ps(p, s)
        if region == DiagramRegion.R1:
            Ts = Region1.T1_ps(p, s)
            return self._unit_converter.fromSIunit_u(Region1.u1_pT(p, Ts))

        if region == DiagramRegion.R2:
            Ts = Region2.T2_ps(p, s)
            return self._unit_converter.fromSIunit_u(Region2.u2_pT(p, Ts))

        if region == DiagramRegion.R3:
            rhos = 1 / Region3.v3_ps(p, s)
            Ts = Region3.T3_ps(p, s)
            return self._unit_converter.fromSIunit_u(Region3.u3_rhoT(rhos, Ts))

        if region == DiagramRegion.R4:
            if p < 16.529:
                uLp = Region1.u1_pT(p, Region4.T4_p(p))
                uVp = Region2.u2_pT(p, Region4.T4_p(p))
            else:
                uLp = Region3.u3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p))
                uVp = Region3.u3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p))
            xs = Region4.x4_ps(p, s)
            return self._unit_converter.fromSIunit_u((xs * uVp + (1 - xs) * uLp))

        if region == DiagramRegion.R5:
            Ts = Region5.T5_ps(p, s)
            return self._unit_converter.fromSIunit_u(Region5.u5_pT(p, Ts))

        self.logger.debug(
            "Region switch u_ps returned unknown value %d for input p %f and s %f",
            region,
            p,
            s,
        )
        return float("NaN")

    def CpV_p(self, p: float) -> float:
        """
        Section 1.9 Specific isobaric heat capacity
        Saturated vapor heat capacity as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: saturated vapor heat capacity or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            if p < 16.529:
                return self._unit_converter.fromSIunit_Cp(Region2.Cp2_pT(p, Region4.T4_p(p)))
            return self._unit_converter.fromSIunit_Cp(Region3.Cp3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p)))
        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def CpL_p(self, p: float) -> float:
        """
        Section 1.9 Specific isobaric heat capacity
        Saturated liquid heat capacity as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: saturated liquid heat capacity or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            if p < 16.529:
                return self._unit_converter.fromSIunit_Cp(Region1.Cp1_pT(p, Region4.T4_p(p)))
            return self._unit_converter.fromSIunit_Cp(Region3.Cp3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p)))
        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def CpV_t(self, t: float) -> float:
        """
        Section 1.9 Specific isobaric heat capacity
        Saturated vapor heat capacity as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturated vapor heat capacity or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if FREEZING_TEMPERATURE_H2O < T < CRITICAL_TEMPERATURE:
            if T <= 623.15:
                return self._unit_converter.fromSIunit_Cp(Region2.Cp2_pT(Region4.p4_T(T), T))
            return self._unit_converter.fromSIunit_Cp(
                Region3.Cp3_rhoT(
                    1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T)))),
                    T,
                )
            )
        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def CpL_t(self, t: float) -> float:
        """
        Section 1.9 Specific isobaric heat capacity
        Saturated liquid heat capacity as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturated liquid heat capacity or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if FREEZING_TEMPERATURE_H2O < T < CRITICAL_TEMPERATURE:
            if T <= 623.15:
                Cp = Region1.Cp1_pT(Region4.p4_T(T), T)
                return self._unit_converter.fromSIunit_Cp(Cp)

            p = Region4.p4_T(T)
            rho = 1 / Region3.v3_ph(p, Region4.h4L_p(p))
            Cp = Region3.Cp3_rhoT(rho, T)
            return self._unit_converter.fromSIunit_Cp(Cp)

        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def Cp_pt(self, p: float, t: float) -> float:
        """
        Section 1.9 Specific isobaric heat capacity
        Specific isobaric heat capacity as a function of pressure and temperature

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: specific isobaric heat capacity or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        T = self._unit_converter.toSIunit_T(t)
        region = select_region_pT(p, T)

        if region == DiagramRegion.R1:
            return self._unit_converter.fromSIunit_Cp(Region1.Cp1_pT(p, T))

        if region == DiagramRegion.R2:
            return self._unit_converter.fromSIunit_Cp(Region2.Cp2_pT(p, T))

        if region == DiagramRegion.R3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self._unit_converter.fromSIunit_Cp(Region3.Cp3_rhoT(rhos, T))

        if region == DiagramRegion.R4:
            self.logger.debug(
                "function Cp_pt is not available in region 4 for input p %f and T %f",
                p,
                T,
            )
            return float("NaN")

        if region == DiagramRegion.R5:
            return self._unit_converter.fromSIunit_Cp(Region5.Cp5_pT(p, T))

        self.logger.debug(
            "Region switch Cp_pt returned unknown value %d for input p %f and T %f",
            region,
            p,
            T,
        )
        return float("NaN")

    def Cp_ph(self, p: float, h: float) -> float:
        """
        Section 1.9 Specific isobaric heat capacity
        Specific isobaric heat capacity as a function of pressure and enthalpy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)

        :return: specific isobaric heat capacity or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        h = self._unit_converter.toSIunit_h(h)
        region = select_region_ph(p, h)

        if region == DiagramRegion.R1:
            Ts = Region1.T1_ph(p, h)
            return self._unit_converter.fromSIunit_Cp(Region1.Cp1_pT(p, Ts))

        if region == DiagramRegion.R2:
            Ts = Region2.T2_ph(p, h)
            return self._unit_converter.fromSIunit_Cp(Region2.Cp2_pT(p, Ts))

        if region == DiagramRegion.R3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self._unit_converter.fromSIunit_Cp(Region3.Cp3_rhoT(rhos, Ts))

        if region == DiagramRegion.R4:
            self.logger.debug(
                "function Cp_ph is not available in region 4  for input p %f and h %f",
                p,
                h,
            )
            return float("NaN")

        if region == DiagramRegion.R5:
            Ts = Region5.T5_ph(p, h)
            return self._unit_converter.fromSIunit_Cp(Region5.Cp5_pT(p, Ts))

        self.logger.debug(
            "Region switch Cp_ph returned unknown value %d for input p %f and h %f",
            region,
            p,
            h,
        )
        return float("NaN")

    def Cp_ps(self, p: float, s: float) -> float:
        """
        Section 1.9 Specific isobaric heat capacity
        Specific isobaric heat capacity as a function of pressure and entropy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: specific isobaric heat capacity or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        s = self._unit_converter.toSIunit_s(s)
        region = select_region_ps(p, s)

        if region == DiagramRegion.R1:
            Ts = Region1.T1_ps(p, s)
            return self._unit_converter.fromSIunit_Cp(Region1.Cp1_pT(p, Ts))

        if region == DiagramRegion.R2:
            Ts = Region2.T2_ps(p, s)
            return self._unit_converter.fromSIunit_Cp(Region2.Cp2_pT(p, Ts))

        if region == DiagramRegion.R3:
            rhos = 1 / Region3.v3_ps(p, s)
            Ts = Region3.T3_ps(p, s)
            return self._unit_converter.fromSIunit_Cp(Region3.Cp3_rhoT(rhos, Ts))

        if region == DiagramRegion.R4:
            self.logger.debug(
                "function Cp_ps is not available in region 4 for input p %f and s %f",
                p,
                s,
            )
            return float("NaN")

        if region == DiagramRegion.R5:
            Ts = Region5.T5_ps(p, s)
            return self._unit_converter.fromSIunit_Cp(Region5.Cp5_pT(p, Ts))

        self.logger.debug(
            "Region switch Cp_ps returned unknown value %d for input p %f and s %f",
            region,
            p,
            s,
        )
        return float("NaN")

    def CvV_p(self, p: float) -> float:
        """
        Section 1.10 Specific isochoric heat capacity
        Saturated vapor isochoric heat capacity as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: saturated vapor isochoric heat capacity or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            if p < 16.529:
                Cv = Region2.Cv2_pT(p, Region4.T4_p(p))
                return self._unit_converter.fromSIunit_Cv(Cv)

            rho = 1 / Region3.v3_ph(p, Region4.h4V_p(p))
            Cv = Region3.Cv3_rhoT(rho, Region4.T4_p(p))
            return self._unit_converter.fromSIunit_Cv(Cv)

        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def CvL_p(self, p: float) -> float:
        """
        Section 1.10 Specific isochoric heat capacity
        Saturated liquid isochoric heat capacity as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: saturated liquid isochoric heat capacity or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            if p < 16.529:
                Cv = Region1.Cv1_pT(p, Region4.T4_p(p))
                return self._unit_converter.fromSIunit_Cv(Cv)

            rho = 1 / Region3.v3_ph(p, Region4.h4L_p(p))
            Cv = Region3.Cv3_rhoT(rho, Region4.T4_p(p))
            return self._unit_converter.fromSIunit_Cv(Cv)

        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def CvV_t(self, t: float) -> float:
        """
        Section 1.10 Specific isochoric heat capacity
        Saturated vapor isochoric heat capacity as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturated vapor isochoric heat capacity or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if FREEZING_TEMPERATURE_H2O < T < CRITICAL_TEMPERATURE:
            if T <= 623.15:
                Cv = Region2.Cv2_pT(Region4.p4_T(T), T)
                return self._unit_converter.fromSIunit_Cv(Cv)

            p = Region4.p4_T(T)
            rho = 1 / Region3.v3_ph(p, Region4.h4V_p(p))
            Cv = Region3.Cv3_rhoT(rho, T)
            return self._unit_converter.fromSIunit_Cv(Cv)

        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def CvL_t(self, t: float) -> float:
        """
        Section 1.10 Specific isochoric heat capacity
        Saturated liquid isochoric heat capacity as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturated liquid isochoric heat capacity or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if FREEZING_TEMPERATURE_H2O < T < CRITICAL_TEMPERATURE:
            if T <= 623.15:
                Cv = Region1.Cv1_pT(Region4.p4_T(T), T)
                return self._unit_converter.fromSIunit_Cv(Cv)

            p = Region4.p4_T(T)
            rho = 1 / Region3.v3_ph(p, Region4.h4L_p(p))
            Cv = Region3.Cv3_rhoT(rho, T)
            return self._unit_converter.fromSIunit_Cv(Cv)

        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def Cv_pt(self, p: float, t: float) -> float:
        """
        Section 1.10 Specific isochoric heat capacity
        Specific isochoric heat capacity as a function of pressure and temperature

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: specific isochoric heat capacity or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        T = self._unit_converter.toSIunit_T(t)
        region = select_region_pT(p, T)

        if region == DiagramRegion.R1:
            return self._unit_converter.fromSIunit_Cv(Region1.Cv1_pT(p, T))

        if region == DiagramRegion.R2:
            return self._unit_converter.fromSIunit_Cv(Region2.Cv2_pT(p, T))

        if region == DiagramRegion.R3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self._unit_converter.fromSIunit_Cv(Region3.Cv3_rhoT(rhos, T))

        if region == DiagramRegion.R4:
            self.logger.debug(
                "function Cv_pt is not available in region 4 for input p %f and T %f",
                p,
                T,
            )
            return float("NaN")

        if region == DiagramRegion.R5:
            return self._unit_converter.fromSIunit_Cv(Region5.Cv5_pT(p, T))

        self.logger.debug(
            "Region switch Cv_pt returned unknown value %d for input p %f and T %f",
            region,
            p,
            T,
        )
        return float("NaN")

    def Cv_ph(self, p: float, h: float) -> float:
        """
        Section 1.10 Specific isochoric heat capacity
        Specific isochoric heat capacity as a function of pressure and enthalpy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)

        :return: specific isochoric heat capacity or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        h = self._unit_converter.toSIunit_h(h)
        region = select_region_ph(p, h)

        if region == DiagramRegion.R1:
            Ts = Region1.T1_ph(p, h)
            return self._unit_converter.fromSIunit_Cv(Region1.Cv1_pT(p, Ts))

        if region == DiagramRegion.R2:
            Ts = Region2.T2_ph(p, h)
            return self._unit_converter.fromSIunit_Cv(Region2.Cv2_pT(p, Ts))

        if region == DiagramRegion.R3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self._unit_converter.fromSIunit_Cv(Region3.Cv3_rhoT(rhos, Ts))

        if region == DiagramRegion.R4:
            self.logger.debug(
                "function Cv_ph is not available in region 4 for input p %f and h %f",
                p,
                h,
            )
            return float("NaN")

        if region == DiagramRegion.R5:
            Ts = Region5.T5_ph(p, h)
            return self._unit_converter.fromSIunit_Cv(Region5.Cv5_pT(p, Ts))

        self.logger.debug(
            "Region switch Cv_ph returned unknown value %d for input p %f and h %f",
            region,
            p,
            h,
        )
        return float("NaN")

    def Cv_ps(self, p: float, s: float) -> float:
        """
        Section 1.10 Specific isochoric heat capacity
        Specific isochoric heat capacity as a function of pressure and entropy
        In the two-phase region, returns a quality-weighted interpolation of
        the saturated liquid and vapor values; this is an estimate, not an
        equilibrium-mixture derivative.

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: specific isochoric heat capacity or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        s = self._unit_converter.toSIunit_s(s)
        region = select_region_ps(p, s)

        if region == DiagramRegion.R1:
            Ts = Region1.T1_ps(p, s)
            return self._unit_converter.fromSIunit_Cv(Region1.Cv1_pT(p, Ts))

        if region == DiagramRegion.R2:
            Ts = Region2.T2_ps(p, s)
            return self._unit_converter.fromSIunit_Cv(Region2.Cv2_pT(p, Ts))

        if region == DiagramRegion.R3:
            rhos = 1 / Region3.v3_ps(p, s)
            Ts = Region3.T3_ps(p, s)
            return self._unit_converter.fromSIunit_Cv(Region3.Cv3_rhoT(rhos, Ts))

        if region == DiagramRegion.R4:
            quality = Region4.x4_ps(p, s)
            p_user = self._unit_converter.fromSIunit_p(p)
            cv_vapor = self.CvV_p(p_user)
            cv_liquid = self.CvL_p(p_user)
            return quality * cv_vapor + (1 - quality) * cv_liquid

        if region == DiagramRegion.R5:
            Ts = Region5.T5_ps(p, s)
            return self._unit_converter.fromSIunit_Cv(Region5.Cv5_pT(p, Ts))

        self.logger.debug(
            "Region switch Cv_ps returned unknown value %d for input p %f and s %f",
            region,
            p,
            s,
        )
        return float("NaN")

    def wV_p(self, p: float) -> float:
        """
        Section 1.11 Speed of sound
        Saturated vapor speed of sound as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: speed of sound in saturated vapor or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            if p < 16.529:
                w = Region2.w2_pT(p, Region4.T4_p(p))
                return self._unit_converter.fromSIunit_w(w)

            rho = 1 / Region3.v3_ph(p, Region4.h4V_p(p))
            w = Region3.w3_rhoT(rho, Region4.T4_p(p))
            return self._unit_converter.fromSIunit_w(w)

        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def wL_p(self, p: float) -> float:
        """
        Section 1.11 Speed of sound
        Saturated liquid speed of sound as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: speed of sound in saturated liquid or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            if p < 16.529:
                w = Region1.w1_pT(p, Region4.T4_p(p))
                return self._unit_converter.fromSIunit_w(w)

            rho = 1 / Region3.v3_ph(p, Region4.h4L_p(p))
            w = Region3.w3_rhoT(rho, Region4.T4_p(p))
            return self._unit_converter.fromSIunit_w(w)

        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def wV_t(self, t: float) -> float:
        """
        Section 1.11 Speed of sound
        Saturated vapor speed of sound as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: speed of sound in saturated vapor or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if FREEZING_TEMPERATURE_H2O < T < CRITICAL_TEMPERATURE:
            if T <= 623.15:
                w = Region2.w2_pT(Region4.p4_T(T), T)
                return self._unit_converter.fromSIunit_w(w)

            p = Region4.p4_T(T)
            rho = 1 / Region3.v3_ph(p, Region4.h4V_p(p))
            w = Region3.w3_rhoT(rho, T)
            return self._unit_converter.fromSIunit_w(w)

        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def wL_t(self, t: float) -> float:
        """
        Section 1.11 Speed of sound
        Saturated liquid speed of sound as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: speed of sound in saturated liquid or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if FREEZING_TEMPERATURE_H2O < T < CRITICAL_TEMPERATURE:
            if T <= 623.15:
                w = Region1.w1_pT(Region4.p4_T(T), T)
                return self._unit_converter.fromSIunit_w(w)

            p = Region4.p4_T(T)
            roh = 1 / Region3.v3_ph(p, Region4.h4L_p(p))
            w = Region3.w3_rhoT(roh, T)
            return self._unit_converter.fromSIunit_w(w)

        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def w_pt(self, p: float, t: float) -> float:
        """
        Section 1.11 Speed of sound
        Speed of sound as a function of pressure and temperature

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: speed of sound or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        T = self._unit_converter.toSIunit_T(t)
        region = select_region_pT(p, T)

        if region == DiagramRegion.R1:
            return self._unit_converter.fromSIunit_w(Region1.w1_pT(p, T))

        if region == DiagramRegion.R2:
            return self._unit_converter.fromSIunit_w(Region2.w2_pT(p, T))

        if region == DiagramRegion.R3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self._unit_converter.fromSIunit_w(Region3.w3_rhoT(rhos, T))

        if region == DiagramRegion.R4:
            self.logger.debug(
                "function w_pt is not available in region 4 for input p %f and T %f",
                p,
                T,
            )
            return float("NaN")

        if region == DiagramRegion.R5:
            return self._unit_converter.fromSIunit_w(Region5.w5_pT(p, T))

        self.logger.debug(
            "Region switch w_pt returned unknown value %d for input p %f and T %f",
            region,
            p,
            T,
        )
        return float("NaN")

    def w_ph(self, p: float, h: float) -> float:
        """
        Section 1.11 Speed of sound
        Speed of sound as a function of pressure and enthalpy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)

        :return: speed of sound or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        h = self._unit_converter.toSIunit_h(h)
        region = select_region_ph(p, h)

        if region == DiagramRegion.R1:
            Ts = Region1.T1_ph(p, h)
            return self._unit_converter.fromSIunit_w(Region1.w1_pT(p, Ts))

        if region == DiagramRegion.R2:
            Ts = Region2.T2_ph(p, h)
            return self._unit_converter.fromSIunit_w(Region2.w2_pT(p, Ts))

        if region == DiagramRegion.R3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self._unit_converter.fromSIunit_w(Region3.w3_rhoT(rhos, Ts))

        if region == DiagramRegion.R4:
            self.logger.debug("function w_ph is not available in region 4 p %f and h %f", p, h)
            return float("NaN")

        if region == DiagramRegion.R5:
            Ts = Region5.T5_ph(p, h)
            return self._unit_converter.fromSIunit_w(Region5.w5_pT(p, Ts))

        self.logger.debug(
            "Region switch w_ph returned unknown value %d for input p %f and h %f",
            region,
            p,
            h,
        )
        return float("NaN")

    def w_ps(self, p: float, s: float) -> float:
        """
        Section 1.11 Speed of sound
        Speed of sound as a function of pressure and entropy
        In the two-phase region, returns a quality-weighted interpolation of
        the saturated liquid and vapor values; this is an estimate, not an
        equilibrium-mixture speed of sound.

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: speed of sound or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        s = self._unit_converter.toSIunit_s(s)
        region = select_region_ps(p, s)

        if region == DiagramRegion.R1:
            Ts = Region1.T1_ps(p, s)
            return self._unit_converter.fromSIunit_w(Region1.w1_pT(p, Ts))

        if region == DiagramRegion.R2:
            Ts = Region2.T2_ps(p, s)
            return self._unit_converter.fromSIunit_w(Region2.w2_pT(p, Ts))

        if region == DiagramRegion.R3:
            rhos = 1 / Region3.v3_ps(p, s)
            Ts = Region3.T3_ps(p, s)
            return self._unit_converter.fromSIunit_w(Region3.w3_rhoT(rhos, Ts))

        if region == DiagramRegion.R4:
            quality = Region4.x4_ps(p, s)
            p_user = self._unit_converter.fromSIunit_p(p)
            w_vapor = self.wV_p(p_user)
            w_liquid = self.wL_p(p_user)
            return quality * w_vapor + (1 - quality) * w_liquid

        if region == DiagramRegion.R5:
            Ts = Region5.T5_ps(p, s)
            return self._unit_converter.fromSIunit_w(Region5.w5_pT(p, Ts))

        self.logger.debug(
            "Region switch w_ps returned unknown value %d for input p %f and s %f",
            region,
            p,
            s,
        )
        return float("NaN")

    def my_pt(self, p: float, t: float) -> float:
        """
        Section 1.12 Viscosity
        Note: Viscosity is not part of IAPWS Steam IF97. Equations from "Revised Release on the IAPWS Formulation 1985 for the
        Viscosity of Ordinary Water Substance", 2003 are used. Viscosity in the mixed region (4) is interpolated according to the
        density. This is not true since it will be two phases.
        Viscosity as a function of pressure and temperature

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: viscosity or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        T = self._unit_converter.toSIunit_T(t)
        region = select_region_pT(p, T)
        if region == DiagramRegion.R4:
            self.logger.debug(
                "function my_pt is not available in region 4 for input p %f and T %f",
                p,
                T,
            )
            return float("NaN")

        if region == DiagramRegion.NILL:
            self.logger.debug(
                "Region switch my_pt returned unknown value %d for input p %f and T %f",
                region,
                p,
                T,
            )
            return float("NaN")

        return self._unit_converter.fromSIunit_my(my_AllRegions_pT(p, T))

    def my_ph(self, p: float, h: float) -> float:
        """
        Section 1.12 Viscosity
        Note: Viscosity is not part of IAPWS Steam IF97. Equations from "Revised
        Release on the IAPWS Formulation 1985 for the Viscosity of Ordinary Water
        Substance", 2003 are used. Viscosity in the mixed region (4) is interpolated
        according to the density. This is not true since it will be two phases.
        Viscosity as a function of pressure and enthalpy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)

        :return: viscosity or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        h = self._unit_converter.toSIunit_h(h)
        region = select_region_ph(p, h)

        if region == DiagramRegion.R4:
            self.logger.debug(
                "function my_pt is not available in region 4 for input p %f and h %f",
                p,
                h,
            )
            return float("NaN")

        if region == DiagramRegion.NILL:
            self.logger.debug(
                "Region switch my_ph returned unknown value %d for input p %f and h %f",
                region,
                p,
                h,
            )
            return float("NaN")

        return self._unit_converter.fromSIunit_my(my_AllRegions_ph(p, h))

    def my_ps(self, p: float, s: float) -> float:
        """
        Section 1.12 Viscosity
        Note: Viscosity is not part of IAPWS Steam IF97. Equations
        from "Revised Release on the IAPWS Formulation 1985 for the
        Viscosity of Ordinary Water Substance", 2003 are used. Viscosity in
        the mixed region (4) is interpolated according to the density. This
        is not true since it will be two phases.
        Viscosity as a function of pressure and entropy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: viscosity
        """
        h = self.h_ps(p, s)
        return self.my_ph(p, h)

    def pr_pt(self, p: float, t: float) -> float:
        """
        Section 1.13 Prandtl
        Prandtl number as a function of preasuere and temperature

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: prandtl number
        """
        Cp = self._unit_converter.toSIunit_Cp(self.Cp_pt(p, t))
        my = self._unit_converter.toSIunit_my(self.my_pt(p, t))
        tc = self._unit_converter.toSIunit_tc(self.tc_pt(p, t))
        return Cp * 1000 * my / tc

    def pr_ph(self, p: float, h: float) -> float:
        """
        Section 1.13 Prandtl
        Prandtl number as a function of preasuere and enthalpy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)

        :return: prandtl number
        """
        Cp = self._unit_converter.toSIunit_Cp(self.Cp_ph(p, h))
        my = self._unit_converter.toSIunit_my(self.my_ph(p, h))
        tc = self._unit_converter.toSIunit_tc(self.tc_ph(p, h))
        return Cp * 1000 * my / tc

    def st_t(self, t: float) -> float:
        """
        Section 1.15 Surface tension
        Surface tension for two phase water/steam as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: surface tension
        """
        T = self._unit_converter.toSIunit_T(t)
        return self._unit_converter.fromSIunit_st(surface_tension_T(T))

    def st_p(self, p: float) -> float:
        """
        Section 1.15 Surface tension
        Surface tension for two phase water/steam as a function of preasure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: surface tension
        """
        T = self.tsat_p(p)
        T = self._unit_converter.toSIunit_T(T)
        return self._unit_converter.fromSIunit_st(surface_tension_T(T))

    def tcL_p(self, p: float) -> float:
        """
        Section 1.16 Thermal conductivity
        Note: Revised release on the IAPS Formulation 1985 for the Thermal
        Conductivity of ordinary water substance (IAPWS 1998)
        Saturated liquid thermal conductivity as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: saturated liquid thermal conductivity
        """
        t = self.tsat_p(p)
        v = self.vL_p(p)
        p = self._unit_converter.toSIunit_p(p)
        T = self._unit_converter.toSIunit_T(t)
        v = self._unit_converter.toSIunit_v(v)
        rho = 1.0 / v
        return self._unit_converter.fromSIunit_tc(tc_ptrho(p, T, rho))

    def tcV_p(self, p: float) -> float:
        """
        Section 1.16 Thermal conductivity
        Note: Revised release on the IAPS Formulation 1985 for the Thermal
        Conductivity of ordinary water substance (IAPWS 1998)
        Saturated vapor thermal conductivity as a function of pressure

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)

        :return: saturated vapor thermal conductivity
        """
        ps = p
        T = self.tsat_p(p)
        v = self.vV_p(ps)
        p = self._unit_converter.toSIunit_p(p)
        T = self._unit_converter.toSIunit_T(T)
        v = self._unit_converter.toSIunit_v(v)
        rho = 1.0 / v
        return self._unit_converter.fromSIunit_tc(tc_ptrho(p, T, rho))

    def tcL_t(self, t: float) -> float:
        """
        Section 1.16 Thermal conductivity
        Note: Revised release on the IAPS Formulation 1985 for the Thermal
        Conductivity of ordinary water substance (IAPWS 1998)
        Saturated vapor thermal conductivity as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturated liquid thermal conductivity
        """
        Ts = t
        p = self.psat_t(Ts)
        v = self.vL_t(Ts)
        p = self._unit_converter.toSIunit_p(p)
        T = self._unit_converter.toSIunit_T(Ts)
        v = self._unit_converter.toSIunit_v(v)
        rho = 1 / v
        return self._unit_converter.fromSIunit_tc(tc_ptrho(p, T, rho))

    def tcV_t(self, t: float) -> float:
        """
        Section 1.16 Thermal conductivity
        Note: Revised release on the IAPS Formulation 1985 for the Thermal
        Conductivity of ordinary water substance (IAPWS 1998)
        Saturated liquid thermal conductivity as a function of temperature

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: saturated vapor thermal conductivity
        """
        Ts = t
        p = self.psat_t(Ts)
        v = self.vV_t(Ts)
        p = self._unit_converter.toSIunit_p(p)
        T = self._unit_converter.toSIunit_T(Ts)
        v = self._unit_converter.toSIunit_v(v)
        rho = 1 / v
        return self._unit_converter.fromSIunit_tc(tc_ptrho(p, T, rho))

    def tc_pt(self, p: float, t: float) -> float:
        """
        Section 1.16 Thermal conductivity
        Note: Revised release on the IAPS Formulation 1985 for the Thermal
        Conductivity of ordinary water substance (IAPWS 1998)
        Thermal conductivity as a function of pressure and temperature

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: thermal conductivity
        """

        Ts = t
        ps = p
        v = self.v_pt(ps, Ts)
        p = self._unit_converter.toSIunit_p(ps)
        T = self._unit_converter.toSIunit_T(Ts)
        v = self._unit_converter.toSIunit_v(v)
        rho = 1 / v

        return self._unit_converter.fromSIunit_tc(tc_ptrho(p, T, rho))

    def tc_ph(self, p: float, h: float) -> float:
        """
        Section 1.16 Thermal conductivity
        Note: Revised release on the IAPS Formulation 1985 for the Thermal
        Conductivity of ordinary water substance (IAPWS 1998)
        Thermal conductivity as a function of pressure and enthalpy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)

        :return: thermal conductivity
        """
        hs = h
        ps = p
        v = self.v_ph(ps, hs)
        T = self.t_ph(ps, hs)
        p = self._unit_converter.toSIunit_p(ps)
        T = self._unit_converter.toSIunit_T(T)
        v = self._unit_converter.toSIunit_v(v)
        rho = 1 / v
        return self._unit_converter.fromSIunit_tc(tc_ptrho(p, T, rho))

    def tc_hs(self, h: float, s: float) -> float:
        """
        Section 1.16 Thermal conductivity
        Note: Revised release on the IAPS Formulation 1985 for the Thermal
        Conductivity of ordinary water substance (IAPWS 1998)
        Thermal conductivity as a function of enthalpy and entropy

        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)
        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: thermal conductivity
        """
        hs = h
        p = self.p_hs(hs, s)
        ps = p
        v = self.v_ph(ps, hs)
        T = self.t_ph(ps, hs)
        p = self._unit_converter.toSIunit_p(p)
        T = self._unit_converter.toSIunit_T(T)
        v = self._unit_converter.toSIunit_v(v)
        rho = 1 / v
        return self._unit_converter.fromSIunit_tc(tc_ptrho(p, T, rho))

    def x_ph(self, p: float, h: float) -> float:
        """
        Section 1.17 Vapour fraction
        Vapour fraction as a function of pressure and enthalpy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)

        :return: vapor fraction or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        h = self._unit_converter.toSIunit_h(h)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            return self._unit_converter.fromSIunit_x(Region4.x4_ph(p, h))
        self.logger.debug("pressure out of range")
        return float("NaN")

    def x_ps(self, p: float, s: float) -> float:
        """
        Section 1.17 Vapour fraction
        Vapour fraction as a function of pressure and entropy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: vapor fraction or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        s = self._unit_converter.toSIunit_s(s)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            return self._unit_converter.fromSIunit_x(Region4.x4_ps(p, s))
        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def vx_ph(self, p: float, h: float) -> float:
        """
        Section 1.18 Vapour Volume Fraction
        Vapour volume fraction as a function of pressure and enthalpy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param h: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)

        :return: vapor volume fraction or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        h = self._unit_converter.toSIunit_h(h)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            if p < 16.529:
                vL = Region1.v1_pT(p, Region4.T4_p(p))
                vV = Region2.v2_pT(p, Region4.T4_p(p))
            else:
                vL = Region3.v3_ph(p, Region4.h4L_p(p))
                vV = Region3.v3_ph(p, Region4.h4V_p(p))
            xs = Region4.x4_ph(p, h)
            vx = xs * vV / (xs * vV + (1 - xs) * vL)
            return self._unit_converter.fromSIunit_vx(vx)

        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def vx_ps(self, p: float, s: float) -> float:
        """
        Section 1.18 Vapour Volume Fraction
        Vapour volume fraction as a function of pressure and entropy

        :param p: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        :param s: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)

        :return: vapor volume fraction or NaN if arguments are out of range
        """
        p = self._unit_converter.toSIunit_p(p)
        s = self._unit_converter.toSIunit_s(s)
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            if p < 16.529:
                vL = Region1.v1_pT(p, Region4.T4_p(p))
                vV = Region2.v2_pT(p, Region4.T4_p(p))
            else:
                vL = Region3.v3_ph(p, Region4.h4L_p(p))
                vV = Region3.v3_ph(p, Region4.h4V_p(p))
            xs = Region4.x4_ps(p, s)
            vx = xs * vV / (xs * vV + (1 - xs) * vL)
            return self._unit_converter.fromSIunit_vx(vx)

        self.logger.debug("pressure %f out of range", p)
        return float("NaN")

    def pmelt_t(self, t: float, hint: IceType = IceType.NONE) -> float:
        """R14-08(2011) Pressure along the melting curve as a function of temperature
        http://www.iapws.org/relguide/MeltSub2011.pdf

        Because of the shape of the meltin curve it is not possible to automatically
        select the correct region automatically. Therefore the optional
        hint-parameter is used to tell the function which area you are interested in.
        The hint-values are namend after the ice types.
        `IceType.Ih` = 1
        `IceType.III` = 3
        `IceType.V` = 5
        `IceType.VI` = 6
        `IceType.VII` = 7
        If the hint is not one of the values above or None(Default), an Exception is raised

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :param hint: optional :class:`IceType` selecting the melting-curve branch; use NONE for automatic selection where unambiguous

        :raises ValueError: unknown value for hint

        :return: pressure or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)

        if hint is None or hint == IceType.NONE:
            if 251.165 <= T < 256.164:
                self.logger.debug("can't select ice type based on temperature %f, hint required", T)
                return float("NaN")

            if 256.164 <= T < 273.31:
                self.logger.debug("chose ice type V based on temperature")
                return self._unit_converter.fromSIunit_p(pmelt_T_iceV(T))

            if 273.31 <= T < 355:
                self.logger.debug("chose ice type VI based on temperature")
                return self._unit_converter.fromSIunit_p(pmelt_T_iceVI(T))

            if 355 <= T < 751:
                self.logger.debug("chose ice type VII based on temperature")
                return self._unit_converter.fromSIunit_p(pmelt_T_iceVII(T))

            self.logger.debug("temperature %f out of range", T)
            return float("NaN")

        if hint == IceType.Ih:
            if 251.165 <= T < 273.16:
                return self._unit_converter.fromSIunit_p(pmelt_T_iceIh(T))
            self.logger.debug("temperature %f out of range", T)
            return float("NaN")

        if hint == IceType.III:
            if 251.165 <= T < 256.164:
                return self._unit_converter.fromSIunit_p(pmelt_T_iceIII(T))
            self.logger.debug("temperature %f out of range", T)
            return float("NaN")

        if hint == IceType.V:
            if 256.164 <= T < 273.31:
                return self._unit_converter.fromSIunit_p(pmelt_T_iceV(T))
            self.logger.debug("temperature %f out of range", T)
            return float("NaN")

        if hint == IceType.VI:
            if 273.31 <= T < 355:
                return self._unit_converter.fromSIunit_p(pmelt_T_iceVI(T))
            self.logger.debug("temperature %f out of range", T)
            return float("NaN")

        if hint == IceType.VII:
            if 355 <= T < 751:
                return self._unit_converter.fromSIunit_p(pmelt_T_iceVII(T))
            self.logger.debug("temperature %f out of range", T)
            return float("NaN")

        raise ValueError("unknown value for parameter 'hint'")

    def psubl_t(self, t: float) -> float:
        """R14-08(2011) Pressure along the sublimation curve as a function of temperature
        http://www.iapws.org/relguide/MeltSub2011.pdf

        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)

        :return: pressure or NaN if arguments are out of range
        """
        T = self._unit_converter.toSIunit_T(t)
        if 50 <= T < 273.16:
            return self._unit_converter.fromSIunit_p(psubl_T(T))
        self.logger.debug("temperature %f out of range", T)
        return float("NaN")

    def R12_my_rhot(self, rho: float, t: float, industrial_application: bool = True) -> float:
        """Return opt-in IAPWS R12-08 viscosity for density and temperature.

        Unlike the legacy ``my_pt``/``my_ph`` methods, this method uses the
        IAPWS 2008 viscosity formulation. It is opt-in to preserve existing
        viscosity results. ``rho`` and ``t`` use the selected unit system
        (kg/m³ and °C for MKS; lb/ft³ and °F for FLS); the returned viscosity
        follows the selected system (Pa·s for MKS; lbm/(ft·hr) for FLS).

        :param rho: density in the selected system (kg/m³ for MKS/BARE; lb/ft³ for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :param industrial_application: omit critical enhancement when True
        :return: dynamic viscosity in the selected unit system
        """
        self.logger.warning("this function is still experimental, use at your own risk!")

        density_si = self._unit_converter.toSIunit_rho(rho)
        temperature_si = self._unit_converter.toSIunit_T(t)
        viscosity_si = my_rhoT(density_si, temperature_si, industrial_application)
        return self._unit_converter.fromSIunit_my(viscosity_si)

    def p_rhoT(self, rho: float, t: float) -> float:
        """Calculate IAPWS-95 pressure from density and temperature.

        :param rho: density in the selected system (kg/m³ for MKS/BARE; lb/ft³ for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :return: pressure in the selected system (bar for MKS/BARE; psi for FLS)
        """
        return self._unit_converter.fromSIunit_p(
            _R6_p_rhoT(
                self._unit_converter.toSIunit_rho(rho),
                self._unit_converter.toSIunit_T(t),
            )
        )

    def u_rhoT(self, rho: float, t: float) -> float:
        """Calculate IAPWS-95 internal energy from density and temperature.

        :param rho: density in the selected system (kg/m³ for MKS/BARE; lb/ft³ for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :return: specific internal energy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)
        """
        return self._unit_converter.fromSIunit_u(
            _R6_u_rhoT(
                self._unit_converter.toSIunit_rho(rho),
                self._unit_converter.toSIunit_T(t),
            )
        )

    def s_rhoT(self, rho: float, t: float) -> float:
        """Calculate IAPWS-95 entropy from density and temperature.

        :param rho: density in the selected system (kg/m³ for MKS/BARE; lb/ft³ for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :return: specific entropy in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)
        """
        return self._unit_converter.fromSIunit_s(
            _R6_s_rhoT(
                self._unit_converter.toSIunit_rho(rho),
                self._unit_converter.toSIunit_T(t),
            )
        )

    def h_rhoT(self, rho: float, t: float) -> float:
        """Calculate IAPWS-95 enthalpy from density and temperature.

        :param rho: density in the selected system (kg/m³ for MKS/BARE; lb/ft³ for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :return: specific enthalpy in the selected system (kJ/kg for MKS/BARE; Btu/lb for FLS)
        """
        return self._unit_converter.fromSIunit_h(
            _R6_h_rhoT(
                self._unit_converter.toSIunit_rho(rho),
                self._unit_converter.toSIunit_T(t),
            )
        )

    def cv_rhoT(self, rho: float, t: float) -> float:
        """Calculate IAPWS-95 isochoric heat capacity from density and temperature.

        :param rho: density in the selected system (kg/m³ for MKS/BARE; lb/ft³ for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :return: specific heat capacity in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)
        """
        return self._unit_converter.fromSIunit_Cv(
            _R6_cv_rhoT(
                self._unit_converter.toSIunit_rho(rho),
                self._unit_converter.toSIunit_T(t),
            )
        )

    def cp_rhoT(self, rho: float, t: float) -> float:
        """Calculate IAPWS-95 isobaric heat capacity from density and temperature.

        :param rho: density in the selected system (kg/m³ for MKS/BARE; lb/ft³ for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :return: specific heat capacity in the selected system (kJ/(kg·K) for MKS/BARE; Btu/(lb·°F) for FLS)
        """
        return self._unit_converter.fromSIunit_Cp(
            _R6_cp_rhoT(
                self._unit_converter.toSIunit_rho(rho),
                self._unit_converter.toSIunit_T(t),
            )
        )

    def w_rhoT(self, rho: float, t: float) -> float:
        """Calculate IAPWS-95 speed of sound from density and temperature.

        :param rho: density in the selected system (kg/m³ for MKS/BARE; lb/ft³ for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :return: speed in the selected system (m/s for MKS/BARE; ft/s for FLS)
        """
        return self._unit_converter.fromSIunit_w(
            _R6_w_rhoT(
                self._unit_converter.toSIunit_rho(rho),
                self._unit_converter.toSIunit_T(t),
            )
        )

    def joulethomson_rhoT(self, rho: float, t: float) -> float:
        """Calculate IAPWS-95 Joule-Thomson coefficient from density and temperature.

        :param rho: density in the selected system (kg/m³ for MKS/BARE; lb/ft³ for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :return: Joule-Thomson temperature-pressure coefficient in selected temperature and pressure units
        """
        return self._unit_converter.fromSIunit_pressure_derivative(
            _R6_joulethomson_rhoT(
                self._unit_converter.toSIunit_rho(rho),
                self._unit_converter.toSIunit_T(t),
            )
        )

    def delta_T_rhoT(self, rho: float, t: float) -> float:
        """Calculate IAPWS-95 isothermal throttling coefficient from density and temperature.

        :param rho: density in the selected system (kg/m³ for MKS/BARE; lb/ft³ for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :return: specific-volume units in the selected system (m³/kg for MKS/BARE; ft³/lb for FLS)
        """
        return self._unit_converter.fromSIunit_v(
            _R6_delta_T_rhoT(
                self._unit_converter.toSIunit_rho(rho),
                self._unit_converter.toSIunit_T(t),
            )
        )

    def isen_t_p_rhoT(self, rho: float, t: float) -> float:
        """Calculate the IAPWS-95 isentropic temperature-pressure coefficient.

        :param rho: density in the selected system (kg/m³ for MKS/BARE; lb/ft³ for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :return: coefficient in selected temperature units per selected pressure unit
        """
        return self._unit_converter.fromSIunit_pressure_derivative(
            _R6_isen_t_p_rhoT(
                self._unit_converter.toSIunit_rho(rho),
                self._unit_converter.toSIunit_T(t),
            )
        )

    def sec_virial_coefficients(self, rho: float, t: float) -> tuple[float, float]:
        """Calculate IAPWS-95 second virial coefficient B and dB/dT.

        :param rho: density in the selected system; retained for API consistency
            but not used because virial coefficients are zero-density limits
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :return: (B, dB/dT), in selected specific-volume and temperature units
        """
        coefficient, derivative = _R6_sec_virial_coefficients(
            self._unit_converter.toSIunit_rho(rho),
            self._unit_converter.toSIunit_T(t),
        )
        return (
            self._unit_converter.fromSIunit_v(coefficient),
            self._unit_converter.fromSIunit_v(derivative) * self._unit_converter.toSIunit_T_interval(1),
        )

    def third_virial_coefficients(self, rho: float, t: float) -> tuple[float, float]:
        """Calculate IAPWS-95 third virial coefficient C and dC/dT.

        :param rho: density in the selected system; retained for API consistency
            but not used because virial coefficients are zero-density limits
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :return: (C, dC/dT), in selected squared-specific-volume and temperature units
        """
        coefficient, derivative = _R6_third_virial_coefficients(
            self._unit_converter.toSIunit_rho(rho),
            self._unit_converter.toSIunit_T(t),
        )
        return (
            self._unit_converter.fromSIunit_v_squared(coefficient),
            self._unit_converter.fromSIunit_v_squared(derivative) * self._unit_converter.toSIunit_T_interval(1),
        )

    def phase_equilibrium_condition(self, rho_liquid: float, rho_vapor: float, t: float) -> tuple[float, float]:
        """Return IAPWS-95 pressure-equality and Maxwell residuals.

        :param rho_liquid: liquid density in the selected system (kg/m³ for
            MKS/BARE; lb/ft³ for FLS)
        :param rho_vapor: vapor density in the selected system (kg/m³ for
            MKS/BARE; lb/ft³ for FLS)
        :param t: temperature in the selected system (°C for MKS/BARE; °F for FLS)
        :return: pressure residual (bar for MKS/BARE; psi for FLS) and
            dimensionless Maxwell residual
        """
        pressure_residual, maxwell_residual = _R6_phase_equilibrium_condition(
            self._unit_converter.toSIunit_rho(rho_liquid),
            self._unit_converter.toSIunit_rho(rho_vapor),
            self._unit_converter.toSIunit_T(t),
        )
        return (
            self._unit_converter.fromSIunit_p(pressure_residual),
            maxwell_residual,
        )
