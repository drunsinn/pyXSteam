#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Main module for the heavy water parts of pyXSteam"""
import logging
from .Constants import (
    CRITICAL_TEMPERATURE_D20_1992,
    CRITICAL_PRESSURE_D20_1992,
    CRITICAL_DENSITY_D20_1992,
    UnitSystem,
)
from .UnitConverter import UnitConverter
from .IAPWS_R4 import myHW_rhoT_R4, tcHW_rhoT_R4
from .IAPWS_R5 import surface_tension_T


class XSteam_HW:
    """Main pyXSteam for Heavy Water object. Abstract of all other functions to allow auto selection of
    the correct region for each set of parameters.

    :param unitSystem: unit system used for input and output values. For supported values
        see the enum UnitSystem.
    """

    UNIT_SYSTEM_BARE = UnitSystem.BARE
    UNIT_SYSTEM_MKS = UnitSystem.MKS
    UNIT_SYSTEM_FLS = UnitSystem.FLS

    def __init__(self, unitSystem: UnitSystem = UnitSystem.BARE):
        """
        Constructor method
        """
        self.logger = logging.getLogger("pyXSteam:Heavy Water")
        self._unit_converter = UnitConverter(unitSystem)
        self.logger.info(
            "initialised pyXSteam for Heavy Water with Unit System %s",
            self._unit_converter,
        )

    def criticalTemperatur(self):
        """
        :return: specific temperature with conversion to the selected unit system
        """
        return self._unit_converter.fromSIunit_T(CRITICAL_TEMPERATURE_D20_1992)

    def criticalPressure(self):
        """
        :return: specific pressure with conversion to the selected unit system
        """
        return self._unit_converter.fromSIunit_p(CRITICAL_PRESSURE_D20_1992)

    def criticalDensity(self):
        """
        :return: specific density with conversion to the selected unit system
        """
        return self._unit_converter.fromSIunit_p(CRITICAL_DENSITY_D20_1992)

    def my_rhoT(self, rho: float, T: float) -> float:
        """Viscosity as a function of density and temperature for heavy water
        substance

        Source: IAPWS R4-84(2007)
        Revised Release on Viscosity and Thermal Conductivity of Heavy Water Substance
        http://www.iapws.org/relguide/TransD2O-2007.pdf
        Appendix A

        :param rho: density for heavy water
        :param T: temperature for heavy water

        :raises ValueError: value of density is zero or negative

        :return: viscosity or NaN if arguments are out of range
        """
        if rho <= 0.0:
            self.logger.error("negative values for density rho not allowed %f", rho)
            raise ValueError("rho out of range")

        rho = self._unit_converter.toSIunit_p(rho)
        T = self._unit_converter.toSIunit_T(T)

        if T < 277.0 or T > 775.0:
            self.logger.error("temperature out of range")
            return float("NaN")

        return self._unit_converter.fromSIunit_T(myHW_rhoT_R4(rho, T))

    def tc_rhoT(self, rho: float, T: float) -> float:
        """Thermal conductivity as a function of density and temperature for heavy water
        substance

        Source: IAPWS R4-84(2007)
        Revised Release on Viscosity and Thermal Conductivity of Heavy Water Substance
        http://www.iapws.org/relguide/TransD2O-2007.pdf
        Appendix B

        :param rho: density for heavy water
        :param T: temperature for heavy water

        :raises ValueError: value of density is zero or negative

        :return: thermal conductivity or NaN if arguments are out of range
        """
        if rho <= 0.0:
            self.logger.error("negative values for density rho not allowed %f", rho)
            raise ValueError("rho out of range")

        rho = self._unit_converter.toSIunit_p(rho)
        T = self._unit_converter.toSIunit_T(T)

        if T < 277.0 or T > 825.0:
            self.logger.error("temperature out of range")
            return float("NaN")

        return self._unit_converter.fromSIunit_tc(tcHW_rhoT_R4(rho, T))

    def st_t(self, t: float) -> float:
        """
        Surface tension for two phase water/steam as a function of temperature

        :param t: temperature

        :return: surface tension
        """
        T = self._unit_converter.toSIunit_T(t)
        return self._unit_converter.fromSIunit_st(surface_tension_T(T))
