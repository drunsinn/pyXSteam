#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import logging
from . import Constants
from .UnitConverter import UnitConverter
from . import IAPWS_R4

class XSteam_HW(object):
    """Main pyXSteam for Heavy Water object. Abstract of all other functions to allow auto selection of
    the correct region for each set of parameters.

    Args:
        unitSystem (int): set the unit system used for input and output values.
            Can be eather 0 (UNIT_SYSTEM_BARE), 1 (UNIT_SYSTEM_MKS) or 2 (UNIT_SYSTEM_FLS).
    """

    # Copy constant Values to expose them to the User
    UNIT_SYSTEM_BARE = UnitConverter.__UNIT_SYSTEM_BARE__
    UNIT_SYSTEM_MKS = UnitConverter.__UNIT_SYSTEM_MKS__
    UNIT_SYSTEM_FLS = UnitConverter.__UNIT_SYSTEM_FLS__

    def __init__(self, unitSystem=UnitConverter.__UNIT_SYSTEM_BARE__):
        self.logger = logging.getLogger(__name__)
        self.unitConverter = UnitConverter(unitSystem)
        self.logger.info('initialised pyXSteam for Heavy Water with Unit System "{}"'.format(self.unitConverter))

    def criticalTemperatur(self):
        """returns the specific temperature with conversion to the selected unit system"""
        return self.unitConverter.fromSIunit_T(Constants.__CRITICAL_TEMPERATURE_D20_1992__)

    def criticalPressure(self):
        """returns the specific pressure with conversion to the selected unit system"""
        return self.unitConverter.fromSIunit_p(Constants.__CRITICAL_PRESSURE_D20_1992__)

    def criticalDensity(self):
        """returns the specific density with conversion to the selected unit system"""
        return self.unitConverter.fromSIunit_p(Constants.__CRITICAL_DENSITY_D20_1992__)

    def my_rhoT(self, rho, T):
        """Viscosity as a function of density and temperature for heavy water
        substance

        Source: IAPWS R4-84(2007)
        Revised Release on Viscosity and Thermal Conductivity of Heavy Water Substance
        http://www.iapws.org/relguide/TransD2O-2007.pdf
        Appendix A

        Args:
            rho (float): density value for heavy water
            T (float): temperature value for heavy water

        Returns:
            my (float): viscosity or NaN if arguments are out of range
        """
        rho = self.unitConverter.toSIunit_p(rho)
        T = self.unitConverter.toSIunit_T(T)

        if T < 277.0 or T > 775.0:
            self.logger.error('temperature out of range')
            return float("NaN")

        self.logger.warning('input for desity wasnt checked!')

        return self.unitConverter.fromSIunit_T(IAPWS_R4.myHW_rhoT_R4(rho, T))

    def tc_rhoT(self, rho, T):
        """Thermal conductivity as a function of density and temperature for heavy water
        substance

        Source: IAPWS R4-84(2007)
        Revised Release on Viscosity and Thermal Conductivity of Heavy Water Substance
        http://www.iapws.org/relguide/TransD2O-2007.pdf
        Appendix B

        Args:
            rho (float): density value for heavy water
            T (float): temperature value for heavy water

        Returns:
            λ (float): thermal conductivity or NaN if arguments are out of range
        """
        rho = self.unitConverter.toSIunit_p(rho)
        T = self.unitConverter.toSIunit_T(T)

        if T < 277.0 or T > 825.0:
            self.logger.error('temperature out of range')
            return float("NaN")

        self.logger.warning('input for desity wasnt checked!')

        return self.unitConverter.fromSIunit_tc(IAPWS_R4.tcHW_rhoT_R4(rho, T))
