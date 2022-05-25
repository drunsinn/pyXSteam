#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Constants for the calculation of water steam properties

Sources:

* IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997

* IAWPS Release on Vaues of Temperature, Pressure and Density of Ordinary and Heavy Water Substances at their Respective Critical Points Released September 1992, Revision of the Release of 1992

"""
from enum import IntEnum

SPECIFIC_GAS_CONSTANT = 0.461526  # kJ / (kg K)
CRITICAL_TEMPERATURE = 647.096  # K
CRITICAL_PRESSURE = 22.06395  # MPa
CRITICAL_DENSITY = 322  # kg / m
TRIPLE_POINT_TEMPERATURE = 273.16  # K (Eq9 Page 7)
TRIPLE_POINT_PRESSURE = 0.000611657  # MPa (Eq9 Page 7)
TRIPLE_POINT_SPECIFIC_ENTHALPY = 0.611783e-3  # kJ / kg (Eq10 Page 7)
FREEZING_TEMPERATURE_H2O = 273.15  # K


# IAWPS Release on Vaues of Temperature, Pressure and Density of Ordinary and
# Heavy Water Substances at their Respective Critical Points
# Released September 1992, Revision of the Release of 1992
CRITICAL_TEMPERATURE_H20_1992 = 647.096  # +-0.1 K
CRITICAL_PRESSURE_H20_1992 = 22.067  # +0.27*(+-0.1)+-0.005 MPa
CRITICAL_DENSITY_H20_1992 = 322  # +-3 kg / m

CRITICAL_TEMPERATURE_D20_1992 = 643.847  # +-0.2 K
CRITICAL_PRESSURE_D20_1992 = 21.671  # +0.27*(+-0.2)+-0.01 MPa
CRITICAL_DENSITY_D20_1992 = 356  # +-5 kg / m


# Other common constants used in calculations
ABSOLUTE_ZERO_CELSIUS = -273.15  # °C
ABSOLUTE_ZERO_FAHRENHEIT = -459.67  # °F


# IAPWS R15-11
__SPECIFIC_GAS_CONSTANT_IAPWS_R15_11__ = 0.46151805  # kJ kg^-1 K^-1
__CRITICAL_TEMPERATURE_IAPWS_R15_11__ = 647.096  # K
__CRITICAL_DENSITY_IAPWS_R15_11__ = 322.0  # kg m^-1


class UnitSystem(IntEnum):
    BARE = 1  # m/kg/sec/K/MPa/W
    MKS = 1  # m/kg/sec/°C/bar/W
    FLS = 2  # ft/lb/sec/°F/psi/btu

    @classmethod
    def has_value(cls, value):
        return value in cls._value2member_map_


class IceType(IntEnum):
    Ih = 1
    III = 3
    V = 5
    VI = 6
    VII = 7

    @classmethod
    def has_value(cls, value):
        return value in cls._value2member_map_


class DiagramRegion(IntEnum):
    NILL = 0  # Error, Outside valid area
    R1 = 1
    R2 = 2
    R3 = 3
    R4 = 4
    R5 = 5
