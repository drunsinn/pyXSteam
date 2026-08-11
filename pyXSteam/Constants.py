#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Constants for the calculation of water steam properties

Sources:
* IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Waterand Steam, September 1997
* IAWPS R2-83(1992) Values of Temperature, Pressure and Density of Ordinary and Heavy Water Substances at their Respective Critical Points
* IAPWS R15-11
* IAPWS R17 and IAPWS R18
* https://doi.org/10.1063/1.5053993
"""

from enum import IntEnum

SPECIFIC_GAS_CONSTANT = 0.461526  # [kJ / (kg K)]
CRITICAL_TEMPERATURE = 647.096  # [K]
CRITICAL_PRESSURE = 22.06395  # [MPa]
CRITICAL_DENSITY = 322.0  # [kg / m³]
TRIPLE_POINT_TEMPERATURE = 273.16  # [K] (Eq9 Page 7)
TRIPLE_POINT_PRESSURE = 0.000611657  # [MPa] (Eq9 Page 7)
TRIPLE_POINT_SPECIFIC_ENTHALPY = 0.611783e-3  # [kJ / kg] (Eq10 Page 7)
FREEZING_TEMPERATURE_H2O = 273.15  # [K]


# IAWPS R2-83(1992) Release on Values of Temperature, Pressure and Density of Ordinary and Heavy Water Substances at their
# Respective Critical Points, Released September 1992, Revision of the Release of 1983
# http://www.iapws.org/relguide/crits.pdf
CRITICAL_TEMPERATURE_H20_1992 = 647.096  # [K] ±0.1
CRITICAL_PRESSURE_H20_1992 = 22.067  # [MPa] ±0.005
CRITICAL_DENSITY_H20_1992 = 322.0  # [kg / m³] ±3

CRITICAL_TEMPERATURE_D20_1992 = 643.847  # [K] ±0.2
CRITICAL_PRESSURE_D20_1992 = 21.671  # [MPa] ±0.01
CRITICAL_DENSITY_D20_1992 = 356.0  # [kg / m³] ±5

# Other common constants used in calculations
ABSOLUTE_ZERO_CELSIUS = -273.15  # [°C]
ABSOLUTE_ZERO_FAHRENHEIT = -459.67  # [°F]

# IAPWS R15-11
__SPECIFIC_GAS_CONSTANT_IAPWS_R15_11__ = 0.46151805  # [kJ kg^-1 K^-1]
__CRITICAL_TEMPERATURE_IAPWS_R15_11__ = 647.096  # [K]
__CRITICAL_DENSITY_IAPWS_R15_11__ = 322.0  # [kg / m³]

# IAPWS R17 and IAPWS R18
__REFERENCE_TEMPERATURE_D20_R17_R18__ = 643.847  # T* in [K]
__REFERENCE_PREASSURE_D20_R17_R18__ = 21.6618  # p* in [MPa]
__REFERENCE_DENSITY_D20_R17_R18__ = 356.0  # ρ* in [kg / m³]
__REFERENCE_VISCOSITY_D20_R17_R18__ = 1.00e-6  # μ* in [Pa s]
__REFERENCE_THERMAL_CONDUCTIVITY_D20_R18__ = 1e-3  # λ* in [W m^-1 K^-1]
__SPECIFIC_GAS_CONSTANT_D20_R17_R18__ = 0.41515199  # R in [kJ kg^-1 K^-1]


# https://doi.org/10.1063/1.5053993
__TRIPLE_POINT_TEMPERATURE_D2O_RESHW_2018__ = 276.969  # [K]
__TRIPLE_POINT_PRESSURE_D20_RESHW_2018__ = 0.00061159  # [MPa]


class UnitSystem(IntEnum):
    """enum for supported unit systems"""

    BARE = 1  # [m/kg/sec/K/MPa/W]
    MKS = 1  # [m/kg/sec/°C/bar/W]
    FLS = 2  # [ft/lb/sec/°F/psi/btu]

    @classmethod
    def has_value(cls, value: int):
        """check if value is member of enum"""
        return value in cls._value2member_map_


class IceType(IntEnum):
    """enum for the types of ice"""

    Ih = 1
    III = 3
    V = 5
    VI = 6
    VII = 7
    NONE = -1

    @classmethod
    def has_value(cls, value: int):
        """check if value is member of enum"""
        return value in cls._value2member_map_


class DiagramRegion(IntEnum):
    """enum for the regions"""

    NILL = 0  # Error, Outside valid area
    R1 = 1
    R2 = 2
    R3 = 3
    R4 = 4
    R5 = 5
