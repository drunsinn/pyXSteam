#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Constants for the calculation of water steam properties

Sources:

* IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997

* IAWPS Releas on Vaues of Temperature, Pressure and Density of Ordinary and Heavy Water Substances at their Respective Critical Points Released September 1992, Revision of the Release of 1992

"""

__SPECIFIC_GAS_CONSTANT__ = 0.461526  # kJ kg^-1 K^-1
__CRITICAL_TEMPERATURE__ = 647.096  # K
__CRITICAL_PRESSURE__ = 22.06395  # MPa
__CRITICAL_DENSITY__ = 322  # kg m^-1
__TRIPLE_POINT_TEMPERATURE__ = 273.16  # K (Eq9 Page 7)
__TRIPLE_POINT_PRESSURE__ = 0.000611657  # MPa (Eq9 Page 7)
__TRIPLE_POINT_SPECIFIC_ENTHALPY__ = 0.611783E-3  # kJ kg^-1 (Eq10 Page 7)


# IAWPS Releas on Vaues of Temperature, Pressure and Density of Ordinary and
# Heavy Water Substances at their Respective Critical Points
# Released September 1992, Revision of the Release of 1992
__CRITICAL_TEMPERATURE_H20_1992__ = 647.096  # +-0.1 K
__CRITICAL_PRESSURE_H20_1992__ = 22.067  # + 0.27*(+-0.1)+-0.005 MPa
__CRITICAL_DENSITY_H20_1992__ = 322  # +-3 kg m^-1

__CRITICAL_TEMPERATURE_D20_1992__ = 643.847  # +-0.2 K
__CRITICAL_PRESSURE_D20_1992__ = 21.671  # + 0.27*(+-0.2)+-0.01 MPa
__CRITICAL_DENSITY_D20_1992__ = 356  # +-5 kg m^-1


# Other common constants used in calculations
__ABSOLUTE_ZERO_CELSIUS__ = -273.15  # °C
__ABSOLUTE_ZERO_FAHRENHEIT__ = -459.67  # °F


# IAPWS R15-11
__SPECIFIC_GAS_CONSTANT_IAPWS_R15_11__ = 0.46151805 # kJ kg^-1 K^-1
__CRITICAL_TEMPERATURE_IAPWS_R15_11__ = 647.096  # K
__CRITICAL_DENSITY_IAPWS_R15_11__ = 322.0  # kg m^-1
