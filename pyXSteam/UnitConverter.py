#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class to convert between the unit system used by pyXSteam and the ones
a enduser might use.
"""
import logging
from .Constants import ABSOLUTE_ZERO_CELSIUS, UnitSystem


class UnitConverter(object):
    """
    Helper class to convert user units to SI-units and back

    :param unitSystem: unit system used for input and output values. For supported values see the enum UnitSystem.
    """

    def __init__(self, unitSystem: UnitSystem = UnitSystem.BARE):
        """
        Constructor method
        """
        self.logger = logging.getLogger(__name__)
        self.set_unitSystem(unitSystem)
        self.logger.debug("set unit converter to %s", self.__str__())

    def set_unitSystem(self, unitSystem: UnitSystem):
        """
        change unit system

        :param unitSystem: new unit system to use for input and output values

        :raises ValueError: unknown value for unit system
        """
        if UnitSystem.has_value(unitSystem):
            self._unit_system = unitSystem
        else:
            self.logger.error("Unknown Unit System selected")
            raise ValueError("Unknown Unit System")

    def toSIunit_p(self, ins: float) -> float:
        """
        convert preasure from user selected unit system to SI units

        :param ins: preasure in [bar] or [psi]
        """
        if self._unit_system is UnitSystem.MKS:
            return float(ins / 10)  # bar to MPa
        elif self._unit_system is UnitSystem.FLS:
            return float(ins * 0.00689475729)  # psi to MPa
        return float(ins)

    def fromSIunit_p(self, ins: float) -> float:
        """
        convert preasure from SI units to user selected unit system

        :param ins: preasure in [MPa]
        """
        if self._unit_system is UnitSystem.MKS:
            return float(ins * 10)  # bar to MPa
        elif self._unit_system is UnitSystem.FLS:
            return float(ins / 0.00689475729)  # MPa to psi
        return float(ins)

    def toSIunit_T(self, ins: float) -> float:
        """
        convert temperature from user selected unit system to SI units

        :param ins: temperature in [°C] or [°F]
        """
        if self._unit_system is UnitSystem.MKS:
            # degC to Kelvin
            return float(ins - ABSOLUTE_ZERO_CELSIUS)
        elif self._unit_system is UnitSystem.FLS:
            return float((5 / 9) * (ins - 32) - ABSOLUTE_ZERO_CELSIUS)  # degF to Kelvin
        return float(ins)

    def fromSIunit_T(self, ins: float) -> float:
        """
        convert temperature from SI units to user selected unit system

        :param ins: temperature in [K]
        """
        if self._unit_system is UnitSystem.MKS:
            # Kelvin to degC
            return float(ins + ABSOLUTE_ZERO_CELSIUS)
        elif self._unit_system is UnitSystem.FLS:
            return float((ins + ABSOLUTE_ZERO_CELSIUS) * (9 / 5) + 32)  # Kelvin to degF
        return float(ins)

    def toSIunit_h(self, ins: float) -> float:
        """
        convert enthalpy from user selected unit system to SI units

        :param ins: enthalpy [kJ / kg] or [btu / lb]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(2.32600 * ins)  # btu/lb to kJ/kg
        return float(ins)

    def fromSIunit_h(self, ins: float) -> float:
        """
        convert enthalpy from SI units to user selected unit system

        :param ins: enthalpy [kJ / kg]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 2.32600)  # kJ/kg to btu/lb
        return float(ins)

    def toSIunit_v(self, ins: float) -> float:
        """
        convert specific volume from user selected unit system to SI units

        :param ins: specific volume in [m³ / kg] or [ft³ / lb]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 0.0624279606)  # ft³/lb to m³/kg
        return float(ins)

    def fromSIunit_v(self, ins: float) -> float:
        """
        convert specific volume from SI units to user selected unit system

        :param ins: specific volume in [m³ / kg]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 0.0624279606)  # m³/kg to ft³/lb
        return float(ins)

    def toSIunit_s(self, ins: float) -> float:
        """
        convert specific entropy from user selected unit system to SI units

        :param ins: specific volume in [kJ / (kg °C)] or [btu / (lb °F)]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 0.238845896627)  # btu/(lb degF) to kJ/(kg degC)
        return float(ins)

    def fromSIunit_s(self, ins: float) -> float:
        """
        convert specific entropy from SI units to user selected unit system

        :param ins: specific entropy in [kJ / (kg °C)]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 0.238845896627)  # kJ/(kg degC) to btu/(lb degF)
        return float(ins)

    def toSIunit_u(self, ins: float) -> float:
        """
        convert specific internal energy from user selected unit system to SI units

        :param ins: specific internal energy in [kJ / kg] or [btu / lb]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 2.32600)  # btu/lb to kJ/kg
        return float(ins)

    def fromSIunit_u(self, ins: float) -> float:
        """
        convert specific internal energy from SI units to user selected unit system

        :param ins: specific internal energy in [kJ / kg]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 2.32600)  # kJ/kg to btu/lb
        return float(ins)

    def toSIunit_Cp(self, ins: float) -> float:
        """
        convert specific isobaric heat capacity from user selected unit system to SI units

        :param ins: specific isobaric heat capacity in [kJ / (kg °C)] or [btu / (lb °F)]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 0.238846)  # btu/(lb degF) to kJ/(kg degC)
        return float(ins)

    def fromSIunit_Cp(self, ins: float) -> float:
        """
        convert specific isobaric heat capacity from SI units to user selected unit system

        :param ins: specific isobaric heat capacity in [kJ / (kg °C)]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 0.238846)  # kJ/(kg degC) to btu/(lb degF)
        return float(ins)

    def toSIunit_Cv(self, ins: float) -> float:
        """
        convert specific isochoric heat capacity from user selected unit system to SI units

        :param ins: specific isochoric heat capacity in [kJ / (kg °C)] or [btu / (lb °F)]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 0.238846)  # btu/(lb degF) to kJ/(kg degC)
        return float(ins)

    def fromSIunit_Cv(self, ins: float) -> float:
        """
        convert specific isochoric heat capacity from SI units to user selected unit system

        :param ins: specific isochoric heat capacity in [kJ / (kg °C)]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 0.238846)  # kJ/(kg degC) to btu/(lb degF)
        return float(ins)

    def toSIunit_w(self, ins: float) -> float:
        """
        convert speed of sound from user selected unit system to SI units

        :param ins: speed of sound in [m / s] or [ft / s]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 0.3048)  # ft/s to m/s
        return float(ins)

    def fromSIunit_w(self, ins: float) -> float:
        """
        convert speed of sound from SI units to user selected unit system

        :param ins: speed of sound in [m / s]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 0.3048)  # m/s to ft/s
        return float(ins)

    def toSIunit_tc(self, ins: float) -> float:
        """
        convert thermal conductivity from user selected unit system to SI units

        :param ins: thermal conductivity in [W / (m °C)] or [btu / (h ft °F)]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 0.577789)  # btu/(h*ft*degF) to W/(m*degC)
        return float(ins)

    def fromSIunit_tc(self, ins: float) -> float:
        """
        convert thermal conductivity from SI units to user selected unit system

        :param ins: thermal conductivity in [W / (m °C)]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 0.577789)  # W/(m*degC) to btu/(h*ft*degF)
        return float(ins)

    def toSIunit_st(self, ins: float) -> float:
        """
        convert surface tension from user selected unit system to SI units

        :param ins: surface tension in [N / m] or [lb / ft]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 0.068521766)  # lb/ft to N/m
        return float(ins)

    def fromSIunit_st(self, ins: float) -> float:
        """
        convert surface tension from SI units to user selected unit system

        :param ins: surface tension in [N / m]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 0.068521766)  # N/m to lb/ft
        return float(ins)

    def toSIunit_x(self, ins: float) -> float:
        """
        convert vapor fraction from user selected unit system to SI units

        :param ins: vapor fraction

        :raises ValueError: value of vapour fraction out of range
        """
        # TODO: Check if <= should be <
        if 0.0 <= ins <= 1.0:
            return float(ins)
        self.logger.error("value of vapour fraction out of range: 0 < x < 1")
        raise ValueError("Vapour fraction out of Range")

    def fromSIunit_x(self, ins: float) -> float:
        """
        convert vapor fraction from SI units to user selected unit system

        :param ins: vapor fraction

        :raises ValueError: value of vapour fraction out of range
        """
        # TODO: Check if <= should be <
        if 0.0 <= ins <= 1.0:
            return float(ins)
        self.logger.error("value of vapour fraction out of range: 0 < x < 1")
        raise ValueError("Vapour fraction out of Range")

    def toSIunit_vx(self, ins: float) -> float:
        """
        convert vapor volume fraction from user selected unit system to SI units

        :param ins: vapor volume fraction

        :raises ValueError: value of vapour volume fraction out of range
        """
        # TODO: Check if <= should be <
        if 0.0 <= ins <= 1.0:
            return float(ins)
        self.logger.error("value of vapour volume fraction out of range: 0 < x < 1")
        raise ValueError("Vapour volume fraction out of Range")

    def fromSIunit_vx(self, ins: float) -> float:
        """
        convert vapor volume fraction from SI units to user selected unit system

        :param ins: vapor volume fraction

        :raises ValueError: value of vapour volume fraction out of range
        """
        # TODO: Check if <= should be <
        if 0.0 <= ins <= 1.0:
            return float(ins)
        self.logger.error("value of vapour volume fraction out of range: 0 < x < 1")
        raise ValueError("Vapour volume fraction out of Range")

    def toSIunit_my(self, ins: float) -> float:
        """
        convert viscosity from user selected unit system to SI units

        :param ins: viscosity in [Pa s], [N s / m²] or [lbm / ft / hr]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 2419.088311)  # lbm/ft/hr to PaS (N*s/m²)
        return float(ins)

    def fromSIunit_my(self, ins: float) -> float:
        """
        convert viscosity from SI units to user selected unit system

        :param ins: viscosity in [PaS] or [N*s/m²]
        """
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 2419.088311)  # PaS (N*s/m²) to lbm/ft/hr
        return float(ins)

    def __str__(self):
        """
        :return: string representation of the selected unit system
        """
        result = ""
        if self._unit_system is UnitSystem.FLS:
            result = "FLS (ft/lb/sec/°F/psi/btu)"
        elif self._unit_system is UnitSystem.MKS:
            result = "MKS (m/kg/sec/°C/bar/W)"
        else:
            result = "BARE (m/kg/sec/K/MPa/W)"
        return result
