#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class to convert between the unit system used by pyXSteam and the ones
a enduser might use.
"""
import logging
from .Constants import __ABSOLUTE_ZERO_CELSIUS__, UnitSystem


class UnitConverter(object):
    """Helper class to convert user units to SI-units and back"""

    def __init__(self, unitSystem: UnitSystem = UnitSystem.BARE):
        """Initialise the unit converter. Parameter is the unit system used by the application"""
        self.logger = logging.getLogger("pyXSteam-UnitConverter")
        self.set_unitSystem(unitSystem)
        self.logger.debug("set unit converter to %s", self.__str__())

    def set_unitSystem(self, unitSystem):
        """change unit system"""
        if UnitSystem.has_value(unitSystem):
            self._unit_system = unitSystem
        else:
            self.logger.error("Unknown Unit System selected")
            raise ValueError("Unknown Unit System")

    def toSIunit_p(self, ins: float) -> float:
        """function toSIunit_p = toSIunit_p( ins )"""
        if self._unit_system is UnitSystem.MKS:
            return float(ins / 10)  # bar to MPa
        elif self._unit_system is UnitSystem.FLS:
            return float(ins * 0.00689475729)  # psi to MPa
        return float(ins)

    def fromSIunit_p(self, ins: float) -> float:
        """function fromSIunit_p = fromSIunit_p( ins )"""
        if self._unit_system is UnitSystem.MKS:
            return float(ins * 10)  # bar to MPa
        elif self._unit_system is UnitSystem.FLS:
            return float(ins / 0.00689475729)  # MPa to psi
        return float(ins)

    def toSIunit_T(self, ins: float) -> float:
        """function toSIunit_T = toSIunit_T( ins )"""
        if self._unit_system is UnitSystem.MKS:
            # degC to Kelvin
            return float(ins - __ABSOLUTE_ZERO_CELSIUS__)
        elif self._unit_system is UnitSystem.FLS:
            return float(
                (5 / 9) * (ins - 32) - __ABSOLUTE_ZERO_CELSIUS__
            )  # degF to Kelvin
        return float(ins)

    def fromSIunit_T(self, ins: float) -> float:
        """function fromSIunit_T = fromSIunit_T( ins )"""
        if self._unit_system is UnitSystem.MKS:
            # Kelvin to degC
            return float(ins + __ABSOLUTE_ZERO_CELSIUS__)
        elif self._unit_system is UnitSystem.FLS:
            return float(
                (ins + __ABSOLUTE_ZERO_CELSIUS__) * (9 / 5) + 32
            )  # Kelvin to degF
        return float(ins)

    def toSIunit_h(self, ins: float) -> float:
        """function toSIunit_h = toSIunit_h( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(2.32600 * ins)  # btu/lb to kJ/kg
        return float(ins)

    def fromSIunit_h(self, ins: float) -> float:
        """function  fromSIunit_h = fromSIunit_h( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 2.32600)  # kJ/kg to btu/lb
        return float(ins)

    def toSIunit_v(self, ins: float) -> float:
        """function toSIunit_v = toSIunit_v( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 0.0624279606)  # ft³/lb to m³/kg
        return float(ins)

    def fromSIunit_v(self, ins: float) -> float:
        """function fromSIunit_v = fromSIunit_v( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 0.0624279606)  # m³/kg to ft³/lb
        return float(ins)

    def toSIunit_s(self, ins: float) -> float:
        """function toSIunit_s = toSIunit_s( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 0.238845896627)  # btu/(lb degF) to kJ/(kg degC)
        return float(ins)

    def fromSIunit_s(self, ins: float) -> float:
        """function fromSIunit_s = fromSIunit_s( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 0.238845896627)  # kJ/(kg degC) to btu/(lb degF)
        return float(ins)

    def toSIunit_u(self, ins: float) -> float:
        """function toSIunit_u = toSIunit_u( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 2.32600)  # btu/lb to kJ/kg
        return float(ins)

    def fromSIunit_u(self, ins: float) -> float:
        """function fromSIunit_u = fromSIunit_u( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 2.32600)  # kJ/kg to btu/lb
        return float(ins)

    def toSIunit_Cp(self, ins: float) -> float:
        """function toSIunit_Cp = toSIunit_Cp( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 0.238846)  # btu/(lb degF) to kJ/(kg degC)
        return float(ins)

    def fromSIunit_Cp(self, ins: float) -> float:
        """function fromSIunit_Cp = fromSIunit_Cp( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 0.238846)  # kJ/(kg degC) to btu/(lb degF)
        return float(ins)

    def toSIunit_Cv(self, ins: float) -> float:
        """function toSIunit_Cv = toSIunit_Cv( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 0.238846)  # btu/(lb degF) to kJ/(kg degC)
        return float(ins)

    def fromSIunit_Cv(self, ins: float) -> float:
        """function fromSIunit_Cv = fromSIunit_Cv( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 0.238846)  # kJ/(kg degC) to btu/(lb degF)
        return float(ins)

    def toSIunit_w(self, ins: float) -> float:
        """function toSIunit_w = toSIunit_w( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 0.3048)  # ft/s to m/s
        return float(ins)

    def fromSIunit_w(self, ins: float) -> float:
        """function fromSIunit_w = fromSIunit_w( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 0.3048)  # m/s to ft/s
        return float(ins)

    def toSIunit_tc(self, ins: float) -> float:
        """function toSIunit_tc = toSIunit_tc( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 0.577789)  # btu/(h*ft*degF) to W/(m*degC)
        return float(ins)

    def fromSIunit_tc(self, ins: float) -> float:
        """function fromSIunit_tc = fromSIunit_tc( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 0.577789)  # W/(m*degC) to btu/(h*ft*degF)
        return float(ins)

    def toSIunit_st(self, ins: float) -> float:
        """function toSIunit_st = toSIunit_st( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 0.068521766)  # lb/ft to N/m
        return float(ins)

    def fromSIunit_st(self, ins: float) -> float:
        """function fromSIunit_st = fromSIunit_st( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 0.068521766)  # N/m to lb/ft
        return float(ins)

    def toSIunit_x(self, ins: float) -> float:
        """function toSIunit_x = toSIunit_x( ins )"""
        if 0.0 <= ins <= 1.0:
            return float(ins)
        self.logger.error("value of vapour fraction out of range: 0 < x < 1")
        raise ValueError("Vapour fraction out of Range")

    def fromSIunit_x(self, ins: float) -> float:
        """function fromSIunit_x = fromSIunit_x( ins )"""
        if 0.0 <= ins <= 1.0:
            return float(ins)
        self.logger.error("value of vapour fraction out of range: 0 < x < 1")
        raise ValueError("Vapour fraction out of Range")

    def toSIunit_vx(self, ins: float) -> float:
        """function toSIunit_vx = toSIunit_vx( ins )"""
        if 0.0 <= ins <= 1.0:
            return float(ins)
        self.logger.error("value of vapour volume fraction out of range: 0 < x < 1")
        raise ValueError("Vapour volume fraction out of Range")

    def fromSIunit_vx(self, ins: float) -> float:
        """function fromSIunit_vx = fromSIunit_vx( ins )"""
        if 0.0 <= ins <= 1.0:
            return float(ins)
        self.logger.error("value of vapour volume fraction out of range: 0 < x < 1")
        raise ValueError("Vapour volume fraction out of Range")

    def toSIunit_my(self, ins: float) -> float:
        """function toSIunit_my = toSIunit_my( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins / 2419.088311)  # lbm/ft/hr to PaS (N*s/m²)
        return float(ins)

    def fromSIunit_my(self, ins: float) -> float:
        """function fromSIunit_my = fromSIunit_my( ins )"""
        if self._unit_system is UnitSystem.FLS:
            return float(ins * 2419.088311)  # PaS (N*s/m²) to lbm/ft/hr
        return float(ins)

    def __str__(self):
        """returns string representation of the selected unit system"""
        result = ""
        if self._unit_system is UnitSystem.FLS:
            result = "FLS (ft/lb/sec/°F/psi/btu)"
        elif self._unit_system is UnitSystem.MKS:
            result = "MKS (m/kg/sec/°C/bar/W)"
        else:
            result = "BARE (m/kg/sec/K/MPa/W)"
        return result
