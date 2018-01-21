# -*- coding: utf-8 -*-
import logging
logger = logging.getLogger('pyXSteam-UnitConverter')

class UnitConverter(object):
    """Helper class to convert user Units to SI-Units and back"""
    __UNIT_SYSTEM_BARE__ = 0
    __UNIT_SYSTEM_MKS__ = 1
    __UNIT_SYSTEM_FLS__ = 2


    def __init__(self, unitSystem):
        self.logger = logging.getLogger(__name__)
        if unitSystem is self.__UNIT_SYSTEM_BARE__ or unitSystem is self.__UNIT_SYSTEM_MKS__ or unitSystem is self.__UNIT_SYSTEM_FLS__:
            self.unitSystem = unitSystem
        else:
            self.logger.critical('Unknown Unit System')
            raise ValueError('Unknown Unit System')


    def toSIunit_p(self, ins):
        """ function toSIunit_p = toSIunit_p( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins / 10)  # bar to MPa
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins * 0.00689475729)  # psi to MPa
        else:
            return float(ins)


    def fromSIunit_p(self, ins):
        """ function fromSIunit_p = fromSIunit_p( ins ) """
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins * 10)  # bar to MPa
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins / 0.00689475729)  # MPa to psi
        else:
            return float(ins)


    def toSIunit_T(self, ins):
        """ function toSIunit_T = toSIunit_T( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins + 273.15)  # degC to Kelvin
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float((5 / 9) * (ins - 32) + 273.15)  # degF to Kelvin
        else:
            return float(ins)


    def fromSIunit_T(self, ins):
        """ function fromSIunit_T = fromSIunit_T( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins - 273.15)  # Kelvin to degC
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float((ins - 273.15) * (9 / 5) + 32)  # Kelvin to degF
        else:
            return float(ins)


    def toSIunit_h(self, ins):
        """ function toSIunit_h = toSIunit_h( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(2.32600 * ins)  # btu/lb to kJ/kg
        else:
            return float(ins)


    def fromSIunit_h(self, ins):
        """ function  fromSIunit_h = fromSIunit_h( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins / 2.32600)  # kJ/kg to btu/lb
        else:
            return float(ins)


    def toSIunit_v(self, ins):
        """ function toSIunit_v = toSIunit_v( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins * 0.0624279606)  # ft^3/lb to m^3/kg
        else:
            return float(ins)


    def fromSIunit_v(self, ins):
        """ function fromSIunit_v = fromSIunit_v( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins / 0.0624279606)  # m^3/kg to ft^3/lb
        else:
            return float(ins)


    def toSIunit_s(self, ins):
        """ function toSIunit_s = toSIunit_s( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins / 0.238845896627)  # btu/(lb degF) to kJ/(kg degC)
        else:
            return float(ins)


    def fromSIunit_s(self, ins):
        """ function fromSIunit_s = fromSIunit_s( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins * 0.238845896627)  # kJ/(kg degC) to btu/(lb degF)
        else:
            return float(ins)


    def toSIunit_u(self, ins):
        """ function toSIunit_u = toSIunit_u( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins * 2.32600)  # btu/lb to kJ/kg
        else:
            return float(ins)


    def fromSIunit_u(self, ins):
        """ function fromSIunit_u = fromSIunit_u( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins / 2.32600)  # kJ/kg to btu/lb
        else:
            return float(ins)


    def toSIunit_Cp(self, ins):
        """ function toSIunit_Cp = toSIunit_Cp( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins / 0.238846)  # btu/(lb degF) to kJ/(kg degC)
        else:
            return float(ins)


    def fromSIunit_Cp(self, ins):
        """ function fromSIunit_Cp = fromSIunit_Cp( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins * 0.238846)  # kJ/(kg degC) to btu/(lb degF)
        else:
            return float(ins)


    def toSIunit_Cv(self, ins):
        """ function toSIunit_Cv = toSIunit_Cv( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins / 0.238846)  # btu/(lb degF) to kJ/(kg degC)
        else:
            return float(ins)


    def fromSIunit_Cv(self, ins):
        """ function fromSIunit_Cv = fromSIunit_Cv( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins * 0.238846)  # kJ/(kg degC) to btu/(lb degF)
        else:
            return float(ins)


    def toSIunit_w(self, ins):
        """ function toSIunit_w = toSIunit_w( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins * 0.3048)  # ft/s to m/s
        else:
            return float(ins)


    def fromSIunit_w(self, ins):
        """ function fromSIunit_w = fromSIunit_w( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins / 0.3048)  # m/s to ft/s
        else:
            return float(ins)


    def toSIunit_tc(self, ins):
        """ function toSIunit_tc = toSIunit_tc( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins / 0.577789)  # btu/(h*ft*degF) to W/(m*degC)
        else:
            return float(ins)


    def fromSIunit_tc(self, ins):
        """ function fromSIunit_tc = fromSIunit_tc( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins * 0.577789)  # W/(m*degC) to btu/(h*ft*degF)
        else:
            return float(ins)


    def toSIunit_st(self, ins):
        """ function toSIunit_st = toSIunit_st( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins / 0.068521766)  # lb/ft to N/m
        else:
            return float(ins)


    def fromSIunit_st(self, ins):
        """ function fromSIunit_st = fromSIunit_st( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins * 0.068521766)  # N/m to lb/ft
        else:
            return float(ins)


    def toSIunit_x(self, ins):
        """ function toSIunit_x = toSIunit_x( ins )"""
        if ins >= 0.0 and ins <= 1.0:
            return float(ins)
        else:
            raise Exception('Vapour fraction out of Range')


    def fromSIunit_x(self, ins):
        """ function fromSIunit_x = fromSIunit_x( ins )"""
        if ins >= 0.0 and ins <= 1.0:
            return float(ins)
        else:
            raise Exception('Vapour fraction out of Range')


    def toSIunit_vx(self, ins):
        """ function toSIunit_vx = toSIunit_vx( ins )"""
        if ins >= 0.0 and ins <= 1.0:
            return float(ins)
        else:
            raise Exception('Vapour volume fraction out of Range')


    def fromSIunit_vx(self, ins):
        """ function fromSIunit_vx = fromSIunit_vx( ins )"""
        if ins >= 0.0 and ins <= 1.0:
            return float(ins)
        else:
            raise Exception('Vapour volume fraction out of Range')


    def toSIunit_my(self, ins):
        """ function toSIunit_my = toSIunit_my( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins / 2419.088311)  # lbm/ft/hr to PaS (N*s/m^2)
        else:
            return float(ins)


    def fromSIunit_my(self, ins):
        """ function fromSIunit_my = fromSIunit_my( ins )"""
        if self.unitSystem is self.__UNIT_SYSTEM_MKS__:
            return float(ins)
        elif self.unitSystem is self.__UNIT_SYSTEM_FLS__:
            return float(ins * 2419.088311)  # PaS (N*s/m^2) to lbm/ft/hr
        else:
            return float(ins)
