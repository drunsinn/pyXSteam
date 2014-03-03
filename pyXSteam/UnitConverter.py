# -*- coding: utf-8 -*-

class UnitConverter(object):
    """Helper class to convert user Units to SI-Units and back"""


    def __init__(self, mksSystem = True):
        self.mksSystem = mksSystem


    def toSIunit_p(self, ins):
        """ function toSIunit_p = toSIunit_p( ins )"""
        if self.mksSystem:
            # bar to MPa
            return float(ins / 10)
        else:
            # psi to MPa
            return float(ins * 0.00689475729)


    def fromSIunit_p(self, ins):
        """ function fromSIunit_p = fromSIunit_p( ins ) """
        if self.mksSystem:
            # bar to MPa
            return float(ins * 10)
        else:
            # MPa to psi
            return float(ins / 0.00689475729)


    def toSIunit_T(self, ins):
        """ function toSIunit_T = toSIunit_T( ins )"""
        if self.mksSystem:
            # degC to Kelvin
            return float(ins + 273.15)
        else:
            # degF to Kelvin
            return float((5 / 9) * (ins - 32) + 273.15)


    def fromSIunit_T(self, ins):
        """ function fromSIunit_T = fromSIunit_T( ins )"""
        if self.mksSystem:
            # Kelvin to degC
            return float(ins - 273.15)
        else:
            # Kelvin to degF
            return float((ins - 273.15) * (9 / 5) + 32)


    def toSIunit_h(self, ins):
        """ function toSIunit_h = toSIunit_h( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # btu/lb to kJ/kg
            return float(2.32600 * ins)


    def fromSIunit_h(self, ins):
        """ function  fromSIunit_h = fromSIunit_h( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # kJ/kg to btu/lb
            return float(ins / 2.32600)


    def toSIunit_v(self, ins):
        """ function toSIunit_v = toSIunit_v( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # ft^3/lb to m^3/kg
            return float(ins * 0.0624279606)


    def fromSIunit_v(self, ins):
        """ function fromSIunit_v = fromSIunit_v( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # m^3/kg to ft^3/lb
            return float(ins / 0.0624279606)


    def toSIunit_s(self, ins):
        """ function toSIunit_s = toSIunit_s( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # btu/(lb degF) to kJ/(kg degC)
            return float(ins / 0.238845896627)


    def fromSIunit_s(self, ins):
        """ function fromSIunit_s = fromSIunit_s( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # kJ/(kg degC) to btu/(lb degF)
            return float(ins * 0.238845896627)


    def toSIunit_u(self, ins):
        """ function toSIunit_u = toSIunit_u( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # btu/lb to kJ/kg
            return float(ins * 2.32600)


    def fromSIunit_u(self, ins):
        """ function fromSIunit_u = fromSIunit_u( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # kJ/kg to btu/lb
            return float(ins / 2.32600)


    def toSIunit_Cp(self, ins):
        """ function toSIunit_Cp = toSIunit_Cp( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # btu/(lb degF) to kJ/(kg degC)
            return float(ins / 0.238846)


    def fromSIunit_Cp(self, ins):
        """ function fromSIunit_Cp = fromSIunit_Cp( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # kJ/(kg degC) to btu/(lb degF)
            return float(ins * 0.238846)


    def toSIunit_Cv(self, ins):
        """ function toSIunit_Cv = toSIunit_Cv( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # btu/(lb degF) to kJ/(kg degC)
            return float(ins / 0.238846)


    def fromSIunit_Cv(self, ins):
        """ function fromSIunit_Cv = fromSIunit_Cv( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # kJ/(kg degC) to btu/(lb degF)
            return float(ins * 0.238846)


    def toSIunit_w(self, ins):
        """ function toSIunit_w = toSIunit_w( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # ft/s to m/s
            return float(ins * 0.3048)


    def fromSIunit_w(self, ins):
        """ function fromSIunit_w = fromSIunit_w( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # m/s to ft/s
            return float(ins / 0.3048)


    def toSIunit_tc(self, ins):
        """ function toSIunit_tc = toSIunit_tc( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # btu/(h*ft*degF) to W/(m*degC)
            return float(ins / 0.577789)


    def fromSIunit_tc(self, ins):
        """ function fromSIunit_tc = fromSIunit_tc( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # W/(m*degC) to btu/(h*ft*degF)
            return float(ins * 0.577789)


    def toSIunit_st(self, ins):
        """ function toSIunit_st = toSIunit_st( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # lb/ft to N/m
            return float(ins / 0.068521766)


    def fromSIunit_st(self, ins):
        """ function fromSIunit_st = fromSIunit_st( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # N/m to lb/ft
            return float(ins * 0.068521766)


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
        if self.mksSystem:
            return float(ins)
        else:
            #  lbm/ft/hr to PaS (N*s/m^2)
            return float(ins / 2419.088311)


    def fromSIunit_my(self, ins):
        """ function fromSIunit_my = fromSIunit_my( ins )"""
        if self.mksSystem:
            return float(ins)
        else:
            # PaS (N*s/m^2) to lbm/ft/hr
            return float(ins * 2419.088311)
