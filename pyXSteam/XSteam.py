#!/usr/bin/python
# -*- coding: UTF-8 -*-
import math
# import Region1
# import Region2
# import Region3
# import Region4
# import Region5
from Regions import Region1, Region2, Region3, Region4, Region5
# from pyXSteam import UnitConverter
import RegionSelection
import TransportProperties

class XSteam(object):
    ''' Library for calculation properties of H2O accordint to IF-97
        This Python Library is based on the original XSteam Library for 
        Matlab and Excel from Magnus Holmgren, www.x-eng.com.
        Es take no responsibilities for any errors in the code or damage thereby!
        See Readme for examples
    '''
    def __init__(self, mksSystem = True):
        self.unitConverter = UnitConverter(mksSystem)

#%***********************************************************************************************************
# %*1.2 Temperature
    def tsat_p(self, p):
        '''Saturation-Temperature as a function of preasure'''
        # case 'tsat_p'
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            return self.unitConverter.fromSIunit_T(Region4.T4_p(p))
        else:
            return float("NaN")

    def tsat_s(self, s):
        '''Saturation-Temperature as a function of entropy'''
        # case 'tsat_s'
        s = self.unitConverter.toSIunit_s(s)
        if (s > -0.0001545495919) and (s < 9.155759395):
            ps = Region4.p4_s(s)
            return self.unitConverter.fromSIunit_T(Region4.T4_p(ps))
        else:
            return float("NaN")

    def t_ph(self, p, h):
        '''Temperture as a function of pressure and enthalpy'''
        # case 't_ph'
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        Region = RegionSelection.region_ph(p, h)
        if Region == 1:
            return self.unitConverter.fromSIunit_T(Region1.T1_ph(p, h))
        elif Region == 2:
            return self.unitConverter.fromSIunit_T(Region2.T2_ph(p, h))
        elif Region == 3:
            return self.unitConverter.fromSIunit_T(Region3.T3_ph(p, h))
        elif Region == 4:
            return self.unitConverter.fromSIunit_T(Region4.T4_p(p))
        elif Region == 5:
            return self.unitConverter.fromSIunit_T(Region5.T5_ph(p, h))
        else:
            return float("NaN")

    def t_ps(self, p, s):
        '''Temperture as a function of pressure and entropy'''
        # case 't_ps'
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        Region = RegionSelection.region_ps(p, s)
        if Region == 1:
            return self.unitConverter.fromSIunit_T(Region1.T1_ps(p, s))
        elif Region == 2:
            return self.unitConverter.fromSIunit_T(Region2.T2_ps(p, s))
        elif Region == 3:
            return self.unitConverter.fromSIunit_T(Region3.T3_ps(p, s))
        elif Region == 4:
            return self.unitConverter.fromSIunit_T(Region4.T4_p(p))
        elif Region == 5:
            return self.unitConverter.fromSIunit_T(Region5.T5_ps(p, s))
        else:
            return float("NaN")

    def t_hs(self, h, s):
        '''Temperture as a function of enthalpy and entropy'''
        # case 't_hs'
        h = self.unitConverter.toSIunit_h(h)
        s = self.unitConverter.toSIunit_s(s)
        Region = RegionSelection.region_hs(h, s)
        if Region == 1:
            p1 = Region1.p1_hs(h, s)
            return self.unitConverter.fromSIunit_T(Region1.T1_ph(p1, h))
        elif Region == 2:
            p2 = Region2.p2_hs(h, s)
            return self.unitConverter.fromSIunit_T(Region2.T2_ph(p2, h))
        elif Region == 3:
            p3 = Region3.p3_hs(h, s)
            return self.unitConverter.fromSIunit_T(Region3.T3_ph(p3, h))
        elif Region == 4:
            return self.unitConverter.fromSIunit_T(Region4.T4_hs(h, s))
        elif Region == 5:
            raise Exception('functions of hs is not avlaible in region 5')
        else:
            return float("NaN")

#%***********************************************************************************************************
# %*1.3 Pressure (p)
    def psat_s(self, s):
        '''Saturation-Pressure as a function of entropy'''
        # case 'psat_s'
        s = self.unitConverter.toSIunit_s(s)
        if (s > -0.0001545495919) and (s < 9.155759395):
            return self.unitConverter.fromSIunit_p(Region4.p4_s(s))
        else:
            return float("NaN")

    def psat_t(self, t):
        '''Saturation-Pressure as a function of temperature'''
        # case 'psat_t'
        T = self.unitConverter.toSIunit_T(t)
        if (T < 647.096) and (T > 273.1):
            return self.unitConverter.fromSIunit_p(Region4.p4_T(T))
        else:
            return float("NaN")

    def p_hs(self, h, s):
        '''Pressure  as a function of enthalpy and entropy'''
        # case 'p_hs'
        h = self.unitConverter.toSIunit_h(h)
        s = self.unitConverter.toSIunit_s(s)
        Region = RegionSelection.region_hs(h, s)
        if Region == 1:
            return self.unitConverter.fromSIunit_p(Region1.p1_hs(h, s))
        elif Region == 2:
            return self.unitConverter.fromSIunit_p(Region2.p2_hs(h, s))
        elif Region == 3:
            return self.unitConverter.fromSIunit_p(Region3.p3_hs(h, s))
        elif Region == 4:
            tSat = Region4.T4_hs(h, s)
            return self.unitConverter.fromSIunit_p(Region4.p4_T(tSat))
        elif Region == 5:
            raise Exception('functions of hs is not avlaible in region 5')
        else:
            return float("NaN")


    def p_hrho(self, h, rho):
        '''Pressure as a function of h and rho. Very unaccurate for solid water region since it's almost incompressible!'''
        '''Not valid for water or sumpercritical since water rho does not change very much with p. Uses iteration to find p.'''
        # case 'p_hrho'
        h = self.unitConverter.toSIunit_h(h)
        High_Bound = self.unitConverter.fromSIunit_p(100)
        Low_Bound = self.unitConverter.fromSIunit_p(0.000611657)
        ps = self.unitConverter.fromSIunit_p(10)
        rhos = 1 / self.v_ph(ps, h)
        while math.fabs(rho - rhos) > 0.0000001:
            # rhos = 1 / XSteam('v_ph', ps, h)
            rhos = 1 / self.v_ph(ps, h)
            if rhos >= rho:
                High_Bound = ps
            else:
                Low_Bound = ps
            ps = (Low_Bound + High_Bound) / 2
        return ps;


#%***********************************************************************************************************
# %*1.4 Enthalpy (h)
    def hV_p(self, p):
        '''Saturated vapour enthalpy as a function of preasure'''
        # case 'hv_p'
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            return self.unitConverter.fromSIunit_h(Region4.h4V_p(p))
        else:
            return float("NaN")

    def hL_p(self, p):
        '''Saturated liquid enthalpy as a function of preasure'''
        # case 'hl_p'
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            return self.unitConverter.fromSIunit_h(Region4.h4L_p(p))
        else:
            return float("NaN")

    def hV_t(self, t):
        '''Saturated vapour enthalpy as a function of temperature'''
        # case 'hv_t'
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            p = Region4.p4_T(T)
            return self.unitConverter.fromSIunit_h(Region4.h4V_p(p))
        else:
            return float("NaN")

    def hL_t(self, t):
        '''Saturated liquid enthalpy as a function of temperature'''
        # case 'hl_t'
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            p = Region4.p4_T(T)
            return self.unitConverter.fromSIunit_h(Region4.h4L_p(p))
        else:
            return float("NaN")

    def h_pt(self, p, t):
        '''Entalpy as a function of pressure and temperature'''
        # case 'h_pt'
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        Region = RegionSelection.region_pT(p, T)
        if Region == 1:
            return self.unitConverter.fromSIunit_h(Region1.h1_pT(p, T))
        elif Region == 2:
            return self.unitConverter.fromSIunit_h(Region2.h2_pT(p, T))
        elif Region == 3:
            return self.unitConverter.fromSIunit_h(Region3.h3_pT(p, T))
        elif Region == 4:
            return float("NaN")
        elif Region == 5:
            return self.unitConverter.fromSIunit_h(Region5.h5_pT(p, T))
        else:
            return float("NaN")

    def h_ps(self, p, s):
        '''Entalpy as a function of pressure and entropy'''
        # case 'h_ps'
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        Region = RegionSelection.region_ps(p, s)
        if Region == 1:
            return self.unitConverter.fromSIunit_h(Region1.h1_pT(p, Region1.T1_ps(p, s)))
        elif Region == 2:
            return self.unitConverter.fromSIunit_h(Region2.h2_pT(p, Region2.T2_ps(p, s)))
        elif Region == 3:
            return self.unitConverter.fromSIunit_h(Region3.h3_rhoT(1 / Region3.v3_ps(p, s), Region3.T3_ps(p, s)))
        elif Region == 4:
            xs = Region4.x4_ps(p, s)
            return self.unitConverter.fromSIunit_h(xs * Region4.h4V_p(p) + (1 - xs) * Region4.h4L_p(p))
        elif Region == 5:
            return self.unitConverter.fromSIunit_h(Region5.h5_pT(p, Region5.T5_ps(p, s)))
        else:
            return float("NaN")

    def h_px(self, p, x):
        '''Entalpy as a function of pressure and vapour fraction'''
        # case 'h_px'
        p = self.unitConverter.toSIunit_p(p)
        x = self.unitConverter.toSIunit_x(x)
        if (x > 1) or (x < 0) or (p >= 22.064):
            return float("NaN")
        hL = Region4.h4L_p(p)
        hV = Region4.h4V_p(p)
        return hL + x * (hV - hL)

    def h_prho(self, p, rho):
        '''Entalpy as a function of pressure and density. Observe for low temperatures (liquid) this equation has 2 solutions'''
        # case 'h_prho'
        p = self.unitConverter.toSIunit_p(p)
        rho = 1 / self.unitConverter.toSIunit_v(1 / float(rho))
        Region = RegionSelection.Region_prho(p, rho)
        if Region == 1:
            return self.unitConverter.fromSIunit_h(Region1.h1_pT(p, Region1.T1_prho(p, rho)))
        elif Region == 2:
            return self.unitConverter.fromSIunit_h(Region2.h2_pT(p, Region2.T2_prho(p, rho)))
        elif Region == 3:
            return self.unitConverter.fromSIunit_h(Region3.h3_rhoT(rho, Region3.T3_prho(p, rho)))
        elif Region == 4:
                if p < 16.529:
                    vV = Region2.v2_pT(p, Region4.T4_p(p))
                    vL = Region1.v1_pT(p, Region4.T4_p(p))
                else:
                    vV = Region3.v3_ph(p, Region4.h4V_p(p))
                    vL = Region3.v3_ph(p, Region4.h4L_p(p))
                hV = Region4.h4V_p(p)
                hL = Region4.h4L_p(p)
                x = (1 / rho - vL) / (vV - vL)
                return self.unitConverter.fromSIunit_h((1 - x) * hL + x * hV)
        elif Region == 5:
            return self.unitConverter.fromSIunit_h(Region5.h5_pT(p, Region5.T5_prho(p, rho)))
        else:
            return float("NaN")

    def h_tx(self, t, x):
        '''Entalpy as a function of temperature and vapour fraction'''
        # case 'h_tx'
        T = self.unitConverter.toSIunit_T(t)
        x = self.unitConverter.toSIunit_x(x)
        if (x > 1) or (x < 0) or (T >= 647.096):
            return float("NaN")
        p = Region4.p4_T(T)
        hL = Region4.h4L_p(p)
        hV = Region4.h4V_p(p)
        return hL + x * (hV - hL)

#%***********************************************************************************************************
# %*1.5 Specific Volume (v)
    def vV_p(self, p):
        '''Saturated vapour volume as a function of preasure'''
        # case {'vv_p','rhov_p'}
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_v(Region2.v2_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_v(Region3.v3_ph(p, Region4.h4V_p(p)))
        else:
            return float("NaN")

    def vL_p(self, p):
        '''Saturated liquid volume as a function of preasure'''
        # case {'vl_p', 'rhol_p'}
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_v(Region1.v1_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_v(Region3.v3_ph(p, Region4.h4L_p(p)))
        else:
            return float("NaN")

    def vV_t(self, t):
        '''Saturated vapour volume as a function of temperature'''
        # case {'vv_t', 'rhov_t'}
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if T <= 623.15:
                return self.unitConverter.fromSIunit_v(Region2.v2_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_v(Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T))))
        else:
            return float("NaN")

    def vL_t(self, t):
        '''Saturated liquid volume as a function of temperature'''
        # case {'vl_t', 'rhol_t'}
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if T <= 623.15:
                return self.unitConverter.fromSIunit_v(Region1.v1_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_v(Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T))))
        else:
            return float("NaN")

    def v_pt(self, p, t):
        '''Specific volume as a function of pressure and temperature'''
        # case {'v_pt', 'rho_pt'}
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        Region = RegionSelection.region_pT(p, T)
        if Region == 1:
            return self.unitConverter.fromSIunit_v(Region1.v1_pT(p, T))
        elif Region == 2:
            return self.unitConverter.fromSIunit_v(Region2.v2_pT(p, T))
        elif Region == 3:
            return self.unitConverter.fromSIunit_v(Region3.v3_ph(p, Region3.h3_pT(p, T)))
        elif Region == 4:
            return float("NaN")
        elif Region == 5:
            return self.unitConverter.fromSIunit_v(Region5.v5_pT(p, T))
        else:
            return float("NaN")

    def v_ph(self, p, h):
        '''Specific volume as a function of pressure and enthalpy'''
        # case {'v_ph', 'rho_ph'}
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        Region = RegionSelection.region_ph(p, h)
        if Region == 1:
            return self.unitConverter.fromSIunit_v(Region1.v1_pT(p, Region1.T1_ph(p, h)))
        elif Region == 2:
            return self.unitConverter.fromSIunit_v(Region2.v2_pT(p, Region2.T2_ph(p, h)))
        elif Region == 3:
            return self.unitConverter.fromSIunit_v(Region3.v3_ph(p, h))
        elif Region == 4:
            xs = Region4.x4_ph(p, h)
            if p < 16.529:
                v4v = Region2.v2_pT(p, Region4.T4_p(p))
                v4L = Region1.v1_pT(p, Region4.T4_p(p))
            else:
                v4v = Region3.v3_ph(p, Region4.h4V_p(p))
                v4L = Region3.v3_ph(p, Region4.h4L_p(p))
            return self.unitConverter.fromSIunit_v((xs * v4v + (1 - xs) * v4L))
        elif Region == 5:
            Ts = Region5.T5_ph(p, h)
            return self.unitConverter.fromSIunit_v(Region5.v5_pT(p, Ts))
        else:
            return float("NaN")

    def v_ps(self, p, s):
        '''Specific volume as a function of pressure and entropy'''
        # case {'v_ps', 'rho_ps'}
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        Region = RegionSelection.region_ps(p, s)
        if Region == 1:
            return self.unitConverter.fromSIunit_v(Region1.v1_pT(p, Region1.T1_ps(p, s)))
        elif Region == 2:
            return self.unitConverter.fromSIunit_v(Region2.v2_pT(p, Region2.T2_ps(p, s)))
        elif Region == 3:
            return self.unitConverter.fromSIunit_v(Region3.v3_ps(p, s))
        elif Region == 4:
            xs = Region4.x4_ps(p, s)
            if p < 16.529:
                v4v = Region2.v2_pT(p, Region4.T4_p(p))
                v4L = Region1.v1_pT(p, Region4.T4_p(p))
            else:
                v4v = Region3.v3_ph(p, Region4.h4V_p(p))
                v4L = Region3.v3_ph(p, Region4.h4L_p(p))
            return self.unitConverter.fromSIunit_v((xs * v4v + (1 - xs) * v4L))
        elif Region == 5:
            Ts = Region5.T5_ps(p, s)
            return self.unitConverter.fromSIunit_v(Region5.v5_pT(p, Ts))
        else:
            return float("NaN")

#%***********************************************************************************************************
# %*1.6 Density (rho)
# % Density is calculated as 1/v. Se section 1.5 Volume
    def rhoV_p(self, p):
        '''Saturated vapour density as a function of preasure'''
        return 1 / self.vV_p(p)

    def rhoL_p(self, p):
        '''Saturated liquid density as a function of preasure'''
        return 1 / self.vL_p(p)

    def rhoV_t(self, t):
        '''Saturated vapour density as a function of temperature'''
        return 1 / self.vV_t(t)

    def rhoL_t(self, t):
        '''Saturated liquid density as a function of temperature'''
        return 1 / self.vL_t(t)

    def rho_pt(self, p, t):
        '''Density as a function of pressure and temperature'''
        return 1 / self.v_pt(p, t)

    def rho_ph(self, p, h):
        '''Density as a function of pressure and enthalpy'''
        return 1 / self.v_pt(p, h)

    def rho_ps(self, p, s):
        '''Density as a function of pressure and entropy'''
        return 1 / self.v_ps(p, s)

#%***********************************************************************************************************
# %*1.7 Specific entropy (s)
    def sV_p(self, p):
        '''Saturated vapour entropy as a function of preasure'''
        # case 'sv_p'
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_s(Region2.s2_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_s(Region3.s3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p)))
        else:
            return float("NaN")

    def sL_p(self, p):
        '''Saturated liquid entropy as a function of preasure'''
        # case 'sl_p'
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_s(Region1.s1_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_s(Region3.s3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p)))
        else:
            return float("NaN")

    def sV_t(self, t):
        '''Saturated vapour entropy as a function of temperature'''
        # case 'sv_t'
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if T <= 623.15:
                return self.unitConverter.fromSIunit_s(Region2.s2_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_s(Region3.s3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T)))), T))
        else:
            return float("NaN")

    def sL_t(self, t):
        '''Saturated liquid entropy as a function of temperature'''
        # case 'sl_t'
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if T <= 623.15:
                return self.unitConverter.fromSIunit_s(Region1.s1_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_s(Region3.s3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T)))), T))
        else:
            return float("NaN")

    def s_pt(self, p, t):
        '''Specific entropy as a function of pressure and temperature (Returns saturated vapour entalpy if mixture)'''
        # case 's_pt'
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        Region = RegionSelection.region_pT(p, T)
        if Region == 1:
            return self.unitConverter.fromSIunit_s(Region1.s1_pT(p, T))
        elif Region == 2:
            return self.unitConverter.fromSIunit_s(Region2.s2_pT(p, T))
        elif Region == 3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self.unitConverter.fromSIunit_s(Region3.s3_rhoT(rhos, T))
        elif Region == 4:
            return float("NaN")
        elif Region == 5:
            return self.unitConverter.fromSIunit_s(Region5.s5_pT(p, T))
        else:
            return float("NaN")

    def s_ph(self, p, h):
        '''Specific entropy as a function of pressure and enthalpy'''
        # case 's_ph'
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        Region = RegionSelection.region_ph(p, h)
        if Region == 1:
            T = Region1.T1_ph(p, h)
            return self.unitConverter.fromSIunit_s(Region1.s1_pT(p, T))
        elif Region == 2:
            T = Region2.T2_ph(p, h)
            return self.unitConverter.fromSIunit_s(Region2.s2_pT(p, T))
        elif Region == 3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self.unitConverter.fromSIunit_s(Region3.s3_rhoT(rhos, Ts))
        elif Region == 4:
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
            return self.unitConverter.fromSIunit_s((xs * s4v + (1 - xs) * s4L))
        elif Region == 5:
            T = Region5.T5_ph(p, h)
            return self.unitConverter.fromSIunit_s(Region5.s5_pT(p, T))
        else:
            return float("NaN")

#%***********************************************************************************************************
# %*1.8 Specific internal energy (u)
    def uV_p(self, p):
        '''Saturated vapour internal energy as a function of preasure'''
        # case 'uv_p'
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_u(Region2.u2_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_u(Region3.u3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p)))
        else:
            return float("NaN")


    def uL_p(self, p):
        '''Saturated liquid internal energy as a function of preasure'''
        # case 'ul_p'
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_u(Region1.u1_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_u(Region3.u3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p)))
        else:
            return float("NaN")

    def uV_t(self, t):
        '''Saturated vapour internal energy as a function of temperature'''
        # case 'uv_t'
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if T <= 623.15:
                return self.unitConverter.fromSIunit_u(Region2.u2_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_u(Region3.u3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T)))), T))
        else:
            return float("NaN")

    def uL_t(self, t):
        '''Saturated liquid internal energy as a function of temperature'''
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if T <= 623.15:
                return self.unitConverter.fromSIunit_u(Region1.u1_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_u(Region3.u3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T)))), T))
        else:
            return float("NaN")

    def u_pt(self, p, t):
        '''Specific internal energy as a function of pressure and temperature'''
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        Region = RegionSelection.region_pT(p, T)
        if Region == 1:
            return self.unitConverter.fromSIunit_u(Region1.u1_pT(p, T))
        elif Region == 2:
            return self.unitConverter.fromSIunit_u(Region2.u2_pT(p, T))
        elif Region == 3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self.unitConverter.fromSIunit_u(Region3.u3_rhoT(rhos, T))
        elif Region == 4:
            return float("NaN")
        elif Region == 5:
            return self.unitConverter.fromSIunit_u(Region5.u5_pT(p, T))
        else:
            return float("NaN")

    def u_ph(self, p, h):
        '''Specific internal energy as a function of pressure and enthalpy'''
        # case 'u_ph'
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        Region = RegionSelection.region_ph(p, h)
        if Region == 1:
            Ts = Region1.T1_ph(p, h)
            return self.unitConverter.fromSIunit_u(Region1.u1_pT(p, Ts))
        elif Region == 2:
            Ts = Region2.T2_ph(p, h)
            return self.unitConverter.fromSIunit_u(Region2.u2_pT(p, Ts))
        elif Region == 3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self.unitConverter.fromSIunit_u(Region3.u3_rhoT(rhos, Ts))
        elif Region == 4:
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
            return self.unitConverter.fromSIunit_u((xs * u4v + (1 - xs) * u4L))
        elif Region == 5:
            Ts = Region5.T5_ph(p, h)
            return self.unitConverter.fromSIunit_u(Region5.u5_pT(p, Ts))
        else:
            return float("NaN")

    def u_ps(self, p, s):
        '''Specific internal energy as a function of pressure and entropy'''
        # case 'u_ps'
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        Region = RegionSelection.region_ps(p, s)
        if Region == 1:
            Ts = Region1.T1_ps(p, s)
            return self.unitConverter.fromSIunit_u(Region1.u1_pT(p, Ts))
        elif Region == 2:
            Ts = Region2.T2_ps(p, s)
            return self.unitConverter.fromSIunit_u(Region2.u2_pT(p, Ts))
        elif Region == 3:
            rhos = 1 / Region3.v3_ps(p, s)
            Ts = Region3.T3_ps(p, s)
            return self.unitConverter.fromSIunit_u(Region3.u3_rhoT(rhos, Ts))
        elif Region == 4:
            if p < 16.529:
                uLp = Region1.u1_pT(p, Region4.T4_p(p))
                uVp = Region2.u2_pT(p, Region4.T4_p(p))
            else:
                uLp = Region3.u3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p))
                uVp = Region3.u3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p))
            xs = Region4.x4_ps(p, s)
            return self.unitConverter.fromSIunit_u((xs * uVp + (1 - xs) * uLp))
        elif Region == 5:
            Ts = Region5.T5_ps(p, s)
            return self.unitConverter.fromSIunit_u(Region5.u5_pT(p, Ts))
        else:
            return float("NaN")

#%***********************************************************************************************************
# %*1.9 Specific isobaric heat capacity (Cp)
    def CpV_p(self, p):
        '''Saturated vapour heat capacity as a function of preasure'''
        # case 'cpv_p'
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_Cp(Region2.Cp2_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_Cp(Region3.Cp3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p)))
        else:
            return float("NaN")

    def CpL_p(self, p):
        '''Saturated liquid heat capacity as a function of preasure'''
        # case 'cpl_p'
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_Cp(Region1.Cp1_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_Cp(Region3.Cp3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p)))
        else:
            return float("NaN")

    def CpV_t(self, t):
        '''Saturated vapour heat capacity as a function of temperature'''
        # case 'cpv_t'
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if (T <= 623.15):
                return self.unitConverter.fromSIunit_Cp(Region2.Cp2_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_Cp(Region3.Cp3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T)))), T))
        else:
            return float("NaN")

    def CpL_t(self, t):
        '''Saturated liquid heat capacity as a function of temperature'''
        # case 'cpl_t'
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if (T <= 623.15):
                return self.unitConverter.fromSIunit_Cp(Region1.Cp1_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_Cp(Region3.Cp3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T)))), T))
        else:
            return float("NaN")

    def Cp_pt(self, p, t):
        '''Specific isobaric heat capacity as a function of pressure and temperature'''
        # case 'cp_pt'
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        Region = RegionSelection.region_pT(p, T)
        if Region == 1:
            return self.unitConverter.fromSIunit_Cp(Region1.Cp1_pT(p, T))
        elif Region == 2:
            return self.unitConverter.fromSIunit_Cp(Region2.Cp2_pT(p, T))
        elif Region == 3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self.unitConverter.fromSIunit_Cp(Region3.Cp3_rhoT(rhos, T))
        elif Region == 4:
            return float("NaN")
        elif Region == 5:
            return self.unitConverter.fromSIunit_Cp(Region5.Cp5_pT(p, T))
        else:
            return float("NaN")

    def Cp_ph(self, p, h):
        '''Specific isobaric heat capacity as a function of pressure and enthalpy'''
        # case 'cp_ph'
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        Region = RegionSelection.region_ph(p, h)
        if Region == 1:
            Ts = Region1.T1_ph(p, h)
            return self.unitConverter.fromSIunit_Cp(Region1.Cp1_pT(p, Ts))
        elif Region == 2:
            Ts = Region2.T2_ph(p, h)
            return self.unitConverter.fromSIunit_Cp(Region2.Cp2_pT(p, Ts))
        elif Region == 3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self.unitConverter.fromSIunit_Cp(Region3.Cp3_rhoT(rhos, Ts))
        elif Region == 4:
            return float("NaN")
        elif Region == 5:
            Ts = Region5.T5_ph(p, h)
            return self.unitConverter.fromSIunit_Cp(Region5.Cp5_pT(p, Ts))
        else:
            return float("NaN")

    def Cp_ps(self, p, s):
        '''Specific isobaric heat capacity as a function of pressure and entropy'''
        # case 'cp_ps'
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        Region = RegionSelection.region_ps(p, s)
        if Region == 1:
            Ts = Region1.T1_ps(p, s)
            return self.unitConverter.fromSIunit_Cp(Region1.Cp1_pT(p, Ts))
        elif Region == 2:
            Ts = Region2.T2_ps(p, s)
            return self.unitConverter.fromSIunit_Cp(Region2.Cp2_pT(p, Ts))
        elif Region == 3:
            rhos = 1 / Region3.v3_ps(p, s)
            Ts = Region3.T3_ps(p, s)
            return self.unitConverter.fromSIunit_Cp(Region3.Cp3_rhoT(rhos, Ts))
        elif Region == 4:
            return float("NaN")
        elif Region == 5:
            Ts = Region5.T5_ps(p, s)
            return self.unitConverter.fromSIunit_Cp(Region5.Cp5_pT(p, Ts))
        else:
            return float("NaN")

#%***********************************************************************************************************
# %*1.10 Specific isochoric heat capacity (Cv)
    def CvV_p(self, p):
        '''Saturated vapour isochoric heat capacity as a function of preasure'''
        # case 'cvv_p'
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_Cv(Region2.Cv2_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_Cv(Region3.Cv3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p)))
        else:
            return float("NaN")

    def CvL_p(self, p):
        '''Saturated liquid isochoric heat capacity as a function of preasure'''
        # case 'cvl_p'
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_Cv(Region1.Cv1_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_Cv(Region3.Cv3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p)))
        else:
            return float("NaN")

    def CvV_t(self, t):
        '''Saturated vapour isochoric heat capacity as a function of temperature'''
        # case 'cvv_t'
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if (T <= 623.15):
                return self.unitConverter.fromSIunit_Cv(Region2.Cv2_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_Cv(Region3.Cv3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T)))), T))
        else:
            return float("NaN")

    def CvL_t(self, t):
        '''Saturated liquid isochoric heat capacity as a function of temperature'''
        # case 'cvl_t'
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if (T <= 623.15):
                return self.unitConverter.fromSIunit_Cv(Region1.Cv1_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_Cv(Region3.Cv3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T)))), T))
        else:
            return float("NaN")

    def Cv_pt(self, p, t):
        '''Specific isochoric heat capacity as a function of pressure and temperature'''
        # case 'cv_pt'
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        Region = RegionSelection.region_pT(p, T)
        if Region == 1:
            return self.unitConverter.fromSIunit_Cv(Region1.Cv1_pT(p, T))
        elif Region == 2:
            return self.unitConverter.fromSIunit_Cv(Region2.Cv2_pT(p, T))
        elif Region == 3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self.unitConverter.fromSIunit_Cv(Region3.Cv3_rhoT(rhos, T))
        elif Region == 4:
            return float("NaN")
        elif Region == 5:
            return self.unitConverter.fromSIunit_Cv(Region5.Cv5_pT(p, T))
        else:
            return float("NaN")

    def Cv_ph(self, p, h):
        '''Specific isochoric heat capacity as a function of pressure and enthalpy'''
        # case 'cv_ph'
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        Region = RegionSelection.region_ph(p, h)
        if Region == 1:
            Ts = Region1.T1_ph(p, h)
            return self.unitConverter.fromSIunit_Cv(Region1.Cv1_pT(p, Ts))
        elif Region == 2:
            Ts = Region2.T2_ph(p, h)
            return self.unitConverter.fromSIunit_Cv(Region2.Cv2_pT(p, Ts))
        elif Region == 3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self.unitConverter.fromSIunit_Cv(Region3.Cv3_rhoT(rhos, Ts))
        elif Region == 4:
            return float("NaN")
        elif Region == 5:
            Ts = Region5.T5_ph(p, h)
            return self.unitConverter.fromSIunit_Cv(Region5.Cv5_pT(p, Ts))
        else:
            return float("NaN")

    def Cv_ps(self, p, s):
        '''Specific isochoric heat capacity as a function of pressure and entropy'''
        # case 'cv_ps'
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        Region = RegionSelection.region_ps(p, s)
        if Region == 1:
            Ts = Region1.T1_ps(p, s)
            return self.unitConverter.fromSIunit_Cv(Region1.Cv1_pT(p, Ts))
        elif Region == 2:
            Ts = Region2.T2_ps(p, s)
            return self.unitConverter.fromSIunit_Cv(Region2.Cv2_pT(p, Ts))
        elif Region == 3:
            rhos = 1 / Region3.v3_ps(p, s)
            Ts = Region3.T3_ps(p, s)
            return self.unitConverter.fromSIunit_Cv(Region3.Cv3_rhoT(rhos, Ts))
        elif Region == 4:
            return float("NaN")  # % (xs * CvVp + (1 - xs) * CvLp) / Cv_scale - Cv_offset
        elif Region == 5:
            Ts = Region5.T5_ps(p, s)
            return self.unitConverter.fromSIunit_Cv(Region5.Cv5_pT(p, Ts))
        else:
            return float("NaN")

#%***********************************************************************************************************
# %*1.11 Speed of sound
    def wV_p(self, p):
        '''Saturated vapour speed of sound as a function of preasure'''
        # case 'wv_p'
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_w(Region2.w2_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_w(Region3.w3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p)))
        else:
            return float("NaN")

    def wL_p(self, p):
        '''Saturated liquid speed of sound as a function of preasure'''
        # case 'wl_p'
        p = self.unitConverter.toSIunit_p(p)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                return self.unitConverter.fromSIunit_w(Region1.w1_pT(p, Region4.T4_p(p)))
            else:
                return self.unitConverter.fromSIunit_w(Region3.w3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p)))
        else:
            return float("NaN")

    def wV_t(self, t):
        '''Saturated vapour speed of sound as a function of temperature'''
        # case 'wv_t'
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if (T <= 623.15):
                return self.unitConverter.fromSIunit_w(Region2.w2_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_w(Region3.w3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4V_p(Region4.p4_T(T)))), T))
        else:
            return float("NaN")

    def wL_t(self, t):
        '''Saturated liquid speed of sound as a function of temperature'''
        # case 'wl_t'
        T = self.unitConverter.toSIunit_T(t)
        if (T > 273.15) and (T < 647.096):
            if (T <= 623.15):
                return self.unitConverter.fromSIunit_w(Region1.w1_pT(Region4.p4_T(T), T))
            else:
                return self.unitConverter.fromSIunit_w(Region3.w3_rhoT(1 / (Region3.v3_ph(Region4.p4_T(T), Region4.h4L_p(Region4.p4_T(T)))), T))
        else:
            return float("NaN")

    def w_pt(self, p, t):
        '''Speed of sound as a function of pressure and temperature'''
        # case 'w_pt'
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        Region = RegionSelection.region_pT(p, T)
        if Region == 1:
            return self.unitConverter.fromSIunit_w(Region1.w1_pT(p, T))
        elif Region == 2:
            return self.unitConverter.fromSIunit_w(Region2.w2_pT(p, T))
        elif Region == 3:
            hs = Region3.h3_pT(p, T)
            rhos = 1 / Region3.v3_ph(p, hs)
            return self.unitConverter.fromSIunit_w(Region3.w3_rhoT(rhos, T))
        elif Region == 4:
            return float("NaN")
        elif Region == 5:
            return self.unitConverter.fromSIunit_w(Region5.w5_pT(p, T))
        else:
            return float("NaN")

    def w_ph(self, p, h):
        '''Speed of sound as a function of pressure and enthalpy'''
        # case 'w_ph'
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        Region = RegionSelection.region_ph(p, h)
        if Region == 1:
            Ts = Region1.T1_ph(p, h)
            return self.unitConverter.fromSIunit_w(Region1.w1_pT(p, Ts))
        elif Region == 2:
            Ts = Region2.T2_ph(p, h)
            return self.unitConverter.fromSIunit_w(Region2.w2_pT(p, Ts))
        elif Region == 3:
            rhos = 1 / Region3.v3_ph(p, h)
            Ts = Region3.T3_ph(p, h)
            return self.unitConverter.fromSIunit_w(Region3.w3_rhoT(rhos, Ts))
        elif Region == 4:
            return float("NaN")
        elif Region == 5:
            Ts = Region5.T5_ph(p, h)
            return self.unitConverter.fromSIunit_w(Region5.w5_pT(p, Ts))
        else:
            return float("NaN")

    def w_ps(self, p, s):
        '''Speed of sound as a function of pressure and entropy'''
        # case 'w_ps'
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        Region = RegionSelection.region_ps(p, s)
        if Region == 1:
            Ts = Region1.T1_ps(p, s)
            return self.unitConverter.fromSIunit_w(Region1.w1_pT(p, Ts))
        elif Region == 2:
            Ts = Region2.T2_ps(p, s)
            return self.unitConverter.fromSIunit_w(Region2.w2_pT(p, Ts))
        elif Region == 3:
            rhos = 1 / Region3.v3_ps(p, s)
            Ts = Region3.T3_ps(p, s)
            return self.unitConverter.fromSIunit_w(Region3.w3_rhoT(rhos, Ts))
        elif Region == 4:
            return float("NaN")  # % (xs * wVp + (1 - xs) * wLp) / w_scale - w_offset
        elif Region == 5:
            Ts = Region5.T5_ps(p, s)
            return self.unitConverter.fromSIunit_w(Region5.w5_pT(p, Ts))
        else:
            return float("NaN")

#%***********************************************************************************************************
# %* 1.12 Viscosity
# %Viscosity is not part of IAPWS Steam IF97. Equations from
# %"Revised Release on the IAPWS Formulation 1985 for the Viscosity of Ordinary Water Substance", 2003 are used.
# %Viscosity in the mixed region (4) is interpolated according to the density. This is not true since it will be two fases.
    def my_pt(self, p, t):
        '''Viscosity as a function of pressure and temperature'''
        # case 'my_pt'
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        Region = RegionSelection.region_pT(p, T)
        if Region == 4:
            return float("NaN")
        elif Region in [1, 2, 3, 5]:
            return self.unitConverter.fromSIunit_my(TransportProperties.my_AllRegions_pT(p, T))
        else:
            return float("NaN")

    def my_ph(self, p, h):
        '''Viscosity as a function of pressure and enthalpy'''
        # case {'my_ph'}
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        Region = RegionSelection.region_ph(p, h)
        if Region in [1, 2, 3, 5]:
            return self.unitConverter.fromSIunit_my(TransportProperties.my_AllRegions_ph(p, h))
        elif Region == 4:
            return float("NaN")
        else:
            return float("NaN")

    def my_ps(self, p, s):
        '''Viscosity as a function of pressure and entropy'''
        # case 'my_ps'
        h = self.h_ps(p, s)
        return self.my_ph(p, h)

#%***********************************************************************************************************
# %*1.13 Prandtl
    def pr_pt(self, p, t):
        # case 'pr_pt'
        Cp = self.unitConverter.toSIunit_Cp(self.Cp_pt(p, t))
        my = self.unitConverter.toSIunit_my(self.my_pt(p, t))
        tc = self.unitConverter.toSIunit_tc(self.tc_pt(p, t))
        return Cp * 1000 * my / tc

    def pr_ph(self, p, h):
        # case 'pr_ph'
        Cp = self.unitConverter.toSIunit_Cp(self.Cp_ph(p, h))
        my = self.unitConverter.toSIunit_my(self.my_ph(p, h))
        tc = self.unitConverter.toSIunit_tc(self.tc_ph(p, h))
        return Cp * 1000 * my / tc

#%***********************************************************************************************************
# %*1.15 Surface tension
    def st_t(self, t):
        '''Surface tension for two phase water/steam as a function of T'''
        # case 'st_t'
        T = self.unitConverter.toSIunit_T(t)
        return self.unitConverter.fromSIunit_st(TransportProperties.Surface_Tension_T(T))

    def st_p(self, p):
        '''Surface tension for two phase water/steam as a function of T'''
        # case 'st_p'
        T = self.tsat_p(p)
        T = self.unitConverter.toSIunit_T(T)
        return self.unitConverter.fromSIunit_st(TransportProperties.Surface_Tension_T(T))

#%***********************************************************************************************************
# %*1.16 Thermal conductivity
# %Revised release on the IAPS Formulation 1985 for the Thermal Conductivity of ordinary water substance (IAPWS 1998)
    def tcL_p(self, p):
        '''Saturated vapour thermal conductivity as a function of preasure'''
        # case 'tcl_p'
        t = self.tsat_p(p)
        v = self.vL_p(p)
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(t)
        v = self.unitConverter.toSIunit_v(v)
        rho = 1 / v
        return self.unitConverter.fromSIunit_tc(TransportProperties.tc_ptrho(p, T, rho))

    def tcV_p(self, p):
        '''Saturated liquid thermal conductivity as a function of preasure'''
        # case 'tcv_p'
        ps = p;
        T = self.tsat_p(p)
        v = self.vV_p(ps)
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(T)
        v = self.unitConverter.toSIunit_v(v)
        rho = 1 / v
        return self.unitConverter.fromSIunit_tc(TransportProperties.tc_ptrho(p, T, rho))

    def tcL_t(self, t):
        '''Saturated vapour thermal conductivity as a function of temperature'''
        # case 'tcl_t'
        Ts = t
        p = self.psat_t(Ts)  # XSteam('psat_T', Ts)
        v = self.vL_t(Ts)  # XSteam('vL_T', Ts)
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(Ts)
        v = self.unitConverter.toSIunit_v(v)
        rho = 1 / v
        return self.unitConverter.fromSIunit_tc(TransportProperties.tc_ptrho(p, T, rho))

    def tcV_t(self, t):
        '''Saturated liquid thermal conductivity as a function of temperature'''
        # case 'tcv_t'
        Ts = t
        p = self.psat_t(Ts)  # XSteam('psat_T', Ts)
        v = self.vV_t(Ts)  # XSteam('vV_T', Ts)
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(Ts)
        v = self.unitConverter.toSIunit_v(v)
        rho = 1 / v
        return self.unitConverter.fromSIunit_tc(TransportProperties.tc_ptrho(p, T, rho))

    def tc_pt(self, p, t):
        '''Thermal conductivity as a function of pressure and temperature'''
        # case 'tc_pt'
        Ts = t
        ps = p
        v = self.v_pt(ps, Ts)  # XSteam('v_pT', ps, Ts)
        p = self.unitConverter.toSIunit_p(ps)
        T = self.unitConverter.toSIunit_T(Ts)
        v = self.unitConverter.toSIunit_v(v)
        rho = 1 / v
        return self.unitConverter.fromSIunit_tc(TransportProperties.tc_ptrho(p, T, rho))

    def tc_ph(self, p, h):
        '''Thermal conductivity as a function of pressure and enthalpy'''
        # case 'tc_ph'
        hs = h
        ps = p
        v = self.v_ph(ps, hs)  # XSteam('v_ph', ps, hs)
        T = self.t_ph(ps, hs)  # XSteam('T_ph', ps, hs)
        p = self.unitConverter.toSIunit_p(ps)
        T = self.unitConverter.toSIunit_T(T)
        v = self.unitConverter.toSIunit_v(v)
        rho = 1 / v
        return self.unitConverter.fromSIunit_tc(TransportProperties.tc_ptrho(p, T, rho))

    def tc_hs(self, h, s):
        '''Thermal conductivity as a function of enthalpy and entropy'''
        # case 'tc_hs'
        hs = h
        p = self.p_hs(hs, s)  # XSteam('p_hs', hs, s)
        ps = p
        v = self.v_ph(ps, hs)  # XSteam('v_ph', ps, hs)
        T = self.t_ph(ps, hs)  # XSteam('T_ph', ps, hs)
        p = self.unitConverter.toSIunit_p(p)
        T = self.unitConverter.toSIunit_T(T)
        v = self.unitConverter.toSIunit_v(v)
        rho = 1 / v
        return self.unitConverter.fromSIunit_tc(TransportProperties.tc_ptrho(p, T, rho))

#%***********************************************************************************************************
# %*1.17 Vapour fraction
    def x_ph(self, p, h):
        '''Vapour fraction as a function of pressure and enthalpy'''
        # case 'x_ph'
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        if (p > 0.000611657) and (p < 22.06395):
            return self.unitConverter.fromSIunit_x(Region4.x4_ph(p, h))
        else:
            return float("NaN")

    def x_ps(self, p, s):
        '''Vapour fraction as a function of pressure and entropy'''
        # case 'x_ps'
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        if (p > 0.000611657) and (p < 22.06395):
            return self.unitConverter.fromSIunit_x(Region4.x4_ps(p, s))
        else:
            return float("NaN")


#%***********************************************************************************************************
# %*1.18 Vapour Volume Fraction
    def vx_ph(self, p, h):
        '''Vapour volume fraction as a function of pressure and enthalpy'''
        # case 'vx_ph'
        p = self.unitConverter.toSIunit_p(p)
        h = self.unitConverter.toSIunit_h(h)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                vL = Region1.v1_pT(p, Region4.T4_p(p))
                vV = Region2.v2_pT(p, Region4.T4_p(p))
            else:
                vL = Region3.v3_ph(p, Region4.h4L_p(p))
                vV = Region3.v3_ph(p, Region4.h4V_p(p))
            xs = Region4.x4_ph(p, h)
            return self.unitConverter.fromSIunit_vx((xs * vV / (xs * vV + (1 - xs) * vL)))
        else:
            return float("NaN")

    def vx_ps(self, p, s):
        '''Vapour volume fraction as a function of pressure and entropy'''
        # case 'vx_ps'
        p = self.unitConverter.toSIunit_p(p)
        s = self.unitConverter.toSIunit_s(s)
        if (p > 0.000611657) and (p < 22.06395):
            if p < 16.529:
                vL = Region1.v1_pT(p, Region4.T4_p(p))
                vV = Region2.v2_pT(p, Region4.T4_p(p))
            else:
                vL = Region3.v3_ph(p, Region4.h4L_p(p))
                vV = Region3.v3_ph(p, Region4.h4V_p(p))
            xs = Region4.x4_ps(p, s)
            return self.unitConverter.fromSIunit_vx((xs * vV / (xs * vV + (1 - xs) * vL)))
        else:
            return float("NaN")

class UnitConverter(object):
    def __init__(self, mksSystem = True):
        self.mksSystem = mksSystem

    def toSIunit_p(self, Ins):
        # function toSIunit_p = toSIunit_p( Ins )
        if self.mksSystem:
            # %Translate bar to MPa
            return float(Ins / 10)
        else:
            # % Translate psi to MPa
            return float(Ins * 0.00689475729)

    def fromSIunit_p(self, Ins):
        # function fromSIunit_p = fromSIunit_p( Ins )
        if self.mksSystem:
            # %Translate bar to MPa
            return float(Ins * 10)
        else:
            # %Translate MPa to psi
            return float(Ins / 0.00689475729)

    def toSIunit_T(self, Ins):
        # function toSIunit_T = toSIunit_T( Ins )
        if self.mksSystem:
            # %Translate degC to Kelvin
            return float(Ins + 273.15)
        else:
            # %Translate degF to Kelvin
            return float((5 / 9) * (Ins - 32) + 273.15)

    def fromSIunit_T(self, Ins):
        # function fromSIunit_T = fromSIunit_T( Ins )
        if self.mksSystem:
            # %Translate Kelvin to degC
            return float(Ins - 273.15)
        else:
            # %Translate Kelvin to degF
            return float((Ins - 273.15) * (9 / 5) + 32)

    def toSIunit_h(self, Ins):
        # function toSIunit_h = toSIunit_h( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %Translate btu/lb to kJ/kg
            return float(2.32600 * Ins)

    def fromSIunit_h(self, Ins):
        # function  fromSIunit_h = fromSIunit_h( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %Translate kJ/kg to btu/lb
            return float(Ins / 2.32600)

    def toSIunit_v(self, Ins):
        # function toSIunit_v = toSIunit_v( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %ft^3/lb to m^3/kg
            return float(Ins * 0.0624279606)

    def fromSIunit_v(self, Ins):
        # function fromSIunit_v = fromSIunit_v( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %m^3/kg to ft^3/lb
            return float(Ins / 0.0624279606)

    def toSIunit_s(self, Ins):
        # function toSIunit_s = toSIunit_s( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %btu/(lb degF) to kJ/(kg degC)
            return float(Ins / 0.238845896627)

    def fromSIunit_s(self, Ins):
        # function fromSIunit_s = fromSIunit_s( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %kJ/(kg degC) to btu/(lb degF)
            return float(Ins * 0.238845896627)

    def toSIunit_u(self, Ins):
        # function toSIunit_u = toSIunit_u( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %Translate btu/lb to kJ/kg
            return float(Ins * 2.32600)

    def fromSIunit_u(self, Ins):
        # function fromSIunit_u = fromSIunit_u( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %Translate kJ/kg to btu/lb
            return float(Ins / 2.32600)

    def toSIunit_Cp(self, Ins):
        # function toSIunit_Cp = toSIunit_Cp( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %btu/(lb degF) to kJ/(kg degC)
            return float(Ins / 0.238846)

    def fromSIunit_Cp(self, Ins):
        # function fromSIunit_Cp = fromSIunit_Cp( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %kJ/(kg degC) to btu/(lb degF)
            return float(Ins * 0.238846)

    def toSIunit_Cv(self, Ins):
        # function toSIunit_Cv = toSIunit_Cv( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %btu/(lb degF) to kJ/(kg degC)
            return float(Ins / 0.238846)

    def fromSIunit_Cv(self, Ins):
        # function fromSIunit_Cv = fromSIunit_Cv( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %kJ/(kg degC) to btu/(lb degF)
            return float(Ins * 0.238846)

    def toSIunit_w(self, Ins):
        # function toSIunit_w = toSIunit_w( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %ft/s to m/s
            return float(Ins * 0.3048)

    def fromSIunit_w(self, Ins):
        # function fromSIunit_w = fromSIunit_w( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %m/s to ft/s
            return float(Ins / 0.3048)

    def toSIunit_tc(self, Ins):
        # function toSIunit_tc = toSIunit_tc( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %btu/(h*ft*degF) to W/(m*degC)
            return float(Ins / 0.577789)

    def fromSIunit_tc(self, Ins):
        # function fromSIunit_tc = fromSIunit_tc( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %W/(m*degC) to btu/(h*ft*degF)
            return float(Ins * 0.577789)

    def toSIunit_st(self, Ins):
        # function toSIunit_st = toSIunit_st( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %lb/ft to N/m
            return float(Ins / 0.068521766)

    def fromSIunit_st(self, Ins):
        # function fromSIunit_st = fromSIunit_st( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %N/m to lb/ft
            return float(Ins * 0.068521766)

    def toSIunit_x(self, Ins):
        # function toSIunit_x = toSIunit_x( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %Fraction 0-1
            return float(Ins)

    def fromSIunit_x(self, Ins):
        # function fromSIunit_x = fromSIunit_x( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %Fraction 0-1
            return float(Ins)

    def toSIunit_vx(self, Ins):
        # function toSIunit_vx = toSIunit_vx( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %Fraction 0-1
            return float(Ins)

    def fromSIunit_vx(self, Ins):
        # function fromSIunit_vx = fromSIunit_vx( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %Fraction 0-1
            return float(Ins)

    def toSIunit_my(self, Ins):
        # function toSIunit_my = toSIunit_my( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # % lbm/ft/hr to PaS (N*s/m^2)
            return float(Ins / 2419.088311)

    def fromSIunit_my(self, Ins):
        # function fromSIunit_my = fromSIunit_my( Ins )
        if self.mksSystem:
            return float(Ins)
        else:
            # %PaS (N*s/m^2) to lbm/ft/hr
            return float(Ins * 2419.088311)
