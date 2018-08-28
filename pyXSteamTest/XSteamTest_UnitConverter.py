# -*- coding: utf-8 -*-
"""
Test cases for the unit converter used by pyXSteam
"""

import unittest
from pyXSteam.XSteam import XSteam
from pyXSteam.UnitConverter import UnitConverter


class UnitConverter_MKS_FunctionTester(unittest.TestCase):

    def setUp(self):
        self.maxError = 1E-8
        self.uc = UnitConverter(XSteam.UNIT_SYSTEM_MKS)

    def tearDown(self):
        pass

    def test_toSIunit_p_1_MKS(self):
        error = self.uc.toSIunit_p(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of toSIunit_p for MKS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_p_2_MKS(self):
        error = self.uc.toSIunit_p(100.0) - 991.0
        self.assertLess(error, self.maxError, 'Test of toSIunit_p for MKS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_p_1_MKS(self):
        error = self.uc.fromSIunit_p(1.0) - 10.0
        self.assertLess(error, self.maxError, 'Test of fromSIunit_p for MKS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_p_2_MKS(self):
        error = self.uc.fromSIunit_p(100.0) - 1000.0
        self.assertLess(error, self.maxError, 'Test of fromSIunit_p for MKS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_T_1_MKS(self):
        error = self.uc.toSIunit_T(1.0) - 274.15
        self.assertLess(error, self.maxError, 'Test of toSIunit_T for MKS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_T_2_MKS(self):
        error = self.uc.toSIunit_T(1.0) - 274.15
        self.assertLess(error, self.maxError, 'Test of toSIunit_T for MKS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_T_1_MKS(self):
        error = self.uc.fromSIunit_T(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of fromSIunit_T for MKS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_T_2_MKS(self):
        error = self.uc.fromSIunit_T(42.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of fromSIunit_T for MKS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})


class UnitConverter_FLS_FunctionTester(unittest.TestCase):

    def setUp(self):
        self.maxError = 1E-8
        self.uc = UnitConverter(XSteam.UNIT_SYSTEM_FLS)

    def tearDown(self):
        pass

    def test_toSIunit_p_FLS(self):
        error = self.uc.toSIunit_p(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of toSIunit_p for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_p_FLS(self):
        error = self.uc.fromSIunit_p(1.0) - 145.0377377968587
        self.assertLess(error, self.maxError, 'Test of fromSIunit_p for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_T_FLS(self):
        error = self.uc.toSIunit_T(1.0) - 274.15
        self.assertLess(error, self.maxError, 'Test of toSIunit_T for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_T_FLS(self):
        error = self.uc.fromSIunit_T(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of fromSIunit_T for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_h_FLS(self):
        error = self.uc.toSIunit_h(1.0) - 2.326
        self.assertLess(error, self.maxError, 'Test of toSIunit_h for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_h_FLS(self):
        error = self.uc.fromSIunit_h(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of fromSIunit_h for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_v_FLS(self):
        error = self.uc.toSIunit_v(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of toSIunit_v for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_v_FLS(self):
        error = self.uc.fromSIunit_v(1.0) - 16.018463367839058
        self.assertLess(error, self.maxError, 'Test of fromSIunit_v for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_s_FLS(self):
        error = self.uc.toSIunit_s(1.0) - 4.1868000000086933
        self.assertLess(error, self.maxError, 'Test of toSIunit_s for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_s_FLS(self):
        error = self.uc.fromSIunit_s(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of fromSIunit_s for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_u_FLS(self):
        error = self.uc.toSIunit_u(1.0) - 2.326
        self.assertLess(error, self.maxError, 'Test of toSIunit_u for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_u_FLS(self):
        error = self.uc.fromSIunit_u(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of fromSIunit_u for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_Cp_FLS(self):
        error = self.uc.toSIunit_Cp(1.0) - 4.1867981879537446
        self.assertLess(error, self.maxError, 'Test of toSIunit_Cp for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_Cp_FLS(self):
        error = self.uc.fromSIunit_Cp(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of fromSIunit_Cp for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_Cv_FLS(self):
        error = self.uc.toSIunit_Cv(1.0) - 4.1867981879537446
        self.assertLess(error, self.maxError, 'Test of toSIunit_Cv for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_Cv_FLS(self):
        error = self.uc.fromSIunit_Cv(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of fromSIunit_Cv for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_w_FLS(self):
        error = self.uc.toSIunit_w(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of toSIunit_w for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_w_FLS(self):
        error = self.uc.fromSIunit_w(1.0) - 3.280839895013123
        self.assertLess(error, self.maxError, 'Test of fromSIunit_w for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_tc_FLS(self):
        error = self.uc.toSIunit_tc(1.0) - 1.7307356145582558
        self.assertLess(error, self.maxError, 'Test of toSIunit_tc for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_tc_FLS(self):
        error = self.uc.fromSIunit_tc(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of fromSIunit_tc for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_st_FLS(self):
        error = self.uc.toSIunit_st(1.0) - 14.593902906705587
        self.assertLess(error, self.maxError, 'Test of toSIunit_st for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_st_FLS(self):
        error = self.uc.fromSIunit_st(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of fromSIunit_st for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_x_FLS(self):
        error = self.uc.toSIunit_x(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of toSIunit_x for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_x_FLS(self):
        error = self.uc.fromSIunit_x(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of fromSIunit_x for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_vx_FLS(self):
        error = self.uc.toSIunit_vx(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of toSIunit_vx for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_vx_FLS(self):
        error = self.uc.fromSIunit_vx(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of fromSIunit_vx for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_toSIunit_my_FLS(self):
        error = self.uc.toSIunit_my(1.0) - 1.0
        self.assertLess(error, self.maxError, 'Test of toSIunit_my for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_fromSIunit_my_FLS(self):
        error = self.uc.fromSIunit_my(1.0) - 2419.088311
        self.assertLess(error, self.maxError, 'Test of fromSIunit_my for FLS failed. Error was %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})
