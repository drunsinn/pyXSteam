# -*- coding: utf-8 -*-

import unittest
from pyXSteam import IAPWS_R14


class R14_FunctionTester(unittest.TestCase):
    """R14-08(2011) computer program verification"""

    def setUp(self):
        self.max_error = 1e-6
        self.max_error_ice_Ih = 0.002
        self.max_error_ice_III = 0.003
        self.max_error_ice_V = 0.003
        self.max_error_ice_VI = 0.0035
        self.max_error_ice_VII = 0.007

    def tearDown(self):
        pass

    def test_pmelt_T_function(self):
        """R14-08(2011) computer program verification for pmelt with values from Table 3"""
        error = abs(IAPWS_R14.pmelt_T_iceIh(260.0) - 138.268)
        self.assertLess(error, self.max_error_ice_Ih, "Error for pmelt_T_iceIh to big")

        error = abs(IAPWS_R14.pmelt_T_iceIII(254.0) - 268.685)
        self.assertLess(error, self.max_error_ice_III, "Error for pmelt_T_iceIII to big")

        error = abs(IAPWS_R14.pmelt_T_iceV(265.0) - 479.640)
        self.assertLess(error, self.max_error_ice_V, "Error for pmelt_t_iceV to big")

        error = abs(IAPWS_R14.pmelt_T_iceVI(320.0) - 1356.76)
        self.assertLess(error, self.max_error_ice_VI, "Error for pmelt_t_iceVI to big")

        error = abs(IAPWS_R14.pmelt_T_iceVII(550.0) - 6308.71)
        self.assertLess(error, self.max_error_ice_VII, "Error for pmelt_t_iceVII to big")

    def test_R14_psubl_T_function(self):
        """R14-08(2011) computer program verification for psubl with values from Table 3"""
        error = abs(IAPWS_R14.psubl_T(230.0) - 8.94735e-6)
        self.assertLess(error, self.max_error, "Error for psubl_t to big")

    def test_pmelt_T_function_custom(self):
        """R14-08(2011) computer program verification for pmelt with custom values"""
        error = abs(IAPWS_R14.pmelt_T_iceIh(251.165) - 208.566)
        self.assertLess(error, self.max_error_ice_Ih, "Error for pmelt_T_iceIh to big")

        error = abs(IAPWS_R14.pmelt_T_iceIII(251.165) - 208.566)
        self.assertLess(error, self.max_error_ice_III, "Error for pmelt_T_iceIII to big")

        error = abs(IAPWS_R14.pmelt_T_iceV(256.164) - 350.1)
        self.assertLess(error, self.max_error_ice_V, "Error for pmelt_t_iceV to big")

        error = abs(IAPWS_R14.pmelt_T_iceVI(273.31) - 632.4)
        self.assertLess(error, self.max_error_ice_VI, "Error for pmelt_t_iceVI to big")

        error = abs(IAPWS_R14.pmelt_T_iceVII(355.0) - 2216)
        self.assertLess(error, self.max_error_ice_VII, "Error for pmelt_t_iceVII to big")
