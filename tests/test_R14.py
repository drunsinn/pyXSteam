# -*- coding: utf-8 -*-

import unittest
from pyXSteam import IAPWS_R14


class R14_FunctionTester(unittest.TestCase):
    def setUp(self):
        self.maxError = 1e-6
        self.maxError_ice_III = 0.003
        self.maxError_ice_V = 0.003
        self.maxError_ice_VI = 0.003
        self.maxError_ice_VII = 0.007
        self.maxError_ice_Ih = 0.002

    def tearDown(self):
        pass

    def test_R14_pmelt_T_function_Ih_1(self):
        error = IAPWS_R14.pmelt_T_iceIh(251.165) - 208.566
        self.assertLess(
            error,
            self.maxError_ice_Ih,
            "pmelt_T_iceIh not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError_ice_Ih},
        )

    def test_R14_pmelt_T_function_Ih_2(self):
        error = IAPWS_R14.pmelt_T_iceIh(254) - 268.685
        self.assertLess(
            error,
            self.maxError_ice_Ih,
            "pmelt_T_iceIh not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError_ice_Ih},
        )

    def test_R14_pmelt_T_function_III_1(self):
        error = IAPWS_R14.pmelt_T_iceIII(251.165) - 208.566
        self.assertLess(
            error,
            self.maxError_ice_III,
            "pmelt_T_iceIII not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError_ice_III},
        )

    def test_R14_pmelt_T_function_III_2(self):
        error = IAPWS_R14.pmelt_T_iceIII(254.0) - 268.685
        self.assertLess(
            error,
            self.maxError_ice_III,
            "pmelt_T_iceIII not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError_ice_III},
        )

    def test_R14_pmelt_T_function_V_1(self):
        error = IAPWS_R14.pmelt_T_iceV(256.164) - 350.1
        self.assertLess(
            error,
            self.maxError_ice_V,
            "pmelt_t not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError_ice_V},
        )

    def test_R14_pmelt_T_function_V_2(self):
        error = IAPWS_R14.pmelt_T_iceV(265) - 479.640
        self.assertLess(
            error,
            self.maxError_ice_V,
            "pmelt_t not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError_ice_V},
        )

    def test_R14_pmelt_T_function_VI_1(self):
        error = IAPWS_R14.pmelt_T_iceVI(273.31) - 632.4
        self.assertLess(
            error,
            self.maxError_ice_VI,
            "pmelt_t not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError_ice_VI},
        )

    def test_R14_pmelt_T_function_VI_2(self):
        error = IAPWS_R14.pmelt_T_iceVI(320) - 1356.76
        self.assertLess(
            error,
            self.maxError_ice_VI,
            "pmelt_t not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError_ice_VI},
        )

    def test_R14_pmelt_T_function_VII_1(self):
        error = IAPWS_R14.pmelt_T_iceVII(355.0) - 2216
        self.assertLess(
            error,
            self.maxError_ice_VII,
            "pmelt_t not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError_ice_VII},
        )

    def test_R14_pmelt_T_function_VII_2(self):
        error = IAPWS_R14.pmelt_T_iceVII(550) - 6308.71
        self.assertLess(
            error,
            self.maxError_ice_VII,
            "pmelt_t not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError_ice_VII},
        )

    def test_R14_psubl_T_function(self):
        error = IAPWS_R14.psubl_T(230.0) - 8.94735e-6
        self.assertLess(
            error,
            self.maxError,
            "psubl_t not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError},
        )
