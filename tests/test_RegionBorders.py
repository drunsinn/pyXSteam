# -*- coding: utf-8 -*-

import unittest

from test_Regions import REGION_ERROR
import pyXSteam.RegionBorders as RB


class BordersTester(unittest.TestCase):
    def setUp(self):
        self.max_error = REGION_ERROR

    def tearDown(self):
        pass

    def _calc_err(
        self,
        val,
        ref,
    ):
        return abs(val - ref) / ref

    def test_R7_97_B23(self):
        """R7-97(2012): computer-program verification for boundary between regions 2 and 3, Eqs. (5) and (6)"""
        error = self._calc_err(RB.pB23_T(0.623150000e3), 0.165291643e2)
        self.assertLess(error, self.max_error, "Error for p_B23(T) to big")

        error = self._calc_err(RB.TB23_p(0.165291643e2), 0.623150000e3)
        self.assertLess(error, self.max_error, "Error for T_B23(p) to big")

    def test_R7_97_B2bc(self):
        """R7-97(2012): computer-program verification for boundary between sub regions 2b and 2c, Eqs. (20) and (21)"""
        error = self._calc_err(RB.pB2bc_h(0.3516004323e4), 0.100000000e3)
        self.assertLess(error, self.max_error, "Error for pB2bc(T) to big")

        error = self._calc_err(RB.hB2bc_p(0.100000000e3), 0.3516004323e4)
        self.assertLess(error, self.max_error, "Error for hB2bc(p) to big")

    def test_SR2_01_B2bc(self):
        """SR2-01(2016): computer-program verification for boundary between sub regions 2b and 2c, Eq. (2)"""
        error = self._calc_err(RB.hB2bc_s(7.0), 3376.437884) / 3376.437884
        self.assertLess(error, self.max_error, "Error for hB2bc(s) to big")

    def test_SR4_04_B13(self):
        """SR4-04(2014): computer-program verification for boundary between regions 1 and 3, Eq. (7)"""
        # Table 24
        error = self._calc_err(RB.hB13_s(3.7), 1.632525047e3)
        self.assertLess(error, self.max_error, "Error for hB13(s) to big")

        error = self._calc_err(RB.hB13_s(3.6), 1.593027215e3)
        self.assertLess(error, self.max_error, "Error for hB13(s) to big")

        error = self._calc_err(RB.hB13_s(3.5), 1.566104611e3)
        self.assertLess(error, self.max_error, "Error for hB13(s) to big")

    def test_SR4_04_B23(self):
        """SR4-04(2014): computer-program verification for boundary between regions 2 and 3, Eq. (8)"""
        error = self._calc_err(RB.TB23_hs(2600.0, 5.1), 7.135259364e2)
        self.assertLess(error, self.max_error, "Error for hB13(s) to big")

        error = self._calc_err(RB.TB23_hs(2700.0, 5.15), 7.685345532e2)
        self.assertLess(error, self.max_error, "Error for hB13(s) to big")

        error = self._calc_err(RB.TB23_hs(2800.0, 5.2), 8.176202120e2)
        self.assertLess(error, self.max_error, "Error for hB13(s) to big")
