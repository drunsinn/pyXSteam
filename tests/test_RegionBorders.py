# -*- coding: utf-8 -*-

import unittest
from pyXSteam.RegionBorders import pB23_T, TB23_p, pB2bc_h, hB2bc_p


class BordersTester(unittest.TestCase):
    def setUp(self):
        self.max_error = 1e-6

    def tearDown(self):
        pass

    def test_R7_97_B23(self):
        """R7-97(2012): computer-program verification for boundary between regions 2 and 3, Eqs. (5) and (6)"""
        pB23_T_error = abs(pB23_T(0.623150000e3) - 0.165291643e2)
        self.assertLess(pB23_T_error, self.max_error, "Error for p_B23(T) to big")

        TB23_p_error = abs(TB23_p(0.165291643e2) - 0.623150000e3)
        self.assertLess(TB23_p_error, self.max_error, "Error for T_B23(p) to big")

    def test_R7_97_B2bc(self):
        """R7-97(2012): computer-program verification for boundary between sub regions 2b and 2c, Eqs. (20) and (21)"""
        pB2bc_T_error = abs(pB2bc_h(0.3516004323e4) - 0.100000000e3)
        self.assertLess(pB2bc_T_error, self.max_error, "Error for pB2bc(T) to big")

        hB2bc_p_error = abs(hB2bc_p(0.100000000e3) - 0.3516004323e4)
        self.assertLess(hB2bc_p_error, self.max_error, "Error for hB2bc(p) to big")
