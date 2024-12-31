# -*- coding: utf-8 -*-

import math
import unittest
import numpy

from pyXSteam import XSteam_HW
from pyXSteam.Constants import TRIPLE_POINT_TEMPERATURE, CRITICAL_TEMPERATURE
from pyXSteam.IAPWS_R5 import surface_tension_T


class HWTester(unittest.TestCase):
    """tester for functions specific to heavy water in XSteam_HW"""

    def setUp(self):
        self.max_error = 1e-6
        self.steam_table = XSteam_HW(XSteam_HW.UNIT_SYSTEM_MKS)

    def tearDown(self):
        pass

    def test_R5_surface_tension(self):
        """R5-85(1994) test calculation of surface tension from T"""
        # selected values from Table 1
        in_T = [3.8, 5, 30, 55, 80, 105, 130, 155, 180, 205, 230, 255, 280, 305, 330, 355, 370]
        in_T = [t + TRIPLE_POINT_TEMPERATURE for t in in_T]  # # quick convert to Klevin
        ref = [74.93, 74.76, 71.09, 67.06, 62.67, 57.96, 52.95, 47.67, 42.16, 36.45, 30.59, 24.65, 18.69, 12.83, 7.24, 2.26, 0.05]
        res = numpy.zeros(len(ref))

        for i, t in enumerate(in_T):
            res[i] = surface_tension_T(t)

        error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(error, 1e-1, "Test of st_D2O(t) failed")

        assert math.isnan(surface_tension_T(TRIPLE_POINT_TEMPERATURE - 0.01))

        assert math.isnan(surface_tension_T(CRITICAL_TEMPERATURE + 0.01))
