# -*- coding: utf-8 -*-

import math
import unittest

from pyXSteam import XSteam_HW
from pyXSteam.Constants import TRIPLE_POINT_TEMPERATURE, CRITICAL_TEMPERATURE
from pyXSteam.TransportProperties_HW import surface_tension_T
from pyXSteam.tables import R5_85
from . import helpers


class HWTester(unittest.TestCase):
    """tester for functions specific to heavy water in XSteam_HW"""

    def setUp(self):
        self.max_error = 1e-1
        self.steam_table = XSteam_HW(XSteam_HW.UNIT_SYSTEM_MKS)

    def tearDown(self):
        pass

    def test_R5_surface_tension(self):
        """R5-85(1994) test calculation of surface tension from T"""
        # selected values from Table 1

        helpers.array_1d_test(
            surface_tension_T,
            [t + TRIPLE_POINT_TEMPERATURE for t in R5_85.Table1_t],  # quick convert to Klevin,
            R5_85.Table1_st_calc,
            self.max_error,
        )

        assert math.isnan(surface_tension_T(TRIPLE_POINT_TEMPERATURE - 0.01))

        assert math.isnan(surface_tension_T(CRITICAL_TEMPERATURE + 0.01))
