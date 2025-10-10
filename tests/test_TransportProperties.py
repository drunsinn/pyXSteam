# -*- coding: utf-8 -*-

import math
import numpy
import unittest

from pyXSteam.TransportProperties import surface_tension_T
from pyXSteam.Constants import TRIPLE_POINT_TEMPERATURE, CRITICAL_TEMPERATURE


class TransportTester(unittest.TestCase):
    """tests for functions in Transport Properties"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_R1_surface_tension(self):
        """R1-76(2014) test calculation of surface tension from T"""
        # selected values from Table 1
        in_T = [0.011, 30.0, 55.0, 80.0, 105.0, 130.0, 155.0, 180.0, 205.0, 230.0, 255.0, 280.0, 305.0, 330.0, 355.0, 370.0]  # in °C
        in_T = [t + TRIPLE_POINT_TEMPERATURE for t in in_T]  # # quick convert to Klevin
        ref = [75.65, 71.19, 67.10, 62.67, 57.94, 52.93, 47.67, 42.19, 36.53, 30.74, 24.87, 18.99, 13.22, 7.70, 2.74, 0.39]
        res = numpy.zeros(len(ref))

        for i, t in enumerate(in_T):
            res[i] = surface_tension_T(t)

        error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(error, 1e-2, "Test of st(t) failed")

        assert math.isnan(surface_tension_T(TRIPLE_POINT_TEMPERATURE - 0.01))

        assert math.isnan(surface_tension_T(CRITICAL_TEMPERATURE + 0.01))
