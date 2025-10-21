# -*- coding: utf-8 -*-

import unittest

import numpy
from pyXSteam import IAPWS_R15


class TransportTester(unittest.TestCase):
    def setUp(self):
        self.max_error = 1e-8
        # self.max_matrix_error = 3e-3

    def tearDown(self):
        pass

    def test_R15_thermal_conductvity(self):
        """R15-11 test correlating equation"""
        self.skipTest("Not implemented yet")

        in_T = [298.15, 298.15, 298.15, 873.15]  # in K
        in_rho = [0.0, 998.0, 1200.0, 0.0]  # in kg/m^3
        ref = [18.4341883, 607.712868, 799.038144, 79.1034659]  # mW / m K
        res = numpy.zeros(len(ref))

        for i, (T, rho) in enumerate(zip(in_T, in_rho)):
            res[i] = IAPWS_R15.R15_11_correlation(T, rho)

        error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(error, self.max_error, "Test of correlation function for T and rho in R15 failed")
