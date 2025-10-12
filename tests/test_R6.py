# -*- coding: utf-8 -*-

import unittest
import numpy
from pyXSteam import IAPWS_R6
from pyXSteam.tables import R6_95


class R6_FunctionTester(unittest.TestCase):

    def setUp(self):
        self.max_error = 1e-8
        # self.max_matrix_error = 3e-3

    def tearDown(self):
        pass

    def test_R6_internal(self):
        """check internal calculations with values from Table 6"""
        delta = 838.025 / R6_95.CRITICAL_DENSITY
        tau = R6_95.CRITICAL_TEMPERATURE / 500.0

        error = abs(IAPWS_R6.eq_phi_o(tau, delta) - 0.204797733e1)
        self.assertLess(error, self.max_error, "Error for _phi_0 to big")

        error = abs(IAPWS_R6.eq_phi_o_delta(delta) - 0.384236747)
        self.assertLess(error, self.max_error, "Error for _phi_0_delta to big")

        error = abs(IAPWS_R6.eq_phi_o_deltadelta(delta) + 0.147637878)
        self.assertLess(error, self.max_error, "Error for _phi_0_delta_delta to big")

        error = abs(IAPWS_R6.eq_phi_o_tau(tau) - 0.904611106e1)
        self.assertLess(error, self.max_error, "Error for _phi_0_tau to big")

        error = abs(IAPWS_R6.eq_phi_o_tautau(tau) + 0.193249185e1)
        self.assertLess(error, self.max_error, "Error for _phi_0_tau_tau to big")

        error = abs(IAPWS_R6.eq_phi_r(tau, delta) + 0.342693206e1)
        self.assertLess(error, self.max_error, "Error for _phi_r to big")

        error = abs(IAPWS_R6.eq_phi_r_delta(tau, delta) + 0.364366650)
        self.assertLess(error, self.max_error, "Error for eq_phi_r_delta to big")

        error = abs(IAPWS_R6.eq_phi_r_deltadelta(tau, delta) - 0.856063701)
        self.assertLess(error, self.max_error, "Error for eq_phi_r_deltadelta to big")

        error = abs(IAPWS_R6.eq_phi_r_tau(tau, delta) + 0.581403435e1)
        self.assertLess(error, self.max_error, "Error for eq_phi_r_tau to big")

        error = abs(IAPWS_R6.eq_phi_r_tautau(tau, delta) + 0.223440737e1)
        self.assertLess(error, self.max_error, "Error for eq_phi_r_tautau to big")

        error = abs(IAPWS_R6.eq_phi_r_deltatau(tau, delta) + 0.112176915e1)
        self.assertLess(error, self.max_error, "Error for eq_phi_r_deltatau to big")

    def test_p_function(self):
        """check functions against values of R6-95 Table 7"""

        res = numpy.zeros(len(R6_95.Table7_p))
        for i, (rho, T) in enumerate(zip(R6_95.Table7_rho, R6_95.Table7_T)):
            res[i] = IAPWS_R6.R6_p_rhoT(rho, T)
        error = numpy.sum(numpy.absolute((res - R6_95.Table7_p) / R6_95.Table7_p))
        self.assertLess(error, self.max_error * 2, "Test of p(rho,T) Function failed")

        res = numpy.zeros(len(R6_95.Table7_cv))
        for i, (rho, T) in enumerate(zip(R6_95.Table7_rho, R6_95.Table7_T)):
            res[i] = IAPWS_R6.R6_cv_rhoT(rho, T)
        error = numpy.sum(numpy.absolute((res - R6_95.Table7_cv) / R6_95.Table7_cv))
        self.assertLess(error, self.max_error * 2, "Test of cv(rho,T) Function failed")

        res = numpy.zeros(len(R6_95.Table7_w))
        for i, (rho, T) in enumerate(zip(R6_95.Table7_rho, R6_95.Table7_T)):
            res[i] = IAPWS_R6.R6_w_rhoT(rho, T)
        # FIXME
        error = numpy.sum(numpy.absolute((res - R6_95.Table7_w) / R6_95.Table7_w))
        self.assertLess(error, 5.3e-3, "Test of w(rho,T) Function failed")

        res = numpy.zeros(len(R6_95.Table7_s))
        for i, (rho, T) in enumerate(zip(R6_95.Table7_rho, R6_95.Table7_T)):
            res[i] = IAPWS_R6.R6_s_rhoT(rho, T)
        error = numpy.sum(numpy.absolute((res - R6_95.Table7_s) / R6_95.Table7_s))
        self.assertLess(error, self.max_error * 2, "Test of s(rho,T) Function failed")
