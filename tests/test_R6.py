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
        in_T = [
            300.0,
            300.0,
            300.0,
            500.0,
            500.0,
            500.0,
            500.0,
            647.0,
            900.0,
            900.0,
            900.0,
        ]  # in K

        in_rho = [
            0.9965560e3,
            0.1005308e4,
            0.1188202e4,
            0.4350000,
            0.4532000e1,
            0.8380250e3,
            0.1084564e4,
            0.3580000e3,
            0.2410000,
            0.5261500e2,
            0.8707690e3,
        ]  # in kg / m^3

        ref = [
            0.992418352e-1,
            0.200022515e2,
            0.700004704e3,
            0.999679423e-1,
            0.999938125,
            0.100003858e2,
            0.700000405e3,
            0.220384756e2,
            0.100062559,
            0.200000690e2,
            0.700000006e3,
        ]  # p in MPa

        res = numpy.zeros(len(ref))
        for i, (rho, T) in enumerate(zip(in_rho, in_T)):
            res[i] = IAPWS_R6.R6_p_rhoT(rho, T)

        error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(error, self.max_error * 2, "Test of p(rho,T) Function failed")

        ref = [
            0.413018112e1,
            0.406798347e1,
            0.346135580e1,
            0.150817541e1,
            0.166991025e1,
            0.322106219e1,
            0.307437693e1,
            0.618315728e1,
            0.175890657e1,
            0.193510526e1,
            0.266422350e1,
        ]  # cv in kJ / kg K

        res = numpy.zeros(len(ref))
        for i, (rho, T) in enumerate(zip(in_rho, in_T)):
            res[i] = IAPWS_R6.R6_cv_rhoT(rho, T)

        error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(error, self.max_error * 2, "Test of cv(rho,T) Function failed")

        ref = [
            0.150151914e4,
            0.153492501e4,
            0.244357992e4,
            0.548314253e3,
            0.535739001e3,
            0.127128441e4,
            0.241200877e4,
            0.252145078e3,
            0.724027147e3,
            0.698445674e3,
            0.201933608e4,
        ]  # w in m / s

        res = numpy.zeros(len(ref))
        for i, (rho, T) in enumerate(zip(in_rho, in_T)):
            res[i] = IAPWS_R6.R6_w_rhoT(rho, T)

        # FIXME
        error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(error, 5.3e-3, "Test of w(rho,T) Function failed")

        ref = [
            0.393062643,
            0.387405401,
            0.132609616,
            0.794488271e1,
            0.682502725e1,
            0.256690919e1,
            0.203237509e1,
            0.432092307e1,
            0.916653194e1,
            0.659070225e1,
            0.417223802e1,
        ]  # s in kJ / kg K

        res = numpy.zeros(len(ref))
        for i, (rho, T) in enumerate(zip(in_rho, in_T)):
            res[i] = IAPWS_R6.R6_s_rhoT(rho, T)

        error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(error, self.max_error * 2, "Test of s(rho,T) Function failed")
