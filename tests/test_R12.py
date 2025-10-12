# -*- coding: utf-8 -*-

import unittest
import numpy

from pyXSteam import IAPWS_R12
from pyXSteam.tables import R12_08


class R12_FunctionTester(unittest.TestCase):
    def setUp(self):
        self.max_error = 1e-7

    def tearDown(self):
        pass

    def test_R12_industrial(self):

        res = numpy.zeros(len(R12_08.Table4_my))
        for i, (T, rho) in enumerate(zip(R12_08.Table4_T, R12_08.Table4_rho)):
            res[i] = IAPWS_R12.my_rhoT(rho, T, industrial_application=True)

        error = numpy.sum(numpy.absolute((res - R12_08.Table4_my) / R12_08.Table4_my))
        self.assertLess(error, self.max_error, "Test of simplifyed viscosity function for rho and T in R12 failed")

    # def test_R12(self):
    #     values = list()
    #     values.append({"T": 647.35, "rho": 122, "result": 25.520677e-6})
    #     values.append({"T": 647.35, "rho": 222, "result": 31.337589e-6})
    #     # values.append({"T": 647.35, "rho": 272, "result": 36.228143e-6})
    #     # values.append({"T": 647.35, "rho": 322, "result": 42.961579e-6})
    #     # values.append({"T": 647.35, "rho": 372, "result": 45.688204e-6})
    #     values.append({"T": 647.35, "rho": 422, "result": 49.436256e-6})

    #     for value in values:
    #         calc = IAPWS_R12.eq10(T=value["T"], rho=value["rho"], industrial=False)
    #         error = calc - value["result"]
    #         self.assertLess(
    #             error,
    #             self.max_error,
    #             "ep10 not passed for values T %(T)f rho %(rho)f: Error is %(error)e allowed: %(max)e"
    #             % {
    #                 "T": value["T"],
    #                 "rho": value["rho"],
    #                 "error": error,
    #                 "max": self.max_error,
    #             },
    #         )

    def test_R12_internal(self):

        res = numpy.zeros(len(R12_08.Table5_xi))
        for i, (T, rho) in enumerate(zip(R12_08.Table5_T, R12_08.Table5_rho)):
            print(IAPWS_R12.R12_xi(rho, T))
            res[i] = IAPWS_R12.R12_xi(rho, T)

        error = numpy.sum(numpy.absolute((res - R12_08.Table5_xi) / R12_08.Table5_xi))
        self.assertLess(error, self.max_error, "Test of internal calculation of correlation length ξ failed")

        # in_T = [647.35] * 6  # K
        # in_rho = [
        #     122.0,
        #     222.0,
        #     272.0,
        #     322.0,
        #     372.0,
        #     422.0,
        # ]  # kg / m^3

        # ref_xi = [
        #     0.309247,
        #     1.571405,
        #     5.266522,
        #     16.590209,
        #     5.603768,
        #     1.876244,
        # ]  # nm

        # res = numpy.zeros(len(ref_xi))
        # for i, (rho, T) in enumerate(zip(in_rho, in_T)):
        #     res[i] = IAPWS_R12.R12_xi(rho, T)

        # error = numpy.sum(numpy.absolute((res - ref_xi) / ref_xi))
        # self.assertLess(error, self.max_error * 2, "Test of internal function R12_xi(rho,T) Function failed")

        # ref_my_dash_2 = [
        #     1.00000289,
        #     1.00375120,
        #     1.03416789,
        #     1.09190440,
        #     1.03665871,
        #     1.00596332,
        # ]

        # res = numpy.zeros(len(ref_my_dash_2))
        # for i, (rho, T) in enumerate(zip(in_rho, in_T)):
        #     res[i] = IAPWS_R12.R12_my_dash_2(rho, T)

        # error = numpy.sum(numpy.absolute((res - ref_my_dash_2) / ref_my_dash_2))
        # self.assertLess(error, self.max_error * 2, "Test of internal function R12_my_dash_2(rho,T) Function failed")
