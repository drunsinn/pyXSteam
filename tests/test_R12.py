# -*- coding: utf-8 -*-

import unittest
import numpy

from pyXSteam import IAPWS_R12
from pyXSteam.tables import R12_08
from . import helpers


class R12_FunctionTester(unittest.TestCase):
    def setUp(self):
        self.max_error = 1e-7

    def tearDown(self):
        pass

    def test_R12_industrial(self):

        helpers.array_2d_test(
            IAPWS_R12.my_rhoT,
            (R12_08.Table4_rho, R12_08.Table4_T),
            R12_08.Table4_my,
            self.max_error,
        )

    def test_R12(self):
        self.skipTest("Not implemented yet")

        values = list()
        values.append({"T": 647.35, "rho": 122, "result": 25.520677e-6})
        values.append({"T": 647.35, "rho": 222, "result": 31.337589e-6})
        values.append({"T": 647.35, "rho": 272, "result": 36.228143e-6})
        values.append({"T": 647.35, "rho": 322, "result": 42.961579e-6})
        values.append({"T": 647.35, "rho": 372, "result": 45.688204e-6})
        values.append({"T": 647.35, "rho": 422, "result": 49.436256e-6})

        for value in values:
            calc = IAPWS_R12.eq10(T=value["T"], rho=value["rho"], industrial=False)
            error = calc - value["result"]
            self.assertLess(
                error,
                self.max_error,
                "ep10 not passed for values T %(T)f rho %(rho)f: Error is %(error)e allowed: %(max)e"
                % {
                    "T": value["T"],
                    "rho": value["rho"],
                    "error": error,
                    "max": self.max_error,
                },
            )

    def test_R12_internal(self):

        self.skipTest("Not implemented yet")

        res = numpy.zeros(len(R12_08.Table5_xi))
        for i, (T, rho) in enumerate(zip(R12_08.Table5_T, R12_08.Table5_rho)):
            print(IAPWS_R12.R12_xi(rho, T))
            res[i] = IAPWS_R12.R12_xi(rho, T)

        error = numpy.sum(numpy.absolute((res - R12_08.Table5_xi) / R12_08.Table5_xi))
        self.assertLess(error, self.max_error, "Test of internal calculation of correlation length ξ failed")
