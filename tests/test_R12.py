# -*- coding: utf-8 -*-

import unittest
from pyXSteam import IAPWS_R12


class R12_FunctionTester(unittest.TestCase):
    def setUp(self):
        self.maxError = 1e-6

    def tearDown(self):
        pass

    def test_R12_industrial(self):
        values = list()
        values.append({"T":298.15, "rho":998, "result":889.735100E-6})
        values.append({"T":298.15, "rho":1200, "result":1437.649467E-6})
        values.append({"T":373.15, "rho":1000, "result":307.883622E-6})
        values.append({"T":433.15, "rho":1, "result":14.538324E-6})
        values.append({"T":433.15, "rho":1000, "result":217.685358E-6})
        values.append({"T":873.15, "rho":1, "result":32.619287E-6})
        values.append({"T":873.15, "rho":100, "result":35.802262E-6})
        values.append({"T":873.15, "rho":600, "result":77.430195E-6})
        values.append({"T":1173.15, "rho":1, "result":44.217245E-6})
        values.append({"T":1173.15, "rho":100, "result":47.640433E-6})
        values.append({"T":1173.15, "rho":400, "result":64.154608E-6})

        for value in values:
            calc = IAPWS_R12.eq10(T=value["T"], rho=value["rho"], industrial=True)
            error = calc - value["result"]
            self.assertLess(
                error,
                self.maxError,
                "ep10 not passed for values T %(T)f rho %(rho)f: Error is %(error)e allowed: %(max)e"
                % {"T":value["T"], "rho":value["rho"], "error": error, "max": self.maxError},
            )

    def test_R12(self):
        values = list()
        values.append({"T":647.35, "rho":122, "result":25.520677E-6})
        values.append({"T":647.35, "rho":222, "result":31.337589E-6})
        #values.append({"T":647.35, "rho":272, "result":36.228143E-6})
        #values.append({"T":647.35, "rho":322, "result":42.961579E-6})
        #values.append({"T":647.35, "rho":372, "result":45.688204E-6})
        values.append({"T":647.35, "rho":422, "result":49.436256E-6})

        for value in values:
            calc = IAPWS_R12.eq10(T=value["T"], rho=value["rho"], industrial=False)
            error = calc - value["result"]
            self.assertLess(
                error,
                self.maxError,
                "ep10 not passed for values T %(T)f rho %(rho)f: Error is %(error)e allowed: %(max)e"
                % {"T":value["T"], "rho":value["rho"], "error": error, "max": self.maxError},
            )