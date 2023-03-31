# -*- coding: utf-8 -*-

import unittest
from pyXSteam import XSteam


class R1_FunctionTester(unittest.TestCase):
    def setUp(self):
        self.maxError = 1e-6
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)

    def tearDown(self):
        pass

    def test_R1(self):
        values = list()
        values.append({"t":0.01, "result":75.65})
        values.append({"t":30, "result":71.19})
        values.append({"t":55, "result":67.10})
        values.append({"t":80, "result":62.67})
        values.append({"t":105, "result":57.94})
        values.append({"t":130, "result":52.93})
        values.append({"t":155, "result":47.67})
        values.append({"t":180, "result":42.19})
        values.append({"t":205, "result":36.53})
        values.append({"t":230, "result":30.74})
        values.append({"t":255, "result":24.87})
        values.append({"t":280, "result":18.99})
        values.append({"t":305, "result":13.22})
        values.append({"t":330, "result":7.70})
        values.append({"t":355, "result":2.74})
        values.append({"t":370, "result":0.39})

        for value in values:
            calc = self.steamTable.st_t(t=value["t"])
            error = calc - value["result"]
            self.assertLess(
                error,
                self.maxError,
                "st_t not passed for values t %(t)f: Error is %(error)e allowed: %(max)e"
                % {"t":value["t"], "error": error, "max": self.maxError},
            )
