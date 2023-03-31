# -*- coding: utf-8 -*-

import unittest
from pyXSteam import XSteam_HW


class R5_FunctionTester(unittest.TestCase):
    def setUp(self):
        self.maxError = 1e-6
        self.steamTable = XSteam_HW(XSteam_HW.UNIT_SYSTEM_MKS)

    def tearDown(self):
        pass

    def test_R5(self):
        values = list()
        values.append({"t":3.8, "result":74.93})
        values.append({"t":5, "result":74.76})
        values.append({"t":30, "result":71.09})
        values.append({"t":55, "result":67.06})
        values.append({"t":80, "result":62.67})
        values.append({"t":105, "result":57.96 })
        values.append({"t":130, "result":52.95})
        values.append({"t":155, "result":47.67})
        values.append({"t":180, "result":42.16})
        values.append({"t":205, "result":36.45})
        values.append({"t":230, "result":30.59})
        values.append({"t":255, "result":24.65})
        values.append({"t":280, "result":18.69})
        values.append({"t":305, "result":12.83})
        values.append({"t":330, "result":7.24})
        values.append({"t":355, "result":2.26})
        values.append({"t":370, "result":0.05})

        for value in values:
            calc = self.steamTable.st_t(t=value["t"])
            error = calc - value["result"]
            self.assertLess(
                error,
                self.maxError,
                "st_t (D2O) not passed for values t %(t)f: Error is %(error)e allowed: %(max)e"
                % {"t":value["t"], "error": error, "max": self.maxError},
            )
