# -*- coding: utf-8 -*-
""" """

import unittest
from pyXSteam.XSteam import XSteam


class FixedBugs_Tester(unittest.TestCase):
    def setUp(self):
        self.maxError = 1e-6
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)

    def tearDown(self):
        pass

    def test_inf_loop_h_pt(self):
        # Bug: calling XSteam('h_pt',800,554) creates an infinite loop.
        # reported by Miika Wallius on 3 Jul 2020
        # https://www.mathworks.com/matlabcentral/fileexchange/9817-x-steam-thermodynamic-properties-of-water-and-steam

        error = self.steamTable.h_pt(800.0, 554.0) - 2733.6986817
        self.assertLess(
            error,
            self.maxError,
            "h_pt not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError},
        )

    def test_inf_loop_rho_pt(self):
        # Bug: pyXSteam hangs
        # reported by @annhak as issue #10
        # https://github.com/drunsinn/pyXSteam/issues/10

        error = self.steamTable.rho_pt(167.11, 351.73) - 113.6761489
        self.assertLess(
            error,
            self.maxError,
            "rho_pt not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError},
        )

    def test_missing_unit_conversion_h_xx(self):
        # Bug: missing unit conversion in h_px for FLS units
        # reported by Grayson Gall on 29 Jan 2022
        # https://www.mathworks.com/matlabcentral/fileexchange/9817-x-steam-thermodynamic-properties-of-water-and-steam

        self.steamTable._unit_converter.set_unitSystem(XSteam.UNIT_SYSTEM_FLS)

        error = self.steamTable.h_px(160, 0.5) - 765.8069264
        self.assertLess(
            error,
            self.maxError,
            " not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError},
        )

        error = self.steamTable.h_tx(160, 0.5) - 628.905695
        self.assertLess(
            error,
            self.maxError,
            " not passed Error %(error)e allowed: %(max)e"
            % {"error": error, "max": self.maxError},
        )
