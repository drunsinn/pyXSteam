import math
import unittest

from pyXSteam import IceType, XSteam


class XSteamPathCoverageTester(unittest.TestCase):
    def setUp(self):
        self.steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)

    def test_standard_constants(self):
        self.assertAlmostEqual(self.steam_table.specificGasConstant(), 0.461526)
        self.assertAlmostEqual(self.steam_table.criticalTemperatur(), 373.946)
        self.assertAlmostEqual(self.steam_table.criticalPressure(), 220.6395)
        self.assertAlmostEqual(self.steam_table.triplePointTemperatur(), 0.01)
        self.assertAlmostEqual(self.steam_table.triplePointPressure(), 0.00611657)
        self.assertAlmostEqual(self.steam_table.zeroPointTemperature(), -273.15)

    def test_region_five_public_properties(self):
        # IAPWS R7-97 Table 42 state: T=1500 K, p=0.5 MPa.
        pressure_bar = 5.0
        temperature_celsius = 1226.85
        enthalpy = 5219.76855
        entropy = 9.65408875

        self.assertAlmostEqual(
            self.steam_table.h_pt(pressure_bar, temperature_celsius),
            5219.76855,
            delta=0.01,
        )
        self.assertAlmostEqual(
            self.steam_table.s_pt(pressure_bar, temperature_celsius),
            entropy,
            delta=1e-5,
        )
        self.assertAlmostEqual(
            self.steam_table.v_pt(pressure_bar, temperature_celsius),
            1.38455090,
            delta=1e-5,
        )
        self.assertAlmostEqual(
            self.steam_table.t_ph(pressure_bar, enthalpy),
            temperature_celsius,
            delta=0.01,
        )
        self.assertAlmostEqual(
            self.steam_table.t_ps(pressure_bar, entropy),
            temperature_celsius,
            delta=0.01,
        )
        self.assertAlmostEqual(
            self.steam_table.h_ps(pressure_bar, entropy),
            enthalpy,
            delta=0.02,
        )
        self.assertTrue(math.isnan(self.steam_table.p_hs(enthalpy, entropy)))
        self.assertTrue(math.isnan(self.steam_table.t_hs(enthalpy, entropy)))

    def test_melting_and_sublimation_curves(self):
        self.assertAlmostEqual(self.steam_table.pmelt_t(-13.15, IceType.Ih), 1382.68, places=2)
        self.assertAlmostEqual(self.steam_table.pmelt_t(-19.15, IceType.III), 2686.85, places=2)
        self.assertAlmostEqual(self.steam_table.pmelt_t(-8.15), 4796.40, places=2)
        self.assertAlmostEqual(self.steam_table.pmelt_t(46.85), 13567.6, places=1)
        self.assertAlmostEqual(self.steam_table.pmelt_t(276.85), 63087.1, places=1)
        self.assertTrue(math.isnan(self.steam_table.pmelt_t(-30.0)))
        self.assertTrue(math.isnan(self.steam_table.pmelt_t(-30.0, IceType.Ih)))
        self.assertTrue(math.isnan(self.steam_table.pmelt_t(-30.0, IceType.III)))
        self.assertTrue(math.isnan(self.steam_table.pmelt_t(20.0, IceType.V)))
        self.assertTrue(math.isnan(self.steam_table.pmelt_t(82.0, IceType.VI)))
        self.assertTrue(math.isnan(self.steam_table.pmelt_t(500.0, IceType.VII)))
        with self.assertRaises(ValueError):
            self.steam_table.pmelt_t(0.0, 2)

        self.assertAlmostEqual(self.steam_table.psubl_t(-43.15), 8.94735e-5, delta=1e-8)
        self.assertTrue(math.isnan(self.steam_table.psubl_t(-300.0)))
