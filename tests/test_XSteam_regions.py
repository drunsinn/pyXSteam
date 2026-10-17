import math
import unittest

from pyXSteam import XSteam


class XSteamRegionRoutingTester(unittest.TestCase):
    def setUp(self):
        self.steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)

    def test_single_phase_pt_verification_tables(self):
        # IAPWS IF97 Tables 5 and 15; pressure is converted from MPa to bar.
        states = (
            # region, pressure [bar], temperature [C], v, h, u, s, cp, w
            ("R1 Table 5", 30.0, 26.85, 0.00100215168, 115.331273, 112.324818, 0.392294792, 4.17301218, 1507.73921),
            ("R1 Table 5", 800.0, 26.85, 0.000971180894, 184.142828, 106.448356, 0.368563852, 4.01008987, 1634.69054),
            ("R1 Table 5", 30.0, 226.85, 0.00120241800, 975.542239, 971.934985, 2.58041912, 4.65580682, 1240.71337),
            ("R2 Table 15", 0.035, 26.85, 39.4913866, 2549.91145, 2411.69160, 8.52238967, 1.91300162, 427.920172),
            ("R2 Table 15", 0.035, 426.85, 92.3015898, 3335.68375, 3012.62819, 10.1749996, 2.08141274, 644.289068),
            ("R2 Table 15", 300.0, 426.85, 0.00542946619, 2631.49474, 2468.61076, 5.17540298, 10.3505092, 480.386523),
            # IAPWS IF97 Table 42, the vector matching the current Region 5 implementation.
            ("R5 Table 42", 5.0, 1226.85, 1.38455090, 5219.76855, 4527.49310, 9.65408875, 2.61609445, 917.068690),
        )
        tolerances = (1e-7, 0.01, 0.01, 1e-6, 1e-5, 0.01)
        region5_tolerances = (5e-6, 0.01, 0.01, 1e-5, 1e-5, 0.01)
        methods = ("v_pt", "h_pt", "u_pt", "s_pt", "Cp_pt", "w_pt")

        for case in states:
            region, pressure, temperature, *expected = case
            with self.subTest(region=region, pressure=pressure, temperature=temperature):
                case_tolerances = region5_tolerances if region.startswith("R5") else tolerances
                for method, value, tolerance in zip(methods, expected, case_tolerances):
                    with self.subTest(property=method):
                        self.assertAlmostEqual(
                            getattr(self.steam_table, method)(pressure, temperature),
                            value,
                            delta=tolerance,
                        )

    def test_region3_ph_and_ps_verification_tables(self):
        # IAPWS Supplementary Release SR3-03 (2014), Tables 5, 8, 12, and 15.
        ph_states = (
            (20.0, 1700.0, 629.3083892, 0.001749903962),
            (50.0, 2000.0, 690.5718338, 0.001908139035),
            (100.0, 2100.0, 733.6163014, 0.001676229776),
            (20.0, 2500.0, 641.8418053, 0.006670547043),
            (50.0, 2400.0, 735.1848618, 0.002801244590),
            (100.0, 2700.0, 842.0460876, 0.002404234998),
        )
        for pressure_mpa, enthalpy, temperature_k, volume in ph_states:
            pressure_bar = pressure_mpa * 10
            with self.subTest(pressure_mpa=pressure_mpa, enthalpy=enthalpy):
                self.assertAlmostEqual(
                    self.steam_table.t_ph(pressure_bar, enthalpy),
                    temperature_k - 273.15,
                    delta=0.001,
                )
                self.assertAlmostEqual(
                    self.steam_table.v_ph(pressure_bar, enthalpy),
                    volume,
                    delta=1e-11,
                )

        ps_states = (
            (20.0, 3.8, 628.2959869, 0.001733791463),
            (50.0, 3.6, 629.7158726, 0.001469680170),
            (100.0, 4.0, 705.6880237, 0.001555893131),
            (20.0, 5.0, 640.1176443, 0.006262101987),
            (50.0, 4.5, 716.3687517, 0.002332634294),
            (100.0, 5.0, 847.4332825, 0.002449610757),
        )
        for pressure_mpa, entropy, temperature_k, volume in ps_states:
            pressure_bar = pressure_mpa * 10
            with self.subTest(pressure_mpa=pressure_mpa, entropy=entropy):
                self.assertAlmostEqual(
                    self.steam_table.t_ps(pressure_bar, entropy),
                    temperature_k - 273.15,
                    delta=0.001,
                )
                self.assertAlmostEqual(
                    self.steam_table.v_ps(pressure_bar, entropy),
                    volume,
                    delta=1e-11,
                )

    def test_region3_hs_routes_against_verification_table(self):
        # IAPWS Supplementary Release SR4-04 (2014), Table 5.
        states = (
            (1700.0, 3.8, 25.55703246),
            (2000.0, 4.2, 45.40873468),
            (2100.0, 4.3, 60.78123340),
            (2600.0, 5.1, 34.34999263),
            (2400.0, 4.7, 63.63924887),
            (2700.0, 5.0, 88.39043281),
        )
        for enthalpy, entropy, pressure_mpa in states:
            pressure_bar = pressure_mpa * 10
            with self.subTest(enthalpy=enthalpy, entropy=entropy):
                calculated_pressure = self.steam_table.p_hs(enthalpy, entropy)
                self.assertAlmostEqual(calculated_pressure, pressure_bar, delta=0.001)
                self.assertAlmostEqual(
                    self.steam_table.s_ph(calculated_pressure, enthalpy),
                    entropy,
                    delta=2e-4,
                )
                self.assertAlmostEqual(
                    self.steam_table.h_ps(calculated_pressure, entropy),
                    enthalpy,
                    delta=0.1,
                )
                self.assertAlmostEqual(
                    self.steam_table.t_hs(enthalpy, entropy),
                    self.steam_table.t_ph(calculated_pressure, enthalpy),
                    delta=0.001,
                )

    def test_region4_mixture_routes(self):
        pressure = 1.0
        quality = 0.5
        enthalpy = self.steam_table.h_px(pressure, quality)
        saturation_temperature = self.steam_table.tsat_p(pressure)

        self.assertAlmostEqual(enthalpy, 1546.193063, delta=1e-5)
        self.assertAlmostEqual(self.steam_table.t_ph(pressure, enthalpy), saturation_temperature, delta=1e-8)
        self.assertAlmostEqual(self.steam_table.x_ph(pressure, enthalpy), quality, delta=1e-10)

        entropy = self.steam_table.s_ph(pressure, enthalpy)
        self.assertTrue(math.isfinite(entropy))
        self.assertAlmostEqual(self.steam_table.h_ps(pressure, entropy), enthalpy, delta=1e-5)
        self.assertAlmostEqual(self.steam_table.t_ps(pressure, entropy), saturation_temperature, delta=1e-8)
        self.assertAlmostEqual(self.steam_table.x_ps(pressure, entropy), quality, delta=1e-10)

        mixture_volume = (1 - quality) * self.steam_table.vL_p(pressure) + quality * self.steam_table.vV_p(pressure)
        self.assertAlmostEqual(self.steam_table.v_ph(pressure, enthalpy), mixture_volume, delta=1e-10)
        self.assertAlmostEqual(self.steam_table.h_prho(pressure, 1 / mixture_volume), enthalpy, delta=1e-5)

    def test_single_phase_pt_routes_reject_saturation_and_invalid_states(self):
        saturation_temperature = self.steam_table.tsat_p(1.0)
        for method in ("h_pt", "u_pt", "s_pt", "Cp_pt", "Cv_pt", "w_pt"):
            with self.subTest(method=method, state="saturation"):
                self.assertTrue(math.isnan(getattr(self.steam_table, method)(1.0, saturation_temperature)))
            with self.subTest(method=method, state="invalid"):
                self.assertTrue(math.isnan(getattr(self.steam_table, method)(0.0, 100.0)))
