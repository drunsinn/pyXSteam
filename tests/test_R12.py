# -*- coding: utf-8 -*-

import unittest
import math
from unittest.mock import patch

from pyXSteam import IAPWS_R12, XSteam, UnitSystem
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

    def test_R12_correlation_length(self):
        """R12-08 Table 5 verification for the correlation length."""
        for temperature, density, reference in zip(R12_08.Table5_T, R12_08.Table5_rho, R12_08.Table5_xi):
            with self.subTest(temperature=temperature, density=density):
                self.assertAlmostEqual(
                    IAPWS_R12.R12_xi(density, temperature),
                    reference,
                    delta=5e-7,
                )

    def test_R12_critical_enhancement(self):
        """R12-08 Table 5 verification for both Y branches and full viscosity."""
        for temperature, density, expected_mu2, expected_viscosity in zip(
            R12_08.Table5_T,
            R12_08.Table5_rho,
            R12_08.Table5_my2dash,
            R12_08.Table5_my,
        ):
            with self.subTest(temperature=temperature, density=density):
                y = IAPWS_R12.R12_Y(density, temperature)
                mu2 = IAPWS_R12.R12_my_dash_2(density, temperature)
                viscosity = IAPWS_R12.my_rhoT(density, temperature, industrial_application=False)
                self.assertAlmostEqual(
                    y,
                    math.log(expected_mu2) / R12_08.Table3_x_u,
                    delta=1e-7,
                )
                self.assertAlmostEqual(mu2, expected_mu2, delta=1e-8)
                self.assertAlmostEqual(viscosity, expected_viscosity, delta=5e-13)

    def test_R12_critical_enhancement_temperature_gate(self):
        cases = (
            (645.91 - 1e-6, 0.0),
            (645.91, 0.0),
            (645.91 + 1e-6, 1.0),
            (650.77 - 1e-6, 1.0),
            (650.77, 0.0),
            (650.77 + 1e-6, 0.0),
        )
        for temperature, expected_calls in cases:
            with self.subTest(temperature=temperature):
                with patch("pyXSteam.IAPWS_R12.R12_Y", return_value=0.5) as y_function:
                    actual = IAPWS_R12.R12_my_dash_2(100.0, temperature)
                self.assertAlmostEqual(actual, math.exp(0.068 * 0.5) if expected_calls else 1.0)
                self.assertEqual(y_function.call_count, int(expected_calls))

    def test_R12_critical_enhancement_density_gate(self):
        cases = (
            (245.8 - 1e-6, 0.0),
            (245.8, 0.0),
            (245.8 + 1e-6, 1.0),
            (405.3 - 1e-6, 1.0),
            (405.3, 0.0),
            (405.3 + 1e-6, 0.0),
        )
        for density, expected_calls in cases:
            with self.subTest(density=density):
                with patch("pyXSteam.IAPWS_R12.R12_Y", return_value=0.5) as y_function:
                    actual = IAPWS_R12.R12_my_dash_2(density, 700.0)
                self.assertAlmostEqual(actual, math.exp(0.068 * 0.5) if expected_calls else 1.0)
                self.assertEqual(y_function.call_count, int(expected_calls))

    def test_R12_industrial_and_full_formulations(self):
        temperature = 647.35
        density = 322.0
        enhancement = IAPWS_R12.R12_my_dash_2(density, temperature)
        with patch(
            "pyXSteam.IAPWS_R12.R12_my_dash_2",
            wraps=IAPWS_R12.R12_my_dash_2,
        ) as enhancement_function:
            industrial_viscosity = IAPWS_R12.my_rhoT(density, temperature)
            self.assertEqual(enhancement_function.call_count, 0)
            full_viscosity = IAPWS_R12.my_rhoT(density, temperature, industrial_application=False)
            self.assertEqual(enhancement_function.call_count, 1)

        self.assertAlmostEqual(
            full_viscosity,
            industrial_viscosity * enhancement,
            delta=1e-18,
        )

    def test_R12_invalid_inputs(self):
        invalid_states = (
            (-1.0, 300.0),
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, -1.0),
            (float("nan"), 300.0),
            (1.0, float("nan")),
            (float("inf"), 300.0),
            (1.0, float("inf")),
        )
        for density, temperature in invalid_states:
            with self.subTest(density=density, temperature=temperature):
                for industrial in (True, False):
                    self.assertTrue(
                        math.isnan(
                            IAPWS_R12.my_rhoT(
                                density,
                                temperature,
                                industrial_application=industrial,
                            )
                        )
                    )
                self.assertTrue(math.isnan(IAPWS_R12.R12_xi(density, temperature)))
                self.assertTrue(math.isnan(IAPWS_R12.R12_Y(density, temperature)))
                self.assertTrue(math.isnan(IAPWS_R12.R12_my_dash_2(density, temperature)))

    def test_R12_zero_density_and_extrapolation(self):
        zero_density_viscosity = IAPWS_R12.my_rhoT(0.0, 300.0)
        self.assertTrue(math.isfinite(zero_density_viscosity))
        self.assertGreater(zero_density_viscosity, 0.0)
        self.assertEqual(IAPWS_R12.R12_my_dash_2(0.0, 647.35), 1.0)

        extrapolated = IAPWS_R12.my_rhoT(1.0, 1174.0)
        self.assertTrue(math.isfinite(extrapolated))
        self.assertGreater(extrapolated, 0.0)

    def test_R12_xsteam_wrapper_units(self):
        density_si = 998.0
        temperature_si = 298.15
        viscosity_si = IAPWS_R12.my_rhoT(density_si, temperature_si)

        mks = XSteam(UnitSystem.MKS)
        self.assertAlmostEqual(mks.R12_my_rhot(density_si, 25.0), viscosity_si, delta=1e-15)
        self.assertAlmostEqual(mks.my_pt(1.0, 100.0), 1.22704e-5, delta=1e-6)
        self.assertAlmostEqual(mks.my_ph(1.0, 100.0), 9.14003770302e-4, delta=1e-6)

        fls = XSteam(UnitSystem.FLS)
        density_fls = fls._unit_converter.fromSIunit_rho(density_si)
        temperature_fls = fls._unit_converter.fromSIunit_T(temperature_si)
        viscosity_fls = fls._unit_converter.fromSIunit_my(viscosity_si)
        self.assertAlmostEqual(fls.R12_my_rhot(density_fls, temperature_fls), viscosity_fls, delta=1e-10)
