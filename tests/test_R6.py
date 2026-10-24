# -*- coding: utf-8 -*-

import unittest
from pyXSteam import IAPWS_R6
from pyXSteam import XSteam
from pyXSteam.tables import R6_95
from . import helpers


class R6_FunctionTester(unittest.TestCase):

    def setUp(self):
        self.max_error = 1e-8
        self.max_matrix_error = self.max_error * 2

    def tearDown(self):
        pass

    def test_R6_internal(self):
        """check internal calculations with values from Table 6"""
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

        self.assertEqual(IAPWS_R6.eq_phi_o_deltatau(), 0.0)

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

    def test_residual_density_second_derivative(self):
        """Check the second density derivative against a finite difference near criticality."""
        temperature = 647.35
        tau = R6_95.CRITICAL_TEMPERATURE / temperature
        step = 1e-5

        for density in (222.0, 272.0, 372.0, 422.0):
            delta = density / R6_95.CRITICAL_DENSITY
            with self.subTest(density=density):
                numerical = (
                    IAPWS_R6.eq_phi_r_delta(tau, delta + step)
                    - IAPWS_R6.eq_phi_r_delta(tau, delta - step)
                ) / (2 * step)
                analytical = IAPWS_R6.eq_phi_r_deltadelta(tau, delta)
                self.assertAlmostEqual(analytical, numerical, delta=1e-8)

    def test_table7(self):
        """check functions against values of R6-95 Table 7"""

        helpers.array_2d_test(
            IAPWS_R6.R6_p_rhoT,
            (R6_95.Table7_rho, R6_95.Table7_T),
            R6_95.Table7_p,
            self.max_matrix_error,
        )

        helpers.array_2d_test(
            IAPWS_R6.R6_cv_rhoT,
            (R6_95.Table7_rho, R6_95.Table7_T),
            R6_95.Table7_cv,
            self.max_matrix_error,
        )

        helpers.array_2d_test(
            IAPWS_R6.R6_w_rhoT,
            (R6_95.Table7_rho, R6_95.Table7_T),
            R6_95.Table7_w,
            5.3e-3,
        )

        helpers.array_2d_test(
            IAPWS_R6.R6_s_rhoT,
            (R6_95.Table7_rho, R6_95.Table7_T),
            R6_95.Table7_s,
            self.max_matrix_error,
        )

    def test_table8(self):
        """Check the saturation properties in R6-95 Table 8."""
        pressure_tolerances = (3e-7, 4e-7, 3e-7)
        table8_max_matrix_error = 5e-8
        for phase, densities in (
            ("liquid", R6_95.Table_rho_dash),
            ("vapor", R6_95.Table_rho_double_dash),
        ):
            for i, (rho, temperature, pressure, tolerance) in enumerate(
                zip(
                    densities,
                    R6_95.Table8_T,
                    R6_95.Table8_p,
                    pressure_tolerances,
                )
            ):
                with self.subTest(phase=phase, temperature=temperature):
                    error = abs(IAPWS_R6.R6_p_rhoT(rho, temperature) - pressure)
                    self.assertLess(
                        error,
                        tolerance,
                        f"Error for {phase} pressure at Table 8 row {i}: {error}",
                    )

        helpers.array_2d_test(
            IAPWS_R6.R6_h_rhoT,
            (R6_95.Table_rho_dash, R6_95.Table8_T),
            R6_95.Table_h_dash,
            table8_max_matrix_error,
        )
        helpers.array_2d_test(
            IAPWS_R6.R6_h_rhoT,
            (R6_95.Table_rho_double_dash, R6_95.Table8_T),
            R6_95.Table_h_double_dash,
            table8_max_matrix_error,
        )
        helpers.array_2d_test(
            IAPWS_R6.R6_s_rhoT,
            (R6_95.Table_rho_dash, R6_95.Table8_T),
            R6_95.Table_s_dash,
            table8_max_matrix_error,
        )
        helpers.array_2d_test(
            IAPWS_R6.R6_s_rhoT,
            (R6_95.Table_rho_double_dash, R6_95.Table8_T),
            R6_95.Table_s_double_dash,
            table8_max_matrix_error,
        )

        for i, temperature in enumerate(R6_95.Table8_T):
            with self.subTest(phase_equilibrium_temperature=temperature):
                pressure_error, maxwell_error = IAPWS_R6.R6_phase_equilibrium_condition(
                    R6_95.Table_rho_dash[i],
                    R6_95.Table_rho_double_dash[i],
                    temperature,
                )
                self.assertAlmostEqual(pressure_error, 0.0, delta=1e-6)
                self.assertAlmostEqual(maxwell_error, 0.0, delta=2e-4)

    def test_second_virial_coefficient(self):
        """Check the published second virial coefficient at 600 K."""
        coefficient, derivative = IAPWS_R6.R6_sec_virial_coefficients(1.0, 600.0)
        self.assertAlmostEqual(coefficient, -0.555366808e-2, delta=5e-12)
        step = 0.01
        upper = IAPWS_R6.R6_sec_virial_coefficients(1.0, 600.0 + step)[0]
        lower = IAPWS_R6.R6_sec_virial_coefficients(1.0, 600.0 - step)[0]
        self.assertAlmostEqual(derivative, (upper - lower) / (2 * step), delta=1e-11)

    def test_third_virial_coefficient(self):
        """Check the published third virial coefficient at 600 K."""
        coefficient, derivative = IAPWS_R6.R6_third_virial_coefficients(1.0, 600.0)
        self.assertAlmostEqual(coefficient, -0.669015050e-5, delta=5e-15)
        step = 0.01
        upper = IAPWS_R6.R6_third_virial_coefficients(1.0, 600.0 + step)[0]
        lower = IAPWS_R6.R6_third_virial_coefficients(1.0, 600.0 - step)[0]
        self.assertAlmostEqual(derivative, (upper - lower) / (2 * step), delta=1e-14)

    def test_xsteam_r6_methods(self):
        """Check XSteam's MKS-unit wrappers against the R6 module."""
        steam = XSteam(XSteam.UNIT_SYSTEM_MKS)
        rho = R6_95.Table7_rho[5]
        temperature = steam._unit_converter.fromSIunit_T(R6_95.Table7_T[5])
        temperature_si = steam._unit_converter.toSIunit_T(temperature)
        reference_properties = {
            "p": steam._unit_converter.fromSIunit_p(IAPWS_R6.R6_p_rhoT(rho, temperature_si)),
            "u": steam._unit_converter.fromSIunit_u(IAPWS_R6.R6_u_rhoT(rho, temperature_si)),
            "s": steam._unit_converter.fromSIunit_s(IAPWS_R6.R6_s_rhoT(rho, temperature_si)),
            "h": steam._unit_converter.fromSIunit_h(IAPWS_R6.R6_h_rhoT(rho, temperature_si)),
            "cv": steam._unit_converter.fromSIunit_Cv(IAPWS_R6.R6_cv_rhoT(rho, temperature_si)),
            "cp": steam._unit_converter.fromSIunit_Cp(IAPWS_R6.R6_cp_rhoT(rho, temperature_si)),
            "w": steam._unit_converter.fromSIunit_w(IAPWS_R6.R6_w_rhoT(rho, temperature_si)),
            "joulethomson": steam._unit_converter.fromSIunit_pressure_derivative(IAPWS_R6.R6_joulethomson_rhoT(rho, temperature_si)),
            "delta_T": steam._unit_converter.fromSIunit_v(IAPWS_R6.R6_delta_T_rhoT(rho, temperature_si)),
            "isen_t_p": steam._unit_converter.fromSIunit_pressure_derivative(IAPWS_R6.R6_isen_t_p_rhoT(rho, temperature_si)),
        }
        for method_name, reference in reference_properties.items():
            with self.subTest(property=method_name):
                method = getattr(steam, f"{method_name}_rhoT")
                self.assertAlmostEqual(method(rho, temperature), reference)

        volume_scale = steam._unit_converter.fromSIunit_v(1)
        temperature_scale = steam._unit_converter.toSIunit_T(1) - steam._unit_converter.toSIunit_T(0)
        for method_name, module_method in (
            ("sec", IAPWS_R6.R6_sec_virial_coefficients),
            ("third", IAPWS_R6.R6_third_virial_coefficients),
        ):
            with self.subTest(virial=method_name):
                coefficient, derivative = module_method(rho, temperature_si)
                if method_name == "sec":
                    steam_result = steam.sec_virial_coefficients(rho, temperature)
                    coefficient_scale = volume_scale
                else:
                    steam_result = steam.third_virial_coefficients(rho, temperature)
                    coefficient_scale = volume_scale**2
                self.assertAlmostEqual(steam_result[0], coefficient * coefficient_scale)
                self.assertAlmostEqual(
                    steam_result[1],
                    derivative * coefficient_scale * temperature_scale,
                )

        self.assertEqual(
            steam.phase_equilibrium_condition(
                R6_95.Table_rho_dash[1],
                R6_95.Table_rho_double_dash[1],
                450 - 273.15,
            ),
            (
                steam._unit_converter.fromSIunit_p(
                    IAPWS_R6.R6_phase_equilibrium_condition(
                        R6_95.Table_rho_dash[1],
                        R6_95.Table_rho_double_dash[1],
                        450,
                    )[0]
                ),
                IAPWS_R6.R6_phase_equilibrium_condition(
                    R6_95.Table_rho_dash[1],
                    R6_95.Table_rho_double_dash[1],
                    450,
                )[1],
            ),
        )

    def test_xsteam_r6_fls_units(self):
        """Check R6 wrappers convert density, temperature, and pressure in FLS."""
        steam = XSteam(XSteam.UNIT_SYSTEM_FLS)
        rho_si = R6_95.Table7_rho[5]
        temperature_si = R6_95.Table7_T[5]
        rho_fls = steam._unit_converter.fromSIunit_rho(rho_si)
        temperature_fls = steam._unit_converter.fromSIunit_T(temperature_si)

        pressure_si = IAPWS_R6.R6_p_rhoT(rho_si, temperature_si)
        enthalpy_si = IAPWS_R6.R6_h_rhoT(rho_si, temperature_si)
        self.assertAlmostEqual(
            steam.p_rhoT(rho_fls, temperature_fls),
            steam._unit_converter.fromSIunit_p(pressure_si),
        )
        self.assertAlmostEqual(
            steam.h_rhoT(rho_fls, temperature_fls),
            steam._unit_converter.fromSIunit_h(enthalpy_si),
        )
        for name, module_method in (
            ("joulethomson_rhoT", IAPWS_R6.R6_joulethomson_rhoT),
            ("isen_t_p_rhoT", IAPWS_R6.R6_isen_t_p_rhoT),
        ):
            with self.subTest(derivative=name):
                expected = steam._unit_converter.fromSIunit_pressure_derivative(module_method(rho_si, temperature_si))
                self.assertAlmostEqual(getattr(steam, name)(rho_fls, temperature_fls), expected)

        coefficient_si = IAPWS_R6.R6_third_virial_coefficients(rho_si, temperature_si)
        coefficient_fls = steam.third_virial_coefficients(rho_fls, temperature_fls)
        volume_scale = steam._unit_converter.fromSIunit_v(1)
        temperature_scale = steam._unit_converter.toSIunit_T_interval(1)
        self.assertAlmostEqual(coefficient_fls[0], coefficient_si[0] * volume_scale**2)
        self.assertAlmostEqual(
            coefficient_fls[1],
            coefficient_si[1] * volume_scale**2 * temperature_scale,
        )
