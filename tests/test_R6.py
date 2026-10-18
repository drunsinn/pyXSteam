# -*- coding: utf-8 -*-

import unittest
from pyXSteam import IAPWS_R6
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
