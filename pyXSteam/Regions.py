#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Section 2: IAPWS IF 97 Calling functions
"""
import math
import logging
from .RegionBorders import TB23_p
from .Constants import (
    SPECIFIC_GAS_CONSTANT,
    CRITICAL_TEMPERATURE,
    CRITICAL_PRESSURE,
    CRITICAL_DENSITY,
    TRIPLE_POINT_PRESSURE,
    FREEZING_TEMPERATURE_H2O,
)
from .IF97_Tables import R1, R2, R3, R4, R5


class Region1:
    """
    Section 2.1: IAPWS IF 97 Calling functions to calculate the properties
    of water in Region 1

    Release on the IAPWS Industrial formulation 1997 for the Thermodynamic
    Properties of Water and Steam, September 1997
    """

    _logger = logging.getLogger(__name__)

    @classmethod
    def v1_pT(cls, p: float, T: float) -> float:
        """function v1_pT = v1_pT(p, T)

        5 Equations for Region 1, Section. 5.1 Basic Equation

        Equation 7, Table 3, Page 6

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific volume in [m³ / kg]
        """
        # I1 = R1.Table2_I
        # J1 = R1.Table2_J
        # n1 = R1.Table2_n
        R = SPECIFIC_GAS_CONSTANT
        Pi = p / 16.53  # 16.53 MPa
        tau = 1386 / T  # 1386 K
        gamma_der_pi = 0
        for I, J, n in zip(R1.Table2_I, R1.Table2_J, R1.Table2_n):
            # for i in range(0, 34):
            gamma_der_pi = gamma_der_pi - n * I * (7.1 - Pi) ** (I - 1) * (tau - 1.222) ** J
        return R * T / p * Pi * gamma_der_pi / 1000

    @classmethod
    def h1_pT(cls, p: float, T: float) -> float:
        """
        calculate enthalpy from preasure and temperature in region 1

        5 Equations for Region 1, Section. 5.1 Basic Equation

        Equation 7, Table 3, Page 6

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: enthalpy in [kJ / kg]
        """
        # I1 = R1.Table2_I
        # J1 = R1.Table2_J
        # n1 = R1.Table2_n
        R = SPECIFIC_GAS_CONSTANT
        Pi = p / 16.53
        tau = 1386 / T
        gamma_der_tau = 0
        for I, J, n in zip(R1.Table2_I, R1.Table2_J, R1.Table2_n):
            # for i in range(0, 34):
            gamma_der_tau = gamma_der_tau + (n * (7.1 - Pi) ** I * J * (tau - 1.222) ** (J - 1))
        return R * T * tau * gamma_der_tau

    @classmethod
    def u1_pT(cls, p: float, T: float) -> float:
        """function u1_pT = u1_pT(p, T)

        5 Equations for Region 1, Section. 5.1 Basic Equation

        Equation 7, Table 3, Page 6

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific internal energy in [kJ / kg]
        """
        # I1 = R1.Table2_I
        # J1 = R1.Table2_J
        # n1 = R1.Table2_n
        R = SPECIFIC_GAS_CONSTANT
        Pi = p / 16.53
        tau = 1386 / T
        gamma_der_tau = 0
        gamma_der_pi = 0
        for I, J, n in zip(R1.Table2_I, R1.Table2_J, R1.Table2_n):
            # for i in range(0, 34):
            gamma_der_pi = gamma_der_pi - n * I * (7.1 - Pi) ** (I - 1) * (tau - 1.222) ** J
            gamma_der_tau = gamma_der_tau + (n * (7.1 - Pi) ** I * J * (tau - 1.222) ** (J - 1))
        return R * T * (tau * gamma_der_tau - Pi * gamma_der_pi)

    @classmethod
    def s1_pT(cls, p: float, T: float) -> float:
        """function s1_pT = s1_pT(p, T)

        5 Equations for Region 1, Section. 5.1 Basic Equation

        Equation 7, Table 3, Page 6

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific entropy in [kJ / (kg K)]
        """
        # I1 = R1.Table2_I
        # J1 = R1.Table2_J
        # n1 = R1.Table2_n
        R = SPECIFIC_GAS_CONSTANT
        Pi = p / 16.53
        tau = 1386 / T
        gamma = 0
        gamma_der_tau = 0
        for I, J, n in zip(R1.Table2_I, R1.Table2_J, R1.Table2_n):
            # for i in range(0, 34):
            gamma_der_tau = gamma_der_tau + (n * (7.1 - Pi) ** I * J * (tau - 1.222) ** (J - 1))
            gamma = gamma + n * (7.1 - Pi) ** I * (tau - 1.222) ** J
        return R * tau * gamma_der_tau - R * gamma

    @classmethod
    def Cp1_pT(cls, p: float, T: float) -> float:
        """function Cp1_pT = Cp1_pT(p, T)

        5 Equations for Region 1, Section. 5.1 Basic Equation

        Equation 7, Table 3, Page 6

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific isobaric heat capacity in [kJ / (kg K)]
        """
        # I1 = R1.Table2_I
        # J1 = R1.Table2_J
        # n1 = R1.Table2_n
        R = SPECIFIC_GAS_CONSTANT
        Pi = p / 16.53
        tau = 1386 / T
        gamma_der_tautau = 0
        for I, J, n in zip(R1.Table2_I, R1.Table2_J, R1.Table2_n):
            # for i in range(0, 34):
            gamma_der_tautau = gamma_der_tautau + (n * (7.1 - Pi) ** I * J * (J - 1) * (tau - 1.222) ** (J - 2))
        return -R * tau**2 * gamma_der_tautau

    @classmethod
    def Cv1_pT(cls, p: float, T: float) -> float:
        """function Cv1_pT = Cv1_pT(p, T)

        5 Equations for Region 1, Section. 5.1 Basic Equation

        Equation 7, Table 3, Page 6

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific isochoric heat capacity in [kJ / (kg K)]
        """
        R = SPECIFIC_GAS_CONSTANT
        Pi = p / 16.53
        tau = 1386 / T
        gamma_der_pi = 0
        gamma_der_pipi = 0
        gamma_der_pitau = 0
        gamma_der_tautau = 0
        for I, J, n in zip(R1.Table2_I, R1.Table2_J, R1.Table2_n):
            # for i in range(0, 34):
            gamma_der_pi = gamma_der_pi - n * I * (7.1 - Pi) ** (I - 1) * (tau - 1.222) ** J
            gamma_der_pipi = gamma_der_pipi + n * I * (I - 1) * (7.1 - Pi) ** (I - 2) * (tau - 1.222) ** J
            gamma_der_pitau = gamma_der_pitau - n * I * (7.1 - Pi) ** (I - 1) * J * (tau - 1.222) ** (J - 1)
            gamma_der_tautau = gamma_der_tautau + n * (7.1 - Pi) ** I * J * (J - 1) * (tau - 1.222) ** (J - 2)
        return R * (-(tau**2) * gamma_der_tautau + (gamma_der_pi - tau * gamma_der_pitau) ** 2 / gamma_der_pipi)

    @classmethod
    def w1_pT(cls, p: float, T: float) -> float:
        """function w1_pT = w1_pT(p, T)

        5 Equations for Region 1, Section. 5.1 Basic Equation

        Equation 7, Table 3, Page 6

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: speed of sound in [m / s]
        """
        R = SPECIFIC_GAS_CONSTANT
        Pi = p / 16.53
        tau = 1386 / T
        gamma_der_pi = 0
        gamma_der_pipi = 0
        gamma_der_pitau = 0
        gamma_der_tautau = 0
        for I, J, n in zip(R1.Table2_I, R1.Table2_J, R1.Table2_n):
            # for i in range(0, 34):
            gamma_der_pi = gamma_der_pi - n * I * (7.1 - Pi) ** (I - 1) * (tau - 1.222) ** J
            gamma_der_pipi = gamma_der_pipi + n * I * (I - 1) * (7.1 - Pi) ** (I - 2) * (tau - 1.222) ** J
            gamma_der_pitau = gamma_der_pitau - n * I * (7.1 - Pi) ** (I - 1) * J * (tau - 1.222) ** (J - 1)
            gamma_der_tautau = gamma_der_tautau + n * (7.1 - Pi) ** I * J * (J - 1) * (tau - 1.222) ** (J - 2)
        return (
            1000 * R * T * gamma_der_pi**2 / ((gamma_der_pi - tau * gamma_der_pitau) ** 2 / (tau**2 * gamma_der_tautau) - gamma_der_pipi)
        ) ** 0.5

    @classmethod
    def T1_ph(cls, p: float, h: float) -> float:
        """function T1_ph = T1_ph(p, h)

        5 Equations for Region 1, Section. 5.1 Basic Equation, 5.2.1 The Backward
        Equation T (p, h)

        Equation 11, Table 6, Page 10

        :param p: preasure in [MPa]
        :param h: enthalpy in [kJ / kg]

        :return: temperature in [K]
        """
        Pi = p / 1
        eta = h / 2500
        T = 0
        for I, J, n in zip(R1.Table6_I, R1.Table6_J, R1.Table6_n):
            # for i in range(0, 20):
            T = T + n * Pi**I * (eta + 1) ** J
        return T

    @classmethod
    def T1_ps(cls, p: float, s: float) -> float:
        """function T1_ps = T1_ps(p, s)

        5 Equations for Region 1, Section. 5.1 Basic Equation, 5.2.2 The Backward
        Equation T (p, s)

        Equation 13, Table 8, Page 11

        :param p: preasure in [MPa]
        :param s: specific entropy in [kJ / (kg K)]

        :return: temperature in [K]
        """
        Pi = p / 1
        Sigma = s / 1
        T = 0
        for I, J, n in zip(R1.Table8_I, R1.Table8_J, R1.Table8_n):
            # for i in range(0, 20):
            T = T + n * Pi**I * (Sigma + 2) ** J
        return T

    @classmethod
    def p1_hs(cls, h: float, s: float) -> float:
        """function p1_hs = p1_hs(h, s)

        Supplementary Release on Backward Equations for Pressure as a Function of
        Enthalpy and Entropy p(h, s) to the IAPWS Industrial formulation 1997 for
        the Thermodynamic Properties of Water and Steam

        5 Backward Equation p(h, s) for Region 1

        Equation 1, Table 2, Page 5

        :param h: enthalpy in [kJ / kg]
        :param s: specific entropy in [kJ / (kg K)]

        :return: preasure in [MPa]
        """
        eta = h / 3400
        Sigma = s / 7.6
        p = 0
        for I, J, n in zip(R1.Sub_psh12_Table2_I, R1.Sub_psh12_Table2_J, R1.Sub_psh12_Table2_n):
            # for i in range(0, 19):
            p = p + n * (eta + 0.05) ** I * (Sigma + 0.05) ** J
        return p * 100

    @classmethod
    def T1_prho(cls, p: float, rho: float) -> float:
        """function T1_prho = T1_prho(p , rho)

        Solve by iteration. Observe that for low temperatures this equation
        has 2 solutions. Solve with half interval method

        :param p: pressure in [MPa]
        :param rho: density in [kg / m³]

        :return: temperature in [K]
        """
        Ts = float("NaN")
        Low_Bound = FREEZING_TEMPERATURE_H2O
        High_Bound = Region4.T4_p(p)
        rhos = -1000
        step_counter = 0
        while math.fabs(rho - rhos) > 0.00001:
            step_counter += 1
            last_rhos = rhos

            Ts = (Low_Bound + High_Bound) / 2
            rhos = 1 / Region1.v1_pT(p, Ts)

            if last_rhos == rhos:
                cls._logger.warning(
                    "T1_prho stopped iterating after %d steps because values did not converge",
                    step_counter,
                )
                break

            if rhos < rho:
                High_Bound = Ts
            else:
                Low_Bound = Ts
        return Ts


class Region2:
    """
    Section 2.2: IAPWS IF 97 Calling functions to calculate the properties
    of water in Region 2
    """

    _logger = logging.getLogger(__name__)

    @classmethod
    def v2_pT(cls, p: float, T: float) -> float:
        """function v2_pT = v2_pT(p, T)

        6 Equations for Region 2, Section. 6.1 Basic Equation

        Table 11 and 12, Page 14 and 15

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific volume in [m³ / kg]
        """
        R = SPECIFIC_GAS_CONSTANT  # Eq 1
        Pi = p  # Eq 1
        tau = 540 / T  # Eq 1
        g0_pi = 1 / Pi  # see table 13
        gr_pi = 0
        for I, J, n in zip(R2.Table11_I, R2.Table11_J, R2.Table11_n):
            # for i in range(0, 43):
            gr_pi = gr_pi + n * I * Pi ** (I - 1) * (tau - 0.5) ** J  # see table 14
        return R * T / p * Pi * (g0_pi + gr_pi) / 1000  # see table 12

    @staticmethod
    def v2_pT_meta(p, T):
        """
        6 Equations for Region 2, Section. 6.2 Supplementary Equation for the Metastable-Vapor Region

        Table 16, Page 18

        specific volume
        """
        R = SPECIFIC_GAS_CONSTANT  # Eq 1
        Pi = p  # Eq 1
        tau = 540 / T  # Eq 1

        # table 13 - dimensionless gibbs free energy - gamma 0 pi
        g0_pi = 1 / Pi

        # table 14 - residual dimensionless gibbs free energy - part r pi
        gr_pi = 0
        for I, J, n in zip(R2.Table16_I, R2.Table16_J, R2.Table16_n):
            # for i, n in enumerate(nr):
            gr_pi = gr_pi + n * I * Pi ** (I - 1) * (tau - 0.5) ** J

        return R * T / p * Pi * (g0_pi + gr_pi) / 1000  # see table 12

    @staticmethod
    def h2_pT(p, T):
        """function h2_pT = h2_pT(p, T)

        6 Equations for Region 2, Section. 6.1 Basic Equation

        Table 11 and 12, Page 14 and 15

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: enthalpy in [kJ / kg]
        """
        R = SPECIFIC_GAS_CONSTANT
        Pi = p
        tau = 540 / T
        g0_tau = 0
        for J, n in zip(R2.Table10_J0, R2.Table10_n0):
            # for i in range(0, 9):
            g0_tau = g0_tau + n * J * tau ** (J - 1)
        gr_tau = 0
        for I, J, n in zip(R2.Table11_I, R2.Table11_J, R2.Table11_n):
            # for i in range(0, 43):
            gr_tau = gr_tau + n * Pi**I * J * (tau - 0.5) ** (J - 1)
        return R * T * tau * (g0_tau + gr_tau)

    @staticmethod
    def h2_pT_meta(p, T):
        """
        6 Equations for Region 2, Section. 6.2 Supplementary Equation for the Metastable-Vapor Region

        Table 16, Page 18

        specific volume
        """
        R = SPECIFIC_GAS_CONSTANT  # Eq 1
        Pi = p  # Eq 1
        tau = 540 / T  # Eq 1

        # table 13 - dimensionless gibbs free energy - gamma 0 tau
        g0_tau = 0
        for J, n in zip(R2.Table10_J0, R2.Table10_n0):
            g0_tau += n * J * tau ** (J - 1)

        # table 14 - residual dimensionless gibbs free energy - part r tau
        gr_tau = 0
        for I, J, n in zip(R2.Table16_I, R2.Table16_J, R2.Table16_n):
            gr_tau += n * Pi**I * J * (tau - 0.5) ** (J - 1)

        return R * T * tau * (g0_tau + gr_tau)  # h2_pT

    @staticmethod
    def u2_pT(p, T):
        """function u2_pT = u2_pT(p, T)

        6 Equations for Region 2, Section. 6.1 Basic Equation

        Table 11 and 12, Page 14 and 15

        :param p: preasure in [MPa]
        :param T: temperature in [K]
        :return: specific internal energy in [kJ / kg]
        """
        R = SPECIFIC_GAS_CONSTANT
        Pi = p
        tau = 540 / T
        g0_pi = 1 / Pi
        g0_tau = 0

        for J, n in zip(R2.Table10_J0, R2.Table10_n0):
            g0_tau = g0_tau + n * J * tau ** (J - 1)

        gr_pi = 0
        gr_tau = 0

        for I, J, n in zip(R2.Table11_I, R2.Table11_J, R2.Table11_n):
            gr_pi = gr_pi + n * I * Pi ** (I - 1) * (tau - 0.5) ** J
            gr_tau = gr_tau + n * Pi**I * J * (tau - 0.5) ** (J - 1)

        return R * T * (tau * (g0_tau + gr_tau) - Pi * (g0_pi + gr_pi))

    @staticmethod
    def u2_pT_meta(p, T):
        """
        6 Equations for Region 2, Section. 6.2 Supplementary Equation for the Metastable-Vapor Region

        Table 16, Page 18

        specific volume
        """
        R = SPECIFIC_GAS_CONSTANT  # Eq 1
        Pi = p  # Eq 1
        tau = 540 / T  # Eq 1

        # table 13 - dimensionless gibbs free energy - gamma 0 pi
        g0_pi = 1 / Pi

        # table 13 - dimensionless gibbs free energy - gamma 0 tau
        g0_tau = 0
        for J, n in zip(R2.Table10_J0, R2.Table10_n0):
            g0_tau += n * J * tau ** (J - 1)

        # table 14 - residual dimensionless gibbs free energy - part r pi
        gr_pi = 0
        for I, J, n in zip(R2.Table16_I, R2.Table16_J, R2.Table16_n):
            gr_pi = gr_pi + n * I * Pi ** (I - 1) * (tau - 0.5) ** J

        # table 14 - residual dimensionless gibbs free energy - part r tau
        gr_tau = 0
        for I, J, n in zip(R2.Table16_I, R2.Table16_J, R2.Table16_n):
            gr_tau += n * Pi**I * J * (tau - 0.5) ** (J - 1)

        return R * T * (tau * (g0_tau + gr_tau) - Pi * (g0_pi + gr_pi))  # u2_pT

    @staticmethod
    def s2_pT(p, T):
        """function s2_pT = s2_pT(p, T)

        6 Equations for Region 2, Section. 6.1 Basic Equation

        Table 11 and 12, Page 14 and 15

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific entropy in [kJ / (kg K)]
        """
        R = SPECIFIC_GAS_CONSTANT
        Pi = p
        tau = 540 / T
        g0 = math.log(Pi)
        g0_tau = 0
        for J, n in zip(R2.Table10_J0, R2.Table10_n0):
            g0 = g0 + n * tau**J
            g0_tau = g0_tau + n * J * tau ** (J - 1)

        gr = 0
        gr_tau = 0
        for I, J, n in zip(R2.Table11_I, R2.Table11_J, R2.Table11_n):
            gr = gr + n * Pi**I * (tau - 0.5) ** J
            gr_tau = gr_tau + n * Pi**I * J * (tau - 0.5) ** (J - 1)
        return R * (tau * (g0_tau + gr_tau) - (g0 + gr))

    @staticmethod
    def s2_pT_meta(p, T):
        """
        6 Equations for Region 2, Section. 6.2 Supplementary Equation for the Metastable-Vapor Region

        Table 16, Page 18

        specific volume
        """
        R = SPECIFIC_GAS_CONSTANT  # Eq 1
        Pi = p  # Eq 1
        tau = 540 / T  # Eq 1

        # table 13 - dimensionless gibbs free energy - gamma 0
        g0 = math.log(Pi)
        for J, n in zip(R2.Table10_J0, R2.Table10_n0):
            g0 += n * tau**J

        # table 13 - dimensionless gibbs free energy - gamma 0 tau
        g0_tau = 0
        for J, n in zip(R2.Table10_J0, R2.Table10_n0):
            g0_tau += n * J * tau ** (J - 1)

        # table 14 - residual dimensionless gibbs free energy - part r
        gr = 0
        for I, J, n in zip(R2.Table16_I, R2.Table16_J, R2.Table16_n):
            gr += n * Pi**I * (tau - 0.5) ** J

        # table 14 - residual dimensionless gibbs free energy - part r tau
        gr_tau = 0
        for I, J, n in zip(R2.Table16_I, R2.Table16_J, R2.Table16_n):
            gr_tau += n * Pi**I * J * (tau - 0.5) ** (J - 1)

        return R * (tau * (g0_tau + gr_tau) - (g0 + gr))  # s2_pT

    @staticmethod
    def Cp2_pT(p, T):
        """function Cp2_pT = Cp2_pT(p, T)

        6 Equations for Region 2, Section. 6.1 Basic Equation

        Table 11 and 12, Page 14 and 15

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific isobaric heat capacity in [kJ / (kg K)]
        """
        R = SPECIFIC_GAS_CONSTANT
        Pi = p
        tau = 540 / T

        g0_tautau = 0
        for J, n in zip(R2.Table10_J0, R2.Table10_n0):
            g0_tautau = g0_tautau + n * J * (J - 1) * tau ** (J - 2)

        gr_tautau = 0
        for I, J, n in zip(R2.Table11_I, R2.Table11_J, R2.Table11_n):
            gr_tautau = gr_tautau + n * Pi**I * J * (J - 1) * (tau - 0.5) ** (J - 2)
        return -R * tau**2 * (g0_tautau + gr_tautau)

    @staticmethod
    def Cp2_pT_meta(p, T):
        """
        6 Equations for Region 2, Section. 6.2 Supplementary Equation for the Metastable-Vapor Region

        Table 16, Page 18

        specific volume
        """
        R = SPECIFIC_GAS_CONSTANT  # Eq 1
        Pi = p  # Eq 1
        tau = 540 / T  # Eq 1

        # table 13 - dimensionless gibbs free energy - gamma 0 tau tau
        g0_tautau = 0
        for J, n in zip(R2.Table10_J0, R2.Table10_n0):
            g0_tautau += n * J * (J - 1) * tau ** (J - 2)

        # table 14 - residual dimensionless gibbs free energy - part r tau tau
        gr_tautau = 0
        for I, J, n in zip(R2.Table16_I, R2.Table16_J, R2.Table16_n):
            gr_tautau += n * Pi**I * J * (J - 1) * (tau - 0.5) ** (J - 2)

        return -R * tau**2 * (g0_tautau + gr_tautau)  # Cp2_pT

    @staticmethod
    def Cv2_pT(p, T):
        """function Cv2_pT = Cv2_pT(p, T)

        6 Equations for Region 2, Section. 6.1 Basic Equation

        Table 11 and 12, Page 14 and 15

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific isochoric heat capacity in [kJ / (kg K)]
        """
        R = SPECIFIC_GAS_CONSTANT
        Pi = p
        tau = 540 / T

        g0_tautau = 0
        for J, n in zip(R2.Table10_J0, R2.Table10_n0):
            g0_tautau = g0_tautau + n * J * (J - 1) * tau ** (J - 2)

        gr_pi = 0
        gr_pitau = 0
        gr_pipi = 0
        gr_tautau = 0
        for I, J, n in zip(R2.Table11_I, R2.Table11_J, R2.Table11_n):
            gr_pi = gr_pi + n * I * Pi ** (I - 1) * (tau - 0.5) ** J
            gr_pipi = gr_pipi + n * I * (I - 1) * Pi ** (I - 2) * (tau - 0.5) ** J
            gr_pitau = gr_pitau + n * I * Pi ** (I - 1) * J * (tau - 0.5) ** (J - 1)
            gr_tautau = gr_tautau + n * Pi**I * J * (J - 1) * (tau - 0.5) ** (J - 2)
        return R * (-(tau**2) * (g0_tautau + gr_tautau) - (1 + Pi * gr_pi - tau * Pi * gr_pitau) ** 2 / (1 - Pi**2 * gr_pipi))

    @classmethod
    def w2_pT(cls, p: float, T: float) -> float:
        """function w2_pT = w2_pT(p, T)

        6 Equations for Region 2, Section. 6.1 Basic Equation

        Table 11 and 12, Page 14 and 15

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: speed of sound in [m / s]
        """
        R = SPECIFIC_GAS_CONSTANT
        Pi = p
        tau = 540 / T

        g0_tautau = 0
        for J, n in zip(R2.Table10_J0, R2.Table10_n0):
            g0_tautau = g0_tautau + n * J * (J - 1) * tau ** (J - 2)

        gr_pi = 0
        gr_pitau = 0
        gr_pipi = 0
        gr_tautau = 0
        for I, J, n in zip(R2.Table11_I, R2.Table11_J, R2.Table11_n):
            gr_pi = gr_pi + n * I * Pi ** (I - 1) * (tau - 0.5) ** J
            gr_pipi = gr_pipi + n * I * (I - 1) * Pi ** (I - 2) * (tau - 0.5) ** J
            gr_pitau = gr_pitau + n * I * Pi ** (I - 1) * J * (tau - 0.5) ** (J - 1)
            gr_tautau = gr_tautau + n * Pi**I * J * (J - 1) * (tau - 0.5) ** (J - 2)
        return (
            1000
            * R
            * T
            * (1 + 2 * Pi * gr_pi + Pi**2 * gr_pi**2)
            / ((1 - Pi**2 * gr_pipi) + (1 + Pi * gr_pi - tau * Pi * gr_pitau) ** 2 / (tau**2 * (g0_tautau + gr_tautau)))
        ) ** 0.5

    @staticmethod
    def w2_pT_meta(p, T):
        """
        6 Equations for Region 2, Section. 6.2 Supplementary Equation for the Metastable-Vapor Region

        Table 16, Page 18

        specific volume
        """
        R = SPECIFIC_GAS_CONSTANT  # Eq 1
        Pi = p  # Eq 1
        tau = 540 / T  # Eq 1

        # table 13 - dimensionless gibbs free energy - gamma 0 tau tau
        g0_tautau = 0
        for J, n in zip(R2.Table10_J0, R2.Table10_n0):
            g0_tautau += n * J * (J - 1) * tau ** (J - 2)

        # table 14 - residual dimensionless gibbs free energy - part r pi
        gr_pi = 0
        for I, J, n in zip(R2.Table16_I, R2.Table16_J, R2.Table16_n):
            gr_pi = gr_pi + n * I * Pi ** (I - 1) * (tau - 0.5) ** J

        # table 14 - residual dimensionless gibbs free energy - part r pi pi
        gr_pipi = 0
        for I, J, n in zip(R2.Table16_I, R2.Table16_J, R2.Table16_n):
            gr_pipi += n * I * (I - 1) * Pi ** (I - 2) * (tau - 0.5) ** J
        if isinstance(gr_pipi, complex):
            if gr_pipi.imag != 0:
                raise Exception()
            gr_pipi = gr_pipi.real

        # table 14 - residual dimensionless gibbs free energy - part r tau tau
        gr_tautau = 0
        for I, J, n in zip(R2.Table16_I, R2.Table16_J, R2.Table16_n):
            gr_tautau += n * Pi**I * J * (J - 1) * (tau - 0.5) ** (J - 2)

        # table 14 - residual dimensionless gibbs free energy - part r pi tau
        gr_pitau = 0
        for I, J, n in zip(R2.Table16_I, R2.Table16_J, R2.Table16_n):
            gr_pitau += n * I * Pi ** (I - 1) * J * (tau - 0.5) ** (J - 1)

        part_1 = 1 + 2 * Pi * gr_pi + Pi**2 * gr_pi**2
        part_2_a = 1 - Pi**2 * gr_pipi
        part_2_b = (1 + Pi * gr_pi - tau * Pi * gr_pitau) ** 2  # x
        part_2_c = tau**2 * (g0_tautau + gr_tautau)
        return math.sqrt(1000 * R * T * part_1 / (part_2_a + part_2_b / part_2_c))

    @staticmethod
    def T2_ph(p, h):
        """function T2_ph = T2_ph(p, h)

        6 Equations for Region 2, 6.3.1 The Backward Equations T(p, h) for
        Subregions 2a, 2b, and 2c

        :param p: preasure in [MPa]
        :param h: enthalpy in [kJ / kg]

        :return: temperature in [K]
        """
        if p < 4:
            sub_reg = 1
        else:
            if p < (905.84278514723 - 0.67955786399241 * h + 1.2809002730136e-04 * h**2):
                sub_reg = 2
            else:
                sub_reg = 3
        if sub_reg == 1:
            Ts = 0
            hs = h / 2000
            for I, J, n in zip(R2.Table20_I, R2.Table20_J, R2.Table20_n):
                # for i in range(0, 34):
                Ts = Ts + n * p**I * (hs - 2.1) ** J
        elif sub_reg == 2:
            # Subregion B
            # Table 21, Eq 23, page 23
            Ts = 0
            hs = h / 2000
            # TODO check list range!
            # for i = 1 : 38
            for I, J, n in zip(R2.Table21_I, R2.Table21_J, R2.Table21_n):
                # for i in range(0, 38):
                Ts = Ts + n * (p - 2) ** I * (hs - 2.6) ** J
        else:
            # Subregion C
            # Table 22, Eq 24, page 24
            Ts = 0
            hs = h / 2000
            for I, J, n in zip(R2.Table22_I, R2.Table22_J, R2.Table22_n):
                # for i in range(0, 23):
                Ts = Ts + n * (p + 25) ** I * (hs - 1.8) ** J
        return Ts

    @classmethod
    def T2_ps(cls, p: float, s: float) -> float:
        """function T2_ps = T2_ps(p, s)

        6 Equations for Region 2,6.3.2 The Backward Equations T( p, s ) for
        Subregions 2a, 2b, and 2c

        Page 26

        :param p: preasure in [MPa]
        :param s: specific entropy in [kJ / (kg K)]

        :return: temperature in [K]
        """
        if p < 4:
            sub_reg = 1
        else:
            if s < 5.85:
                sub_reg = 3
            else:
                sub_reg = 2
        if sub_reg == 1:
            # Subregion A
            # Table 25, Eq 25, page 26
            Pi = p
            Sigma = s / 2
            teta = 0
            # TODO check list range!
            # for i = 1 : 46
            for I, J, n in zip(R2.Table25_I, R2.Table25_J, R2.Table25_n):
                # for i in range(0, 46):
                teta = teta + n * Pi**I * (Sigma - 2) ** J
        elif sub_reg == 2:
            # Subregion B
            # Table 26, Eq 26, page 27
            Pi = p
            Sigma = s / 0.7853
            teta = 0
            for I, J, n in zip(R2.Table26_I, R2.Table26_J, R2.Table26_n):
                # for i in range(0, 44):
                teta = teta + n * Pi**I * (10 - Sigma) ** J
        else:
            # Subregion C
            # Table 27, Eq 27, page 28
            Pi = p
            Sigma = s / 2.9251
            teta = 0
            for I, J, n in zip(R2.Table27_I, R2.Table27_J, R2.Table27_n):
                # for i in range(0, 30):
                teta = teta + n * Pi**I * (2 - Sigma) ** J
        return teta

    @classmethod
    def p2_hs(cls, h: float, s: float) -> float:
        """function p2_hs = p2_hs(h, s)

        Supplementary Release on Backward Equations for Pressure as a function
        of Enthalpy and Entropy p(h,s) to the IAPWS Industrial formulation 1997
        for the Thermodynamic Properties of Water and Steam

        Chapter 6: Backward Equations p(h,s) for Region 2

        :param h: enthalpy in [kJ / kg]
        :param s: specific entropy in [kJ / (kg K)]

        :return: preasure in [MPa]
        """
        if h < (-3498.98083432139 + 2575.60716905876 * s - 421.073558227969 * s**2 + 27.6349063799944 * s**3):
            sub_reg = 1
        else:
            if s < 5.85:
                sub_reg = 3
            else:
                sub_reg = 2
        if sub_reg == 1:
            # Subregion A
            # Table 6, Eq 3, page 8
            eta = h / 4200
            Sigma = s / 12
            Pi = 0
            for I, J, n in zip(R2.Sub_psh12_Table6_I, R2.Sub_psh12_Table6_J, R2.Sub_psh12_Table6_n):
                Pi = Pi + n * (eta - 0.5) ** I * (Sigma - 1.2) ** J
            p2_hs = Pi**4 * 4
        elif sub_reg == 2:
            # Subregion B
            # Table 7, Eq 4, page 9
            eta = h / 4100
            Sigma = s / 7.9
            Pi = 0
            for I, J, n in zip(R2.Sub_psh12_Table7_I, R2.Sub_psh12_Table7_J, R2.Sub_psh12_Table7_n):
                # for i in range(0, 33):
                Pi = Pi + n * (eta - 0.6) ** I * (Sigma - 1.01) ** J
            p2_hs = Pi**4 * 100
        else:
            # Subregion C
            # Table 8, Eq 5, page 10
            eta = h / 3500
            Sigma = s / 5.9
            Pi = 0
            for I, J, n in zip(R2.Sub_psh12_Table8_I, R2.Sub_psh12_Table8_J, R2.Sub_psh12_Table8_n):
                # for i in range(0, 31):
                Pi = Pi + n * (eta - 0.7) ** I * (Sigma - 1.1) ** J
            p2_hs = Pi**4 * 100
        return p2_hs

    @classmethod
    def T2_prho(cls, p: float, rho: float) -> float:
        """function T2_prho=T2_prho(p,rho)
        Solve by iteration. Observe that of low temperatures this equation has 2 solutions. Solve with half interval method

        :param p: preasure in [MPa]
        :param rho: density in [kg / m³]

        :return: temperature in [K]
        """
        Ts = float("NaN")
        if p < 16.5292:
            Low_Bound = Region4.T4_p(p)
        else:
            Low_Bound = TB23_p(p)
        High_Bound = 1073.15
        rhos = -1000
        step_counter = 0
        while math.fabs(rho - rhos) > 0.000001:
            step_counter += 1
            last_rhos = rhos

            Ts = (Low_Bound + High_Bound) / 2
            rhos = 1 / Region2.v2_pT(p, Ts)

            if last_rhos == rhos:
                cls._logger.warning(
                    "T2_prho stopped iterating after %d steps because values did not converge",
                    step_counter,
                )
                break

            if rhos < rho:
                High_Bound = Ts
            else:
                Low_Bound = Ts
        return Ts


class Region3:
    """
    Section 2.3: IAPWS IF 97 Calling functions to calculate the properties
    of water in Region 3
    """

    _logger = logging.getLogger(__name__)

    @classmethod
    def p3_rhoT(cls, rho: float, T: float) -> float:
        """function p3_rhoT = p3_rhoT(rho, T)

        7 Basic Equation for Region 3, Section. 6.1 Basic Equation

        Table 30 and 31, Page 30 and 31

        :param rho: density in [kg / m³]
        :param T: temperature in [K]

        :return: preasure in [MPa]
        """
        R = SPECIFIC_GAS_CONSTANT
        tc = CRITICAL_TEMPERATURE
        rhoc = CRITICAL_DENSITY
        delta = rho / rhoc
        tau = tc / T
        fidelta = 0
        for I, J, n in zip(R3.Table30_I, R3.Table30_J, R3.Table30_n):
            fidelta = fidelta + n * I * delta ** (I - 1) * tau**J

        fidelta = fidelta + (R3.Table30_n[0] / delta)
        return (rho * R * T * delta * fidelta) / 1000.0

    @classmethod
    def u3_rhoT(cls, rho: float, T: float) -> float:
        """function u3_rhoT = u3_rhoT(rho, T)

        7 Basic Equation for Region 3, Section. 6.1 Basic Equation

        Table 30 and 31, Page 30 and 31

        :param rho: density in [kg / m³]
        :param T: temperature in [K]

        :return: specific internal energy in [kJ / kg]
        """
        R = SPECIFIC_GAS_CONSTANT
        tc = CRITICAL_TEMPERATURE
        rhoc = CRITICAL_DENSITY
        delta = rho / rhoc
        tau = tc / T
        fitau = 0
        # TODO check table range
        for I, J, n in zip(R3.Table30_I[1:], R3.Table30_J[1:], R3.Table30_n[1:]):
            fitau = fitau + n * delta**I * J * tau ** (J - 1)
        return R * T * (tau * fitau)

    @classmethod
    def h3_rhoT(cls, rho: float, T: float) -> float:
        """function h3_rhoT = h3_rhoT(rho, T)

        7 Basic Equation for Region 3, Section. 6.1 Basic Equation

        Table 30 and 31, Page 30 and 31

        :param rho: density in [kg / m³]
        :param T: temperature in [K]

        :return: enthalpy in [kJ / kg]
        """
        R = SPECIFIC_GAS_CONSTANT
        tc = CRITICAL_TEMPERATURE
        rhoc = CRITICAL_DENSITY
        delta = rho / rhoc
        tau = tc / T
        fidelta = 0
        fitau = 0
        for I, J, n in zip(R3.Table30_I[1:], R3.Table30_J[1:], R3.Table30_n[1:]):
            fidelta = fidelta + n * I * delta ** (I - 1) * tau**J
            fitau = fitau + n * delta**I * J * tau ** (J - 1)
        fidelta = fidelta + R3.Table30_n[0] / delta
        return R * T * (tau * fitau + delta * fidelta)

    @classmethod
    def s3_rhoT(cls, rho: float, T: float) -> float:
        """function s3_rhoT = s3_rhoT(rho, T)

        7 Basic Equation for Region 3, Section. 6.1 Basic Equation

        Table 30 and 31, Page 30 and 31

        :param rho: density in [kg / m³]
        :param T: temperature in [K]

        :return: specific entropy in [kJ / (kg K)]
        """
        R = SPECIFIC_GAS_CONSTANT
        tc = CRITICAL_TEMPERATURE
        rhoc = CRITICAL_DENSITY
        delta = rho / rhoc
        tau = tc / T
        fi = 0
        fitau = 0
        for I, J, n in zip(R3.Table30_I[1:], R3.Table30_J[1:], R3.Table30_n[1:]):
            fi = fi + n * delta**I * tau**J
            fitau = fitau + n * delta**I * J * tau ** (J - 1)
        fi = fi + R3.Table30_n[0] * math.log(delta)
        return R * (tau * fitau - fi)

    @classmethod
    def Cp3_rhoT(cls, rho: float, T: float) -> float:
        """function Cp3_rhoT = Cp3_rhoT(rho, T)

        7 Basic Equation for Region 3, Section. 6.1 Basic Equation

        Table 30 and 31, Page 30 and 31

        :param rho: density in [kg / m³]
        :param T: temperature in [K]

        :return: specific isobaric heat capacity in [kJ / (kg K)]
        """
        R = SPECIFIC_GAS_CONSTANT
        tc = CRITICAL_TEMPERATURE
        rhoc = CRITICAL_DENSITY
        delta = rho / rhoc
        tau = tc / T
        fitautau = 0
        fidelta = 0
        fideltatau = 0
        fideltadelta = 0
        for I, J, n in zip(R3.Table30_I[1:], R3.Table30_J[1:], R3.Table30_n[1:]):
            fitautau = fitautau + n * delta**I * J * (J - 1) * tau ** (J - 2)
            fidelta = fidelta + n * I * delta ** (I - 1) * tau**J
            fideltatau = fideltatau + n * I * delta ** (I - 1) * J * tau ** (J - 1)
            fideltadelta = fideltadelta + n * I * (I - 1) * delta ** (I - 2) * tau**J
        fidelta = fidelta + R3.Table30_n[0] / delta
        fideltadelta = fideltadelta - R3.Table30_n[0] / (delta**2)
        return R * (
            -(tau**2) * fitautau + (delta * fidelta - delta * tau * fideltatau) ** 2 / (2 * delta * fidelta + delta**2 * fideltadelta)
        )

    @classmethod
    def Cv3_rhoT(cls, rho: float, T: float) -> float:
        """function Cv3_rhoT = Cv3_rhoT(rho, T)

        7 Basic Equation for Region 3, Section. 6.1 Basic Equation

        Table 30 and 31, Page 30 and 31

        :param rho: density in [kg / m³]
        :param T: temperature in [K]

        :return: specific isochoric heat capacity in [kJ / (kg K)]
        """
        R = SPECIFIC_GAS_CONSTANT
        tc = CRITICAL_TEMPERATURE
        rhoc = CRITICAL_DENSITY
        delta = rho / rhoc
        tau = tc / T
        fitautau = 0
        # TODO:vvvv Check for mistake vvvvv
        # for i = 1 : 40
        # IAWPS says i=2..40
        for I, J, n in zip(R3.Table30_I[1:], R3.Table30_J[1:], R3.Table30_n[1:]):
            # for i in range(1, 40):
            fitautau = fitautau + n * delta**I * J * (J - 1) * tau ** (J - 2)
        return R * -(tau * tau * fitautau)

    @classmethod
    def w3_rhoT(cls, rho: float, T: float) -> float:
        """function w3_rhoT = w3_rhoT(rho, T)

        7 Basic Equation for Region 3, Section. 6.1 Basic Equation

        Table 30 and 31, Page 30 and 31

        :param rho: density in [kg / m³]
        :param T: temperature in [K]

        :return: speed of sound in [m / s]
        """
        R = SPECIFIC_GAS_CONSTANT
        tc = CRITICAL_TEMPERATURE
        rhoc = CRITICAL_DENSITY
        delta = rho / rhoc
        tau = tc / T
        fitautau = 0
        fidelta = 0
        fideltatau = 0
        fideltadelta = 0
        for I, J, n in zip(R3.Table30_I[1:], R3.Table30_J[1:], R3.Table30_n[1:]):
            fitautau = fitautau + n * delta**I * J * (J - 1) * tau ** (J - 2)
            fidelta = fidelta + n * I * delta ** (I - 1) * tau**J
            fideltatau = fideltatau + n * I * delta ** (I - 1) * J * tau ** (J - 1)
            fideltadelta = fideltadelta + n * I * (I - 1) * delta ** (I - 2) * tau**J
        fidelta = fidelta + R3.Table30_n[0] / delta
        fideltadelta = fideltadelta - R3.Table30_n[0] / (delta**2)
        return (
            1000
            * R
            * T
            * (2 * delta * fidelta + delta**2 * fideltadelta - (delta * fidelta - delta * tau * fideltatau) ** 2 / (tau**2 * fitautau))
        ) ** 0.5

    @classmethod
    def T3_ph(cls, p: float, h: float) -> float:
        """function T3_ph = T3_ph(p, h)

        Revised Supplementary Release on Backward Equations for the functions T(p,h),
        v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial formulation
        1997 for the Thermodynamic Properties of Water and Steam 2004

        Section 3.3 Backward Equations T(p,h) and v(p,h) for Subregions 3a and 3b

        Boundary equation, Eq 1 Page 5

        :param p: preasure in [MPa]
        :param h: enthalpy in [kJ / kg]

        :return: temperature in [K]
        """
        h3ab = 2014.64004206875 + 3.74696550136983 * p - 2.19921901054187e-02 * p**2 + 8.7513168600995e-05 * p**3
        if h < h3ab:
            # Subregion 3a
            # Eq 2, Table 3, Page 7
            # TODO check table!
            Ii = [
                -12,
                -12,
                -12,
                -12,
                -12,
                -12,
                -12,
                -12,
                -10,
                -10,
                -10,
                -8,
                -8,
                -8,
                -8,
                -5,
                -3,
                -2,
                -2,
                -2,
                -1,
                -1,
                0,
                0,
                1,
                3,
                3,
                4,
                4,
                10,
                12,
            ]
            Ji = [
                0,
                1,
                2,
                6,
                14,
                16,
                20,
                22,
                1,
                5,
                12,
                0,
                2,
                4,
                10,
                2,
                0,
                1,
                3,
                4,
                0,
                2,
                0,
                1,
                1,
                0,
                1,
                0,
                3,
                4,
                5,
            ]
            ni = [
                -1.33645667811215e-07,
                4.55912656802978e-06,
                -1.46294640700979e-05,
                6.3934131297008e-03,
                372.783927268847,
                -7186.54377460447,
                573494.7521034,
                -2675693.29111439,
                -3.34066283302614e-05,
                -2.45479214069597e-02,
                47.8087847764996,
                7.64664131818904e-06,
                1.28350627676972e-03,
                1.71219081377331e-02,
                -8.51007304583213,
                -1.36513461629781e-02,
                -3.84460997596657e-06,
                3.37423807911655e-03,
                -0.551624873066791,
                0.72920227710747,
                -9.92522757376041e-03,
                -0.119308831407288,
                0.793929190615421,
                0.454270731799386,
                0.20999859125991,
                -6.42109823904738e-03,
                -0.023515586860454,
                2.52233108341612e-03,
                -7.64885133368119e-03,
                1.36176427574291e-02,
                -1.33027883575669e-02,
            ]
            ps = p / 100
            hs = h / 2300
            Ts = 0
            for i in range(0, 31):
                Ts = Ts + ni[i] * (ps + 0.24) ** Ii[i] * (hs - 0.615) ** Ji[i]
            T3_ph = Ts * 760
        else:
            # Subregion 3b
            # Eq 3, Table 4, Page 7,8
            # TODO check table!
            Ii = [
                -12,
                -12,
                -10,
                -10,
                -10,
                -10,
                -10,
                -8,
                -8,
                -8,
                -8,
                -8,
                -6,
                -6,
                -6,
                -4,
                -4,
                -3,
                -2,
                -2,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                0,
                0,
                1,
                3,
                5,
                6,
                8,
            ]
            Ji = [
                0,
                1,
                0,
                1,
                5,
                10,
                12,
                0,
                1,
                2,
                4,
                10,
                0,
                1,
                2,
                0,
                1,
                5,
                0,
                4,
                2,
                4,
                6,
                10,
                14,
                16,
                0,
                2,
                1,
                1,
                1,
                1,
                1,
            ]
            ni = [
                3.2325457364492e-05,
                -1.27575556587181e-04,
                -4.75851877356068e-04,
                1.56183014181602e-03,
                0.105724860113781,
                -85.8514221132534,
                724.140095480911,
                2.96475810273257e-03,
                -5.92721983365988e-03,
                -1.26305422818666e-02,
                -0.115716196364853,
                84.9000969739595,
                -1.08602260086615e-02,
                1.54304475328851e-02,
                7.50455441524466e-02,
                2.52520973612982e-02,
                -6.02507901232996e-02,
                -3.07622221350501,
                -5.74011959864879e-02,
                5.03471360939849,
                -0.925081888584834,
                3.91733882917546,
                -77.314600713019,
                9493.08762098587,
                -1410437.19679409,
                8491662.30819026,
                0.861095729446704,
                0.32334644281172,
                0.873281936020439,
                -0.436653048526683,
                0.286596714529479,
                -0.131778331276228,
                6.76682064330275e-03,
            ]
            hs = h / 2800
            ps = p / 100
            Ts = 0
            for i in range(0, 33):
                Ts = Ts + ni[i] * (ps + 0.298) ** Ii[i] * (hs - 0.72) ** Ji[i]
            T3_ph = Ts * 860
        return T3_ph

    @classmethod
    def v3_ph(cls, p: float, h: float) -> float:
        """function v3_ph = v3_ph(p, h)
        Revised Supplementary Release on Backward Equations for the functions T(p, h), v(p, h) and T(p, s), v(p, s) for Region 3
        of the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam 2004

        Section 3.3 Backward Equations T(p, h) and v(p, h) for Subregions 3a and 3b

        Boundary equation, Eq 1 Page 5

        :param p: preasure in [MPa]
        :param h: enthalpy in [kJ / kg]

        :return: specific volume in [m³ / kg]
        """
        h3ab = 2014.64004206875 + 3.74696550136983 * p - 2.19921901054187e-02 * p**2 + 8.7513168600995e-05 * p**3
        if h < h3ab:
            # Subregion 3a
            # Eq 4, Table 6, Page 9
            Ii = [
                -12,
                -12,
                -12,
                -12,
                -10,
                -10,
                -10,
                -8,
                -8,
                -6,
                -6,
                -6,
                -4,
                -4,
                -3,
                -2,
                -2,
                -1,
                -1,
                -1,
                -1,
                0,
                0,
                1,
                1,
                1,
                2,
                2,
                3,
                4,
                5,
                8,
            ]
            Ji = [
                6,
                8,
                12,
                18,
                4,
                7,
                10,
                5,
                12,
                3,
                4,
                22,
                2,
                3,
                7,
                3,
                16,
                0,
                1,
                2,
                3,
                0,
                1,
                0,
                1,
                2,
                0,
                2,
                0,
                2,
                2,
                2,
            ]
            ni = [
                5.29944062966028e-03,
                -0.170099690234461,
                11.1323814312927,
                -2178.98123145125,
                -5.06061827980875e-04,
                0.556495239685324,
                -9.43672726094016,
                -0.297856807561527,
                93.9353943717186,
                1.92944939465981e-02,
                0.421740664704763,
                -3689141.2628233,
                -7.37566847600639e-03,
                -0.354753242424366,
                -1.99768169338727,
                1.15456297059049,
                5683.6687581596,
                8.08169540124668e-03,
                0.172416341519307,
                1.04270175292927,
                -0.297691372792847,
                0.560394465163593,
                0.275234661176914,
                -0.148347894866012,
                -6.51142513478515e-02,
                -2.92468715386302,
                6.64876096952665e-02,
                3.52335014263844,
                -1.46340792313332e-02,
                -2.24503486668184,
                1.10533464706142,
                -4.08757344495612e-02,
            ]
            ps = p / 100
            hs = h / 2100
            vs = 0
            for i in range(0, 32):
                vs = vs + ni[i] * (ps + 0.128) ** Ii[i] * (hs - 0.727) ** Ji[i]
            v3_ph = vs * 0.0028
        else:
            # Subregion 3b
            # Eq 5, Table 7, Page 9
            Ii = [
                -12,
                -12,
                -8,
                -8,
                -8,
                -8,
                -8,
                -8,
                -6,
                -6,
                -6,
                -6,
                -6,
                -6,
                -4,
                -4,
                -4,
                -3,
                -3,
                -2,
                -2,
                -1,
                -1,
                -1,
                -1,
                0,
                1,
                1,
                2,
                2,
            ]
            Ji = [
                0,
                1,
                0,
                1,
                3,
                6,
                7,
                8,
                0,
                1,
                2,
                5,
                6,
                10,
                3,
                6,
                10,
                0,
                2,
                1,
                2,
                0,
                1,
                4,
                5,
                0,
                0,
                1,
                2,
                6,
            ]
            ni = [
                -2.25196934336318e-09,
                1.40674363313486e-08,
                2.3378408528056e-06,
                -3.31833715229001e-05,
                1.07956778514318e-03,
                -0.271382067378863,
                1.07202262490333,
                -0.853821329075382,
                -2.15214194340526e-05,
                7.6965608822273e-04,
                -4.31136580433864e-03,
                0.453342167309331,
                -0.507749535873652,
                -100.475154528389,
                -0.219201924648793,
                -3.21087965668917,
                607.567815637771,
                5.57686450685932e-04,
                0.18749904002955,
                9.05368030448107e-03,
                0.285417173048685,
                3.29924030996098e-02,
                0.239897419685483,
                4.82754995951394,
                -11.8035753702231,
                0.169490044091791,
                -1.79967222507787e-02,
                3.71810116332674e-02,
                -5.36288335065096e-02,
                1.6069710109252,
            ]
            ps = p / 100
            hs = h / 2800
            vs = 0
            for i in range(0, 30):
                vs = vs + ni[i] * (ps + 0.0661) ** Ii[i] * (hs - 0.72) ** Ji[i]
            v3_ph = vs * 0.0088
        return v3_ph

    @classmethod
    def T3_ps(cls, p: float, s: float) -> float:
        """function T3_ps = T3_ps(p, s)

        Revised Supplementary Release on Backward Equations for the functions T(p,h),
        v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial formulation
        1997 for the Thermodynamic Properties of Water and Steam 2004

        3.4 Backward Equations T(p,s) and v(p,s) for Subregions 3a and 3b

        Boundary equation, Eq 6 Page 11

        :param p: preasure in [MPa]
        :param s: specific entropy in [kJ / (kg K)]

        :return: temperature in [K]
        """
        if s <= 4.41202148223476:
            # Subregion 3a
            # Eq 6, Table 10, Page 11
            Ii = [
                -12,
                -12,
                -10,
                -10,
                -10,
                -10,
                -8,
                -8,
                -8,
                -8,
                -6,
                -6,
                -6,
                -5,
                -5,
                -5,
                -4,
                -4,
                -4,
                -2,
                -2,
                -1,
                -1,
                0,
                0,
                0,
                1,
                2,
                2,
                3,
                8,
                8,
                10,
            ]
            Ji = [
                28,
                32,
                4,
                10,
                12,
                14,
                5,
                7,
                8,
                28,
                2,
                6,
                32,
                0,
                14,
                32,
                6,
                10,
                36,
                1,
                4,
                1,
                6,
                0,
                1,
                4,
                0,
                0,
                3,
                2,
                0,
                1,
                2,
            ]
            ni = [
                1500420082.63875,
                -159397258480.424,
                5.02181140217975e-04,
                -67.2057767855466,
                1450.58545404456,
                -8238.8953488889,
                -0.154852214233853,
                11.2305046746695,
                -29.7000213482822,
                43856513263.5495,
                1.37837838635464e-03,
                -2.97478527157462,
                9717779473494.13,
                -5.71527767052398e-05,
                28830.794977842,
                -74442828926270.3,
                12.8017324848921,
                -368.275545889071,
                6.64768904779177e15,
                0.044935925195888,
                -4.22897836099655,
                -0.240614376434179,
                -4.74341365254924,
                0.72409399912611,
                0.923874349695897,
                3.99043655281015,
                3.84066651868009e-02,
                -3.59344365571848e-03,
                -0.735196448821653,
                0.188367048396131,
                1.41064266818704e-04,
                -2.57418501496337e-03,
                1.23220024851555e-03,
            ]
            Sigma = s / 4.4
            Pi = p / 100
            teta = 0
            for i in range(0, 33):
                teta = teta + ni[i] * (Pi + 0.24) ** Ii[i] * (Sigma - 0.703) ** Ji[i]
            T3_ps = teta * 760
        else:
            # Subregion 3b
            # Eq 7, Table 11, Page 11
            Ii = [
                -12,
                -12,
                -12,
                -12,
                -8,
                -8,
                -8,
                -6,
                -6,
                -6,
                -5,
                -5,
                -5,
                -5,
                -5,
                -4,
                -3,
                -3,
                -2,
                0,
                2,
                3,
                4,
                5,
                6,
                8,
                12,
                14,
            ]
            Ji = [
                1,
                3,
                4,
                7,
                0,
                1,
                3,
                0,
                2,
                4,
                0,
                1,
                2,
                4,
                6,
                12,
                1,
                6,
                2,
                0,
                1,
                1,
                0,
                24,
                0,
                3,
                1,
                2,
            ]
            ni = [
                0.52711170160166,
                -40.1317830052742,
                153.020073134484,
                -2247.99398218827,
                -0.193993484669048,
                -1.40467557893768,
                42.6799878114024,
                0.752810643416743,
                22.6657238616417,
                -622.873556909932,
                -0.660823667935396,
                0.841267087271658,
                -25.3717501764397,
                485.708963532948,
                880.531517490555,
                2650155.92794626,
                -0.359287150025783,
                -656.991567673753,
                2.41768149185367,
                0.856873461222588,
                0.655143675313458,
                -0.213535213206406,
                5.62974957606348e-03,
                -316955725450471,
                -6.99997000152457e-04,
                1.19845803210767e-02,
                1.93848122022095e-05,
                -2.15095749182309e-05,
            ]
            Sigma = s / 5.3
            Pi = p / 100
            teta = 0
            for i in range(0, 28):
                teta = teta + ni[i] * (Pi + 0.76) ** Ii[i] * (Sigma - 0.818) ** Ji[i]
            T3_ps = teta * 860
        return T3_ps

    @classmethod
    def v3_ps(cls, p: float, s: float) -> float:
        """function v3_ps = v3_ps(p, s)

        Revised Supplementary Release on Backward Equations for the functions T(p, h),
        v(p, h) and T(p, s), v(p, s) for Region 3 of the IAPWS Industrial formulation
        1997 for the Thermodynamic Properties of Water and Steam 2004

        3.4 Backward Equations T(p, s) and v(p, s) for Subregions 3a and 3b

        Boundary equation, Eq 6 Page 11

        :param p: preasure in [MPa]
        :param s: specific entropy in [kJ / (kg K)]

        :return: specific volume in [m³ / kg]
        """
        if s <= 4.41202148223476:
            # Subregion 3a
            # Eq 8, Table 13, Page 14
            Ii = [
                -12,
                -12,
                -12,
                -10,
                -10,
                -10,
                -10,
                -8,
                -8,
                -8,
                -8,
                -6,
                -5,
                -4,
                -3,
                -3,
                -2,
                -2,
                -1,
                -1,
                0,
                0,
                0,
                1,
                2,
                4,
                5,
                6,
            ]
            Ji = [
                10,
                12,
                14,
                4,
                8,
                10,
                20,
                5,
                6,
                14,
                16,
                28,
                1,
                5,
                2,
                4,
                3,
                8,
                1,
                2,
                0,
                1,
                3,
                0,
                0,
                2,
                2,
                0,
            ]
            ni = [
                79.5544074093975,
                -2382.6124298459,
                17681.3100617787,
                -1.10524727080379e-03,
                -15.3213833655326,
                297.544599376982,
                -35031520.6871242,
                0.277513761062119,
                -0.523964271036888,
                -148011.182995403,
                1600148.99374266,
                1708023226634.27,
                2.46866996006494e-04,
                1.6532608479798,
                -0.118008384666987,
                2.537986423559,
                0.965127704669424,
                -28.2172420532826,
                0.203224612353823,
                1.10648186063513,
                0.52612794845128,
                0.277000018736321,
                1.08153340501132,
                -7.44127885357893e-02,
                1.64094443541384e-02,
                -6.80468275301065e-02,
                0.025798857610164,
                -1.45749861944416e-04,
            ]
            Pi = p / 100
            Sigma = s / 4.4
            omega = 0
            for i in range(0, 28):
                omega = omega + ni[i] * (Pi + 0.187) ** Ii[i] * (Sigma - 0.755) ** Ji[i]
            v3_ps = omega * 0.0028
        else:
            # Subregion 3b
            # Eq 9, Table 14, Page 14
            Ii = [
                -12,
                -12,
                -12,
                -12,
                -12,
                -12,
                -10,
                -10,
                -10,
                -10,
                -8,
                -5,
                -5,
                -5,
                -4,
                -4,
                -4,
                -4,
                -3,
                -2,
                -2,
                -2,
                -2,
                -2,
                -2,
                0,
                0,
                0,
                1,
                1,
                2,
            ]
            Ji = [
                0,
                1,
                2,
                3,
                5,
                6,
                0,
                1,
                2,
                4,
                0,
                1,
                2,
                3,
                0,
                1,
                2,
                3,
                1,
                0,
                1,
                2,
                3,
                4,
                12,
                0,
                1,
                2,
                0,
                2,
                2,
            ]
            ni = [
                5.91599780322238e-05,
                -1.85465997137856e-03,
                1.04190510480013e-02,
                5.9864730203859e-03,
                -0.771391189901699,
                1.72549765557036,
                -4.67076079846526e-04,
                1.34533823384439e-02,
                -8.08094336805495e-02,
                0.508139374365767,
                1.28584643361683e-03,
                -1.63899353915435,
                5.86938199318063,
                -2.92466667918613,
                -6.14076301499537e-03,
                5.76199014049172,
                -12.1613320606788,
                1.67637540957944,
                -7.44135838773463,
                3.78168091437659e-02,
                4.01432203027688,
                16.0279837479185,
                3.17848779347728,
                -3.58362310304853,
                -1159952.60446827,
                0.199256573577909,
                -0.122270624794624,
                -19.1449143716586,
                -1.50448002905284e-02,
                14.6407900162154,
                -3.2747778718823,
            ]
            Pi = p / 100
            Sigma = s / 5.3
            omega = 0
            for i in range(0, 31):
                omega = omega + ni[i] * (Pi + 0.298) ** Ii[i] * (Sigma - 0.816) ** Ji[i]
            v3_ps = omega * 0.0088
        return v3_ps

    @classmethod
    def p3_hs(cls, h: float, s: float) -> float:
        """function p3_hs = p3_hs(h, s)

        Supplementary Release on Backward Equations () , p h s for Region 3,
        Equations as a function of h and s for the Region Boundaries, and an
        Equation sat, T hs for Region 4 of the IAPWS Industrial formulation 1997
        for the Thermodynamic Properties of Water and Steam 2004

        Section 3 Backward functions p(h, s), T(h, s), and v(h, s) for Region 3

        :param h: enthalpy in [kJ / kg]
        :param s: specific entropy in [kJ / (kg K)]

        :return: preasure in [MPa]
        """
        if s < 4.41202148223476:
            # Subregion 3a
            # Eq 1, Table 3, Page 8
            Sigma = s / 4.4
            eta = h / 2300
            Pi = 0
            for I, J, n in zip(R3.Sub_psh3_Table3_I, R3.Sub_psh3_Table3_J, R3.Sub_psh3_Table3_n):
                # for i in range(0, 33):
                Pi = Pi + n * (eta - 1.01) ** I * (Sigma - 0.75) ** J
            p3_hs = Pi * 99
        else:
            # Subregion 3b
            # Eq 2, Table 4, Page 8
            Sigma = s / 5.3
            eta = h / 2800
            Pi = 0
            # TODO check table range
            # for i = 1 : 35
            for I, J, n in zip(R3.Sub_psh3_Table4_I, R3.Sub_psh3_Table4_J, R3.Sub_psh3_Table4_n):
                # for i in range(0, 35):
                Pi = Pi + n * (eta - 0.681) ** I * (Sigma - 0.792) ** J
            p3_hs = 16.6 / Pi
        return p3_hs

    @classmethod
    def h3_pT(cls, p: float, T: float) -> float:
        """function h3_pT = h3_pT(p, T)

        Not available with if 97

        Solve function T3_ph - T = 0 with half interval method.

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: enthalpy in [kJ / kg]
        """
        hs = float("NaN")
        if p < CRITICAL_PRESSURE:  # Below triple point
            Ts = Region4.T4_p(p)  # Saturation temperature
            if T <= Ts:  # Liquid side
                High_Bound = Region4.h4L_p(p)  # Max h ???r liauid h.
                Low_Bound = Region1.h1_pT(p, 623.15)
            else:
                Low_Bound = Region4.h4V_p(p)  # Min h ???r Vapour h.
                High_Bound = Region2.h2_pT(p, TB23_p(p))
        else:  # Above triple point. R3 from R2 till R3.
            Low_Bound = Region1.h1_pT(p, 623.15)
            High_Bound = Region2.h2_pT(p, TB23_p(p))

        Ts = T + 1
        step_counter = 0
        while math.fabs(T - Ts) > 0.00001:
            step_counter += 1
            last_Ts = Ts

            hs = (Low_Bound + High_Bound) / 2
            Ts = Region3.T3_ph(p, hs)

            if last_Ts == Ts:
                cls._logger.warning(
                    "h3_pT stopped iterating after %d steps because values did not converge",
                    step_counter,
                )
                break

            if Ts > T:
                High_Bound = hs
            else:
                Low_Bound = hs
        return hs

    @classmethod
    def T3_prho(cls, p: float, rho: float) -> float:
        """function T3_prho = T3_prho(p, rho)

        Solve by iteration. Observe that of low temperatures this equation has
        2 solutions. Solve with half interval method

        :param p: preasure in [MPa]
        :param rho: density in [kg / m³]

        :return: temperature in [K]
        """
        Ts = float("NaN")
        Low_Bound = 623.15
        High_Bound = 1073.15
        ps = -1000
        step_counter = 0
        while math.fabs(p - ps) > 0.00000001:
            step_counter += 1
            last_ps = ps

            Ts = (Low_Bound + High_Bound) / 2
            ps = Region3.p3_rhoT(rho, Ts)

            if last_ps == ps:
                cls._logger.warning(
                    "T3_prho stopped iterating after %d steps because values did not converge",
                    step_counter,
                )
                break

            if ps > p:
                High_Bound = Ts
            else:
                Low_Bound = Ts
        return Ts

    @staticmethod
    def psat3_h(h: float) -> float:
        """
        calculate saturation preasure from preasure for region 3

        Section 4.2 Region 3. pSat_h  & pSat_s

        Revised Supplementary Release on Backward Equations for the functions T(p,h), v(p,h) s& T(p,s), v(p,s) for Region 3
        of the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water & Steam 2004 Section 4 Boundary
        Equations psat(h) & psat(s) for the Saturation Lines of Region 3

        See pictures Page 17, Eq 10, Table 17, Page 18

        :param h: enthalpy in [kJ / kg]

        :return: saturation preasure in [MPa]
        """
        Ii = [0, 1, 1, 1, 1, 5, 7, 8, 14, 20, 22, 24, 28, 36]
        Ji = [0, 1, 3, 4, 36, 3, 0, 24, 16, 16, 3, 18, 8, 24]
        ni = [
            0.600073641753024,
            -9.36203654849857,
            24.6590798594147,
            -107.014222858224,
            -91582131580576.8,
            -8623.32011700662,
            -23.5837344740032,
            2.52304969384128e17,
            -3.89718771997719e18,
            -3.33775713645296e22,
            35649946963.6328,
            -1.48547544720641e26,
            3.30611514838798e18,
            8.13641294467829e37,
        ]
        hs = h / 2600
        ps = 0
        for i in range(0, 14):
            ps = ps + ni[i] * (hs - 1.02) ** Ii[i] * (hs - 0.608) ** Ji[i]
        return ps * 22

    @staticmethod
    def psat3_s(s: float) -> float:
        """
        calculate saturation preasure from Specific entropy for region 3

        Section 4.2 Region 3. pSat_h  & pSat_s

        :param s: Specific entropy in [kJ / (kg K)]

        :return: saturation preasure in [MPa]
        """
        Ii = [0, 1, 1, 4, 12, 12, 16, 24, 28, 32]
        Ji = [0, 1, 32, 7, 4, 14, 36, 10, 0, 18]
        ni = [
            0.639767553612785,
            -12.9727445396014,
            -2.24595125848403e15,
            1774667.41801846,
            7170793495.71538,
            -3.78829107169011e17,
            -9.55586736431328e34,
            1.87269814676188e23,
            119254746466.473,
            1.10649277244882e36,
        ]
        Sigma = s / 5.2
        Pi = 0
        for i in range(0, 10):
            Pi = Pi + ni[i] * (Sigma - 1.03) ** Ii[i] * (Sigma - 0.699) ** Ji[i]
        return Pi * 22


class Region4:
    """
    Section 2.4: IAPWS IF 97 Calling functions to calculate the properties of
    water in Region 4
    """

    _logger = logging.getLogger(__name__)

    @classmethod
    def p4_T(cls, T: float) -> float:
        """function p4_T = p4_T(T)

        Section 8.1 The Saturation-Pressure Equation

        Eq 30, Page 33

        :param T: temperature in [K]

        :return: preasure in [MPa]
        """
        teta = T - 0.23855557567849 / (T - 650.17534844798)
        a = teta**2 + 1167.0521452767 * teta - 724213.16703206
        B = -17.073846940092 * teta**2 + 12020.82470247 * teta - 3232555.0322333
        C = 14.91510861353 * teta**2 - 4823.2657361591 * teta + 405113.40542057
        return (2 * C / (-B + (B**2 - 4 * a * C) ** 0.5)) ** 4

    @classmethod
    def T4_p(cls, p: float) -> float:
        """function T4_p = T4_p(p)

        Section 8.2 The Saturation-Temperature Equation

        Eq 31, Page 34

        :param p: preasure in [MPa]

        :return: temperature in [K]
        """
        beta = p**0.25
        E = beta**2 - 17.073846940092 * beta + 14.91510861353
        f = 1167.0521452767 * beta**2 + 12020.82470247 * beta - 4823.2657361591
        G = -724213.16703206 * beta**2 - 3232555.0322333 * beta + 405113.40542057
        D = 2 * G / (-f - (f**2 - 4 * E * G) ** 0.5)
        return (650.17534844798 + D - ((650.17534844798 + D) ** 2 - 4 * (-0.23855557567849 + 650.17534844798 * D)) ** 0.5) / 2

    @classmethod
    def h4_s(cls, s: float) -> float:
        """function h4_s = h4_s(s)

        Supplementary Release on Backward Equations () , p h s for Region 3, Equations
        as a function of h and s for the Region Boundaries, and an Equation() sat,
        T hs for Region 4 of the IAPWS Industrial formulation 1997 for the
        Thermodynamic Properties of Water and Steam 4 Equations for Region Boundaries
        Given Enthalpy and Entropy

        See picture page 14

        :param s: specific entropy in [kJ / (kg K)]

        :return: enthalpy in [kJ / kg]
        """
        if -0.0001545495919 < s <= 3.77828134:
            # hL1_s
            # Eq 3, Table 9, Page 16
            Sigma = s / 3.8
            eta = 0
            for I, J, n in zip(R3.Sub_psh3_Table9_I, R3.Sub_psh3_Table9_J, R3.Sub_psh3_Table9_n):
                eta = eta + n * (Sigma - 1.09) ** I * (Sigma + 0.0000366) ** J
            h4_s = eta * 1700
        elif 3.77828134 < s <= 4.41202148223476:
            # hL3_s
            # Eq 4, Table 10, Page 16
            Sigma = s / 3.8
            eta = 0
            for I, J, n in zip(R3.Sub_psh3_Table10_I, R3.Sub_psh3_Table10_J, R3.Sub_psh3_Table10_n):
                eta = eta + n * (Sigma - 1.09) ** I * (Sigma + 0.0000366) ** J
            h4_s = eta * 1700
        elif 4.41202148223476 < s <= 5.85:
            # Section 4.4 Equations () 2ab " h s and ( ) 2c3b "h s for the
            # Saturated Vapor Line
            # Page 19, Eq 5
            # hV2c3b_s(s)
            Sigma = s / 5.9
            eta = 0
            for I, J, n in zip(R3.Sub_psh3_Table17_I, R3.Sub_psh3_Table17_J, R3.Sub_psh3_Table17_n):
                eta = eta + n * (Sigma - 1.02) ** I * (Sigma - 0.726) ** J
            h4_s = eta**4 * 2800
        elif 5.85 < s < 9.155759395:
            # Section 4.4 Equations () 2ab " h s and ( ) 2c3b "h s for the
            # Saturated Vapor Line
            # Page 20, Eq 6
            Sigma1 = s / 5.21
            Sigma2 = s / 9.2
            eta = 0
            for I, J, n in zip(R3.Sub_psh3_Table16_I, R3.Sub_psh3_Table16_J, R3.Sub_psh3_Table16_n):
                eta = eta + n * (1 / Sigma1 - 0.513) ** I * (Sigma2 - 0.524) ** J
            h4_s = math.exp(eta) * 2800
        else:
            h4_s = -99999
        return h4_s

    @classmethod
    def p4_s(cls, s: float) -> float:
        """function p4_s = p4_s(s)

        Uses h4_s and p_hs for the different regions to determine p4_s

        :param s: specific entropy in [kJ / (kg K)]

        :return: preasure in [MPa]
        """
        h_sat = Region4.h4_s(s)
        if -0.0001545495919 < s <= 3.77828134:
            p4_s = Region1.p1_hs(h_sat, s)
        elif 3.77828134 < s <= 5.210887663:
            p4_s = Region3.p3_hs(h_sat, s)
        elif 5.210887663 < s < 9.155759395:
            p4_s = Region2.p2_hs(h_sat, s)
        else:
            p4_s = -99999
        return p4_s

    @classmethod
    def h4L_p(cls, p: float) -> float:
        """function h4L_p = h4L_p(p)

        :param p: preasure in [MPa]

        :return: enthalpy in [kJ / kg]
        """
        hs = float("NaN")
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            Ts = Region4.T4_p(p)
            if p < 16.529:
                h4L_p = Region1.h1_pT(p, Ts)
            else:
                # Iterate to find the the backward solution of p3sat_h
                Low_Bound = 1670.858218
                High_Bound = 2087.23500164864
                ps = -1000
                step_counter = 0
                while math.fabs(p - ps) > 0.00001:
                    step_counter += 1
                    last_ps = ps

                    hs = (Low_Bound + High_Bound) / 2
                    ps = Region3.psat3_h(hs)

                    if last_ps == ps:
                        cls._logger.warning(
                            "h4L_p stopped iterating after %d steps because values did not converge",
                            step_counter,
                        )
                        break

                    if ps > p:
                        High_Bound = hs
                    else:
                        Low_Bound = hs
                h4L_p = hs
        else:
            h4L_p = -99999
        return h4L_p

    @classmethod
    def h4V_p(cls, p: float) -> float:
        """function h4V_p = h4V_p(p)

        :param p: preasure in [MPa]

        :return: enthalpy in [kJ / kg]
        """
        hs = float("NaN")
        if TRIPLE_POINT_PRESSURE < p < CRITICAL_PRESSURE:
            Ts = Region4.T4_p(p)
            if p < 16.529:
                h4V_p = Region2.h2_pT(p, Ts)
            else:
                # Iterate to find the the backward solution of p3sat_h
                Low_Bound = 2087.23500164864
                High_Bound = 2563.592004 + 5
                ps = -1000
                step_counter = 0
                while math.fabs(p - ps) > 0.000001:
                    step_counter += 1
                    last_ps = ps

                    hs = (Low_Bound + High_Bound) / 2
                    ps = Region3.psat3_h(hs)

                    if last_ps == ps:
                        cls._logger.warning(
                            "h4V_p stopped iterating after %d steps because values did not converge",
                            step_counter,
                        )
                        break

                    if ps < p:
                        High_Bound = hs
                    else:
                        Low_Bound = hs
                h4V_p = hs
        else:
            h4V_p = -99999
        return h4V_p

    @classmethod
    def x4_ph(cls, p: float, h: float) -> float:
        """function x4_ph = x4_ph(p, h)

        Calculate vapour fraction from hL and hV for given p

        :param p: preasure in [MPa]
        :param h: enthalpy in [kJ / kg]

        :return: vapor fraction
        """
        h4v = Region4.h4V_p(p)
        h4L = Region4.h4L_p(p)
        if h > h4v:
            x4_ph = 1
        elif h < h4L:
            x4_ph = 0
        else:
            x4_ph = (h - h4L) / (h4v - h4L)
        return x4_ph

    @classmethod
    def x4_ps(cls, p: float, s: float) -> float:
        """function x4_ps = x4_ps(p, s)

        :param p: preasure in [MPa]
        :param s: specific entropy in [kJ / (kg K)]

        :return: vapor fraction
        """
        if p < 16.529:
            ssv = Region2.s2_pT(p, Region4.T4_p(p))
            ssL = Region1.s1_pT(p, Region4.T4_p(p))
        else:
            ssv = Region3.s3_rhoT(1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p))
            ssL = Region3.s3_rhoT(1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p))
        if s < ssL:
            x4_ps = 0
        elif s > ssv:
            x4_ps = 1
        else:
            x4_ps = (s - ssL) / (ssv - ssL)
        return x4_ps

    @classmethod
    def T4_hs(cls, h: float, s: float) -> float:
        """function T4_hs = T4_hs(h, s)

        Supplementary Release on Backward Equations ( ) , p h s for Region 3,

        Chapter 5.3 page 30.

        :param h: enthalpy in [kJ / kg]
        :param s: specific entropy in [kJ / (kg K)]

        :return: temperature in [K]
        """
        PL = float("NaN")
        p = float("NaN")
        Ii = R3.Sub_psh3_Table28_I
        Ji = R3.Sub_psh3_Table28_J
        ni = R3.Sub_psh3_Table28_n
        if 5.210887825 < s < 9.15546555571324:
            Sigma = s / 9.2
            eta = h / 2800
            teta = 0
            for i in range(0, 36):
                teta = teta + ni[i] * (eta - 0.119) ** Ii[i] * (Sigma - 1.07) ** Ji[i]
            T4_hs = teta * 550
        else:
            # function psat_h
            if -0.0001545495919 < s <= 3.77828134:
                Low_Bound = 0.000611
                High_Bound = 165.291642526045
                hL = -1000
                while (math.fabs(hL - h) > 0.00001) and (math.fabs(High_Bound - Low_Bound) > 0.0001):
                    PL = (Low_Bound + High_Bound) / 2
                    Ts = Region4.T4_p(PL)
                    hL = Region1.h1_pT(PL, Ts)
                    if hL > h:
                        High_Bound = PL
                    else:
                        Low_Bound = PL
            elif 3.77828134 < s <= 4.41202148223476:
                PL = Region3.psat3_h(h)
            elif 4.41202148223476 < s <= 5.210887663:
                PL = Region3.psat3_h(h)
            Low_Bound = 0.000611
            High_Bound = PL
            sss = -1000
            while (math.fabs(s - sss) > 0.000001) and (math.fabs(High_Bound - Low_Bound) > 0.0000001):
                p = (Low_Bound + High_Bound) / 2
                # Calculate s4_ph
                Ts = Region4.T4_p(p)
                xs = Region4.x4_ph(p, h)
                if p < 16.529:
                    s4v = Region2.s2_pT(p, Ts)
                    s4L = Region1.s1_pT(p, Ts)
                else:
                    v4v = Region3.v3_ph(p, Region4.h4V_p(p))
                    s4v = Region3.s3_rhoT(1 / v4v, Ts)
                    v4L = Region3.v3_ph(p, Region4.h4L_p(p))
                    s4L = Region3.s3_rhoT(1 / v4L, Ts)
                sss = xs * s4v + (1 - xs) * s4L
                if sss < s:
                    High_Bound = p
                else:
                    Low_Bound = p
            T4_hs = Region4.T4_p(p)
        return T4_hs


class Region5:
    """
    Section 2.5: IAPWS IF 97 Calling functions to calculate the properties
    of water in Region 5
    """

    nir = [
        -1.2563183589592e-04,
        2.1774678714571e-03,
        -0.004594282089991,
        -3.9724828359569e-06,
        1.2919228289784e-07,
    ]

    Iir = [1, 1, 1, 2, 3]
    Jir = [0, 1, 3, 9, 3]

    _logger = logging.getLogger(__name__)

    @classmethod
    def h5_pT(cls, p: float, T: float) -> float:
        """function h5_pT = h5_pT(p, T)

        Basic Equation for Region 5

        Eq 32,33, Page 36, Tables 37-41

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: enthalpy in [kJ / kg]
        """
        Iir = Region5.Iir
        Jir = Region5.Jir
        nir = Region5.nir
        R = SPECIFIC_GAS_CONSTANT
        tau = 1000 / T
        Pi = p

        gamma0_tau = 0
        for J, n in zip(R4.Table37_J, R4.Table37_n):
            gamma0_tau = gamma0_tau + n * J * tau ** (J - 1)

        gammar_tau = 0
        for i in range(0, 5):
            gammar_tau = gammar_tau + nir[i] * Pi ** Iir[i] * Jir[i] * tau ** (Jir[i] - 1)
        return R * T * tau * (gamma0_tau + gammar_tau)

    @classmethod
    def v5_pT(cls, p: float, T: float) -> float:
        """function v5_pT = v5_pT(p, T)

        Basic Equation for Region 5

        Eq 32,33, Page 36, Tables 37-41

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific volume in [m³ / kg]
        """
        Iir = Region5.Iir
        Jir = Region5.Jir
        nir = Region5.nir
        R = SPECIFIC_GAS_CONSTANT
        tau = 1000 / T  #
        Pi = p  #
        gamma0_pi = 1 / Pi  #
        gammar_pi = 0  #
        for i in range(0, 5):
            gammar_pi = gammar_pi + nir[i] * Iir[i] * Pi ** (Iir[i] - 1) * tau ** Jir[i]
        return R * T / p * Pi * (gamma0_pi + gammar_pi) / 1000

    @classmethod
    def u5_pT(cls, p: float, T: float) -> float:
        """function u5_pT = u5_pT(p, T)

        Basic Equation for Region 5

        Eq 32,33, Page 36, Tables 37-41

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific internal energy in [kJ / kg]
        """
        Ji0 = R4.Table37_J
        ni0 = R4.Table37_n
        Iir = Region5.Iir
        Jir = Region5.Jir
        nir = Region5.nir
        R = SPECIFIC_GAS_CONSTANT
        tau = 1000 / T
        Pi = p

        gamma0_pi = 1 / Pi

        gamma0_tau = 0
        for J, n in zip(R4.Table37_J, R4.Table37_n):
            gamma0_tau = gamma0_tau + n * J * tau ** (J - 1)

        gammar_pi = 0
        gammar_tau = 0
        for i in range(0, 5):
            gammar_pi = gammar_pi + nir[i] * Iir[i] * Pi ** (Iir[i] - 1) * tau ** Jir[i]
            gammar_tau = gammar_tau + nir[i] * Pi ** Iir[i] * Jir[i] * tau ** (Jir[i] - 1)
        return R * T * (tau * (gamma0_tau + gammar_tau) - Pi * (gamma0_pi + gammar_pi))

    @classmethod
    def Cp5_pT(cls, p: float, T: float) -> float:
        """function Cp5_pT = Cp5_pT(p, T)

        Basic Equation for Region 5

        Eq 32,33, Page 36, Tables 37-41

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific isobaric heat capacity in [kJ / (kg K)]
        """
        Iir = Region5.Iir
        Jir = Region5.Jir
        nir = Region5.nir
        R = SPECIFIC_GAS_CONSTANT
        tau = 1000 / T
        Pi = p

        gamma0_tautau = 0
        for J, n in zip(R4.Table37_J, R4.Table37_n):
            gamma0_tautau = gamma0_tautau + n * J * (J - 1) * tau ** (J - 2)

        gammar_tautau = 0
        for i in range(0, 5):
            gammar_tautau = gammar_tautau + nir[i] * Pi ** Iir[i] * Jir[i] * (Jir[i] - 1) * tau ** (Jir[i] - 2)
        return -R * tau**2 * (gamma0_tautau + gammar_tautau)

    @classmethod
    def s5_pT(cls, p: float, T: float) -> float:
        """function s5_pT = s5_pT(p, T)

        Basic Equation for Region 5

        Eq 32, 33, Page 36, Tables 37 - 41

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific entropy in [kJ / (kg K)]
        """
        Iir = Region5.Iir
        Jir = Region5.Jir
        nir = Region5.nir
        R = SPECIFIC_GAS_CONSTANT
        tau = 1000 / T
        Pi = p

        gamma0 = math.log(Pi)
        gamma0_tau = 0
        for J, n in zip(R4.Table37_J, R4.Table37_n):
            gamma0 = gamma0 + n * tau**J
            gamma0_tau = gamma0_tau + n * J * tau ** (J - 1)

        gammar = 0
        gammar_tau = 0
        for i in range(0, 5):
            gammar = gammar + nir[i] * Pi ** Iir[i] * tau ** Jir[i]
            gammar_tau = gammar_tau + nir[i] * Pi ** Iir[i] * Jir[i] * tau ** (Jir[i] - 1)
        return R * (tau * (gamma0_tau + gammar_tau) - (gamma0 + gammar))

    @classmethod
    def Cv5_pT(cls, p: float, T: float) -> float:
        """function Cv5_pT = Cv5_pT(p, T)

        Basic Equation for Region 5

        Eq 32, 33, Page 36, Tables 37 - 41

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific isochoric heat capacity in [kJ / (kg K)]
        """
        Ji0 = R4.Table37_J
        ni0 = R4.Table37_n
        Iir = Region5.Iir
        Jir = Region5.Jir
        nir = Region5.nir
        R = SPECIFIC_GAS_CONSTANT
        tau = 1000 / T
        Pi = p

        gamma0_tautau = 0
        for J, n in zip(R4.Table37_J, R4.Table37_n):
            gamma0_tautau = gamma0_tautau + n * (J - 1) * J * tau ** (J - 2)

        gammar_pi = 0
        gammar_pitau = 0
        gammar_pipi = 0
        gammar_tautau = 0
        for i in range(0, 5):
            gammar_pi = gammar_pi + nir[i] * Iir[i] * Pi ** (Iir[i] - 1) * tau ** Jir[i]
            gammar_pitau = gammar_pitau + nir[i] * Iir[i] * Pi ** (Iir[i] - 1) * Jir[i] * tau ** (Jir[i] - 1)
            gammar_pipi = gammar_pipi + nir[i] * Iir[i] * (Iir[i] - 1) * Pi ** (Iir[i] - 2) * tau ** Jir[i]
            gammar_tautau = gammar_tautau + nir[i] * Pi ** Iir[i] * Jir[i] * (Jir[i] - 1) * tau ** (Jir[i] - 2)
        return R * (
            -(tau**2 * (gamma0_tautau + gammar_tautau)) - (1 + Pi * gammar_pi - tau * Pi * gammar_pitau) ** 2 / (1 - Pi**2 * gammar_pipi)
        )

    @classmethod
    def w5_pT(cls, p: float, T: float) -> float:
        """function w5_pT = w5_pT(p, T)

        Basic Equation for Region 5

        Eq 32, 33, Page 36, Tables 37 - 41

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: speed of sound in [m / s]
        """
        Iir = Region5.Iir
        Jir = Region5.Jir
        nir = Region5.nir
        R = SPECIFIC_GAS_CONSTANT
        tau = 1000 / T
        Pi = p

        gamma0_tautau = 0
        for J, n in zip(R4.Table37_J, R4.Table37_n):
            gamma0_tautau = gamma0_tautau + n * (J - 1) * J * tau ** (J - 2)

        gammar_pi = 0
        gammar_pitau = 0
        gammar_pipi = 0
        gammar_tautau = 0
        for i in range(0, 5):
            gammar_pi = gammar_pi + nir[i] * Iir[i] * Pi ** (Iir[i] - 1) * tau ** Jir[i]
            gammar_pitau = gammar_pitau + nir[i] * Iir[i] * Pi ** (Iir[i] - 1) * Jir[i] * tau ** (Jir[i] - 1)
            gammar_pipi = gammar_pipi + nir[i] * Iir[i] * (Iir[i] - 1) * Pi ** (Iir[i] - 2) * tau ** Jir[i]
            gammar_tautau = gammar_tautau + nir[i] * Pi ** Iir[i] * Jir[i] * (Jir[i] - 1) * tau ** (Jir[i] - 2)

        return (
            1000
            * R
            * T
            * (1 + 2 * Pi * gammar_pi + Pi**2 * gammar_pi**2)
            / ((1 - Pi**2 * gammar_pipi) + (1 + Pi * gammar_pi - tau * Pi * gammar_pitau) ** 2 / (tau**2 * (gamma0_tautau + gammar_tautau)))
        ) ** 0.5

    @classmethod
    def T5_ph(cls, p: float, h: float) -> float:
        """function T5_ph = T5_ph(p, h)

        Solve with half interval method

        :param p: preasure in [MPa]
        :param h: enthalpy in [kJ / kg]

        :return: temperature in [K]
        """
        Ts = float("NaN")
        Low_Bound = 1073.15
        High_Bound = 2273.15
        hs = h - 1
        step_counter = 0
        while math.fabs(h - hs) > 0.00001:
            step_counter += 1
            last_hs = hs

            Ts = (Low_Bound + High_Bound) / 2
            hs = Region5.h5_pT(p, Ts)

            if last_hs == hs:
                cls._logger.warning(
                    "T5_ph stopped iterating after %d steps because values did not converge",
                    step_counter,
                )
                break

            if hs > h:
                High_Bound = Ts
            else:
                Low_Bound = Ts
        return Ts

    @classmethod
    def T5_ps(cls, p: float, s: float) -> float:
        """function T5_ps = T5_ps(p, s)

        Solve with half interval method

        :param p: preasure in [MPa]
        :param s: specific entropy in [kJ / (kg K)]

        :return: temperature in [K]
        """
        Ts = float("NaN")
        Low_Bound = 1073.15
        High_Bound = 2273.15
        ss = s - 1
        step_counter = 0
        while math.fabs(s - ss) > 0.00001:
            step_counter += 1
            last_ss = ss

            Ts = (Low_Bound + High_Bound) / 2
            ss = Region5.s5_pT(p, Ts)

            if last_ss == ss:
                cls._logger.warning(
                    "T5_ps stopped iterating after %d steps because values did not converge",
                    step_counter,
                )
                break

            if ss > s:
                High_Bound = Ts
            else:
                Low_Bound = Ts
        return Ts

    @classmethod
    def T5_prho(cls, p: float, rho: float) -> float:
        """function T5_prho = T5_prho(p, rho)

        Solve by iteration. Observe that for low temperatures this equation has
        2 solutions. Solve with half interval method

        :param p: preasure in [MPa]
        :param rho: density in [kg / m³]

        :return: temperature in [K]
        """
        Ts = float("NaN")
        Low_Bound = 1073.15
        High_Bound = 2073.15
        rhos = -1000
        step_counter = 0
        while math.fabs(rho - rhos) > 0.000001:
            step_counter += 1
            last_rhos = rhos

            Ts = (Low_Bound + High_Bound) / 2
            rhos = 1 / Region2.v2_pT(p, Ts)

            if last_rhos == rhos:
                cls._logger.warning(
                    "T5_prho stopped iterating after %d steps because values did not converge",
                    step_counter,
                )
                break

            if rhos < rho:
                High_Bound = Ts
            else:
                Low_Bound = Ts
        return Ts
