#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Section 2: IAPWS IF 97 Calling functions
"""
import math
import logging
from .RegionBorders import B23T_p, p3sat_h
from .Constants import (
    SPECIFIC_GAS_CONSTANT,
    CRITICAL_TEMPERATURE,
    CRITICAL_PRESSURE,
    TRIPLE_POINT_PRESSURE,
    FREEZING_TEMPERATURE_H2O,
)


class Region1:
    """
    Section 2.1: IAPWS IF 97 Calling functions to calculate the properties of water in Region 1

    Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
    """

    @staticmethod
    def v1_pT(p: float, T: float) -> float:
        """function v1_pT = v1_pT(p, T)

        5 Equations for Region 1, Section. 5.1 Basic Equation

        Equation 7, Table 3, Page 6

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific volume in [m³ / kg]
        """
        I1 = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
            8,
            8,
            21,
            23,
            29,
            30,
            31,
            32,
        ]
        J1 = [
            -2,
            -1,
            0,
            1,
            2,
            3,
            4,
            5,
            -9,
            -7,
            -1,
            0,
            1,
            3,
            -3,
            0,
            1,
            3,
            17,
            -4,
            0,
            6,
            -5,
            -2,
            10,
            -8,
            -11,
            -6,
            -29,
            -31,
            -38,
            -39,
            -40,
            -41,
        ]
        n1 = [
            0.14632971213167,
            -0.84548187169114,
            -3.756360367204,
            3.3855169168385,
            -0.95791963387872,
            0.15772038513228,
            -0.016616417199501,
            8.1214629983568e-04,
            2.8319080123804e-04,
            -6.0706301565874e-04,
            -0.018990068218419,
            -0.032529748770505,
            -0.021841717175414,
            -5.283835796993e-05,
            -4.7184321073267e-04,
            -3.0001780793026e-04,
            4.7661393906987e-05,
            -4.4141845330846e-06,
            -7.2694996297594e-16,
            -3.1679644845054e-05,
            -2.8270797985312e-06,
            -8.5205128120103e-10,
            -2.2425281908e-06,
            -6.5171222895601e-07,
            -1.4341729937924e-13,
            -4.0516996860117e-07,
            -1.2734301741641e-09,
            -1.7424871230634e-10,
            -6.8762131295531e-19,
            1.4478307828521e-20,
            2.6335781662795e-23,
            -1.1947622640071e-23,
            1.8228094581404e-24,
            -9.3537087292458e-26,
        ]
        R = SPECIFIC_GAS_CONSTANT
        Pi = p / 16.53  # 16.53 MPa
        tau = 1386 / T  # 1386 K
        gamma_der_pi = 0
        for i in range(0, 34):
            gamma_der_pi = (
                gamma_der_pi
                - n1[i] * I1[i] * (7.1 - Pi) ** (I1[i] - 1) * (tau - 1.222) ** J1[i]
            )
        return R * T / p * Pi * gamma_der_pi / 1000

    @staticmethod
    def h1_pT(p: float, T: float) -> float:
        """
        function h1_pT = h1_pT(p, T)

        5 Equations for Region 1, Section. 5.1 Basic Equation

        Equation 7, Table 3, Page 6

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: enthalpy in [kJ / kg]
        """
        I1 = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
            8,
            8,
            21,
            23,
            29,
            30,
            31,
            32,
        ]
        J1 = [
            -2,
            -1,
            0,
            1,
            2,
            3,
            4,
            5,
            -9,
            -7,
            -1,
            0,
            1,
            3,
            -3,
            0,
            1,
            3,
            17,
            -4,
            0,
            6,
            -5,
            -2,
            10,
            -8,
            -11,
            -6,
            -29,
            -31,
            -38,
            -39,
            -40,
            -41,
        ]
        n1 = [
            0.14632971213167,
            -0.84548187169114,
            -3.756360367204,
            3.3855169168385,
            -0.95791963387872,
            0.15772038513228,
            -0.016616417199501,
            8.1214629983568e-04,
            2.8319080123804e-04,
            -6.0706301565874e-04,
            -0.018990068218419,
            -0.032529748770505,
            -0.021841717175414,
            -5.283835796993e-05,
            -4.7184321073267e-04,
            -3.0001780793026e-04,
            4.7661393906987e-05,
            -4.4141845330846e-06,
            -7.2694996297594e-16,
            -3.1679644845054e-05,
            -2.8270797985312e-06,
            -8.5205128120103e-10,
            -2.2425281908e-06,
            -6.5171222895601e-07,
            -1.4341729937924e-13,
            -4.0516996860117e-07,
            -1.2734301741641e-09,
            -1.7424871230634e-10,
            -6.8762131295531e-19,
            1.4478307828521e-20,
            2.6335781662795e-23,
            -1.1947622640071e-23,
            1.8228094581404e-24,
            -9.3537087292458e-26,
        ]
        R = SPECIFIC_GAS_CONSTANT
        Pi = p / 16.53
        tau = 1386 / T
        gamma_der_tau = 0
        for i in range(0, 34):
            gamma_der_tau = gamma_der_tau + (
                n1[i] * (7.1 - Pi) ** I1[i] * J1[i] * (tau - 1.222) ** (J1[i] - 1)
            )
        return R * T * tau * gamma_der_tau

    @staticmethod
    def u1_pT(p: float, T: float) -> float:
        """function u1_pT = u1_pT(p, T)

        5 Equations for Region 1, Section. 5.1 Basic Equation

        Equation 7, Table 3, Page 6

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific internal energy in [kJ / kg]
        """
        I1 = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
            8,
            8,
            21,
            23,
            29,
            30,
            31,
            32,
        ]
        J1 = [
            -2,
            -1,
            0,
            1,
            2,
            3,
            4,
            5,
            -9,
            -7,
            -1,
            0,
            1,
            3,
            -3,
            0,
            1,
            3,
            17,
            -4,
            0,
            6,
            -5,
            -2,
            10,
            -8,
            -11,
            -6,
            -29,
            -31,
            -38,
            -39,
            -40,
            -41,
        ]
        n1 = [
            0.14632971213167,
            -0.84548187169114,
            -3.756360367204,
            3.3855169168385,
            -0.95791963387872,
            0.15772038513228,
            -0.016616417199501,
            8.1214629983568e-04,
            2.8319080123804e-04,
            -6.0706301565874e-04,
            -0.018990068218419,
            -0.032529748770505,
            -0.021841717175414,
            -5.283835796993e-05,
            -4.7184321073267e-04,
            -3.0001780793026e-04,
            4.7661393906987e-05,
            -4.4141845330846e-06,
            -7.2694996297594e-16,
            -3.1679644845054e-05,
            -2.8270797985312e-06,
            -8.5205128120103e-10,
            -2.2425281908e-06,
            -6.5171222895601e-07,
            -1.4341729937924e-13,
            -4.0516996860117e-07,
            -1.2734301741641e-09,
            -1.7424871230634e-10,
            -6.8762131295531e-19,
            1.4478307828521e-20,
            2.6335781662795e-23,
            -1.1947622640071e-23,
            1.8228094581404e-24,
            -9.3537087292458e-26,
        ]
        R = SPECIFIC_GAS_CONSTANT
        Pi = p / 16.53
        tau = 1386 / T
        gamma_der_tau = 0
        gamma_der_pi = 0
        for i in range(0, 34):
            gamma_der_pi = (
                gamma_der_pi
                - n1[i] * I1[i] * (7.1 - Pi) ** (I1[i] - 1) * (tau - 1.222) ** J1[i]
            )
            gamma_der_tau = gamma_der_tau + (
                n1[i] * (7.1 - Pi) ** I1[i] * J1[i] * (tau - 1.222) ** (J1[i] - 1)
            )
        return R * T * (tau * gamma_der_tau - Pi * gamma_der_pi)

    @staticmethod
    def s1_pT(p: float, T: float) -> float:
        """function s1_pT = s1_pT(p, T)

        5 Equations for Region 1, Section. 5.1 Basic Equation

        Equation 7, Table 3, Page 6

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific entropy in [kJ / (kg K)]
        """
        I1 = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
            8,
            8,
            21,
            23,
            29,
            30,
            31,
            32,
        ]
        J1 = [
            -2,
            -1,
            0,
            1,
            2,
            3,
            4,
            5,
            -9,
            -7,
            -1,
            0,
            1,
            3,
            -3,
            0,
            1,
            3,
            17,
            -4,
            0,
            6,
            -5,
            -2,
            10,
            -8,
            -11,
            -6,
            -29,
            -31,
            -38,
            -39,
            -40,
            -41,
        ]
        n1 = [
            0.14632971213167,
            -0.84548187169114,
            -3.756360367204,
            3.3855169168385,
            -0.95791963387872,
            0.15772038513228,
            -0.016616417199501,
            8.1214629983568e-04,
            2.8319080123804e-04,
            -6.0706301565874e-04,
            -0.018990068218419,
            -0.032529748770505,
            -0.021841717175414,
            -5.283835796993e-05,
            -4.7184321073267e-04,
            -3.0001780793026e-04,
            4.7661393906987e-05,
            -4.4141845330846e-06,
            -7.2694996297594e-16,
            -3.1679644845054e-05,
            -2.8270797985312e-06,
            -8.5205128120103e-10,
            -2.2425281908e-06,
            -6.5171222895601e-07,
            -1.4341729937924e-13,
            -4.0516996860117e-07,
            -1.2734301741641e-09,
            -1.7424871230634e-10,
            -6.8762131295531e-19,
            1.4478307828521e-20,
            2.6335781662795e-23,
            -1.1947622640071e-23,
            1.8228094581404e-24,
            -9.3537087292458e-26,
        ]
        R = SPECIFIC_GAS_CONSTANT
        Pi = p / 16.53
        tau = 1386 / T
        gamma = 0
        gamma_der_tau = 0
        for i in range(0, 34):
            gamma_der_tau = gamma_der_tau + (
                n1[i] * (7.1 - Pi) ** I1[i] * J1[i] * (tau - 1.222) ** (J1[i] - 1)
            )
            gamma = gamma + n1[i] * (7.1 - Pi) ** I1[i] * (tau - 1.222) ** J1[i]
        return R * tau * gamma_der_tau - R * gamma

    @staticmethod
    def Cp1_pT(p: float, T: float) -> float:
        """function Cp1_pT = Cp1_pT(p, T)

        5 Equations for Region 1, Section. 5.1 Basic Equation

        Equation 7, Table 3, Page 6

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific isobaric heat capacity in [kJ / (kg K)]
        """
        I1 = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
            8,
            8,
            21,
            23,
            29,
            30,
            31,
            32,
        ]
        J1 = [
            -2,
            -1,
            0,
            1,
            2,
            3,
            4,
            5,
            -9,
            -7,
            -1,
            0,
            1,
            3,
            -3,
            0,
            1,
            3,
            17,
            -4,
            0,
            6,
            -5,
            -2,
            10,
            -8,
            -11,
            -6,
            -29,
            -31,
            -38,
            -39,
            -40,
            -41,
        ]
        n1 = [
            0.14632971213167,
            -0.84548187169114,
            -3.756360367204,
            3.3855169168385,
            -0.95791963387872,
            0.15772038513228,
            -0.016616417199501,
            8.1214629983568e-04,
            2.8319080123804e-04,
            -6.0706301565874e-04,
            -0.018990068218419,
            -0.032529748770505,
            -0.021841717175414,
            -5.283835796993e-05,
            -4.7184321073267e-04,
            -3.0001780793026e-04,
            4.7661393906987e-05,
            -4.4141845330846e-06,
            -7.2694996297594e-16,
            -3.1679644845054e-05,
            -2.8270797985312e-06,
            -8.5205128120103e-10,
            -2.2425281908e-06,
            -6.5171222895601e-07,
            -1.4341729937924e-13,
            -4.0516996860117e-07,
            -1.2734301741641e-09,
            -1.7424871230634e-10,
            -6.8762131295531e-19,
            1.4478307828521e-20,
            2.6335781662795e-23,
            -1.1947622640071e-23,
            1.8228094581404e-24,
            -9.3537087292458e-26,
        ]
        R = SPECIFIC_GAS_CONSTANT
        Pi = p / 16.53
        tau = 1386 / T
        gamma_der_tautau = 0
        for i in range(0, 34):
            gamma_der_tautau = gamma_der_tautau + (
                n1[i]
                * (7.1 - Pi) ** I1[i]
                * J1[i]
                * (J1[i] - 1)
                * (tau - 1.222) ** (J1[i] - 2)
            )
        return -R * tau**2 * gamma_der_tautau

    @staticmethod
    def Cv1_pT(p: float, T: float) -> float:
        """function Cv1_pT = Cv1_pT(p, T)

        5 Equations for Region 1, Section. 5.1 Basic Equation

        Equation 7, Table 3, Page 6

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific isochoric heat capacity in [kJ / (kg K)]
        """
        I1 = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
            8,
            8,
            21,
            23,
            29,
            30,
            31,
            32,
        ]
        J1 = [
            -2,
            -1,
            0,
            1,
            2,
            3,
            4,
            5,
            -9,
            -7,
            -1,
            0,
            1,
            3,
            -3,
            0,
            1,
            3,
            17,
            -4,
            0,
            6,
            -5,
            -2,
            10,
            -8,
            -11,
            -6,
            -29,
            -31,
            -38,
            -39,
            -40,
            -41,
        ]
        n1 = [
            0.14632971213167,
            -0.84548187169114,
            -3.756360367204,
            3.3855169168385,
            -0.95791963387872,
            0.15772038513228,
            -0.016616417199501,
            8.1214629983568e-04,
            2.8319080123804e-04,
            -6.0706301565874e-04,
            -0.018990068218419,
            -0.032529748770505,
            -0.021841717175414,
            -5.283835796993e-05,
            -4.7184321073267e-04,
            -3.0001780793026e-04,
            4.7661393906987e-05,
            -4.4141845330846e-06,
            -7.2694996297594e-16,
            -3.1679644845054e-05,
            -2.8270797985312e-06,
            -8.5205128120103e-10,
            -2.2425281908e-06,
            -6.5171222895601e-07,
            -1.4341729937924e-13,
            -4.0516996860117e-07,
            -1.2734301741641e-09,
            -1.7424871230634e-10,
            -6.8762131295531e-19,
            1.4478307828521e-20,
            2.6335781662795e-23,
            -1.1947622640071e-23,
            1.8228094581404e-24,
            -9.3537087292458e-26,
        ]
        R = SPECIFIC_GAS_CONSTANT
        Pi = p / 16.53
        tau = 1386 / T
        gamma_der_pi = 0
        gamma_der_pipi = 0
        gamma_der_pitau = 0
        gamma_der_tautau = 0
        for i in range(0, 34):
            gamma_der_pi = (
                gamma_der_pi
                - n1[i] * I1[i] * (7.1 - Pi) ** (I1[i] - 1) * (tau - 1.222) ** J1[i]
            )
            gamma_der_pipi = (
                gamma_der_pipi
                + n1[i]
                * I1[i]
                * (I1[i] - 1)
                * (7.1 - Pi) ** (I1[i] - 2)
                * (tau - 1.222) ** J1[i]
            )
            gamma_der_pitau = gamma_der_pitau - n1[i] * I1[i] * (7.1 - Pi) ** (
                I1[i] - 1
            ) * J1[i] * (tau - 1.222) ** (J1[i] - 1)
            gamma_der_tautau = gamma_der_tautau + n1[i] * (7.1 - Pi) ** I1[i] * J1[
                i
            ] * (J1[i] - 1) * (tau - 1.222) ** (J1[i] - 2)
        return R * (
            -(tau**2) * gamma_der_tautau
            + (gamma_der_pi - tau * gamma_der_pitau) ** 2 / gamma_der_pipi
        )

    @staticmethod
    def w1_pT(p: float, T: float) -> float:
        """function w1_pT = w1_pT(p, T)

        5 Equations for Region 1, Section. 5.1 Basic Equation

        Equation 7, Table 3, Page 6

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: speed of sound in [m / s]
        """
        I1 = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
            8,
            8,
            21,
            23,
            29,
            30,
            31,
            32,
        ]
        J1 = [
            -2,
            -1,
            0,
            1,
            2,
            3,
            4,
            5,
            -9,
            -7,
            -1,
            0,
            1,
            3,
            -3,
            0,
            1,
            3,
            17,
            -4,
            0,
            6,
            -5,
            -2,
            10,
            -8,
            -11,
            -6,
            -29,
            -31,
            -38,
            -39,
            -40,
            -41,
        ]
        n1 = [
            0.14632971213167,
            -0.84548187169114,
            -3.756360367204,
            3.3855169168385,
            -0.95791963387872,
            0.15772038513228,
            -0.016616417199501,
            8.1214629983568e-04,
            2.8319080123804e-04,
            -6.0706301565874e-04,
            -0.018990068218419,
            -0.032529748770505,
            -0.021841717175414,
            -5.283835796993e-05,
            -4.7184321073267e-04,
            -3.0001780793026e-04,
            4.7661393906987e-05,
            -4.4141845330846e-06,
            -7.2694996297594e-16,
            -3.1679644845054e-05,
            -2.8270797985312e-06,
            -8.5205128120103e-10,
            -2.2425281908e-06,
            -6.5171222895601e-07,
            -1.4341729937924e-13,
            -4.0516996860117e-07,
            -1.2734301741641e-09,
            -1.7424871230634e-10,
            -6.8762131295531e-19,
            1.4478307828521e-20,
            2.6335781662795e-23,
            -1.1947622640071e-23,
            1.8228094581404e-24,
            -9.3537087292458e-26,
        ]
        R = SPECIFIC_GAS_CONSTANT
        Pi = p / 16.53
        tau = 1386 / T
        gamma_der_pi = 0
        gamma_der_pipi = 0
        gamma_der_pitau = 0
        gamma_der_tautau = 0
        for i in range(0, 34):
            gamma_der_pi = (
                gamma_der_pi
                - n1[i] * I1[i] * (7.1 - Pi) ** (I1[i] - 1) * (tau - 1.222) ** J1[i]
            )
            gamma_der_pipi = (
                gamma_der_pipi
                + n1[i]
                * I1[i]
                * (I1[i] - 1)
                * (7.1 - Pi) ** (I1[i] - 2)
                * (tau - 1.222) ** J1[i]
            )
            gamma_der_pitau = gamma_der_pitau - n1[i] * I1[i] * (7.1 - Pi) ** (
                I1[i] - 1
            ) * J1[i] * (tau - 1.222) ** (J1[i] - 1)
            gamma_der_tautau = gamma_der_tautau + n1[i] * (7.1 - Pi) ** I1[i] * J1[
                i
            ] * (J1[i] - 1) * (tau - 1.222) ** (J1[i] - 2)
        return (
            1000
            * R
            * T
            * gamma_der_pi**2
            / (
                (gamma_der_pi - tau * gamma_der_pitau) ** 2
                / (tau**2 * gamma_der_tautau)
                - gamma_der_pipi
            )
        ) ** 0.5

    @staticmethod
    def T1_ph(p: float, h: float) -> float:
        """function T1_ph = T1_ph(p, h)

        5 Equations for Region 1, Section. 5.1 Basic Equation, 5.2.1 The Backward Equation T (p, h)

        Equation 11, Table 6, Page 10

        :param p: preasure in [MPa]
        :param h: enthalpy in [kJ / kg]

        :return: temperature in [K]
        """
        I1 = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6]
        J1 = [0, 1, 2, 6, 22, 32, 0, 1, 2, 3, 4, 10, 32, 10, 32, 10, 32, 32, 32, 32]
        n1 = [
            -238.72489924521,
            404.21188637945,
            113.49746881718,
            -5.8457616048039,
            -1.528548241314e-04,
            -1.0866707695377e-06,
            -13.391744872602,
            43.211039183559,
            -54.010067170506,
            30.535892203916,
            -6.5964749423638,
            9.3965400878363e-03,
            1.157364750534e-07,
            -2.5858641282073e-05,
            -4.0644363084799e-09,
            6.6456186191635e-08,
            8.0670734103027e-11,
            -9.3477771213947e-13,
            5.8265442020601e-15,
            -1.5020185953503e-17,
        ]
        Pi = p / 1
        eta = h / 2500
        T = 0
        for i in range(0, 20):
            T = T + n1[i] * Pi ** I1[i] * (eta + 1) ** J1[i]
        return T

    @staticmethod
    def T1_ps(p: float, s: float) -> float:
        """function T1_ps = T1_ps(p, s)

        5 Equations for Region 1, Section. 5.1 Basic Equation, 5.2.2 The Backward Equation T (p, s)

        Equation 13, Table 8, Page 11

        :param p: preasure in [MPa]
        :param s: Specific entropy in [kJ / (kg K)]

        :return: temperature in [K]
        """
        I1 = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4]
        J1 = [0, 1, 2, 3, 11, 31, 0, 1, 2, 3, 12, 31, 0, 1, 2, 9, 31, 10, 32, 32]
        n1 = [
            174.78268058307,
            34.806930892873,
            6.5292584978455,
            0.33039981775489,
            -1.9281382923196e-07,
            -2.4909197244573e-23,
            -0.26107636489332,
            0.22592965981586,
            -0.064256463395226,
            7.8876289270526e-03,
            3.5672110607366e-10,
            1.7332496994895e-24,
            5.6608900654837e-04,
            -3.2635483139717e-04,
            4.4778286690632e-05,
            -5.1322156908507e-10,
            -4.2522657042207e-26,
            2.6400441360689e-13,
            7.8124600459723e-29,
            -3.0732199903668e-31,
        ]
        Pi = p / 1
        Sigma = s / 1
        T = 0
        for i in range(0, 20):
            T = T + n1[i] * Pi ** I1[i] * (Sigma + 2) ** J1[i]
        return T

    @staticmethod
    def p1_hs(h: float, s: float) -> float:
        """function p1_hs = p1_hs(h, s)

        Supplementary Release on Backward Equations for Pressure as a Function of Enthalpy and Entropy p(h, s) to the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam

        5 Backward Equation p(h, s) for Region 1

        Equation 1, Table 2, Page 5

        :param h: enthalpy in [kJ / kg]
        :param s: Specific entropy in [kJ / (kg K)]

        :return: preasure in [MPa]
        """
        I1 = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 5]
        J1 = [0, 1, 2, 4, 5, 6, 8, 14, 0, 1, 4, 6, 0, 1, 10, 4, 1, 4, 0]
        n1 = [
            -0.691997014660582,
            -18.361254878756,
            -9.28332409297335,
            65.9639569909906,
            -16.2060388912024,
            450.620017338667,
            854.68067822417,
            6075.23214001162,
            32.6487682621856,
            -26.9408844582931,
            -319.9478483343,
            -928.35430704332,
            30.3634537455249,
            -65.0540422444146,
            -4309.9131651613,
            -747.512324096068,
            730.000345529245,
            1142.84032569021,
            -436.407041874559,
        ]
        eta = h / 3400
        Sigma = s / 7.6
        p = 0
        for i in range(0, 19):
            p = p + n1[i] * (eta + 0.05) ** I1[i] * (Sigma + 0.05) ** J1[i]
        return p * 100

    @staticmethod
    def T1_prho(p: float, rho: float) -> float:
        """function T1_prho = T1_prho(p , rho)

        Solve by iteration. Observe that for low temperatures this equation has 2 solutions. Solve with half interval method

        :param p: pressure in [MPa]
        :param rho: density in [kg / m³]

        :return: temperature in [K]
        """
        logger = logging.getLogger("pyXSteam")
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
                logger.warning(
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
    Section 2.2: IAPWS IF 97 Calling functions to calculate the properties of water in Region 3
    """

    @staticmethod
    def v2_pT(p: float, T: float) -> float:
        """function v2_pT = v2_pT(p, T)

        6 Equations for Region 2, Section. 6.1 Basic Equation

        Table 11 and 12, Page 14 and 15

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific volume in [m³ / kg]
        """
        Ir = [
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
            6,
            6,
            6,
            7,
            7,
            7,
            8,
            8,
            9,
            10,
            10,
            10,
            16,
            16,
            18,
            20,
            20,
            20,
            21,
            22,
            23,
            24,
            24,
            24,
        ]
        Jr = [
            0,
            1,
            2,
            3,
            6,
            1,
            2,
            4,
            7,
            36,
            0,
            1,
            3,
            6,
            35,
            1,
            2,
            3,
            7,
            3,
            16,
            35,
            0,
            11,
            25,
            8,
            36,
            13,
            4,
            10,
            14,
            29,
            50,
            57,
            20,
            35,
            48,
            21,
            53,
            39,
            26,
            40,
            58,
        ]
        nr = [
            -1.7731742473213e-03,
            -0.017834862292358,
            -0.045996013696365,
            -0.057581259083432,
            -0.05032527872793,
            -3.3032641670203e-05,
            -1.8948987516315e-04,
            -3.9392777243355e-03,
            -0.043797295650573,
            -2.6674547914087e-05,
            2.0481737692309e-08,
            4.3870667284435e-07,
            -3.227767723857e-05,
            -1.5033924542148e-03,
            -0.040668253562649,
            -7.8847309559367e-10,
            1.2790717852285e-08,
            4.8225372718507e-07,
            2.2922076337661e-06,
            -1.6714766451061e-11,
            -2.1171472321355e-03,
            -23.895741934104,
            -5.905956432427e-18,
            -1.2621808899101e-06,
            -0.038946842435739,
            1.1256211360459e-11,
            -8.2311340897998,
            1.9809712802088e-08,
            1.0406965210174e-19,
            -1.0234747095929e-13,
            -1.0018179379511e-09,
            -8.0882908646985e-11,
            0.10693031879409,
            -0.33662250574171,
            8.9185845355421e-25,
            3.0629316876232e-13,
            -4.2002467698208e-06,
            -5.9056029685639e-26,
            3.7826947613457e-06,
            -1.2768608934681e-15,
            7.3087610595061e-29,
            5.5414715350778e-17,
            -9.436970724121e-07,
        ]
        R = SPECIFIC_GAS_CONSTANT
        Pi = p
        tau = 540 / T
        g0_pi = 1 / Pi
        gr_pi = 0
        for i in range(0, 43):
            gr_pi = gr_pi + nr[i] * Ir[i] * Pi ** (Ir[i] - 1) * (tau - 0.5) ** Jr[i]
        return R * T / p * Pi * (g0_pi + gr_pi) / 1000

    @staticmethod
    def h2_pT(p: float, T: float) -> float:
        """function h2_pT = h2_pT(p, T)

        6 Equations for Region 2, Section. 6.1 Basic Equation

        Table 11 and 12, Page 14 and 15

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: enthalpy in [kJ / kg]
        """
        J0 = [0, 1, -5, -4, -3, -2, -1, 2, 3]
        n0 = [
            -9.6927686500217,
            10.086655968018,
            -0.005608791128302,
            0.071452738081455,
            -0.40710498223928,
            1.4240819171444,
            -4.383951131945,
            -0.28408632460772,
            0.021268463753307,
        ]
        Ir = [
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
            6,
            6,
            6,
            7,
            7,
            7,
            8,
            8,
            9,
            10,
            10,
            10,
            16,
            16,
            18,
            20,
            20,
            20,
            21,
            22,
            23,
            24,
            24,
            24,
        ]
        Jr = [
            0,
            1,
            2,
            3,
            6,
            1,
            2,
            4,
            7,
            36,
            0,
            1,
            3,
            6,
            35,
            1,
            2,
            3,
            7,
            3,
            16,
            35,
            0,
            11,
            25,
            8,
            36,
            13,
            4,
            10,
            14,
            29,
            50,
            57,
            20,
            35,
            48,
            21,
            53,
            39,
            26,
            40,
            58,
        ]
        nr = [
            -1.7731742473213e-03,
            -0.017834862292358,
            -0.045996013696365,
            -0.057581259083432,
            -0.05032527872793,
            -3.3032641670203e-05,
            -1.8948987516315e-04,
            -3.9392777243355e-03,
            -0.043797295650573,
            -2.6674547914087e-05,
            2.0481737692309e-08,
            4.3870667284435e-07,
            -3.227767723857e-05,
            -1.5033924542148e-03,
            -0.040668253562649,
            -7.8847309559367e-10,
            1.2790717852285e-08,
            4.8225372718507e-07,
            2.2922076337661e-06,
            -1.6714766451061e-11,
            -2.1171472321355e-03,
            -23.895741934104,
            -5.905956432427e-18,
            -1.2621808899101e-06,
            -0.038946842435739,
            1.1256211360459e-11,
            -8.2311340897998,
            1.9809712802088e-08,
            1.0406965210174e-19,
            -1.0234747095929e-13,
            -1.0018179379511e-09,
            -8.0882908646985e-11,
            0.10693031879409,
            -0.33662250574171,
            8.9185845355421e-25,
            3.0629316876232e-13,
            -4.2002467698208e-06,
            -5.9056029685639e-26,
            3.7826947613457e-06,
            -1.2768608934681e-15,
            7.3087610595061e-29,
            5.5414715350778e-17,
            -9.436970724121e-07,
        ]
        R = SPECIFIC_GAS_CONSTANT
        Pi = p
        tau = 540 / T
        g0_tau = 0
        for i in range(0, 9):
            g0_tau = g0_tau + n0[i] * J0[i] * tau ** (J0[i] - 1)
        gr_tau = 0
        for i in range(0, 43):
            gr_tau = gr_tau + nr[i] * Pi ** Ir[i] * Jr[i] * (tau - 0.5) ** (Jr[i] - 1)
        return R * T * tau * (g0_tau + gr_tau)

    @staticmethod
    def u2_pT(p: float, T: float) -> float:
        """function u2_pT = u2_pT(p, T)

        6 Equations for Region 2, Section. 6.1 Basic Equation

        Table 11 and 12, Page 14 and 15

        :param p: preasure in [MPa]
        :param T: temperature in [K]
        :return: specific internal energy in [kJ / kg]
        """
        J0 = [0, 1, -5, -4, -3, -2, -1, 2, 3]
        n0 = [
            -9.6927686500217,
            10.086655968018,
            -0.005608791128302,
            0.071452738081455,
            -0.40710498223928,
            1.4240819171444,
            -4.383951131945,
            -0.28408632460772,
            0.021268463753307,
        ]
        Ir = [
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
            6,
            6,
            6,
            7,
            7,
            7,
            8,
            8,
            9,
            10,
            10,
            10,
            16,
            16,
            18,
            20,
            20,
            20,
            21,
            22,
            23,
            24,
            24,
            24,
        ]
        Jr = [
            0,
            1,
            2,
            3,
            6,
            1,
            2,
            4,
            7,
            36,
            0,
            1,
            3,
            6,
            35,
            1,
            2,
            3,
            7,
            3,
            16,
            35,
            0,
            11,
            25,
            8,
            36,
            13,
            4,
            10,
            14,
            29,
            50,
            57,
            20,
            35,
            48,
            21,
            53,
            39,
            26,
            40,
            58,
        ]
        nr = [
            -1.7731742473213e-03,
            -0.017834862292358,
            -0.045996013696365,
            -0.057581259083432,
            -0.05032527872793,
            -3.3032641670203e-05,
            -1.8948987516315e-04,
            -3.9392777243355e-03,
            -0.043797295650573,
            -2.6674547914087e-05,
            2.0481737692309e-08,
            4.3870667284435e-07,
            -3.227767723857e-05,
            -1.5033924542148e-03,
            -0.040668253562649,
            -7.8847309559367e-10,
            1.2790717852285e-08,
            4.8225372718507e-07,
            2.2922076337661e-06,
            -1.6714766451061e-11,
            -2.1171472321355e-03,
            -23.895741934104,
            -5.905956432427e-18,
            -1.2621808899101e-06,
            -0.038946842435739,
            1.1256211360459e-11,
            -8.2311340897998,
            1.9809712802088e-08,
            1.0406965210174e-19,
            -1.0234747095929e-13,
            -1.0018179379511e-09,
            -8.0882908646985e-11,
            0.10693031879409,
            -0.33662250574171,
            8.9185845355421e-25,
            3.0629316876232e-13,
            -4.2002467698208e-06,
            -5.9056029685639e-26,
            3.7826947613457e-06,
            -1.2768608934681e-15,
            7.3087610595061e-29,
            5.5414715350778e-17,
            -9.436970724121e-07,
        ]
        R = SPECIFIC_GAS_CONSTANT
        Pi = p
        tau = 540 / T
        g0_pi = 1 / Pi
        g0_tau = 0
        for i in range(0, 9):
            g0_tau = g0_tau + n0[i] * J0[i] * tau ** (J0[i] - 1)
        gr_pi = 0
        gr_tau = 0
        for i in range(0, 43):
            gr_pi = gr_pi + nr[i] * Ir[i] * Pi ** (Ir[i] - 1) * (tau - 0.5) ** Jr[i]
            gr_tau = gr_tau + nr[i] * Pi ** Ir[i] * Jr[i] * (tau - 0.5) ** (Jr[i] - 1)
        return R * T * (tau * (g0_tau + gr_tau) - Pi * (g0_pi + gr_pi))

    @staticmethod
    def s2_pT(p: float, T: float) -> float:
        """function s2_pT = s2_pT(p, T)

        6 Equations for Region 2, Section. 6.1 Basic Equation

        Table 11 and 12, Page 14 and 15

        :param p: preasure in [MPa]
        :param T: temperature in [K]
        :return: specific entropy in [kJ / (kg K)]
        """
        J0 = [0, 1, -5, -4, -3, -2, -1, 2, 3]
        n0 = [
            -9.6927686500217,
            10.086655968018,
            -0.005608791128302,
            0.071452738081455,
            -0.40710498223928,
            1.4240819171444,
            -4.383951131945,
            -0.28408632460772,
            0.021268463753307,
        ]
        Ir = [
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
            6,
            6,
            6,
            7,
            7,
            7,
            8,
            8,
            9,
            10,
            10,
            10,
            16,
            16,
            18,
            20,
            20,
            20,
            21,
            22,
            23,
            24,
            24,
            24,
        ]
        Jr = [
            0,
            1,
            2,
            3,
            6,
            1,
            2,
            4,
            7,
            36,
            0,
            1,
            3,
            6,
            35,
            1,
            2,
            3,
            7,
            3,
            16,
            35,
            0,
            11,
            25,
            8,
            36,
            13,
            4,
            10,
            14,
            29,
            50,
            57,
            20,
            35,
            48,
            21,
            53,
            39,
            26,
            40,
            58,
        ]
        nr = [
            -1.7731742473213e-03,
            -0.017834862292358,
            -0.045996013696365,
            -0.057581259083432,
            -0.05032527872793,
            -3.3032641670203e-05,
            -1.8948987516315e-04,
            -3.9392777243355e-03,
            -0.043797295650573,
            -2.6674547914087e-05,
            2.0481737692309e-08,
            4.3870667284435e-07,
            -3.227767723857e-05,
            -1.5033924542148e-03,
            -0.040668253562649,
            -7.8847309559367e-10,
            1.2790717852285e-08,
            4.8225372718507e-07,
            2.2922076337661e-06,
            -1.6714766451061e-11,
            -2.1171472321355e-03,
            -23.895741934104,
            -5.905956432427e-18,
            -1.2621808899101e-06,
            -0.038946842435739,
            1.1256211360459e-11,
            -8.2311340897998,
            1.9809712802088e-08,
            1.0406965210174e-19,
            -1.0234747095929e-13,
            -1.0018179379511e-09,
            -8.0882908646985e-11,
            0.10693031879409,
            -0.33662250574171,
            8.9185845355421e-25,
            3.0629316876232e-13,
            -4.2002467698208e-06,
            -5.9056029685639e-26,
            3.7826947613457e-06,
            -1.2768608934681e-15,
            7.3087610595061e-29,
            5.5414715350778e-17,
            -9.436970724121e-07,
        ]
        R = SPECIFIC_GAS_CONSTANT
        Pi = p
        tau = 540 / T
        g0 = math.log(Pi)
        g0_tau = 0
        for i in range(0, 9):
            g0 = g0 + n0[i] * tau ** J0[i]
            g0_tau = g0_tau + n0[i] * J0[i] * tau ** (J0[i] - 1)
        gr = 0
        gr_tau = 0
        for i in range(0, 43):
            gr = gr + nr[i] * Pi ** Ir[i] * (tau - 0.5) ** Jr[i]
            gr_tau = gr_tau + nr[i] * Pi ** Ir[i] * Jr[i] * (tau - 0.5) ** (Jr[i] - 1)
        return R * (tau * (g0_tau + gr_tau) - (g0 + gr))

    @staticmethod
    def Cp2_pT(p: float, T: float) -> float:
        """function Cp2_pT = Cp2_pT(p, T)

        6 Equations for Region 2, Section. 6.1 Basic Equation

        Table 11 and 12, Page 14 and 15

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific isobaric heat capacity in [kJ / (kg K)]
        """
        J0 = [0, 1, -5, -4, -3, -2, -1, 2, 3]
        n0 = [
            -9.6927686500217,
            10.086655968018,
            -0.005608791128302,
            0.071452738081455,
            -0.40710498223928,
            1.4240819171444,
            -4.383951131945,
            -0.28408632460772,
            0.021268463753307,
        ]
        Ir = [
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
            6,
            6,
            6,
            7,
            7,
            7,
            8,
            8,
            9,
            10,
            10,
            10,
            16,
            16,
            18,
            20,
            20,
            20,
            21,
            22,
            23,
            24,
            24,
            24,
        ]
        Jr = [
            0,
            1,
            2,
            3,
            6,
            1,
            2,
            4,
            7,
            36,
            0,
            1,
            3,
            6,
            35,
            1,
            2,
            3,
            7,
            3,
            16,
            35,
            0,
            11,
            25,
            8,
            36,
            13,
            4,
            10,
            14,
            29,
            50,
            57,
            20,
            35,
            48,
            21,
            53,
            39,
            26,
            40,
            58,
        ]
        nr = [
            -1.7731742473213e-03,
            -0.017834862292358,
            -0.045996013696365,
            -0.057581259083432,
            -0.05032527872793,
            -3.3032641670203e-05,
            -1.8948987516315e-04,
            -3.9392777243355e-03,
            -0.043797295650573,
            -2.6674547914087e-05,
            2.0481737692309e-08,
            4.3870667284435e-07,
            -3.227767723857e-05,
            -1.5033924542148e-03,
            -0.040668253562649,
            -7.8847309559367e-10,
            1.2790717852285e-08,
            4.8225372718507e-07,
            2.2922076337661e-06,
            -1.6714766451061e-11,
            -2.1171472321355e-03,
            -23.895741934104,
            -5.905956432427e-18,
            -1.2621808899101e-06,
            -0.038946842435739,
            1.1256211360459e-11,
            -8.2311340897998,
            1.9809712802088e-08,
            1.0406965210174e-19,
            -1.0234747095929e-13,
            -1.0018179379511e-09,
            -8.0882908646985e-11,
            0.10693031879409,
            -0.33662250574171,
            8.9185845355421e-25,
            3.0629316876232e-13,
            -4.2002467698208e-06,
            -5.9056029685639e-26,
            3.7826947613457e-06,
            -1.2768608934681e-15,
            7.3087610595061e-29,
            5.5414715350778e-17,
            -9.436970724121e-07,
        ]
        R = SPECIFIC_GAS_CONSTANT
        Pi = p
        tau = 540 / T
        g0_tautau = 0
        for i in range(0, 9):
            g0_tautau = g0_tautau + n0[i] * J0[i] * (J0[i] - 1) * tau ** (J0[i] - 2)
        gr_tautau = 0
        for i in range(0, 43):
            gr_tautau = gr_tautau + nr[i] * Pi ** Ir[i] * Jr[i] * (Jr[i] - 1) * (
                tau - 0.5
            ) ** (Jr[i] - 2)
        return -R * tau**2 * (g0_tautau + gr_tautau)

    @staticmethod
    def Cv2_pT(p: float, T: float) -> float:
        """function Cv2_pT = Cv2_pT(p, T)

        6 Equations for Region 2, Section. 6.1 Basic Equation

        Table 11 and 12, Page 14 and 15

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific isochoric heat capacity in [kJ / (kg K)]
        """
        J0 = [0, 1, -5, -4, -3, -2, -1, 2, 3]
        n0 = [
            -9.6927686500217,
            10.086655968018,
            -0.005608791128302,
            0.071452738081455,
            -0.40710498223928,
            1.4240819171444,
            -4.383951131945,
            -0.28408632460772,
            0.021268463753307,
        ]
        Ir = [
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
            6,
            6,
            6,
            7,
            7,
            7,
            8,
            8,
            9,
            10,
            10,
            10,
            16,
            16,
            18,
            20,
            20,
            20,
            21,
            22,
            23,
            24,
            24,
            24,
        ]
        Jr = [
            0,
            1,
            2,
            3,
            6,
            1,
            2,
            4,
            7,
            36,
            0,
            1,
            3,
            6,
            35,
            1,
            2,
            3,
            7,
            3,
            16,
            35,
            0,
            11,
            25,
            8,
            36,
            13,
            4,
            10,
            14,
            29,
            50,
            57,
            20,
            35,
            48,
            21,
            53,
            39,
            26,
            40,
            58,
        ]
        nr = [
            -1.7731742473213e-03,
            -0.017834862292358,
            -0.045996013696365,
            -0.057581259083432,
            -0.05032527872793,
            -3.3032641670203e-05,
            -1.8948987516315e-04,
            -3.9392777243355e-03,
            -0.043797295650573,
            -2.6674547914087e-05,
            2.0481737692309e-08,
            4.3870667284435e-07,
            -3.227767723857e-05,
            -1.5033924542148e-03,
            -0.040668253562649,
            -7.8847309559367e-10,
            1.2790717852285e-08,
            4.8225372718507e-07,
            2.2922076337661e-06,
            -1.6714766451061e-11,
            -2.1171472321355e-03,
            -23.895741934104,
            -5.905956432427e-18,
            -1.2621808899101e-06,
            -0.038946842435739,
            1.1256211360459e-11,
            -8.2311340897998,
            1.9809712802088e-08,
            1.0406965210174e-19,
            -1.0234747095929e-13,
            -1.0018179379511e-09,
            -8.0882908646985e-11,
            0.10693031879409,
            -0.33662250574171,
            8.9185845355421e-25,
            3.0629316876232e-13,
            -4.2002467698208e-06,
            -5.9056029685639e-26,
            3.7826947613457e-06,
            -1.2768608934681e-15,
            7.3087610595061e-29,
            5.5414715350778e-17,
            -9.436970724121e-07,
        ]
        R = SPECIFIC_GAS_CONSTANT
        Pi = p
        tau = 540 / T
        g0_tautau = 0
        for i in range(0, 9):
            g0_tautau = g0_tautau + n0[i] * J0[i] * (J0[i] - 1) * tau ** (J0[i] - 2)
        gr_pi = 0
        gr_pitau = 0
        gr_pipi = 0
        gr_tautau = 0
        for i in range(0, 43):
            gr_pi = gr_pi + nr[i] * Ir[i] * Pi ** (Ir[i] - 1) * (tau - 0.5) ** Jr[i]
            gr_pipi = (
                gr_pipi
                + nr[i] * Ir[i] * (Ir[i] - 1) * Pi ** (Ir[i] - 2) * (tau - 0.5) ** Jr[i]
            )
            gr_pitau = gr_pitau + nr[i] * Ir[i] * Pi ** (Ir[i] - 1) * Jr[i] * (
                tau - 0.5
            ) ** (Jr[i] - 1)
            gr_tautau = gr_tautau + nr[i] * Pi ** Ir[i] * Jr[i] * (Jr[i] - 1) * (
                tau - 0.5
            ) ** (Jr[i] - 2)
        return R * (
            -(tau**2) * (g0_tautau + gr_tautau)
            - (1 + Pi * gr_pi - tau * Pi * gr_pitau) ** 2 / (1 - Pi**2 * gr_pipi)
        )

    @staticmethod
    def w2_pT(p: float, T: float) -> float:
        """function w2_pT = w2_pT(p, T)

        6 Equations for Region 2, Section. 6.1 Basic Equation

        Table 11 and 12, Page 14 and 15

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: speed of sound in [m / s]
        """
        J0 = [0, 1, -5, -4, -3, -2, -1, 2, 3]
        n0 = [
            -9.6927686500217,
            10.086655968018,
            -0.005608791128302,
            0.071452738081455,
            -0.40710498223928,
            1.4240819171444,
            -4.383951131945,
            -0.28408632460772,
            0.021268463753307,
        ]
        Ir = [
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
            6,
            6,
            6,
            7,
            7,
            7,
            8,
            8,
            9,
            10,
            10,
            10,
            16,
            16,
            18,
            20,
            20,
            20,
            21,
            22,
            23,
            24,
            24,
            24,
        ]
        Jr = [
            0,
            1,
            2,
            3,
            6,
            1,
            2,
            4,
            7,
            36,
            0,
            1,
            3,
            6,
            35,
            1,
            2,
            3,
            7,
            3,
            16,
            35,
            0,
            11,
            25,
            8,
            36,
            13,
            4,
            10,
            14,
            29,
            50,
            57,
            20,
            35,
            48,
            21,
            53,
            39,
            26,
            40,
            58,
        ]
        nr = [
            -1.7731742473213e-03,
            -0.017834862292358,
            -0.045996013696365,
            -0.057581259083432,
            -0.05032527872793,
            -3.3032641670203e-05,
            -1.8948987516315e-04,
            -3.9392777243355e-03,
            -0.043797295650573,
            -2.6674547914087e-05,
            2.0481737692309e-08,
            4.3870667284435e-07,
            -3.227767723857e-05,
            -1.5033924542148e-03,
            -0.040668253562649,
            -7.8847309559367e-10,
            1.2790717852285e-08,
            4.8225372718507e-07,
            2.2922076337661e-06,
            -1.6714766451061e-11,
            -2.1171472321355e-03,
            -23.895741934104,
            -5.905956432427e-18,
            -1.2621808899101e-06,
            -0.038946842435739,
            1.1256211360459e-11,
            -8.2311340897998,
            1.9809712802088e-08,
            1.0406965210174e-19,
            -1.0234747095929e-13,
            -1.0018179379511e-09,
            -8.0882908646985e-11,
            0.10693031879409,
            -0.33662250574171,
            8.9185845355421e-25,
            3.0629316876232e-13,
            -4.2002467698208e-06,
            -5.9056029685639e-26,
            3.7826947613457e-06,
            -1.2768608934681e-15,
            7.3087610595061e-29,
            5.5414715350778e-17,
            -9.436970724121e-07,
        ]
        R = SPECIFIC_GAS_CONSTANT
        Pi = p
        tau = 540 / T
        g0_tautau = 0
        for i in range(0, 9):
            g0_tautau = g0_tautau + n0[i] * J0[i] * (J0[i] - 1) * tau ** (J0[i] - 2)
        gr_pi = 0
        gr_pitau = 0
        gr_pipi = 0
        gr_tautau = 0
        for i in range(0, 43):
            gr_pi = gr_pi + nr[i] * Ir[i] * Pi ** (Ir[i] - 1) * (tau - 0.5) ** Jr[i]
            gr_pipi = (
                gr_pipi
                + nr[i] * Ir[i] * (Ir[i] - 1) * Pi ** (Ir[i] - 2) * (tau - 0.5) ** Jr[i]
            )
            gr_pitau = gr_pitau + nr[i] * Ir[i] * Pi ** (Ir[i] - 1) * Jr[i] * (
                tau - 0.5
            ) ** (Jr[i] - 1)
            gr_tautau = gr_tautau + nr[i] * Pi ** Ir[i] * Jr[i] * (Jr[i] - 1) * (
                tau - 0.5
            ) ** (Jr[i] - 2)
        return (
            1000
            * R
            * T
            * (1 + 2 * Pi * gr_pi + Pi**2 * gr_pi**2)
            / (
                (1 - Pi**2 * gr_pipi)
                + (1 + Pi * gr_pi - tau * Pi * gr_pitau) ** 2
                / (tau**2 * (g0_tautau + gr_tautau))
            )
        ) ** 0.5

    @staticmethod
    def T2_ph(p: float, h: float) -> float:
        """function T2_ph = T2_ph(p, h)

        6 Equations for Region 2, 6.3.1 The Backward Equations T(p, h) for Subregions 2a, 2b, and 2c

        :param p: preasure in [MPa]
        :param h: enthalpy in [kJ / kg]

        :return: temperature in [K]
        """
        if p < 4:
            sub_reg = 1
        else:
            if p < (
                905.84278514723 - 0.67955786399241 * h + 1.2809002730136e-04 * h**2
            ):
                sub_reg = 2
            else:
                sub_reg = 3
        if sub_reg == 1:
            Ji = [
                0,
                1,
                2,
                3,
                7,
                20,
                0,
                1,
                2,
                3,
                7,
                9,
                11,
                18,
                44,
                0,
                2,
                7,
                36,
                38,
                40,
                42,
                44,
                24,
                44,
                12,
                32,
                44,
                32,
                36,
                42,
                34,
                44,
                28,
            ]
            Ii = [
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                3,
                3,
                4,
                4,
                4,
                5,
                5,
                5,
                6,
                6,
                7,
            ]
            ni = [
                1089.8952318288,
                849.51654495535,
                -107.81748091826,
                33.153654801263,
                -7.4232016790248,
                11.765048724356,
                1.844574935579,
                -4.1792700549624,
                6.2478196935812,
                -17.344563108114,
                -200.58176862096,
                271.96065473796,
                -455.11318285818,
                3091.9688604755,
                252266.40357872,
                -6.1707422868339e-03,
                -0.31078046629583,
                11.670873077107,
                128127984.04046,
                -985549096.23276,
                2822454697.3002,
                -3594897141.0703,
                1722734991.3197,
                -13551.334240775,
                12848734.66465,
                1.3865724283226,
                235988.32556514,
                -13105236.545054,
                7399.9835474766,
                -551966.9703006,
                3715408.5996233,
                19127.72923966,
                -415351.64835634,
                -62.459855192507,
            ]
            Ts = 0
            hs = h / 2000
            for i in range(0, 34):
                Ts = Ts + ni[i] * p ** (Ii[i]) * (hs - 2.1) ** Ji[i]
        elif sub_reg == 2:
            # Subregion B
            # Table 21, Eq 23, page 23
            Ji = [
                0,
                1,
                2,
                12,
                18,
                24,
                28,
                40,
                0,
                2,
                6,
                12,
                18,
                24,
                28,
                40,
                2,
                8,
                18,
                40,
                1,
                2,
                12,
                24,
                2,
                12,
                18,
                24,
                28,
                40,
                18,
                24,
                40,
                28,
                2,
                28,
                1,
                40,
            ]
            Ii = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                2,
                2,
                2,
                2,
                3,
                3,
                3,
                3,
                4,
                4,
                4,
                4,
                4,
                4,
                5,
                5,
                5,
                6,
                7,
                7,
                9,
                9,
            ]
            ni = [
                1489.5041079516,
                743.07798314034,
                -97.708318797837,
                2.4742464705674,
                -0.63281320016026,
                1.1385952129658,
                -0.47811863648625,
                8.5208123431544e-03,
                0.93747147377932,
                3.3593118604916,
                3.3809355601454,
                0.16844539671904,
                0.73875745236695,
                -0.47128737436186,
                0.15020273139707,
                -0.002176411421975,
                -0.021810755324761,
                -0.10829784403677,
                -0.046333324635812,
                7.1280351959551e-05,
                1.1032831789999e-04,
                1.8955248387902e-04,
                3.0891541160537e-03,
                1.3555504554949e-03,
                2.8640237477456e-07,
                -1.0779857357512e-05,
                -7.6462712454814e-05,
                1.4052392818316e-05,
                -3.1083814331434e-05,
                -1.0302738212103e-06,
                2.821728163504e-07,
                1.2704902271945e-06,
                7.3803353468292e-08,
                -1.1030139238909e-08,
                -8.1456365207833e-14,
                -2.5180545682962e-11,
                -1.7565233969407e-18,
                8.6934156344163e-15,
            ]
            Ts = 0
            hs = h / 2000
            # for i = 1 : 38
            for i in range(0, 38):
                Ts = Ts + ni[i] * (p - 2) ** (Ii[i]) * (hs - 2.6) ** Ji[i]
        else:
            # Subregion C
            # Table 22, Eq 24, page 24
            Ji = [
                0,
                4,
                0,
                2,
                0,
                2,
                0,
                1,
                0,
                2,
                0,
                1,
                4,
                8,
                4,
                0,
                1,
                4,
                10,
                12,
                16,
                20,
                22,
            ]
            Ii = [
                -7,
                -7,
                -6,
                -6,
                -5,
                -5,
                -2,
                -2,
                -1,
                -1,
                0,
                0,
                1,
                1,
                2,
                6,
                6,
                6,
                6,
                6,
                6,
                6,
                6,
            ]
            ni = [
                -3236839855524.2,
                7326335090218.1,
                358250899454.47,
                -583401318515.9,
                -10783068217.47,
                20825544563.171,
                610747.83564516,
                859777.2253558,
                -25745.72360417,
                31081.088422714,
                1208.2315865936,
                482.19755109255,
                3.7966001272486,
                -10.842984880077,
                -0.04536417267666,
                1.4559115658698e-13,
                1.126159740723e-12,
                -1.7804982240686e-11,
                1.2324579690832e-07,
                -1.1606921130984e-06,
                2.7846367088554e-05,
                -5.9270038474176e-04,
                1.2918582991878e-03,
            ]
            Ts = 0
            hs = h / 2000
            for i in range(0, 23):
                Ts = Ts + ni[i] * (p + 25) ** (Ii[i]) * (hs - 1.8) ** Ji[i]
        return Ts

    @staticmethod
    def T2_ps(p: float, s: float) -> float:
        """function T2_ps = T2_ps(p, s)

        6 Equations for Region 2,6.3.2 The Backward Equations T( p, s ) for Subregions 2a, 2b, and 2c

        Page 26

        :param p: preasure in [MPa]
        :param s: Specific entropy in [kJ / (kg K)]

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
            Ii = [
                -1.5,
                -1.5,
                -1.5,
                -1.5,
                -1.5,
                -1.5,
                -1.25,
                -1.25,
                -1.25,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -0.75,
                -0.75,
                -0.5,
                -0.5,
                -0.5,
                -0.5,
                -0.25,
                -0.25,
                -0.25,
                -0.25,
                0.25,
                0.25,
                0.25,
                0.25,
                0.5,
                0.5,
                0.5,
                0.5,
                0.5,
                0.5,
                0.5,
                0.75,
                0.75,
                0.75,
                0.75,
                1,
                1,
                1.25,
                1.25,
                1.5,
                1.5,
            ]
            Ji = [
                -24,
                -23,
                -19,
                -13,
                -11,
                -10,
                -19,
                -15,
                -6,
                -26,
                -21,
                -17,
                -16,
                -9,
                -8,
                -15,
                -14,
                -26,
                -13,
                -9,
                -7,
                -27,
                -25,
                -11,
                -6,
                1,
                4,
                8,
                11,
                0,
                1,
                5,
                6,
                10,
                14,
                16,
                0,
                4,
                9,
                17,
                7,
                18,
                3,
                15,
                5,
                18,
            ]
            ni = [
                -392359.83861984,
                515265.7382727,
                40482.443161048,
                -321.93790923902,
                96.961424218694,
                -22.867846371773,
                -449429.14124357,
                -5011.8336020166,
                0.35684463560015,
                44235.33584819,
                -13673.388811708,
                421632.60207864,
                22516.925837475,
                474.42144865646,
                -149.31130797647,
                -197811.26320452,
                -23554.39947076,
                -19070.616302076,
                55375.669883164,
                3829.3691437363,
                -603.91860580567,
                1936.3102620331,
                4266.064369861,
                -5978.0638872718,
                -704.01463926862,
                338.36784107553,
                20.862786635187,
                0.033834172656196,
                -4.3124428414893e-05,
                166.53791356412,
                -139.86292055898,
                -0.78849547999872,
                0.072132411753872,
                -5.9754839398283e-03,
                -1.2141358953904e-05,
                2.3227096733871e-07,
                -10.538463566194,
                2.0718925496502,
                -0.072193155260427,
                2.074988708112e-07,
                -0.018340657911379,
                2.9036272348696e-07,
                0.21037527893619,
                2.5681239729999e-04,
                -0.012799002933781,
                -8.2198102652018e-06,
            ]
            Pi = p
            Sigma = s / 2
            teta = 0
            # for i = 1 : 46
            for i in range(0, 46):
                teta = teta + ni[i] * Pi ** Ii[i] * (Sigma - 2) ** Ji[i]
        elif sub_reg == 2:
            # Subregion B
            # Table 26, Eq 26, page 27
            Ii = [
                -6,
                -6,
                -5,
                -5,
                -4,
                -4,
                -4,
                -3,
                -3,
                -3,
                -3,
                -2,
                -2,
                -2,
                -2,
                -1,
                -1,
                -1,
                -1,
                -1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                1,
                1,
                1,
                1,
                1,
                2,
                2,
                2,
                3,
                3,
                3,
                4,
                4,
                5,
                5,
                5,
            ]
            Ji = [
                0,
                11,
                0,
                11,
                0,
                1,
                11,
                0,
                1,
                11,
                12,
                0,
                1,
                6,
                10,
                0,
                1,
                5,
                8,
                9,
                0,
                1,
                2,
                4,
                5,
                6,
                9,
                0,
                1,
                2,
                3,
                7,
                8,
                0,
                1,
                5,
                0,
                1,
                3,
                0,
                1,
                0,
                1,
                2,
            ]
            ni = [
                316876.65083497,
                20.864175881858,
                -398593.99803599,
                -21.816058518877,
                223697.85194242,
                -2784.1703445817,
                9.920743607148,
                -75197.512299157,
                2970.8605951158,
                -3.4406878548526,
                0.38815564249115,
                17511.29508575,
                -1423.7112854449,
                1.0943803364167,
                0.89971619308495,
                -3375.9740098958,
                471.62885818355,
                -1.9188241993679,
                0.41078580492196,
                -0.33465378172097,
                1387.0034777505,
                -406.63326195838,
                41.72734715961,
                2.1932549434532,
                -1.0320050009077,
                0.35882943516703,
                5.2511453726066e-03,
                12.838916450705,
                -2.8642437219381,
                0.56912683664855,
                -0.099962954584931,
                -3.2632037778459e-03,
                2.3320922576723e-04,
                -0.1533480985745,
                0.029072288239902,
                3.7534702741167e-04,
                1.7296691702411e-03,
                -3.8556050844504e-04,
                -3.5017712292608e-05,
                -1.4566393631492e-05,
                5.6420857267269e-06,
                4.1286150074605e-08,
                -2.0684671118824e-08,
                1.6409393674725e-09,
            ]
            Pi = p
            Sigma = s / 0.7853
            teta = 0
            for i in range(0, 44):
                teta = teta + ni[i] * Pi ** Ii[i] * (10 - Sigma) ** Ji[i]
        else:
            # Subregion C
            # Table 27, Eq 27, page 28
            Ii = [
                -2,
                -2,
                -1,
                0,
                0,
                0,
                0,
                1,
                1,
                1,
                1,
                2,
                2,
                2,
                3,
                3,
                3,
                4,
                4,
                4,
                5,
                5,
                5,
                6,
                6,
                7,
                7,
                7,
                7,
                7,
            ]
            Ji = [
                0,
                1,
                0,
                0,
                1,
                2,
                3,
                0,
                1,
                3,
                4,
                0,
                1,
                2,
                0,
                1,
                5,
                0,
                1,
                4,
                0,
                1,
                2,
                0,
                1,
                0,
                1,
                3,
                4,
                5,
            ]
            ni = [
                909.68501005365,
                2404.566708842,
                -591.6232638713,
                541.45404128074,
                -270.98308411192,
                979.76525097926,
                -469.66772959435,
                14.399274604723,
                -19.104204230429,
                5.3299167111971,
                -21.252975375934,
                -0.3114733441376,
                0.60334840894623,
                -0.042764839702509,
                5.8185597255259e-03,
                -0.014597008284753,
                5.6631175631027e-03,
                -7.6155864584577e-05,
                2.2440342919332e-04,
                -1.2561095013413e-05,
                6.3323132660934e-07,
                -2.0541989675375e-06,
                3.6405370390082e-08,
                -2.9759897789215e-09,
                1.0136618529763e-08,
                5.9925719692351e-12,
                -2.0677870105164e-11,
                -2.0874278181886e-11,
                1.0162166825089e-10,
                -1.6429828281347e-10,
            ]
            Pi = p
            Sigma = s / 2.9251
            teta = 0
            for i in range(0, 30):
                teta = teta + ni[i] * Pi ** Ii[i] * (2 - Sigma) ** Ji[i]
        return teta

    @staticmethod
    def p2_hs(h: float, s: float) -> float:
        """function p2_hs = p2_hs(h, s)

        Supplementary Release on Backward Equations for Pressure as a function of Enthalpy and Entropy p(h,s) to the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam

        Chapter 6: Backward Equations p(h,s) for Region 2

        :param h: enthalpy in [kJ / kg]
        :param s: Specific entropy in [kJ / (kg K)]

        :return: preasure in [MPa]
        """
        if h < (
            -3498.98083432139
            + 2575.60716905876 * s
            - 421.073558227969 * s**2
            + 27.6349063799944 * s**3
        ):
            sub_reg = 1
        else:
            if s < 5.85:
                sub_reg = 3
            else:
                sub_reg = 2
        if sub_reg == 1:
            # Subregion A
            # Table 6, Eq 3, page 8
            Ii = [
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                2,
                2,
                2,
                3,
                3,
                3,
                3,
                3,
                4,
                5,
                5,
                6,
                7,
            ]
            Ji = [
                1,
                3,
                6,
                16,
                20,
                22,
                0,
                1,
                2,
                3,
                5,
                6,
                10,
                16,
                20,
                22,
                3,
                16,
                20,
                0,
                2,
                3,
                6,
                16,
                16,
                3,
                16,
                3,
                1,
            ]
            ni = [
                -1.82575361923032e-02,
                -0.125229548799536,
                0.592290437320145,
                6.04769706185122,
                238.624965444474,
                -298.639090222922,
                0.051225081304075,
                -0.437266515606486,
                0.413336902999504,
                -5.16468254574773,
                -5.57014838445711,
                12.8555037824478,
                11.414410895329,
                -119.504225652714,
                -2847.7798596156,
                4317.57846408006,
                1.1289404080265,
                1974.09186206319,
                1516.12444706087,
                1.41324451421235e-02,
                0.585501282219601,
                -2.97258075863012,
                5.94567314847319,
                -6236.56565798905,
                9659.86235133332,
                6.81500934948134,
                -6332.07286824489,
                -5.5891922446576,
                4.00645798472063e-02,
            ]
            eta = h / 4200
            Sigma = s / 12
            Pi = 0
            for i in range(0, 29):
                Pi = Pi + ni[i] * (eta - 0.5) ** Ii[i] * (Sigma - 1.2) ** Ji[i]
            p2_hs = Pi**4 * 4
        elif sub_reg == 2:
            # Subregion B
            # Table 7, Eq 4, page 9
            Ii = [
                0,
                0,
                0,
                0,
                0,
                1,
                1,
                1,
                1,
                1,
                1,
                2,
                2,
                2,
                3,
                3,
                3,
                3,
                4,
                4,
                5,
                5,
                6,
                6,
                6,
                7,
                7,
                8,
                8,
                8,
                8,
                12,
                14,
            ]
            Ji = [
                0,
                1,
                2,
                4,
                8,
                0,
                1,
                2,
                3,
                5,
                12,
                1,
                6,
                18,
                0,
                1,
                7,
                12,
                1,
                16,
                1,
                12,
                1,
                8,
                18,
                1,
                16,
                1,
                3,
                14,
                18,
                10,
                16,
            ]
            ni = [
                8.01496989929495e-02,
                -0.543862807146111,
                0.337455597421283,
                8.9055545115745,
                313.840736431485,
                0.797367065977789,
                -1.2161697355624,
                8.72803386937477,
                -16.9769781757602,
                -186.552827328416,
                95115.9274344237,
                -18.9168510120494,
                -4334.0703719484,
                543212633.012715,
                0.144793408386013,
                128.024559637516,
                -67230.9534071268,
                33697238.0095287,
                -586.63419676272,
                -22140322476.9889,
                1716.06668708389,
                -570817595.806302,
                -3121.09693178482,
                -2078413.8463301,
                3056059461577.86,
                3221.57004314333,
                326810259797.295,
                -1441.04158934487,
                410.694867802691,
                109077066873.024,
                -24796465425889.3,
                1888019068.65134,
                -123651009018773,
            ]
            eta = h / 4100
            Sigma = s / 7.9
            Pi = 0
            for i in range(0, 33):
                Pi = Pi + ni[i] * (eta - 0.6) ** Ii[i] * (Sigma - 1.01) ** Ji[i]
            p2_hs = Pi**4 * 100
        else:
            # Subregion C
            # Table 8, Eq 5, page 10
            Ii = [
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                1,
                1,
                1,
                1,
                2,
                2,
                2,
                2,
                2,
                3,
                3,
                3,
                3,
                3,
                4,
                5,
                5,
                5,
                5,
                6,
                6,
                10,
                12,
                16,
            ]
            Ji = [
                0,
                1,
                2,
                3,
                4,
                8,
                0,
                2,
                5,
                8,
                14,
                2,
                3,
                7,
                10,
                18,
                0,
                5,
                8,
                16,
                18,
                18,
                1,
                4,
                6,
                14,
                8,
                18,
                7,
                7,
                10,
            ]
            ni = [
                0.112225607199012,
                -3.39005953606712,
                -32.0503911730094,
                -197.5973051049,
                -407.693861553446,
                13294.3775222331,
                1.70846839774007,
                37.3694198142245,
                3581.44365815434,
                423014.446424664,
                -751071025.760063,
                52.3446127607898,
                -228.351290812417,
                -960652.417056937,
                -80705929.2526074,
                1626980172256.69,
                0.772465073604171,
                46392.9973837746,
                -13731788.5134128,
                1704703926305.12,
                -25110462818730.8,
                31774883083552,
                53.8685623675312,
                -55308.9094625169,
                -1028615.22421405,
                2042494187562.34,
                273918446.626977,
                -2.63963146312685e15,
                -1078908541.08088,
                -29649262098.0124,
                -1.11754907323424e15,
            ]
            eta = h / 3500
            Sigma = s / 5.9
            Pi = 0
            for i in range(0, 31):
                Pi = Pi + ni[i] * (eta - 0.7) ** Ii[i] * (Sigma - 1.1) ** Ji[i]
            p2_hs = Pi**4 * 100
        return p2_hs

    @staticmethod
    def T2_prho(p: float, rho: float) -> float:
        """function T2_prho=T2_prho(p,rho)
        Solve by iteration. Observe that of low temperatures this equation has 2 solutions. Solve with half interval method

        :param p: preasure in [MPa]
        :param rho: density in [kg / m³]

        :return: temperature in [K]
        """
        logger = logging.getLogger("pyXSteam")
        if p < 16.5292:
            Low_Bound = Region4.T4_p(p)
        else:
            Low_Bound = B23T_p(p)
        High_Bound = 1073.15
        rhos = -1000
        step_counter = 0
        while math.fabs(rho - rhos) > 0.000001:
            step_counter += 1
            last_rhos = rhos

            Ts = (Low_Bound + High_Bound) / 2
            rhos = 1 / Region2.v2_pT(p, Ts)

            if last_rhos == rhos:
                logger.warning(
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
    Section 2.3: IAPWS IF 97 Calling functions to calculate the properties of water in Region 3
    """

    @staticmethod
    def p3_rhoT(rho: float, T: float) -> float:
        """function p3_rhoT = p3_rhoT(rho, T)

        7 Basic Equation for Region 3, Section. 6.1 Basic Equation

        Table 30 and 31, Page 30 and 31

        :param rho: density in [kg / m³]
        :param T: temperature in [K]

        :return: preasure in [MPa]
        """
        Ii = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            4,
            5,
            5,
            5,
            6,
            6,
            6,
            7,
            8,
            9,
            9,
            10,
            10,
            11,
        ]
        Ji = [
            0,
            0,
            1,
            2,
            7,
            10,
            12,
            23,
            2,
            6,
            15,
            17,
            0,
            2,
            6,
            7,
            22,
            26,
            0,
            2,
            4,
            16,
            26,
            0,
            2,
            4,
            26,
            1,
            3,
            26,
            0,
            2,
            26,
            2,
            26,
            2,
            26,
            0,
            1,
            26,
        ]
        ni = [
            1.0658070028513,
            -15.732845290239,
            20.944396974307,
            -7.6867707878716,
            2.6185947787954,
            -2.808078114862,
            1.2053369696517,
            -8.4566812812502e-03,
            -1.2654315477714,
            -1.1524407806681,
            0.88521043984318,
            -0.64207765181607,
            0.38493460186671,
            -0.85214708824206,
            4.8972281541877,
            -3.0502617256965,
            0.039420536879154,
            0.12558408424308,
            -0.2799932969871,
            1.389979956946,
            -2.018991502357,
            -8.2147637173963e-03,
            -0.47596035734923,
            0.0439840744735,
            -0.44476435428739,
            0.90572070719733,
            0.70522450087967,
            0.10770512626332,
            -0.32913623258954,
            -0.50871062041158,
            -0.022175400873096,
            0.094260751665092,
            0.16436278447961,
            -0.013503372241348,
            -0.014834345352472,
            5.7922953628084e-04,
            3.2308904703711e-03,
            8.0964802996215e-05,
            -1.6557679795037e-04,
            -4.4923899061815e-05,
        ]
        R = SPECIFIC_GAS_CONSTANT
        tc = CRITICAL_TEMPERATURE
        rhoc = 322.0  # kg/m3
        delta = rho / rhoc
        tau = tc / T
        fidelta = 0
        for i in range(1, 40):
            fidelta = fidelta + ni[i] * Ii[i] * delta ** (Ii[i] - 1) * tau ** Ji[i]

        fidelta = fidelta + (ni[0] / delta)
        return (rho * R * T * delta * fidelta) / 1000.0

    @staticmethod
    def u3_rhoT(rho: float, T: float) -> float:
        """function u3_rhoT = u3_rhoT(rho, T)

        7 Basic Equation for Region 3, Section. 6.1 Basic Equation

        Table 30 and 31, Page 30 and 31

        :param rho: density in [kg / m³]
        :param T: temperature in [K]

        :return: specific internal energy in [kJ / kg]
        """
        Ii = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            4,
            5,
            5,
            5,
            6,
            6,
            6,
            7,
            8,
            9,
            9,
            10,
            10,
            11,
        ]
        Ji = [
            0,
            0,
            1,
            2,
            7,
            10,
            12,
            23,
            2,
            6,
            15,
            17,
            0,
            2,
            6,
            7,
            22,
            26,
            0,
            2,
            4,
            16,
            26,
            0,
            2,
            4,
            26,
            1,
            3,
            26,
            0,
            2,
            26,
            2,
            26,
            2,
            26,
            0,
            1,
            26,
        ]
        ni = [
            1.0658070028513,
            -15.732845290239,
            20.944396974307,
            -7.6867707878716,
            2.6185947787954,
            -2.808078114862,
            1.2053369696517,
            -8.4566812812502e-03,
            -1.2654315477714,
            -1.1524407806681,
            0.88521043984318,
            -0.64207765181607,
            0.38493460186671,
            -0.85214708824206,
            4.8972281541877,
            -3.0502617256965,
            0.039420536879154,
            0.12558408424308,
            -0.2799932969871,
            1.389979956946,
            -2.018991502357,
            -8.2147637173963e-03,
            -0.47596035734923,
            0.0439840744735,
            -0.44476435428739,
            0.90572070719733,
            0.70522450087967,
            0.10770512626332,
            -0.32913623258954,
            -0.50871062041158,
            -0.022175400873096,
            0.094260751665092,
            0.16436278447961,
            -0.013503372241348,
            -0.014834345352472,
            5.7922953628084e-04,
            3.2308904703711e-03,
            8.0964802996215e-05,
            -1.6557679795037e-04,
            -4.4923899061815e-05,
        ]
        R = SPECIFIC_GAS_CONSTANT
        tc = CRITICAL_TEMPERATURE
        rhoc = 322.0  # kg/m3
        delta = rho / rhoc
        tau = tc / T
        fitau = 0
        for i in range(1, 40):
            fitau = fitau + ni[i] * delta ** Ii[i] * Ji[i] * tau ** (Ji[i] - 1)
        return R * T * (tau * fitau)

    @staticmethod
    def h3_rhoT(rho: float, T: float) -> float:
        """function h3_rhoT = h3_rhoT(rho, T)

        7 Basic Equation for Region 3, Section. 6.1 Basic Equation

        Table 30 and 31, Page 30 and 31

        :param rho: density in [kg / m³]
        :param T: temperature in [K]

        :return: enthalpy in [kJ / kg]
        """
        Ii = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            4,
            5,
            5,
            5,
            6,
            6,
            6,
            7,
            8,
            9,
            9,
            10,
            10,
            11,
        ]
        Ji = [
            0,
            0,
            1,
            2,
            7,
            10,
            12,
            23,
            2,
            6,
            15,
            17,
            0,
            2,
            6,
            7,
            22,
            26,
            0,
            2,
            4,
            16,
            26,
            0,
            2,
            4,
            26,
            1,
            3,
            26,
            0,
            2,
            26,
            2,
            26,
            2,
            26,
            0,
            1,
            26,
        ]
        ni = [
            1.0658070028513,
            -15.732845290239,
            20.944396974307,
            -7.6867707878716,
            2.6185947787954,
            -2.808078114862,
            1.2053369696517,
            -8.4566812812502e-03,
            -1.2654315477714,
            -1.1524407806681,
            0.88521043984318,
            -0.64207765181607,
            0.38493460186671,
            -0.85214708824206,
            4.8972281541877,
            -3.0502617256965,
            0.039420536879154,
            0.12558408424308,
            -0.2799932969871,
            1.389979956946,
            -2.018991502357,
            -8.2147637173963e-03,
            -0.47596035734923,
            0.0439840744735,
            -0.44476435428739,
            0.90572070719733,
            0.70522450087967,
            0.10770512626332,
            -0.32913623258954,
            -0.50871062041158,
            -0.022175400873096,
            0.094260751665092,
            0.16436278447961,
            -0.013503372241348,
            -0.014834345352472,
            5.7922953628084e-04,
            3.2308904703711e-03,
            8.0964802996215e-05,
            -1.6557679795037e-04,
            -4.4923899061815e-05,
        ]
        R = SPECIFIC_GAS_CONSTANT
        tc = CRITICAL_TEMPERATURE
        rhoc = 322.0  # kg/m3
        delta = rho / rhoc
        tau = tc / T
        fidelta = 0
        fitau = 0
        for i in range(1, 40):
            fidelta = fidelta + ni[i] * Ii[i] * delta ** (Ii[i] - 1) * tau ** Ji[i]
            fitau = fitau + ni[i] * delta ** Ii[i] * Ji[i] * tau ** (Ji[i] - 1)
        fidelta = fidelta + ni[0] / delta
        return R * T * (tau * fitau + delta * fidelta)

    @staticmethod
    def s3_rhoT(rho: float, T: float) -> float:
        """function s3_rhoT = s3_rhoT(rho, T)

        7 Basic Equation for Region 3, Section. 6.1 Basic Equation

        Table 30 and 31, Page 30 and 31

        :param rho: density in [kg / m³]
        :param T: temperature in [K]

        :return: specific entropy in [kJ / (kg K)]
        """
        Ii = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            4,
            5,
            5,
            5,
            6,
            6,
            6,
            7,
            8,
            9,
            9,
            10,
            10,
            11,
        ]
        Ji = [
            0,
            0,
            1,
            2,
            7,
            10,
            12,
            23,
            2,
            6,
            15,
            17,
            0,
            2,
            6,
            7,
            22,
            26,
            0,
            2,
            4,
            16,
            26,
            0,
            2,
            4,
            26,
            1,
            3,
            26,
            0,
            2,
            26,
            2,
            26,
            2,
            26,
            0,
            1,
            26,
        ]
        ni = [
            1.0658070028513,
            -15.732845290239,
            20.944396974307,
            -7.6867707878716,
            2.6185947787954,
            -2.808078114862,
            1.2053369696517,
            -8.4566812812502e-03,
            -1.2654315477714,
            -1.1524407806681,
            0.88521043984318,
            -0.64207765181607,
            0.38493460186671,
            -0.85214708824206,
            4.8972281541877,
            -3.0502617256965,
            0.039420536879154,
            0.12558408424308,
            -0.2799932969871,
            1.389979956946,
            -2.018991502357,
            -8.2147637173963e-03,
            -0.47596035734923,
            0.0439840744735,
            -0.44476435428739,
            0.90572070719733,
            0.70522450087967,
            0.10770512626332,
            -0.32913623258954,
            -0.50871062041158,
            -0.022175400873096,
            0.094260751665092,
            0.16436278447961,
            -0.013503372241348,
            -0.014834345352472,
            5.7922953628084e-04,
            3.2308904703711e-03,
            8.0964802996215e-05,
            -1.6557679795037e-04,
            -4.4923899061815e-05,
        ]
        R = SPECIFIC_GAS_CONSTANT
        tc = CRITICAL_TEMPERATURE
        rhoc = 322.0  # kg/m3
        delta = rho / rhoc
        tau = tc / T
        fi = 0
        fitau = 0
        for i in range(1, 40):
            fi = fi + ni[i] * delta ** Ii[i] * tau ** Ji[i]
            fitau = fitau + ni[i] * delta ** Ii[i] * Ji[i] * tau ** (Ji[i] - 1)
        fi = fi + ni[0] * math.log(delta)
        return R * (tau * fitau - fi)

    @staticmethod
    def Cp3_rhoT(rho: float, T: float) -> float:
        """function Cp3_rhoT = Cp3_rhoT(rho, T)

        7 Basic Equation for Region 3, Section. 6.1 Basic Equation

        Table 30 and 31, Page 30 and 31

        :param rho: density in [kg / m³]
        :param T: temperature in [K]

        :return: specific isobaric heat capacity in [kJ / (kg K)]
        """
        Ii = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            4,
            5,
            5,
            5,
            6,
            6,
            6,
            7,
            8,
            9,
            9,
            10,
            10,
            11,
        ]
        Ji = [
            0,
            0,
            1,
            2,
            7,
            10,
            12,
            23,
            2,
            6,
            15,
            17,
            0,
            2,
            6,
            7,
            22,
            26,
            0,
            2,
            4,
            16,
            26,
            0,
            2,
            4,
            26,
            1,
            3,
            26,
            0,
            2,
            26,
            2,
            26,
            2,
            26,
            0,
            1,
            26,
        ]
        ni = [
            1.0658070028513,
            -15.732845290239,
            20.944396974307,
            -7.6867707878716,
            2.6185947787954,
            -2.808078114862,
            1.2053369696517,
            -8.4566812812502e-03,
            -1.2654315477714,
            -1.1524407806681,
            0.88521043984318,
            -0.64207765181607,
            0.38493460186671,
            -0.85214708824206,
            4.8972281541877,
            -3.0502617256965,
            0.039420536879154,
            0.12558408424308,
            -0.2799932969871,
            1.389979956946,
            -2.018991502357,
            -8.2147637173963e-03,
            -0.47596035734923,
            0.0439840744735,
            -0.44476435428739,
            0.90572070719733,
            0.70522450087967,
            0.10770512626332,
            -0.32913623258954,
            -0.50871062041158,
            -0.022175400873096,
            0.094260751665092,
            0.16436278447961,
            -0.013503372241348,
            -0.014834345352472,
            5.7922953628084e-04,
            3.2308904703711e-03,
            8.0964802996215e-05,
            -1.6557679795037e-04,
            -4.4923899061815e-05,
        ]
        R = SPECIFIC_GAS_CONSTANT
        tc = CRITICAL_TEMPERATURE
        rhoc = 322.0  # kg/m3
        delta = rho / rhoc
        tau = tc / T
        fitautau = 0
        fidelta = 0
        fideltatau = 0
        fideltadelta = 0
        for i in range(1, 40):
            fitautau = fitautau + ni[i] * delta ** Ii[i] * Ji[i] * (
                Ji[i] - 1
            ) * tau ** (Ji[i] - 2)
            fidelta = fidelta + ni[i] * Ii[i] * delta ** (Ii[i] - 1) * tau ** Ji[i]
            fideltatau = fideltatau + ni[i] * Ii[i] * delta ** (Ii[i] - 1) * Ji[
                i
            ] * tau ** (Ji[i] - 1)
            fideltadelta = (
                fideltadelta
                + ni[i] * Ii[i] * (Ii[i] - 1) * delta ** (Ii[i] - 2) * tau ** Ji[i]
            )
        fidelta = fidelta + ni[0] / delta
        fideltadelta = fideltadelta - ni[0] / (delta**2)
        return R * (
            -(tau**2) * fitautau
            + (delta * fidelta - delta * tau * fideltatau) ** 2
            / (2 * delta * fidelta + delta**2 * fideltadelta)
        )

    @staticmethod
    def Cv3_rhoT(rho: float, T: float) -> float:
        """function Cv3_rhoT = Cv3_rhoT(rho, T)

        7 Basic Equation for Region 3, Section. 6.1 Basic Equation

        Table 30 and 31, Page 30 and 31

        :param rho: density in [kg / m³]
        :param T: temperature in [K]

        :return: specific isochoric heat capacity in [kJ / (kg K)]
        """
        Ii = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            4,
            5,
            5,
            5,
            6,
            6,
            6,
            7,
            8,
            9,
            9,
            10,
            10,
            11,
        ]
        Ji = [
            0,
            0,
            1,
            2,
            7,
            10,
            12,
            23,
            2,
            6,
            15,
            17,
            0,
            2,
            6,
            7,
            22,
            26,
            0,
            2,
            4,
            16,
            26,
            0,
            2,
            4,
            26,
            1,
            3,
            26,
            0,
            2,
            26,
            2,
            26,
            2,
            26,
            0,
            1,
            26,
        ]
        ni = [
            1.0658070028513,
            -15.732845290239,
            20.944396974307,
            -7.6867707878716,
            2.6185947787954,
            -2.808078114862,
            1.2053369696517,
            -8.4566812812502e-03,
            -1.2654315477714,
            -1.1524407806681,
            0.88521043984318,
            -0.64207765181607,
            0.38493460186671,
            -0.85214708824206,
            4.8972281541877,
            -3.0502617256965,
            0.039420536879154,
            0.12558408424308,
            -0.2799932969871,
            1.389979956946,
            -2.018991502357,
            -8.2147637173963e-03,
            -0.47596035734923,
            0.0439840744735,
            -0.44476435428739,
            0.90572070719733,
            0.70522450087967,
            0.10770512626332,
            -0.32913623258954,
            -0.50871062041158,
            -0.022175400873096,
            0.094260751665092,
            0.16436278447961,
            -0.013503372241348,
            -0.014834345352472,
            5.7922953628084e-04,
            3.2308904703711e-03,
            8.0964802996215e-05,
            -1.6557679795037e-04,
            -4.4923899061815e-05,
        ]
        R = SPECIFIC_GAS_CONSTANT
        tc = CRITICAL_TEMPERATURE
        rhoc = 322.0  # kg/m3
        delta = rho / rhoc
        tau = tc / T
        fitautau = 0
        # TODO:vvvv Check for mistake vvvvv
        # for i = 1 : 40
        # IAWPS says i=2..40
        for i in range(1, 40):
            fitautau = fitautau + ni[i] * delta ** Ii[i] * Ji[i] * (
                Ji[i] - 1
            ) * tau ** (Ji[i] - 2)
        return R * -(tau * tau * fitautau)

    @staticmethod
    def w3_rhoT(rho: float, T: float) -> float:
        """function w3_rhoT = w3_rhoT(rho, T)

        7 Basic Equation for Region 3, Section. 6.1 Basic Equation

        Table 30 and 31, Page 30 and 31

        :param rho: density in [kg / m³]
        :param T: temperature in [K]

        :return: speed of sound in [m / s]
        """
        Ii = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            4,
            5,
            5,
            5,
            6,
            6,
            6,
            7,
            8,
            9,
            9,
            10,
            10,
            11,
        ]
        Ji = [
            0,
            0,
            1,
            2,
            7,
            10,
            12,
            23,
            2,
            6,
            15,
            17,
            0,
            2,
            6,
            7,
            22,
            26,
            0,
            2,
            4,
            16,
            26,
            0,
            2,
            4,
            26,
            1,
            3,
            26,
            0,
            2,
            26,
            2,
            26,
            2,
            26,
            0,
            1,
            26,
        ]
        ni = [
            1.0658070028513,
            -15.732845290239,
            20.944396974307,
            -7.6867707878716,
            2.6185947787954,
            -2.808078114862,
            1.2053369696517,
            -8.4566812812502e-03,
            -1.2654315477714,
            -1.1524407806681,
            0.88521043984318,
            -0.64207765181607,
            0.38493460186671,
            -0.85214708824206,
            4.8972281541877,
            -3.0502617256965,
            0.039420536879154,
            0.12558408424308,
            -0.2799932969871,
            1.389979956946,
            -2.018991502357,
            -8.2147637173963e-03,
            -0.47596035734923,
            0.0439840744735,
            -0.44476435428739,
            0.90572070719733,
            0.70522450087967,
            0.10770512626332,
            -0.32913623258954,
            -0.50871062041158,
            -0.022175400873096,
            0.094260751665092,
            0.16436278447961,
            -0.013503372241348,
            -0.014834345352472,
            5.7922953628084e-04,
            3.2308904703711e-03,
            8.0964802996215e-05,
            -1.6557679795037e-04,
            -4.4923899061815e-05,
        ]
        R = SPECIFIC_GAS_CONSTANT
        tc = CRITICAL_TEMPERATURE
        rhoc = 322.0  # kg/m3
        delta = rho / rhoc
        tau = tc / T
        fitautau = 0
        fidelta = 0
        fideltatau = 0
        fideltadelta = 0
        for i in range(1, 40):
            fitautau = fitautau + ni[i] * delta ** Ii[i] * Ji[i] * (
                Ji[i] - 1
            ) * tau ** (Ji[i] - 2)
            fidelta = fidelta + ni[i] * Ii[i] * delta ** (Ii[i] - 1) * tau ** Ji[i]
            fideltatau = fideltatau + ni[i] * Ii[i] * delta ** (Ii[i] - 1) * Ji[
                i
            ] * tau ** (Ji[i] - 1)
            fideltadelta = (
                fideltadelta
                + ni[i] * Ii[i] * (Ii[i] - 1) * delta ** (Ii[i] - 2) * tau ** Ji[i]
            )
        fidelta = fidelta + ni[0] / delta
        fideltadelta = fideltadelta - ni[0] / (delta**2)
        return (
            1000
            * R
            * T
            * (
                2 * delta * fidelta
                + delta**2 * fideltadelta
                - (delta * fidelta - delta * tau * fideltatau) ** 2
                / (tau**2 * fitautau)
            )
        ) ** 0.5

    @staticmethod
    def T3_ph(p: float, h: float) -> float:
        """function T3_ph = T3_ph(p, h)

        Revised Supplementary Release on Backward Equations for the functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam 2004

        Section 3.3 Backward Equations T(p,h) and v(p,h) for Subregions 3a and 3b

        Boundary equation, Eq 1 Page 5

        :param p: preasure in [MPa]
        :param h: enthalpy in [kJ / kg]

        :return: temperature in [K]
        """
        h3ab = (
            2014.64004206875
            + 3.74696550136983 * p
            - 2.19921901054187e-02 * p**2
            + 8.7513168600995e-05 * p**3
        )
        if h < h3ab:
            # Subregion 3a
            # Eq 2, Table 3, Page 7
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

    @staticmethod
    def v3_ph(p: float, h: float) -> float:
        """function v3_ph = v3_ph(p, h)
        Revised Supplementary Release on Backward Equations for the functions T(p, h), v(p, h) and T(p, s), v(p, s) for Region 3 of the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam 2004

        Section 3.3 Backward Equations T(p, h) and v(p, h) for Subregions 3a and 3b

        Boundary equation, Eq 1 Page 5

        :param p: preasure in [MPa]
        :param h: enthalpy in [kJ / kg]

        :return: specific volume in [m³ / kg]
        """
        h3ab = (
            2014.64004206875
            + 3.74696550136983 * p
            - 2.19921901054187e-02 * p**2
            + 8.7513168600995e-05 * p**3
        )
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

    @staticmethod
    def T3_ps(p: float, s: float) -> float:
        """function T3_ps = T3_ps(p, s)

        Revised Supplementary Release on Backward Equations for the functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam 2004

        3.4 Backward Equations T(p,s) and v(p,s) for Subregions 3a and 3b

        Boundary equation, Eq 6 Page 11

        :param p: preasure in [MPa]
        :param s: Specific entropy in [kJ / (kg K)]

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

    @staticmethod
    def v3_ps(p: float, s: float) -> float:
        """function v3_ps = v3_ps(p, s)

        Revised Supplementary Release on Backward Equations for the functions T(p, h), v(p, h) and T(p, s), v(p, s) for Region 3 of the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam 2004

        3.4 Backward Equations T(p, s) and v(p, s) for Subregions 3a and 3b

        Boundary equation, Eq 6 Page 11

        :param p: preasure in [MPa]
        :param s: Specific entropy in [kJ / (kg K)]

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

    @staticmethod
    def p3_hs(h: float, s: float) -> float:
        """function p3_hs = p3_hs(h, s)

        Supplementary Release on Backward Equations () , p h s for Region 3, Equations as a function of h and s for the Region Boundaries, and an Equation sat , T hs for Region 4 of the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam 2004

        Section 3 Backward functions p(h, s), T(h, s), and v(h, s) for Region 3

        :param h: enthalpy in [kJ / kg]
        :param s: Specific entropy in [kJ / (kg K)]

        :return: preasure in [MPa]
        """
        if s < 4.41202148223476:
            # Subregion 3a
            # Eq 1, Table 3, Page 8
            Ii = [
                0,
                0,
                0,
                1,
                1,
                1,
                1,
                1,
                2,
                2,
                3,
                3,
                3,
                4,
                4,
                4,
                4,
                5,
                6,
                7,
                8,
                10,
                10,
                14,
                18,
                20,
                22,
                22,
                24,
                28,
                28,
                32,
                32,
            ]
            Ji = [
                0,
                1,
                5,
                0,
                3,
                4,
                8,
                14,
                6,
                16,
                0,
                2,
                3,
                0,
                1,
                4,
                5,
                28,
                28,
                24,
                1,
                32,
                36,
                22,
                28,
                36,
                16,
                28,
                36,
                16,
                36,
                10,
                28,
            ]
            ni = [
                7.70889828326934,
                -26.0835009128688,
                267.416218930389,
                17.2221089496844,
                -293.54233214597,
                614.135601882478,
                -61056.2757725674,
                -65127225.1118219,
                73591.9313521937,
                -11664650591.4191,
                35.5267086434461,
                -596.144543825955,
                -475.842430145708,
                69.6781965359503,
                335.674250377312,
                25052.6809130882,
                146997.380630766,
                5.38069315091534e19,
                1.43619827291346e21,
                3.64985866165994e19,
                -2547.41561156775,
                2.40120197096563e27,
                -3.93847464679496e29,
                1.47073407024852e24,
                -4.26391250432059e31,
                1.94509340621077e38,
                6.66212132114896e23,
                7.06777016552858e33,
                1.75563621975576e41,
                1.08408607429124e28,
                7.30872705175151e43,
                1.5914584739887e24,
                3.77121605943324e40,
            ]
            Sigma = s / 4.4
            eta = h / 2300
            Pi = 0
            for i in range(0, 33):
                Pi = Pi + ni[i] * (eta - 1.01) ** Ii[i] * (Sigma - 0.75) ** Ji[i]
            p3_hs = Pi * 99
        else:
            # Subregion 3b
            # Eq 2, Table 4, Page 8
            Ii = [
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
                -8,
                -6,
                -6,
                -6,
                -6,
                -5,
                -4,
                -4,
                -4,
                -3,
                -3,
                -3,
                -3,
                -2,
                -2,
                -1,
                0,
                2,
                2,
                5,
                6,
                8,
                10,
                14,
                14,
            ]
            Ji = [
                2,
                10,
                12,
                14,
                20,
                2,
                10,
                14,
                18,
                2,
                8,
                2,
                6,
                7,
                8,
                10,
                4,
                5,
                8,
                1,
                3,
                5,
                6,
                0,
                1,
                0,
                3,
                0,
                1,
                0,
                1,
                1,
                1,
                3,
                7,
            ]
            ni = [
                1.25244360717979e-13,
                -1.26599322553713e-02,
                5.06878030140626,
                31.7847171154202,
                -391041.161399932,
                -9.75733406392044e-11,
                -18.6312419488279,
                510.973543414101,
                373847.005822362,
                2.99804024666572e-08,
                20.0544393820342,
                -4.98030487662829e-06,
                -10.230180636003,
                55.2819126990325,
                -206.211367510878,
                -7940.12232324823,
                7.82248472028153,
                -58.6544326902468,
                3550.73647696481,
                -1.15303107290162e-04,
                -1.75092403171802,
                257.98168774816,
                -727.048374179467,
                1.21644822609198e-04,
                3.93137871762692e-02,
                7.04181005909296e-03,
                -82.910820069811,
                -0.26517881813125,
                13.7531682453991,
                -52.2394090753046,
                2405.56298941048,
                -22736.1631268929,
                89074.6343932567,
                -23923456.5822486,
                5687958081.29714,
            ]
            Sigma = s / 5.3
            eta = h / 2800
            Pi = 0
            # for i = 1 : 35
            for i in range(0, 35):
                Pi = Pi + ni[i] * (eta - 0.681) ** Ii[i] * (Sigma - 0.792) ** Ji[i]
            p3_hs = 16.6 / Pi
        return p3_hs

    @staticmethod
    def h3_pT(p: float, T: float) -> float:
        """function h3_pT = h3_pT(p, T)

        Not available with if 97

        Solve function T3_ph - T = 0 with half interval method.

        ver2.6 Start corrected bug

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: enthalpy in [kJ / kg]
        """
        logger = logging.getLogger("pyXSteam")
        if p < CRITICAL_PRESSURE:  # Below triple point
            Ts = Region4.T4_p(p)  # Saturation temperature
            if T <= Ts:  # Liquid side
                High_Bound = Region4.h4L_p(p)  # Max h ???r liauid h.
                Low_Bound = Region1.h1_pT(p, 623.15)
            else:
                Low_Bound = Region4.h4V_p(p)  # Min h ???r Vapour h.
                High_Bound = Region2.h2_pT(p, B23T_p(p))
        else:  # Above triple point. R3 from R2 till R3.
            Low_Bound = Region1.h1_pT(p, 623.15)
            High_Bound = Region2.h2_pT(p, B23T_p(p))

        Ts = T + 1
        step_counter = 0
        while math.fabs(T - Ts) > 0.00001:
            step_counter += 1
            last_Ts = Ts

            hs = (Low_Bound + High_Bound) / 2
            Ts = Region3.T3_ph(p, hs)

            if last_Ts == Ts:
                logger.warning(
                    "h3_pT stopped iterating after %d steps because values did not converge",
                    step_counter,
                )
                break

            if Ts > T:
                High_Bound = hs
            else:
                Low_Bound = hs
        return hs

    @staticmethod
    def T3_prho(p: float, rho: float) -> float:
        """function T3_prho = T3_prho(p, rho)

        Solve by iteration. Observe that of low temperatures this equation has 2 solutions. Solve with half interval method

        :param p: preasure in [MPa]
        :param rho: density in [kg / m³]

        :return: temperature in [K]
        """
        logger = logging.getLogger("pyXSteam")
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
                logger.warning(
                    "T3_prho stopped iterating after %d steps because values did not converge",
                    step_counter,
                )
                break

            if ps > p:
                High_Bound = Ts
            else:
                Low_Bound = Ts
        return Ts


class Region4:
    """
    Section 2.4: IAPWS IF 97 Calling functions to calculate the properties of water in Region 4
    """

    @staticmethod
    def p4_T(T: float) -> float:
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

    @staticmethod
    def T4_p(p: float) -> float:
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
        return (
            650.17534844798
            + D
            - (
                (650.17534844798 + D) ** 2
                - 4 * (-0.23855557567849 + 650.17534844798 * D)
            )
            ** 0.5
        ) / 2

    @staticmethod
    def h4_s(s: float) -> float:
        """function h4_s = h4_s(s)

        Supplementary Release on Backward Equations () , p h s for Region 3, Equations as a function of h and s for the Region Boundaries, and an Equation() sat , T hs for Region 4 of the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam 4 Equations for Region Boundaries Given Enthalpy and Entropy

        See picture page 14

        :param s: Specific entropy in [kJ / (kg K)]

        :return: enthalpy in [kJ / kg]
        """
        if (s > -0.0001545495919) and (s <= 3.77828134):
            # hL1_s
            # Eq 3, Table 9, Page 16
            Ii = [
                0,
                0,
                1,
                1,
                2,
                2,
                3,
                3,
                4,
                4,
                4,
                5,
                5,
                7,
                8,
                12,
                12,
                14,
                14,
                16,
                20,
                20,
                22,
                24,
                28,
                32,
                32,
            ]
            Ji = [
                14,
                36,
                3,
                16,
                0,
                5,
                4,
                36,
                4,
                16,
                24,
                18,
                24,
                1,
                4,
                2,
                4,
                1,
                22,
                10,
                12,
                28,
                8,
                3,
                0,
                6,
                8,
            ]
            ni = [
                0.332171191705237,
                6.11217706323496e-04,
                -8.82092478906822,
                -0.45562819254325,
                -2.63483840850452e-05,
                -22.3949661148062,
                -4.28398660164013,
                -0.616679338856916,
                -14.682303110404,
                284.523138727299,
                -113.398503195444,
                1156.71380760859,
                395.551267359325,
                -1.54891257229285,
                19.4486637751291,
                -3.57915139457043,
                -3.35369414148819,
                -0.66442679633246,
                32332.1885383934,
                3317.66744667084,
                -22350.1257931087,
                5739538.75852936,
                173.226193407919,
                -3.63968822121321e-02,
                8.34596332878346e-07,
                5.03611916682674,
                65.5444787064505,
            ]
            Sigma = s / 3.8
            eta = 0
            for i in range(0, 27):
                eta = (
                    eta + ni[i] * (Sigma - 1.09) ** Ii[i] * (Sigma + 0.0000366) ** Ji[i]
                )
            h4_s = eta * 1700
        elif (s > 3.77828134) and (s <= 4.41202148223476):
            # hL3_s
            # Eq 4, Table 10, Page 16
            Ii = [0, 0, 0, 0, 2, 3, 4, 4, 5, 5, 6, 7, 7, 7, 10, 10, 10, 32, 32]
            Ji = [1, 4, 10, 16, 1, 36, 3, 16, 20, 36, 4, 2, 28, 32, 14, 32, 36, 0, 6]
            ni = [
                0.822673364673336,
                0.181977213534479,
                -0.011200026031362,
                -7.46778287048033e-04,
                -0.179046263257381,
                4.24220110836657e-02,
                -0.341355823438768,
                -2.09881740853565,
                -8.22477343323596,
                -4.99684082076008,
                0.191413958471069,
                5.81062241093136e-02,
                -1655.05498701029,
                1588.70443421201,
                -85.0623535172818,
                -31771.4386511207,
                -94589.0406632871,
                -1.3927384708869e-06,
                0.63105253224098,
            ]
            Sigma = s / 3.8
            eta = 0
            for i in range(0, 19):
                eta = (
                    eta + ni[i] * (Sigma - 1.09) ** Ii[i] * (Sigma + 0.0000366) ** Ji[i]
                )
            h4_s = eta * 1700
        elif (s > 4.41202148223476) and (s <= 5.85):
            # Section 4.4 Equations () 2ab " h s and ( ) 2c3b "h s for the Saturated Vapor Line
            # Page 19, Eq 5
            # hV2c3b_s(s)
            Ii = [0, 0, 0, 1, 1, 5, 6, 7, 8, 8, 12, 16, 22, 22, 24, 36]
            Ji = [0, 3, 4, 0, 12, 36, 12, 16, 2, 20, 32, 36, 2, 32, 7, 20]
            ni = [
                1.04351280732769,
                -2.27807912708513,
                1.80535256723202,
                0.420440834792042,
                -105721.24483466,
                4.36911607493884e24,
                -328032702839.753,
                -6.7868676080427e15,
                7439.57464645363,
                -3.56896445355761e19,
                1.67590585186801e31,
                -3.55028625419105e37,
                396611982166.538,
                -4.14716268484468e40,
                3.59080103867382e18,
                -1.16994334851995e40,
            ]
            Sigma = s / 5.9
            eta = 0
            for i in range(0, 16):
                eta = eta + ni[i] * (Sigma - 1.02) ** Ii[i] * (Sigma - 0.726) ** Ji[i]
            h4_s = eta**4 * 2800
        elif (s > 5.85) and (s < 9.155759395):
            # Section 4.4 Equations () 2ab " h s and ( ) 2c3b "h s for the Saturated Vapor Line
            # Page 20, Eq 6
            Ii = [
                1,
                1,
                2,
                2,
                4,
                4,
                7,
                8,
                8,
                10,
                12,
                12,
                18,
                20,
                24,
                28,
                28,
                28,
                28,
                28,
                32,
                32,
                32,
                32,
                32,
                36,
                36,
                36,
                36,
                36,
            ]
            Ji = [
                8,
                24,
                4,
                32,
                1,
                2,
                7,
                5,
                12,
                1,
                0,
                7,
                10,
                12,
                32,
                8,
                12,
                20,
                22,
                24,
                2,
                7,
                12,
                14,
                24,
                10,
                12,
                20,
                22,
                28,
            ]
            ni = [
                -524.581170928788,
                -9269472.18142218,
                -237.385107491666,
                21077015581.2776,
                -23.9494562010986,
                221.802480294197,
                -5104725.33393438,
                1249813.96109147,
                2000084369.96201,
                -815.158509791035,
                -157.612685637523,
                -11420042233.2791,
                6.62364680776872e15,
                -2.27622818296144e18,
                -1.71048081348406e31,
                6.60788766938091e15,
                1.66320055886021e22,
                -2.18003784381501e29,
                -7.87276140295618e29,
                1.51062329700346e31,
                7957321.70300541,
                1.31957647355347e15,
                -3.2509706829914e23,
                -4.18600611419248e25,
                2.97478906557467e34,
                -9.53588761745473e19,
                1.66957699620939e24,
                -1.75407764869978e32,
                3.47581490626396e34,
                -7.10971318427851e38,
            ]
            Sigma1 = s / 5.21
            Sigma2 = s / 9.2
            eta = 0
            for i in range(0, 30):
                eta = (
                    eta
                    + ni[i] * (1 / Sigma1 - 0.513) ** Ii[i] * (Sigma2 - 0.524) ** Ji[i]
                )
            h4_s = math.exp(eta) * 2800
        else:
            h4_s = -99999
        return h4_s

    @staticmethod
    def p4_s(s: float) -> float:
        """function p4_s = p4_s(s)

        Uses h4_s and p_hs for the different regions to determine p4_s

        :param s: Specific entropy in [kJ / (kg K)]

        :return: preasure in [MPa]
        """
        h_sat = Region4.h4_s(s)
        if s > -0.0001545495919 and s <= 3.77828134:
            p4_s = Region1.p1_hs(h_sat, s)
        elif s > 3.77828134 and s <= 5.210887663:
            p4_s = Region3.p3_hs(h_sat, s)
        elif s > 5.210887663 and s < 9.155759395:
            p4_s = Region2.p2_hs(h_sat, s)
        else:
            p4_s = -99999
        return p4_s

    @staticmethod
    def h4L_p(p: float) -> float:
        """function h4L_p = h4L_p(p)

        :param p: preasure in [MPa]

        :return: enthalpy in [kJ / kg]
        """
        logger = logging.getLogger("pyXSteam")
        if (p > TRIPLE_POINT_PRESSURE) and (p < CRITICAL_PRESSURE):
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
                    ps = p3sat_h(hs)

                    if last_ps == ps:
                        logger.warning(
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

    @staticmethod
    def h4V_p(p: float) -> float:
        """function h4V_p = h4V_p(p)

        :param p: preasure in [MPa]

        :return: enthalpy in [kJ / kg]
        """
        logger = logging.getLogger("pyXSteam")
        if (p > TRIPLE_POINT_PRESSURE) and (p < CRITICAL_PRESSURE):
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
                    ps = p3sat_h(hs)

                    if last_ps == ps:
                        logger.warning(
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

    @staticmethod
    def x4_ph(p: float, h: float) -> float:
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

    @staticmethod
    def x4_ps(p: float, s: float) -> float:
        """function x4_ps = x4_ps(p, s)

        :param p: preasure in [MPa]
        :param s: Specific entropy in [kJ / (kg K)]

        :return: vapor fraction
        """
        if p < 16.529:
            ssv = Region2.s2_pT(p, Region4.T4_p(p))
            ssL = Region1.s1_pT(p, Region4.T4_p(p))
        else:
            ssv = Region3.s3_rhoT(
                1 / (Region3.v3_ph(p, Region4.h4V_p(p))), Region4.T4_p(p)
            )
            ssL = Region3.s3_rhoT(
                1 / (Region3.v3_ph(p, Region4.h4L_p(p))), Region4.T4_p(p)
            )
        if s < ssL:
            x4_ps = 0
        elif s > ssv:
            x4_ps = 1
        else:
            x4_ps = (s - ssL) / (ssv - ssL)
        return x4_ps

    @staticmethod
    def T4_hs(h: float, s: float) -> float:
        """function T4_hs = T4_hs(h, s)

        Supplementary Release on Backward Equations ( ) , p h s for Region 3,

        Chapter 5.3 page 30.

        The if 97 function is only valid for part of region4. Use iteration outside.


        :param h: enthalpy in [kJ / kg]
        :param s: Specific entropy in [kJ / (kg K)]
        :return: temperature in [K]
        """
        Ii = [
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            3,
            3,
            3,
            3,
            4,
            4,
            5,
            5,
            5,
            5,
            6,
            6,
            6,
            8,
            10,
            10,
            12,
            14,
            14,
            16,
            16,
            18,
            18,
            18,
            20,
            28,
        ]
        Ji = [
            0,
            3,
            12,
            0,
            1,
            2,
            5,
            0,
            5,
            8,
            0,
            2,
            3,
            4,
            0,
            1,
            1,
            2,
            4,
            16,
            6,
            8,
            22,
            1,
            20,
            36,
            24,
            1,
            28,
            12,
            32,
            14,
            22,
            36,
            24,
            36,
        ]
        ni = [
            0.179882673606601,
            -0.267507455199603,
            1.162767226126,
            0.147545428713616,
            -0.512871635973248,
            0.421333567697984,
            0.56374952218987,
            0.429274443819153,
            -3.3570455214214,
            10.8890916499278,
            -0.248483390456012,
            0.30415322190639,
            -0.494819763939905,
            1.07551674933261,
            7.33888415457688e-02,
            1.40170545411085e-02,
            -0.106110975998808,
            1.68324361811875e-02,
            1.25028363714877,
            1013.16840309509,
            -1.51791558000712,
            52.4277865990866,
            23049.5545563912,
            2.49459806365456e-02,
            2107964.67412137,
            366836848.613065,
            -144814105.365163,
            -1.7927637300359e-03,
            4899556021.00459,
            471.262212070518,
            -82929439019.8652,
            -1715.45662263191,
            3557776.82973575,
            586062760258.436,
            -12988763.5078195,
            31724744937.1057,
        ]
        if (s > 5.210887825) and (s < 9.15546555571324):
            Sigma = s / 9.2
            eta = h / 2800
            teta = 0
            for i in range(0, 36):
                teta = teta + ni[i] * (eta - 0.119) ** Ii[i] * (Sigma - 1.07) ** Ji[i]
            T4_hs = teta * 550
        else:
            # function psat_h
            if (s > -0.0001545495919) and (s <= 3.77828134):
                Low_Bound = 0.000611
                High_Bound = 165.291642526045
                hL = -1000
                while (math.fabs(hL - h) > 0.00001) and (
                    math.fabs(High_Bound - Low_Bound) > 0.0001
                ):
                    PL = (Low_Bound + High_Bound) / 2
                    Ts = Region4.T4_p(PL)
                    hL = Region1.h1_pT(PL, Ts)
                    if hL > h:
                        High_Bound = PL
                    else:
                        Low_Bound = PL
            elif (s > 3.77828134) and (s <= 4.41202148223476):
                PL = p3sat_h(h)
            elif (s > 4.41202148223476) and (s <= 5.210887663):
                PL = p3sat_h(h)
            Low_Bound = 0.000611
            High_Bound = PL
            sss = -1000
            while (math.fabs(s - sss) > 0.000001) and (
                math.fabs(High_Bound - Low_Bound) > 0.0000001
            ):
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
    Section 2.5: IAPWS IF 97 Calling functions to calculate the properties of water in Region 5
    """

    @staticmethod
    def h5_pT(p: float, T: float) -> float:
        """function h5_pT = h5_pT(p, T)

        Basic Equation for Region 5

        Eq 32,33, Page 36, Tables 37-41

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: enthalpy in [kJ / kg]
        """
        Ji0 = [0, 1, -3, -2, -1, 2]
        ni0 = [
            -13.179983674201,
            6.8540841634434,
            -0.024805148933466,
            0.36901534980333,
            -3.1161318213925,
            -0.32961626538917,
        ]
        Iir = [1, 1, 1, 2, 3]
        Jir = [0, 1, 3, 9, 3]
        nir = [
            -1.2563183589592e-04,
            2.1774678714571e-03,
            -0.004594282089991,
            -3.9724828359569e-06,
            1.2919228289784e-07,
        ]
        R = SPECIFIC_GAS_CONSTANT
        tau = 1000 / T
        Pi = p
        gamma0_tau = 0
        for i in range(0, 6):
            gamma0_tau = gamma0_tau + ni0[i] * Ji0[i] * tau ** (Ji0[i] - 1)
        gammar_tau = 0
        for i in range(0, 5):
            gammar_tau = gammar_tau + nir[i] * Pi ** Iir[i] * Jir[i] * tau ** (
                Jir[i] - 1
            )
        return R * T * tau * (gamma0_tau + gammar_tau)

    @staticmethod
    def v5_pT(p: float, T: float) -> float:
        """function v5_pT = v5_pT(p, T)

        Basic Equation for Region 5

        Eq 32,33, Page 36, Tables 37-41

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific volume in [m³ / kg]
        """
        Iir = [1, 1, 1, 2, 3]
        Jir = [0, 1, 3, 9, 3]
        nir = [
            -1.2563183589592e-04,
            2.1774678714571e-03,
            -0.004594282089991,
            -3.9724828359569e-06,
            1.2919228289784e-07,
        ]
        R = SPECIFIC_GAS_CONSTANT
        tau = 1000 / T  #
        Pi = p  #
        gamma0_pi = 1 / Pi  #
        gammar_pi = 0  #
        for i in range(0, 5):
            gammar_pi = gammar_pi + nir[i] * Iir[i] * Pi ** (Iir[i] - 1) * tau ** Jir[i]
        return R * T / p * Pi * (gamma0_pi + gammar_pi) / 1000

    @staticmethod
    def u5_pT(p: float, T: float) -> float:
        """function u5_pT = u5_pT(p, T)

        Basic Equation for Region 5

        Eq 32,33, Page 36, Tables 37-41

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific internal energy in [kJ / kg]
        """
        Ji0 = [0, 1, -3, -2, -1, 2]
        ni0 = [
            -13.179983674201,
            6.8540841634434,
            -0.024805148933466,
            0.36901534980333,
            -3.1161318213925,
            -0.32961626538917,
        ]
        Iir = [1, 1, 1, 2, 3]
        Jir = [0, 1, 3, 9, 3]
        nir = [
            -1.2563183589592e-04,
            2.1774678714571e-03,
            -0.004594282089991,
            -3.9724828359569e-06,
            1.2919228289784e-07,
        ]
        R = SPECIFIC_GAS_CONSTANT
        tau = 1000 / T
        Pi = p
        gamma0_pi = 1 / Pi
        gamma0_tau = 0
        for i in range(0, 6):
            gamma0_tau = gamma0_tau + ni0[i] * Ji0[i] * tau ** (Ji0[i] - 1)
        gammar_pi = 0
        gammar_tau = 0
        for i in range(0, 5):
            gammar_pi = gammar_pi + nir[i] * Iir[i] * Pi ** (Iir[i] - 1) * tau ** Jir[i]
            gammar_tau = gammar_tau + nir[i] * Pi ** Iir[i] * Jir[i] * tau ** (
                Jir[i] - 1
            )
        return R * T * (tau * (gamma0_tau + gammar_tau) - Pi * (gamma0_pi + gammar_pi))

    @staticmethod
    def Cp5_pT(p: float, T: float) -> float:
        """function Cp5_pT = Cp5_pT(p, T)

        Basic Equation for Region 5

        Eq 32,33, Page 36, Tables 37-41

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific isobaric heat capacity in [kJ / (kg K)]
        """
        Ji0 = [0, 1, -3, -2, -1, 2]
        ni0 = [
            -13.179983674201,
            6.8540841634434,
            -0.024805148933466,
            0.36901534980333,
            -3.1161318213925,
            -0.32961626538917,
        ]
        Iir = [1, 1, 1, 2, 3]
        Jir = [0, 1, 3, 9, 3]
        nir = [
            -1.2563183589592e-04,
            2.1774678714571e-03,
            -0.004594282089991,
            -3.9724828359569e-06,
            1.2919228289784e-07,
        ]
        R = SPECIFIC_GAS_CONSTANT
        tau = 1000 / T
        Pi = p
        gamma0_tautau = 0
        for i in range(0, 6):
            gamma0_tautau = gamma0_tautau + ni0[i] * Ji0[i] * (Ji0[i] - 1) * tau ** (
                Ji0[i] - 2
            )
        gammar_tautau = 0
        for i in range(0, 5):
            gammar_tautau = gammar_tautau + nir[i] * Pi ** Iir[i] * Jir[i] * (
                Jir[i] - 1
            ) * tau ** (Jir[i] - 2)
        return -R * tau**2 * (gamma0_tautau + gammar_tautau)

    @staticmethod
    def s5_pT(p: float, T: float) -> float:
        """function s5_pT = s5_pT(p, T)

        Basic Equation for Region 5

        Eq 32, 33, Page 36, Tables 37 - 41

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific entropy in [kJ / (kg K)]
        """
        Ji0 = [0, 1, -3, -2, -1, 2]
        ni0 = [
            -13.179983674201,
            6.8540841634434,
            -0.024805148933466,
            0.36901534980333,
            -3.1161318213925,
            -0.32961626538917,
        ]
        Iir = [1, 1, 1, 2, 3]
        Jir = [0, 1, 3, 9, 3]
        nir = [
            -1.2563183589592e-04,
            2.1774678714571e-03,
            -0.004594282089991,
            -3.9724828359569e-06,
            1.2919228289784e-07,
        ]
        R = SPECIFIC_GAS_CONSTANT
        tau = 1000 / T
        Pi = p
        gamma0 = math.log(Pi)
        gamma0_tau = 0
        for i in range(0, 6):
            gamma0_tau = gamma0_tau + ni0[i] * Ji0[i] * tau ** (Ji0[i] - 1)
            gamma0 = gamma0 + ni0[i] * tau ** Ji0[i]
        gammar = 0
        gammar_tau = 0
        for i in range(0, 5):
            gammar = gammar + nir[i] * Pi ** Iir[i] * tau ** Jir[i]
            gammar_tau = gammar_tau + nir[i] * Pi ** Iir[i] * Jir[i] * tau ** (
                Jir[i] - 1
            )
        return R * (tau * (gamma0_tau + gammar_tau) - (gamma0 + gammar))

    @staticmethod
    def Cv5_pT(p: float, T: float) -> float:
        """function Cv5_pT = Cv5_pT(p, T)

        Basic Equation for Region 5

        Eq 32, 33, Page 36, Tables 37 - 41

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: specific isochoric heat capacity in [kJ / (kg K)]
        """
        Ji0 = [0, 1, -3, -2, -1, 2]
        ni0 = [
            -13.179983674201,
            6.8540841634434,
            -0.024805148933466,
            0.36901534980333,
            -3.1161318213925,
            -0.32961626538917,
        ]
        Iir = [1, 1, 1, 2, 3]
        Jir = [0, 1, 3, 9, 3]
        nir = [
            -1.2563183589592e-04,
            2.1774678714571e-03,
            -0.004594282089991,
            -3.9724828359569e-06,
            1.2919228289784e-07,
        ]
        R = SPECIFIC_GAS_CONSTANT
        tau = 1000 / T
        Pi = p
        gamma0_tautau = 0
        for i in range(0, 6):
            gamma0_tautau = gamma0_tautau + ni0[i] * (Ji0[i] - 1) * Ji0[i] * tau ** (
                Ji0[i] - 2
            )
        gammar_pi = 0
        gammar_pitau = 0
        gammar_pipi = 0
        gammar_tautau = 0
        for i in range(0, 5):
            gammar_pi = gammar_pi + nir[i] * Iir[i] * Pi ** (Iir[i] - 1) * tau ** Jir[i]
            gammar_pitau = gammar_pitau + nir[i] * Iir[i] * Pi ** (Iir[i] - 1) * Jir[
                i
            ] * tau ** (Jir[i] - 1)
            gammar_pipi = (
                gammar_pipi
                + nir[i] * Iir[i] * (Iir[i] - 1) * Pi ** (Iir[i] - 2) * tau ** Jir[i]
            )
            gammar_tautau = gammar_tautau + nir[i] * Pi ** Iir[i] * Jir[i] * (
                Jir[i] - 1
            ) * tau ** (Jir[i] - 2)
        return R * (
            -(tau**2 * (gamma0_tautau + gammar_tautau))
            - (1 + Pi * gammar_pi - tau * Pi * gammar_pitau) ** 2
            / (1 - Pi**2 * gammar_pipi)
        )

    @staticmethod
    def w5_pT(p: float, T: float) -> float:
        """function w5_pT = w5_pT(p, T)

        Basic Equation for Region 5

        Eq 32, 33, Page 36, Tables 37 - 41

        :param p: preasure in [MPa]
        :param T: temperature in [K]

        :return: speed of sound in [m / s]
        """
        Ji0 = [0, 1, -3, -2, -1, 2]
        ni0 = [
            -13.179983674201,
            6.8540841634434,
            -0.024805148933466,
            0.36901534980333,
            -3.1161318213925,
            -0.32961626538917,
        ]
        Iir = [1, 1, 1, 2, 3]
        Jir = [0, 1, 3, 9, 3]
        nir = [
            -1.2563183589592e-04,
            2.1774678714571e-03,
            -0.004594282089991,
            -3.9724828359569e-06,
            1.2919228289784e-07,
        ]
        R = SPECIFIC_GAS_CONSTANT
        tau = 1000 / T
        Pi = p
        gamma0_tautau = 0
        for i in range(0, 6):
            gamma0_tautau = gamma0_tautau + ni0[i] * (Ji0[i] - 1) * Ji0[i] * tau ** (
                Ji0[i] - 2
            )
        gammar_pi = 0
        gammar_pitau = 0
        gammar_pipi = 0
        gammar_tautau = 0
        for i in range(0, 5):
            gammar_pi = gammar_pi + nir[i] * Iir[i] * Pi ** (Iir[i] - 1) * tau ** Jir[i]
            gammar_pitau = gammar_pitau + nir[i] * Iir[i] * Pi ** (Iir[i] - 1) * Jir[
                i
            ] * tau ** (Jir[i] - 1)
            gammar_pipi = (
                gammar_pipi
                + nir[i] * Iir[i] * (Iir[i] - 1) * Pi ** (Iir[i] - 2) * tau ** Jir[i]
            )
            gammar_tautau = gammar_tautau + nir[i] * Pi ** Iir[i] * Jir[i] * (
                Jir[i] - 1
            ) * tau ** (Jir[i] - 2)
        return (
            1000
            * R
            * T
            * (1 + 2 * Pi * gammar_pi + Pi**2 * gammar_pi**2)
            / (
                (1 - Pi**2 * gammar_pipi)
                + (1 + Pi * gammar_pi - tau * Pi * gammar_pitau) ** 2
                / (tau**2 * (gamma0_tautau + gammar_tautau))
            )
        ) ** 0.5

    @staticmethod
    def T5_ph(p: float, h: float) -> float:
        """function T5_ph = T5_ph(p, h)

        Solve with half interval method

        :param p: preasure in [MPa]
        :param h: enthalpy in [kJ / kg]

        :return: temperature in [K]
        """
        logger = logging.getLogger("pyXSteam")
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
                logger.warning(
                    "T5_ph stopped iterating after %d steps because values did not converge",
                    step_counter,
                )
                break

            if hs > h:
                High_Bound = Ts
            else:
                Low_Bound = Ts
        return Ts

    @staticmethod
    def T5_ps(p: float, s: float) -> float:
        """function T5_ps = T5_ps(p, s)

        Solve with half interval method

        :param p: preasure in [MPa]
        :param s: Specific entropy in [kJ / (kg K)]

        :return: temperature in [K]
        """
        logger = logging.getLogger("pyXSteam")
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
                logger.warning(
                    "T5_ps stopped iterating after %d steps because values did not converge",
                    step_counter,
                )
                break

            if ss > s:
                High_Bound = Ts
            else:
                Low_Bound = Ts
        return Ts

    @staticmethod
    def T5_prho(p: float, rho: float) -> float:
        """function T5_prho = T5_prho(p, rho)

        Solve by iteration. Observe that for low temperatures this equation has 2 solutions. Solve with half interval method

        :param p: preasure in [MPa]
        :param rho: density in [kg / m³]

        :return: temperature in [K]
        """
        logger = logging.getLogger("pyXSteam")
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
                logger.warning(
                    "T5_prho stopped iterating after %d steps because values did not converge",
                    step_counter,
                )
                break

            if rhos < rho:
                High_Bound = Ts
            else:
                Low_Bound = Ts
        return Ts
