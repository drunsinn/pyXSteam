#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IAPWS R6-95(2018)
Revised Release on the IAPWS Formulation 1995 for the Thermodynamic
Properties of Ordinary Water Substance for General and Scientific Use
"""

import sys
import math
import logging

from .tables import R6_95

logger = logging.getLogger(__name__)

# TODO: make functions available via XSteam
# TODO: verify functions R6_* for which no test values are available in R6
# TODO: add missing tests


def eq_Theta(tau: float, delta: float, A: float, beta: float) -> float:
    """line 2 of eq 6 or table 5 - function not to be used separately

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param A: coefficient A_i
    :param beta: coefficient β_i

    :return: θ
    """
    return (1 - tau) + A * math.pow(math.pow(delta - 1, 2), (1 / (2 * beta)))


def eq_Delta(tau: float, delta: float, a: float, A: float, B: float, beta: float) -> float:
    """line 1 of eq 6 or table 5 - function not to be used separately

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param a: coefficient a_i
    :param A: coefficient A_i (for θ)
    :param B: coefficient B_i
    :param beta: coefficient β_i (for θ)

    :return: ∆
    """
    Theta = eq_Theta(tau, delta, A, beta)
    return math.pow(Theta, 2) + B * math.pow(math.pow(delta - 1, 2), a)


def eq_Psi(tau: float, delta: float, C: float, D: float) -> float:
    """line 3 of eq 6 or table 5 - function not to be used separately

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param C: coefficient C_i
    :param D: coefficient D_i

    :return: ψ
    """
    return math.exp(-1 * C * math.pow(delta - 1, 2) - D * math.pow(tau - 1, 2))


def eq_dDelta_ddelta(tau: float, delta: float, a: float, A: float, B: float, beta: float) -> float:
    """line 6 in table 5 continued - function not to be used separately
    Derivatives of the distance function ∆

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param a: coefficient a_i
    :param A: coefficient A_i
    :param B: coefficient B_i
    :param beta: coefficient β_i

    :return: ∂∆ / ∂δ
    """
    Theta = eq_Theta(tau, delta, A, beta)
    part_1 = A * Theta * (2 / beta) * math.pow(math.pow(delta - 1, 2), (1 / (2 * beta)) - 1)
    part_2 = 2 * B * a * math.pow(math.pow(delta - 1, 2), a - 1)
    return (delta - 1) * (part_1 + part_2)


def eq_ddDelta_ddeltadelta(tau: float, delta: float, a: float, A: float, B: float, beta: float):
    """line 7 in table 5 continued - function not to be used separately
    Derivatives of the distance function ∆

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param a: coefficient a_i
    :param A: coefficient A_i
    :param B: coefficient B_i
    :param beta: coefficient β_i

    :return: ∂^2∆ / ∂δ^2
    """

    Theta = eq_Theta(tau, delta, A, beta)
    part_1 = (1 / (delta - 1)) * eq_dDelta_ddelta(tau, delta, a, A, B, beta)
    part_2_1 = 4 * B * a * (a - 1) * math.pow(math.pow(delta - 1, 2), a - 2)
    part_2_2 = 2 * math.pow(A, 2) * math.pow(1 / beta, 2) * math.pow(math.pow(math.pow(delta - 1, 2), (1 / (2 * beta)) - 1), 2)
    part_2_3 = A * Theta * (4 / beta) * ((1 / (2 * beta)) - 1) * math.pow(math.pow(delta - 1, 2), (1 / (2 * beta)) - 2)
    return part_1 + math.pow((delta - 1), 2) * (part_2_1 + part_2_2 + part_2_3)


def eq_dDelta_b_ddelta(tau: float, delta: float, a: float, b: float, A: float, B: float, beta: float) -> float:
    """line 1 of table 5 continued - function not to be used separately
    Derivatives of the distance function ∆^b

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param a: coefficient a_i
    :param b: coefficient b_i
    :param A: coefficient A_i
    :param B: coefficient B_i
    :param beta: coefficient β_i

    :return: ∂∆^b / ∂δ
    """
    Delta = eq_Delta(tau, delta, a, A, B, beta)
    dDelta_ddelta = eq_dDelta_ddelta(tau, delta, a, A, B, beta)
    return b * math.pow(Delta, b - 1) * dDelta_ddelta


def eq_ddDelta_b_ddeltadelta(tau: float, delta: float, a: float, b: float, A: float, B: float, beta: float) -> float:
    """line 2 in table 5 continued - function not to be used separately
    Derivatives of the distance function ∆^b

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param a: coefficient a_i
    :param b: coefficient b_i
    :param A: coefficient A_i
    :param B: coefficient B_i
    :param beta: coefficient β_i

    :return: ∂^2∆^b / ∂δ^2
    """
    Delta = eq_Delta(tau, delta, a, A, B, beta)
    ddDelta_ddeltadelta = eq_ddDelta_ddeltadelta(tau, delta, a, A, B, beta)
    dDelta_ddelta = eq_dDelta_ddelta(tau, delta, a, A, B, beta)
    return b * (math.pow(Delta, b - 1) * ddDelta_ddeltadelta + (b - 1) * math.pow(Delta, b - 2) * math.pow(dDelta_ddelta, 2))


def eq_dDelta_b_dtau(tau: float, delta: float, a: float, b: float, A: float, B: float, beta: float) -> float:
    """line 3 in table 5 continued - function not to be used separately
    Derivatives of the distance function ∆^b

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param a: coefficient a_i
    :param b: coefficient b_i
    :param A: coefficient A_i
    :param B: coefficient B_i
    :param beta: coefficient β_i

    :return: ∂∆^b / ∂τ
    """
    Theta = eq_Theta(tau, delta, A, beta)
    Delta = eq_Delta(tau, delta, a, A, B, beta)
    return -2 * Theta * b * math.pow(Delta, b - 1)


def eq_ddDelta_b_dtautau(tau: float, delta: float, a: float, b: float, A: float, B: float, beta: float) -> float:
    """line 4 of table 5 continued - function not to be used separately
    Derivatives of the distance function ∆^b

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param a: coefficient a_i
    :param b: coefficient b_i
    :param A: coefficient A_i
    :param B: coefficient B_i
    :param beta: coefficient β_i

    :return: ∂^2∆^b / ∂τ^2
    """
    Theta = eq_Theta(tau, delta, A, beta)
    Delta = eq_Delta(tau, delta, a, A, B, beta)
    return 2 * b * math.pow(Delta, b - 1) + 4 * math.pow(Theta, 2) * b * (b - 1) * math.pow(Delta, b - 2)


def eq_ddDelta_b_ddeltadtau(tau: float, delta: float, a: float, b: float, A: float, B: float, beta: float) -> float:
    """line 5 of table 5 continued - function not to be used separately
    Derivatives of the distance function ∆^b

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param a: coefficient a_i
    :param b: coefficient b_i
    :param A: coefficient A_i
    :param B: coefficient B_i
    :param beta: coefficient β_i

    :return: ∂^2∆^b / ∂δ∂τ
    """

    Theta = eq_Theta(tau, delta, A, beta)
    Delta = eq_Delta(tau, delta, a, A, B, beta)
    dDelta_ddelta = eq_dDelta_ddelta(tau, delta, a, A, B, beta)
    part_1 = -1 * A * b * (2 / beta) * math.pow(Delta, b - 1) * (delta - 1) * math.pow(math.pow(delta - 1, 2), (1 / (2 * beta)) - 1)
    part_2 = -2 * Theta * b * (b - 1) * math.pow(Delta, b - 2) * dDelta_ddelta
    return part_1 + part_2


def eq_dPsi_ddelta(tau: float, delta: float, C: float, D: float) -> float:
    """line 1 in column 2 of table 5 continued - function not to be used separately
    Derivatives of the exponential function ψ

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param C: coefficient C_i
    :param D: coefficient D_i

    :return: ∂ψ / ∂δ
    """
    Psi = eq_Psi(tau, delta, C, D)
    return -2 * C * (delta - 1) * Psi


def eq_ddPsi_ddeltadelta(tau: float, delta: float, C: float, D: float) -> float:
    """line 2 in column 2 of table 5 continued - function not to be used separately
    Derivatives of the exponential function ψ

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param C: coefficient C_i
    :param D: coefficient D_i

    :return: ∂^2ψ / ∂δ^2
    """
    Psi = eq_Psi(tau, delta, C, D)
    return (2 * C * math.pow(delta - 1, 2) - 1) * 2 * C * Psi


def eq_dPsi_dtau(tau: float, delta: float, C: float, D: float) -> float:
    """line 3 in column 2 of table 5 continued - function not to be used separately
    Derivatives of the exponential function ψ

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param C: coefficient C_i
    :param D: coefficient D_i

    :return: ∂ψ / ∂τ
    """
    Psi = eq_Psi(tau, delta, C, D)
    return -2 * D * (tau - 1) * Psi


def eq_ddPsi_dtautau(tau: float, delta: float, C: float, D: float) -> float:
    """line 4 in column 2 of table 5 continued - function not to be used separately
    Derivatives of the exponential function ψ

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param C: coefficient C_i
    :param D: coefficient D_i

    :return: ∂^2ψ / ∂τ^2
    """
    Psi = eq_Psi(tau, delta, C, D)
    return (2 * D * math.pow(tau - 1, 2) - 1) * 2 * D * Psi


def eq_ddPsi_ddeltadtau(tau: float, delta: float, C: float, D: float) -> float:
    """line 5 in column 2 of table 5 continued - function not to be used separately
    Derivatives of the exponential function ψ

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ
    :param C: coefficient C_i
    :param D: coefficient D_i

    :return: ∂^2ψ / ∂δ∂τ
    """
    Psi = eq_Psi(tau, delta, C, D)
    return 4 * C * D * (delta - 1) * (tau - 1) * Psi


def eq_phi_o(tau: float, delta: float) -> float:
    """Ep 5 or line 1 in table 4 - function not to be used separately
    The ideal-gas part φo of the dimensionless Helmholtz free energy

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ

    :return: φ^o
    """
    res = math.log(delta)
    res += R6_95.Table1_n[0]
    res += R6_95.Table1_n[1] * tau
    res += R6_95.Table1_n[2] * math.log(tau)
    for n, gamma in zip(R6_95.Table1_n[3:], R6_95.Table1_gamma[3:]):
        res += n * math.log(1 - math.exp(-1 * gamma * tau))
    return res


def eq_phi_o_delta(delta: float) -> float:
    """line 2 in table 4 - function not to be used separately
    The ideal-gas part φo of the dimensionless Helmholtz free energy

    :param delta: density coefficient δ

    :return: φ^o_δ
    """
    return +1 / delta


def eq_phi_o_deltadelta(delta: float) -> float:
    """line 3 in table 4 - function not to be used separately
    The ideal-gas part φo of the dimensionless Helmholtz free energy

    :param delta: density coefficient δ

    :return: φ^o_δδ
    """
    return -1 / math.pow(delta, 2)


def eq_phi_o_tau(tau: float) -> float:
    """line 5 in table 4 - function not to be used separately
    The ideal-gas part φo of the dimensionless Helmholtz free energy

    :param tau: temperature coefficient τ

    :return: φ^o_τ
    """
    res = R6_95.Table1_n[1]
    res += R6_95.Table1_n[2] / tau
    for n, gamma in zip(R6_95.Table1_n[3:], R6_95.Table1_gamma[3:]):
        res += n * gamma * (math.pow(1 - math.exp(-1 * gamma * tau), -1) - 1)
    return res


def eq_phi_o_tautau(tau: float) -> float:
    """line 6 in table 4 - function not to be used separately
    The ideal-gas part φo of the dimensionless Helmholtz free energy

    :param tau: temperature coefficient τ

    :return: φ^o_ττ
    """
    sum_4 = R6_95.Table1_n[2] / math.pow(tau, 2)
    sum_5 = 0
    for n, gamma in zip(R6_95.Table1_n[3:], R6_95.Table1_gamma[3:]):
        # Here, the typesetting in R6 is difficult to read. correct form determined by testing
        sum_5 += n * math.pow(gamma, 2) * math.exp(-1 * gamma * tau) * math.pow(1 - math.exp(-1 * gamma * tau), -2)
    return 0 - sum_4 - sum_5


def eq_phi_o_deltatau() -> float:
    """Table 4, Eq 6 - function not to be used separately
    The ideal-gas part φo of the dimensionless Helmholtz free energy

    :return: φ^o_δτ
    """
    return 0.0


def eq_phi_r(tau: float, delta: float) -> float:
    """Eq 6 or line 1 in table 5 - function not to be used separately
    The residual part φr of the dimensionless Helmholtz free energy

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ

    :return: φ^r
    """
    res = 0
    for n, d, t in R6_95.tab2_sec1():
        res += n * math.pow(delta, d) * math.pow(tau, t)

    for n, d, t, c in R6_95.tab2_sec2():
        res += n * math.pow(delta, d) * math.pow(tau, t) * math.exp(-1 * math.pow(delta, c))

    for n, d, t, alpha, beta, gamma, epsilon in R6_95.tab2_sec3():
        exp_1 = -1 * alpha * math.pow((delta - epsilon), 2) - beta * math.pow((tau - gamma), 2)
        res += n * math.pow(delta, d) * math.pow(tau, t) * math.exp(exp_1)

    for n, a, b, A, B, C, D, beta in R6_95.tab2_sec4():
        Delta = eq_Delta(tau, delta, a, A, B, beta)
        Psi = eq_Psi(tau, delta, C, D)
        res += n * math.pow(Delta, b) * delta * Psi

    return res


def _exp(tau, delta, alpha, beta, gamma, epsilon) -> float:
    return -1 * alpha * math.pow((delta - epsilon), 2) - beta * math.pow((tau - gamma), 2)


def eq_phi_r_delta(tau: float, delta: float) -> float:
    """line 2 in table 5 - function not to be used separately
    The residual part φr of the dimensionless Helmholtz free energy

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ

    :return: φ^r_δ
    """
    res = 0
    for n, d, t in R6_95.tab2_sec1():
        res += n * d * math.pow(delta, d - 1) * math.pow(tau, t)

    for n, d, t, c in R6_95.tab2_sec2():
        res += n * math.exp(-1 * math.pow(delta, c)) * (math.pow(delta, d - 1) * math.pow(tau, t) * (d - c * math.pow(delta, c)))

    for n, d, t, alpha, beta, gamma, epsilon in R6_95.tab2_sec3():
        exp_1 = _exp(tau, delta, alpha, beta, gamma, epsilon)
        res += n * math.pow(delta, d) * math.pow(tau, t) * math.exp(exp_1) * ((d / delta) - 2.0 * alpha * (delta - epsilon))

    for n, a, b, A, B, C, D, beta in R6_95.tab2_sec4():
        Delta = eq_Delta(tau, delta, a, A, B, beta)
        dDelta_b_ddelta = eq_dDelta_b_ddelta(tau, delta, a, b, A, B, beta)
        Psi = eq_Psi(tau, delta, C, D)
        dPsi_ddelta = eq_dPsi_ddelta(tau, delta, C, D)

        res += n * (math.pow(Delta, b) * (Psi + delta * dPsi_ddelta) + dDelta_b_ddelta * delta * Psi)

    return res


def eq_phi_r_deltadelta(tau: float, delta: float) -> float:
    """line 3 in table 5 - function not to be used separately
    The residual part φr of the dimensionless Helmholtz free energy

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ

    :return: φ^r_δδ
    """
    sum_1 = 0
    for n, d, t in R6_95.tab2_sec1():
        sum_1 += n * d * (d - 1) * math.pow(delta, d - 2) * math.pow(tau, t)

    sum_2 = 0
    for n, d, t, c in R6_95.tab2_sec2():
        part_1 = math.pow(delta, d - 2) * math.pow(tau, t)
        part_2_1 = d - c * math.pow(delta, c)
        part_2_2 = d - 1 - c * math.pow(delta, c)
        part_2_3 = math.pow(c, 2) * math.pow(delta, c)
        sum_2 += n * math.exp(-1 * math.pow(delta, c)) * (part_1 * (part_2_1 * part_2_2 - part_2_3))

    sum_3 = 0
    for n, d, t, alpha, beta, gamma, epsilon in R6_95.tab2_sec3():
        exp_1 = _exp(tau, delta, alpha, beta, gamma, epsilon)
        part_1 = -2 * alpha * math.pow(delta, d)
        part_2 = 4 * alpha * math.pow(delta, d) * math.pow(delta - epsilon, 2)
        part_3 = -4 * d * alpha * math.pow(delta, d - 1) * (delta - epsilon)
        part_4 = d * (d - 1) * math.pow(delta, d - 2)
        sum_3 += n * math.pow(tau, t) * math.exp(exp_1) * (part_1 + part_2 + part_3 + part_4)

    sum_4 = 0
    for n, a, b, A, B, C, D, beta in R6_95.tab2_sec4():
        Delta = eq_Delta(tau, delta, a, A, B, beta)
        Psi = eq_Psi(tau, delta, C, D)
        dPsi_ddelta = eq_dPsi_ddelta(tau, delta, C, D)

        ddPsi_ddeltadelta = eq_ddPsi_ddeltadelta(tau, delta, C, D)
        dDelta_b_ddelta = eq_dDelta_b_ddelta(tau, delta, a, b, A, B, beta)
        ddDelta_b_ddeltadelta = eq_ddDelta_b_ddeltadelta(tau, delta, a, b, A, B, beta)

        sum_4 += n * (
            math.pow(Delta, b) * (2 * dPsi_ddelta + delta * ddPsi_ddeltadelta)
            + 2 * dDelta_b_ddelta * (Psi + delta * dPsi_ddelta)
            + ddDelta_b_ddeltadelta * delta * Psi
        )

    return sum_1 + sum_2 + sum_3 + sum_4


def eq_phi_r_tau(tau: float, delta: float) -> float:
    """line 4 in table 5 - function not to be used separately
    The residual part φr of the dimensionless Helmholtz free energy

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ

    :return: φ^r_τ
    """
    sum_1 = 0
    for n, d, t in R6_95.tab2_sec1():
        sum_1 += n * t * math.pow(delta, d) * math.pow(tau, t - 1)

    sum_2 = 0
    for n, d, t, c in R6_95.tab2_sec2():
        sum_2 += n * t * math.pow(delta, d) * math.pow(tau, t - 1) * math.exp(-1 * math.pow(delta, c))

    sum_3 = 0
    for n, d, t, alpha, beta, gamma, epsilon in R6_95.tab2_sec3():
        exp_1 = _exp(tau, delta, alpha, beta, gamma, epsilon)
        sum_3 += n * math.pow(delta, d) * math.pow(tau, t) * math.exp(exp_1) * ((t / tau) - 2 * beta * (tau - gamma))

    sum_4 = 0
    for n, a, b, A, B, C, D, beta in R6_95.tab2_sec4():
        Delta = eq_Delta(tau, delta, a, A, B, beta)
        dDelta_b_dtau = eq_dDelta_b_dtau(tau, delta, a, b, A, B, beta)
        Psi = eq_Psi(tau, delta, C, D)
        dPsi_dtau = eq_dPsi_dtau(tau, delta, C, D)
        sum_4 += n * delta * (dDelta_b_dtau * Psi + math.pow(Delta, b) * dPsi_dtau)
    return sum_1 + sum_2 + sum_3 + sum_4


def eq_phi_r_tautau(tau: float, delta: float) -> float:
    """line 5 in table 5 - function not to be used separately
    The residual part φr of the dimensionless Helmholtz free energy

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ

    :return: φ^r_ττ
    """

    sum_1 = 0
    for n, d, t in R6_95.tab2_sec1():
        sum_1 += n * t * (t - 1) * math.pow(delta, d) * math.pow(tau, t - 2)

    sum_2 = 0
    for n, d, t, c in R6_95.tab2_sec2():
        sum_2 += n * t * (t - 1) * math.pow(delta, d) * math.pow(tau, t - 2) * math.exp(-1 * math.pow(delta, c))

    sum_3 = 0
    for n, d, t, alpha, beta, gamma, epsilon in R6_95.tab2_sec3():
        exp_1 = _exp(tau, delta, alpha, beta, gamma, epsilon)
        part_1 = t / tau - 2 * beta * (tau - gamma)
        part_2 = t / math.pow(tau, 2)
        sum_3 += n * math.pow(delta, d) * math.pow(tau, t) * math.exp(exp_1) * (math.pow(part_1, 2) - part_2 - 2 * beta)

    sum_4 = 0
    for n, a, b, A, B, C, D, beta in R6_95.tab2_sec4():
        Delta = eq_Delta(tau, delta, a, A, B, beta)
        dDelta_b_dtau = eq_dDelta_b_dtau(tau, delta, a, b, A, B, beta)
        ddDelta_b_dtautau = eq_ddDelta_b_dtautau(tau, delta, a, b, A, B, beta)
        Psi = eq_Psi(tau, delta, C, D)
        dPsi_dtau = eq_dPsi_dtau(tau, delta, C, D)
        ddPsi_dtautau = eq_ddPsi_dtautau(tau, delta, C, D)

        sum_4 += n * delta * (ddDelta_b_dtautau * Psi + 2 * dDelta_b_dtau * dPsi_dtau + math.pow(Delta, b) * ddPsi_dtautau)

    return sum_1 + sum_2 + sum_3 + sum_4


def eq_phi_r_deltatau(tau: float, delta: float) -> float:
    """line 6 in table 5 - function not to be used separately
    The residual part φr of the dimensionless Helmholtz free energy

    :param tau: temperature coefficient τ
    :param delta: density coefficient δ

    :return: φ^r_δτ
    """

    sum_1 = 0
    for n, d, t in R6_95.tab2_sec1():
        sum_1 += n * d * t * math.pow(delta, d - 1) * math.pow(tau, t - 1)

    sum_2 = 0
    for n, d, t, c in R6_95.tab2_sec2():
        sum_2 += n * t * math.pow(delta, d - 1) * math.pow(tau, t - 1) * (d - c * math.pow(delta, c)) * math.exp(-1 * math.pow(delta, c))

    sum_3 = 0
    for n, d, t, alpha, beta, gamma, epsilon in R6_95.tab2_sec3():
        exp_1 = _exp(tau, delta, alpha, beta, gamma, epsilon)
        part_1 = (d / delta) - 2 * alpha * (delta - epsilon)
        part_2 = (t / tau) - 2 * beta * (tau - gamma)
        sum_3 += n * math.pow(delta, d) * math.pow(tau, t) * math.exp(exp_1) * part_1 * part_2

    sum_4 = 0
    for n, a, b, A, B, C, D, beta in R6_95.tab2_sec4():
        Delta = eq_Delta(tau, delta, a, A, B, beta)
        dDelta_b_ddelta = eq_dDelta_b_ddelta(tau, delta, a, b, A, B, beta)
        ddDelta_b_ddeltadtau = eq_ddDelta_b_ddeltadtau(tau, delta, a, b, A, B, beta)
        dDelta_b_dtau = eq_dDelta_b_dtau(tau, delta, a, b, A, B, beta)
        Psi = eq_Psi(tau, delta, C, D)
        dPsi_dtau = eq_dPsi_dtau(tau, delta, C, D)
        ddPsi_ddeltadtau = eq_ddPsi_ddeltadtau(tau, delta, C, D)
        dPsi_ddelta = eq_dPsi_ddelta(tau, delta, C, D)

        sum_4 += n * (
            math.pow(Delta, b) * (dPsi_dtau + delta * ddPsi_ddeltadtau)
            + delta * dDelta_b_ddelta * dPsi_dtau
            + dDelta_b_dtau * (Psi + delta * dPsi_ddelta)
            + ddDelta_b_ddeltadtau * delta * Psi
        )

    return sum_1 + sum_2 + sum_3 + sum_4


def R6_p_rhoT(rho: float, T: float) -> float:
    """
    line 1 in table 3, Relations of thermodynamic properties to the ideal-gas part o φ and
    the residual part r φ of the dimensionless Helmholtz free energy and their derivatives

    :param rho: density [kg / m³]
    :param T: temperature [k]

    :return: p preasure in [MPa]
    """
    delta = rho / R6_95.CRITICAL_DENSITY
    tau = R6_95.CRITICAL_TEMPERATURE / T

    if delta == 1:
        logger.info("change delta slightly to avoid div 0 at delta = 1.0")
        delta += sys.float_info.epsilon

    phi_r_delta = eq_phi_r_delta(tau, delta)

    org = 1
    org += delta * phi_r_delta

    return (org * rho * R6_95.SPECIFIC_GAS_CONSTANT * T) / 1000


def R6_u_rhoT(rho: float, T: float):
    """
    line 2 in table 3, Relations of thermodynamic properties to the ideal-gas part o φ and
    the residual part r φ of the dimensionless Helmholtz free energy and their derivatives

    :param rho: density [kg / m³]
    :param T: temperature [K]

    :return: u internal energy
    """
    delta = rho / R6_95.CRITICAL_DENSITY
    tau = R6_95.CRITICAL_TEMPERATURE / T

    if delta == 1:
        logger.info("change delta slightly to avoid div 0 at delta = 1.0")
        delta += sys.float_info.epsilon

    return tau * (eq_phi_o_tau(tau) + eq_phi_r_tau(tau, delta)) * R6_95.SPECIFIC_GAS_CONSTANT * T


def R6_s_rhoT(rho: float, T: float) -> float:
    """
    line 3 in table 3, Relations of thermodynamic properties to the ideal-gas part o φ and
    the residual part r φ of the dimensionless Helmholtz free energy and their derivatives

    :param rho: density [kg / m³]
    :param T: temperature [k]

    :return: s entropy [kJ / kg K]
    """
    delta = rho / R6_95.CRITICAL_DENSITY
    tau = R6_95.CRITICAL_TEMPERATURE / T

    if delta == 1:
        logger.info("change delta slightly to avoid div 0 at delta = 1.0")
        delta += sys.float_info.epsilon

    phi_o_tau = eq_phi_o_tau(tau)
    phi_r_tau = eq_phi_r_tau(tau, delta)
    phi_o = eq_phi_o(tau, delta)
    phi_r = eq_phi_r(tau, delta)

    org = tau * (phi_o_tau + phi_r_tau)
    org += -1 * phi_o
    org += -1 * phi_r

    return org * R6_95.SPECIFIC_GAS_CONSTANT


def R6_h_rhoT(rho: float, T: float) -> float:
    """
    line 4 in table 3, Relations of thermodynamic properties to the ideal-gas part o φ and
    the residual part r φ of the dimensionless Helmholtz free energy and their derivatives

    :param rho: density [kg / m³]
    :param T: temperature [k]

    :return: h enthalpy
    """
    delta = rho / R6_95.CRITICAL_DENSITY
    tau = R6_95.CRITICAL_TEMPERATURE / T

    if delta == 1:
        logger.info("change delta slightly to avoid div 0 at delta = 1.0")
        delta += sys.float_info.epsilon

    phi_o_tau = eq_phi_o_tau(tau)
    phi_r_tau = eq_phi_r_tau(tau, delta)
    phi_r_delta = eq_phi_r_delta(tau, delta)

    org = 1 + tau * (phi_o_tau + phi_r_tau) + delta * phi_r_delta

    return org * R6_95.SPECIFIC_GAS_CONSTANT * T


def R6_cv_rhoT(rho: float, T: float) -> float:
    """
    line 5 in table 3, Relations of thermodynamic properties to the ideal-gas part o φ and
    the residual part r φ of the dimensionless Helmholtz free energy and their derivatives

    :param rho: density [kg / m³]
    :param T: temperature [k]

    :return: cv isochoric heat capacity [kJ / kg K]
    """

    delta = rho / R6_95.CRITICAL_DENSITY
    tau = R6_95.CRITICAL_TEMPERATURE / T

    if delta == 1:
        logger.info("change delta slightly to avoid div 0 at delta = 1.0")
        delta += sys.float_info.epsilon

    phi_o_tautau = eq_phi_o_tautau(tau)
    phi_r_tautau = eq_phi_r_tautau(tau, delta)

    org = -1 * math.pow(tau, 2) * (phi_o_tautau + phi_r_tautau)

    return org * R6_95.SPECIFIC_GAS_CONSTANT


def R6_cp_rhoT(rho: float, T: float) -> float:
    """
    line 6 in table 3, Relations of thermodynamic properties to the ideal-gas part o φ and
    the residual part r φ of the dimensionless Helmholtz free energy and their derivatives

    :param rho: density [kg / m³]
    :param T: temperature [k]

    :return: cp isobaric heat capacity [kJ / kg K]
    """

    delta = rho / R6_95.CRITICAL_DENSITY
    tau = R6_95.CRITICAL_TEMPERATURE / T

    if delta == 1:
        logger.info("change delta slightly to avoid div 0 at delta = 1.0")
        delta += sys.float_info.epsilon

    phi_r_delta = eq_phi_r_delta(tau, delta)
    phi_o_tautau = eq_phi_o_tautau(tau)
    phi_r_tautau = eq_phi_r_tautau(tau, delta)
    phi_r_deltatau = eq_phi_r_deltatau(tau, delta)
    phi_r_deltadelta = eq_phi_r_deltadelta(tau, delta)

    part_1 = math.pow(1 + delta * phi_r_delta - delta * tau * phi_r_deltatau, 2)
    part_2 = 1 + 2 * delta * phi_r_delta + math.pow(delta, 2) * phi_r_deltadelta
    org = -1 * math.pow(tau, 2) * (phi_o_tautau + phi_r_tautau) + part_1 / part_2

    return org * R6_95.SPECIFIC_GAS_CONSTANT


def R6_w_rhoT(rho: float, T: float) -> float:
    """
    line 7 in table 3, Relations of thermodynamic properties to the ideal-gas part o φ and
    the residual part r φ of the dimensionless Helmholtz free energy and their derivatives


    :param rho: density [kg / m³]
    :param T: temperature [K]

    :return: w speed of sound [m / s]
    """
    delta = rho / R6_95.CRITICAL_DENSITY
    tau = R6_95.CRITICAL_TEMPERATURE / T

    if delta == 1:
        logger.info("change delta slightly to avoid div 0 at delta = 1.0")
        delta += sys.float_info.epsilon

    phi_r_delta = eq_phi_r_delta(tau, delta)
    phi_r_deltadelta = eq_phi_r_deltadelta(tau, delta)
    phi_r_tautau = eq_phi_r_tautau(tau, delta)
    phi_r_deltatau = eq_phi_r_deltatau(tau, delta)
    phi_o_tautau = eq_phi_o_tautau(tau)

    org = 1.0
    org += 2 * delta * phi_r_delta
    org += math.pow(delta, 2) * phi_r_deltadelta
    part_1 = math.pow(1 + delta * phi_r_delta - delta * tau * phi_r_deltatau, 2)
    part_2 = math.pow(tau, 2) * (phi_o_tautau + phi_r_tautau)
    org += -1 * (part_1 / part_2)

    return math.sqrt(1000 * org * R6_95.SPECIFIC_GAS_CONSTANT * T)


def R6_joulethomson_rhoT(rho: float, T: float) -> float:
    """
    line 8 in table 3, Relations of thermodynamic properties to the ideal-gas part o φ and
    the residual part r φ of the dimensionless Helmholtz free energy and their derivatives


    :param rho: density [kg / m³]
    :param T: temperature [K]

    :return: μ Joule-Thomson coefficient
    """
    delta = rho / R6_95.CRITICAL_DENSITY
    tau = R6_95.CRITICAL_TEMPERATURE / T

    if delta == 1:
        logger.info("change delta slightly to avoid div 0 at delta = 1.0")
        delta += sys.float_info.epsilon

    phi_r_delta = eq_phi_r_delta(tau, delta)
    phi_r_deltadelta = eq_phi_r_deltadelta(tau, delta)
    phi_r_tautau = eq_phi_r_tautau(tau, delta)
    phi_r_deltatau = eq_phi_r_deltatau(tau, delta)
    phi_o_tautau = eq_phi_o_tautau(tau)

    part_1 = -1 * (delta * phi_r_delta + math.pow(delta, 2) * phi_r_deltadelta + delta * tau * phi_r_deltatau)
    part_2_1 = math.pow(1 + delta * phi_r_delta - delta * tau * phi_r_deltatau, 2)
    part_2_2 = -1 * math.pow(tau, 2) * (phi_o_tautau + phi_r_tautau) * (1 + 2 * delta * phi_r_delta + math.pow(delta, 2) * phi_r_deltadelta)
    org = part_1 / (part_2_1 + part_2_2)
    return org / (R6_95.SPECIFIC_GAS_CONSTANT * rho)


def R6_delta_T_rhoT(rho: float, T: float) -> float:
    """
    line 9 in table 3, Relations of thermodynamic properties to the ideal-gas part o φ and
    the residual part r φ of the dimensionless Helmholtz free energy and their derivatives


    :param rho: density [kg / m³]
    :param T: temperature [K]

    :return: δT Isothermal throttling coefficient
    """
    delta = rho / R6_95.CRITICAL_DENSITY
    tau = R6_95.CRITICAL_TEMPERATURE / T

    if delta == 1:
        logger.info("change delta slightly to avoid div 0 at delta = 1.0")
        delta += sys.float_info.epsilon

    phi_r_delta = eq_phi_r_delta(tau, delta)
    phi_r_deltadelta = eq_phi_r_deltadelta(tau, delta)
    phi_r_deltatau = eq_phi_r_deltatau(tau, delta)

    part_1 = 1 + delta * phi_r_delta - delta * tau * phi_r_deltatau
    part_2 = 1 + 2 * delta * phi_r_delta + math.pow(delta, 2) * phi_r_deltadelta
    org = 1 - (part_1 / part_2)

    return org / rho


def R6_isen_t_p_rhoT(rho: float, T: float) -> float:
    """
    line 10 in table 3, Relations of thermodynamic properties to the ideal-gas part o φ and
    the residual part r φ of the dimensionless Helmholtz free energy and their derivatives


    :param rho: density [kg / m³]
    :param T: temperature [K]

    :return: βs  Isentropic temperature-pressure coefficient
    """
    delta = rho / R6_95.CRITICAL_DENSITY
    tau = R6_95.CRITICAL_TEMPERATURE / T

    if delta == 1:
        logger.info("change delta slightly to avoid div 0 at delta = 1.0")
        delta += sys.float_info.epsilon

    phi_r_delta = eq_phi_r_delta(tau, delta)
    phi_r_deltadelta = eq_phi_r_deltadelta(tau, delta)
    phi_r_tautau = eq_phi_r_tautau(tau, delta)
    phi_r_deltatau = eq_phi_r_deltatau(tau, delta)
    phi_o_tautau = eq_phi_o_tautau(tau)

    part_1 = 1 + delta * phi_r_delta - delta * tau * phi_r_deltatau
    part_2_1 = math.pow(1 + delta * phi_r_delta - delta * tau * phi_r_deltatau, 2)
    part_2_2 = -1 * math.pow(tau, 2) * (phi_o_tautau + phi_r_tautau) * (1 + 2 * delta * phi_r_delta + math.pow(delta, 2) * phi_r_deltadelta)

    org = part_1 / (part_2_1 + part_2_2)

    return org / (rho * R6_95.SPECIFIC_GAS_CONSTANT)


def R6_sec_virial_coefficients(rho: float, T: float) -> float:
    # """
    # Calculate the second virial coefficient and its temperature derivative at a given temperature.

    # :param T: Temperature in [K]

    # :return: A tuple containing the second virial coefficient B [m³/kg] and its temperature derivative dB/dT [m³/(kg K)]
    # """
    # tau = R6_95.CRITICAL_TEMPERATURE / T

    # B = 0.0
    # dB_dT = 0.0

    # for n, d, t in R6_95.tab2_sec1():
    #     B += n * R6_95.SPECIFIC_GAS_CONSTANT * R6_95.CRITICAL_TEMPERATURE * math.pow(R6_95.CRITICAL_DENSITY, -1) * math.pow(tau, t)
    #     dB_dT += -1 * n * t * R6_95.SPECIFIC_GAS_CONSTANT * R6_95.CRITICAL_TEMPERATURE * math.pow(R6_95.CRITICAL_DENSITY, -1) * math.pow(tau, t + 1) / T

    # return B, dB_dT
    raise NotImplementedError("Function R6_sec_virial_coefficients is not yet implemented.")


def R6_third_virial_coefficients(rho: float, T: float) -> float:
    # """
    # Calculate the third virial coefficient and its temperature derivative at a given temperature.

    # :param T: Temperature in [K]

    # :return: A tuple containing the third virial coefficient C [m^6/kg^2] and its temperature derivative dC/dT [m^6/(kg^2 K)]
    # """
    # tau = R6_95.CRITICAL_TEMPERATURE / T

    # C = 0.0
    # dC_dT = 0.0

    # for n, d, t in R6_95.tab2_sec2():
    #     C += n * R6_95.SPECIFIC_GAS_CONSTANT * R6_95.CRITICAL_TEMPERATURE * math.pow(R6_95.CRITICAL_DENSITY, -2) * math.pow(tau, t) * math.exp(-1 * math.pow(R6_95.CRITICAL_DENSITY, d))
    #     dC_dT += -1 * n * t * R6_95.SPECIFIC_GAS_CONSTANT * R6_95.CRITICAL_TEMPERATURE * math.pow(R6_95.CRITICAL_DENSITY, -2) * math.pow(tau, t + 1) * math.exp(-1 * math.pow(R6_95.CRITICAL_DENSITY, d)) / T

    # return C, dC_dT
    raise NotImplementedError("Function R6_third_virial_coefficients is not yet implemented.")


def R6_phase_equilibrium_condition(rho: float, T: float) -> float:
    """Phase-equilibrium condition (Maxwell criterion)"""

    # p_sigma / R6_95.SPECIFIC_GAS_CONSTANT * T * rho_dash = 1 + delta_dash * eq_phi_r_delta(delta_dash, tau)
    # p_sigma / R6_95.SPECIFIC_GAS_CONSTANT * T * rho_double_dash = 1 + delta_double_dash * eq_phi_r_delta(delta_double_dash, tau)
    # (p_sigma / R6_95.SPECIFIC_GAS_CONSTANT * T) * (1/ rho_double_dash - 1/ rho_dash) - math.log(delta_dash / delta_double_dash) = eq_phi_r(delta_dash, tau) - eq_phi_r(delta_double_dash, tau)

    raise NotADirectoryError()
