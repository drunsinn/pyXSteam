#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""collection of demos presenting the functionality of pyXSteam"""
import time
import logging
import math
from scipy.optimize import fsolve, fmin
from pyXSteam.XSteam import XSteam


def demo_simple_cycle():
    """
    converted example from Stu Blair
    found at https://github.com/stu314159/xsteam/blob/42767648a05c6759ad11aea95256cb24e4fc9499/Examples/SimpleRankineCycle.m
    """
    print("Purpose: test xsteam functionality for simple Rankine Cycle")
    print(
        """Problem Description:
    Advanced Boiling Water Reactor (ABWR) produces saturated steam at 7.17 MPa (71.7 bar) to a set of turbine generators.
    The turbines exhaust to condensers maintained at 8 kPa (0.08 bar).
    Assume the turbines and pumps have isentropic efficiencies of 100%
    Assume fluid leaving the condenser is saturated liquid."""
    )
    print(
        """Goal:
    Calculate enthalpy, entropy, and temperature for all state points.
    Calculate net work and thermal efficiency for this cycle"""
    )

    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)

    # state point 1: condenser outlet / feed pump inlet
    # state point 2: feed pump discharge / Reactor core inlet
    # state point 3: reactor core outlet / turbine inlet
    # state point 4: turbine exhaust / condenser inlet

    # %% State Point 1
    # condenser outlet / feed pump inlet
    state_point_1 = {}
    state_point_1["p"] = 0.08  # bar - given
    state_point_1["x"] = 0.0  # given
    state_point_1["t"] = steam_table.tsat_p(state_point_1["p"])
    state_point_1["h"] = steam_table.hL_p(state_point_1["p"])
    state_point_1["s"] = steam_table.sL_p(state_point_1["p"])

    # %% State Point 2
    # feed pump discharge / Reactor core inlet
    # process: isentropic compression
    state_point_2 = {}
    state_point_2["s"] = state_point_1["s"]
    state_point_2["p"] = 71.1  # bar - given
    state_point_2["h"] = steam_table.t_ps(state_point_2["p"], state_point_2["s"])
    state_point_2["t"] = steam_table.t_ph(state_point_2["p"], state_point_2["h"])
    # state_point_2["x"] = float("NaN")

    pump_work = state_point_1["h"] - state_point_2["h"]

    # %% State Point 3
    # reactor core outlet / turbine inlet
    # % process: isobaric heat addition
    state_point_3 = {}
    state_point_3["p"] = state_point_2["p"]
    state_point_3["x"] = 1.0
    state_point_3["t"] = steam_table.tsat_p(state_point_3["p"])
    state_point_3["h"] = steam_table.hV_p(state_point_3["p"])
    state_point_3["s"] = steam_table.sV_p(state_point_3["p"])

    heat_in = state_point_3["h"] - state_point_2["h"]

    # %% State Point 4
    # turbine exhaust / condenser inlet
    # % process: isentropic expansion
    state_point_4 = {}
    state_point_4["s"] = state_point_3["s"]
    state_point_4["p"] = state_point_1["p"]
    # % for isobaric heat rejection in next step
    state_point_4["h"] = steam_table.h_ps(state_point_4["p"], state_point_4["s"])
    state_point_4["x"] = steam_table.x_ph(state_point_4["p"], state_point_4["h"])

    turbine_work = state_point_3["h"] - state_point_4["h"]

    # % % State Point 1
    # % process: isobaric heat rejection

    heat_out = state_point_1["h"] - state_point_4["h"]

    # % % Energy Balance
    net_heat = heat_in + heat_out
    net_work = turbine_work + pump_work

    eta_th = net_work / heat_in

    print(f"Net heat: {net_heat:.3f} Net work: {net_work:.3f}")
    print(f"Thermal efficiency: {eta_th*100:.2f} percent")


def demo_ms_and_ofwh():
    """
    converted example from Stu Blair, found at https://github.com/stu314159/xsteam/blob/42767648a05c6759ad11aea95256cb24e4fc9499/Examples/Rankine_MS_and_OFWH.m
    """
    print("converted example from Stu Blair / https://github.com/stu314159")
    print("Purpose: test XSteam functionality with slightly more complex Rankine. Uses scipy.optimize.fsolve")
    print(
        """Problem Description:
    A Pressurized Water Reactor transfers heat to a Rankine cycle with the following properties:
    * Steam Generator Outlet Pressure: 820 psia, quality = 100%
    * High Pressure turbine: outlet pressure 164 psia, isentropic efficiency of 94%.
    * Moisture separator draining to OFWH (Open Feed Water Heater)
    * Flow extraction downstream of M/S set to make OFWH outlet temperature equal to saturation temp at its pressure
    * LP Turbine with outlet pressure of 3 psia, isentropic efficiency of 94%.
    * Condenser outlet quality = 0.0
    * Main condensate pump (efficiency = 84%) outlet pressure 164 psia
    * Main Feed pump (efficiency = 84%) outlet pressure 820 psia"""
    )

    steam_table = XSteam(XSteam.UNIT_SYSTEM_FLS)

    # % Calculations
    # %% state point 1 - condenser outlet
    state_point_1 = {}
    state_point_1["p"] = 3.0  # % psia - given
    state_point_1["x"] = 0.0  # % quality - given
    state_point_1["t"] = steam_table.hL_p(state_point_1["p"])
    state_point_1["h"] = steam_table.tsat_p(state_point_1["p"])
    state_point_1["s"] = steam_table.sL_p(state_point_1["p"])

    # %% state point 2
    # % compression in main condensate pump
    eta_mcp = 0.84  # % pump isentropic efficiency - given
    state_point_2 = {}
    state_point_2["s_s"] = state_point_1["s"]
    state_point_2["p"] = 164.0  # % psia - given
    state_point_2["h_s"] = steam_table.h_ps(state_point_2["p"], state_point_2["s_s"])
    state_point_2["h"] = state_point_1["h"] - (state_point_1["h"] - state_point_2["h_s"]) / eta_mcp

    # %% state point 3 OFWH exit
    state_point_3 = {}
    state_point_3["p"] = state_point_2["p"]  # % constant pressure in OFWH
    state_point_3["x"] = 0.0
    state_point_3["h"] = steam_table.hL_p(state_point_3["p"])
    state_point_3["s"] = steam_table.sL_p(state_point_3["p"])

    # %% State point 4 MFP exit
    eta_mfp = 0.84
    state_point_4 = {}
    state_point_4["p"] = 820.0  # % psia - given
    state_point_4["s_s"] = state_point_3["s"]
    state_point_4["h_s"] = steam_table.h_ps(state_point_4["p"], state_point_4["s_s"])
    state_point_4["h"] = state_point_3["h"] - (state_point_3["h"] - state_point_4["h_s"]) / eta_mfp

    # %% State point 5 S/G Exit
    state_point_5 = {}
    state_point_5["p"] = state_point_4["p"]  # % assume isobaric in S/G
    state_point_5["x"] = 1.0  # % saturated steam; given
    state_point_5["h"] = steam_table.hV_p(state_point_5["p"])
    state_point_5["s"] = steam_table.sV_p(state_point_5["p"])

    # %% State point 6 HP Turbine Exhaust
    eta_hpt = 0.94  # % hp turbine isentropic efficiency; given
    state_point_6 = {}
    state_point_6["p"] = 164.0  # % psia - given
    state_point_6["s_s"] = state_point_5["s"]
    state_point_6["h_s"] = steam_table.h_ps(state_point_6["p"], state_point_6["s_s"])
    state_point_6["h"] = state_point_5["h"] - eta_hpt * (state_point_5["h"] - state_point_6["h_s"])
    state_point_6["x"] = steam_table.x_ph(state_point_6["p"], state_point_6["h"])

    # %% State point 7 Moisture Separator vapor exit
    state_point_7 = {}
    state_point_7["p"] = state_point_6["p"]  # assume isobaric process in M/S
    state_point_7["x"] = 1.0  # quality - given
    state_point_7["h"] = steam_table.hV_p(state_point_7["p"])
    state_point_7["s"] = steam_table.sV_p(state_point_7["p"])

    # %% State point 8 LP Turbine exhaust
    eta_lpt = 0.94  # % lp turbine isentropic efficiency; given
    state_point_8 = {}
    state_point_8["p"] = state_point_1["p"]
    state_point_8["s_s"] = state_point_7["s"]
    state_point_8["h_s"] = steam_table.h_ps(state_point_8["p"], state_point_8["s_s"])
    state_point_8["h"] = state_point_7["h"] - eta_lpt * (state_point_7["h"] - state_point_8["h_s"])
    state_point_8["x"] = steam_table.x_ph(state_point_8["p"], state_point_8["h"])

    # %% State point 9 Moisture Separator liquid drain to OFWH
    state_point_9 = {}
    # % same pressure as HP Turbine exhaust
    state_point_9["p"] = state_point_6["p"]
    state_point_9["h"] = steam_table.hL_p(state_point_9["p"])

    # http://mathesaurus.sourceforge.net/matlab-numpy.html

    # %% Energy balance on OFWH to find flow fraction f1 at extraction point
    # OFWH_Ebal = @(f1) x(6)*(1-f1)*h(2)+x(6)*f1*h(7)+(1-x(6))*h(9) - h(3)
    def OFWH_Ebal(f1):
        # https://de.mathworks.com/help/matlab/matlab_prog/matlab-operators-and-special-characters.html
        # x(6)*(1-f1)*h(2)+x(6)*f1*h(7)+(1-x(6))*h(9) - h(3)
        return (
            state_point_6["x"] * (1 - f1) * state_point_2["h"]
            + state_point_6["x"] * f1 * state_point_7["h"]
            + (1 - state_point_6["x"]) * state_point_9["h"]
            - state_point_3["h"]
        )

    f1 = fsolve(OFWH_Ebal, 0.05)[0]

    # %% Specific Work and Energy Balance
    w_mcp = (state_point_1["h"] - state_point_2["h"]) * f1 * state_point_6["x"]
    w_mfp = state_point_3["h"] - state_point_4["h"]
    w_hpt = state_point_5["h"] - state_point_6["h"]
    w_lpt = (state_point_7["h"] - state_point_8["h"]) * (1 - f1) * state_point_6["x"]

    w_net = w_mcp + w_mfp + w_hpt + w_lpt

    q_cond = (state_point_1["h"] - state_point_8["h"]) * (1 - f1) * state_point_6["x"]
    q_sg = state_point_5["h"] - state_point_4["h"]

    q_net = q_cond + q_sg
    eta_th = w_net / q_sg

    print(f"Net heat: {q_net:.3f} BTU/lbm Net work: {w_net:.3f} BTU/lbm")
    print(f"Thermal efficiency: {eta_th*100:.2f} percent")


def demo_reheat_ms_ofwh():
    """
    converted example from Stu Blair, found at https://github.com/stu314159/xsteam/blob/42767648a05c6759ad11aea95256cb24e4fc9499/Examples/Rankine_Reheat_MS_OFWH.m
    """
    print("converted example from Stu Blair / https://github.com/stu314159")
    print("Purpose: use fminsearch / scipy.optimize.fmin")
    print(
        """Problem Description:
    (a picture would be better....)
    """
    )

    steam_table = XSteam(XSteam.UNIT_SYSTEM_FLS)

    # %% State point 1 - condenser exit
    state_point_1 = {}
    state_point_1["p"] = 1.5  # % psia
    state_point_1["h"] = steam_table.hL_p(state_point_1["p"])
    state_point_1["s"] = steam_table.sL_p(state_point_1["p"])

    # %% State point 1 -> 2, main condensate pump
    eta_mcp = 0.84
    state_point_2 = {}
    state_point_2["p"] = 164.0  # % psia
    state_point_2["s_s"] = state_point_1["s"]
    state_point_2["h_s"] = steam_table.h_ps(state_point_2["p"], state_point_2["s_s"])
    state_point_2["h"] = state_point_1["h"] - (state_point_1["h"] - state_point_2["h_s"]) / eta_mcp

    # %% State point 3, OFWH exit, saturated liquid
    state_point_3 = {}
    state_point_3["p"] = state_point_2["p"]
    state_point_3["h"] = steam_table.hL_p(state_point_3["p"])
    state_point_3["s"] = steam_table.sL_p(state_point_3["p"])

    # %% State point 3 -> 4, main feed pump
    eta_mfp = 0.84
    state_point_4 = {}
    state_point_4["p"] = 820.0  # % psia
    state_point_4["s_s"] = state_point_3["s"]
    state_point_4["h_s"] = steam_table.h_ps(state_point_4["p"], state_point_4["s_s"])
    state_point_4["h"] = state_point_3["h"] - (state_point_3["h"] - state_point_4["h_s"]) / eta_mfp

    # %% State point 5 - Steam generator exit
    state_point_5 = {}
    state_point_5["p"] = state_point_4["p"]
    state_point_5["h"] = steam_table.hV_p(state_point_5["p"])
    state_point_5["s"] = steam_table.sV_p(state_point_5["p"])

    # %% State point 6 - HP Turbine Exhaust
    eta_hpt = 0.94
    state_point_6 = {}
    state_point_6["p"] = 164.0  # % psia
    state_point_6["s_s"] = state_point_5["s"]
    state_point_6["h_s"] = steam_table.h_ps(state_point_6["p"], state_point_6["s_s"])
    state_point_6["h"] = state_point_5["h"] - eta_hpt * (state_point_5["h"] - state_point_6["h_s"])
    state_point_6["x"] = steam_table.x_ph(state_point_6["p"], state_point_6["h"])

    # %% State point 7 - Moisture Separator Exit
    state_point_7 = {}
    state_point_7["p"] = state_point_6["p"]
    state_point_7["h"] = steam_table.hV_p(state_point_7["p"])
    state_point_7["s"] = steam_table.sV_p(state_point_7["p"])

    # %% State point 8 - Reheater Mid-Pressure Steam exit
    state_point_8 = {}
    state_point_8["p"] = state_point_7["p"]
    state_point_8["t"] = 490.0  # % degrees F
    state_point_8["h"] = steam_table.h_pt(state_point_8["p"], state_point_8["t"])
    state_point_8["s"] = steam_table.s_pt(state_point_8["p"], state_point_8["t"])

    # %% State point 9 - LP Turbine Exhaust
    state_point_9 = {}
    state_point_9["p"] = state_point_1["p"]
    eta_lpt = 0.94
    state_point_9["s_s"] = state_point_8["s"]
    state_point_9["h_s"] = steam_table.h_ps(state_point_9["p"], state_point_9["s_s"])
    state_point_9["h"] = state_point_8["h"] - eta_lpt * (state_point_8["h"] - state_point_9["h_s"])

    # %% State point 10 - Reheater HP Steam exit
    state_point_10 = {}
    state_point_10["p"] = state_point_5["p"]
    # % assume steam exits as a saturated liquid.
    state_point_10["h"] = steam_table.hL_p(state_point_10["p"])

    # %% State point 11 - pressure trap exit to OFWH
    state_point_11 = {}
    state_point_11["p"] = state_point_2["p"]
    # % assume isenthalpic expansion in the trap.
    state_point_11["h"] = state_point_10["h"]

    # %% State point 12 - Moisture Separator liquid drain to OFWH
    state_point_12 = {}
    state_point_12["p"] = state_point_6["p"]
    state_point_12["h"] = steam_table.hL_p(state_point_12["p"])

    # %% Heat Balance - find the flow fractions
    def RH_heatBalance(f):
        return (f[0] * state_point_5["h"] + (1 - f[0]) * state_point_6["x"] * (1 - f[1]) * state_point_7["h"]) - (
            f[0] * state_point_10["h"] + (1 - f[0]) * state_point_6["x"] * (1 - f[1]) * state_point_8["h"]
        )

    def OFWH_heatBalance(f):
        return (
            (1 - f[0]) * (1 - f[1]) * state_point_6["x"] * state_point_2["h"]
            + f[0] * state_point_11["h"]
            + (1 - f[0]) * state_point_6["x"] * f[1] * state_point_7["h"]
            + (1 - f[0]) * (1 - state_point_6["x"]) * state_point_12["h"]
        ) - state_point_3["h"]

    def totalFunctional(f):
        return math.fabs(RH_heatBalance(f)) + math.fabs(OFWH_heatBalance(f))

    # % the strategy is to minimize the total functional.  The minimum value is
    # % when they are both equal to zero.
    initialGuess = [0.1, 0.1]
    f = fmin(func=totalFunctional, x0=initialGuess)

    # % % calculate heat and energy balances
    w_mcp = (state_point_1["h"] - state_point_2["h"]) * (1 - f[0]) * (1 - f[1]) * state_point_6["x"]
    w_mfp = state_point_3["h"] - state_point_4["h"]
    w_hpt = (state_point_5["h"] - state_point_6["h"]) * (1 - f[0])
    w_lpt = (state_point_8["h"] - state_point_9["h"]) * (1 - f[0]) * (1 - f[1]) * state_point_6["x"]

    w_net = w_mcp + w_mfp + w_hpt + w_lpt

    q_cond = (state_point_1["h"] - state_point_9["h"]) * (1 - f[0]) * (1 - f[1]) * state_point_6["x"]
    q_sg = state_point_5["h"] - state_point_4["h"]

    q_net = q_cond + q_sg

    eta_th = w_net / q_sg

    # % % report the results:
    print(f"Net heat: {q_net:.3f} BTU/lbm Net work: {w_net:.3f} BTU/lbm")
    print(f"Thermal efficiency: {eta_th*100:.2f} percent")


def main():
    logger = logging.getLogger("pyXSteam")
    logger.setLevel(logging.ERROR)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter("%(name)s - %(levelname)s - %(message)s"))
    logger.addHandler(sh)

    print("Demos on how to use pyXSteam on the example of the rankine cycle")
    print("converted examples from Stu Blair")
    print("https://github.com/stu314159/xsteam")
    print("requires scipy")
    print("------------------------------")
    print("Select which demo to run:")
    print("1. Simple Rankine Cycle")
    print("2. Rankine cycle with water separator and open feed water heater")
    print("3. Reheat Rankine cycle with water separator and open feed water heater")

    selection = str(input("Please enter selection [1-3]: "))
    print("You selected " + selection)
    print("------------------------------")

    if selection == "1":
        start = time.process_time()
        demo_simple_cycle()
        duration = time.process_time() - start
        print(f"Demo took {duration:.4f} seconds to complete")
    elif selection == "2":
        start = time.process_time()
        demo_ms_and_ofwh()
        duration = time.process_time() - start
        print(f"Demo took {duration:.4f} seconds to complete")
    elif selection == "3":
        start = time.process_time()
        demo_reheat_ms_ofwh()
        duration = time.process_time() - start
        print(f"Demo took {duration:.4f} seconds to complete")
    else:
        print("Unknown selection")

    print("------------------------------")


if __name__ == "__main__":
    main()
