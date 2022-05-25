#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""collection of demos presenting the functionality of pyXSteam"""
import time
import logging
import matplotlib.pyplot as pyplot
import numpy as np
from pyXSteam.XSteam import XSteam
from pyXSteam.XSteam_HW import XSteam_HW


def demo_simpel_values():
    """calculate values and print the results"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    # get saturated liquid enthalpy for a preasure of 220 bar
    print("hV_p(220.0) =", steam_table.hL_p(220.0))
    # get saturated vapour enthalpy for a preasure of 220 bar
    print("hV_p(220.0) =", steam_table.hV_p(220.0))
    print("tcL_p(1.0) =", steam_table.tcL_p(1.0))
    print("tcL_t(25.0) =", steam_table.tcL_t(25.0))
    print("tcV_p(1.0) =", steam_table.tcV_p(1.0))
    print("tcL_t(25.0) =", steam_table.tcV_t(25.0))
    print("tc_hs(100.0, 0.34) =", steam_table.tc_hs(100.0, 0.34))
    print("tc_ph(1.0, 100.0) =", steam_table.tc_ph(1.0, 100.0))
    print("tc_pt(1.0, 25.0) =", steam_table.tc_pt(1.0, 25.0))
    print("w_ps(1.0, 1.0) =", steam_table.w_ps(1.0, 1.0))


def demo_generate_ph_diagramm(path=None, precision=1.0):
    """Generate a p(h) Diagram showing the Saturation Line"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    p_krit = (
        steam_table.criticalPressure() - 0.0001
    )  # minus 0.0001 or else hL_V returns NaN
    h_krit = steam_table.hL_p(p_krit)
    p = np.arange(0.0, 1000, precision)
    p2 = np.arange(0.5, p_krit, precision)
    vapor_fraction = np.arange(0.1, 1.0, 0.1)
    h = np.arange(200.0, 4500.0, 100.0)
    rho = np.arange(1, 15.0, precision * 2)
    nph_px = np.frompyfunc(steam_table.h_px, 2, 1)
    nph_pt = np.frompyfunc(steam_table.h_pt, 2, 1)
    nphL_p = np.frompyfunc(steam_table.hL_p, 1, 1)
    nphV_p = np.frompyfunc(steam_table.hV_p, 1, 1)
    npp_hrho = np.frompyfunc(steam_table.p_hrho, 2, 1)

    # boiling and dew lines
    hL = nphL_p(p)
    hV = nphV_p(p)

    # vapor fraction
    for vf in vapor_fraction:
        h_px = nph_px(p2, vf)
        pyplot.plot(h_px, p2, linewidth=1, color="g")

    # temperature
    for temp in range(0, 900, 30):
        h_pt = nph_pt(p, temp)
        pyplot.plot(h_pt, p, linewidth=1, color="r")

    # density
    for r in rho:
        p_hrho = npp_hrho(h, r)
        pyplot.plot(h, p_hrho, linewidth=1, color="y")

    # critical point
    pyplot.plot([h_krit], [p_krit], marker="s", mfc="k", ms=8)
    (line1,) = pyplot.plot(hL, p, linewidth=2, color="b")
    (line2,) = pyplot.plot(hV, p, linewidth=2, color="r")

    pyplot.xlabel("h in [kJ/kg]")
    pyplot.ylabel("p in [bar]")
    pyplot.yscale("log")
    pyplot.grid()
    if path is None:
        pyplot.show()
    else:
        pyplot.savefig(path, bbox_inches="tight")


def demo_generate_Tp_diagramm():
    """Generate a T(p) Diagram showing the Saturation Curve"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    p = np.arange(-100.0, 250.0, 1.0)
    ntsat_p = np.frompyfunc(steam_table.tsat_p, 1, 1)
    tsat = ntsat_p(p)
    (line1,) = pyplot.plot(tsat, p)
    pyplot.xlabel("t")
    pyplot.ylabel("p")
    pyplot.setp(line1, linewidth=1, color="b")
    pyplot.show()


def demo_generate_pvT_diagramm():
    """Generate a Diagram showing the v(p,T) as a 3D survace"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    fig = pyplot.figure()
    ax = pyplot.axes(projection="3d")

    p = np.arange(0.0, 300.0, 5.0)
    t = np.arange(120, 400.0, 5.0)
    p, t = np.meshgrid(p, t)
    npv_pt = np.frompyfunc(steam_table.v_pt, 2, 1)
    v = npv_pt(p, t)

    colour_map = pyplot.get_cmap("hot")
    ax.plot_surface(
        p, t, v, cmap=colour_map, rstride=1, cstride=1, linewidth=0, shade=True
    )

    ax.set_xlabel("p")
    ax.set_ylabel("t")
    ax.set_zlabel("v")

    pyplot.show()


def demo_moillier_diagramm():
    """Generate a moillier diagram"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    s = np.arange(2.0, 10.0, 0.01)
    pSteps = [0.006117, 0.01, 0.02, 1.0, 2.0, 3.0, 10, 100, 1000]
    nph_ps = np.frompyfunc(steam_table.h_ps, 2, 1)
    for pstep in pSteps:
        h = nph_ps(pstep, s)
        (hline,) = pyplot.plot(s, h)
        pyplot.setp(hline, linewidth=1, color="b")
    pyplot.xlabel("s in [kJ/(kg K)]")
    pyplot.ylabel("h in [kJ/kg]")
    pyplot.show()


def demo_ice_diagramm():
    """Generate a diagram showing the sublimation and melting preasure"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_BARE)
    t_subl = np.arange(50.0, 273.16, 2.0)
    t_melt_Ih = np.arange(251.165, 273.16, 2.0)
    t_melt_III = np.arange(251.165, 256.164, 2.0)
    t_melt_V = np.arange(256.164, 273.31, 2.0)
    t_melt_VI = np.arange(273.31, 355.0, 2.0)
    t_melt_VII = np.arange(355.0, 751.0, 2.0)

    psubl_func = np.frompyfunc(steam_table.psubl_t, 1, 1)
    pmelt_func = np.frompyfunc(steam_table.pmelt_t, 2, 1)

    pyplot.plot(t_subl, psubl_func(t_subl), linewidth=1, color="b", label="t_subl")
    pyplot.plot(
        t_melt_Ih,
        pmelt_func(t_melt_Ih, steam_table.TYPE_ICE_Ih),
        linewidth=2,
        color="g",
        label="t_melt_Ih",
    )
    pyplot.plot(
        t_melt_III,
        pmelt_func(t_melt_III, steam_table.TYPE_ICE_III),
        linewidth=1,
        color="r",
        label="t_melt_III",
    )
    pyplot.plot(
        t_melt_V,
        pmelt_func(t_melt_V, steam_table.TYPE_ICE_V),
        linewidth=2,
        color="y",
        label="t_melt_V",
    )
    pyplot.plot(
        t_melt_VI,
        pmelt_func(t_melt_VI, steam_table.TYPE_ICE_VI),
        linewidth=1,
        color="g",
        label="t_melt_VI",
    )
    pyplot.plot(
        t_melt_VII,
        pmelt_func(t_melt_VII, steam_table.TYPE_ICE_VII),
        linewidth=2,
        color="r",
        label="t_melt_VII",
    )

    pyplot.legend(loc="upper left")

    pyplot.xlabel("T in [K]")
    pyplot.ylabel("p in [MPa]")

    pyplot.show()


def demo_simpel_values_heavy_water():
    """calculate values for heavy water and print the results"""
    steam_table_hw = XSteam_HW(XSteam_HW.UNIT_SYSTEM_MKS)
    print("my_rhoT(1.0, 320.0) =", steam_table_hw.my_rhoT(1.0, 320.0))
    print("my_rhoT(1.0, 320.0) =", steam_table_hw.tc_rhoT(1.0, 320.0))


if __name__ == "__main__":
    logger = logging.getLogger("pyXSteam")
    logger.setLevel(logging.ERROR)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter("%(name)s - %(levelname)s - %(message)s"))
    logger.addHandler(sh)

    print("Collection of simple demos on how to use pyXSteam")
    print("requires matplotlib and numpy")
    print("------------------------------")
    print("Select which demo to run:")
    print("1. Run simple calculations")
    print("2. generate ph diagram")
    print("3. generate Tp diagram")
    print("4. generate pvT diagram")
    print("5. generate moillier diagram")
    print("6. generate ice metling and sublimation diagram")
    print("7. Run sinple calculations for heavy water")

    selection = str(input("Please enter selection [1-7]: "))
    print("You selected " + selection)
    print("------------------------------")

    if selection == "1":
        start = time.process_time()
        demo_simpel_values()
        duration = time.process_time() - start
        print("Demo took %d seconds to complete", duration)
    elif selection == "2":
        start = time.process_time()
        demo_generate_ph_diagramm()
        duration = time.process_time() - start
        print("Demo took %d seconds to complete", duration)
    elif selection == "3":
        start = time.process_time()
        demo_generate_Tp_diagramm()
        duration = time.process_time() - start
        print("Demo took %d seconds to complete", duration)
    elif selection == "4":
        start = time.process_time()
        demo_generate_pvT_diagramm()
        duration = time.process_time() - start
        print("Demo took %d seconds to complete", duration)
    elif selection == "5":
        start = time.process_time()
        demo_moillier_diagramm()
        duration = time.process_time() - start
        print("Demo took %d seconds to complete", duration)
    elif selection == "6":
        start = time.process_time()
        demo_ice_diagramm()
        duration = time.process_time() - start
        print("Demo took %d seconds to complete", duration)
    elif selection == "7":
        start = time.process_time()
        demo_simpel_values_heavy_water()
        duration = time.process_time() - start
        print("Demo took %d seconds to complete", duration)
    else:
        print("Unknown selection")

    print("------------------------------")
