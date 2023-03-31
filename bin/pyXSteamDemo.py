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
    print("saturated liquid enthalpy @ 220 bar")
    print(f" hV_p(220.0) = {steam_table.hL_p(220.0)}")
    print("saturated vapour enthalpy @ 220 bar")
    print(f" hV_p(220.0) = {steam_table.hV_p(220.0)}")
    print("saturated liquid thermal conductivity @ 1 bar")
    print(f" tcL_p(1.0) = {steam_table.tcL_p(1.0)}")
    print("saturated liquid thermal conductivity @ 25 °C")
    print(f" tcL_t(25.0) = {steam_table.tcL_t(25.0)}")
    print("saturated vapour thermal conductivity @ 1 bar")
    print(f" tcV_p(1.0) = {steam_table.tcV_p(1.0)}")
    print("saturated vapour thermal conductivity @ 25 °C")
    print(f" tcL_t(25.0) = {steam_table.tcV_t(25.0)}")
    print("thermal conductivity @ enthalpy 100 kJ / kg and specific entropy 0.34 kJ / (kg °C)")
    print(f" tc_hs(100.0, 0.34) = {steam_table.tc_hs(100.0, 0.34)}")
    print("thermal conductivity @ 1 bar and enthalpy 100 kJ / kg")
    print(f" tc_ph(1.0, 100.0) = {steam_table.tc_ph(1.0, 100.0)}")
    print("thermal conductivity @ 1 bar and 25 °C")
    print(f" tc_pt(1.0, 25.0) = {steam_table.tc_pt(1.0, 25.0)}")
    print("speef of sound @ 1 bar and specific entropy 1.0 kJ / (kg °C)")
    print(f" w_ps(1.0, 1.0) = {steam_table.w_ps(1.0, 1.0)}")


def demo_generate_ph_diagramm(precision=1.0):
    """Generate a p(h) Diagram showing the Saturation Line"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)

    p_krit = (
        steam_table.criticalPressure() - 0.0001
    )  # minus 0.0001 or else hL_V returns NaN

    h_krit = steam_table.hL_p(p_krit)

    p_range = np.arange(0.0, 1000, precision)
    p2_range = np.arange(0.5, p_krit, precision)
    vf_range = np.arange(0.1, 1.0, 0.1)
    h_range = np.arange(200.0, 4500.0, 100.0)
    temp_range = np.arange(0, 900, 30)

    nph_px = np.frompyfunc(steam_table.h_px, 2, 1)
    nph_pt = np.frompyfunc(steam_table.h_pt, 2, 1)
    nphL_p = np.frompyfunc(steam_table.hL_p, 1, 1)
    nphV_p = np.frompyfunc(steam_table.hV_p, 1, 1)

    # boiling and dew lines
    hL = nphL_p(p_range)
    hV = nphV_p(p_range)

    # vapor fraction
    for i, vf in enumerate(vf_range):
        h_px = nph_px(p2_range, vf)
        if i == 0:
            pyplot.plot(h_px, p2_range, linewidth=1, color="g",
                        label="vapour fraction")
        else:
            pyplot.plot(h_px, p2_range, linewidth=1, color="g")

    # temperature
    for i, temp in enumerate(temp_range):
        h_pt = nph_pt(p_range, temp)
        if i == 0:
            pyplot.plot(h_pt, p_range, linewidth=1,
                        color="r", label="temperature")
        else:
            pyplot.plot(h_pt, p_range, linewidth=1, color="r")

    # critical point
    pyplot.plot([h_krit], [p_krit], marker="s",
                mfc="k", ms=8, label="critical point")
    (line1,) = pyplot.plot(hL, p_range, linewidth=2, color="b", label="liquide line")
    (line2,) = pyplot.plot(hV, p_range, linewidth=2, color="r", label="vapour line")

    pyplot.title("saturation Line p(h)")

    pyplot.legend(loc="upper right")

    pyplot.xlabel("h in [kJ/kg]")
    pyplot.ylabel("p in [bar]")

    pyplot.yscale("log")
    pyplot.grid()

    pyplot.show()


def demo_generate_Tp_diagramm():
    """Generate a T(p) Diagram showing the Saturation Curve"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    p = np.arange(-100.0, 250.0, 1.0)
    ntsat_p = np.frompyfunc(steam_table.tsat_p, 1, 1)
    tsat = ntsat_p(p)
    (line1,) = pyplot.plot(tsat, p)
    pyplot.setp(line1, linewidth=1, color="b")

    pyplot.title("Saturtion curve t(p)")
    pyplot.xlabel("t in [°C]")
    pyplot.ylabel("p in [bar]")
    pyplot.show()


def demo_generate_pvT_diagramm():
    """Generate a Diagram showing the v(p,T) as a 3D survace"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    #fig = pyplot.figure()
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

    pyplot.title("v(p,t)")
    ax.set_xlabel("p in [bar]")
    ax.set_ylabel("t in [°C]")
    ax.set_zlabel("v in [m³ / kg]")
    pyplot.show()


def demo_moillier_diagramm():
    """Generate a moillier diagram"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    s = np.arange(2.0, 15.0, 0.01)
    pSteps = [0.006117, 0.01, 0.02, 1.0, 2.0, 3.0, 10, 100, 1000]
    nph_ps = np.frompyfunc(steam_table.h_ps, 2, 1)
    for pstep in pSteps:
        h = nph_ps(pstep, s)
        (hline,) = pyplot.plot(s, h)
        pyplot.setp(hline, linewidth=1, label=f"{pstep} bar")

    pyplot.title("Moillier Diagram")
    pyplot.legend(loc="upper left")
    pyplot.xlabel("s in [kJ/(kg K)]")
    pyplot.ylabel("h in [kJ/kg]")
    pyplot.show()


def demo_ice_diagramm():
    """Generate a diagram showing the sublimation and melting preasure of ice"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_BARE)
    t_subl = np.arange(50.0, 273.16, 2.0)
    t_melt_Ih = np.arange(251.165, 273.16, 2.0)
    t_melt_III = np.arange(251.165, 256.164, 2.0)
    t_melt_V = np.arange(256.164, 273.31, 2.0)
    t_melt_VI = np.arange(273.31, 355.0, 2.0)
    t_melt_VII = np.arange(355.0, 751.0, 2.0)

    psubl_func = np.frompyfunc(steam_table.psubl_t, 1, 1)
    pmelt_func = np.frompyfunc(steam_table.pmelt_t, 2, 1)

    pyplot.plot(t_subl, psubl_func(t_subl),
                linewidth=1, color="b", label="sublimation temperature")
    pyplot.plot(
        t_melt_Ih,
        pmelt_func(t_melt_Ih, steam_table.TYPE_ICE_Ih),
        linewidth=2,
        color="g",
        label="melting temperature for ice type Ih",
    )
    pyplot.plot(
        t_melt_III,
        pmelt_func(t_melt_III, steam_table.TYPE_ICE_III),
        linewidth=1,
        color="r",
        label="melting temperature for ice type III",
    )
    pyplot.plot(
        t_melt_V,
        pmelt_func(t_melt_V, steam_table.TYPE_ICE_V),
        linewidth=2,
        color="y",
        label="melting temperature for ice type V",
    )
    pyplot.plot(
        t_melt_VI,
        pmelt_func(t_melt_VI, steam_table.TYPE_ICE_VI),
        linewidth=1,
        color="g",
        label="melting temperature for ice type VI",
    )
    pyplot.plot(
        t_melt_VII,
        pmelt_func(t_melt_VII, steam_table.TYPE_ICE_VII),
        linewidth=2,
        color="r",
        label="melting temperature for ice type VII",
    )

    pyplot.title("Sublimation and melting preasure of ice")
    pyplot.legend(loc="upper left")
    pyplot.xlabel("T in [K]")
    pyplot.ylabel("p in [MPa]")

    pyplot.show()


def demo_simpel_values_heavy_water():
    """calculate values for heavy water and print the results"""
    steam_table_hw = XSteam_HW(XSteam_HW.UNIT_SYSTEM_MKS)
    print("calculate values for heavy water D2O")
    print("* viscolity @ 320 °C:")
    print(f" my_rhoT(1.0, 320.0) = {steam_table_hw.my_rhoT(1.0, 320.0)}")
    print("* thermal conductivity @ 320 °C:")
    print(f" tc_rhoT(1.0, 320.0) = {steam_table_hw.tc_rhoT(1.0, 320.0)}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.ERROR)

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
        print(f"Demo took {duration:.4f} seconds to complete")
    elif selection == "2":
        start = time.process_time()
        demo_generate_ph_diagramm()
        duration = time.process_time() - start
        print(f"Demo took {duration:.4f} seconds to complete")
    elif selection == "3":
        start = time.process_time()
        demo_generate_Tp_diagramm()
        duration = time.process_time() - start
        print(f"Demo took {duration:.4f} seconds to complete")
    elif selection == "4":
        start = time.process_time()
        demo_generate_pvT_diagramm()
        duration = time.process_time() - start
        print(f"Demo took {duration:.4f} seconds to complete")
    elif selection == "5":
        start = time.process_time()
        demo_moillier_diagramm()
        duration = time.process_time() - start
        print(f"Demo took {duration:.4f} seconds to complete")
    elif selection == "6":
        start = time.process_time()
        demo_ice_diagramm()
        duration = time.process_time() - start
        print(f"Demo took {duration:.4f} seconds to complete")
    elif selection == "7":
        start = time.process_time()
        demo_simpel_values_heavy_water()
        duration = time.process_time() - start
        print(f"Demo took {duration:.4f} seconds to complete")
    else:
        print("Unknown selection")

    print("------------------------------")
