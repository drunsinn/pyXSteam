#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""collection of demos presenting the functionality of pyXSteam"""
import time
import logging
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
import numpy as np
from pyXSteam.XSteam import XSteam
from pyXSteam.XSteam_HW import XSteam_HW


def demo_simpel_values():
    """calculate values and print the results"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    # get saturated liquid enthalpy for a preasure of 220 bar
    print('hV_p(220.0) =', steam_table.hL_p(220.0))
    # get saturated vapour enthalpy for a preasure of 220 bar
    print('hV_p(220.0) =', steam_table.hV_p(220.0))
    print('tcL_p(1.0) =', steam_table.tcL_p(1.0))
    print('tcL_t(25.0) =', steam_table.tcL_t(25.0))
    print('tcV_p(1.0) =', steam_table.tcV_p(1.0))
    print('tcL_t(25.0) =', steam_table.tcV_t(25.0))
    print('tc_hs(100.0, 0.34) =', steam_table.tc_hs(100.0, 0.34))
    print('tc_ph(1.0, 100.0) =', steam_table.tc_ph(1.0, 100.0))
    print('tc_pt(1.0, 25.0) =', steam_table.tc_pt(1.0, 25.0))
    print('w_ps(1.0, 1.0) =', steam_table.w_ps(1.0, 1.0))


def demo_generate_ph_diagramm(path=None, precision=1.0):
    """Generate a p(h) Diagramm showing the Saturation Line"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    p_krit = steam_table.criticalPressure() - 0.0001  # minus 0.0001 or else hL_V returns NaN
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
    # Siede und Taulinie
    hL = nphL_p(p)
    hV = nphV_p(p)
    # Dampfgehalt
    for vf in vapor_fraction:
        h_px = nph_px(p2, vf)
        line, = pyplot.plot(h_px, p2)
        pyplot.setp(line, linewidth=1, color='g')
    # Temperatur
    for temp in range(0, 900, 30):
        h_pt = nph_pt(p, temp)
        line, = pyplot.plot(h_pt, p)
        pyplot.setp(line, linewidth=1, color='r')
    # Dichte
    for r in rho:
        p_hrho = npp_hrho(h, r)
        line, = pyplot.plot(h, p_hrho)
        pyplot.setp(line, linewidth=1, color='y')
    # Kritischer Punkt
    pyplot.plot([h_krit], [p_krit], marker='s', mfc='k', ms=8)
    line1, = pyplot.plot(hL, p)
    line2, = pyplot.plot(hV, p)
    pyplot.xlabel("h in [kJ/kg]")
    pyplot.ylabel("p in [bar]")
    pyplot.setp(line1, linewidth=2, color='b')
    pyplot.setp(line2, linewidth=2, color='r')
    pyplot.yscale('log')
    pyplot.grid()
    if path is None:
        pyplot.show()
    else:
        pyplot.savefig(path, bbox_inches='tight')


def demo_generate_Tp_diagramm():
    """Generate a T(p) Diagramm showing the Saturation Curve"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    p = np.arange(-100.0, 250.0, 1.0)
    ntsat_p = np.frompyfunc(steam_table.tsat_p, 1, 1)
    tsat = ntsat_p(p)
    line1, = pyplot.plot(tsat, p)
    pyplot.xlabel("t")
    pyplot.ylabel("p")
    pyplot.setp(line1, linewidth=1, color='b')
    pyplot.show()


def demo_generate_pvT_diagramm():
    """Generate a Diagramm showing the v(p,T) as a 3D survace"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    fig = pyplot.figure()
    ax = Axes3D(fig)
    p = np.arange(-10.0, 300.0, 5.0)
    t = np.arange(-50.0, 400.0, 5.0)
    p, t = np.meshgrid(p, t)
    npv_pt = np.frompyfunc(steam_table.v_pt, 2, 1)
    v = npv_pt(p, t)
    ax.plot_surface(v, p, t, rstride=1, cstride=1, linewidth=0, shade=True)
    ax.set_xlabel("v")
    ax.set_ylabel("p")
    ax.set_zlabel("t")
    pyplot.show()


def demo_moillier_diagramm():
    """Generate a moillier diagramm"""
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    s = np.arange(2.0, 10.0, 0.01)
    pSteps = [0.006117, 0.01, 0.02, 1.0, 2.0, 3.0, 10, 100, 1000]
    nph_ps = np.frompyfunc(steam_table.h_ps, 2, 1)
    for pstep in pSteps:
        h = nph_ps(pstep, s)
        hline, = pyplot.plot(s, h)
        pyplot.setp(hline, linewidth=1, color='b')
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

    line1, = pyplot.plot(t_subl, psubl_func(t_subl))
    line2, = pyplot.plot(t_melt_Ih, pmelt_func(t_melt_Ih, steam_table.TYPE_ICE_Ih))
    line3, = pyplot.plot(t_melt_III, pmelt_func(t_melt_III, steam_table.TYPE_ICE_III))
    line4, = pyplot.plot(t_melt_V, pmelt_func(t_melt_V, steam_table.TYPE_ICE_V))
    line5, = pyplot.plot(t_melt_VI, pmelt_func(t_melt_VI, steam_table.TYPE_ICE_VI))
    line6, = pyplot.plot(t_melt_VII, pmelt_func(t_melt_VII, steam_table.TYPE_ICE_VII))

    pyplot.xlabel("T in [K]")
    pyplot.ylabel("p in [MPa]")
    pyplot.setp(line1, linewidth=1, color='b')
    pyplot.setp(line2, linewidth=1, color='g')
    pyplot.setp(line3, linewidth=1, color='r')
    pyplot.setp(line4, linewidth=1, color='y')
    pyplot.setp(line5, linewidth=1, color='g')
    pyplot.setp(line6, linewidth=1, color='r')
    pyplot.show()

def demo_simpel_values_heavy_water():
    """calculate values for heavy water and print the results"""
    steam_table_hw = XSteam_HW(XSteam_HW.UNIT_SYSTEM_MKS)
    print('my_rhoT(1.0, 320.0) =', steam_table_hw.my_rhoT(1.0, 320.0))
    print('my_rhoT(1.0, 320.0) =', steam_table_hw.tc_rhoT(1.0, 320.0))

if __name__ == '__main__':
    logger = logging.getLogger('pyXSteam')
    logger.setLevel(logging.DEBUG)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(sh)

    print("Select which demo to run:")
    print("1. Run simple calculations")
    print("2. generate ph diagram")
    print("3. generate Tp diagram")
    print("4. generate pvT diagramm")
    print("5. generate moillier diagram")
    print("6. generate ice metling and sublimation diagram")
    print("7. Run sinple calculations for heavy water")

    selection = str(input("Please enter selection [1-7]: "))
    print("You selected " + selection)

    if selection == '1':
        start = time.process_time()
        demo_simpel_values()
        print("Demo took", time.process_time() - start, "seconds to complete")
    elif selection == '2':
        start = time.process_time()
        demo_generate_ph_diagramm()
        print("Demo took", time.process_time() - start, "seconds to complete")
    elif selection == '3':
        start = time.process_time()
        demo_generate_Tp_diagramm()
        print("Demo took", time.process_time() - start, "seconds to complete")
    elif selection == '4':
        start = time.process_time()
        demo_generate_pvT_diagramm()
        print("Demo took", time.process_time() - start, "seconds to complete")
    elif selection == '5':
        start = time.process_time()
        demo_moillier_diagramm()
        print("Demo took", time.process_time() - start, "seconds to complete")
    elif selection == '6':
        start = time.process_time()
        demo_ice_diagramm()
        print("Demo took", time.process_time() - start, "seconds to complete")
    elif selection == '7':
        start = time.process_time()
        demo_simpel_values_heavy_water()
        print("Demo took", time.process_time() - start, "seconds to complete")
    else:
        print("Unknown selection")
