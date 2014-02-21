# -*- coding: utf-8 -*-
from pyXSteam.XSteam import XSteam
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
import numpy as np

def demo_simpel_Values():
    steamTable = XSteam(mksSystem = True)
    print steamTable.hL_p(220.0)  # get saturated liquid enthalpy for a preasure of 220 bar
    print steamTable.hV_p(220.0)  # get saturated vapour enthalpy for a preasure of 220 bar

def demo_generate_ph_Diagramm(path = None, precision = 1.0):
    '''Generate a p(h) Diagramm showing the Saturation Line'''
    steamTable = XSteam(mksSystem = True)
    p_krit = steamTable.critPreasure()
    h_krit = steamTable.hL_p(p_krit)

    p = np.arange(0.0, 1000, precision)
    p2 = np.arange(0.5, p_krit, precision)
    vaporFrac = np.arange(0.1, 1.0, 0.1)
    # t = np.arange(0, 800, 1.6)

    nph_px = np.frompyfunc(steamTable.h_px, 2, 1)
    nph_pt = np.frompyfunc(steamTable.h_pt, 2, 1)
    nphL_p = np.frompyfunc(steamTable.hL_p, 1, 1)
    nphV_p = np.frompyfunc(steamTable.hV_p, 1, 1)

    hL = nphL_p(p)
    hV = nphV_p(p)

    # Dampfgehalt
    for vf in vaporFrac:
        h_px = nph_px(p2, vf)
        line, = pyplot.plot(h_px, p2)
        pyplot.setp(line, linewidth = 1, color = 'g')

    # Temperatur
    for temp in range(0, 800, 50):
        h_pt = nph_pt(p, temp)
        line, = pyplot.plot(h_pt, p)
        pyplot.setp(line, linewidth = 1, color = 'r')

    pyplot.plot([h_krit], [ p_krit], marker = 's', mfc = 'k', ms = 4)

    line1, = pyplot.plot(hL, p)
    line2, = pyplot.plot(hV, p)
    pyplot.xlabel("h")
    pyplot.ylabel("p")
    pyplot.setp(line1, linewidth = 2, color = 'b')
    pyplot.setp(line2, linewidth = 2, color = 'r')
    pyplot.yscale('log')

    if path == None:
        pyplot.show()
    else:
        pyplot.savefig(path, bbox_inches = 'tight')


def demo_generate_Tp_Diagramm():
    '''Generate a T(p) Diagramm showing the Saturation Curve'''
    steamTable = XSteam(mksSystem = True)
    p = np.arange(0.0, 250.0, 1.0)
    ntsat_p = np.frompyfunc(steamTable.tsat_p, 1, 1)
    tsat = ntsat_p(p)

    line1, = pyplot.plot(tsat, p)
    pyplot.xlabel("T")
    pyplot.ylabel("p")
    pyplot.setp(line1, linewidth = 1, color = 'b')
    pyplot.show()

def demo_generate_pvT_Diagramm():
    steamTable = XSteam(mksSystem = True)
    fig = pyplot.figure()
    ax = Axes3D(fig)

    p = np.arange(-10.0, 300.0, 5.0)
    t = np.arange(-50.0, 400.0, 5.0)
    p, t = np.meshgrid(p, t)

    npv_pt = np.frompyfunc(steamTable.v_pt, 2, 1)
    v = npv_pt(p, t)

    ax.plot_surface(v, p, t, rstride = 1, cstride = 1, linewidth = 0, shade = True)
    ax.set_xlabel("v")
    ax.set_ylabel("p")
    ax.set_zlabel("t")
    pyplot.show()

def demo_Moillier_Diagramm():
    steamTable = XSteam(mksSystem = True)
    pyplot.xlabel("s")
    pyplot.ylabel("h")
    s = np.arange(5.5 , 9.0, 0.01)
    # h = np.arange(1800 , 4200, 50)

    pSteps = [0.006117, 0.01, 0.02, 1.0, 2.0, 3.0, 10, 100, 1000]
    nph_ps = np.frompyfunc(steamTable.h_ps, 2, 1)

    for pstep in pSteps:
        print pstep
        h = nph_ps(pstep, s)
        hline, = pyplot.plot(s, h)
        pyplot.setp(hline, linewidth = 1, color = 'b')

#     h1 = nph_ps(0.006117, s)
#     line1, = pyplot.plot(s, h1)
#     h2 = nph_ps(0.01, s)
#     line2, = pyplot.plot(s, h2)
#     h3 = nph_ps(0.02, s)
#     line3, = pyplot.plot(s, h3)
#     h4 = nph_ps(1, s)
#     line4, = pyplot.plot(s, h4)
#
#     h5 = nph_ps(2, s)
#     line5, = pyplot.plot(s, h5)
#     h6 = nph_ps(3, s)
#     line6, = pyplot.plot(s, h6)
#     h7 = nph_ps(10, s)
#     line7, = pyplot.plot(s, h7)
#     h8 = nph_ps(100, s)
#     line8, = pyplot.plot(s, h8)
#     h9 = nph_ps(1000, s)
#     line9, = pyplot.plot(s, h9)
#     pyplot.setp(line2, linewidth = 1, color = 'b')
#     pyplot.setp(line3, linewidth = 1, color = 'b')
#     pyplot.setp(line4, linewidth = 1, color = 'b')
#     pyplot.setp(line5, linewidth = 1, color = 'b')
#     pyplot.setp(line6, linewidth = 1, color = 'b')
#     pyplot.setp(line7, linewidth = 1, color = 'b')
#     pyplot.setp(line8, linewidth = 1, color = 'b')
#     pyplot.setp(line9, linewidth = 1, color = 'b')
    pyplot.show()

if __name__ == '__main__':
    demo_generate_ph_Diagramm()
    # demo_generate_Tp_Diagramm()
    # demo_generate_pvT_Diagramm()
    # demo_Moillier_Diagramm()

