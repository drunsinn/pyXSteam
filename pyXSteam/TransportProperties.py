# -*- coding: utf-8 -*-
"""
Section 5: Transport properties
"""
import math
import logging
from . import RegionSelection
from .Regions import Region1, Region2, Region3, Region4, Region5
from . import Constants

logger = logging.getLogger(__name__)

# Section 5.1 Viscosity (IAPWS formulation 1985, Revised 2003)


def my_AllRegions_pT(p, T):
    """function my_AllRegions_pT = my_AllRegions_pT(p, T)"""
    h0 = [0.5132047, 0.3205656, 0, 0, -0.7782567, 0.1885447]
    h1 = [0.2151778, 0.7317883, 1.241044, 1.476783, 0, 0]
    h2 = [-0.2818107, -1.070786, -1.263184, 0, 0, 0]
    h3 = [0.1778064, 0.460504, 0.2340379, -0.4924179, 0, 0]
    h4 = [-0.0417661, 0, 0, 0.1600435, 0, 0]
    h5 = [0, -0.01578386, 0, 0, 0, 0]
    h6 = [0, 0, 0, -0.003629481, 0, 0]

    # %Calcualte density.
    # switch region_pT(p, T)
    if RegionSelection.region_pT(p, T) == 1:
        rho = 1 / Region1.v1_pT(p, T)
    elif RegionSelection.region_pT(p, T) == 2:
        rho = 1 / Region2.v2_pT(p, T)
    elif RegionSelection.region_pT(p, T) == 3:
        hs = Region3.h3_pT(p, T)
        rho = 1 / Region3.v3_ph(p, hs)
    elif RegionSelection.region_pT(p, T) == 4:
        logger.warning('Region switch returned unknown value')
        return float('NaN')
    elif RegionSelection.region_pT(p, T) == 5:
        rho = 1 / Region5.v5_pT(p, T)
    else:
        logger.warning('Region switch returned unknown value')
        return float("NaN")

    rhos = rho / 317.763
    Ts = T / 647.226
    # ps = p / 22.115

    # % Check valid area
    if (T > (900 + 273.15)) or ((T > (600 + 273.15)) and (p > 300)) or ((T > (150 + 273.15)) and (p > 350)) or (p > 500):
        logger.warning('Temperature and/or preasure out of range')
        return float("NaN")

    my0 = Ts ** 0.5 / (1 + 0.978197 / Ts + 0.579829 / (Ts ** 2) - 0.202354 / (Ts ** 3))
    sum = 0
    # TODO:vvvv Check for mistake? vvvvv
    # Original Code: for i = 0 : 5
    # Matlab: Index of first Element is 1
    # Python: Index of first Element is 0 -> Pythonindex = Matlabindex - 1
    # Matlab: For-loop: for i = 0 : 5 -> 0, 1, 2, 3, 4, 5
    # Python: For-loop: for i in range(0, 5): ->  0, 1, 2, 3, 4
    # The IAPWS Document says i=0..5 and j=0..6 , so range(0,6) is correct....
    for i in range(0, 6):
        sum = sum + h0[i] * (((1 / Ts) - 1) ** i) + \
                    h1[i] * (((1 / Ts) - 1) ** i) * ((rhos - 1) ** 1) + \
                    h2[i] * (((1 / Ts) - 1) ** i) * ((rhos - 1) ** 2) + \
                    h3[i] * (((1 / Ts) - 1) ** i) * ((rhos - 1) ** 3) + \
                    h4[i] * (((1 / Ts) - 1) ** i) * ((rhos - 1) ** 4) + \
                    h5[i] * (((1 / Ts) - 1) ** i) * ((rhos - 1) ** 5) + \
                    h6[i] * (((1 / Ts) - 1) ** i) * ((rhos - 1) ** 6)
    my1 = math.exp(rhos * sum)
    mys = my0 * my1
    return mys * 0.000055071


def my_AllRegions_ph(p, h):
    """ function my_AllRegions_ph = my_AllRegions_ph(p, h) """
    h0 = [0.5132047, 0.3205656, 0, 0, -0.7782567, 0.1885447]
    h1 = [0.2151778, 0.7317883, 1.241044, 1.476783, 0, 0]
    h2 = [-0.2818107, -1.070786, -1.263184, 0, 0, 0]
    h3 = [0.1778064, 0.460504, 0.2340379, -0.4924179, 0, 0]
    h4 = [-0.0417661, 0, 0, 0.1600435, 0, 0]
    h5 = [0, -0.01578386, 0, 0, 0, 0]
    h6 = [0, 0, 0, -0.003629481, 0, 0]

    # % Calcualte density.
    # switch region_ph(p, h)
    if RegionSelection.region_ph(p, h) == 1:
        Ts = Region1.T1_ph(p, h)
        T = Ts
        rho = 1 / Region1.v1_pT(p, Ts)
    elif RegionSelection.region_ph(p, h) == 2:
        Ts = Region2.T2_ph(p, h)
        T = Ts
        rho = 1 / Region2.v2_pT(p, Ts)
    elif RegionSelection.region_ph(p, h) == 3:
        rho = 1 / Region3.v3_ph(p, h)
        T = Region3.T3_ph(p, h)
    elif RegionSelection.region_ph(p, h) == 4:
        xs = Region4.x4_ph(p, h)
        if p < 16.529:
            v4v = Region2.v2_pT(p, Region4.T4_p(p))
            v4L = Region1.v1_pT(p, Region4.T4_p(p))
        else:
            v4v = Region3.v3_ph(p, Region4.h4V_p(p))
            v4L = Region3.v3_ph(p, Region4.h4L_p(p))
        rho = 1 / (xs * v4v + (1 - xs) * v4L)
        T = Region4.T4_p(p)
    elif RegionSelection.region_ph(p, h) == 5:
        Ts = Region5.T5_ph(p, h)
        T = Ts
        rho = 1 / Region5.v5_pT(p, Ts)
    else:
        logger.warning('Region switch returned unknown value')
        return float("NaN")

    rhos = rho / 317.763
    Ts = T / 647.226
    # ps = p / 22.115
    # % Check valid area
    if (T > (900 + 273.15)) or (T > (600 + 273.15) and (p > 300)) or (T > (150 + 273.15) and (p > 350)) or (p > 500):
        # my_AllRegions_ph = NaN;
        return float("NaN")

    my0 = Ts ** 0.5 / (1 + 0.978197 / Ts + 0.579829 / (Ts ** 2) - 0.202354 / (Ts ** 3))

    sum = 0
    # TODO:vvvv Check for mistake vvvvv
    # Original Code: for i = 0 : 5
    # Same Problem as in my_AllRegions_pT, see there for explanation
    for i in range(0, 6):
        sum = sum + h0[i] * (1 / Ts - 1) ** i + \
            h1[i] * (1 / Ts - 1) ** i * (rhos - 1) ** 1 + \
            h2[i] * (1 / Ts - 1) ** i * (rhos - 1) ** 2 + \
            h3[i] * (1 / Ts - 1) ** i * (rhos - 1) ** 3 + \
            h4[i] * (1 / Ts - 1) ** i * (rhos - 1) ** 4 + \
            h5[i] * (1 / Ts - 1) ** i * (rhos - 1) ** 5 + \
            h6[i] * (1 / Ts - 1) ** i * (rhos - 1) ** 6

    my1 = math.exp(rhos * sum)
    mys = my0 * my1
    return mys * 0.000055071


def tc_ptrho(p, T, rho):
    """
    function tc_ptrho = tc_ptrho(p, T, rho)

    Section  5.2 Thermal Conductivity (IAPWS formulation 1985)

    Revised release on the IAPWS formulation 1985 for the Thermal Conductivity of ordinary water IAPWS, September 1998
    Page 8
     - ver2.6 Start corrected bug
    """

    if T < 273.15:
        # tc_ptrho = NaN; % Out of range of validity (para. B4)
        logger.warning('Temperature out of range')
        return float("NaN")
    elif T < 500 + 273.15:
        if p > 100:
            # tc_ptrho = NaN; % Out of range of validity (para. B4)
            logger.warning('Preasure out of range')
            return float("NaN")
    elif T <= 650 + 273.15:
        if p > 70:
            # tc_ptrho = NaN; % Out of range of validity (para. B4)
            logger.warning('Preasure out of range')
            return float("NaN")
    else:  # T <= 800 + 273.15:
        if p > 40:
            # tc_ptrho = NaN; % Out of range of validity (para. B4)
            logger.warning('Preasure out of range')
            return float("NaN")

    # % ver2.6 End corrected bug
    T = T / 647.26  # Page 8, Eq 4
    rho = rho / 317.7  # Page 8, Eq 5

    tc0 = T ** 0.5 * (0.0102811 + 0.0299621 * T + 0.0156146 * (T ** 2) - 0.00422464 * (T ** 3))  # Page 9, Eq 9

    tc1 = -0.397070 + 0.400302 * rho + 1.06 * math.exp(-0.171587 * ((rho + 2.392190) ** 2))  # Page 9, Eq 10

    dT = abs(T - 1) + 0.00308976  # Page 9, Eq 12
    Q = 2 + 0.0822994 / (dT ** (3 / 5))  # Page 10, Eq 13
    if T >= 1:  # Page 10, Eq 14
        s = 1 / dT
    else:
        s = 10.0932 / (dT ** (3 / 5))

    tc2 = (0.0701309 / (T ** 10) + 0.0118520) * (rho ** (9 / 5)) * math.exp(0.642857 * (1 - rho ** (14 / 5))) + 0.00169937 * s * (rho ** Q) * math.exp((Q / (1 + Q)) * (1 - rho ** (1 + Q))) - 1.02 * math.exp(-4.11717 * (T ** (3 / 2)) - 6.17937 / (rho ** 5))  # Page 9, Eq 11
    return tc0 + tc1 + tc2  # Page 9, Eq 8


def Surface_Tension_T(T):
    """
    function Surface_Tension_T = Surface_Tension_T(T)

    Section 5.3 Surface Tension

    IAPWS Release on Surface Tension of Ordinary Water Substance, September 1994
    """
    # tc = 647.096  # % K
    tc = Constants.__CRITICAL_TEMPERATURE__
    B = 0.2358  # % N / m
    bb = -0.625  #
    my = 1.256  #
    if (T < 0.01) or (T > tc):
        # Surface_Tension_T = NaN# % "Out of valid region"
        logger.warning('Temperature out of range')
        return float("NaN")

    tau = 1 - T / tc
    return B * tau ** my * (1 + bb * tau)
