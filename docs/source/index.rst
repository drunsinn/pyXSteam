.. pyXSteam documentation master file, created by
   sphinx-quickstart on Fri Mar  8 16:32:11 2019.

Welcome to pyXSteam's documentation!
####################################

Original Released by Magnus Holmgren for Matlab and Excel:
<http://xsteam.sourceforge.net> and/or <http://www.x-eng.com>

XSteam provides (mostly) accurate steam and water properties from 0 -
1000 bar and from 0 - 2000 °C according to the [IAPWS release IF-97](http://www.iapws.org/relguide/IF97-Rev.pdf). For
accuracy of the functions in different regions see IF-97 Page 4

Also includes thermal conductivity and viscosity, which are not part of
the IF97 release.
* `Thermal Conductivity: (IAPWS 1998) <http://www.iapws.org/relguide/ThCond.pdf>`_

This Python Library is based on the original XSteam Library for Matlab and Excel
from Magnus Holmgren, www.x-eng.com.
We take no responsibilities for any errors in the code or damage thereby!
See README.md for examples

Some effort has been made to include the refined function of more recent releases
and also functions for calculations on heavy water. This includes:
* IAPWS R4
* IAPWS R14
* IAPWS R15-11

Notes
*****
Density (rho)
=============
Density is calculated as 1/v. See section 1.5 Volume

Viscosity
=========
Viscosity is not part of IAPWS Steam IF97. Equations from "Revised Release
on the IAPWS Formulation 1985 for the Viscosity of Ordinary Water
Substance", 2003 are used. Viscosity in the mixed region (4) is interpolated
according to the density. This is not true since it will be two phases.

Thermal conductivity
====================
Revised release on the IAPS Formulation 1985 for the Thermal Conductivity
of ordinary water substance (IAPWS 1998)

Nomenclature
============
All Functions follow the same naming schema: First the wanted property,
then a underscore `_`, then the wanted input properties Example:
`t_ph` is temperature as a function of pressure and enthalpy. For a list
of valid functions se bellow:

========   ============================================================
Property   Description
========   ============================================================
t          Temperature (°C or °F)
p          Pressure (bar or psi)
h          Enthalpy (kJ/kg or btu/lb)
v          Specific volume (m3/kg or ft\^3/lb)|
rho        Density (kg/m3 or lb/ft\^3)
s          Specific entropy (kJ/(kg °C) or btu/(lb °F))
u          Specific internal energy (kJ/kg or btu/lb)
Cp         Specific isobaric heat capacity (kJ/(kg °C) or btu/(lb °F))
Cv         Specific isochoric heat capacity (kJ/(kg °C) or btu/(lb °F))
w          Speed of sound (m/s or ft/s)
my         Viscosity (N s/m\^2 or lbm/ft/hr)
tc         Thermal Conductivity (W/(m °C) or btu/(h ft °F))
st         Surface Tension (N/m or lb/ft)
x          Vapor fraction
vx         Vapor Volume Fraction
========   ============================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   pyXSteamAvailibleFunctions
   pyXSteam
   pyXSteamDemo

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
