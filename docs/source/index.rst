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

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   pyXSteam
   pyXSteamDemo

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
