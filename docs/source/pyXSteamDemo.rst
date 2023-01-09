Tutorial and Demos
##################

Usage
*****

Simple Example::

    from pyXSteam.XSteam import XSteam
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    print steam_table.hL_p(220.0)

By using the unitSystem Parameter, you can tell XSteam witch Unit System you are using.::

    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS) # m/kg/sec/°C/bar/W
    steam_table = XSteam(XSteam.UNIT_SYSTEM_FLS) # ft/lb/sec/°F/psi/btu
    steam_table = XSteam(XSteam.UNIT_SYSTEM_BARE) # m/kg/sec/K/MPa/W

To enable logging, add the following lines to your code::

    import logging
    logger = logging.getLogger('pyXSteam')
    logger.setLevel(logging.DEBUG) sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(sh)

Calculate single values
=======================
This is a simple example::

    >>> from pyXSteam.XSteam import XSteam
    >>> steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    >>> steam_table.hL_p(220.0)
    2021.909286172027
    >>> steam_table.tcV_p(1)
    0.02475366759235046

Usage with numpy
================
By converting one to the functions to a NumPy universal function it is easy to use pyXSteam on NumPy arrays::

    >>> npv_pt = np.frompyfunc(steam_table.v_pt, 2, 1)

* `NumPy documentation for frombyfunc <https://numpy.org/doc/stable/reference/generated/numpy.frompyfunc.html>`_
* `pyXSteamDemo.py <https://github.com/drunsinn/pyXSteam/blob/master/bin/pyXSteamDemo.py>`_

Usage with pandas
=================
Using pandas together with pyXSteam is very easy as well. Using
the `apply <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.apply.html>`_ function and a
`lambda <https://docs.python.org/3/reference/expressions.html#lambda>`_ you can add a new column based on the
values of other columns::

    >>> df['h'] = df.apply(lambda x: steamTable.h_pt(x.P, x.T), axis=1)

* `Pandas documentation on DataFrame.apply <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.apply.html>`_
* `Python documentation on lambdas <https://docs.python.org/3/reference/expressions.html#lambda>`_

Generate Diagrams
=================
Diagrams based on the calculated values can easily be created using numpy and matplotlib.

Example: To draw a T(p) diagram showing the saturation curve::

    from pyXSteam.XSteam import XSteam
    import matplotlib.pyplot as pyplot
    import numpy as np
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    p = np.arange(-100.0, 250.0, 1.0)
    ntsat_p = np.frompyfunc(steam_table.tsat_p, 1, 1)
    tsat = ntsat_p(p)
    line1, = pyplot.plot(tsat, p)
    pyplot.xlabel("t")
    pyplot.ylabel("p")
    pyplot.setp(line1, linewidth=1, color='b')
    pyplot.show()

For more demos, see `pyXSteamDemo.py <https://github.com/drunsinn/pyXSteam/blob/master/bin/pyXSteamDemo.py>`_

Heavy Water functions
=====================
The functions to calculate values for heavy water are available through the
class XSteamHW

    >>> from pyXSteam.XSteamHW import XSteamHW
    >>> steam_table = XSteamHW(XSteam.UNIT_SYSTEM_MKS)
    >>> steam_table.my_rhoT(1.2, 300.0)

Content of the demo files
*************************

Main demo file pyXSteamDemo.py
==============================
.. literalinclude:: ../../bin/pyXSteamDemo.py
    :language: python

Example on how to calculate the values for a rankine cycle
==========================================================
Matlab example from converted example from Stu Blair converted to python
Original can be found at `his github page <https://github.com/stu314159/xsteam/blob/master/Examples/SimpleRankineCycle.m>`_

.. literalinclude:: ../../bin/pyXSteamRankineDemo.py
    :language: python
