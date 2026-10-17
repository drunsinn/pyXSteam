Tutorial and demos
##################

Usage
*****

The main classes and unit-system constants are exported from the package:

.. code-block:: python

    from pyXSteam import XSteam

    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    print(steam_table.hL_p(220.0))

The unit system controls the input and output units:

* ``UNIT_SYSTEM_MKS`` uses pressure in bar and temperature in degrees Celsius.
  Enthalpy is in kJ/kg, density in kg/m³, and thermal conductivity in W/(m K).
* ``UNIT_SYSTEM_FLS`` uses pressure in psi and temperature in degrees Fahrenheit.
* ``UNIT_SYSTEM_BARE`` is currently an alias of ``UNIT_SYSTEM_MKS`` in this
  version; it does **not** select MPa and kelvin. Use ``UNIT_SYSTEM_MKS`` for
  the documented MKS conversions.

Calculate single values
=======================

Methods are named after the property and input properties. For example,
``hL_p`` returns saturated-liquid enthalpy as a function of pressure, while
``t_ph`` returns temperature from pressure and enthalpy:

.. code-block:: pycon

    >>> from pyXSteam import XSteam
    >>> steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    >>> steam_table.hL_p(220.0)
    2021.909286172027
    >>> steam_table.t_ph(1.0, 500.0)
    99.60591861133764

Methods that cannot calculate a result for their input commonly return ``NaN``.
For example, ``tsat_p`` expects a pressure between the triple-point and critical
pressures. Check results with ``math.isnan`` or ``numpy.isnan`` when inputs may
be outside a method's valid range.

Usage with NumPy
================

NumPy is an optional dependency. ``numpy.frompyfunc`` can apply a pyXSteam
method element-wise to arrays:

.. code-block:: python

    import numpy as np
    from pyXSteam import XSteam

    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    pressures_bar = np.arange(1.0, 221.0, 1.0)
    tsat = np.frompyfunc(steam_table.tsat_p, 1, 1)(pressures_bar)
    tsat = np.asarray(tsat, dtype=float)

See the `NumPy frompyfunc documentation
<https://numpy.org/doc/stable/reference/generated/numpy.frompyfunc.html>`_ and
the `pyXSteam demo source
<https://github.com/drunsinn/pyXSteam/blob/master/pyXSteam/scripts/pyXSteamDemo.py>`_.

Usage with pandas
=================

If pandas is installed, ``DataFrame.apply`` can calculate a result from columns
of input data. This example expects pressure in bar and temperature in degrees
Celsius:

.. code-block:: python

    import pandas as pd
    from pyXSteam import XSteam

    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    states = pd.DataFrame({"p": [1.0, 5.0], "t": [100.0, 200.0]})
    states["h"] = states.apply(
        lambda row: steam_table.h_pt(row["p"], row["t"]),
        axis=1,
    )

See the `pandas DataFrame.apply documentation
<https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.apply.html>`_.

Generate diagrams
=================

The following example plots the saturation temperature against pressure. The
selected pressures are within the saturation curve's valid range:

.. code-block:: python

    import matplotlib.pyplot as plt
    import numpy as np
    from pyXSteam import XSteam

    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    pressure_bar = np.arange(1.0, 221.0, 1.0)
    tsat = np.frompyfunc(steam_table.tsat_p, 1, 1)(pressure_bar)
    tsat = np.asarray(tsat, dtype=float)

    plt.plot(pressure_bar, tsat)
    plt.xlabel("Pressure [bar]")
    plt.ylabel("Saturation temperature [°C]")
    plt.grid(True)
    plt.show()

For more examples, see the `pyXSteam demo source
<https://github.com/drunsinn/pyXSteam/blob/master/pyXSteam/scripts/pyXSteamDemo.py>`_.
Install its optional dependencies with ``pip install "pyXSteam[demo]"``.

More diagrams
-------------

For additional examples of diagrams generated with pyXSteam, see the
repositories by Richard J Smith:

* `steam-properties <https://github.com/Richard66NZ/steam-properties>`_
* `steam-cycles-Ideal <https://github.com/Richard66NZ/steam-cycles-Ideal>`_
* `steam-cycles-non-ideal <https://github.com/Richard66NZ/steam-cycles-non-ideal>`_
* `steam-cycles-regenerative <https://github.com/Richard66NZ/steam-cycles-regenerative>`_

Heavy water
===========

The heavy-water API is exposed as ``XSteam_HW``:

.. code-block:: python

    from pyXSteam import XSteam_HW

    heavy_water = XSteam_HW(XSteam_HW.UNIT_SYSTEM_MKS)

Heavy-water support is still in progress and should not be relied on for
production calculations. In particular, the viscosity and thermal-conductivity
methods use the superseded IAPWS R4-84 formulation.

Logging
-------

pyXSteam uses Python's standard logging hierarchy and does not configure the
root logger. Configure logging once in the application. Routine out-of-range
diagnostics are logged at ``DEBUG`` because many methods return ``NaN`` for
unsupported inputs:

.. code-block:: python

    import logging

    logging.basicConfig(level=logging.WARNING)
    logging.getLogger("pyXSteam").setLevel(logging.DEBUG)

Demo scripts
************

The interactive demo commands are installed with the optional dependencies:

.. code-block:: console

    $ pip install "pyXSteam[demo]"
    $ pyXSteamDemo
    $ pyXSteamRankineDemo

Main demo file ``pyXSteamDemo.py``
==================================

.. literalinclude:: ../../pyXSteam/scripts/pyXSteamDemo.py
    :language: python

Rankine cycle demo
==================

This example implements Rankine cycles based on examples by Stu Blair. The
original Matlab example is available on `GitHub
<https://github.com/stu314159/xsteam/blob/master/Examples/SimpleRankineCycle.m>`_.

.. literalinclude:: ../../pyXSteam/scripts/pyXSteamRankineDemo.py
    :language: python
