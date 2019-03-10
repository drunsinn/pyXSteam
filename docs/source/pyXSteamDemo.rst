Tutorial and Demos
##################

Calculate single values
+++++++++++++++++++++++
This is a simple example::

    >>> from pyXSteam.XSteam import XSteam
    >>> steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
    >>> steamTable.hL_p(220.0)
    2021.909286172027
    >>> steamTable.tcV_p(1)
    0.02475366759235046

Generate Diagrams
+++++++++++++++++
Diagrams based on the calculated values can easily be created using numpy and matplotlib.

Example: To draw a T(p) diagramm showing the saturation curve::

    from pyXSteam.XSteam import XSteam
    import matplotlib.pyplot as pyplot
    import numpy as np
    steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
    p = np.arange(-100.0, 250.0, 1.0)
    ntsat_p = np.frompyfunc(steamTable.tsat_p, 1, 1)
    tsat = ntsat_p(p)
    line1, = pyplot.plot(tsat, p)
    pyplot.xlabel("t")
    pyplot.ylabel("p")
    pyplot.setp(line1, linewidth=1, color='b')
    pyplot.show()

For more demos, see pyXSteamDemo.py

.. literalinclude:: ../../bin/pyXSteamDemo.py
    :language: python
