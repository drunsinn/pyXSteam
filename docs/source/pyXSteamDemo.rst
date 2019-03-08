pyXSteam Demos
##############

Calculate Values
++++++++++++++++
This is a simple example::

    from pyXSteam.XSteam import XSteam
    steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
    steamTable.hL_p(220.0)
    >>>
    steamTable.tcV_p(1.0)
    >>>

Generate Diagrams
+++++++++++++++++

.. literalinclude:: ../../bin/pyXSteamDemo.py
    :linenos:
    :language: python
    :lines: 77-86
