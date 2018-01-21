#!/usr/bin/python
# -*- coding: UTF-8 -*-

import unittest

from pyXSteam.XSteam import XSteam
from pyXSteamTest import XSteamTest

def suite():
    # import unittest
    # import doctest
    suite = unittest.TestSuite()
    suite.addTests(XSteamTest.suite())
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity = 1).run(suite())
