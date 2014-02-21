#!/usr/bin/python
# -*- coding: UTF-8 -*-
'''
Created on 03.02.2014

@author: max
'''
import unittest
from pyXSteam import XSteam
import test_pyXSteam

def suite():
    # import unittest
    # import doctest
    suite = unittest.TestSuite()
    # suite.addTests(doctest.DocTestSuite(helloworld))
    suite.addTests(test_pyXSteam.suite())
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity = 2).run(suite())
