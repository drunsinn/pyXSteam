#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test runner for pyXSteam
"""
import unittest
import XSteamTest_Regions
import XSteamTest_MKS
import XSteamTest_UnitConverter


def suite():
    """Run tests described in IAPWS release IF-97."""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTest(loader.loadTestsFromTestCase(XSteamTest_Regions.Region1Tester))
    suite.addTest(loader.loadTestsFromTestCase(XSteamTest_Regions.Region2Tester))
    suite.addTest(loader.loadTestsFromTestCase(XSteamTest_Regions.Region3Tester))
    suite.addTest(loader.loadTestsFromTestCase(XSteamTest_Regions.Region4Tester))
    suite.addTest(loader.loadTestsFromTestCase(XSteamTest_Regions.Region5Tester))
    suite.addTest(loader.loadTestsFromTestCase(XSteamTest_MKS.MKS_FunctionTester))
    suite.addTest(loader.loadTestsFromTestCase(XSteamTest_UnitConverter.UnitConverter_MKS_FunctionTester))
    suite.addTest(loader.loadTestsFromTestCase(XSteamTest_UnitConverter.UnitConverter_FLS_FunctionTester))
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
