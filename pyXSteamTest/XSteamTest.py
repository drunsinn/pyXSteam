# -*- coding: UTF-8 -*-
import unittest
import XSteamTest_Regions
import XSteamTest_MKS
# import XSteamTest_FLS

def suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTest(loader.loadTestsFromTestCase(XSteamTest_Regions.Region1Tester))
    suite.addTest(loader.loadTestsFromTestCase(XSteamTest_Regions.Region2Tester))
    suite.addTest(loader.loadTestsFromTestCase(XSteamTest_Regions.Region3Tester))
    suite.addTest(loader.loadTestsFromTestCase(XSteamTest_Regions.Region4Tester))
    suite.addTest(loader.loadTestsFromTestCase(XSteamTest_Regions.Region5Tester))
    suite.addTest(loader.loadTestsFromTestCase(XSteamTest_MKS.MKS_FunctionTester))
    # suite.addTest(loader.loadTestsFromTestCase(XSteamTest_FLS.FLS_FunctionTester))
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity = 2).run(suite())
