# -*- coding: UTF-8 -*-
import unittest
import TestXSteam_Regions
import TestXSteam_MKS

def suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTest(loader.loadTestsFromTestCase(TestXSteam_Regions.Region1Tester))
    suite.addTest(loader.loadTestsFromTestCase(TestXSteam_Regions.Region2Tester))
    suite.addTest(loader.loadTestsFromTestCase(TestXSteam_Regions.Region3Tester))
    suite.addTest(loader.loadTestsFromTestCase(TestXSteam_Regions.Region4Tester))
    suite.addTest(loader.loadTestsFromTestCase(TestXSteam_Regions.Region5Tester))
    suite.addTest(loader.loadTestsFromTestCase(TestXSteam_MKS.MKS_FunctionTester))
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity = 2).run(suite())

# if __name__ == "__main__":
#    # import sys;sys.argv = ['', 'Test.testName']
#    unittest.main()