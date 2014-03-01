# -*- coding: utf-8 -*-
import unittest
import numpy
from pyXSteam.XSteam import XSteam
from pyXSteam.Regions import Region1, Region2, Region3, Region4, Region5

class Region1Tester(unittest.TestCase):

    def setUp(self):
        self.maxError = 1E-8

    def tearDown(self):
        pass

    def test_pT_functions(self):
        # '''Tests to verify all functions with the Parameters p and T of Region 1 by comparing the Results to IF-97 Page 9 Table 5'''
        # %* 7.1 Verifiy region 1
        # %IF-97 Table 5, Page 9
        p = [3.0, 80.0, 3.0]
        T = [300.0, 300.0, 500.0]
        IF97 = [[0.00100215168, 0.000971180894, 0.001202418],
                [115.331273, 184.142828, 975.542239],
                [112.324818, 106.448356, 971.934985],
                [0.392294792, 0.368563852, 2.58041912],
                [4.17301218, 4.01008987, 4.65580682],
                [1507.73921, 1634.69054, 1240.71337]]
        # IF97 = [[0.100215168e-2, 0.971180894e-3, 0.120241800e-2],
        #        [0.115331273e3, 0.184142828e3, 0.975542239e3],
        #        [0.112324818e3, 0.106448356e3, 0.971934985e3],
        #        [0.392294792e0, 0.368563852e0, 0.258041912e1],
        #        [0.417301218e1, 0.401008987e1, 0.465580682e1],
        #        [0.150773921e4, 0.163469054e4, 0.124071337e4]]
        R1 = numpy.zeros((6, 3))

        for i in range(0, 3):
            R1[0][i] = Region1.v1_pT(p[i], T[i])
            R1[1][i] = Region1.h1_pT(p[i], T[i])
            R1[2][i] = Region1.u1_pT(p[i], T[i])
            R1[3][i] = Region1.s1_pT(p[i], T[i])
            R1[4][i] = Region1.Cp1_pT(p[i], T[i])
            R1[5][i] = Region1.w1_pT(p[i], T[i])

        print IF97
        print R1
        Region1_error = numpy.sum(numpy.absolute((R1 - IF97) / IF97))
        self.assertLess(Region1_error, self.maxError, 'Test of p,T Functions for Region 1 failed. Error was %(error)e allowed: %(max)e' % {'error' : Region1_error, 'max' : self.maxError})

    def test_ph_function(self):
        # '''Tests to verify all functions with the Parameters p and h of Region 1 by comparing the Results to IF-97 Page 11 Table 7'''
        # % IF - 97 Table 7, Page 11
        p = [3.0, 80.0, 80.0]
        h = [500.0, 500.0, 1500.0]
        IF97 = [391.798509, 378.108626, 611.041229]
        R1 = numpy.zeros(3)
        for i in range(0, 3):
            R1[i] = Region1.T1_ph(p[i], h[i])

        T1_ph_error = numpy.sum(numpy.absolute((R1 - IF97) / IF97))
        self.assertLess(T1_ph_error, self.maxError, 'Test of T(p,h) Function for Region 1 failed. Error was %(error)e allowed: %(max)e' % {'error' : T1_ph_error, 'max' : self.maxError})
        # if T1_ph_error > 1E-8:
        #    if self.verbose:
        #        print'\tRegion 1 Test of ph functions failed. Error is:' , T1_ph_error
        #    return False
        # else:
        #    return True

    def test_ps_function(self):
        # '''Tests to verify all functions with the Parameters p and s of Region 1 by comparing the Results to IF-97 Page 12 Table 9'''
        # % IF - 97 Table 9, Page 12
        p = [3.0, 80.0, 80.0]
        s = [0.5, 0.5, 3.0];
        IF97 = [307.842258, 309.979785, 565.899909]
        R1 = numpy.zeros(3)
        for i in range(0, 3):
            R1[i] = Region1.T1_ps(p[i], s[i])

        T1_ps_error = numpy.sum(numpy.absolute((R1 - IF97) / IF97))
        self.assertLess(T1_ps_error, self.maxError, 'Test of T(p,s) Function for Region 1 failed. Error was %(error)e allowed: %(max)e' % {'error' : T1_ps_error, 'max' : self.maxError})


    def test_hs_function(self):
        # '''Tests to verify all functions with the Parameters h and s of Region 1 by comparing the Results to IF-97 Page 6 Table 3'''
        # % Supplementary Release on Backward Equations
        # % for Pressure as a Function of Enthalpy and Entropy p(h, s)
        # % Table 3, Page 6
        h = [0.001, 90.0, 1500.0]
        s = [0.0, 0.0, 3.4]
        IF97 = [0.0009800980612, 91.929547272, 58.68294423]
        R1 = numpy.zeros(3)
        for i in range(0, 3):
            R1[i] = Region1.p1_hs(h[i], s[i])

        p1_hs_error = numpy.sum(numpy.absolute((R1 - IF97) / IF97))
        self.assertLess(p1_hs_error, self.maxError, 'Test of p(h,s) Function for Region 1 failed. Error was %(error)e allowed: %(max)e' % {'error' : p1_hs_error, 'max' : self.maxError})
        # if p1_hs_error > 1E-8:
        #    if self.verbose:
        #        print '\tRegion 1 Test of hs functions failed. Error is:' , p1_hs_error
        #    return False
        # else:
        #    return True
class Region2Tester(unittest.TestCase):

    def setUp(self):
        self.maxError = 1E-8

    def tearDown(self):
        pass

    def test_pT_function(self):
        # '''Tests to verify all functions with the Parameters p and T of Region 2 by comparing the Results to IF-97 Page 17 Table 15'''
        # %* 7.2 Verifiy region 2
        # % IF-97 Table 15, Page 17
        p = [0.0035 , 0.0035 , 30.0]
        T = [300.0, 700.0, 700.0]
        # Fun = {'v2_pT', 'h2_pT', 'u2_pT', 's2_pT', 'Cp2_pT', 'w2_pT'};
        IF97 = [[39.4913866, 92.3015898, 0.00542946619],
                [2549.91145, 3335.68375, 2631.49474],
                [2411.6916, 3012.62819, 2468.61076],
                [8.52238967, 10.1749996, 5.17540298],
                [1.91300162, 2.08141274, 10.3505092],
                [427.920172, 644.289068, 480.386523]]
        R2 = numpy.zeros((6, 3))

        for i in range(0, 3):
            R2[0][i] = Region2.v2_pT(p[i], T[i])
            R2[1][i] = Region2.h2_pT(p[i], T[i])
            R2[2][i] = Region2.u2_pT(p[i], T[i])
            R2[3][i] = Region2.s2_pT(p[i], T[i])
            R2[4][i] = Region2.Cp2_pT(p[i], T[i])
            R2[5][i] = Region2.w2_pT(p[i], T[i])

        Region2_error = numpy.sum(numpy.absolute((R2 - IF97) / IF97))
        self.assertLess(Region2_error, self.maxError, 'Test of p,T Functions for Region 2 failed. Error was %(error)e allowed: %(max)e' % {'error' : Region2_error, 'max' : self.maxError})
        # if Region2_error > 1E-8:
        #    if self.verbose:
        #        print '\tRegion 2 Test of pT functions failed. Error is:', Region2_error
        #    return False
        # else:
        #    return True

    def test_ph_function(self):
        # '''Tests to verify all functions with the Parameters p and h of Region 2 by comparing the Results to IF-97 Page 25 Table 24'''
        # % IF - 97 Table 24, Page 25
        p = [0.01 / 10, 30 / 10, 30 / 10, 50 / 10, 50 / 10, 250 / 10, 400 / 10, 600 / 10, 600 / 10]
        h = [3000.0, 3000.0, 4000.0, 3500.0, 4000.0, 3500.0, 2700.0, 2700.0, 3200.0]
        IF97 = [534.433241, 575.37337, 1010.77577, 801.299102, 1015.31583, 875.279054, 743.056411, 791.137067, 882.75686]
        R2 = numpy.zeros(9)
        for i in range(0, 9):
            R2[i] = Region2.T2_ph(p[i], h[i])
        T2_ph_error = numpy.sum(numpy.absolute((R2 - IF97) / IF97))
        self.assertLess(T2_ph_error, self.maxError, 'Test of ph Function for Region 2 failed. Error was %(error)e allowed: %(max)e' % {'error' : T2_ph_error, 'max' : self.maxError})
        # if T2_ph_error > 1E-8:
        #    if self.verbose:
        #        print '\tRegion 2 Test of ph functions failed. Error is:', T2_ph_error
        #    return False
        # else:
        #    return True

    def test_ps_function(self):
        # '''Tests to verify all functions with the Parameters p and s of Region 2 by comparing the Results to IF-97 Page 29 Table 29'''
        # % IF - 97 Table 29, Page 29
        p = [0.1, 0.1, 2.5, 8.0, 8.0, 90.0, 20.0, 80.0, 80.0]
        s = [7.5, 8.0, 8.0, 6.0, 7.5, 6.0, 5.75, 5.25, 5.75]
        IF97 = [399.517097, 514.127081, 1039.84917, 600.48404, 1064.95556, 1038.01126, 697.992849, 854.011484, 949.017998]
        R2 = numpy.zeros(9)
        for i in range(0, 9):
            R2[i] = Region2.T2_ps(p[i], s[i])
        T2_ps_error = numpy.sum(numpy.absolute((R2 - IF97) / IF97))
        self.assertLess(T2_ps_error, self.maxError, 'Test of ps Function for Region 2 failed. Error was %(error)e allowed: %(max)e' % {'error' : T2_ps_error, 'max' : self.maxError})
        # if T2_ps_error > 1E-8:
        #    if self.verbose:
        #        print '\tRegion 2 Test of ps functions failed. Error is:', T2_ps_error
        #    return False
        # else:
        #    return True

    def test_hs_function(self):
        # '''Tests to verify all functions with the Parameters h and s of Region 2 by comparing the Results to IF-97 Page 6 Table 3'''
        # % Supplementary Release on Backward Equations for Pressure as a Function of Enthalpy and Entropy p(h, s)
        # % Table 3, Page 6
        h = [2800.0, 2800.0, 4100.0, 2800.0, 3600.0, 3600.0, 2800.0, 2800.0, 3400.0]
        s = [6.5, 9.5, 9.5, 6, 6, 7, 5.1, 5.8, 5.8]
        IF97 = [1.371012767, 0.001879743844, 0.1024788997, 4.793911442, 83.95519209, 7.527161441, 94.3920206, 8.414574124, 83.76903879]
        R2 = numpy.zeros(9)
        for i in range(0, 9):
            R2[i] = Region2.p2_hs(h[i], s[i])
        p2_hs_error = numpy.sum(numpy.absolute((R2 - IF97) / IF97))
        self.assertLess(p2_hs_error, self.maxError, 'Test of hs Function for Region 2 failed. Error was %(error)e allowed: %(max)e' % {'error' : p2_hs_error, 'max' : self.maxError})
        # if p2_hs_error > 1E-8:
        #    if self.verbose:
        #        print '\tRegion 2 Test of hs functions failed. Error is:', p2_hs_error
        #    return False
        # else:
        #    return True

class Region3Tester(unittest.TestCase):

    def setUp(self):
        self.maxError = 1E-8

    def tearDown(self):
        pass

    def test_rhoT_function(self):
        # '''Tests to verify all functions with the Parameters rho and T of Region 3 by comparing the Results to IF-97 Page 32 Table 33'''
        # % IF-97 Table 33, Page 32
        T = [650.0, 650.0, 750.0]
        rho = [500.0, 200.0, 500.0]
        # Fun={'p3_rhoT','h3_rhoT','u3_rhoT','s3_rhoT','Cp3_rhoT','w3_rhoT'};
        IF97 = [[25.5837018, 22.2930643, 78.3095639],
                [1863.43019, 2375.12401, 2258.68845],
                [1812.26279, 2263.65868, 2102.06932],
                [4.05427273, 4.85438792, 4.46971906],
                [13.8935717, 44.6579342, 6.34165359],
                [502.005554, 383.444594, 760.696041]]

        R3 = numpy.zeros((6, 3))
        for i in range(0, 3):
            R3[0][i] = Region3.p3_rhoT(rho[i], T[i])
            R3[1][i] = Region3.h3_rhoT(rho[i], T[i])
            R3[2][i] = Region3.u3_rhoT(rho[i], T[i])
            R3[3][i] = Region3.s3_rhoT(rho[i], T[i])
            R3[4][i] = Region3.Cp3_rhoT(rho[i], T[i])
            R3[5][i] = Region3.w3_rhoT(rho[i], T[i])

        # print R3
        # print IF97
        Region3_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(Region3_error, self.maxError, 'Test of rhoT Function for Region 3 failed. Error was %(error)e allowed: %(max)e' % {'error' : Region3_error, 'max' : self.maxError})
        # if Region3_error > 1E-8:
        #    if self.verbose:
        #        print '\tRegion 3 Test of rhoT functions failed. Error is:', Region3_error
        #    return False
        # else:
        #    return True

    def test_T_ph_function(self):
        # '''Tests to verify all temperature functions with the Parameters p and h of Region 3'''
        # % T3_ph
        p = [20.0, 50.0, 100.0, 20.0, 50.0, 100.0]
        h = [1700.0, 2000.0, 2100.0, 2500.0, 2400.0, 2700.0]
        IF97 = [629.3083892, 690.5718338, 733.6163014, 641.8418053, 735.1848618, 842.0460876]
        R3 = numpy.zeros(6)
        for i in range(0, 6):
            R3[i] = Region3.T3_ph(p[i], h[i]);
        T3_ph_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(T3_ph_error, self.maxError, 'Test of T(p,h) Function for Region 3 failed. Error was %(error)e allowed: %(max)e' % {'error' : T3_ph_error, 'max' : self.maxError})
        # if T3_ph_error > 1E-8:
        #    if self.verbose:
        #        print '\tRegion 3 Test of T(p,h) functions failed. Error is:', T3_ph_error
        #    return False
        # else:
        #    return True

    def test_v_ph_function(self):
        # '''Tests to verify all v functions with the Parameters p and h of Region 3'''
        # % v3_ph
        p = [20.0, 50.0, 100.0, 20.0, 50.0, 100.0]
        h = [1700.0, 2000.0, 2100.0, 2500.0, 2400.0, 2700.0]
        IF97 = [0.001749903962, 0.001908139035, 0.001676229776, 0.006670547043, 0.0028012445, 0.002404234998]
        R3 = numpy.zeros(6)
        for i in range(0, 6):
            R3[i] = Region3.v3_ph(p[i], h[i])
        v3_ph_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(v3_ph_error, 1E-7, 'Test of v(p,h) Function for Region 3 failed. Error was %(error)e allowed: %(max)e' % {'error' : v3_ph_error, 'max' : 1E-7})
        # if v3_ph_error > 1E-7:
        #    if self.verbose:
        #        print '\tRegion 3 Test of v(p,h) functions failed. Error is:', v3_ph_error
        #    return False
        # else:
        #    return True

    def test_T_ps_function(self):
        # '''Tests to verify all T functions with the Parameters p and s of Region 3'''
        # % T3_ps
        p = [20.0, 50.0, 100.0, 20.0, 50.0, 100.0]
        s = [3.7, 3.5, 4, 5, 4.5, 5.0]
        IF97 = [620.8841563, 618.1549029, 705.6880237, 640.1176443, 716.3687517, 847.4332825]
        R3 = numpy.zeros(6)
        for i in range(0, 6):
            R3[i] = Region3.T3_ps(p[i], s[i])
        T3_ps_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(T3_ps_error, self.maxError, 'Test of T(p,s) Function for Region 3 failed. Error was %(error)e allowed: %(max)e' % {'error' : T3_ps_error, 'max' : self.maxError})
        # if T3_ps_error > 1E-8:
        #    if self.verbose:
        #        print '\tRegion 3 Test of ps functions failed. Error is:', T3_ps_error
        #    return False
        # else:
        #    return True

    def test_v_ps_function(self):
        # '''Tests to verify all v functions with the Parameters p and s of Region 3'''
        # % v3_ps
        p = [20.0, 50.0, 100.0, 20.0, 50.0, 100.0]
        s = [3.7, 3.5, 4.0, 5.0, 4.5, 5.0]
        IF97 = [0.001639890984, 0.001423030205, 0.001555893131, 0.006262101987, 0.002332634294, 0.002449610757]
        R3 = numpy.zeros(6)
        for i in range(0, 6):
            R3[i] = Region3.v3_ps(p[i], s[i])
        v3_ps_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(v3_ps_error, self.maxError, 'Test of v(p,s) Function for Region 3 failed. Error was %(error)e allowed: %(max)e' % {'error' : v3_ps_error, 'max' : self.maxError})
        # if v3_ps_error > 1E-8:
        #    if self.verbose:
        #        print '\tRegion 3 Test of v(p,s) functions failed. Error is:', v3_ps_error
        #    return False
        # else:
        #    return True

    def test_hs_function(self):
        # '''Tests to verify all functions with the Parameters h and s of Region 3'''
        # % p3_hs
        h = [1700.0, 2000.0, 2100.0, 2500.0, 2400.0, 2700.0]
        s = [3.8, 4.2, 4.3, 5.1, 4.7, 5.0]
        IF97 = [25.55703246, 45.40873468, 60.7812334, 17.20612413, 63.63924887, 88.39043281]
        R3 = numpy.zeros(6)
        for i in range(0, 6):
            R3[i] = Region3.p3_hs(h[i], s[i])
        p3_hs_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(p3_hs_error, self.maxError, 'Test of p(h,s) Function for Region 3 failed. Error was %(error)e allowed: %(max)e' % {'error' : p3_hs_error, 'max' : self.maxError})
        # if p3_hs_error > 1E-8:
        #    if self.verbose:
        #        print '\tRegion 3 Test of hs functions failed. Error is:', p3_hs_error
        #    return False
        # else:
        #    return True

    def test_pT_function(self):
        # '''Tests to verify all functions with the Parameters p and T of Region 3'''
        # % h3_pT (Iteration)
        p = [25.583702, 22.293064, 78.309564]
        T = [650.0, 650.0, 750.0]
        IF97 = [1863.271389, 2375.696155, 2258.626582]
        R3 = numpy.zeros(3)
        for i in range(0, 3):
            R3[i] = Region3.h3_pT(p[i], T[i])
        h3_pT_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(h3_pT_error, 1E-6, 'Test of h(p,T) Function for Region 3 failed. Error was %(error)e allowed: %(max)e' % {'error' : h3_pT_error, 'max' :  1E-6})
        # if h3_pT_error > 1E-6:  # % Decimals in IF97
        #    if self.verbose:
        #        print '\tRegion 3 Test of pT functions failed. Error is:', h3_pT_error
        #    return False
        # else:
        #    return True

class Region4Tester(unittest.TestCase):

    def setUp(self):
        self.maxError = 1E-7

    def tearDown(self):
        pass

    def test_T_function(self):
        # '''Tests to verify all functions with the Parameter T of Region 4 by comparing the Results to IF-97 Page 34 Table 35'''
        # %Saturation pressure, If97, Table 35, Page 34
        T = [300.0, 500.0, 600.0]
        IF97 = [0.00353658941, 2.63889776, 12.3443146]
        R4 = numpy.zeros(3)
        for i in range(0, 3):
            R4[i] = Region4.p4_T(T[i])
        p4_t_error = numpy.sum(numpy.absolute((R4 - IF97) / IF97))
        self.assertLess(p4_t_error, self.maxError, 'Test of p(T) Function for Region 4 failed. Error was %(error)e allowed: %(max)e' % {'error' : p4_t_error, 'max' : self.maxError})
        # if p4_t_error > 1E-7:
        #    if self.verbose:
        #        print '\tRegion 4 Test of T functions failed. Error is:', p4_t_error
        #    return False
        # else:
        #    return True

    def test_p_functions(self):
        # '''Tests to verify all functions with the Parameters p of Region 4'''
        p = [0.1, 1.0, 10.0]
        IF97 = [372.755919, 453.035632, 584.149488]
        R4 = numpy.zeros(3)
        for i in range(0, 3):
            R4[i] = Region4.T4_p(p[i])
        T4_p_error = numpy.sum(numpy.absolute((R4 - IF97) / IF97))
        self.assertLess(T4_p_error, self.maxError, 'Test of T(p) Function for Region 4 failed. Error was %(error)e allowed: %(max)e' % {'error' : T4_p_error, 'max' : self.maxError})
        # if T4_p_error > 1E-7:
        #    if self.verbose:
        #        print '\tRegion 4 Test of p functions failed. Error is:', T4_p_error
        #    return False
        # else:
        #    return True

    def test_s_functions(self):
        # '''Tests to verify all functions with the Parameters s of Region 4'''
        s = [1.0, 2.0, 3.0, 3.8, 4.0, 4.2, 7.0, 8.0, 9.0, 5.5, 5.0, 4.5]
        IF97 = [308.5509647, 700.6304472, 1198.359754, 1685.025565, 1816.891476, 1949.352563, 2723.729985, 2599.04721, 2511.861477, 2687.69385, 2451.623609, 2144.360448]
        R4 = numpy.zeros(12)
        for i in range(0, 12):
            R4[i] = Region4.h4_s(s[i])
        h4_s_error = numpy.sum(numpy.absolute((R4 - IF97) / IF97))
        self.assertLess(h4_s_error, self.maxError, 'Test of h(s) Function for Region 4 failed. Error was %(error)e allowed: %(max)e' % {'error' : h4_s_error, 'max' : self.maxError})
        # if h4_s_error > 1E-7:
        #    if self.verbose:
        #        print '\tRegion 4 Test of s functions failed. Error is:', h4_s_error
        #    return False
        # else:
        #    return True

class Region5Tester(unittest.TestCase):

    def setUp(self):
        self.maxError = 1E-8

    def tearDown(self):
        pass

    def test_pT_function(self):
        # '''Tests to verify all functions with the Parameters p and T of Region 5 by comparing the Results to IF-97 Page 39 Table 42'''
        # % IF-97 Table 42, Page 39
        T = [1500.0, 1500.0, 2000.0]
        p = [0.5, 8.0, 8.0]
        IF97 = [[1.38455354, 0.0865156616, 0.115743146],
                [5219.76332, 5206.09634, 6583.80291],
                [4527.48654, 4513.97105, 5657.85774],
                [9.65408431, 8.36546724, 9.15671044],
                [2.61610228, 2.64453866, 2.8530675],
                [917.071933, 919.708859, 1054.35806]]
        R5 = numpy.zeros((6, 3))
        for i in range(0, 3):
            R5[0][i] = Region5.v5_pT(p[i], T[i])
            R5[1][i] = Region5.h5_pT(p[i], T[i])
            R5[2][i] = Region5.u5_pT(p[i], T[i])
            R5[3][i] = Region5.s5_pT(p[i], T[i])
            R5[4][i] = Region5.Cp5_pT(p[i], T[i])
            R5[5][i] = Region5.w5_pT(p[i], T[i])

        Region5_error = numpy.sum(numpy.absolute((R5 - IF97) / IF97))
        self.assertLess(Region5_error, self.maxError, 'Test of p,T Function for Region 5 failed. Error was %(error)e allowed: %(max)e' % {'error' : Region5_error, 'max' : self.maxError})
        # if Region5_error > 1E-8:
        #    if self.verbose:
        #        print '\tRegion 5 Test of pT functions failed. Error is:', Region5_error
        #    return False
        # else:
        #    return True

    def test_ph_function(self):
        # '''Tests to verify all functions with the Parameters p and h of Region 5'''
        # %T5_ph (Iteration)
        p = [0.5, 8.0, 8.0]
        h = [5219.76331549428, 5206.09634477373, 6583.80290533381]
        IF97 = [1500.0, 1500.0, 2000.0];
        R5 = numpy.zeros(3)
        for i in range(0, 3):
            R5[i] = Region5.T5_ph(p[i], h[i])
        T5_ph_error = numpy.sum(numpy.absolute((R5 - IF97) / IF97))
        self.assertLess(T5_ph_error, self.maxError, 'Test of T(p,h) Function for Region 5 failed. Error was %(error)e allowed: %(max)e' % {'error' : T5_ph_error, 'max' : self.maxError})
        # if T5_ph_error > 1E-7:  # % Decimals in IF97
        #    if self.verbose:
        #        print '\tRegion 5 Test of ph functions failed. Error is:', T5_ph_error
        #    return False
        # else:
        #    return True

    def test_ps_function(self):
        # '''Tests to verify all functions with the Parameters p and s of Region 5'''
        # %T5_ps (Iteration)
        p = [0.5, 8.0, 8.0]
        s = [9.65408430982588, 8.36546724495503, 9.15671044273249]
        IF97 = [1500.0, 1500.0, 2000.0]
        R5 = numpy.zeros(3)
        for i in range(0, 3):
            R5[i] = Region5.T5_ps(p[i], s[i]);
        T5_ps_error = numpy.sum(numpy.absolute((R5 - IF97) / IF97))
        self.assertLess(T5_ps_error, 1E-4, 'Test of T(p,s) Function for Region 5 failed. Error was %(error)e allowed: %(max)e' % {'error' : T5_ps_error, 'max' : 1E-4})
        # if T5_ps_error > 1E-4:  # % Decimals in IF97
        #    if self.verbose:
        #        print '\tRegion 5 Test of ps functions failed. Error is:', T5_ps_error
        #    return False
        # else:
        #    return True

class CallFunctionTester(unittest.TestCase):

    def setUp(self):
        self.maxError = 1E-8
        self.steamTable = XSteam(mksSystem = True)

    def tearDown(self):
        pass

    def test_tsat_p(self):
        error = self.steamTable.tsat_p(1.0) - 99.60591861
        self.assertLess(error, self.maxError, 'tsat_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_t_ph(self):
        error = self.steamTable.t_ph(1.0, 100.0) - 23.84481908
        self.assertLess(error, self.maxError, 't_ph not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_t_ps(self):
        error = self.steamTable.t_ps(1.0, 1.0) - 73.70859421
        self.assertLess(error, self.maxError, 't_ps not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_t_hs(self):
        error = self.steamTable.t_hs(100.0, 2.0) - 13.84933511
        self.assertLess(error, self.maxError, 't_hs not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_psat_t(self):
        error = self.steamTable.psat_t(100.0) - 1.014179779
        self.assertLess(error, self.maxError, 'psat_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_p_hs(self):
        error = self.steamTable.p_hs(84.0, 0.296) - 2.295498269
        self.assertLess(error, self.maxError, 'p_hs not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_hV_p(self):
        error = self.steamTable.hV_p(1.0) - 2674.949641
        self.assertLess(error, self.maxError, 'hV_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_hL_p(self):
        error = self.steamTable.hL_p(1.0) - 417.4364858
        self.assertLess(error, self.maxError, 'hL_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_hV_t(self):
        error = self.steamTable.hV_t(100.0) - 2675.572029
        self.assertLess(error, self.maxError, 'hV_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_hL_t(self):
        error = self.steamTable.hL_t(100.0) - 419.099155
        self.assertLess(error, self.maxError, 'hL_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_h_pt(self):
        error = self.steamTable.h_pt(1.0, 20.0) - 84.01181117
        self.assertLess(error, self.maxError, 'h_pt not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_h_ps(self):
        error = self.steamTable.h_ps(1.0, 1.0) - 308.6107171
        self.assertLess(error, self.maxError, 'h_ps not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_h_px(self):
        error = self.steamTable.h_px(1.0, 0.5) - 1546.193063
        self.assertLess(error, self.maxError, 'h_px not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_h_prho(self):
        error = self.steamTable.h_prho(1.0, 2.0) - 1082.773391
        self.assertLess(error, self.maxError, 'h_prho not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_h_tx(self):
        error = self.steamTable.h_tx(100.0, 0.5) - 1547.33559211
        self.assertLess(error, self.maxError, 'h_tx not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_vV_p(self):
        error = self.steamTable.vV_p(1.0) - 1.694022523
        self.assertLess(error, self.maxError, 'vV_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_vL_p(self):
        error = self.steamTable.vL_p(1.0) - 0.001043148
        self.assertLess(error, self.maxError, 'vL_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_vV_t(self):
        error = self.steamTable.vV_t(100.0) - 1.671860601
        self.assertLess(error, self.maxError, 'vV_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_vL_t(self):
        error = self.steamTable.vL_t(100.0) - 0.001043455
        self.assertLess(error, self.maxError, 'vL_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_v_pt(self):
        error = self.steamTable.v_pt(1.0, 100.0) - 1.695959407
        self.assertLess(error, self.maxError, 'v_pt not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_v_ph(self):
        error = self.steamTable.v_ph(1.0, 1000.0) - 0.437925658
        self.assertLess(error, self.maxError, 'v_ph not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_v_ps(self):
        error = self.steamTable.v_ps(1.0, 5.0) - 1.03463539
        self.assertLess(error, self.maxError, 'v_ps not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_rhoV_p(self):
        error = self.steamTable.rhoV_p(1.0) - 0.590310924
        self.assertLess(error, self.maxError, 'rhoV_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_rhoL_p(self):
        error = self.steamTable.rhoL_p(1.0) - 958.6368897
        self.assertLess(error, self.maxError, 'rhoL_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_rhoV_t(self):
        error = self.steamTable.rhoV_t(100.0) - 0.598135993
        self.assertLess(error, self.maxError, 'rhoV_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_rhoL_t(self):
        error = self.steamTable.rhoL_t(100.0) - 958.3542773
        self.assertLess(error, self.maxError, 'rhoL_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_rho_pt(self):
        error = self.steamTable.rho_pt(1.0, 100.0) - 0.589636754
        self.assertLess(error, self.maxError, 'rho_pt not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_rho_ph(self):
        error = self.steamTable.rho_ph(1.0, 1000.0) - 2.283492601
        self.assertLess(error, self.maxError, 'rho_ph not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_rho_ps(self):
        error = self.steamTable.rho_ps(1.0, 1.0) - 975.6236788
        self.assertLess(error, self.maxError, 'rho_ps not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_sV_p(self):
        error = self.steamTable.sV_p(0.006117) - 9.155465556
        self.assertLess(error, self.maxError, 'sV_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_sL_p(self):
        error = self.steamTable.sL_p(0.0061171) - 1.8359e-05
        self.assertLess(error, self.maxError, 'sL_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_sV_t(self):
        error = self.steamTable.sV_t(0.0001) - 9.155756716
        self.assertLess(error, self.maxError, 'sV_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_sL_t(self):
        error = self.steamTable.sL_t(100.0) - 1.307014328
        self.assertLess(error, self.maxError, 'sL_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})
    def test_s_pt(self):
        error = self.steamTable.s_pt(1.0, 20.0) - 0.296482921
        self.assertLess(error, self.maxError, 's_pt not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_s_ph(self):
        error = self.steamTable.s_ph(1.0, 84.01181117) - 0.296813845
        self.assertLess(error, self.maxError, 's_ph not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_uV_p(self):
        error = self.steamTable.uV_p(1.0) - 2505.547389
        self.assertLess(error, self.maxError, 'uV_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_uL_p(self):
        error = self.steamTable.uL_p(1.0) - 417.332171
        self.assertLess(error, self.maxError, 'uL_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_uV_t(self):
        error = self.steamTable.uV_t(100.0) - 2506.015308
        self.assertLess(error, self.maxError, 'uV_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_uL_t(self):
        error = self.steamTable.uL_t(100.0) - 418.9933299
        self.assertLess(error, self.maxError, 'uL_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_u_pt(self):
        error = self.steamTable.u_pt(1.0, 100.0) - 2506.171426
        self.assertLess(error, self.maxError, 'u_pt not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_u_ph(self):
        error = self.steamTable.u_ph(1.0, 1000.0) - 956.2074342
        self.assertLess(error, self.maxError, 'u_ph not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_u_ps(self):
        error = self.steamTable.u_ps(1.0, 1.0) - 308.5082185
        self.assertLess(error, self.maxError, 'u_ps not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_CpV_p(self):
        error = self.steamTable.CpV_p(1.0) - 2.075938025
        self.assertLess(error, self.maxError, 'CpV_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_CpL_p(self):
        error = self.steamTable.CpL_p(1.0) - 4.216149431
        self.assertLess(error, self.maxError, 'CpL_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_CpV_t(self):
        error = self.steamTable.CpV_t(100.0) - 2.077491868
        self.assertLess(error, self.maxError, 'CpV_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_CpL_t(self):
        error = self.steamTable.CpL_t(100.0) - 4.216645119
        self.assertLess(error, self.maxError, 'CpL_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})
    def test_Cp_pt(self):
        error = self.steamTable.Cp_pt(1.0, 100.0) - 2.074108555
        self.assertLess(error, self.maxError, 'Cp_pt not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_Cp_ph(self):
        error = self.steamTable.Cp_ph(1.0, 200.0) - 4.17913573169
        self.assertLess(error, self.maxError, 'Cp_ph not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_Cp_ps(self):
        error = self.steamTable.Cp_ps(1.0, 1.0) - 4.190607038
        self.assertLess(error, self.maxError, 'Cp_ps not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_CvV_p(self):
        error = self.steamTable.CvV_p(1.0) - 1.552696979
        self.assertLess(error, self.maxError, 'CvV_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_CvL_p(self):
        error = self.steamTable.CvL_p(1.0) - 3.769699683
        self.assertLess(error, self.maxError, 'CvL_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_CvV_t(self):
        error = self.steamTable.CvV_t(100.0) - 1.553698696
        self.assertLess(error, self.maxError, 'CvV_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_CvL_t(self):
        error = self.steamTable.CvL_t(100.0) - 3.76770022
        self.assertLess(error, self.maxError, 'CvL_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_Cv_pt(self):
        error = self.steamTable.Cv_pt(1.0, 100.0) - 1.551397249
        self.assertLess(error, self.maxError, 'Cv_pt not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_Cv_ph(self):
        error = self.steamTable.Cv_ph(1.0, 200.0) - 4.035176364
        self.assertLess(error, self.maxError, 'Cv_ph not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_Cv_ps(self):
        error = self.steamTable.Cv_ps(1.0, 1.0) - 3.902919468
        self.assertLess(error, self.maxError, 'Cv_ps not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_wV_p(self):
        error = self.steamTable.wV_p(1.0) - 472.0541571
        self.assertLess(error, self.maxError, 'wV_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_wL_p(self):
        error = self.steamTable.wL_p(1.0) - 1545.451948
        self.assertLess(error, self.maxError, 'wL_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_wV_t(self):
        error = self.steamTable.wV_t(100.0) - 472.2559492
        self.assertLess(error, self.maxError, 'wV_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_wL_t(self):
        error = self.steamTable.wL_t(100.0) - 1545.092249
        self.assertLess(error, self.maxError, 'wL_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_w_pt(self):
        error = self.steamTable.w_pt(1.0, 100.0) - 472.3375235
        self.assertLess(error, self.maxError, 'w_pt not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_w_ph(self):
        error = self.steamTable.w_ph(1.0, 200.0) - 1542.682475
        self.assertLess(error, self.maxError, 'w_ph not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_w_ps(self):
        error = self.steamTable.w_ps(1.0, 1.0) - 1557.8585
        self.assertLess(error, self.maxError, 'w_ps not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})


    def test_my_pt(self):
        error = self.steamTable.my_pt(1.0, 100.0) - 1.22704e-05
        self.assertLess(error, self.maxError, 'my_pt not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})


    def test_my_ph(self):
        error = self.steamTable.my_ph(1.0, 100.0) - 0.000914003770302
        self.assertLess(error, self.maxError, 'my_ph not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})


    def test_my_ps(self):
        error = self.steamTable.my_ps(1.0, 1.0) - 0.000384222
        self.assertLess(error, self.maxError, 'my_ps not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_tcL_p(self):
        error = self.steamTable.tcL_p(1.0) - 0.677593822
        self.assertLess(error, self.maxError, 'tcL_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_tcV_p(self):
        error = self.steamTable.tcV_p(1.0) - 0.024753668
        self.assertLess(error, self.maxError, 'tcV_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_tcL_t(self):
        error = self.steamTable.tcL_t(25.0) - 0.607458162
        self.assertLess(error, self.maxError, 'tcL_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_tcV_t(self):
        error = self.steamTable.tcV_t(25.0) - 0.018326723
        self.assertLess(error, self.maxError, 'tcV_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})


    def test_tc_pt(self):
        error = self.steamTable.tc_pt(1.0, 25.0) - 0.607509806
        self.assertLess(error, self.maxError, 'tc_pt not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})


    def test_tc_ph(self):
        error = self.steamTable.tc_ph(1.0, 100.0) - 0.605710062
        self.assertLess(error, self.maxError, 'tc_ph not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})


    def test_tc_hs(self):
        error = self.steamTable.tc_hs(100.0, 0.34) - 0.606283124
        self.assertLess(error, self.maxError, 'tc_hs not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_st_t(self):
        error = self.steamTable.st_t(100.0) - 0.0589118685877
        self.assertLess(error, self.maxError, 'st_t not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_st_p(self):
        error = self.steamTable.st_p(1.0) - 0.058987784
        self.assertLess(error, self.maxError, 'st_p not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_x_ph(self):
        error = self.steamTable.x_ph(1.0, 1000.0) - 0.258055424
        self.assertLess(error, self.maxError, 'x_ph not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_x_ps(self):
        error = self.steamTable.x_ps(1.0, 4.0) - 0.445397961
        self.assertLess(error, self.maxError, 'x_ps not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_vx_ph(self):
        error = self.steamTable.vx_ph(1.0, 418.0) - 0.288493093
        self.assertLess(error, self.maxError, 'vx_ph not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

    def test_vx_ps(self):
        error = self.steamTable.vx_ps(1.0, 4.0) - 0.999233827
        self.assertLess(error, self.maxError, 'vx_ps not passed Error %(error)e allowed: %(max)e' % {'error' : error, 'max' : self.maxError})

def suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTest(loader.loadTestsFromTestCase(Region1Tester))
    suite.addTest(loader.loadTestsFromTestCase(Region2Tester))
    suite.addTest(loader.loadTestsFromTestCase(Region3Tester))
    suite.addTest(loader.loadTestsFromTestCase(Region4Tester))
    suite.addTest(loader.loadTestsFromTestCase(Region5Tester))
    suite.addTest(loader.loadTestsFromTestCase(CallFunctionTester))
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity = 2).run(suite())

# if __name__ == "__main__":
#    # import sys;sys.argv = ['', 'Test.testName']
#    unittest.main()
