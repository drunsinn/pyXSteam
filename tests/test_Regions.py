# -*- coding: utf-8 -*-
"""
These Tests are taken form the original XSteam Matlab Script.
Some Errors are calculated with the help of numpy matrix functions.
Due to some rounding Errors the max. allowedError had to be increased....
"""

import unittest
import numpy
from pyXSteam.Regions import Region1, Region2, Region3, Region4, Region5


class Region1Tester(unittest.TestCase):

    def setUp(self):
        self.maxError = 1E-8
        self.maxMatrixError = 2E-8  # The accumulated Error is bigger than the error of each single function

    def tearDown(self):
        pass

    def test_pT_functions(self):
        """Tests to verify all functions with the Parameters p and T of Region 1 by comparing the Results to IF-97 Page 9 Table 5"""
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
        R1 = numpy.zeros((6, 3))

        for i in range(0, 3):
            R1[0][i] = Region1.v1_pT(p[i], T[i])
            R1[1][i] = Region1.h1_pT(p[i], T[i])
            R1[2][i] = Region1.u1_pT(p[i], T[i])
            R1[3][i] = Region1.s1_pT(p[i], T[i])
            R1[4][i] = Region1.Cp1_pT(p[i], T[i])
            R1[5][i] = Region1.w1_pT(p[i], T[i])

        Region1_error = numpy.sum(numpy.absolute((R1 - IF97) / IF97))
        self.assertLess(Region1_error, self.maxMatrixError, 'Test of p,T Functions for Region 1 failed. Error was %(error)e allowed: %(max)e' % {'error': Region1_error, 'max': self.maxMatrixError})

    def test_ph_function(self):
        """Tests to verify all functions with the Parameters p and h of Region 1 by comparing the Results to IF-97 Page 11 Table 7"""
        # % IF - 97 Table 7, Page 11
        p = [3.0, 80.0, 80.0]
        h = [500.0, 500.0, 1500.0]
        IF97 = [391.798509, 378.108626, 611.041229]
        R1 = numpy.zeros(3)
        for i in range(0, 3):
            R1[i] = Region1.T1_ph(p[i], h[i])

        T1_ph_error = numpy.sum(numpy.absolute((R1 - IF97) / IF97))
        self.assertLess(T1_ph_error, self.maxError, 'Test of T(p,h) Function for Region 1 failed. Error was %(error)e allowed: %(max)e' % {'error': T1_ph_error, 'max': self.maxError})

    def test_ps_function(self):
        """Tests to verify all functions with the Parameters p and s of Region 1 by comparing the Results to IF-97 Page 12 Table 9"""
        # % IF - 97 Table 9, Page 12
        p = [3.0, 80.0, 80.0]
        s = [0.5, 0.5, 3.0]
        IF97 = [307.842258, 309.979785, 565.899909]
        R1 = numpy.zeros(3)
        for i in range(0, 3):
            R1[i] = Region1.T1_ps(p[i], s[i])

        T1_ps_error = numpy.sum(numpy.absolute((R1 - IF97) / IF97))
        self.assertLess(T1_ps_error, self.maxError, 'Test of T(p,s) Function for Region 1 failed. Error was %(error)e allowed: %(max)e' % {'error': T1_ps_error, 'max': self.maxError})

    def test_hs_function(self):
        """Tests to verify all functions with the Parameters h and s of Region 1 by comparing the Results to IF-97 Page 6 Table 3"""
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
        self.assertLess(p1_hs_error, self.maxError, 'Test of p(h,s) Function for Region 1 failed. Error was %(error)e allowed: %(max)e' % {'error': p1_hs_error, 'max': self.maxError})


class Region2Tester(unittest.TestCase):

    def setUp(self):
        self.maxError = 1E-8
        self.maxMatrixError = 2E-8  # The accumulated Error is bigger than the error of each single function

    def tearDown(self):
        pass

    def test_pT_function(self):
        """ Tests to verify all functions with the Parameters p and T of Region 2 by comparing the Results to IF-97 Page 17 Table 15"""
        # %* 7.2 Verifiy region 2
        # % IF-97 Table 15, Page 17
        p = [0.0035, 0.0035, 30.0]
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
        self.assertLess(Region2_error, self.maxMatrixError, 'Test of p,T Functions for Region 2 failed. Error was %(error)e allowed: %(max)e' % {'error': Region2_error, 'max': self.maxMatrixError})

    def test_ph_function(self):
        """ Tests to verify all functions with the Parameters p and h of Region 2 by comparing the Results to IF-97 Page 25 Table 24"""
        # % IF - 97 Table 24, Page 25
        p = [0.01 / 10, 30 / 10, 30 / 10, 50 / 10, 50 / 10, 250 / 10, 400 / 10, 600 / 10, 600 / 10]
        h = [3000.0, 3000.0, 4000.0, 3500.0, 4000.0, 3500.0, 2700.0, 2700.0, 3200.0]
        IF97 = [534.433241, 575.37337, 1010.77577, 801.299102, 1015.31583, 875.279054, 743.056411, 791.137067, 882.75686]
        R2 = numpy.zeros(9)
        for i in range(0, 9):
            R2[i] = Region2.T2_ph(p[i], h[i])
        T2_ph_error = numpy.sum(numpy.absolute((R2 - IF97) / IF97))
        self.assertLess(T2_ph_error, self.maxMatrixError, 'Test of ph Function for Region 2 failed. Error was %(error)e allowed: %(max)e' % {'error': T2_ph_error, 'max': self.maxMatrixError})

    def test_ps_function(self):
        """ Tests to verify all functions with the Parameters p and s of Region 2 by comparing the Results to IF-97 Page 29 Table 29"""
        # % IF - 97 Table 29, Page 29
        p = [0.1, 0.1, 2.5, 8.0, 8.0, 90.0, 20.0, 80.0, 80.0]
        s = [7.5, 8.0, 8.0, 6.0, 7.5, 6.0, 5.75, 5.25, 5.75]
        IF97 = [399.517097, 514.127081, 1039.84917, 600.48404, 1064.95556, 1038.01126, 697.992849, 854.011484, 949.017998]
        R2 = numpy.zeros(9)
        for i in range(0, 9):
            R2[i] = Region2.T2_ps(p[i], s[i])
        T2_ps_error = numpy.sum(numpy.absolute((R2 - IF97) / IF97))
        self.assertLess(T2_ps_error, self.maxMatrixError, 'Test of ps Function for Region 2 failed. Error was %(error)e allowed: %(max)e' % {'error': T2_ps_error, 'max': self.maxMatrixError})

    def test_hs_function(self):
        """ Tests to verify all functions with the Parameters h and s of Region 2 by comparing the Results to IF-97 Page 6 Table 3"""
        # % Supplementary Release on Backward Equations for Pressure as a Function of Enthalpy and Entropy p(h, s)
        # % Table 3, Page 6
        h = [2800.0, 2800.0, 4100.0, 2800.0, 3600.0, 3600.0, 2800.0, 2800.0, 3400.0]
        s = [6.5, 9.5, 9.5, 6, 6, 7, 5.1, 5.8, 5.8]
        IF97 = [1.371012767, 0.001879743844, 0.1024788997, 4.793911442, 83.95519209, 7.527161441, 94.3920206, 8.414574124, 83.76903879]
        R2 = numpy.zeros(9)
        for i in range(0, 9):
            R2[i] = Region2.p2_hs(h[i], s[i])
        p2_hs_error = numpy.sum(numpy.absolute((R2 - IF97) / IF97))
        self.assertLess(p2_hs_error, self.maxError, 'Test of hs Function for Region 2 failed. Error was %(error)e allowed: %(max)e' % {'error': p2_hs_error, 'max': self.maxError})


class Region3Tester(unittest.TestCase):

    def setUp(self):
        self.maxError = 1E-8
        self.maxMatrixError = 2E-8  # The accumulated Error is bigger than the error of each single function

    def tearDown(self):
        pass

    def test_rhoT_function(self):
        """ Tests to verify all functions with the Parameters rho and T of Region 3 by comparing the Results to IF-97 Page 32 Table 33"""
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
        self.assertLess(Region3_error, self.maxMatrixError, 'Test of rhoT Function for Region 3 failed. Error was %(error)e allowed: %(max)e' % {'error': Region3_error, 'max': self.maxMatrixError})

    def test_T_ph_function(self):
        """ Tests to verify all temperature functions with the Parameters p and h of Region 3"""
        # % T3_ph
        p = [20.0, 50.0, 100.0, 20.0, 50.0, 100.0]
        h = [1700.0, 2000.0, 2100.0, 2500.0, 2400.0, 2700.0]
        IF97 = [629.3083892, 690.5718338, 733.6163014, 641.8418053, 735.1848618, 842.0460876]
        R3 = numpy.zeros(6)
        for i in range(0, 6):
            R3[i] = Region3.T3_ph(p[i], h[i])
        T3_ph_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(T3_ph_error, self.maxError, 'Test of T(p,h) Function for Region 3 failed. Error was %(error)e allowed: %(max)e' % {'error': T3_ph_error, 'max': self.maxError})

    def test_v_ph_function(self):
        """ Tests to verify all v functions with the Parameters p and h of Region 3"""
        # % v3_ph
        p = [20.0, 50.0, 100.0, 20.0, 50.0, 100.0]
        h = [1700.0, 2000.0, 2100.0, 2500.0, 2400.0, 2700.0]
        IF97 = [0.001749903962, 0.001908139035, 0.001676229776, 0.006670547043, 0.0028012445, 0.002404234998]
        R3 = numpy.zeros(6)
        for i in range(0, 6):
            R3[i] = Region3.v3_ph(p[i], h[i])
        v3_ph_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(v3_ph_error, 1E-7, 'Test of v(p,h) Function for Region 3 failed. Error was %(error)e allowed: %(max)e' % {'error': v3_ph_error, 'max': 1E-7})

    def test_T_ps_function(self):
        """ Tests to verify all T functions with the Parameters p and s of Region 3"""
        # % T3_ps
        p = [20.0, 50.0, 100.0, 20.0, 50.0, 100.0]
        s = [3.7, 3.5, 4, 5, 4.5, 5.0]
        IF97 = [620.8841563, 618.1549029, 705.6880237, 640.1176443, 716.3687517, 847.4332825]
        R3 = numpy.zeros(6)
        for i in range(0, 6):
            R3[i] = Region3.T3_ps(p[i], s[i])
        T3_ps_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(T3_ps_error, self.maxError, 'Test of T(p,s) Function for Region 3 failed. Error was %(error)e allowed: %(max)e' % {'error': T3_ps_error, 'max': self.maxError})

    def test_v_ps_function(self):
        """ Tests to verify all v functions with the Parameters p and s of Region 3"""
        # % v3_ps
        p = [20.0, 50.0, 100.0, 20.0, 50.0, 100.0]
        s = [3.7, 3.5, 4.0, 5.0, 4.5, 5.0]
        IF97 = [0.001639890984, 0.001423030205, 0.001555893131, 0.006262101987, 0.002332634294, 0.002449610757]
        R3 = numpy.zeros(6)
        for i in range(0, 6):
            R3[i] = Region3.v3_ps(p[i], s[i])
        v3_ps_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(v3_ps_error, self.maxError, 'Test of v(p,s) Function for Region 3 failed. Error was %(error)e allowed: %(max)e' % {'error': v3_ps_error, 'max': self.maxError})

    def test_hs_function(self):
        """ Tests to verify all functions with the Parameters h and s of Region 3"""
        # % p3_hs
        h = [1700.0, 2000.0, 2100.0, 2500.0, 2400.0, 2700.0]
        s = [3.8, 4.2, 4.3, 5.1, 4.7, 5.0]
        IF97 = [25.55703246, 45.40873468, 60.7812334, 17.20612413, 63.63924887, 88.39043281]
        R3 = numpy.zeros(6)
        for i in range(0, 6):
            R3[i] = Region3.p3_hs(h[i], s[i])
        p3_hs_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(p3_hs_error, self.maxError, 'Test of p(h,s) Function for Region 3 failed. Error was %(error)e allowed: %(max)e' % {'error': p3_hs_error, 'max': self.maxError})

    def test_pT_function(self):
        """ Tests to verify all functions with the Parameters p and T of Region 3"""
        # % h3_pT (Iteration)
        p = [25.583702, 22.293064, 78.309564]
        T = [650.0, 650.0, 750.0]
        IF97 = [1863.271389, 2375.696155, 2258.626582]
        R3 = numpy.zeros(3)
        for i in range(0, 3):
            R3[i] = Region3.h3_pT(p[i], T[i])
        h3_pT_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(h3_pT_error, 1E-6, 'Test of h(p,T) Function for Region 3 failed. Error was %(error)e allowed: %(max)e' % {'error': h3_pT_error, 'max':  1E-6})


class Region4Tester(unittest.TestCase):

    def setUp(self):
        self.maxError = 1E-7
        self.maxMatrixError = 2E-8  # The accumulated Error is bigger than the error of each single function

    def tearDown(self):
        pass

    def test_T_function(self):
        """ Tests to verify all functions with the Parameter T of Region 4 by comparing the Results to IF-97 Page 34 Table 35"""
        # %Saturation pressure, If97, Table 35, Page 34
        T = [300.0, 500.0, 600.0]
        IF97 = [0.00353658941, 2.63889776, 12.3443146]
        R4 = numpy.zeros(3)
        for i in range(0, 3):
            R4[i] = Region4.p4_T(T[i])
        p4_t_error = numpy.sum(numpy.absolute((R4 - IF97) / IF97))
        self.assertLess(p4_t_error, self.maxError, 'Test of p(T) Function for Region 4 failed. Error was %(error)e allowed: %(max)e' % {'error': p4_t_error, 'max': self.maxError})

    def test_p_functions(self):
        """ Tests to verify all functions with the Parameters p of Region 4"""
        p = [0.1, 1.0, 10.0]
        IF97 = [372.755919, 453.035632, 584.149488]
        R4 = numpy.zeros(3)
        for i in range(0, 3):
            R4[i] = Region4.T4_p(p[i])
        T4_p_error = numpy.sum(numpy.absolute((R4 - IF97) / IF97))
        self.assertLess(T4_p_error, self.maxError, 'Test of T(p) Function for Region 4 failed. Error was %(error)e allowed: %(max)e' % {'error': T4_p_error, 'max': self.maxError})

    def test_s_functions(self):
        """ Tests to verify all functions with the Parameters s of Region 4"""
        s = [1.0, 2.0, 3.0, 3.8, 4.0, 4.2, 7.0, 8.0, 9.0, 5.5, 5.0, 4.5]
        IF97 = [308.5509647, 700.6304472, 1198.359754, 1685.025565, 1816.891476, 1949.352563, 2723.729985, 2599.04721, 2511.861477, 2687.69385, 2451.623609, 2144.360448]
        R4 = numpy.zeros(12)
        for i in range(0, 12):
            R4[i] = Region4.h4_s(s[i])
        h4_s_error = numpy.sum(numpy.absolute((R4 - IF97) / IF97))
        self.assertLess(h4_s_error, self.maxError, 'Test of h(s) Function for Region 4 failed. Error was %(error)e allowed: %(max)e' % {'error': h4_s_error, 'max': self.maxError})


class Region5Tester(unittest.TestCase):

    def setUp(self):
        self.maxError = 1E-8
        self.maxMatrixError = 2E-8  # The accumulated Error is bigger than the error of each single function

    def tearDown(self):
        pass

    def test_pT_function(self):
        """ Tests to verify all functions with the Parameters p and T of Region 5 by comparing the Results to IF-97 Page 39 Table 42"""
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
        self.assertLess(Region5_error, self.maxMatrixError, 'Test of p,T Function for Region 5 failed. Error was %(error)e allowed: %(max)e' % {'error': Region5_error, 'max': self.maxMatrixError})

    def test_ph_function(self):
        """ Tests to verify all functions with the Parameters p and h of Region 5"""
        # %T5_ph (Iteration)
        p = [0.5, 8.0, 8.0]
        h = [5219.76331549428, 5206.09634477373, 6583.80290533381]
        IF97 = [1500.0, 1500.0, 2000.0]
        R5 = numpy.zeros(3)
        for i in range(0, 3):
            R5[i] = Region5.T5_ph(p[i], h[i])
        T5_ph_error = numpy.sum(numpy.absolute((R5 - IF97) / IF97))
        self.assertLess(T5_ph_error, self.maxError, 'Test of T(p,h) Function for Region 5 failed. Error was %(error)e allowed: %(max)e' % {'error': T5_ph_error, 'max': self.maxError})

    def test_ps_function(self):
        """ Tests to verify all functions with the Parameters p and s of Region 5"""
        # %T5_ps (Iteration)
        p = [0.5, 8.0, 8.0]
        s = [9.65408430982588, 8.36546724495503, 9.15671044273249]
        IF97 = [1500.0, 1500.0, 2000.0]
        R5 = numpy.zeros(3)
        for i in range(0, 3):
            R5[i] = Region5.T5_ps(p[i], s[i])
        T5_ps_error = numpy.sum(numpy.absolute((R5 - IF97) / IF97))
        self.assertLess(T5_ps_error, 1E-4, 'Test of T(p,s) Function for Region 5 failed. Error was %(error)e allowed: %(max)e' % {'error': T5_ps_error, 'max': 1E-4})
