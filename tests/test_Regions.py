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
    """tests for functions in region 1"""

    def setUp(self):
        self.max_error = 1e-8
        # The accumulated Error is bigger than the error of each single function
        self.max_matrix_error = 2e-8

    def tearDown(self):
        pass

    def test_pT_functions(self):
        """R7-97(2012) Tests to verify all functions with the Parameters p and T of Region 1"""
        # Table 5
        in_p = [3.0, 80.0, 3.0]
        in_T = [300.0, 300.0, 500.0]
        ref = [
            [0.100215168e-2, 0.971180894e-3, 0.120241800e-2],  # v
            [0.115331273e3, 0.184142828e3, 0.975542239e3],  # h
            [0.112324818e3, 0.106448356e3, 0.971934985e3],  # u
            [0.392294792, 0.368563852, 0.258041912e1],  # s
            [0.417301218e1, 0.401008987e1, 0.465580682e1],  # cp
            [0.150773921e4, 0.163469054e4, 0.124071337e4],  # w
        ]
        res = numpy.zeros((6, 3))
        for i, (p, T) in enumerate(zip(in_p, in_T)):
            res[0][i] = Region1.v1_pT(p, T)
            res[1][i] = Region1.h1_pT(p, T)
            res[2][i] = Region1.u1_pT(p, T)
            res[3][i] = Region1.s1_pT(p, T)
            res[4][i] = Region1.Cp1_pT(p, T)
            res[5][i] = Region1.w1_pT(p, T)

        Region1_error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(Region1_error, self.max_matrix_error, "Test of *(p,T) functions for Region 1 failed.")

    def test_ph_function(self):
        """R7-97(2012) Tests to verify all functions with the Parameters p and h of Region 1"""
        # Table 7
        in_p = [3.0, 80.0, 80.0]
        in_h = [500.0, 500.0, 1500.0]
        ref = [0.391798509e3, 0.378108626e3, 0.611041229e3]
        res = numpy.zeros(3)
        for i, (p, h) in enumerate(zip(in_p, in_h)):
            res[i] = Region1.T1_ph(p, h)

        T1_ph_error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(T1_ph_error, self.max_error, "Test of T(p,h) Function for Region 1 failed")

    def test_ps_function(self):
        """R7-97(2012) Tests to verify all functions with the Parameters p and s of Region 1"""
        # Table 9
        in_p = [3.0, 80.0, 80.0]
        in_s = [0.5, 0.5, 3.0]
        ref = [0.307842258e3, 0.309979785e3, 0.565899909e3]
        res = numpy.zeros(3)
        for i, (p, s) in enumerate(zip(in_p, in_s)):
            res[i] = Region1.T1_ps(p, s)

        T1_ps_error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(T1_ps_error, self.max_error, "Test of T(p,s) Function for Region 1 failed.")

    def test_hs_function(self):
        """Tests to verify all functions with the Parameters h and s of Region 1 by comparing the Results to IF-97 Page 6 Table 3"""
        # TODO
        # % Supplementary Release on Backward Equations
        # % for Pressure as a Function of Enthalpy and Entropy p(h, s)
        # % Table 3, Page 6
        in_h = [0.001, 90.0, 1500.0]
        in_s = [0.0, 0.0, 3.4]
        ref = [0.0009800980612, 91.929547272, 58.68294423]
        res = numpy.zeros(3)
        for i, (h, s) in enumerate(zip(in_h, in_s)):
            res[i] = Region1.p1_hs(h, s)

        p1_hs_error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(p1_hs_error, self.max_error, "Test of p(h,s) Function for Region 1 failed.")


class Region2Tester(unittest.TestCase):
    """tests for functions in region 2"""

    def setUp(self):
        self.max_error = 1e-8
        # The accumulated Error is bigger than the error of each single function
        self.max_matrix_error = 2e-8

    def tearDown(self):
        pass

    def test_pT_function(self):
        """R7-97(2012) Tests to verify all functions with the Parameters p and T of Region 2"""
        # Table 15
        in_p = [0.0035, 0.0035, 30.0]
        in_T = [300.0, 700.0, 700.0]
        ref = [
            [0.394913866e2, 0.923015898e2, 0.542946619e-2],  # v
            [0.254991145e4, 0.333568375e4, 0.263149474e4],  # h
            [0.241169160e4, 0.301262819e4, 0.246861076e4],  # u
            [0.852238967e1, 0.101749996e2, 0.517540298e1],  # s
            [0.191300162e1, 0.208141274e1, 0.103505092e2],  # cp
            [0.427920172e3, 0.644289068e3, 0.480386523e3],  # w
        ]
        res = numpy.zeros((6, 3))

        for i, (p, T) in enumerate(zip(in_p, in_T)):
            res[0][i] = Region2.v2_pT(p, T)
            res[1][i] = Region2.h2_pT(p, T)
            res[2][i] = Region2.u2_pT(p, T)
            res[3][i] = Region2.s2_pT(p, T)
            res[4][i] = Region2.Cp2_pT(p, T)
            res[5][i] = Region2.w2_pT(p, T)

        Region2_error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(Region2_error, self.max_matrix_error, "Test of *(p,T) Functions for Region 2 failed")

    def test_pT_meta_function(self):
        """R7-97(2012) Tests to verify all functions with the Parameters p and T of Metastable-Vapor Region 2"""
        # Table 18
        in_p = [1.0, 1.0, 1.5]
        in_T = [450.0, 440.0, 450.0]

        ref = [
            [0.192516540, 0.186212297, 0.121685206],  # v
            [0.276881115e4, 0.274015123e4, 0.272134539e4],  # h
            [0.257629461e4, 0.255393894e4, 0.253881758e4],  # u
            [0.656660377e1, 0.650218759e1, 0.629170440e1],  # s
            [0.276349265e1, 0.298166443e1, 0.362795578e1],  # cp
            [0.498408101e3, 0.489363295e3, 0.481941819e3],  # w
        ]
        res = numpy.zeros((6, 3))

        for i, (p, T) in enumerate(zip(in_p, in_T)):
            res[0][i] = Region2.v2_pT_meta(p, T)
            res[1][i] = Region2.h2_pT_meta(p, T)
            res[2][i] = Region2.u2_pT_meta(p, T)
            res[3][i] = Region2.s2_pT_meta(p, T)
            res[4][i] = Region2.Cp2_pT_meta(p, T)
            res[5][i] = Region2.w2_pT_meta(p, T)

        Region2_error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(Region2_error, self.max_matrix_error, "Test of *(p,T) Functions for metastable-vapor Region 2 failed.")

    def test_ph_function(self):
        """R7-97(2012) Tests to verify all functions with the Parameters p and h of subregion 2a, 2b and 2c"""
        # Table 24
        in_p = [0.001, 3.0, 3.0, 5.0, 5.0, 25.0, 40.0, 60.0, 60.0]
        in_h = [3000.0, 3000.0, 4000.0, 3500.0, 4000.0, 3500.0, 2700.0, 2700.0, 3200.0]
        ref = [
            0.534433241e3,
            0.575373370e3,
            0.101077577e4,
            0.801299102e3,
            0.101531583e4,
            0.875279054e3,
            0.743056411e3,
            0.791137067e3,
            0.882756860e3,
        ]
        res = numpy.zeros(9)
        for i, (p, h) in enumerate(zip(in_p, in_h)):
            res[i] = Region2.T2_ph(p, h)

        T2_ph_error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(T2_ph_error, self.max_matrix_error, "Test of T(p,h) Function for Region 2 failed")

    def test_ps_function(self):
        """R7-97(2012) Tests to verify all functions with the Parameters p and s of subregion 2a, 2b and 2c"""
        # Table 29
        in_p = [0.1, 0.1, 2.5, 8.0, 8.0, 90.0, 20.0, 80.0, 80.0]
        in_s = [7.5, 8.0, 8.0, 6.0, 7.5, 6.0, 5.75, 5.25, 5.75]
        ref = [
            0.399517097e3,
            0.514127081e3,
            0.103984917e4,
            0.600484040e3,
            0.106495556e4,
            0.103801126e4,
            0.697992849e3,
            0.854011484e3,
            0.949017998e3,
        ]
        res = numpy.zeros(9)
        for i, (p, s) in enumerate(zip(in_p, in_s)):
            res[i] = Region2.T2_ps(p, s)

        T2_ps_error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(T2_ps_error, self.max_matrix_error, "Test of T(p,s) Function for Region 2 failed")

    def test_hs_function(self):
        """Tests to verify all functions with the Parameters h and s of Region 2 by comparing the Results to IF-97 Page 6 Table 3"""
        # TODO
        # % Supplementary Release on Backward Equations for Pressure as a Function of Enthalpy and Entropy p(h, s)
        # % Table 3, Page 6
        h = [2800.0, 2800.0, 4100.0, 2800.0, 3600.0, 3600.0, 2800.0, 2800.0, 3400.0]
        s = [6.5, 9.5, 9.5, 6, 6, 7, 5.1, 5.8, 5.8]
        IF97 = [
            1.371012767,
            0.001879743844,
            0.1024788997,
            4.793911442,
            83.95519209,
            7.527161441,
            94.3920206,
            8.414574124,
            83.76903879,
        ]
        R2 = numpy.zeros(9)
        for i in range(9):
            R2[i] = Region2.p2_hs(h[i], s[i])

        p2_hs_error = numpy.sum(numpy.absolute((R2 - IF97) / IF97))
        self.assertLess(p2_hs_error, self.max_error, "Test of hs Function for Region 2 failed")


class Region3Tester(unittest.TestCase):
    """tests for functions in region 3"""

    def setUp(self):
        self.max_error = 1e-8
        # The accumulated Error is bigger than the error of each single function
        self.max_matrix_error = 2e-8

    def tearDown(self):
        pass

    def test_rhoT_function(self):
        """R7-97(2012) Tests to verify all functions with the Parameters rho and T of region 3"""
        # Table 33
        in_T = [650.0, 650.0, 750.0]
        in_rho = [500.0, 200.0, 500.0]
        ref = [
            [0.255837018e2, 0.222930643e2, 0.783095639e2],  # p
            [0.186343019e4, 0.237512401e4, 0.225868845e4],  # h
            [0.181226279e4, 0.226365868e4, 0.210206932e4],  # u
            [0.405427273e1, 0.485438792e1, 0.446971906e1],  # s
            [0.138935717e2, 0.446579342e2, 0.634165359e1],  # cp
            [0.502005554e3, 0.383444594e3, 0.760696041e3],  # w
        ]
        res = numpy.zeros((6, 3))
        for i, (rho, T) in enumerate(zip(in_rho, in_T)):
            res[0][i] = Region3.p3_rhoT(rho, T)
            res[1][i] = Region3.h3_rhoT(rho, T)
            res[2][i] = Region3.u3_rhoT(rho, T)
            res[3][i] = Region3.s3_rhoT(rho, T)
            res[4][i] = Region3.Cp3_rhoT(rho, T)
            res[5][i] = Region3.w3_rhoT(rho, T)

        Region3_error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(Region3_error, self.max_matrix_error, "Test of rhoT Function for Region 3 failed")

    def test_T_ph_function(self):
        """Tests to verify all temperature functions with the Parameters p and h of Region 3"""
        # TODO
        # % T3_ph
        p = [20.0, 50.0, 100.0, 20.0, 50.0, 100.0]
        h = [1700.0, 2000.0, 2100.0, 2500.0, 2400.0, 2700.0]
        IF97 = [
            629.3083892,
            690.5718338,
            733.6163014,
            641.8418053,
            735.1848618,
            842.0460876,
        ]
        R3 = numpy.zeros(6)
        for i in range(6):
            R3[i] = Region3.T3_ph(p[i], h[i])
        T3_ph_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(
            T3_ph_error,
            self.max_error,
            "Test of T(p,h) Function for Region 3 failed. Error was %(error)e allowed:"
            " %(max)e" % {"error": T3_ph_error, "max": self.max_error},
        )

    def test_v_ph_function(self):
        """Tests to verify all v functions with the Parameters p and h of Region 3"""
        # TODO
        # % v3_ph
        p = [20.0, 50.0, 100.0, 20.0, 50.0, 100.0]
        h = [1700.0, 2000.0, 2100.0, 2500.0, 2400.0, 2700.0]
        IF97 = [
            0.001749903962,
            0.001908139035,
            0.001676229776,
            0.006670547043,
            0.0028012445,
            0.002404234998,
        ]
        R3 = numpy.zeros(6)
        for i in range(6):
            R3[i] = Region3.v3_ph(p[i], h[i])
        v3_ph_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(
            v3_ph_error,
            1e-7,
            "Test of v(p,h) Function for Region 3 failed. Error was %(error)e allowed:" " %(max)e" % {"error": v3_ph_error, "max": 1e-7},
        )

    def test_T_ps_function(self):
        """Tests to verify all T functions with the Parameters p and s of Region 3"""
        # TODO
        # % T3_ps
        p = [20.0, 50.0, 100.0, 20.0, 50.0, 100.0]
        s = [3.7, 3.5, 4, 5, 4.5, 5.0]
        IF97 = [
            620.8841563,
            618.1549029,
            705.6880237,
            640.1176443,
            716.3687517,
            847.4332825,
        ]
        R3 = numpy.zeros(6)
        for i in range(6):
            R3[i] = Region3.T3_ps(p[i], s[i])
        T3_ps_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(
            T3_ps_error,
            self.max_error,
            "Test of T(p,s) Function for Region 3 failed. Error was %(error)e allowed:"
            " %(max)e" % {"error": T3_ps_error, "max": self.max_error},
        )

    def test_v_ps_function(self):
        """Tests to verify all v functions with the Parameters p and s of Region 3"""
        # TODO
        # % v3_ps
        p = [20.0, 50.0, 100.0, 20.0, 50.0, 100.0]
        s = [3.7, 3.5, 4.0, 5.0, 4.5, 5.0]
        IF97 = [
            0.001639890984,
            0.001423030205,
            0.001555893131,
            0.006262101987,
            0.002332634294,
            0.002449610757,
        ]
        R3 = numpy.zeros(6)
        for i in range(6):
            R3[i] = Region3.v3_ps(p[i], s[i])
        v3_ps_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(
            v3_ps_error,
            self.max_error,
            "Test of v(p,s) Function for Region 3 failed. Error was %(error)e allowed:"
            " %(max)e" % {"error": v3_ps_error, "max": self.max_error},
        )

    def test_hs_function(self):
        """Tests to verify all functions with the Parameters h and s of Region 3"""
        # TODO
        # % p3_hs
        h = [1700.0, 2000.0, 2100.0, 2500.0, 2400.0, 2700.0]
        s = [3.8, 4.2, 4.3, 5.1, 4.7, 5.0]
        IF97 = [
            25.55703246,
            45.40873468,
            60.7812334,
            17.20612413,
            63.63924887,
            88.39043281,
        ]
        R3 = numpy.zeros(6)
        for i in range(6):
            R3[i] = Region3.p3_hs(h[i], s[i])
        p3_hs_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(
            p3_hs_error,
            self.max_error,
            "Test of p(h,s) Function for Region 3 failed. Error was %(error)e allowed:"
            " %(max)e" % {"error": p3_hs_error, "max": self.max_error},
        )

    def test_pT_function(self):
        """Tests to verify all functions with the Parameters p and T of Region 3"""
        # TODO
        # % h3_pT (Iteration)
        p = [25.583702, 22.293064, 78.309564]
        T = [650.0, 650.0, 750.0]
        IF97 = [1863.271389, 2375.696155, 2258.626582]
        R3 = numpy.zeros(3)
        for i in range(3):
            R3[i] = Region3.h3_pT(p[i], T[i])
        h3_pT_error = numpy.sum(numpy.absolute((R3 - IF97) / IF97))
        self.assertLess(
            h3_pT_error,
            1e-6,
            "Test of h(p,T) Function for Region 3 failed. Error was %(error)e allowed:" " %(max)e" % {"error": h3_pT_error, "max": 1e-6},
        )


class Region4Tester(unittest.TestCase):
    """tests for functions in region 4"""

    def setUp(self):
        self.max_error = 1e-7
        # The accumulated Error is bigger than the error of each single function
        self.max_matrix_error = 2e-8

    def tearDown(self):
        pass

    def test_T_function(self):
        """R7-97(2012) Tests to verify all functions with the Parameters T of region 4"""
        # Table 35
        in_T = [300.0, 500.0, 600.0]
        ref = [0.353658941e-2, 0.263889776e1, 0.123443146e2]
        res = numpy.zeros(3)
        for i, T in enumerate(in_T):
            res[i] = Region4.p4_T(T)

        p4_t_error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(p4_t_error, self.max_error, "Test of p(T) Function for Region 4 failed")

    def test_p_functions(self):
        """R7-97(2012) Tests to verify all functions with the Parameters p of region 4"""
        in_p = [0.1, 1.0, 10.0]
        ref = [0.372755919e3, 0.453035632e3, 0.584149488e3]
        res = numpy.zeros(3)
        for i, p in enumerate(in_p):
            res[i] = Region4.T4_p(p)

        T4_p_error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(T4_p_error, self.max_error, "Test of T(p) Function for Region 4 failed")

    def test_s_functions(self):
        """Tests to verify all functions with the Parameters s of Region 4"""
        # TODO
        s = [1.0, 2.0, 3.0, 3.8, 4.0, 4.2, 7.0, 8.0, 9.0, 5.5, 5.0, 4.5]
        IF97 = [
            308.5509647,
            700.6304472,
            1198.359754,
            1685.025565,
            1816.891476,
            1949.352563,
            2723.729985,
            2599.04721,
            2511.861477,
            2687.69385,
            2451.623609,
            2144.360448,
        ]
        R4 = numpy.zeros(12)
        for i in range(12):
            R4[i] = Region4.h4_s(s[i])
        h4_s_error = numpy.sum(numpy.absolute((R4 - IF97) / IF97))
        self.assertLess(
            h4_s_error,
            self.max_error,
            "Test of h(s) Function for Region 4 failed. Error was %(error)e allowed:"
            " %(max)e" % {"error": h4_s_error, "max": self.max_error},
        )


class Region5Tester(unittest.TestCase):
    """tests for functions in region 5"""

    def setUp(self):
        self.max_error = 1e-8
        # The accumulated Error is bigger than the error of each single function
        self.max_matrix_error = 2e-8

    def tearDown(self):
        pass

    def test_pT_function(self):
        """R7-97(2012) Tests to verify all functions with the Parameters rhpo and T of region 5"""
        # Table 42
        in_p = [0.5, 8.0, 8.0]
        in_T = [1500.0, 1500.0, 2000.0]
        ref = [
            [0.138455090e1, 0.230761299e-1, 0.311385219e-1],  # v
            [0.521976855e4, 0.516723514e4, 0.657122604e4],  # h
            [0.452749310e4, 0.447495124e4, 0.563707038e4],  # u
            [0.965408875e1, 0.772970133e1, 0.853640523e1],  # s
            [0.261609445e1, 0.272724317e1, 0.288569882e1],  # cp
            [0.917068690e3, 0.928548002e3, 0.106736948e4],  # w
        ]
        res = numpy.zeros((6, 3))
        for i, (p, T) in enumerate(zip(in_p, in_T)):
            res[0][i] = Region5.v5_pT(p, T)
            res[1][i] = Region5.h5_pT(p, T)
            res[2][i] = Region5.u5_pT(p, T)
            res[3][i] = Region5.s5_pT(p, T)
            res[4][i] = Region5.Cp5_pT(p, T)
            res[5][i] = Region5.w5_pT(p, T)

        Region5_error = numpy.sum(numpy.absolute((res - ref) / ref))
        self.assertLess(Region5_error, self.max_matrix_error, "Test of p,T Function for Region 5 failed")

    def test_ph_function(self):
        """Tests to verify all functions with the Parameters p and h of Region 5"""
        # TODO
        # %T5_ph (Iteration)
        p = [0.5, 8.0, 8.0]
        h = [5219.76331549428, 5206.09634477373, 6583.80290533381]
        IF97 = [1500.0, 1500.0, 2000.0]
        R5 = numpy.zeros(3)
        for i in range(3):
            R5[i] = Region5.T5_ph(p[i], h[i])
        T5_ph_error = numpy.sum(numpy.absolute((R5 - IF97) / IF97))
        self.assertLess(
            T5_ph_error,
            self.max_error,
            "Test of T(p,h) Function for Region 5 failed. Error was %(error)e allowed:"
            " %(max)e" % {"error": T5_ph_error, "max": self.max_error},
        )

    def test_ps_function(self):
        """Tests to verify all functions with the Parameters p and s of Region 5"""
        # TODO
        # %T5_ps (Iteration)
        p = [0.5, 8.0, 8.0]
        s = [9.65408430982588, 8.36546724495503, 9.15671044273249]
        IF97 = [1500.0, 1500.0, 2000.0]
        R5 = numpy.zeros(3)
        for i in range(3):
            R5[i] = Region5.T5_ps(p[i], s[i])
        T5_ps_error = numpy.sum(numpy.absolute((R5 - IF97) / IF97))
        self.assertLess(
            T5_ps_error,
            1e-4,
            "Test of T(p,s) Function for Region 5 failed. Error was %(error)e allowed:" " %(max)e" % {"error": T5_ps_error, "max": 1e-4},
        )
