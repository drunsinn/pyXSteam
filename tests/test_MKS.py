# -*- coding: utf-8 -*-
"""
These Tests are taken form the original XSteam Matlab Script.
Some Errors are calculated with the help of numpy matrix functions.
"""

import unittest
from pyXSteam import XSteam, UnitSystem


class MKSFunctionTester(unittest.TestCase):
    def setUp(self):
        self.max_error = 1e-6
        self.steam_table = XSteam(UnitSystem.MKS)

    def tearDown(self):
        pass

    def test_tsat_p(self):
        error = abs(self.steam_table.tsat_p(1.0) - 99.60591861)
        self.assertLess(error, self.max_error, "test of tsat_p not passed")

    def test_t_ph(self):
        error = abs(self.steam_table.t_ph(1.0, 100.0) - 23.84481908)
        self.assertLess(error, self.max_error, "test of t_ph not passed")

    def test_t_ps(self):
        error = abs(self.steam_table.t_ps(1.0, 1.0) - 73.70859421)
        self.assertLess(error, self.max_error, "test of t_ps not passed")

    def test_t_hs(self):
        error = abs(self.steam_table.t_hs(100.0, 0.2) - 13.84933511)
        self.assertLess(error, self.max_error, "test of t_hs not passed")

    def test_psat_t(self):
        error = abs(self.steam_table.psat_t(100.0) - 1.014179779)
        self.assertLess(error, self.max_error, "test of psat_t not passed")

    def test_p_hs(self):
        error = abs(self.steam_table.p_hs(84.0, 0.296) - 2.295498269)
        self.assertLess(error, self.max_error, "test of p_hs not passed")

    def test_hV_p(self):
        error = abs(self.steam_table.hV_p(1.0) - 2674.949641)
        self.assertLess(error, self.max_error, "test of hV_p not passed")

    def test_hL_p(self):
        error = abs(self.steam_table.hL_p(1.0) - 417.4364858)
        self.assertLess(error, self.max_error, "test of hL_p not passed")

    def test_hV_t(self):
        error = abs(self.steam_table.hV_t(100.0) - 2675.572029)
        self.assertLess(error, self.max_error, "test of hV_t not passed")

    def test_hL_t(self):
        error = abs(self.steam_table.hL_t(100.0) - 419.099155)
        self.assertLess(error, self.max_error, "test of hL_t not passed")

    def test_h_pt(self):
        error = abs(self.steam_table.h_pt(1.0, 20.0) - 84.01181117)
        self.assertLess(error, self.max_error, "test of h_pt not passed")

    def test_h_ps(self):
        error = abs(self.steam_table.h_ps(1.0, 1.0) - 308.6107171)
        self.assertLess(error, self.max_error, "test of h_ps not passed")

    def test_h_px(self):
        error = abs(self.steam_table.h_px(1.0, 0.5) - 1546.193063)
        self.assertLess(error, self.max_error, "test of h_px not passed")

    def test_h_prho(self):
        error = abs(self.steam_table.h_prho(1.0, 2.0) - 1082.773391)
        self.assertLess(error, self.max_error, "test of h_prho not passed")

    def test_h_tx(self):
        error = abs(self.steam_table.h_tx(100.0, 0.5) - 1547.33559211)
        self.assertLess(error, self.max_error, "test of h_tx not passed")

    def test_vV_p(self):
        error = abs(self.steam_table.vV_p(1.0) - 1.694022523)
        self.assertLess(error, self.max_error, "test of vV_p not passed")

    def test_vL_p(self):
        error = abs(self.steam_table.vL_p(1.0) - 0.001043148)
        self.assertLess(error, self.max_error, "test of vL_p not passed")

    def test_vV_t(self):
        error = abs(self.steam_table.vV_t(100.0) - 1.671860601)
        self.assertLess(error, self.max_error, "test of vV_t not passed")

    def test_vL_t(self):
        error = abs(self.steam_table.vL_t(100.0) - 0.001043455)
        self.assertLess(error, self.max_error, "test of vL_t not passed")

    def test_v_pt(self):
        error = abs(self.steam_table.v_pt(1.0, 100.0) - 1.695959407)
        self.assertLess(error, self.max_error, "test of v_pt not passed")

    def test_v_ph(self):
        error = abs(self.steam_table.v_ph(1.0, 1000.0) - 0.437925658)
        self.assertLess(error, self.max_error, "test of v_ph not passed")

    def test_v_ps(self):
        error = abs(self.steam_table.v_ps(1.0, 5.0) - 1.03463539)
        self.assertLess(error, self.max_error, "test of v_ps not passed")

    def test_rhoV_p(self):
        error = abs(self.steam_table.rhoV_p(1.0) - 0.590310924)
        self.assertLess(error, self.max_error, "test of rhoV_p not passed")

    def test_rhoL_p(self):
        error = abs(self.steam_table.rhoL_p(1.0) - 958.6368897)
        self.assertLess(error, self.max_error, "test of rhoL_p not passed")

    def test_rhoV_t(self):
        error = abs(self.steam_table.rhoV_t(100.0) - 0.598135993)
        self.assertLess(error, self.max_error, "test of rhoV_t not passed")

    def test_rhoL_t(self):
        error = abs(self.steam_table.rhoL_t(100.0) - 958.3542773)
        self.assertLess(error, self.max_error, "test of rhoL_t not passed")

    def test_rho_pt(self):
        error = abs(self.steam_table.rho_pt(1.0, 100.0) - 0.589636754)
        self.assertLess(error, self.max_error, "test of rho_pt not passed")

    def test_rho_ph(self):
        error = abs(self.steam_table.rho_ph(1.0, 1000.0) - 2.283492601)
        self.assertLess(error, self.max_error, "test of rho_ph not passed")

    def test_rho_ps(self):
        error = abs(self.steam_table.rho_ps(1.0, 1.0) - 975.6236788)
        self.assertLess(error, self.max_error, "test of rho_ps not passed")

    def test_sV_p(self):
        error = abs(self.steam_table.sV_p(0.006117) - 9.155465556)
        self.assertLess(error, self.max_error, "test of sV_p not passed")

    def test_sL_p(self):
        error = abs(self.steam_table.sL_p(0.0061171) - 1.8359e-05)
        self.assertLess(error, self.max_error, "test of sL_p not passed")

    def test_sV_t(self):
        error = abs(self.steam_table.sV_t(0.0001) - 9.155756716)
        self.assertLess(error, self.max_error, "test of sV_t not passed")

    def test_sL_t(self):
        error = abs(self.steam_table.sL_t(100.0) - 1.307014328)
        self.assertLess(error, self.max_error, "test of sL_t not passed")

    def test_s_pt(self):
        error = abs(self.steam_table.s_pt(1.0, 20.0) - 0.296482921)
        self.assertLess(error, self.max_error, "test of s_pt not passed")

    def test_s_ph(self):
        error = abs(self.steam_table.s_ph(1.0, 84.01181117) - 0.296813845)
        self.assertLess(error, self.max_error, "test of s_ph not passed")

    def test_uV_p(self):
        error = abs(self.steam_table.uV_p(1.0) - 2505.547389)
        self.assertLess(error, self.max_error, "test of uV_p not passed")

    def test_uL_p(self):
        error = abs(self.steam_table.uL_p(1.0) - 417.332171)
        self.assertLess(error, self.max_error, "test of uL_p not passed")

    def test_uV_t(self):
        error = abs(self.steam_table.uV_t(100.0) - 2506.015308)
        self.assertLess(error, self.max_error, "test of uV_t not passed")

    def test_uL_t(self):
        error = abs(self.steam_table.uL_t(100.0) - 418.9933299)
        self.assertLess(error, self.max_error, "test of uL_t not passed")

    def test_u_pt(self):
        error = abs(self.steam_table.u_pt(1.0, 100.0) - 2506.171426)
        self.assertLess(error, self.max_error, "test of u_pt not passed")

    def test_u_ph(self):
        error = abs(self.steam_table.u_ph(1.0, 1000.0) - 956.2074342)
        self.assertLess(error, self.max_error, "test of u_ph not passed")

    def test_u_ps(self):
        error = abs(self.steam_table.u_ps(1.0, 1.0) - 308.5082185)
        self.assertLess(error, self.max_error, "test of u_ps not passed")

    def test_CpV_p(self):
        error = abs(self.steam_table.CpV_p(1.0) - 2.075938025)
        self.assertLess(error, self.max_error, "test of cpV_p not passed")

    def test_CpL_p(self):
        error = abs(self.steam_table.CpL_p(1.0) - 4.216149431)
        self.assertLess(error, self.max_error, "test of cpL_p not passed")

    def test_CpV_t(self):
        error = abs(self.steam_table.CpV_t(100.0) - 2.077491868)
        self.assertLess(error, self.max_error, "test of cpV_t not passed")

    def test_CpL_t(self):
        error = abs(self.steam_table.CpL_t(100.0) - 4.216645119)
        self.assertLess(error, self.max_error, "test of cpL_t not passed")

    def test_Cp_pt(self):
        error = abs(self.steam_table.Cp_pt(1.0, 100.0) - 2.074108555)
        self.assertLess(error, self.max_error, "test of cp_pt not passed")

    def test_Cp_ph(self):
        error = abs(self.steam_table.Cp_ph(1.0, 200.0) - 4.17913573169)
        self.assertLess(error, self.max_error, "test of Cp_ph not passed")

    def test_Cp_ps(self):
        error = abs(self.steam_table.Cp_ps(1.0, 1.0) - 4.190607038)
        self.assertLess(error, self.max_error, "test of Cp_ps not passed")

    def test_CvV_p(self):
        error = abs(self.steam_table.CvV_p(1.0) - 1.552696979)
        self.assertLess(error, self.max_error, "test of CvV_p not passed")

    def test_CvL_p(self):
        error = abs(self.steam_table.CvL_p(1.0) - 3.769699683)
        self.assertLess(error, self.max_error, "test of CvL_p not passed")

    def test_CvV_t(self):
        error = abs(self.steam_table.CvV_t(100.0) - 1.553698696)
        self.assertLess(error, self.max_error, "test of CvV_t not passed")

    def test_CvL_t(self):
        error = abs(self.steam_table.CvL_t(100.0) - 3.76770022)
        self.assertLess(error, self.max_error, "test of CvL_t not passed")

    def test_Cv_pt(self):
        error = abs(self.steam_table.Cv_pt(1.0, 100.0) - 1.551397249)
        self.assertLess(error, self.max_error, "test of Cv_pt not passed")

    def test_Cv_ph(self):
        error = abs(self.steam_table.Cv_ph(1.0, 200.0) - 4.035176364)
        self.assertLess(error, self.max_error, "test of Cv_ph not passed")

    def test_Cv_ps(self):
        error = abs(self.steam_table.Cv_ps(1.0, 1.0) - 3.902919468)
        self.assertLess(error, self.max_error, "test of Cv_ps not passed")

    def test_wV_p(self):
        error = abs(self.steam_table.wV_p(1.0) - 472.0541571)
        self.assertLess(error, self.max_error, "test of wV_p not passed")

    def test_wL_p(self):
        error = abs(self.steam_table.wL_p(1.0) - 1545.451948)
        self.assertLess(error, self.max_error, "test of wL_p not passed")

    def test_wV_t(self):
        error = abs(self.steam_table.wV_t(100.0) - 472.2559492)
        self.assertLess(error, self.max_error, "test of wV_t not passed")

    def test_wL_t(self):
        error = abs(self.steam_table.wL_t(100.0) - 1545.092249)
        self.assertLess(error, self.max_error, "test of wL_t not passed")

    def test_w_pt(self):
        error = abs(self.steam_table.w_pt(1.0, 100.0) - 472.3375235)
        self.assertLess(error, self.max_error, "test of w_pt not passed")

    def test_w_ph(self):
        error = abs(self.steam_table.w_ph(1.0, 200.0) - 1542.682475)
        self.assertLess(error, self.max_error, "test of w_ph not passed")

    def test_w_ps(self):
        error = abs(self.steam_table.w_ps(1.0, 1.0) - 1557.858535)
        self.assertLess(error, self.max_error, "test of w_ps not passed")

    def test_my_pt(self):
        error = abs(self.steam_table.my_pt(1.0, 100.0) - 1.22704e-05)
        self.assertLess(error, self.max_error, "test of my_pt not passed")

    def test_my_ph(self):
        error = abs(self.steam_table.my_ph(1.0, 100.0) - 0.000914003770302)
        self.assertLess(error, self.max_error, "test of my_ph not passed")

    def test_my_ps(self):
        error = abs(self.steam_table.my_ps(1.0, 1.0) - 0.000384222)
        self.assertLess(error, self.max_error, "test of my_ps not passed")

    def test_tcL_p(self):
        error = abs(self.steam_table.tcL_p(1.0) - 0.677593822)
        self.assertLess(error, self.max_error, "test of tcL_p not passed")

    def test_tcV_p(self):
        error = abs(self.steam_table.tcV_p(1.0) - 0.024753668)
        self.assertLess(error, self.max_error, "test of tcV_p not passed")

    def test_tcL_t(self):
        error = abs(self.steam_table.tcL_t(25.0) - 0.607458162)
        self.assertLess(error, self.max_error, "test of tcL_t not passed")

    def test_tcV_t(self):
        error = abs(self.steam_table.tcV_t(25.0) - 0.018326723)
        self.assertLess(error, self.max_error, "test of tcV_t not passed")

    def test_tc_pt(self):
        error = abs(self.steam_table.tc_pt(1.0, 25.0) - 0.607509806)
        self.assertLess(error, self.max_error, "test of tc_pt not passed")

    def test_tc_ph(self):
        error = abs(self.steam_table.tc_ph(1.0, 100.0) - 0.605710062)
        self.assertLess(error, self.max_error, "test of tc_ph not passed")

    def test_tc_hs(self):
        error = abs(self.steam_table.tc_hs(100.0, 0.34) - 0.606283124)
        self.assertLess(error, self.max_error, "test of tc_hs not passed")

    def test_st_t(self):
        error = abs(self.steam_table.st_t(100.0) - 58.9118685877)
        self.assertLess(error, self.max_error, "test of st_t not passed")

    def test_st_p(self):
        error = abs(self.steam_table.st_p(1.0) - 58.987784)
        self.assertLess(error, self.max_error, "test of st_p not passed")

    def test_x_ph(self):
        error = abs(self.steam_table.x_ph(1.0, 1000.0) - 0.258055424)
        self.assertLess(error, self.max_error, "test of x_ph not passed")

    def test_x_ps(self):
        error = abs(self.steam_table.x_ps(1.0, 4.0) - 0.445397961)
        self.assertLess(error, self.max_error, "test of x_ps not passed")

    def test_vx_ph(self):
        error = abs(self.steam_table.vx_ph(1.0, 418.0) - 0.288493093)
        self.assertLess(error, self.max_error, "test of vx_ph not passed")

    def test_vx_ps(self):
        error = abs(self.steam_table.vx_ps(1.0, 4.0) - 0.999233827)
        self.assertLess(error, self.max_error, "test of vx_ps not passed")
