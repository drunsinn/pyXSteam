# -*- coding: utf-8 -*-
"""
These Tests are taken form the original XSteam Matlab Script.
Some Errors are calculated with the help of numpy matrix functions.
"""

import unittest
from pyXSteam.XSteam import XSteam


class MKS_FunctionTester(unittest.TestCase):

    def setUp(self):
        self.maxError = 1E-6
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)

    def tearDown(self):
        pass

    def test_tsat_p(self):
        error = self.steamTable.tsat_p(1.0) - 99.60591861
        self.assertLess(error, self.maxError, 'tsat_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_t_ph(self):
        error = self.steamTable.t_ph(1.0, 100.0) - 23.84481908
        self.assertLess(error, self.maxError, 't_ph not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_t_ps(self):
        error = self.steamTable.t_ps(1.0, 1.0) - 73.70859421
        self.assertLess(error, self.maxError, 't_ps not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_t_hs(self):
        error = self.steamTable.t_hs(100.0, 0.2) - 13.84933511
        self.assertLess(error, self.maxError, 't_hs not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_psat_t(self):
        error = self.steamTable.psat_t(100.0) - 1.014179779
        self.assertLess(error, self.maxError, 'psat_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_p_hs(self):
        error = self.steamTable.p_hs(84.0, 0.296) - 2.295498269
        self.assertLess(error, self.maxError, 'p_hs not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_hV_p(self):
        error = self.steamTable.hV_p(1.0) - 2674.949641
        self.assertLess(error, self.maxError, 'hV_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_hL_p(self):
        error = self.steamTable.hL_p(1.0) - 417.4364858
        self.assertLess(error, self.maxError, 'hL_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_hV_t(self):
        error = self.steamTable.hV_t(100.0) - 2675.572029
        self.assertLess(error, self.maxError, 'hV_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_hL_t(self):
        error = self.steamTable.hL_t(100.0) - 419.099155
        self.assertLess(error, self.maxError, 'hL_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_h_pt(self):
        error = self.steamTable.h_pt(1.0, 20.0) - 84.01181117
        self.assertLess(error, self.maxError, 'h_pt not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_h_ps(self):
        error = self.steamTable.h_ps(1.0, 1.0) - 308.6107171
        self.assertLess(error, self.maxError, 'h_ps not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_h_px(self):
        error = self.steamTable.h_px(1.0, 0.5) - 1546.193063
        self.assertLess(error, self.maxError, 'h_px not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_h_prho(self):
        error = self.steamTable.h_prho(1.0, 2.0) - 1082.773391
        self.assertLess(error, self.maxError, 'h_prho not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_h_tx(self):
        error = self.steamTable.h_tx(100.0, 0.5) - 1547.33559211
        self.assertLess(error, self.maxError, 'h_tx not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_vV_p(self):
        error = self.steamTable.vV_p(1.0) - 1.694022523
        self.assertLess(error, self.maxError, 'vV_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_vL_p(self):
        error = self.steamTable.vL_p(1.0) - 0.001043148
        self.assertLess(error, self.maxError, 'vL_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_vV_t(self):
        error = self.steamTable.vV_t(100.0) - 1.671860601
        self.assertLess(error, self.maxError, 'vV_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_vL_t(self):
        error = self.steamTable.vL_t(100.0) - 0.001043455
        self.assertLess(error, self.maxError, 'vL_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_v_pt(self):
        error = self.steamTable.v_pt(1.0, 100.0) - 1.695959407
        self.assertLess(error, self.maxError, 'v_pt not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_v_ph(self):
        error = self.steamTable.v_ph(1.0, 1000.0) - 0.437925658
        self.assertLess(error, self.maxError, 'v_ph not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_v_ps(self):
        error = self.steamTable.v_ps(1.0, 5.0) - 1.03463539
        self.assertLess(error, self.maxError, 'v_ps not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_rhoV_p(self):
        error = self.steamTable.rhoV_p(1.0) - 0.590310924
        self.assertLess(error, self.maxError, 'rhoV_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_rhoL_p(self):
        error = self.steamTable.rhoL_p(1.0) - 958.6368897
        self.assertLess(error, self.maxError, 'rhoL_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_rhoV_t(self):
        error = self.steamTable.rhoV_t(100.0) - 0.598135993
        self.assertLess(error, self.maxError, 'rhoV_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_rhoL_t(self):
        error = self.steamTable.rhoL_t(100.0) - 958.3542773
        self.assertLess(error, self.maxError, 'rhoL_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_rho_pt(self):
        error = self.steamTable.rho_pt(1.0, 100.0) - 0.589636754
        self.assertLess(error, self.maxError, 'rho_pt not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_rho_ph(self):
        error = self.steamTable.rho_ph(1.0, 1000.0) - 2.283492601
        self.assertLess(error, self.maxError, 'rho_ph not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_rho_ps(self):
        error = self.steamTable.rho_ps(1.0, 1.0) - 975.6236788
        self.assertLess(error, self.maxError, 'rho_ps not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_sV_p(self):
        error = self.steamTable.sV_p(0.006117) - 9.155465556
        self.assertLess(error, self.maxError, 'sV_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_sL_p(self):
        error = self.steamTable.sL_p(0.0061171) - 1.8359e-05
        self.assertLess(error, self.maxError, 'sL_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_sV_t(self):
        error = self.steamTable.sV_t(0.0001) - 9.155756716
        self.assertLess(error, self.maxError, 'sV_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_sL_t(self):
        error = self.steamTable.sL_t(100.0) - 1.307014328
        self.assertLess(error, self.maxError, 'sL_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_s_pt(self):
        error = self.steamTable.s_pt(1.0, 20.0) - 0.296482921
        self.assertLess(error, self.maxError, 's_pt not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_s_ph(self):
        error = self.steamTable.s_ph(1.0, 84.01181117) - 0.296813845
        self.assertLess(error, self.maxError, 's_ph not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_uV_p(self):
        error = self.steamTable.uV_p(1.0) - 2505.547389
        self.assertLess(error, self.maxError, 'uV_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_uL_p(self):
        error = self.steamTable.uL_p(1.0) - 417.332171
        self.assertLess(error, self.maxError, 'uL_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_uV_t(self):
        error = self.steamTable.uV_t(100.0) - 2506.015308
        self.assertLess(error, self.maxError, 'uV_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_uL_t(self):
        error = self.steamTable.uL_t(100.0) - 418.9933299
        self.assertLess(error, self.maxError, 'uL_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_u_pt(self):
        error = self.steamTable.u_pt(1.0, 100.0) - 2506.171426
        self.assertLess(error, self.maxError, 'u_pt not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_u_ph(self):
        error = self.steamTable.u_ph(1.0, 1000.0) - 956.2074342
        self.assertLess(error, self.maxError, 'u_ph not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_u_ps(self):
        error = self.steamTable.u_ps(1.0, 1.0) - 308.5082185
        self.assertLess(error, self.maxError, 'u_ps not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_CpV_p(self):
        error = self.steamTable.CpV_p(1.0) - 2.075938025
        self.assertLess(error, self.maxError, 'cpV_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_CpL_p(self):
        error = self.steamTable.CpL_p(1.0) - 4.216149431
        self.assertLess(error, self.maxError, 'cpL_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_CpV_t(self):
        error = self.steamTable.CpV_t(100.0) - 2.077491868
        self.assertLess(error, self.maxError, 'cpV_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_CpL_t(self):
        error = self.steamTable.CpL_t(100.0) - 4.216645119
        self.assertLess(error, self.maxError, 'cpL_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_Cp_pt(self):
        error = self.steamTable.Cp_pt(1.0, 100.0) - 2.074108555
        self.assertLess(error, self.maxError, 'cp_pt not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_Cp_ph(self):
        error = self.steamTable.Cp_ph(1.0, 200.0) - 4.17913573169
        self.assertLess(error, self.maxError, 'Cp_ph not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_Cp_ps(self):
        error = self.steamTable.Cp_ps(1.0, 1.0) - 4.190607038
        self.assertLess(error, self.maxError, 'Cp_ps not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_CvV_p(self):
        error = self.steamTable.CvV_p(1.0) - 1.552696979
        self.assertLess(error, self.maxError, 'CvV_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_CvL_p(self):
        error = self.steamTable.CvL_p(1.0) - 3.769699683
        self.assertLess(error, self.maxError, 'CvL_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_CvV_t(self):
        error = self.steamTable.CvV_t(100.0) - 1.553698696
        self.assertLess(error, self.maxError, 'CvV_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_CvL_t(self):
        error = self.steamTable.CvL_t(100.0) - 3.76770022
        self.assertLess(error, self.maxError, 'CvL_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_Cv_pt(self):
        error = self.steamTable.Cv_pt(1.0, 100.0) - 1.551397249
        self.assertLess(error, self.maxError, 'Cv_pt not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_Cv_ph(self):
        error = self.steamTable.Cv_ph(1.0, 200.0) - 4.035176364
        self.assertLess(error, self.maxError, 'Cv_ph not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_Cv_ps(self):
        error = self.steamTable.Cv_ps(1.0, 1.0) - 3.902919468
        self.assertLess(error, self.maxError, 'Cv_ps not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_wV_p(self):
        error = self.steamTable.wV_p(1.0) - 472.0541571
        self.assertLess(error, self.maxError, 'wV_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_wL_p(self):
        error = self.steamTable.wL_p(1.0) - 1545.451948
        self.assertLess(error, self.maxError, 'wL_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_wV_t(self):
        error = self.steamTable.wV_t(100.0) - 472.2559492
        self.assertLess(error, self.maxError, 'wV_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_wL_t(self):
        error = self.steamTable.wL_t(100.0) - 1545.092249
        self.assertLess(error, self.maxError, 'wL_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_w_pt(self):
        error = self.steamTable.w_pt(1.0, 100.0) - 472.3375235
        self.assertLess(error, self.maxError, 'w_pt not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_w_ph(self):
        error = self.steamTable.w_ph(1.0, 200.0) - 1542.682475
        self.assertLess(error, self.maxError, 'w_ph not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_w_ps(self):  # TODO: Check values and calculation
        error = self.steamTable.w_ps(1.0, 1.0) - 1557.8585
        self.assertLess(error, self.maxError, 'w_ps not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_my_pt(self):
        error = self.steamTable.my_pt(1.0, 100.0) - 1.22704e-05
        self.assertLess(error, self.maxError, 'my_pt not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_my_ph(self):
        error = self.steamTable.my_ph(1.0, 100.0) - 0.000914003770302
        self.assertLess(error, self.maxError, 'my_ph not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_my_ps(self):
        error = self.steamTable.my_ps(1.0, 1.0) - 0.000384222
        self.assertLess(error, self.maxError, 'my_ps not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_tcL_p(self):
        error = self.steamTable.tcL_p(1.0) - 0.677593822
        self.assertLess(error, self.maxError, 'tcL_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_tcV_p(self):
        error = self.steamTable.tcV_p(1.0) - 0.024753668
        self.assertLess(error, self.maxError, 'tcV_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_tcL_t(self):
        error = self.steamTable.tcL_t(25.0) - 0.607458162
        self.assertLess(error, self.maxError, 'tcL_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_tcV_t(self):
        error = self.steamTable.tcV_t(25.0) - 0.018326723
        self.assertLess(error, self.maxError, 'tcV_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_tc_pt(self):
        error = self.steamTable.tc_pt(1.0, 25.0) - 0.607509806
        self.assertLess(error, self.maxError, 'tc_pt not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_tc_ph(self):
        error = self.steamTable.tc_ph(1.0, 100.0) - 0.605710062
        self.assertLess(error, self.maxError, 'tc_ph not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_tc_hs(self):
        error = self.steamTable.tc_hs(100.0, 0.34) - 0.606283124
        self.assertLess(error, self.maxError, 'tc_hs not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_st_t(self):
        error = self.steamTable.st_t(100.0) - 0.0589118685877
        self.assertLess(error, self.maxError, 'st_t not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_st_p(self):
        error = self.steamTable.st_p(1.0) - 0.058987784
        self.assertLess(error, self.maxError, 'st_p not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_x_ph(self):
        error = self.steamTable.x_ph(1.0, 1000.0) - 0.258055424
        self.assertLess(error, self.maxError, 'x_ph not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_x_ps(self):
        error = self.steamTable.x_ps(1.0, 4.0) - 0.445397961
        self.assertLess(error, self.maxError, 'x_ps not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_vx_ph(self):
        error = self.steamTable.vx_ph(1.0, 418.0) - 0.288493093
        self.assertLess(error, self.maxError, 'vx_ph not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})

    def test_vx_ps(self):
        error = self.steamTable.vx_ps(1.0, 4.0) - 0.999233827
        self.assertLess(error, self.maxError, 'vx_ps not passed Error %(error)e allowed: %(max)e' % {'error': error, 'max': self.maxError})
