#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Content of the Tables from related documents

Sources:

* IAPWS SR5-05(2016) Revised Supplementary Release on Backward Equations for Specific Volume
as a Function of Pressure and Temperature v(p,T)
for Region 3 of the IAPWS Industrial Formulation 1997 for the
Thermodynamic Properties of Water and Steam

"""


class SR5_05:
    """class with table data from IAPWS SR5-05(2016)"""

    # TODO

    # Table 1: Numerical values of the coefficients of the equations for subregion boundaries
    Table1_T3ab_I = [0, 1, 2, -1, -2]
    Table1_T3ab_n = [0.154793643129415e4, -0.817661219490113e3, 0.213144632222113e2, -0.191887498864292e4, 0.918419702359447e3]

    Table1_T3cd_I = [0, 1, 2, 3]
    Table1_T3cd_n = [0.585276966696349e3, 0.278233532206915e1, -0.127283549295878e-1, 0.159090746562729e-3]

    Table1_T3gh_I = [0, 1, 2, 3, 4]
    Table1_T3gh_n = [-0.249284240900418e5, 0.428143584791546e4, -0.269029173140130e3, 0.751608051114157e1, -0.787105249910383e-1]

    Table1_T3ij_I = [0, 1, 2, 3, 4]
    Table1_T3ij_n = [0.584814781649163e3, -0.616179320924617, 0.260763050899562, -0.587071076864459e-2, 0.515308185433082e-4]

    Table1_T3jk_I = [0, 1, 2, 3, 4]
    Table1_T3jk_n = [0.617229772068439e3, -0.770600270141675e1, 0.697972596851896, -0.157391839848015e-1, 0.137897492684194e-3]

    Table1_T3mn_I = [0, 1, 2, 3]
    Table1_T3mn_n = [0.535339483742384e3, 0.761978122720128e1, -0.158365725441648, 0.192871054508108e-2]

    Table1_T3op_I = [0, 1, 2, -1, -2]
    Table1_T3op_n = [0.969461372400213e3, -0.332500170441278, 0.642859598466067e2, 0.773845935768222e3, -0.152313732937084e4]

    Table1_T3qu_I = [0, 1, 2, 3]
    Table1_T3qu_n = [0.565603648239126e3, 0.529062258221222e1, -0.102020639611016, 0.122240301070145e-2]

    Table1_T3rx_I = [0, 1, 2, 3]
    Table1_T3rx_n = [0.584561202520006e3, -0.102961025163669e1, 0.243293362700452, -0.249905044740799e-2]
