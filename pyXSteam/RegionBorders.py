# -*- coding: utf-8 -*-
"""
Section 4: Region Borders
"""

"""Section 4.1 Boundary between region 2 and 3."""


def B23p_T(T):
    """
    function B23p_T = B23p_T(T)
    # %Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam
    # %1997
    # %Section 4 Auxiliary Equation for the Boundary between Regions 2 and 3
    # %Eq 5, Page 5
    """
    return 348.05185628969 - 1.1671859879975 * T + 1.0192970039326E-03 * (T ** 2)


def B23T_p(p):
    """function B23T_p = B23T_p(p)
    # %Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam
    # %1997
    # %Section 4 Auxiliary Equation for the Boundary between Regions 2 and 3
    # %Eq 6, Page 6
    """
    return 572.54459862746 + ((p - 13.91883977887) / 1.0192970039326E-03) ** 0.5

# Section 4.2 Region 3. pSat_h  & pSat_s


def p3sat_h(h):
    """
    function p3sat_h = p3sat_h(h)
    # %Revised Supplementary Release on Backward Equations for the functions T(p,h), v(p,h)  & T(p,s), v(p,s) for Region 3 of the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water  & Steam
    # %2004
    # %Section 4 Boundary Equations psat(h)  & psat(s) for the Saturation Lines of Region 3
    # %Se pictures Page 17, Eq 10, Table 17, Page 18
    """
    Ii = [0, 1, 1, 1, 1, 5, 7, 8, 14, 20, 22, 24, 28, 36]
    Ji = [0, 1, 3, 4, 36, 3, 0, 24, 16, 16, 3, 18, 8, 24]
    ni = [0.600073641753024, -9.36203654849857, 24.6590798594147, -107.014222858224, -91582131580576.8, -8623.32011700662, -23.5837344740032, 2.52304969384128E+17, -3.89718771997719E+18, -3.33775713645296E+22, 35649946963.6328, -1.48547544720641E+26, 3.30611514838798E+18, 8.13641294467829E+37]
    hs = h / 2600
    ps = 0
    for i in range(0, 14):
        ps = ps + ni[i] * (hs - 1.02) ** Ii[i] * (hs - 0.608) ** Ji[i]
    return ps * 22


def p3sat_s(s):
    """ function p3sat_s = p3sat_s(s)"""
    Ii = [0, 1, 1, 4, 12, 12, 16, 24, 28, 32]
    Ji = [0, 1, 32, 7, 4, 14, 36, 10, 0, 18]
    ni = [0.639767553612785, -12.9727445396014, -2.24595125848403E+15, 1774667.41801846, 7170793495.71538, -3.78829107169011E+17, -9.55586736431328E+34, 1.87269814676188E+23, 119254746466.473, 1.10649277244882E+36]
    Sigma = s / 5.2
    Pi = 0
    for i in range(0, 10):
        Pi = Pi + ni[i] * (Sigma - 1.03) ** Ii[i] * (Sigma - 0.699) ** Ji[i]
    return Pi * 22

# 4.3 Region boundary 1to3  & 3to2 as a functions of s


def hB13_s(s):
    """
    function hB13_s = hB13_s(s)
    # %Supplementary Release on Backward Equations ( ) , p h s for Region 3,
    # %Chapter 4.5 page 23.
    """
    Ii = [0, 1, 1, 3, 5, 6]
    Ji = [0, -2, 2, -12, -4, -3]
    ni = [0.913965547600543, -4.30944856041991E-05, 60.3235694765419, 1.17518273082168E-18, 0.220000904781292, -69.0815545851641]
    Sigma = s / 3.8
    eta = 0
    for i in range(0, 6):
        eta = eta + ni[i] * (Sigma - 0.884) ** Ii[i] * (Sigma - 0.864) ** Ji[i]
    return eta * 1700


def TB23_hs(h, s):
    """
    function TB23_hs = TB23_hs(h, s)
    # % Supplementary Release on Backward Equations () , p h s for Region 3,
    # % Chapter 4.6 page 25.
    """
    Ii = [-12, -10, -8, -4, -3, -2, -2, -2, -2, 0, 1, 1, 1, 3, 3, 5, 6, 6, 8, 8, 8, 12, 12, 14, 14]
    Ji = [10, 8, 3, 4, 3, -6, 2, 3, 4, 0, -3, -2, 10, -2, -1, -5, -6, -3, -8, -2, -1, -12, -1, -12, 1]
    ni = [6.2909626082981E-04, -8.23453502583165E-04, 5.15446951519474E-08, -1.17565945784945, 3.48519684726192, -5.07837382408313E-12, -2.84637670005479, -2.36092263939673, 6.01492324973779, 1.48039650824546, 3.60075182221907E-04, -1.26700045009952E-02, -1221843.32521413, 0.149276502463272, 0.698733471798484, -2.52207040114321E-02, 1.47151930985213E-02, -1.08618917681849, -9.36875039816322E-04, 81.9877897570217, -182.041861521835, 2.61907376402688E-06, -29162.6417025961, 1.40660774926165E-05, 7832370.62349385]
    Sigma = s / 5.3
    eta = h / 3000
    teta = 0
    for i in range(0, 25):
        teta = teta + ni[i] * (eta - 0.727) ** Ii[i] * (Sigma - 0.864) ** Ji[i]
    return teta * 900
