#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Content of the Tables from related documents

Sources:

* IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water
and Steam, September 1997

* Revised Supplementary Release on Backward Equations for Pressure as a 
Function of Enthalpy and Entropy p(h,s) for Regions 1 and 2 of the IAPWS 
Industrial Formulation 1997 for the Thermodynamic Properties of Water and 
Steam

* Revised Supplementary Release on Backward Equations p(h,s) for Region 3, 
Equations as a Function of h and s for the Region Boundaries, and an Equation 
Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997 for the 
Thermodynamic Properties of Water and Steam 


"""


# Region 1
class R1:
    Table2_I = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32]
    Table2_J = [
        -2,
        -1,
        0,
        1,
        2,
        3,
        4,
        5,
        -9,
        -7,
        -1,
        0,
        1,
        3,
        -3,
        0,
        1,
        3,
        17,
        -4,
        0,
        6,
        -5,
        -2,
        10,
        -8,
        -11,
        -6,
        -29,
        -31,
        -38,
        -39,
        -40,
        -41,
    ]
    Table2_n = [
        0.14632971213167,
        -0.84548187169114,
        -3.756360367204,
        3.3855169168385,
        -0.95791963387872,
        0.15772038513228,
        -0.016616417199501,
        8.1214629983568e-04,
        2.8319080123804e-04,
        -6.0706301565874e-04,
        -0.018990068218419,
        -0.032529748770505,
        -0.021841717175414,
        -5.283835796993e-05,
        -4.7184321073267e-04,
        -3.0001780793026e-04,
        4.7661393906987e-05,
        -4.4141845330846e-06,
        -7.2694996297594e-16,
        -3.1679644845054e-05,
        -2.8270797985312e-06,
        -8.5205128120103e-10,
        -2.2425281908e-06,
        -6.5171222895601e-07,
        -1.4341729937924e-13,
        -4.0516996860117e-07,
        -1.2734301741641e-09,
        -1.7424871230634e-10,
        -6.8762131295531e-19,
        1.4478307828521e-20,
        2.6335781662795e-23,
        -1.1947622640071e-23,
        1.8228094581404e-24,
        -9.3537087292458e-26,
    ]

    Table6_I = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6]
    Table6_J = [0, 1, 2, 6, 22, 32, 0, 1, 2, 3, 4, 10, 32, 10, 32, 10, 32, 32, 32, 32]
    Table6_n = [
        -238.72489924521,
        404.21188637945,
        113.49746881718,
        -5.8457616048039,
        -1.528548241314e-04,
        -1.0866707695377e-06,
        -13.391744872602,
        43.211039183559,
        -54.010067170506,
        30.535892203916,
        -6.5964749423638,
        9.3965400878363e-03,
        1.157364750534e-07,
        -2.5858641282073e-05,
        -4.0644363084799e-09,
        6.6456186191635e-08,
        8.0670734103027e-11,
        -9.3477771213947e-13,
        5.8265442020601e-15,
        -1.5020185953503e-17,
    ]

    Table8_I = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4]
    Table8_J = [0, 1, 2, 3, 11, 31, 0, 1, 2, 3, 12, 31, 0, 1, 2, 9, 31, 10, 32, 32]
    Table8_n = [
        174.78268058307,
        34.806930892873,
        6.5292584978455,
        0.33039981775489,
        -1.9281382923196e-07,
        -2.4909197244573e-23,
        -0.26107636489332,
        0.22592965981586,
        -0.064256463395226,
        7.8876289270526e-03,
        3.5672110607366e-10,
        1.7332496994895e-24,
        5.6608900654837e-04,
        -3.2635483139717e-04,
        4.4778286690632e-05,
        -5.1322156908507e-10,
        -4.2522657042207e-26,
        2.6400441360689e-13,
        7.8124600459723e-29,
        -3.0732199903668e-31,
    ]

    Sub_psh12_Table2_I = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 5]
    Sub_psh12_Table2_J = [0, 1, 2, 4, 5, 6, 8, 14, 0, 1, 4, 6, 0, 1, 10, 4, 1, 4, 0]
    Sub_psh12_Table2_n = [
        -0.691997014660582,
        -18.361254878756,
        -9.28332409297335,
        65.9639569909906,
        -16.2060388912024,
        450.620017338667,
        854.68067822417,
        6075.23214001162,
        32.6487682621856,
        -26.9408844582931,
        -319.9478483343,
        -928.35430704332,
        30.3634537455249,
        -65.0540422444146,
        -4309.9131651613,
        -747.512324096068,
        730.000345529245,
        1142.84032569021,
        -436.407041874559,
    ]


# Region2
class R2:
    Table10_J0 = [0, 1, -5, -4, -3, -2, -1, 2, 3]
    Table10_n0 = [
        -9.6927686500217,
        10.086655968018,
        -0.005608791128302,
        0.071452738081455,
        -0.40710498223928,
        1.4240819171444,
        -4.383951131945,
        -0.28408632460772,
        0.021268463753307,
    ]
    Table11_I = [
        1,
        1,
        1,
        1,
        1,
        2,
        2,
        2,
        2,
        2,
        3,
        3,
        3,
        3,
        3,
        4,
        4,
        4,
        5,
        6,
        6,
        6,
        7,
        7,
        7,
        8,
        8,
        9,
        10,
        10,
        10,
        16,
        16,
        18,
        20,
        20,
        20,
        21,
        22,
        23,
        24,
        24,
        24,
    ]
    Table11_J = [
        0,
        1,
        2,
        3,
        6,
        1,
        2,
        4,
        7,
        36,
        0,
        1,
        3,
        6,
        35,
        1,
        2,
        3,
        7,
        3,
        16,
        35,
        0,
        11,
        25,
        8,
        36,
        13,
        4,
        10,
        14,
        29,
        50,
        57,
        20,
        35,
        48,
        21,
        53,
        39,
        26,
        40,
        58,
    ]
    Table11_n = [
        -1.7731742473213e-03,
        -0.017834862292358,
        -0.045996013696365,
        -0.057581259083432,
        -0.05032527872793,
        -3.3032641670203e-05,
        -1.8948987516315e-04,
        -3.9392777243355e-03,
        -0.043797295650573,
        -2.6674547914087e-05,
        2.0481737692309e-08,
        4.3870667284435e-07,
        -3.227767723857e-05,
        -1.5033924542148e-03,
        -0.040668253562649,
        -7.8847309559367e-10,
        1.2790717852285e-08,
        4.8225372718507e-07,
        2.2922076337661e-06,
        -1.6714766451061e-11,
        -2.1171472321355e-03,
        -23.895741934104,
        -5.905956432427e-18,
        -1.2621808899101e-06,
        -0.038946842435739,
        1.1256211360459e-11,
        -8.2311340897998,
        1.9809712802088e-08,
        1.0406965210174e-19,
        -1.0234747095929e-13,
        -1.0018179379511e-09,
        -8.0882908646985e-11,
        0.10693031879409,
        -0.33662250574171,
        8.9185845355421e-25,
        3.0629316876232e-13,
        -4.2002467698208e-06,
        -5.9056029685639e-26,
        3.7826947613457e-06,
        -1.2768608934681e-15,
        7.3087610595061e-29,
        5.5414715350778e-17,
        -9.436970724121e-07,
    ]
    Table16_I = [1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5]
    Table16_J = [0, 2, 5, 11, 1, 7, 16, 4, 16, 7, 10, 9, 10]
    Table16_n = [
        -0.73362260186506e-2,
        -0.88223831943146e-1,
        -0.72334555213245e-1,
        -0.40813178534455e-2,
        0.20097803380207e-2,
        -0.53045921898642e-1,
        -0.76190409086970e-2,
        -0.63498037657313e-2,
        -0.86043093028588e-1,
        0.75321581522770e-2,
        -0.79238375446139e-2,
        -0.22888160778447e-3,
        -0.26456501482810e-2,
    ]
    Table20_I = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7]
    Table20_J = [0, 1, 2, 3, 7, 20, 0, 1, 2, 3, 7, 9, 11, 18, 44, 0, 2, 7, 36, 38, 40, 42, 44, 24, 44, 12, 32, 44, 32, 36, 42, 34, 44, 28]
    Table20_n = [
        1089.8952318288,
        849.51654495535,
        -107.81748091826,
        33.153654801263,
        -7.4232016790248,
        11.765048724356,
        1.844574935579,
        -4.1792700549624,
        6.2478196935812,
        -17.344563108114,
        -200.58176862096,
        271.96065473796,
        -455.11318285818,
        3091.9688604755,
        252266.40357872,
        -6.1707422868339e-03,
        -0.31078046629583,
        11.670873077107,
        128127984.04046,
        -985549096.23276,
        2822454697.3002,
        -3594897141.0703,
        1722734991.3197,
        -13551.334240775,
        12848734.66465,
        1.3865724283226,
        235988.32556514,
        -13105236.545054,
        7399.9835474766,
        -551966.9703006,
        3715408.5996233,
        19127.72923966,
        -415351.64835634,
        -62.459855192507,
    ]
    Table21_I = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 7, 7, 9, 9]
    Table21_J = [
        0,
        1,
        2,
        12,
        18,
        24,
        28,
        40,
        0,
        2,
        6,
        12,
        18,
        24,
        28,
        40,
        2,
        8,
        18,
        40,
        1,
        2,
        12,
        24,
        2,
        12,
        18,
        24,
        28,
        40,
        18,
        24,
        40,
        28,
        2,
        28,
        1,
        40,
    ]
    Table21_n = [
        1489.5041079516,
        743.07798314034,
        -97.708318797837,
        2.4742464705674,
        -0.63281320016026,
        1.1385952129658,
        -0.47811863648625,
        8.5208123431544e-03,
        0.93747147377932,
        3.3593118604916,
        3.3809355601454,
        0.16844539671904,
        0.73875745236695,
        -0.47128737436186,
        0.15020273139707,
        -0.002176411421975,
        -0.021810755324761,
        -0.10829784403677,
        -0.046333324635812,
        7.1280351959551e-05,
        1.1032831789999e-04,
        1.8955248387902e-04,
        3.0891541160537e-03,
        1.3555504554949e-03,
        2.8640237477456e-07,
        -1.0779857357512e-05,
        -7.6462712454814e-05,
        1.4052392818316e-05,
        -3.1083814331434e-05,
        -1.0302738212103e-06,
        2.821728163504e-07,
        1.2704902271945e-06,
        7.3803353468292e-08,
        -1.1030139238909e-08,
        -8.1456365207833e-14,
        -2.5180545682962e-11,
        -1.7565233969407e-18,
        8.6934156344163e-15,
    ]
    Table22_I = [-7, -7, -6, -6, -5, -5, -2, -2, -1, -1, 0, 0, 1, 1, 2, 6, 6, 6, 6, 6, 6, 6, 6]
    Table22_J = [0, 4, 0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 4, 8, 4, 0, 1, 4, 10, 12, 16, 20, 22]
    Table22_n = [
        -3236839855524.2,
        7326335090218.1,
        358250899454.47,
        -583401318515.9,
        -10783068217.47,
        20825544563.171,
        610747.83564516,
        859777.2253558,
        -25745.72360417,
        31081.088422714,
        1208.2315865936,
        482.19755109255,
        3.7966001272486,
        -10.842984880077,
        -0.04536417267666,
        1.4559115658698e-13,
        1.126159740723e-12,
        -1.7804982240686e-11,
        1.2324579690832e-07,
        -1.1606921130984e-06,
        2.7846367088554e-05,
        -5.9270038474176e-04,
        1.2918582991878e-03,
    ]
    Table25_I = [
        -1.5,
        -1.5,
        -1.5,
        -1.5,
        -1.5,
        -1.5,
        -1.25,
        -1.25,
        -1.25,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -0.75,
        -0.75,
        -0.5,
        -0.5,
        -0.5,
        -0.5,
        -0.25,
        -0.25,
        -0.25,
        -0.25,
        0.25,
        0.25,
        0.25,
        0.25,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.75,
        0.75,
        0.75,
        0.75,
        1,
        1,
        1.25,
        1.25,
        1.5,
        1.5,
    ]
    Table25_J = [
        -24,
        -23,
        -19,
        -13,
        -11,
        -10,
        -19,
        -15,
        -6,
        -26,
        -21,
        -17,
        -16,
        -9,
        -8,
        -15,
        -14,
        -26,
        -13,
        -9,
        -7,
        -27,
        -25,
        -11,
        -6,
        1,
        4,
        8,
        11,
        0,
        1,
        5,
        6,
        10,
        14,
        16,
        0,
        4,
        9,
        17,
        7,
        18,
        3,
        15,
        5,
        18,
    ]
    Table25_n = [
        -392359.83861984,
        515265.7382727,
        40482.443161048,
        -321.93790923902,
        96.961424218694,
        -22.867846371773,
        -449429.14124357,
        -5011.8336020166,
        0.35684463560015,
        44235.33584819,
        -13673.388811708,
        421632.60207864,
        22516.925837475,
        474.42144865646,
        -149.31130797647,
        -197811.26320452,
        -23554.39947076,
        -19070.616302076,
        55375.669883164,
        3829.3691437363,
        -603.91860580567,
        1936.3102620331,
        4266.064369861,
        -5978.0638872718,
        -704.01463926862,
        338.36784107553,
        20.862786635187,
        0.033834172656196,
        -4.3124428414893e-05,
        166.53791356412,
        -139.86292055898,
        -0.78849547999872,
        0.072132411753872,
        -5.9754839398283e-03,
        -1.2141358953904e-05,
        2.3227096733871e-07,
        -10.538463566194,
        2.0718925496502,
        -0.072193155260427,
        2.074988708112e-07,
        -0.018340657911379,
        2.9036272348696e-07,
        0.21037527893619,
        2.5681239729999e-04,
        -0.012799002933781,
        -8.2198102652018e-06,
    ]
    Table26_I = [
        -6,
        -6,
        -5,
        -5,
        -4,
        -4,
        -4,
        -3,
        -3,
        -3,
        -3,
        -2,
        -2,
        -2,
        -2,
        -1,
        -1,
        -1,
        -1,
        -1,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        1,
        1,
        1,
        1,
        1,
        2,
        2,
        2,
        3,
        3,
        3,
        4,
        4,
        5,
        5,
        5,
    ]
    Table26_J = [
        0,
        11,
        0,
        11,
        0,
        1,
        11,
        0,
        1,
        11,
        12,
        0,
        1,
        6,
        10,
        0,
        1,
        5,
        8,
        9,
        0,
        1,
        2,
        4,
        5,
        6,
        9,
        0,
        1,
        2,
        3,
        7,
        8,
        0,
        1,
        5,
        0,
        1,
        3,
        0,
        1,
        0,
        1,
        2,
    ]
    Table26_n = [
        316876.65083497,
        20.864175881858,
        -398593.99803599,
        -21.816058518877,
        223697.85194242,
        -2784.1703445817,
        9.920743607148,
        -75197.512299157,
        2970.8605951158,
        -3.4406878548526,
        0.38815564249115,
        17511.29508575,
        -1423.7112854449,
        1.0943803364167,
        0.89971619308495,
        -3375.9740098958,
        471.62885818355,
        -1.9188241993679,
        0.41078580492196,
        -0.33465378172097,
        1387.0034777505,
        -406.63326195838,
        41.72734715961,
        2.1932549434532,
        -1.0320050009077,
        0.35882943516703,
        5.2511453726066e-03,
        12.838916450705,
        -2.8642437219381,
        0.56912683664855,
        -0.099962954584931,
        -3.2632037778459e-03,
        2.3320922576723e-04,
        -0.1533480985745,
        0.029072288239902,
        3.7534702741167e-04,
        1.7296691702411e-03,
        -3.8556050844504e-04,
        -3.5017712292608e-05,
        -1.4566393631492e-05,
        5.6420857267269e-06,
        4.1286150074605e-08,
        -2.0684671118824e-08,
        1.6409393674725e-09,
    ]
    Table27_I = [-2, -2, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7]
    Table27_J = [0, 1, 0, 0, 1, 2, 3, 0, 1, 3, 4, 0, 1, 2, 0, 1, 5, 0, 1, 4, 0, 1, 2, 0, 1, 0, 1, 3, 4, 5]
    Table27_n = [
        909.68501005365,
        2404.566708842,
        -591.6232638713,
        541.45404128074,
        -270.98308411192,
        979.76525097926,
        -469.66772959435,
        14.399274604723,
        -19.104204230429,
        5.3299167111971,
        -21.252975375934,
        -0.3114733441376,
        0.60334840894623,
        -0.042764839702509,
        5.8185597255259e-03,
        -0.014597008284753,
        5.6631175631027e-03,
        -7.6155864584577e-05,
        2.2440342919332e-04,
        -1.2561095013413e-05,
        6.3323132660934e-07,
        -2.0541989675375e-06,
        3.6405370390082e-08,
        -2.9759897789215e-09,
        1.0136618529763e-08,
        5.9925719692351e-12,
        -2.0677870105164e-11,
        -2.0874278181886e-11,
        1.0162166825089e-10,
        -1.6429828281347e-10,
    ]
    Sub_psh12_Table6_I = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 4, 5, 5, 6, 7]
    Sub_psh12_Table6_J = [1, 3, 6, 16, 20, 22, 0, 1, 2, 3, 5, 6, 10, 16, 20, 22, 3, 16, 20, 0, 2, 3, 6, 16, 16, 3, 16, 3, 1]
    Sub_psh12_Table6_n = [
        -1.82575361923032e-02,
        -0.125229548799536,
        0.592290437320145,
        6.04769706185122,
        238.624965444474,
        -298.639090222922,
        0.051225081304075,
        -0.437266515606486,
        0.413336902999504,
        -5.16468254574773,
        -5.57014838445711,
        12.8555037824478,
        11.414410895329,
        -119.504225652714,
        -2847.7798596156,
        4317.57846408006,
        1.1289404080265,
        1974.09186206319,
        1516.12444706087,
        1.41324451421235e-02,
        0.585501282219601,
        -2.97258075863012,
        5.94567314847319,
        -6236.56565798905,
        9659.86235133332,
        6.81500934948134,
        -6332.07286824489,
        -5.5891922446576,
        4.00645798472063e-02,
    ]
    Sub_psh12_Table7_I = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8, 8, 12, 14]
    Sub_psh12_Table7_J = [0, 1, 2, 4, 8, 0, 1, 2, 3, 5, 12, 1, 6, 18, 0, 1, 7, 12, 1, 16, 1, 12, 1, 8, 18, 1, 16, 1, 3, 14, 18, 10, 16]
    Sub_psh12_Table7_n = [
        8.01496989929495e-02,
        -0.543862807146111,
        0.337455597421283,
        8.9055545115745,
        313.840736431485,
        0.797367065977789,
        -1.2161697355624,
        8.72803386937477,
        -16.9769781757602,
        -186.552827328416,
        95115.9274344237,
        -18.9168510120494,
        -4334.0703719484,
        543212633.012715,
        0.144793408386013,
        128.024559637516,
        -67230.9534071268,
        33697238.0095287,
        -586.63419676272,
        -22140322476.9889,
        1716.06668708389,
        -570817595.806302,
        -3121.09693178482,
        -2078413.8463301,
        3056059461577.86,
        3221.57004314333,
        326810259797.295,
        -1441.04158934487,
        410.694867802691,
        109077066873.024,
        -24796465425889.3,
        1888019068.65134,
        -123651009018773,
    ]
    Sub_psh12_Table8_I = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 6, 6, 10, 12, 16]
    Sub_psh12_Table8_J = [0, 1, 2, 3, 4, 8, 0, 2, 5, 8, 14, 2, 3, 7, 10, 18, 0, 5, 8, 16, 18, 18, 1, 4, 6, 14, 8, 18, 7, 7, 10]
    Sub_psh12_Table8_n = [
        0.112225607199012,
        -3.39005953606712,
        -32.0503911730094,
        -197.5973051049,
        -407.693861553446,
        13294.3775222331,
        1.70846839774007,
        37.3694198142245,
        3581.44365815434,
        423014.446424664,
        -751071025.760063,
        52.3446127607898,
        -228.351290812417,
        -960652.417056937,
        -80705929.2526074,
        1626980172256.69,
        0.772465073604171,
        46392.9973837746,
        -13731788.5134128,
        1704703926305.12,
        -25110462818730.8,
        31774883083552,
        53.8685623675312,
        -55308.9094625169,
        -1028615.22421405,
        2042494187562.34,
        273918446.626977,
        -2.63963146312685e15,
        -1078908541.08088,
        -29649262098.0124,
        -1.11754907323424e15,
    ]


# Region3
class R3:
    Table30_I = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11]
    Table30_J = [
        0,
        0,
        1,
        2,
        7,
        10,
        12,
        23,
        2,
        6,
        15,
        17,
        0,
        2,
        6,
        7,
        22,
        26,
        0,
        2,
        4,
        16,
        26,
        0,
        2,
        4,
        26,
        1,
        3,
        26,
        0,
        2,
        26,
        2,
        26,
        2,
        26,
        0,
        1,
        26,
    ]
    Table30_n = [
        1.0658070028513,
        -15.732845290239,
        20.944396974307,
        -7.6867707878716,
        2.6185947787954,
        -2.808078114862,
        1.2053369696517,
        -8.4566812812502e-03,
        -1.2654315477714,
        -1.1524407806681,
        0.88521043984318,
        -0.64207765181607,
        0.38493460186671,
        -0.85214708824206,
        4.8972281541877,
        -3.0502617256965,
        0.039420536879154,
        0.12558408424308,
        -0.2799932969871,
        1.389979956946,
        -2.018991502357,
        -8.2147637173963e-03,
        -0.47596035734923,
        0.0439840744735,
        -0.44476435428739,
        0.90572070719733,
        0.70522450087967,
        0.10770512626332,
        -0.32913623258954,
        -0.50871062041158,
        -0.022175400873096,
        0.094260751665092,
        0.16436278447961,
        -0.013503372241348,
        -0.014834345352472,
        5.7922953628084e-04,
        3.2308904703711e-03,
        8.0964802996215e-05,
        -1.6557679795037e-04,
        -4.4923899061815e-05,
    ]
    # TODO check table 3!
    Sub_psh3_Table3_I = [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 6, 7, 8, 10, 10, 14, 18, 20, 22, 22, 24, 28, 28, 32, 32]
    Sub_psh3_Table3_J = [0, 1, 5, 0, 3, 4, 8, 14, 6, 16, 0, 2, 3, 0, 1, 4, 5, 28, 28, 24, 1, 32, 36, 22, 28, 36, 16, 28, 36, 16, 36, 10, 28]
    Sub_psh3_Table3_n = [
        7.70889828326934,
        -26.0835009128688,
        267.416218930389,
        17.2221089496844,
        -293.54233214597,
        614.135601882478,
        -61056.2757725674,
        -65127225.1118219,
        73591.9313521937,
        -11664650591.4191,
        35.5267086434461,
        -596.144543825955,
        -475.842430145708,
        69.6781965359503,
        335.674250377312,
        25052.6809130882,
        146997.380630766,
        5.38069315091534e19,
        1.43619827291346e21,
        3.64985866165994e19,
        -2547.41561156775,
        2.40120197096563e27,
        -3.93847464679496e29,
        1.47073407024852e24,
        -4.26391250432059e31,
        1.94509340621077e38,
        6.66212132114896e23,
        7.06777016552858e33,
        1.75563621975576e41,
        1.08408607429124e28,
        7.30872705175151e43,
        1.5914584739887e24,
        3.77121605943324e40,
    ]
    # TODO check table 4!
    Sub_psh3_Table4_I = [
        -12,
        -12,
        -12,
        -12,
        -12,
        -10,
        -10,
        -10,
        -10,
        -8,
        -8,
        -6,
        -6,
        -6,
        -6,
        -5,
        -4,
        -4,
        -4,
        -3,
        -3,
        -3,
        -3,
        -2,
        -2,
        -1,
        0,
        2,
        2,
        5,
        6,
        8,
        10,
        14,
        14,
    ]
    Sub_psh3_Table4_J = [2, 10, 12, 14, 20, 2, 10, 14, 18, 2, 8, 2, 6, 7, 8, 10, 4, 5, 8, 1, 3, 5, 6, 0, 1, 0, 3, 0, 1, 0, 1, 1, 1, 3, 7]
    Sub_psh3_Table4_n = [
        1.25244360717979e-13,
        -1.26599322553713e-02,
        5.06878030140626,
        31.7847171154202,
        -391041.161399932,
        -9.75733406392044e-11,
        -18.6312419488279,
        510.973543414101,
        373847.005822362,
        2.99804024666572e-08,
        20.0544393820342,
        -4.98030487662829e-06,
        -10.230180636003,
        55.2819126990325,
        -206.211367510878,
        -7940.12232324823,
        7.82248472028153,
        -58.6544326902468,
        3550.73647696481,
        -1.15303107290162e-04,
        -1.75092403171802,
        257.98168774816,
        -727.048374179467,
        1.21644822609198e-04,
        3.93137871762692e-02,
        7.04181005909296e-03,
        -82.910820069811,
        -0.26517881813125,
        13.7531682453991,
        -52.2394090753046,
        2405.56298941048,
        -22736.1631268929,
        89074.6343932567,
        -23923456.5822486,
        5687958081.29714,
    ]
    Sub_psh3_Table9_I = [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 7, 8, 12, 12, 14, 14, 16, 20, 20, 22, 24, 28, 32, 32]
    Sub_psh3_Table9_J = [14, 36, 3, 16, 0, 5, 4, 36, 4, 16, 24, 18, 24, 1, 4, 2, 4, 1, 22, 10, 12, 28, 8, 3, 0, 6, 8]
    Sub_psh3_Table9_n = [
        0.332171191705237,
        6.11217706323496e-04,
        -8.82092478906822,
        -0.45562819254325,
        -2.63483840850452e-05,
        -22.3949661148062,
        -4.28398660164013,
        -0.616679338856916,
        -14.682303110404,
        284.523138727299,
        -113.398503195444,
        1156.71380760859,
        395.551267359325,
        -1.54891257229285,
        19.4486637751291,
        -3.57915139457043,
        -3.35369414148819,
        -0.66442679633246,
        32332.1885383934,
        3317.66744667084,
        -22350.1257931087,
        5739538.75852936,
        173.226193407919,
        -3.63968822121321e-02,
        8.34596332878346e-07,
        5.03611916682674,
        65.5444787064505,
    ]
    Sub_psh3_Table10_I = [0, 0, 0, 0, 2, 3, 4, 4, 5, 5, 6, 7, 7, 7, 10, 10, 10, 32, 32]
    Sub_psh3_Table10_J = [1, 4, 10, 16, 1, 36, 3, 16, 20, 36, 4, 2, 28, 32, 14, 32, 36, 0, 6]
    Sub_psh3_Table10_n = [
        0.822673364673336,
        0.181977213534479,
        -0.011200026031362,
        -7.46778287048033e-04,
        -0.179046263257381,
        4.24220110836657e-02,
        -0.341355823438768,
        -2.09881740853565,
        -8.22477343323596,
        -4.99684082076008,
        0.191413958471069,
        5.81062241093136e-02,
        -1655.05498701029,
        1588.70443421201,
        -85.0623535172818,
        -31771.4386511207,
        -94589.0406632871,
        -1.3927384708869e-06,
        0.63105253224098,
    ]
    Sub_psh3_Table16_I = [1, 1, 2, 2, 4, 4, 7, 8, 8, 10, 12, 12, 18, 20, 24, 28, 28, 28, 28, 28, 32, 32, 32, 32, 32, 36, 36, 36, 36, 36]
    Sub_psh3_Table16_J = [8, 24, 4, 32, 1, 2, 7, 5, 12, 1, 0, 7, 10, 12, 32, 8, 12, 20, 22, 24, 2, 7, 12, 14, 24, 10, 12, 20, 22, 28]
    Sub_psh3_Table16_n = [
        -524.581170928788,
        -9269472.18142218,
        -237.385107491666,
        21077015581.2776,
        -23.9494562010986,
        221.802480294197,
        -5104725.33393438,
        1249813.96109147,
        2000084369.96201,
        -815.158509791035,
        -157.612685637523,
        -11420042233.2791,
        6.62364680776872e15,
        -2.27622818296144e18,
        -1.71048081348406e31,
        6.60788766938091e15,
        1.66320055886021e22,
        -2.18003784381501e29,
        -7.87276140295618e29,
        1.51062329700346e31,
        7957321.70300541,
        1.31957647355347e15,
        -3.2509706829914e23,
        -4.18600611419248e25,
        2.97478906557467e34,
        -9.53588761745473e19,
        1.66957699620939e24,
        -1.75407764869978e32,
        3.47581490626396e34,
        -7.10971318427851e38,
    ]
    Sub_psh3_Table17_I = [0, 0, 0, 1, 1, 5, 6, 7, 8, 8, 12, 16, 22, 22, 24, 36]
    Sub_psh3_Table17_J = [0, 3, 4, 0, 12, 36, 12, 16, 2, 20, 32, 36, 2, 32, 7, 20]
    Sub_psh3_Table17_n = [
        1.04351280732769,
        -2.27807912708513,
        1.80535256723202,
        0.420440834792042,
        -105721.24483466,
        4.36911607493884e24,
        -328032702839.753,
        -6.7868676080427e15,
        7439.57464645363,
        -3.56896445355761e19,
        1.67590585186801e31,
        -3.55028625419105e37,
        396611982166.538,
        -4.14716268484468e40,
        3.59080103867382e18,
        -1.16994334851995e40,
    ]
    Sub_psh3_Table23_I = [0, 1, 1, 3, 5, 6]
    Sub_psh3_Table23_J = [0, -2, 2, -12, -4, -3]
    Sub_psh3_Table23_n = [
        0.913965547600543,
        -4.30944856041991e-05,
        60.3235694765419,
        1.17518273082168e-18,
        0.220000904781292,
        -69.0815545851641,
    ]
    Sub_psh3_Table25_I = [-12, -10, -8, -4, -3, -2, -2, -2, -2, 0, 1, 1, 1, 3, 3, 5, 6, 6, 8, 8, 8, 12, 12, 14, 14]
    Sub_psh3_Table25_J = [10, 8, 3, 4, 3, -6, 2, 3, 4, 0, -3, -2, 10, -2, -1, -5, -6, -3, -8, -2, -1, -12, -1, -12, 1]
    Sub_psh3_Table25_n = [
        6.2909626082981e-04,
        -8.23453502583165e-04,
        5.15446951519474e-08,
        -1.17565945784945,
        3.48519684726192,
        -5.07837382408313e-12,
        -2.84637670005479,
        -2.36092263939673,
        6.01492324973779,
        1.48039650824546,
        3.60075182221907e-04,
        -1.26700045009952e-02,
        -1221843.32521413,
        0.149276502463272,
        0.698733471798484,
        -2.52207040114321e-02,
        1.47151930985213e-02,
        -1.08618917681849,
        -9.36875039816322e-04,
        81.9877897570217,
        -182.041861521835,
        2.61907376402688e-06,
        -29162.6417025961,
        1.40660774926165e-05,
        7832370.62349385,
    ]
    Sub_psh3_Table28_I = [
        0,
        0,
        0,
        1,
        1,
        1,
        1,
        2,
        2,
        2,
        3,
        3,
        3,
        3,
        4,
        4,
        5,
        5,
        5,
        5,
        6,
        6,
        6,
        8,
        10,
        10,
        12,
        14,
        14,
        16,
        16,
        18,
        18,
        18,
        20,
        28,
    ]
    Sub_psh3_Table28_J = [
        0,
        3,
        12,
        0,
        1,
        2,
        5,
        0,
        5,
        8,
        0,
        2,
        3,
        4,
        0,
        1,
        1,
        2,
        4,
        16,
        6,
        8,
        22,
        1,
        20,
        36,
        24,
        1,
        28,
        12,
        32,
        14,
        22,
        36,
        24,
        36,
    ]
    Sub_psh3_Table28_n = [
        0.179882673606601,
        -0.267507455199603,
        1.162767226126,
        0.147545428713616,
        -0.512871635973248,
        0.421333567697984,
        0.56374952218987,
        0.429274443819153,
        -3.3570455214214,
        10.8890916499278,
        -0.248483390456012,
        0.30415322190639,
        -0.494819763939905,
        1.07551674933261,
        7.33888415457688e-02,
        1.40170545411085e-02,
        -0.106110975998808,
        1.68324361811875e-02,
        1.25028363714877,
        1013.16840309509,
        -1.51791558000712,
        52.4277865990866,
        23049.5545563912,
        2.49459806365456e-02,
        2107964.67412137,
        366836848.613065,
        -144814105.365163,
        -1.7927637300359e-03,
        4899556021.00459,
        471.262212070518,
        -82929439019.8652,
        -1715.45662263191,
        3557776.82973575,
        586062760258.436,
        -12988763.5078195,
        31724744937.1057,
    ]


# Region4
class R4:
    Table37_J = [0, 1, -3, -2, -1, 2]
    Table34_n = [
        0.11670521452767e4,
        -0.72421316703206e6,
        -0.17073846940092e2,
        0.12020824702470e5,
        -0.32325550322333e7,
        0.14915108613530e2,
        -0.48232657361591e4,
        0.40511340542057e6,
        -0.23855557567849,
        0.65017534844798e3,
    ]

    Table37_J = [0, 1, -3, -2, -1, 2]
    Table37_n = [-13.179983674201, 6.8540841634434, -0.024805148933466, 0.36901534980333, -3.1161318213925, -0.32961626538917]


# Region5
class R5:
    Table38_I = [1, 1, 1, 2, 2, 3]
    Table38_J = [1, 2, 3, 3, 9, 8]
    Table38_n = [
        0.15736404855259e-2,
        0.90153761673944e-3,
        -0.50270077677648e-2,
        0.22440037409485e-5,
        -0.41163275453471e-5,
        0.37919454822955e-7,
    ]
