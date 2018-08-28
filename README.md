# XSteam

Original Released by Magnus Holmgren for Matlab and Excel:
<http://xsteam.sourceforge.net> and/or <http://www.x-eng.com>

XSteam provides (mostly) accurate steam and water properties from 0 -
1000 bar and from 0 - 2000 °C according to the [IAPWS release IF-97](http://www.iapws.org/relguide/IF97-Rev.pdf). For
accuracy of the functions in different regions see IF-97 Page 4

Also includes thermal conductivity and viscosity, which are not part of
the IF97 release.
* Thermal Conductivity: (IAPWS 1998)
<http://www.iapws.org/relguide/ThCond.pdf>
* Viscosity: (2003)

## Requirements

There are no requirements for installing pyXSteam with Python 3.6 and up.
Tests require numpy, demos require numpy and matplotlib

## Install

run `python3 setup.py install`

To run unittests: `python setup.py test` but make sure numpy is installed

To test if setup was successful, run `python3 bin/pyXSteamDemo.py`. This will require numpy and matplotlib to be installed.
There are still (as of v0.4.0) some Errors in Thermal Conductivity and Speed of sound functions, so be warned that they exceed the Error Range. Apart form accuracy errors, there should be no warnings.


## Nomenclature

All Functions follow the same naming schema: First the wanted property,
then a underscore `_`, then the wanted input properties Example:
`t_ph` is temperature as a function of pressure and enthalpy. For a list
of valid functions se bellow:

| Property | Description                                                  |
|----------|--------------------------------------------------------------|
| t        | Temperature (°C or °F)                                       |
| p        | Pressure (bar or psi)                                        |
| h        | Enthalpy (kJ/kg or btu/lb)                                   |
| v        | Specific volume (m3/kg or ft\^3/lb)                          |
| rho      | Density (kg/m3 or lb/ft\^3)                                  |
| s        | Specific entropy (kJ/(kg °C) or btu/(lb °F))                 |
| u        | Specific internal energy (kJ/kg or btu/lb)                   |
| Cp       | Specific isobaric heat capacity (kJ/(kg °C) or btu/(lb °F))  |
| Cv       | Specific isochoric heat capacity (kJ/(kg °C) or btu/(lb °F)) |
| w        | Speed of sound (m/s or ft/s)                                 |
| my       | Viscosity (N s/m\^2 or lbm/ft/hr)                            |
| tc       | Thermal Conductivity (W/(m °C) or btu/(h ft °F))             |
| st       | Surface Tension (N/m or lb/ft)                               |
| x        | Vapor fraction                                               |
| vx       | Vapor Volume Fraction                                        |

## Usage

Simple Example:

    from pyXSteam.XSteam import XSteam
    steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS) print steamTable.hL_p(220.0)

By using the unitSystem Parameter, you can tell XSteam witch Unit System you are using.

    steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS) # m/kg/sec/°C/bar/W
    steamTable = XSteam(XSteam.UNIT_SYSTEM_FLS) # ft/lb/sec/°F/psi/btu
    steamTable = XSteam(XSteam.UNIT_SYSTEM_BARE) # m/kg/sec/K/MPa/W

To enable logging, add the following lines to your code:

    import logging logger = logging.getLogger(\'pyXSteam\')
    logger.setLevel(logging.DEBUG) sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter(\'%(name)s - %(levelname)s - %(message)s\'))
    logger.addHandler(sh)

## Available Functions

### Temperature
| Function | Description                                        |
|----------|----------------------------------------------------|
| tsat_p   | Saturation temperature                             |
| t_ph     | Temperature as a function of pressure and enthalpy |
| t_ps     | Temperature as a function of pressure and entropy  |
| t_hs     | Temperature as a function of enthalpy and entropy  |

### Pressure
| Function | Description                                                                                                    |
|----------|----------------------------------------------------------------------------------------------------------------|
| psat_t   | Saturation pressure                                                                                            |
| p_hs     | Pressure as a function of h and s.                                                                             |
| p_hrho   | Pressure as a function of h and rho. Very inaccurate for solid water region since it's almost incompressible! |

### Enthalpy
| Function | Description                                                                                                         |
|----------|---------------------------------------------------------------------------------------------------------------------|
| hV_p     | Saturated vapor enthalpy                                                                                            |
| hL_p     | Saturated liquid enthalpy                                                                                           |
| hV_t     | Saturated vapor enthalpy                                                                                            |
| hL_t     | Saturated liquid enthalpy                                                                                           |
| h_pt     | Enthalpy as a function of pressure and temperature                                                                  |
| h_ps     | Enthalpy as a function of pressure and entropy                                                                      |
| h_px     | Enthalpy as a function of pressure and vapor fraction                                                               |
| h_prho   | Enthalpy as a function of pressure and density. Observe for low temperatures (liquid) this equation has 2 solutions |
| h_tx     | Enthalpy as a function of temperature and vapor fraction                                                            |

### Specific volume
| Function | Description                                               |
|----------|-----------------------------------------------------------|
| vV_p     | Saturated vapor volume                                    |
| vL_p     | Saturated liquid volume                                   |
| vV_t     | Saturated vapor volume                                    |
| vL_t     | Saturated liquid volume                                   |
| v_pt     | Specific volume as a function of pressure and temperature |
| v_ph     | Specific volume as a function of pressure and enthalpy    |
| v_ps     | Specific volume as a function of pressure and entropy     |

### Density
| Function | Description                                       |
|----------|---------------------------------------------------|
| rhoV_p   | Saturated vapor density                           |
| rhoL_p   | Saturated liquid density                          |
| rhoV_t   | Saturated vapor density                           |
| rhoL_t   | Saturated liquid density                          |
| rho_pt   | Density as a function of pressure and temperature |
| rho_ph   | Density as a function of pressure and enthalpy    |
| rho_ps   | Density as a function of pressure and entropy     |

### Specific entropy
| Function | Description                                                                                              |
|----------|----------------------------------------------------------------------------------------------------------|
| sV_p     | Saturated vapor entropy                                                                                  |
| sL_p     | Saturated liquid entropy                                                                                 |
| sV_t     | Saturated vapor entropy                                                                                  |
| sL_t     | Saturated liquid entropy                                                                                 |
| s_pt     | Specific entropy as a function of pressure and temperature (Returns saturated vapor enthalpy if mixture) |
| s_ph     | Specific entropy as a function of pressure and enthalpy                                                  |

### Specific internal energy
| Function | Description                                                        |
|----------|--------------------------------------------------------------------|
| uV_p     | Saturated vapor internal energy                                    |
| uL_p     | Saturated liquid internal energy                                   |
| uV_t     | Saturated vapor internal energy                                    |
| uL_t     | Saturated liquid internal energy                                   |
| u_pt     | Specific internal energy as a function of pressure and temperature |
| u_ph     | Specific internal energy as a function of pressure and enthalpy    |
| u_ps     | Specific internal energy as a function of pressure and entropy     |

### Specific isobaric heat capacity
| Function | Description                                                               |
|----------|---------------------------------------------------------------------------|
| CpV_p    | Saturated vapor heat capacity                                             |
| CpL_p    | Saturated liquid heat capacity                                            |
| CpV_t    | Saturated vapor heat capacity                                             |
| CpL_t    | Saturated liquid heat capacity                                            |
| Cp_pt    | Specific isobaric heat capacity as a function of pressure and temperature |
| Cp_ph    | Specific isobaric heat capacity as a function of pressure and enthalpy    |
| Cp_ps    | Specific isobaric heat capacity as a function of pressure and entropy     |

### Specific isochoric heat capacity
| Function | Description                                                                |
|----------|----------------------------------------------------------------------------|
| CvV_p    | Saturated vapor isochoric heat capacity                                    |
| CvL_p    | Saturated liquid isochoric heat capacity                                   |
| CvV_t    | Saturated vapor isochoric heat capacity                                    |
| CvL_t    | Saturated liquid isochoric heat capacity                                   |
| Cv_pt    | Specific isochoric heat capacity as a function of pressure and temperature |
| Cv_ph    | Specific isochoric heat capacity as a function of pressure and enthalpy    |
| Cv_ps    | Specific isochoric heat capacity as a function of pressure and entropy     |

### Speed of sound
| Function | Description                                              |
|----------|----------------------------------------------------------|
| wV_p     | Saturated vapor speed of sound                           |
| wL_p     | Saturated liquid speed of sound                          |
| wV_t     | Saturated vapor speed of sound                           |
| wL_t     | Saturated liquid speed of sound                          |
| w_pt     | Speed of sound as a function of pressure and temperature |
| w_ph     | Speed of sound as a function of pressure and enthalpy    |
| w_ps     | Speed of sound as a function of pressure and entropy     |

### Viscosity
| Function | Description                                         |
|----------|-----------------------------------------------------|
| my_pt    | Viscosity as a function of pressure and temperature |
| my_ph    | Viscosity as a function of pressure and enthalpy    |
| my_ps    | Viscosity as a function of pressure and entropy     |

### Thermal Conductivity
| Function | Description                                                    |
|----------|----------------------------------------------------------------|
| tcL_p    | Saturated vapor thermal conductivity                           |
| tcV_p    | Saturated liquid thermal conductivity                          |
| tcL_t    | Saturated vapor thermal conductivity                           |
| tcV_t    | Saturated liquid thermal conductivity                          |
| tc_pt    | Thermal conductivity as a function of pressure and temperature |
| tc_ph    | Thermal conductivity as a function of pressure and enthalpy    |
| tc_hs    | Thermal conductivity as a function of enthalpy and entropy     |

### Surface tension
| Function | Description                                                  |
|----------|--------------------------------------------------------------|
| st_t     | Surface tension for two phase water/steam as a function of T |
| st_p     | Surface tension for two phase water/steam as a function of p |

### vapor fraction
| Function | Description                                           |
|----------|-------------------------------------------------------|
| x_ph     | vapor fraction as a function of pressure and enthalpy |
| x_ps     | vapor fraction as a function of pressure and entropy  |

## vapor volume fraction
| Function | Description                                                  |
|----------|--------------------------------------------------------------|
| vx_ph    | vapor volume fraction as a function of pressure and enthalpy |
| vx_ps    | vapor volume fraction as a function of pressure and entropy  |
