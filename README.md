# pyXSteam

pyXSteam is a Python library designed for calculating thermodynamic properties of water and steam based on the [IAPWS](http://www.iapws.org) (International Association for the Properties of Water and Steam) formulations.

This library was ported from the original Matlab released by Magnus Holmgren. The original can be founbd at [Sourceforge](http://xsteam.sourceforge.net) and/or <http://www.x-eng.com>.

A .NET implementation is available thanks to the work of rogerlew here: [XSteamNET](https://github.com/rogerlew/XSteamNET).

XSteam provides (mostly) accurate steam and water properties from 0 -
1000 bar and from 0 - 2000 °C according to the [IAPWS release IF-97](http://www.iapws.org/relguide/IF97-Rev.pdf). For
accuracy of the functions in different regions see IF-97 Page 4

Also included are functions for thermal conductivity and viscosity which are not part of R7-97(2012).

 - [IAPWS R15-11: Thermal Conductivity of Ordinary Water Substance](http://www.iapws.org/relguide/ThCond.pdf)
 - [IAPWS R12-08: Viscosity of Ordinary Water Substance](http://www.iapws.org/relguide/visc.pdf)
 - [IAPWS R1-76(2014): Surface Tension of Ordinary Water Substance](https://www.iapws.org/relguide/Surf-H2O-2014.pdf)

Some effort has been made to include the refined function of more recent releases and also functions for calculations on heavy water. This implementation is still in progress, these function should not be used. This includes:

 - IAPWS R4 (deprecated)
 - IAPWS R5-85(1994)
 - IAPWS R14


## Contributors to pyXSteam
In chronological order:
- IAPWS
- Magnus Holmgren
- drunsinn
- mohit0749
- gonmolina
- xjtu-blacksmith
- sebastiantuinstra

## Documentation

Check out [pyXSteam at readthedocs](https://pyxsteam.readthedocs.io)

## Requirements

There are no runtime requirements for installing pyXSteam with Python 3.11 and
up. Install the optional demo dependencies with `pip install "pyXSteam[demo]"`
(includes NumPy, Matplotlib, and SciPy).

### Tests and coverage

Tests require `pytest`, `pytest-cov`, and `numpy`. Run `pytest` (or `make test`)
to see uncovered lines and branches in the terminal and generate an HTML report
at `dist/reports/coverage.html`. The test command fails if total branch-aware
coverage falls below 60%. Demo scripts under `pyXSteam/scripts/` are excluded
from the coverage calculation because they are user-facing examples rather than
core library code.

## Install

run `pip install pyXSteam`

To run unittests you need two additional packages: `pytest` and `numpy`. After installing both, just run pytest in the base directory. The tests use some fixed values from the documentation to make sure that every coefficient and formula works as intended.

To test if setup was successful, install the demo dependencies and run
`pyXSteamDemo`.

## Usage

Simple Example:

```
    from pyXSteam import XSteam
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
    print(steam_table.hL_p(220.0))
```

The `unitSystem` parameter selects input and output units:

```
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)  # pressure: bar; temperature: °C
    steam_table = XSteam(XSteam.UNIT_SYSTEM_FLS)  # pressure: psi; temperature: °F
```

`UNIT_SYSTEM_BARE` is currently an alias of `UNIT_SYSTEM_MKS`, so it also uses
bar and °C rather than MPa and kelvin. Use `UNIT_SYSTEM_MKS` for the documented
MKS conversions.

pyXSteam uses Python's standard `logging` package. The library does not configure
the root logger or attach output handlers; applications should configure logging
once at startup. To display warnings and higher-severity messages:

```
    import logging
    logging.basicConfig(level=logging.WARNING)
```

To see pyXSteam's diagnostic messages, including expected out-of-range results
and region selection details, enable debug logging for the package:

```
    import logging
    logging.basicConfig(level=logging.WARNING)
    logging.getLogger("pyXSteam").setLevel(logging.DEBUG)
```

Package loggers follow the `pyXSteam.*` hierarchy and propagate to application
handlers. Out-of-range inputs that are handled by returning `NaN` are logged at
`DEBUG` to avoid flooding applications that calculate many values; iteration
failures remain warnings. Exceptions are raised without also logging them, so
applications can decide how to report or handle them.

## Related projects

Implementations of the IAPWS releases are also available for other programming languages:

- Matlab and Excel: [xsteam](http://xsteam.sourceforge.net)
- Matlab: [xsteam](https://github.com/stu314159/xsteam)
- .NET: [XSteamNET](https://github.com/rogerlew/XSteamNET)


## Notes

### Density (rho)

Density is calculated as 1/v. See section for Specific volume

### Viscosity

The legacy `my_pt`, `my_ph`, and `my_ps` methods use the revised 1985 IAPWS viscosity formulation. Their behavior is unchanged for compatibility.

The IAPWS R12-08 formulation is available separately as the opt-in `XSteam.R12_my_rhot(rho, t, industrial_application=True)` method. It takes density and temperature in the selected unit system and returns viscosity in that system: kg/m³, °C, and Pa·s for MKS; lb/ft³, °F, and lbm/(ft·hr) for FLS. Setting `industrial_application=False` includes the near-critical enhancement; the default industrial form omits it. This method is marked experimental, so validate its suitability for your application.

### Thermal conductivity

Revised release on the IAPS Formulation 1985 for the Thermal Conductivity of ordinary water substance (IAPWS 1998)

## Nomenclature

All Functions follow the same naming schema: First the wanted property, then a underscore `_`, then the wanted input properties Example: `t_ph` is temperature as a function of pressure and enthalpy. For a list of valid functions se below:

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
| p_hrho   | Pressure as a function of h and rho. Very inaccurate for solid water region since it's almost incompressible!  |
| pmelt_t  | Pressure as a function of temperature along the melting curve. Optional parameter to select ice region         |
| psubl_t  | Pressure as a function of temperature along the sublimation curve.                                             |


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

### vapor volume fraction
| Function | Description                                                  |
|----------|--------------------------------------------------------------|
| vx_ph    | vapor volume fraction as a function of pressure and enthalpy |
| vx_ps    | vapor volume fraction as a function of pressure and entropy  |

### Pressure along the Melting and Sublimation Curves
| Function | Description                                                        |
|----------|--------------------------------------------------------------------|
| pmelt_t  | Pressure along the melting curve as a function of temperature      |
| psubl_t  | Pressure along the sublimation curve as a function of temperature  |

## Development
- pull requests are always wellcome!
- code style is enforced by black
