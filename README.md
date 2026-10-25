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
coverage falls below 70%. Demo scripts under `pyXSteam/scripts/` are excluded
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

```python
from pyXSteam import UnitSystem, XSteam

steam_table = XSteam(UnitSystem.MKS)  # pressure: bar; temperature: °C
steam_table_fls = XSteam(UnitSystem.FLS)  # pressure: psi; temperature: °F
```

For the full upgrade checklist when moving from `master` 0.4.x to the planned
1.0.x release, see the [migration guide](./MIGRATION.md). The key breaking
change is that `BARE` is now an alias for MKS (bar and °C), not the former
MPa-and-kelvin unit convention.

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

Most property methods follow `property_inputs`: the property symbol is followed
by underscore-separated independent variables. For example, `t_ph` returns
temperature from pressure and specific enthalpy, while `h_tx` returns specific
enthalpy from saturation temperature and vapor quality. Arguments use the
selected unit system described above; quality and Prandtl number are
unitless.

| Symbol | Meaning | MKS/BARE units | FLS units |
|--------|---------|----------------|-----------|
| `p` | Pressure | bar | psi |
| `t` | Temperature | °C | °F |
| `h` | Specific enthalpy | kJ/kg | Btu/lb |
| `u` | Specific internal energy | kJ/kg | Btu/lb |
| `s` | Specific entropy | kJ/(kg·K) | Btu/(lb·°F) |
| `v` | Specific volume | m³/kg | ft³/lb |
| `rho` | Density | kg/m³ | lb/ft³ |
| `x` | Vapor quality (mass fraction) | dimensionless, 0–1 | dimensionless, 0–1 |
| `vx` | Vapor volume fraction | dimensionless, 0–1 | dimensionless, 0–1 |
| `Cp` | Isobaric specific heat capacity | kJ/(kg·K) | Btu/(lb·°F) |
| `Cv` | Isochoric specific heat capacity | kJ/(kg·K) | Btu/(lb·°F) |
| `w` | Speed of sound | m/s | ft/s |
| `my` | Dynamic viscosity | Pa·s | lbm/(ft·hr) |
| `tc` | Thermal conductivity | W/(m·K) | Btu/(h·ft·°F) |
| `st` | Surface tension | N/m | lb/ft |
| `pr` | Prandtl number | dimensionless | dimensionless |
| `B`, `C` | Second and third virial coefficients | m³/kg; (m³/kg)² | ft³/lb; (ft³/lb)² |
| `dB/dT`, `dC/dT` | Temperature derivatives of virial coefficients | coefficient units/K | coefficient units/°F |
| `Joule-Thomson coefficient` | Temperature change per pressure change at constant enthalpy | K/MPa internally; K/bar as returned by MKS/BARE | °F/psi |
| `delta_T` | Isothermal throttling coefficient | Specific-volume units | Specific-volume units |

In saturation method names, uppercase `V` means saturated vapor and uppercase
`L` means saturated liquid (`hV_p`, `vL_t`). The final lowercase suffixes
`p`, `t`, `h`, and `s` identify pressure, temperature, enthalpy, and entropy
inputs. `rho` means density, `x` means vapor quality, and `hs` means the pair
of enthalpy and entropy inputs. The case-sensitive suffix `rhoT` identifies
IAPWS-95 methods taking density and temperature. The virial-coefficient
methods return a `(coefficient, temperature_derivative)` pair; their `rho`
argument is accepted for API consistency but the zero-density coefficients do
not depend on it.

`pmelt_t` also accepts the optional `hint` argument (`IceType`) to select an
ice melting-curve branch; `industrial_application` in `R12_my_rhot` selects
whether to include the near-critical viscosity enhancement. These options are
not encoded in the method-name suffixes. Methods such as `specificGasConstant`
and `criticalPressure` return constants and do not follow the property-input
naming pattern.


## Available Functions

### Temperature
| Function | Description                                        |
|----------|----------------------------------------------------|
| tsat_p   | Saturation temperature as a function of pressure   |
| tsat_s   | Saturation temperature as a function of entropy    |
| t_ph     | Temperature as a function of pressure and enthalpy |
| t_ps     | Temperature as a function of pressure and entropy  |
| t_hs     | Temperature as a function of enthalpy and entropy  |

### Pressure
| Function | Description                                                                                                    |
|----------|----------------------------------------------------------------------------------------------------------------|
| psat_t   | Saturation pressure as a function of temperature                                                               |
| psat_s   | Saturation pressure as a function of entropy                                                                   |
| p_hs     | Pressure as a function of enthalpy and entropy                                                                 |
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
| s_pt     | Specific entropy as a function of pressure and temperature                                            |
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
| tcL_p    | Saturated liquid thermal conductivity as a function of pressure |
| tcV_p    | Saturated vapor thermal conductivity as a function of pressure  |
| tcL_t    | Saturated liquid thermal conductivity as a function of temperature |
| tcV_t    | Saturated vapor thermal conductivity as a function of temperature |
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
| pmelt_t  | Pressure along the melting curve; optional ice-type hint selects branch |
| psubl_t  | Pressure along the sublimation curve as a function of temperature      |

### Prandtl number
| Function | Description                                           |
|----------|-------------------------------------------------------|
| pr_pt    | Prandtl number as a function of pressure and temperature |
| pr_ph    | Prandtl number as a function of pressure and enthalpy    |

### Opt-in IAPWS R12-08 viscosity
| Function | Description |
|----------|-------------|
| R12_my_rhot | Dynamic viscosity as a function of density and temperature; optionally includes near-critical enhancement |

### IAPWS-95 properties from density and temperature
| Function | Description |
|----------|-------------|
| p_rhoT | Pressure |
| u_rhoT | Specific internal energy |
| s_rhoT | Specific entropy |
| h_rhoT | Specific enthalpy |
| cv_rhoT | Isochoric specific heat capacity |
| cp_rhoT | Isobaric specific heat capacity |
| w_rhoT | Speed of sound |
| joulethomson_rhoT | Joule-Thomson coefficient |
| delta_T_rhoT | Isothermal throttling coefficient |
| isen_t_p_rhoT | Isentropic temperature-pressure coefficient |

### IAPWS-95 virial and phase-equilibrium utilities
| Function | Description |
|----------|-------------|
| sec_virial_coefficients | Second virial coefficient and temperature derivative |
| third_virial_coefficients | Third virial coefficient and temperature derivative |
| phase_equilibrium_condition | Pressure-equality and Maxwell residuals for liquid/vapor densities at a temperature |

### Constants
| Function | Description |
|----------|-------------|
| specificGasConstant | Specific gas constant R |
| criticalTemperatur | Critical temperature |
| criticalPressure | Critical pressure |
| criticalDensity | Critical density |
| triplePointTemperatur | Triple-point temperature |
| triplePointPressure | Triple-point pressure |
| zeroPointTemperature | Absolute-zero temperature |

## Development
- pull requests are always wellcome!
- code style is enforced by black
