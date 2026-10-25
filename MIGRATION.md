# Migration guide: `master` (0.4.x) to `modernisation` (1.x.x)

This guide summarizes developer-facing changes between the current `master`
branch and `modernisation`. It is intended for applications upgrading across
the planned major release. The ordinary `XSteam` property methods remain
available with their existing names; the most important migration work is
selecting and updating the unit system explicitly. Some calculations were also
corrected or gained values where earlier versions returned `NaN`; see
[Behavior changes](#behavior-changes).

## Required: choose a unit system and update inputs/outputs

The old default `BARE` system used SI base values such as MPa and kelvin
(`XSteam()` or `UnitConverter()` with no argument). The new `BARE` is an alias
of `MKS`, which uses bar and degrees Celsius. There is no longer a separate
MPa/kelvin mode.

Consequences:

- Existing calls that relied on the no-argument default may now interpret the
  same numeric arguments differently.
- Passing the old integer value `0` for `BARE` is no longer supported and raises
  `ValueError`.
- The old integer values `1` and `2` continue to identify MKS and FLS,
  respectively, but use the named constants or enum instead of bare integers.

Choose a unit system explicitly and convert values at application boundaries:

```python
from pyXSteam import UnitSystem, XSteam

steam = XSteam(UnitSystem.MKS)  # pressure: bar; temperature: °C
steam_fls = XSteam(UnitSystem.FLS)  # pressure: psi; temperature: °F
```

`XSteam.UNIT_SYSTEM_MKS`, `XSteam.UNIT_SYSTEM_FLS`, and the corresponding
`XSteam_HW` constants remain available. `XSteam.UNIT_SYSTEM_BARE` remains an
alias of MKS; it does not restore the old MPa/kelvin behavior. The former
`UnitConverter.__UNIT_SYSTEM_BARE__`, `__UNIT_SYSTEM_MKS__`, and
`__UNIT_SYSTEM_FLS__` class attributes were removed; use the `UnitSystem`
enum. Apply the same unit-system update to direct `UnitConverter` construction
and to `XSteam_HW`.

## Public `XSteam` API additions

No existing high-level `XSteam` property method was removed or renamed. The
following methods were added:

- IAPWS-95 density/temperature property methods: `p_rhoT`, `u_rhoT`,
  `s_rhoT`, `h_rhoT`, `cv_rhoT`, `cp_rhoT`, `w_rhoT`, `joulethomson_rhoT`,
  `delta_T_rhoT`, and `isen_t_p_rhoT`.
- IAPWS-95 helpers: `sec_virial_coefficients`,
  `third_virial_coefficients`, and `phase_equilibrium_condition`.
- `XSteam_HW.st_t` for heavy-water surface tension.
- In `Regions`, metastable Region 2 forward properties add `_meta` methods
  (`v2_pT_meta`, `h2_pT_meta`, `u2_pT_meta`, `s2_pT_meta`, `Cp2_pT_meta`,
  `w2_pT_meta`); Region 3 adds `psat3_h` and `psat3_s`.
- `IAPWS_R14` adds inverse melting-temperature functions for ice III, V, and
  VI; `IAPWS_R4` adds `myHW_rhoT`.

The IAPWS R12-08 viscosity is available as opt-in
`XSteam.R12_my_rhot(rho, t, industrial_application=True)`. It does not replace
the legacy viscosity formulation used by `my_pt`, `my_ph`, and `my_ps`.

```python
mu = steam.R12_my_rhot(rho_kg_m3, temperature_celsius)
mu_full = steam.R12_my_rhot(
    rho_kg_m3,
    temperature_celsius,
    industrial_application=False,
)
```

The wrapper converts both density and temperature according to the selected
unit system. The full option includes near-critical enhancement; the default
industrial option omits it.

## Enum and constant changes

Import and use the public enums rather than internal integer constants:

```python
from pyXSteam import IceType, UnitSystem

steam = XSteam(UnitSystem.MKS)
pressure = steam.pmelt_t(-13.15, IceType.Ih)
```

- Unit systems are represented by `Constants.UnitSystem`; ice hints use
  `Constants.IceType`; region selection uses `Constants.DiagramRegion`.
  `UnitSystem` and `IceType` are also exported from the top-level `pyXSteam`
  package.
- `XSteam.TYPE_ICE_*` and `XSteam.UNIT_SYSTEM_*` aliases remain available.
- Former internal constants such as `UnitConverter.__UNIT_SYSTEM_BARE__` and
  `IAPWS_R14.__TYPE_ICE_Ih__` have been replaced by the public enum names
  `UnitSystem.BARE` and `IceType.Ih`. Common constants formerly exposed only
  with double-underscore names (for example `Constants.__CRITICAL_TEMPERATURE__`)
  now use public names such as `Constants.CRITICAL_TEMPERATURE`. Update code
  that imported the double-underscore implementation details.
- Region selectors now return `DiagramRegion` (`IntEnum`) values rather than
  plain integers. They compare equal to their former numeric values; prefer
  comparisons against `DiagramRegion.R1` etc.
- The `pyXSteam.name` attribute (`name = "pyXSteam"`) has been removed from
  `pyXSteam/__init__.py`. Use `pyXSteam.__name__` or the literal `"pyXSteam"`
  instead. The `__version__`, `__author__`, `__license__`, and `__email__`
  attributes remain available.

## Behavior changes

These do not rename the high-level methods, but may affect numerical results
or cases that were previously unsupported:

- `criticalDensity()` now returns density converted as density (including FLS),
  rather than applying a pressure conversion to density.
- `h_ps` now calculates and returns a value for Regions 1–5, including the
  two-phase mixture, instead of returning `NaN` for some previously
  unimplemented region cases.
- `Cv_ps` and `w_ps` now return legacy-style quality-weighted estimates for
  Region 4 mixtures instead of `NaN`. These are estimates, especially for
  sound speed, rather than derivative-based equilibrium-mixture properties.
- The Region 1 lower boundary in `select_region_ps` now accepts valid states
  whose entropy is slightly negative; unclassified states now return NILL
  rather than falling through to Region 1.
- Region 3 backward-equation coefficients and Region 5 reference data have
  been reconciled with their IAPWS source values. Results affected by those
  previously incorrect coefficients or verification-data errors can change.

## Low-level module renames and moves

Applications using internal calculation modules directly should update imports
and call sites. The high-level `XSteam` methods are not affected by these
renames.

| Old symbol | New symbol |
|------------|------------|
| `RegionSelection.region_pT` | `RegionSelection.select_region_pT` |
| `RegionSelection.region_ph` | `RegionSelection.select_region_ph` |
| `RegionSelection.region_ps` | `RegionSelection.select_region_ps` |
| `RegionSelection.region_hs` | `RegionSelection.select_region_hs` |
| `RegionSelection.region_prho` | `RegionSelection.select_region_prho` |
| `RegionBorders.B23p_T` | `RegionBorders.pB23_T` |
| `RegionBorders.B23T_p` | `RegionBorders.TB23_p` |
| `RegionBorders.p3sat_h` | `Regions.Region3.psat3_h` |
| `RegionBorders.p3sat_s` | `Regions.Region3.psat3_s` |

Additional IAPWS supplementary equations and region helpers have been added.
Direct imports from `Regions`, `RegionBorders`, or `RegionSelection` are
lower-level interfaces and may require reviewing the updated reference
releases and function signatures.

## Packaging and runtime requirements

- Minimum supported Python is now **3.11**. Upgrade the interpreter before
  installing this release.
- Packaging metadata and build configuration moved from `setup.py` and
  `MANIFEST.in` to `pyproject.toml`. Use a modern PEP 517 installer/build
  frontend; downstream packaging scripts should not invoke `setup.py` directly.
- Core runtime dependencies are empty. Install demo extras only when needed:
  `pip install "pyXSteam[demo]"`.
- Demo scripts moved from `bin/` into `pyXSteam/scripts/` and are exposed as
  the `pyXSteamDemo` and `pyXSteamRankineDemo` console commands.

## Developer upgrade checklist

- [ ] Replace implicit/default or integer `0` unit-system selection with an
      explicit `UnitSystem.MKS` or `UnitSystem.FLS` choice. For applications
      that intentionally relied on the former MPa/K default, add an explicit
      conversion layer; there is no built-in MPa/K mode in the new release.
- [ ] Update pressure, temperature, density, and output-value conversions to
      match the selected unit system; do not assume BARE means MPa and K.
- [ ] Replace imports of double-underscore constants with `UnitSystem`,
      `IceType`, `DiagramRegion`, or the corresponding public constants;
      update `UnitConverter.__UNIT_SYSTEM_*__` references.
- [ ] Update low-level region-selection and boundary-function imports if used.
- [ ] Keep existing `my_pt`/`my_ph`/`my_ps` calls for legacy viscosity; use
      `R12_my_rhot` only when intentionally opting into IAPWS R12-08.
- [ ] Install/test under Python 3.11 or newer, and update build automation to
      use `pyproject.toml`.
- [ ] Run the application's tests with representative MKS and FLS states,
      including tests for out-of-range states and any direct low-level imports.
