Open tasks and ideas
===================

Use this file as the central index for open work. Source comments are kept
in place; references below point back to their locations.

Merge PR #13 and release v1.0.0
-------------------------------

Track the work required to merge the `modernisation` branch into `master` and
publish the 1.0.0 release. Status markers: `[ ]` open, `[~]` in progress,
`[x]` done, `[!]` blocked.

Blockers (must be resolved before the merge):

* [x] Resolve the merge conflicts with `master`. PR #13 reports
  `mergeable=false` / `dirty`. Three files conflict:
  `.github/workflows/python-publish.yml` (take the `master` version:
  `pip install build twine` plus `python -m build`; do not regress to
  `pip wheel -w dist .`), `.gitignore` (union of both; `master` added the
  misspelled `docs/souce/README.md`, the branch added
  `docs/source/README.md`), and `docs/source/README.md` (take the branch
  version). Resolving these also restores the `pull_request` CI trigger.
* [x] Commit the `1.0.0` version bump. `pyXSteam/__init__.py` is committed as
  `__version__ = "0.5.0"` while the working tree holds an uncommitted bump to
  `1.0.0`; the release tag must match the package metadata.
* [x] Finalize `NEWS.txt`. The heading is still the placeholder
  `## v1.x.x, xx.yyyy`; replace it with the real version and date.
* [x] Commit the remaining uncommitted local changes
  (`docs/source/README.md`, `pyXSteam/__init__.py`).

Merge and release:

* [ ] Merge PR #13 into `master`. The repository default branch is `master`,
  not `main`.
* [ ] Verify the CI matrix after the merge. The test workflow
  (`python-package.yml`, Python 3.11-3.14) has not run on this code since
  2025-10-10 because the conflicted PR suppressed `pull_request` runs.
* [ ] Tag `v1.0.0` and create the GitHub release. `python-publish.yml` uploads
  to PyPI on `release: created`; confirm both the sdist and the wheel are
  produced and uploaded.

Follow-up issues found while reviewing the PR:

* [x] Sync `docs/source/README.md` with the root `README.md`. The committed
  copy lacks the new `UnitSystem` example and the migration-guide link.
* [x] Fix the broken `MIGRATION.md` documentation link. The README links
  `./MIGRATION.md`, but Sphinx builds from `docs/source/`, where the file does
  not exist, producing a `myst` cross-reference warning. Also add
  `MIGRATION.md` to the sdist and the documentation toctree.
* [x] Untrack the generated `docs/source/README.md`. It is tracked but listed
  in `.gitignore`, and `docs/Makefile` overwrites it via
  `cp -f ../README.md ./source/`, leaving the tree dirty after every
  `make html`. Both generated copies are now ignored and produced on demand.
* [x] Fix the coverage threshold mismatch in the README. The README states
  coverage fails below 60%, but `pyproject.toml` sets `--cov-fail-under=70`.
* [x] Correct the `main`/`master` references in `MIGRATION.md` and the README.
* [x] Document the removal of the `pyXSteam.name` attribute in `MIGRATION.md`.
* [x] Fix the dead `BARE` branch in `UnitConverter.__str__`. `UnitSystem.BARE`
  is an alias of `MKS` (both value 1), so the branch is unreachable and
  `str(UnitConverter(UnitSystem.BARE))` returns `"MKS (...)"`.
* [x] Add `ruff` linting to CI. `flake8` was removed but `ruff` was never
  added, so CI runs no linter; rename the workflow from
  "Build and Lint with Flake8".
* [ ] Resolve the license inconsistency. The `pyXSteam/__init__.py` docstring
  says GPLv2, `__license__` says MIT, and `LICENSE.txt` is LGPLv2.
* [ ] Update the stale PR #13 description. It claims Python 3.7 is the likely
  minimum (it is 3.11) and one checklist item is still unchecked.

Table data access
-----------------

Decide how to make tabular data in pyXSteam.tables easier to access:

* Include selected, readable tables in the Sphinx documentation, or
* Provide a script that generates CSV files from the data classes on demand.

Choose a direction before implementing either approach. Keep the data classes
as the source of truth to avoid divergent values.

Implementation and behavior
----------------------------

* Complete: `select_region_ps` now explicitly accepts the Region 1 lower
  entropy boundary at 273.15 K, including the slightly negative entropy at
  triple-point pressure. Inputs below the pressure-specific boundary and any
  otherwise unclassified states return NILL instead of defaulting to Region 1.
  (pyXSteam/RegionSelection.py)
* Complete: `Cv_ps` and `w_ps` in Region 4 now return quality-weighted
  interpolations of their saturated liquid and vapour endpoint values, using
  entropy-derived vapour quality. These are estimates, not derivative-based
  equilibrium-mixture properties; both formulas preserve the legacy XSteam
  interpolation convention. Tests cover MKS and FLS at low and high pressures.
  (pyXSteam/XSteam.py)
* Decide whether to use or remove the unused hB2bc_p function.
  (pyXSteam/RegionBorders.py:100)

IAPWS R12-08 viscosity formulation
----------------------------------

The industrial Eq. (10) calculation passes the Table 4 verification points, but
the critical-enhancement path is incomplete and its Table 5 test is skipped.
Complete the work in this order:

* Complete: R12-08 Table 5 exposed a missing square on coefficient `alpha` in
  the IAPWS-95 Gaussian residual term's second density derivative. This
  affected the susceptibility used by Eq. (20); correcting the derivative
  brings the calculated correlation lengths into agreement with all six
  Table 5 values. Added a near-critical finite-difference regression test.
  (pyXSteam/IAPWS_R6.py; tests/test_R6.py; tests/test_R12.py)
* Complete: corrected `R12_Y` using R12-08 Eqs. (15)–(19). The correlation
  length is in nanometers, so `q_C` and `q_D` now use reciprocal nanometers;
  corrected the high-length `w` denominator and Eq. (16)'s final term. Added
  Table 5 checks of `Y`, µ₂, and full viscosity covering the low- and high-
  length branches. (pyXSteam/IAPWS_R12.py; pyXSteam/tables/R12_08.py;
  tests/test_R12.py)
* Complete: restored R12-08 Table 5 verification for correlation length ξ,
  enhancement factor µ₂, and total viscosity µ at all six states. The obsolete
  skipped `eq10` checks were replaced with assertions using the current
  `my_rhoT` API and tolerances matched to the published precision. Coverage
  exercises both `R12_Y` branches. (tests/test_R12.py)
* Complete: retained Table 4 as the industrial Eq. (10) verification. Added
  tests immediately below, at, and above both Eq. (13) temperature and density
  bounds; boundaries are excluded as specified. Verified industrial viscosity
  skips µ₂ and full viscosity applies it. (pyXSteam/IAPWS_R12.py;
  tests/test_R12.py)
* Complete: R12 functions return NaN for non-finite inputs, nonpositive
  temperature, or negative density; `my_rhoT` accepts zero density as the
  dilute-gas limit. Positive finite states outside R12-08's recommended range
  remain extrapolatable, as the release describes. Guarded critical helper
  paths against arithmetic/domain exceptions. Added invalid-input and finite
  extrapolation coverage. (pyXSteam/IAPWS_R12.py; tests/test_R12.py)
* Complete: kept `XSteam.R12_my_rhot` opt-in; legacy `my_pt`, `my_ph`, and
  `my_ps` remain on the existing viscosity formulation. Documented R12-08
  availability and selected-unit-system conventions. Fixed the wrapper to
  convert density to SI for FLS input and convert the result back; added MKS/FLS
  wrapper checks. R15 verification passes with the corrected shared
  susceptibility calculation. (pyXSteam/XSteam.py; README.md;
  tests/test_R12.py; tests/test_R15.py)

Table and calculation checks
----------------------------

Tests and reference data (completed)
------------------------------------

* Complete: reconciled Region 3 T(p,s) and v(p,s) tests with SR3-03 (2014)
  Tables 12 and 15. The uncommented inputs and expected outputs already
  matched the official tables; stale commented legacy values did not.
  Tests now consume the canonical `SR3_03` table arrays directly.
  (tests/test_Regions.py)
* Complete: the Region 3 h(p,T) fixture is copied from the legacy XSteam
  Scilab verification block (`h3_pT`, `IF97` values), not an official R7-97
  Table 33 test. `h3_pT` solves by iterating the supplementary backward
  equation `T3_ph`; at Table 33 enthalpies, `T3_ph` differs from the tabulated
  temperatures by about 0.01 K, explaining why the h(p,T) fixture does not
  match Table 33. The separate `test_rhoT_function` verifies the fundamental
  Region 3 equation against Table 33. Retain the h(p,T) case as a legacy
  regression, not an independent accuracy test.
  ([XSteam Scilab source](https://github.com/KurtJacobson/XSteam/blob/c8c168883e561cfa249a6ebe3d1c4dfb7c85f61c/Scilab/XSteam.sci#L4532-L4541);
  tests/test_Regions.py)
* Assess whether to replace the `h3_pT` inversion through approximate
  `T3_ph` with a direct numerical solve of the Region 3 fundamental equation.
  Keep behavior unchanged until accuracy, convergence, and compatibility are
  evaluated. (pyXSteam/Regions.py:1326)
* Complete: the Region 4 T(p) verification values (0.1, 1, and 10 MPa)
  are from R7-97 (2012), Table 36, "Saturation temperatures calculated from
  Eq. (31) for selected values of p." (tests/test_Regions.py)
* Split the Region 4 h(s) test data by its source tables.
  (tests/test_Regions.py:401)
* Complete: corrected the first Region 4 T(h,s) verification value to the
  SR4-04 (2014) Table 29 value (346.8475498 K); the former expected value had
  a digit transposition. All three Table 29 cases now pass at the common
  relative-error threshold. (tests/test_Regions.py)
* Complete: corrected the Region 5 Table 38 coefficient transcription errors
  (the third coefficient exponent and final temperature exponent), and fixed
  Table 42's second pressure and specific-volume exponents. Region 5 forward
  calculations now use the official six-term coefficient table; forward and
  inverse tests use the Table 42 reference states and values.
  ([IAPWS R7-97](https://iapws.org/technical-guidance/release/IF97-Rev);
  tests/test_Regions.py:450)

Python 3.11 optimization backlog (completed)
---------------------------------------------

The project now requires Python 3.11 or newer. The profiling, optimization
review, and targeted changes below are complete; revisit only if representative
application workloads identify new bottlenecks.

* Baseline profile completed on Python 3.11.17 with `cProfile`: 1,000
  repetitions of 17 public calculations across Regions 1–5 and saturation
  paths produced 257,009 calls in 1.412 profiled seconds. The largest
  exclusive costs were `Region2.s2_pT` (0.287 s / 8,000 calls),
  `Region3.T3_ph` (0.205 s / 16,000 calls), and `Region1.s1_pT`
  (0.189 s / 7,000 calls). The Region 3 `h3_pT` iterative calculation took
  0.361 s cumulatively over 1,000 calls; `select_region_ps` took 0.551 s
  cumulatively over 4,000 calls, including its nested region-boundary
  calculations. These are synthetic-workload baselines, not production
  benchmarks. A profile of the full test run was dominated by pytest-cov
  analysis/reporting rather than library calculations. No algorithms were
  changed; use representative application workloads before optimizing.
* Complete: reviewed the profiled region-selector paths and reused duplicate
  per-call saturation-temperature boundary calculations in the pressure/entropy
  and pressure/density selectors. Further optimization is deferred until
  representative workloads justify it.
* Complete: table lookups are direct reads from static coefficient arrays, and
  unit conversions are inexpensive scalar arithmetic, so neither needs caching.
  A Python 3.11.17 microbenchmark of 100,000 calls found an `lru_cache(512)`
  slowed unique-input calls versus direct calls for `Region2.s2_pT` (2.52 s vs
  2.19 s), `pB23_T` (0.076 s vs 0.031 s), and `UnitConverter.toSIunit_p`
  (0.077 s vs 0.031 s). Exact repeat hits were much faster, especially for
  `Region2.s2_pT` (0.011 s), but those results use a repeated-input synthetic
  workload. Avoid global caching until a representative application workload
  demonstrates enough repeated inputs to offset miss overhead and retained
  cache entries.
* Complete: `typing.Self` is not applicable because the public classes expose
  no fluent or self-returning methods. A `TypeAlias` would only rename existing
  scalar or fixed-length tuple types, so it would not materially clarify the
  current API. Added missing `None`, `float`, and `bool` return annotations to
  constructors, heavy-water accessors, unit-system setters, and enum validators;
  these changes do not alter runtime behavior.
* Complete: no current boundary calculation walks adjacent points, so
  `itertools.pairwise` would add indirection without simplifying the logic.
  Rewrote the viscosity polynomial loops in `TransportProperties.py` to iterate
  over zipped coefficient tuples instead of indexing seven parallel arrays;
  the summation order is preserved.
* Complete: retain the existing `if`/`elif` region dispatch. The `XSteam`
  methods have many region switches, but their branches perform distinct
  property-specific calculations, conversions, and unsupported-region handling;
  replacing them with `match`/`case` would change syntax without reducing
  duplication or improving the flow enough to justify a broad mechanical edit.
* Complete: confirmed the existing package exports (`XSteam`, `XSteam_HW`,
  `UnitConverter`, `UnitSystem`, and `IceType`) and core IAPWS/table imports.
  Kept those exports unchanged; demo scripts appropriately require the optional
  `demo` dependencies. Replaced the hand-maintained package list with setuptools
  package discovery for `pyXSteam*`, removed redundant setuptools-scm/version
  and readme configuration, and configured the existing `LICENSE.txt` via the
  modern `license-files` metadata field without inferring an SPDX identifier.
