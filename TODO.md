Open tasks and ideas
===================

Use this file as the central index for open work. Source comments are kept
in place; references below point back to their locations.

Table data access
-----------------

Decide how to make tabular data in pyXSteam.tables easier to access:

* Include selected, readable tables in the Sphinx documentation, or
* Provide a script that generates CSV files from the data classes on demand.

Choose a direction before implementing either approach. Keep the data classes
as the source of truth to avoid divergent values.

Implementation and behavior
----------------------------

* Confirm whether defaulting to Region 1 is correct for this region-selection
  case. (pyXSteam/RegionSelection.py:170)
* Investigate the commented-out mixture formulas for Cv_ps and w_ps in
  Region 4 and decide whether to implement them or retain the NaN behavior.
  (pyXSteam/XSteam.py:1665, 1891)
* Decide whether to use or remove the unused hB2bc_p function.
  (pyXSteam/RegionBorders.py:100)

Table and calculation checks
----------------------------

* Resolve the provenance of the five Region 5 residual Gibbs coefficients.
  They appear in the XSteam MATLAB implementation, whose comment cites R7-97
  generally, but do not match the six rows of official R7-97 Table 38. Find an
  older/alternate source or establish whether they are erroneous before
  changing them. (pyXSteam/tables/R7_97.py:1038,
  pyXSteam/Regions.py:1803; [XSteam MATLAB source](https://github.com/KurtJacobson/XSteam/blob/c8c168883e561cfa249a6ebe3d1c4dfb7c85f61c/MATLAB/XSteam.m#L2606),
  [official IAPWS R7-97](https://iapws.org/technical-guidance/release/IF97-Rev))

Tests and reference data
------------------------

* Reconcile the Region 3 T(p,s) and v(p,s) test inputs with SR3-03 Tables 12
  and 15; the current inputs differ from the values noted in the comments.
  (tests/test_Regions.py:285, 301)
* Identify the source of the Region 3 p(T) verification values.
  (tests/test_Regions.py:353)
* Find the correct table number for the Region 4 T(p) verification values.
  (tests/test_Regions.py:388)
* Split the Region 4 h(s) test data by its source tables.
  (tests/test_Regions.py:401)
* Clarify and resolve the FIXME in the Region 4 T(h,s) verification test.
  (tests/test_Regions.py:428)
* Resolve the Region 5 p(T) test mismatch with R7-97 Table 42; the test
  currently uses modified pressures and expected values.
  (tests/test_Regions.py:450)
* Identify or verify the source for the Region 5 T(p,h) and T(p,s) test
  values. (tests/test_Regions.py:485, 498)

Python 3.11 optimization backlog
----------------------------------

Upgrade the project floor to Python 3.11 and evaluate the high-value wins
available in the newer runtime before doing broader refactors.

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
* Continue reviewing the region selectors for other repeated calculations;
  make changes only when profiling shows meaningful cost and tests preserve the
  boundary behavior.
* Review opportunities to use `functools.cache`/`lru_cache` for memoized table
  lookups or repeated unit-conversion paths.
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
