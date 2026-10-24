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

* Audit the hot paths for pure-Python bottlenecks and profile them on Python 3.11
  before changing algorithms.
* Replace repeated region/index selection logic with Python 3.11-friendly typing,
  caching, and clearer control flow where it reduces repeated work.
* Review opportunities to use `functools.cache`/`lru_cache` for memoized table
  lookups or repeated unit-conversion paths.
* Evaluate whether `typing.Self`, `TypeAlias`, and stricter annotations can make
  the public API clearer without changing behavior.
* Consider using `itertools.pairwise` and tuple-based iteration for boundary and
  interpolation code where it improves readability and performance.
* Decide whether a small `match`/`case` refactor is worthwhile for the region
  selection and dispatch logic.
* Re-check `__init__.py` exports and module imports for Python 3.11 packaging
  simplifications and import hygiene.
