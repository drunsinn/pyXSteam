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

Legacy roadmap
--------------

Review the old v0.5.0 TODO list, confirm which items are complete, and migrate
only the remaining work here. The list includes: minimum Python version,
function type hints and docstrings, enumerations for magic numbers, import
restructuring, pyproject.toml, demo-script location, dependency groups,
switching to Ruff, the uv workflow, a Makefile, and moving constants into the
tables module. (NEWS.txt:1-13)
