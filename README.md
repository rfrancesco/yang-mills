# Fork of claudio-bonati/yang-mills
## Current project: Extending yang\_mills\_tracedef to 2 (or higher) compactified dimensions
## TESTING Branch

This fork is made specifically for the SU(3) gauge theory, with 0, 1 or 2 compactified dimensions.
To make it work for arbitrary SU(N) or SO(N) gauge groups, we should
- Find a better way to represent the trace-deformation coefficient matrix
- Modify the update routine, in order to calculate all the mixed terms

- /!\ OpenMP Parallelization is probably broken /!\

This branch hosts, currently, a few commits that should fix the parallelization issue.
This is still in testing.

## Changelog

- Introduced a new parameter (gparam.c/h), *param.tracedef\_dim*, to store the number of dimensions included in the trace deformation.
- Renamed *param.htracedef* to *param.htracedef_d* and added *param.htracedef_mixed* to store the coefficients for the case (tracedefdim == 1, 2 && N\_c == 3)
- Introduced a new indicization of the lattice sites, (siorth, par, axis), where sites orthogonal to a given axis are indicized with lexicographical even/odd ordering. This generalizes the (sisp, t) labeling, is equal to it for axis=0, and ensures that for every orthogonal slice, consecutive siorth indices do not refer to neighboring lattice sites. Therefore, within a single orthogonal slice, metropolis computations can be parallelized by splitting the for cycle in 0-\> orth\_volume/2; orth\_volume/2-\>orth\_volume, as in the original code by _claudio-bonati_.
- Wrote a test script for this new indicization: _tests/test\_lexicographic\_orth.c_.
- Generalized Polyakov calculation functions, and tracedef staples functions, to an arbitrary (trace-deformed) axis.
- Wrote functions for Mixed Polyakov loops calculation
- Generalized the Metropolis with tracedef to the case of tracedef\_dim == 0, 1, 2, N\_c = 3.  Sadly, this does not scale elegantly to higher N\_c/ higher tracedef\_dim.

