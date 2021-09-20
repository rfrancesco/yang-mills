# Fork of claudio-bonati/yang-mills
## Current project: Extending yang\_mills\_tracedef to 2 (or higher) compactified dimensions

This fork is made specifically for the SU(3) gauge theory, with 0, 1 or 2 compactified dimensions.
To make it work for arbitrary SU(N) or SO(N) gauge groups, we should
- Find a better way to represent the trace-deformation coefficient matrix
- Modify the update routine, in order to calculate all the mixed terms

## Changelog

- Introduced a new parameter (gparam.c/h), *param.tracedef\_dim*, to store the number of dimensions included in the trace deformation.
- Renamed *param.htracedef* to *param.htracedef_d* and added *param.htracedef_mixed* to store the coefficients for the case (tracedefdim == 1, 2 && N\_c == 3)
- Introduced a new indicization of the lattice sites, (siorth, par, axis), where sites orthogonal to a given axis are indicized with lexicographical even/odd ordering. This generalizes the (sisp, t) labeling, is equal to it for axis=0, and ensures that for every orthogonal slice, consecutive siorth indices do not refer to neighboring lattice sites. Therefore, within a single orthogonal slice, metropolis computations can be parallelized by splitting the for cycle in 0-\> orth\_volume/2; orth\_volume/2-\>orth\_volume, as in the original code by _claudio-bonati_.
- Wrote a test script for this new indicization: _tests/test\_lexicographic\_orth.c_. It automatically checks for consistency (i.e. f(x) <-> f^-1(y), and sisp = siorth when axis=0). It prints to stdout the mapping (siorth, par, t) <-> (cartesian coord.), which can be checked visually.
- Generalized Polyakov calculation functions, and tracedef staples functions, to an arbitrary (trace-deformed) axis.
- Wrote functions for Mixed Polyakov loops calculation (only ~PiPj and PiPj^\dag for now).
- Generalized the Metropolis with tracedef to the case of tracedef\_dim == 0, 1, 2, N\_c = 3.  Sadly, this does not scale to higher N\_c/ higher tracedef\_dim in an elegant way.
- Wrote (lib/gauge_conf_upd.c) _polyakov\_loop_ (Polyakov loop from r to r as a matrix) and _linear\_parallel\_transport_ (Product of links from r, counting steps) to support the Metropolis algorithm
