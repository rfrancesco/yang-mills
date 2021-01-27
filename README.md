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
