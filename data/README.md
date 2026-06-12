# Data files

This directory contains selected output files from Goldbach Lab v03.

## Included 10B scan

- `goldbach_to_10B_v03_report.txt`  
  Text report for the exhaustive first-interior scan up to `10^10`.

- `goldbach_to_10B_v03_first_interior_records.csv`  
  Record cases for first interior prime coordinate, offset above `sqrt(N)`, and ratio `first_x / floor(sqrt(N))`.

- `goldbach_to_10B_v03_failures.csv`  
  Boundary/interior failure records. In the included scan, this contains only `N = 4`, the expected boundary-only case.

## Interpretation

The scan verifies complete-sieve interior positivity for all even `4 < N <= 10^10`. It does not prove the Goldbach conjecture.
