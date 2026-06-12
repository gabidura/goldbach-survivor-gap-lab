# Goldbach Survivor-Gap Lab

This repository accompanies the preprint:

**Goldbach Survivor-Gap Problem: Complete-Sieve Interior, Two-Cloud Survivor Masks, and Computational Validation**

The repository provides reproducible C++ code and selected output files for finite complete-sieve Goldbach survivor-window validation.

The computations do **not** prove the Goldbach conjecture. They validate the complete-sieve interior survivor formulation over tested finite ranges and record worst-case diagnostics for the survivor-gap program.

## Mathematical purpose

For an even integer `N`, the complete-sieve interior test sets

```text
z = floor(sqrt(N))
```

and searches for a prime pair

```text
x + (N - x) = N
```

with

```text
x > z and N - x > z.
```

Such a pair is a complete-sieve interior Goldbach representation: both coordinates survive divisibility by all primes up to `z`, and because both coordinates are greater than `z` and at most `N`, they must be prime.

## Main validation included

The included `goldbach_to_10B_v03_report.txt` records an exhaustive first-interior scan using the dense backend:

```text
Start N: 4
Stop N: 10000000000
Even N tested: 4999999999
Interior positive: 4999999998
Boundary-only: 1
No representation found: 0
```

The unique boundary-only case is `N = 4`, corresponding to `2 + 2`.

The scan verifies complete-sieve interior positivity for every even `4 < N <= 10^10`.

## Repository structure

```text
src/
  goldbach_lab_v03.cpp

scripts/
  build.sh
  run_first_interior_10B.sh
  run_gap_sample_10B.sh
  run_exact_gap_1M.sh
  run_single_example.sh

data/
  goldbach_to_10B_v03_report.txt
  goldbach_to_10B_v03_first_interior_records.csv
  goldbach_to_10B_v03_failures.csv
  README.md

paper/
  Paper_IV_Goldbach_Survivor_Gap_Problem_v04_10B_validation.pdf
  Paper_IV_Goldbach_Survivor_Gap_Problem_v04_10B_validation.tex

docs/
  method_notes.md
  reproducibility.md
  file_manifest.md
```

## Build

Linux/macOS/MinGW with a C++17 compiler:

```bash
g++ -O3 -std=c++17 src/goldbach_lab_v03.cpp -o goldbach_lab_v03
```

or use:

```bash
bash scripts/build.sh
```

## Reproduce the included 10B scan

The full 10B scan is computationally large. On a 32 GB i7-class machine, the dense backend uses about 625 MB for the prime bitset and may run for several hours.

```bash
./goldbach_lab_v03 --mode first-interior --maxN 10000000000 \
  --backend dense --output-prefix goldbach_to_10B_v03 \
  --progress-every 10000000
```

or use:

```bash
bash scripts/run_first_interior_10B.sh
```

## Other useful runs

Sampled complete-window gap scan:

```bash
./goldbach_lab_v03 --mode gap-sample --maxN 10000000000 \
  --step 1000000 --backend dense \
  --output-prefix goldbach_gap_sample_to_10B_v03 \
  --progress-every 1000
```

Exact complete-window gap scan over a smaller range:

```bash
./goldbach_lab_v03 --mode exact-gap --maxN 1000000 \
  --backend dense --output-prefix goldbach_exact_gap_to_1M_v03
```

Single large diagnostic using deterministic 64-bit Miller--Rabin:

```bash
./goldbach_lab_v03 --mode single --N 1000000000000 \
  --backend miller --output-prefix goldbach_single_1T_v03
```

## Caution

This repository is computational and diagnostic. It does not prove Goldbach. It tests and records the complete-sieve interior version of the Goldbach survivor-window problem over finite ranges.

## Related preprints

- Additive Encoding of Primality Constraints and Local Obstruction Structure in a Symmetric Affine System  
  https://zenodo.org/records/20259817

- Exact Finite Packet Sieves and Survivor Convolutions for Symmetric Diagonal Affine Systems  
  https://zenodo.org/records/20573926

- Global, Windowed, and Layered Survivor Envelopes: Primorial Survivor Spaces, Complete-Sieve Prime Interior, and Residual-Prime Attrition  
  https://zenodo.org/records/20617961

- Goldbach Survivor-Gap Problem: Complete-Sieve Interior, Two-Cloud Survivor Masks, and Computational Validation  
  https://zenodo.org/records/20660792
