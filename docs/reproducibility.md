# Reproducibility

## Compiler

The source file is standard C++17.

Recommended build:

```bash
g++ -O3 -std=c++17 src/goldbach_lab_v03.cpp -o goldbach_lab_v03
```

## Reproducing the included result

Run:

```bash
./goldbach_lab_v03 --mode first-interior --maxN 10000000000 \
  --backend dense --output-prefix goldbach_to_10B_v03 \
  --progress-every 10000000
```

Expected output files:

```text
goldbach_to_10B_v03_report.txt
goldbach_to_10B_v03_first_interior_records.csv
goldbach_to_10B_v03_failures.csv
```

The exact running time depends on hardware. The included report recorded approximately:

```text
sieve/build: 85.2197 seconds
scan: 11538.8 seconds
total: 11624 seconds
```

## Validation summary

The included report records:

```text
Even N tested: 4999999999
Interior positive: 4999999998
Boundary-only: 1
No representation found: 0
```

The unique boundary-only case is `N = 4`.
