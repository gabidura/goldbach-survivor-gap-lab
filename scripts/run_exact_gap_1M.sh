#!/usr/bin/env bash
set -euo pipefail

./goldbach_lab_v03 --mode exact-gap --maxN 1000000 \
  --backend dense --output-prefix goldbach_exact_gap_to_1M_v03
