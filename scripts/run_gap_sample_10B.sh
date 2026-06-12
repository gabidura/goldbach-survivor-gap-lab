#!/usr/bin/env bash
set -euo pipefail

./goldbach_lab_v03 --mode gap-sample --maxN 10000000000 \
  --step 1000000 --backend dense \
  --output-prefix goldbach_gap_sample_to_10B_v03 \
  --progress-every 1000
