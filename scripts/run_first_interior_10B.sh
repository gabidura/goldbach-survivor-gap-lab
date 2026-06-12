#!/usr/bin/env bash
set -euo pipefail

./goldbach_lab_v03 --mode first-interior --maxN 10000000000 \
  --backend dense --output-prefix goldbach_to_10B_v03 \
  --progress-every 10000000
