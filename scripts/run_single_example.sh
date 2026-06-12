#!/usr/bin/env bash
set -euo pipefail

./goldbach_lab_v03 --mode single --N 1000000000000 \
  --backend miller --output-prefix goldbach_single_1T_v03
