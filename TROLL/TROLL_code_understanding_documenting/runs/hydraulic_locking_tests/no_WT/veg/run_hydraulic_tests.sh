#!/bin/bash

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: ./run_hydraulic_tests.sh RUN_NAME"
  exit 1
fi

RUN="$1"

BASE_DIR="/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting"
RUN_DIR="$BASE_DIR/runs/hydraulic_locking_tests/no_WT/veg/$RUN"
OUTPUT_DIR="$RUN_DIR/output"
INPUT_DIR="$BASE_DIR/runs/hydraulic_locking_tests/no_WT/veg/common_inputs"

EXE_FILE="$BASE_DIR/TROLL.out"

CLIMATE_FILE="$INPUT_DIR/Paracou_input_climate.txt"
DAILY_FILE="$INPUT_DIR/Paracou_input_daily.txt"
SPECIES_FILE="$INPUT_DIR/Paracou_input_species.txt"
GLOBAL_FILE="$INPUT_DIR/Paracou_input_global.txt"
PEDOLOGY_FILE="$INPUT_DIR/Paracou_input_pedology.txt"

mkdir -p "$RUN_DIR"
mkdir -p "$OUTPUT_DIR"

echo "Running experiment: $RUN"

GIT_COMMIT=$(git -C "$BASE_DIR" rev-parse HEAD)
GIT_COMMIT_SHORT=$(git -C "$BASE_DIR" rev-parse --short HEAD)
GIT_BRANCH=$(git -C "$BASE_DIR" branch --show-current)

cat > "$RUN_DIR/config_used.txt" <<EOF
run_name: $RUN
date: $(date)

git_commit: $GIT_COMMIT
git_commit_short: $GIT_COMMIT_SHORT
git_branch: $GIT_BRANCH
EOF

echo "Compiling..."

cd "$BASE_DIR"
g++ mainTROLL4.0_WTD.cpp -O3 -o TROLL.out -lgsl -lgslcblas -Wall

echo "Running model..."

cd "$OUTPUT_DIR"
"$EXE_FILE" \
  -m"$CLIMATE_FILE" \
  -d"$DAILY_FILE" \
  -s"$SPECIES_FILE" \
  -i"$GLOBAL_FILE" \
  -p"$PEDOLOGY_FILE" \

echo "Run completed: $RUN"