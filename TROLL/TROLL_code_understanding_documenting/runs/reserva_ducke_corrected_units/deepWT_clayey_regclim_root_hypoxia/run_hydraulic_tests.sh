#!/bin/bash

set -euo pipefail

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Find project root by walking up until mainTROLL4.0_WTD.cpp is found
find_project_root() {
  local dir="$1"
  while [ "$dir" != "/" ]; do
    if [ -f "$dir/mainTROLL4.0_WTD.cpp" ]; then
      echo "$dir"
      return 0
    fi
    dir="$(dirname "$dir")"
  done
  return 1
}

BASE_DIR="$(find_project_root "$SCRIPT_DIR")" || {
  echo "Error: could not find project root containing mainTROLL4.0_WTD.cpp"
  exit 1
}

# Assume the script is inside the experiment branch folder
# Example:
# .../runs/theta_w_tests/veg/deepWT/sandy/
BRANCH_DIR="$SCRIPT_DIR"
INPUT_DIR="$BRANCH_DIR/common_inputs"

# Ask for run name if not provided as argument
if [ $# -ge 1 ]; then
  RUN="$1"
else
  read -r -p "Enter run name: " RUN
fi

if [ -z "$RUN" ]; then
  echo "Error: run name cannot be empty."
  exit 1
fi

RUN_DIR="$BRANCH_DIR/$RUN"
OUTPUT_DIR="$RUN_DIR/output"

EXE_FILE="$BASE_DIR/TROLL.out"

CLIMATE_FILE="$INPUT_DIR/Paracou_input_climate.txt"
DAILY_FILE="$INPUT_DIR/Paracou_input_daily.txt"
SPECIES_FILE="$INPUT_DIR/Paracou_input_species.txt"
GLOBAL_FILE="$INPUT_DIR/Paracou_input_global.txt"
PEDOLOGY_FILE="$INPUT_DIR/Paracou_input_pedology.txt"

# Check required input files
for f in \
  "$CLIMATE_FILE" \
  "$DAILY_FILE" \
  "$SPECIES_FILE" \
  "$GLOBAL_FILE" \
  "$PEDOLOGY_FILE"
do
  if [ ! -f "$f" ]; then
    echo "Error: required input file not found: $f"
    exit 1
  fi
done

mkdir -p "$RUN_DIR"
mkdir -p "$OUTPUT_DIR"

echo "Running experiment: $RUN"
echo "Branch directory: $BRANCH_DIR"
echo "Project root: $BASE_DIR"

GIT_COMMIT=$(git -C "$BASE_DIR" rev-parse HEAD)
GIT_COMMIT_SHORT=$(git -C "$BASE_DIR" rev-parse --short HEAD)
GIT_BRANCH=$(git -C "$BASE_DIR" branch --show-current)

{
  echo "run_name: $RUN"
  echo "date: $(date)"
  echo
  echo "branch_dir: $BRANCH_DIR"
  echo "project_root: $BASE_DIR"
  echo
  echo "git_commit: $GIT_COMMIT"
  echo "git_commit_short: $GIT_COMMIT_SHORT"
  echo "git_branch: $GIT_BRANCH"
  echo
  echo "climate_file: $CLIMATE_FILE"
  echo "daily_file: $DAILY_FILE"
  echo "species_file: $SPECIES_FILE"
  echo "global_file: $GLOBAL_FILE"
  echo "pedology_file: $PEDOLOGY_FILE"
  echo
  echo "========================"
  echo "GLOBAL FILE CONTENT"
  echo "========================"
  cat "$GLOBAL_FILE"
  echo
  echo "=========================="
  echo "PEDOLOGY FILE CONTENT"
  echo "=========================="
  cat "$PEDOLOGY_FILE"
} > "$RUN_DIR/config_used.txt"

echo "Compiling..."

cd "$BASE_DIR"
GSL_PREFIX="$(brew --prefix gsl)"

g++ mainTROLL4.0_WTD.cpp -O3 -o TROLL.out \
  -I"$GSL_PREFIX/include" \
  -L"$GSL_PREFIX/lib" \
  -lgsl -lgslcblas \
  -Wall 

echo "Running model..."

cd "$OUTPUT_DIR"
"$EXE_FILE" \
  -m"$CLIMATE_FILE" \
  -d"$DAILY_FILE" \
  -s"$SPECIES_FILE" \
  -i"$GLOBAL_FILE" \
  -p"$PEDOLOGY_FILE" \
  2>&1 | tee "$RUN_DIR/run.log"

echo "Run completed: $RUN"