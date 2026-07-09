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

# ----------------------------------------------------------------------
# Arguments
# ----------------------------------------------------------------------
# Example:
# ./run_model.sh clayey clayey_WTD_4
# ./run_model.sh sandy sandy_WTD_4
# ----------------------------------------------------------------------

if [ $# -ge 1 ]; then
  SOIL="$1"
else
  read -r -p "Enter soil type: " SOIL
fi

if [ -z "$SOIL" ]; then
  echo "Error: soil type cannot be empty."
  exit 1
fi

SOIL_DIR="$SCRIPT_DIR/$SOIL"

if [ ! -d "$SOIL_DIR" ]; then
  echo "Error: soil directory not found: $SOIL_DIR"
  echo
  echo "Available soil directories:"
  find "$SCRIPT_DIR" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort
  exit 1
fi

if [ $# -ge 2 ]; then
  RUN="$2"
else
  echo "Available WTD scenario folders inside $SOIL_DIR:"
  find "$SOIL_DIR" -mindepth 1 -maxdepth 1 -type d ! -name "common_inputs" -exec basename {} \; | sort
  echo
  read -r -p "Enter WTD scenario folder name: " RUN
fi

if [ -z "$RUN" ]; then
  echo "Error: run name cannot be empty."
  exit 1
fi

RUN_DIR="$SOIL_DIR/$RUN"
INPUT_DIR="$SOIL_DIR/common_inputs"
OUTPUT_DIR="$RUN_DIR/output"

if [ ! -d "$RUN_DIR" ]; then
  echo "Error: run directory not found: $RUN_DIR"
  exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
  echo "Error: common input directory not found: $INPUT_DIR"
  exit 1
fi

EXE_FILE="$BASE_DIR/TROLL.out"

CLIMATE_FILE="$INPUT_DIR/Paracou_input_climate.txt"
DAILY_FILE="$INPUT_DIR/Paracou_input_daily.txt"
SPECIES_FILE="$INPUT_DIR/Paracou_input_species.txt"
PEDOLOGY_FILE="$INPUT_DIR/Paracou_input_pedology.txt"

# IMPORTANT:
# The global file is specific to each WTD scenario
GLOBAL_FILE="$RUN_DIR/Paracou_input_global.txt"

# ----------------------------------------------------------------------
# Check required input files
# ----------------------------------------------------------------------

for f in \
  "$CLIMATE_FILE" \
  "$DAILY_FILE" \
  "$SPECIES_FILE" \
  "$PEDOLOGY_FILE" \
  "$GLOBAL_FILE"
do
  if [ ! -f "$f" ]; then
    echo "Error: required input file not found: $f"
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

echo "Running experiment: $RUN"
echo "Soil: $SOIL"
echo "Soil directory: $SOIL_DIR"
echo "Run directory: $RUN_DIR"
echo "Common input directory: $INPUT_DIR"
echo "Project root: $BASE_DIR"

# ----------------------------------------------------------------------
# Git information
# ----------------------------------------------------------------------

GIT_COMMIT=$(git -C "$BASE_DIR" rev-parse HEAD)
GIT_COMMIT_SHORT=$(git -C "$BASE_DIR" rev-parse --short HEAD)
GIT_BRANCH=$(git -C "$BASE_DIR" branch --show-current)

# ----------------------------------------------------------------------
# Save configuration used
# ----------------------------------------------------------------------

{
  echo "run_name: $RUN"
  echo "soil: $SOIL"
  echo "date: $(date)"
  echo
  echo "soil_dir: $SOIL_DIR"
  echo "run_dir: $RUN_DIR"
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

# ----------------------------------------------------------------------
# Compile
# ----------------------------------------------------------------------

echo "Compiling..."

cd "$BASE_DIR"
GSL_PREFIX="$(brew --prefix gsl)"

g++ mainTROLL4.0_WTD.cpp -O3 -o TROLL.out \
  -I"$GSL_PREFIX/include" \
  -L"$GSL_PREFIX/lib" \
  -lgsl -lgslcblas \
  -Wall

# ----------------------------------------------------------------------
# Run model
# ----------------------------------------------------------------------

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