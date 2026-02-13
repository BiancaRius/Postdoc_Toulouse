#!/usr/bin/env zsh
# set -euo pipefail

# Absolute path to the TROLL code directory (where the .cpp lives and where we build the executable)
CODEDIR="/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting"

# Base directory containing the experiment folders
BASE="$CODEDIR/runs/sensitivity_tests"

# Directory containing inputs shared across all experiments (climate, daily, species, etc.)
COMMON="$BASE/common_inputs"

# Executable path (single binary reused across runs). You can rename it if you want.
EXE="$CODEDIR/TROLL_WTD.out"

# Main C++ source file to compile
SRC="$CODEDIR/mainTROLL4.0_WTD.cpp"

# List of sandy experiments to run (folder names under $BASE)
exps=(
  shallow_variable_sandy
  shallow_intermediate_sandy
  shallow_thin_sandy
  deep_variable_sandy
  deep_intermediate_sandy
  deep_thin_sandy
)

# ------------------------------------------------------------
# Compile only if the executable does not exist.
# (If you want to force recompilation, delete $EXE before running.)
# ------------------------------------------------------------
if [[ ! -f "$EXE" ]]; then
  echo "== Compiling: $EXE =="
  g++ "$SRC" -O3 -o "$EXE" -lgsl -lgslcblas -Wall
fi

# Collect useful metadata (safe defaults if not available)
NOW=$(date +"%Y-%m-%d %H:%M:%S")

# Git commit hash (if CODEDIR is a git repo)
GITHASH="NA"
if command -v git >/dev/null 2>&1 && git -C "$CODEDIR" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  GITHASH=$(git -C "$CODEDIR" rev-parse HEAD 2>/dev/null || echo "NA")
fi

# Compiler version (first line of `g++ --version`)
COMPILER="NA"
if command -v g++ >/dev/null 2>&1; then
  COMPILER=$(g++ --version | head -n 1)
fi

# ------------------------------------------------------------
# Loop over experiments and run TROLL for each one
# ------------------------------------------------------------
for EXP in "${exps[@]}"; do
  # Experiment directory (contains the experiment-specific inputs)
  EXPDIR="$BASE/$EXP"

  # Output directory for this experiment
  OUTDIR="$EXPDIR/outputs"
  mkdir -p "$OUTDIR"

  # Output prefix: controls BOTH where outputs go and the filename prefix
  # (TROLL will append rank and file suffixes to this prefix)
  OUTPREFIX="$OUTDIR/$EXP"

  # Build the command as an array (safer than a single string).
  # IMPORTANT: This TROLL parser expects flags in the form -mPATH (no space).
  CMD=(
    "$EXE"
    -m"$COMMON/Paracou_input_climate.txt"
    -d"$COMMON/Paracou_input_daily.txt"
    -s"$COMMON/Paracou_input_species.txt"
    -i"$EXPDIR/Paracou_input_global.txt"
    -p"$EXPDIR/Paracou_input_pedology.txt"
    -o"$OUTPREFIX"
  )

  echo "== Running $EXP =="
  echo "Output prefix: $OUTPREFIX"

  # Paths for run documentation files written inside the experiment folder
  README="$EXPDIR/README.md"
  MANIFEST="$EXPDIR/run_manifest.txt"

  # Write a human-readable README (overwrites each time; last run only)
  {
    echo "# Experiment: $EXP"
    echo ""
    echo "## Last run"
    echo "- Date (local): $NOW"
    echo "- Executable: $EXE"
    echo "- Git commit: $GITHASH"
    echo "- Compiler: $COMPILER"
    echo ""
    echo "## Inputs"
    echo "- Climate: $COMMON/Paracou_input_climate.txt"
    echo "- Daily: $COMMON/Paracou_input_daily.txt"
    echo "- Species: $COMMON/Paracou_input_species.txt"
    echo "- Global: $EXPDIR/Paracou_input_global.txt"
    echo "- Pedology: $EXPDIR/Paracou_input_pedology.txt"
    echo ""
    echo "## Outputs"
    echo "- Output directory: $OUTDIR"
    echo "- Output prefix: $OUTPREFIX"
    echo ""
    echo "## Command"
    echo '```'
    printf "%q " "${CMD[@]}"
    echo ""
    echo '```'
    echo ""
    echo "## Notes"
    echo "- This experiment differs from others only by its pedology (and possibly global) inputs."
    echo "- Outputs are written using the prefix above (TROLL appends rank + file suffix)."
  } > "$README"

  # Write a machine-friendly manifest (easy to parse later)
  {
    echo "EXP=$EXP"
    echo "DATE_LOCAL=$NOW"
    echo "EXECUTABLE=$EXE"
    echo "GIT_COMMIT=$GITHASH"
    echo "COMPILER=$COMPILER"
    echo "CLIMATE=$COMMON/Paracou_input_climate.txt"
    echo "DAILY=$COMMON/Paracou_input_daily.txt"
    echo "SPECIES=$COMMON/Paracou_input_species.txt"
    echo "GLOBAL=$EXPDIR/Paracou_input_global.txt"
    echo "PEDOLOGY=$EXPDIR/Paracou_input_pedology.txt"
    echo "OUTDIR=$OUTDIR"
    echo "OUTPREFIX=$OUTPREFIX"
    echo -n "CMD="
    printf "%q " "${CMD[@]}"
    echo ""
  } > "$MANIFEST"

  # Run TROLL
  "${CMD[@]}"

done

echo "All sandy experiments finished."
