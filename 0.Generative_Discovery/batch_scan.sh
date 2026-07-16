#!/bin/bash
# Batch scan multiple chemical systems through MatterGen
# All generation parameters (batch-size, num-batches, energy-above-hull)
# are controlled by predict_mattergen.py's argparse defaults.
#
# Pipeline per system:
#   1. predict_mattergen.py  → generate + evaluate
#   2. generation_stats.py   → print structural statistics
#   3. post_filter.py        → filter stable & novel candidates to CIF

SYSTEMS=(
    "Ca-P-O-H"
    "Ca-P-O-H-C"
    "Ca-P-O-H-F"
    "Ca-P-O-H-Cl"
    "Ca-Sr-P-O-H"
    "Ca-Mg-P-O-H"
    "Ca-P-Si-O-H"
)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TOTAL=${#SYSTEMS[@]}
SUCCESS=0
FAIL=0

echo "================================================="
echo "Starting MatterGen Batch Scan"
echo "Systems: ${SYSTEMS[*]}"
echo "Total systems: $TOTAL"
echo "================================================="

for SYS in "${SYSTEMS[@]}"; do
    echo ""
    echo "=========================================="
    echo "=== [$((SUCCESS + FAIL + 1))/$TOTAL] System: $SYS ==="
    echo "=========================================="

    RESULTS_DIR="$SCRIPT_DIR/results/$SYS"

    # Step 1: Generate + Evaluate
    echo "[Step 1/3] Generating and evaluating structures..."
    python "$SCRIPT_DIR/predict_mattergen.py" \
        --chemical-system="$SYS" \
        --model-name=chemical_system_energy_above_hull \
        --output-path="$RESULTS_DIR" \
        --evaluate \
        --repo-root=../../mattergen

    if [ $? -ne 0 ]; then
        echo "WARNING: Generation failed for $SYS, skipping stats and filter."
        FAIL=$((FAIL + 1))
        continue
    fi

    # Step 2: Print generation statistics
    EXTXYZ_FILE="$RESULTS_DIR/generated_crystals.extxyz"
    if [ -f "$EXTXYZ_FILE" ]; then
        echo ""
        echo "[Step 2/3] Structural statistics for $SYS:"
        python "$SCRIPT_DIR/generation_stats.py" "$EXTXYZ_FILE"
    else
        echo "WARNING: $EXTXYZ_FILE not found, skipping stats."
    fi

    # Step 3: Filter stable & novel candidates
    echo ""
    echo "[Step 3/3] Filtering candidates for $SYS..."
    python "$SCRIPT_DIR/post_filter.py" \
        --results-dir="$RESULTS_DIR"

    if [ $? -eq 0 ]; then
        # Count filtered candidates
        CAND_DIR="$RESULTS_DIR/filtered_candidates"
        if [ -d "$CAND_DIR" ]; then
            NUM_CIFS=$(find "$CAND_DIR" -name "*.cif" | wc -l | tr -d ' ')
            echo "→ $NUM_CIFS candidate CIFs saved to $CAND_DIR"
        fi
    fi

    SUCCESS=$((SUCCESS + 1))
done

echo ""
echo "================================================="
echo "Batch scan completed."
echo "Success: $SUCCESS / $TOTAL"
if [ $FAIL -gt 0 ]; then
    echo "Failed:  $FAIL / $TOTAL"
fi
echo "================================================="
