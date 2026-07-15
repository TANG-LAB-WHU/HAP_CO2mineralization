#!/bin/bash
# Batch scan multiple chemical systems through MatterGen

SYSTEMS=(
    "Ca-P-O-H"
    "Ca-P-O-H-C"
    "Ca-P-O-H-F"
    "Ca-P-O-H-Cl"
    "Ca-Sr-P-O-H"
    "Ca-Mg-P-O-H"
    "Ca-P-Si-O-H"
)

# Number of batches to generate (64 structures per batch)
BATCHES=${1:-5}
# Energy above hull threshold
EHULL=${2:-0.05}

echo "================================================="
echo "Starting MatterGen Batch Scan"
echo "Batches per system: $BATCHES (x64 = $((BATCHES*64)) structures)"
echo "Energy above hull threshold: $EHULL"
echo "================================================="

for SYS in "${SYSTEMS[@]}"; do
    echo ""
    echo "=== Scanning System: $SYS ==="
    
    python predict_mattergen.py \
        --chemical-system="$SYS" \
        --model-name=chemical_system_energy_above_hull \
        --energy-above-hull=$EHULL \
        --batch-size=64 \
        --num-batches=$BATCHES \
        --evaluate \
        --repo-root=../../mattergen
        
    if [ $? -ne 0 ]; then
        echo "Warning: Scanning failed or interrupted for $SYS"
        # Continue to next system
    fi
done

echo ""
echo "Batch scan completed."
