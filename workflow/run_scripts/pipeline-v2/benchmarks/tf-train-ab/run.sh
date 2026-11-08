#!/bin/bash
# A/B: the iteration-0 model-0 configuration (same data, seed, sel and 300k steps) trained with the TensorFlow backend
cd /work/.research-runs/pipeline-v2/benchmarks/tf-train-ab/model-0 || exit 1
t0=$(date +%s)
dp --tf train input.json > train.log 2>&1 || { echo "train failed"; exit 1; }
echo "train wall $(( $(date +%s) - t0 )) s" > timing.txt
dp --tf freeze -o graph.pb > freeze.log 2>&1 || { echo "freeze failed"; exit 1; }
for s in /work/.research-runs/pipeline-v2/sealed-test/data/train/*__*; do
  dp --tf test -m graph.pb -s $s -n 100000 > test_$(basename $s).log 2>&1
done
echo done > done
