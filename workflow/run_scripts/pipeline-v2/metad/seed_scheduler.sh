#!/bin/bash
# Keeps seed-metad (CP2K labels of the biased frames) off the CPU while DeePMD trains, and lets it run while
# LAMMPS explores. Training slows by about a third with 20 CP2K ranks beside it (PIPELINE_V2.md, "DeePMD training
# and CPU load"); exploration is the CPU window used before (seed-fs). During its own label steps the orchestrator
# pauses and resumes seed containers itself, so this loop never resumes one outside an exploration step.
# 2026-11-07: any dp_explore.py container counted as exploration, so the NPT smoke test (run by hand on 2026-11-06
# 17:35, during label-4) resumed seed-metad beside the label step; label-4 then ran at 7-9 frames per hour instead
# of 16-23 until train-4. Another CP2K container (a label or relabel step) now keeps seed-metad paused.
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*"; }
while true; do
  st=$(docker inspect -f '{{.State.Status}}' seed-metad 2>/dev/null) || { log "seed-metad gone"; exit 0; }
  [ "$st" = exited ] && { log "seed-metad exited"; exit 0; }
  training=$(docker ps --format '{{.Names}}' | grep -c '^dp-train-')
  exploring=$(docker ps --filter ancestor=hap-mlip:2026.10 --no-trunc --format '{{.Command}}' | grep -c 'dp_explore.py')
  labelling=$(docker ps --filter ancestor=cp2k/cp2k:2026.2_mpich_generic_psmp --format '{{.Names}}' | grep -vc '^seed-')
  if [ "$training" -gt 0 ] && [ "$st" = running ]; then
    docker pause seed-metad >/dev/null && log "paused (training)"
  elif [ "$labelling" -gt 0 ] && [ "$st" = running ]; then
    docker pause seed-metad >/dev/null && log "paused (labelling)"
  elif [ "$training" -eq 0 ] && [ "$labelling" -eq 0 ] && [ "$exploring" -gt 0 ] && [ "$st" = paused ]; then
    docker unpause seed-metad >/dev/null && log "resumed (exploration)"
  fi
  sleep 60
done
