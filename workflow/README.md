# Compute handoff contract

The current workflow supports local writing, literature management, validation, and rendering. Existing scientific scripts remain in the numbered project directories. This folder defines how future local or remote computations report provenance back to the paper.

## Profiles

- `profiles/local.yaml`: enabled for lightweight local checks and analysis.
- `profiles/whu-hpc.example.yaml`: disabled Slurm template; contains no credentials.
- `profiles/whu-hpc.local.yaml`: future machine-local configuration; ignored by Git.

Profiles are descriptive in v1. Nothing in the writing workflow automatically connects to or submits work on a remote machine.

## Run manifests

Every result used by the paper should eventually have a compact JSON manifest based on `schemas/run_manifest.example.json`. Store manifests under a stage-appropriate tracked directory or a future `workflow/runs/` index; keep large outputs externally and reference them by stable path plus checksum.

The required chain is:

```text
research question → configuration → Git commit → compute run → manifest → evidence row → manuscript claim
```

