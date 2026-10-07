# Final local validation handoff

Authoritative checkout: `D:\project\HAP_CO2mineralization`, on Shawn-PC.
Branch: `integration/shawn-4090-reproducible`. The local final commit, patch, file hashes
and verification log index are in `deliverables/DELIVERY_MANIFEST.json` and
`deliverables/README_local_results.md`. Nothing was pushed or uploaded.

## Deliverables

All files below are under `deliverables/`:

- `HAP_local_validation_draft.pdf`: 23-page diagnostic working draft, all pages visually reviewed.
- `HAP_local_validation_overleaf.zip`: self-contained XeLaTeX/BibTeX project with 19 cited references and 11 figures; clean offline build tested.
- `HAP_stable_250fs_demonstration.gif`: 51 actual frames; explicitly limited to the stable short window.
- `HAP_recovered_actual_1ps.gif`: separate 101-frame failure animation; failure begins near 395 fs.
- `HAP_LAMMPS_actual_trajectory.gif`: retained original failed 250 fs experiment.
- Matching OVITO dump/extxyz files, `HAP_validation_diagnostics.zip`, and disabled future-sampling plan.
- `HAP_local_integration.patch`: committed integration diff from the fetched `shawn_dev` base.

All 203 GIF frames were decoded; first/middle/final frames of each animation were
actually inspected. These are saved trajectories, with no interpolation. An OVITO
GUI session is not claimed. The clean ZIP rebuild has identical text and page count;
22 pages are pixel-identical, while one bibliography URL line on page 19 moves 0.2864 pt.
That page was separately reviewed. Pixel identity of all pages is not claimed.

## Verify existing calculations without rerunning them

```powershell
Set-Location D:\project\HAP_CO2mineralization
.\workflow\run-refinement.ps1 -Stage all -Run refine-001
python workflow/capture_refinement.py --run .research-runs/refine-001 --verify
python workflow/capture_recovery.py --run .research-runs/recovery-001 --verify
.\workflow\run-local.ps1 -Stage check
```

Captured runs are immutable. The first command verifies and returns; it does not
recalculate physics or rewrite figures. The two captures contain 658 and 244 files.
The aggregate command passes 170 tests and the evidence checker; static checks do
not override failed scientific gates. Full outputs and checksums remain under
`.research-runs/final-validation/`. Use a fresh run name for a new experiment.

## Rebuild the document on this configured machine

```powershell
Set-Location D:\project\HAP_CO2mineralization
python workflow/build_validation_paper.py --run .research-runs/refine-001 --recovery .research-runs/recovery-001
$env:TECTONIC_CACHE_DIR="$PWD\.research-runs\tools\tex-cache"
& '.\.research-runs\tools\tectonic\tectonic.exe' --only-cached --keep-logs --keep-intermediates --outdir drafts/local-validation drafts/local-validation/main.tex
```

Review a rebuilt PDF and bind `visual_qa.json` to its new SHA256 before packaging:

```powershell
python workflow/package_validation_paper.py --run .research-runs/refine-001 --recovery .research-runs/recovery-001
```

For a clean ZIP build, extract it to a fresh local folder and run the same cached
Tectonic executable on its `main.tex`; no raw calculation data are required. Overleaf
should use `main.tex`, XeLaTeX, and BibTeX. No Overleaf login/upload was performed.

## Scientific limits and next use

Optimizer-only recovery lowers one-frame force RMSE to 0.016476 eV/A, passing its
unchanged 0.020 gate. The validation-selected candidate's reused-test force RMSE is
0.107883 eV/A, failing 0.100; this is retrospective/exploratory, not a new holdout.
The first 250 fs has no failed geometry frames and passes strict restart, but the
1 ps runs fail at both timesteps. At the first failure, DFT raises energy by 5.965761
eV while the model lowers its own energy by 2.289111 eV; force RMSE is 1.283081 eV/A.
Each energy change uses its own method's initial energy. No production potential or
mineralization result is claimed.

`RECOVERY.md` records the full verification matrix, optimizer budgets and an exact
fresh-folder replay command. `PARAMETERS.md`, `REFINEMENT.md` and `REPRODUCE_EXTENDED.md`
separate smoke settings from candidate production studies. `NEXT_SAMPLING_PLAN.md`
defines an unexecuted maximum of 159 new DFT attempts, fixed thresholds, targeted
distortions, at most two acquisition rounds and a new sealed holdout. No more DFT
was started for final delivery.

`HPC_HANDOFF.md` and `hpc/` provide disabled job, restart and dependency templates.
Remaining site details are scheduler/account/partition/QoS, modules or Apptainer/image
paths, allocated hardware and resource limits, and shared/scratch paths. No cluster
connection, submission or cluster validation occurred. A pristine workstation bootstrap
was not repeated; the installed pinned local runtimes were actually exercised.
