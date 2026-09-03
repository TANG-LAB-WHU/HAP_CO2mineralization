# Research-agent contract

This repository combines scientific code, evidence, and a living manuscript. Agents must preserve scientific provenance.

## Default scope

- Work on `shawn_dev` unless the user explicitly names another branch.
- Treat `research/research_state.yml` as the current project state.
- Treat `literature/evidence.csv` as the claim-to-evidence ledger.
- Treat `paper/references.bib` as the only manuscript bibliography.
- Keep large trajectories, checkpoints, generated datasets, PDFs, credentials, and machine-local profiles out of Git.

## Scientific integrity

- Never invent numerical results, completed experiments, citations, DOIs, author lists, or file provenance.
- A result may enter the manuscript only when its source path or run manifest is recorded in `literature/evidence.csv`.
- Mark unsupported prose as `TODO-EVIDENCE` and proposed work as `PLANNED`; do not phrase either as a finding.
- Verify bibliographic metadata against a primary publisher, DOI registry, or official proceedings before adding a citation.
- Preserve units, uncertainty, sample counts, seeds, convergence criteria, and negative results.

## Editing and validation

- Prefer small reviewable changes. Do not reorganize the existing numbered computation directories without approval.
- Run `make check` after changing the manuscript, bibliography, or evidence ledger.
- Run `make paper` before proposing a writing-related commit.
- Put publication figures under `paper/figures/`; record the generating script and source data in its README.

## Compute boundary

- The default workflow is local writing and lightweight validation on macOS.
- HPC integration is disabled. Never store VPN passwords, SSH passwords, private keys, tokens, or secrets in this repository.
- Do not run `ssh`, `scp`, `rsync`, `sbatch`, `srun`, `scancel`, or destructive/large computations unless the user explicitly authorizes that action.
- Future remote runs must emit a manifest compatible with `workflow/schemas/run_manifest.example.json`.

