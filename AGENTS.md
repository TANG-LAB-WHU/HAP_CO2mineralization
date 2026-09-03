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
- Never generate synthetic or "realistic sample" scientific data as a fallback for missing, malformed, or all-zero inputs. Fail closed and record the missing evidence instead.
- Treat every legacy computational artifact as manuscript-ineligible until `research/artifact_inventory.csv` records adequate provenance and an explicit eligibility decision.
- Normal program termination is not evidence of scientific validity. Check convergence, stability, completeness, and predefined acceptance criteria independently.
- A literature claim is `supported` only after its citation metadata and a claim-specific source locator have been verified and recorded. An existing bibliography entry alone is not verification.

## Editing and validation

- Prefer small reviewable changes. Do not reorganize the existing numbered computation directories without approval.
- Run `make check` after changing the manuscript, bibliography, or evidence ledger.
- Run `make paper` before proposing a writing-related commit.
- Put publication figures under `paper/figures/`; record the generating script and source data in its README.
- Keep `TODO-EVIDENCE <claim-id>` markers synchronized with non-supported rows in `literature/evidence.csv`.
- Do not copy a figure into `paper/figures/` unless every source artifact is marked `manuscript_eligible=yes`.
- Preserve failed, partial, and negative runs in the inventory; never relabel or hide them to improve the apparent result set.

## Compute boundary

- The default workflow is local writing and lightweight validation on macOS.
- HPC integration is disabled. Never store VPN passwords, SSH passwords, private keys, tokens, or secrets in this repository.
- Do not run `ssh`, `scp`, `rsync`, `sbatch`, `srun`, `scancel`, or destructive/large computations unless the user explicitly authorizes that action.
- Future remote runs must emit a manifest compatible with `workflow/schemas/run_manifest.example.json`.
