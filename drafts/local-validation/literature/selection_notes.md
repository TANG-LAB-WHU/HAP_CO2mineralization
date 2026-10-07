# Method selection and access limits

Parent-supplied primary-source/official-project audit dated 2026-10-07, condensed locally.
The method pool is not a recommendation to cite every entry.

- CP2K: Kuhne2020CP2K (10.1063/5.0007045), publisher/author manuscript arXiv:2003.03868 and official Quickstep bibliography. Software attribution does not establish cutoff, basis, functional or surface accuracy. Cite actual functional/pseudopotential/basis and dispersion separately where used.
- MACE architecture: Batatia2022MACE (10.52202/068431-0830), official NeurIPS proceedings and arXiv:2206.07697. No transfer of molecular benchmark accuracy to HAP; version wording differs, so specify actual model rather than universal layer-count claims.
- MACE-MP: Batatia2025MACEMP (10.1063/5.0297006), now published JCP163,184110; arXiv:2401.00096v3 limitations section. Parent record uses standard 'and others'. Identify exact checkpoint; do not relabel every MACE model as MACE-MP. Reference PBE(+U), absent explicit long-range/spin degrees and intermolecular limitations require material-specific validation. No dispersion correction was enabled for this MACE checkpoint in the local tests.
- LAMMPS: Thompson2022LAMMPS (10.1016/j.cpc.2021.108171), volume publication2022; official doc.lammps.org/Intro_citing.html. Engine citation is not potential/reactivity validation.
- DeepMD: software citations and original local-frame DPMD method are distinct. Actual local image is v3.1.0 PyTorch, so use the v3 JCTC2025 paper (10.1021/acs.jctc.5c00340), plus the smooth descriptor method if discussed. Exact official version citation file inspected locally during authoring: https://raw.githubusercontent.com/deepmodeling/deepmd-kit/v3.1.0/CITATIONS.bib . Do not imply long-range features are enabled simply by citing software with that capability.
- MatterGen: Zeni2025MatterGen (10.1038/s41586-025-08628-5), Nature639,624-632. Generative crystals, not a force field. Official project warns MatterSim evaluation differs from paper DFT. Most source training structures have at most20atoms/unitcell; generated chemistry must be checked. Optional here; rejected sample is not a HAP candidate.
- MatterSim: Yang2024MatterSim, arXiv:2405.04967 and official MODEL_CARD.md. No journal publication was established in supplied audit. Model card flags surfaces/interfaces/long-range settings as potentially requiring fine-tuning. Not a crystal generator or universal validation.

No general methods paper independently supports HAP capture capacity, an adsorption
free energy, a carbonation pathway, a reaction barrier, or durable storage. Cite only
actual methods and only the source needed for the attached claim. Never use nocite{*}.
