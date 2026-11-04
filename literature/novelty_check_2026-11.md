# Novelty check of candidate contributions (2026-11-03)

Requested by Shawn in chat on 2026-11-03. This check covers the four contributions proposed for pipeline v2. It asks only whether published or preprint work already covers each one.

## Scope and limits

- **Services:**
  - Web search, standard and extended modes (20 queries).
  - The DeePMD-kit publication list, read in full.
  - Crossref and arXiv metadata previews.

  Queries and screening decisions are in `search_log.md` (2026-11-03).
- **Not searched:** Web of Science, Scopus, Google Scholar, ChemRxiv directly, and CNKI. A "no close precedent found" verdict below means "none found in this search". It does not establish novelty.
- **Basis of judgement:** abstracts and search snippets only; no full text was read.
- **Status of the references:** the 12 new references are `candidate` in `references.csv`. No `evidence.csv` relationship changed, and no manuscript text changed.

## Verdicts

| # | Candidate contribution | Verdict | Closest precedents |
| --- | --- | --- | --- |
| 1 | DFT-level MLIP for hydrated HAP surfaces with CO2 | Partial precedent | `wang2023hapmlp`, `wang2024hapmethane`, `wang2025capnucleation`, `piaggi2025caco3`, `raman2025mica` |
| 2 | Training set that also covers carbonation products: MatterGen candidates that pass a DFT hull join the active-learning training set and the sealed test | No close precedent found | `podryabinkin2019csp` (conceptual only) |
| 3 | Quantified attrition of a generative + foundation-MLIP screen under DFT | Partial precedent (general finding known) | `szymanski2025baselines`, `han2025phononbench` |
| 4 | CO2 → HCO3⁻/CO3²⁻ mechanism at hydrated HAP surfaces, with biased active learning | System not found; methods established | `zhu2025portlandite`, `yang2022metadynamics`, `bittencourt2021happrobe`, `lee2015scco2`, `piaggi2025caco3` |

### 1. MLIP for hydrated HAP + CO2

- **HAP MLIPs exist.** One Sichuan University group (Xu Dingguo and co-workers) has published:
  - an active-learning MLIP for bulk HAP OH reorientation (`wang2023hapmlp`);
  - an MLIP-assisted study of methane activation on HAP surfaces with defects (`wang2024hapmethane`);
  - an MLIP study of calcium phosphate nucleation in water (`wang2025capnucleation`).

  None of these involves CO2 or a HAP surface under water.
- **Other related potentials:**
  - A Deep Potential for lead apatite exists (He et al., Phys. Rev. B 108, 224114, 2023, DOI 10.1103/PhysRevB.108.224114). It is listed in `search_log.md` but was not added because the round cap was reached.
  - A SCAN Deep Potential covers CaCO3 from solution to solid, including bicarbonate chemistry (`piaggi2025caco3`).
  - Deep Potentials for mineral–water interfaces are established practice (`raman2025mica`).
- **Not found:** an MLIP trained on HAP surfaces with both water and CO2.
- **Honest phrasing:** "extends HAP machine-learning potentials [wang2023hapmlp; wang2024hapmethane] to hydrated HAP surfaces with CO2 and to the candidate carbonation products". Do not write "first" until a database search confirms it. That search should cover the same group's 2025–2026 output specifically.

### 2. Training set that also covers carbonation products

- **Not found:** work that feeds generative-model candidate phases, after DFT hull checks, into the training set and sealed test of an interfacial reaction potential.
- **Closest idea:** active learning coupled to crystal-structure search (`podryabinkin2019csp`). That work aims to accelerate structure prediction, not to cover reaction products.
- **Honest phrasing:** present it as a design choice with a measured consequence, not as a first.
  - The consequence is in the per-system sealed-test errors and the stage-7 stability runs for the S3 phases.
  - Without those numbers it is only a workflow detail.

### 3. Attrition of the generative screen

- **Already known in general:**
  - Stability rates of generated novel materials fall below the published S.U.N. rates (`szymanski2025baselines`).
  - Most generated crystals are dynamically unstable. PhononBench (`han2025phononbench`, preprint) reports 32 % stable on average and 45 % for MatterGen, from MatterSim phonons at −0.1 THz.
  - The MatterGen repository README itself recommends DFT confirmation of results evaluated with an MLIP. That README is not a paper and is not registered.
- **Specific to this project:**
  - the scan is conditioned on the Ca–P–O–H–C system;
  - the hull uses the downstream level of theory, CP2K PBE-D3(BJ) with 34 computed references;
  - three MLIPs agree within 2.5–3.3 meV/atom (median pairwise), yet differ from DFT by −12.8 to +20.7 meV/atom;
  - the possible systematic shift for phases whose decomposition includes ice is a hypothesis, not a finding.
- **Not directly comparable:** our S3 pass rate (10/149) used MACE phonons at −0.3 THz on an 8³ mesh plus 5 ps NPT. PhononBench used MatterSim at −0.1 THz.
- **Honest phrasing:** "consistent with [szymanski2025baselines; han2025phononbench], a chemistry-conditioned scan loses ... under a same-level DFT hull". This belongs in a supporting section, not the headline.

### 4. Mechanism of CO2 conversion at hydrated HAP

- **Static and experimental work on CO2 at HAP:**
  - static DFT of CO2 and other probe molecules on HAP(0001) (`bittencourt2021happrobe`);
  - experimental FTIR of CO2 on HAP (`cheng1998co2`, already registered).
- **DFT-level carbonation on another mineral class:** AIMD of water/scCO2/silicate interfaces, with carbonate formation near cation vacancies (`lee2015scco2`).
- **Closest framing:** `zhu2025portlandite` maps CO2 mineralization steps on portlandite under a nanometric water film. It uses reactive and non-reactive simulations with well-tempered metadynamics; the abstract does not name the force field.
- **Methods already published:**
  - metadynamics-driven active learning of reactive neural-network potentials (`yang2022metadynamics`);
  - a DFT-level carbonate potential that captures bicarbonate chemistry (`piaggi2025caco3`).
- **Not found:** an AIMD or MLIP study of CO2 conversion at a hydrated HAP surface.
- **Honest phrasing:** the novelty would be the system and a DFT-level free-energy mechanism, for example HAP facets and substitutions. It would not be the method. `zhu2025portlandite` is the natural comparison (portlandite against HAP).

## Consequences for the plan

1. Before any mechanism claim, add biased exploration (metadynamics along C–O coordination or proton-transfer coordinates) to the active learning, following `yang2022metadynamics`. Unbiased exploration has not produced a single CO2 conversion.
2. Report the sealed test per system. The water + CO2 interface (energy RMSE 2.8–3.9 meV/atom in iteration 2) fails the 2 meV/atom gate on its own, although the overall gate passes.
3. The defensible headline is a validated DFT-level potential for hydrated HAP–CO2 interfaces that also covers carbonation products, together with its first mechanistic result. Claim 1 alone is incremental after `wang2023hapmlp` and `wang2024hapmethane`.
4. Small-cell MatterGen phases cannot represent A/B-type carbonated apatite. Bulk carbonated-apatite DFT papers were deferred (see `search_log.md`), and they become relevant if that enumeration starts.
5. Before writing "first", run a database search: Web of Science or Scopus, Google Scholar forward citations of `wang2023hapmlp`, ChemRxiv, and CNKI.

## Note on visibility

Web search returned this project's GitHub repository, `TANG-LAB-WHU/HAP_CO2mineralization`, with its README text. The README describes the MatterGen + DeepMD/MACE plan for HAP CO2 mineralization, so it appears to be publicly indexed. Its visibility was not checked (no `gh` CLI on this machine). If the repository is public, pushed branches expose the unpublished workflow.
