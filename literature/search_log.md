# Literature search log

## 2026-09-03 — seed literature for manuscript structure

Scope: mineral-carbonation mechanisms and the core software/methods already named by the repository.

| Topic | Source selected | Identifier | Screening decision |
|---|---|---|---|
| Interfacial CO₂ mineralization | Qomi et al., *Nature Reviews Chemistry* (2022) | 10.1038/s41570-022-00418-1 | Include as mechanism review; not HAP-specific evidence |
| Generative inorganic materials | Zeni et al., *Nature* (2025) | 10.1038/s41586-025-08628-5 | Include for MatterGen method |
| Equivariant MLIP | Batatia et al., NeurIPS (2022) | arXiv:2206.07697 | Include for MACE architecture |
| Deep Potential | Zhang et al., *Physical Review Letters* (2018) | 10.1103/PhysRevLett.120.143001 | Include for Deep Potential method |
| First-principles MD | Kühne et al., *Journal of Chemical Physics* (2020) | 10.1063/5.0007045 | Include for CP2K |
| Surface construction | Ong et al., *Computational Materials Science* (2013) | 10.1016/j.commatsci.2012.10.028 | Include for pymatgen |
| Classical MD engine | Plimpton, *Journal of Computational Physics* (1995) | 10.1006/jcph.1995.1039 | Include for LAMMPS |
| Enhanced sampling | Bonomi et al., *Nature Methods* (2019) | 10.1038/s41592-019-0506-8 | Include for PLUMED reproducibility |

## Next focused searches

- Hydroxyapatite surface terminations and hydration by facet.
- CO₂ adsorption, carbonate substitution, and carbonate formation on calcium-phosphate surfaces.
- Dopant effects in apatites: F, Cl, Sr, and Mg.
- Transferability and uncertainty of universal MLIPs for ionic aqueous interfaces.
- Collective variables and convergence tests for interfacial carbonate formation.

For each search, record databases, full query strings, date range, inclusion criteria, and rejected near-neighbor literature.

## 2026-09-03 — Phase 2A seed metadata verification

Sources checked:

- Crossref REST API (`api.crossref.org/works/{doi}`) for DOI registration metadata.
- DOI resolution and the corresponding Nature, APS, AIP, Elsevier, and official NeurIPS landing pages.
- Publisher citation metadata embedded in the Nature and APS pages when available.

Outcome: all eight seeds moved from `candidate` to `metadata_partial`. This is a machine check only; `verified_by` and `verified_at` remain blank pending human review. No claim was promoted from `pending`, because this pass did not establish a legal full-text, claim-specific locator. The MACE seed gained the DOI `10.52202/068431-0830` from the official NeurIPS proceedings page and Crossref. The PLUMED creator requires human review because Nature presents the group author as “The PLUMED consortium,” its embedded data lists consortium contributors, and the Crossref creator record is empty.

## 2026-09-03 — Phase 2A outline-gap search, round 1

Databases and services:

- Codex web search (`web.run` with `search_query`), restricted during screening to publisher, DOI, official proceedings, and open full-text landing pages. The underlying index provider was not exposed by the interface.
- Crossref REST API for DOI-based de-duplication and candidate metadata preview; this preview is not bibliographic verification.

### G01 — HAP facets, terminations, and hydration

Exact queries:

1. `hydroxyapatite surface termination hydration first principles DOI`
2. `hydroxyapatite water adsorption surface DFT DOI`

Inclusion rule: retain 2–4 primary studies directly treating HAP surface structure, termination, or water adsorption at identified facets. Exclude biomolecule adsorption, collagen-interface studies, non-HAP apatites, and papers without a sufficiently direct surface/hydration scope.

Selected candidates:

| Citation key | DOI | Screening reason |
|---|---|---|
| `astala2008hapwater` | 10.1103/PhysRevB.78.075427 | Primary first-principles study directly pairing HAP surfaces with water adsorption. |
| `corno2009water` | 10.1021/la803253k | Primary periodic-DFT study of water on named HAP facets. |
| `wang2018hapsurface` | 10.1039/C7RA13121F | Primary first-principles study focused on hydroxyl-dependent HAP surface structure. |
| `peccati2018carbonated` | 10.1021/acs.jpcc.7b12738 | Primary DFT comparison of water reactivity on HAP and carbonated-apatite surfaces. |

### G02 — CO2/carbonate interaction with HAP or calcium-phosphate interfaces

Exact queries:

1. `CO2 adsorption hydroxyapatite surface DOI`
2. `carbonate formation hydroxyapatite calcium phosphate CO2 DOI`
3. `site:pubs.acs.org hydroxyapatite carbon dioxide adsorption`
4. `site:sciencedirect.com hydroxyapatite CO2 adsorption`
5. `hydroxyapatite CO2 capture adsorption DOI`
6. `apatite surface carbon dioxide adsorption first principles`

Inclusion rule: retain 2–4 primary studies directly involving CO2 adsorption, carbonation, or carbonate response in hydroxyapatite/apatite. High-temperature capture studies may be retained only as boundary candidates and cannot support hydrated-interface claims without full-text review.

Selected candidates:

| Citation key | DOI | Screening reason |
|---|---|---|
| `cheng1998co2` | 10.1021/la980339n | Primary FTIR study directly addressing CO2 adsorption on nonstoichiometric calcium HAP. |
| `bouharras2025apatiteco2` | 10.1016/j.jece.2025.115450 | Primary experimental/computational study of CO2 adsorption across apatite materials. |
| `nowicki2024capture` | 10.1039/D3MA00909B | Primary carbonation–regeneration study retained as a high-temperature comparison boundary. |
| `mekhemer2019co2` | 10.1016/j.matchemphys.2018.09.007 | Primary spectroscopy/microscopy study of CO2 interfacial interaction with apatites. |

### De-duplication and exclusions

- Eight unique new candidate DOIs were retained; none duplicates the eight seed DOIs.
- `astala2008hapwater` appeared in more than one result set and was collapsed to one registry row by DOI.
- Title normalization plus first-author/year comparison found no additional duplicate among the selected candidates.
- Excluded near-neighbors included HAP–collagen water adsorption (interface scope differs), fluorapatite-only hydration, HAP gas-sensor studies, transition-metal photocatalytic CO2 reduction, and catalyst-conversion studies whose main endpoint was not adsorption/carbonation of the scoped HAP interface.
- Item-level exclusion records were not retained in the initial search. The category-level summary above is the only surviving exclusion information and must not be reconstructed or treated as an item-level audit trail.
- The first-round cap of eight new candidates is reached. No additional search is authorized until these candidates are screened or a new outline gap is approved.

### Required screening-record template for future searches

Record one row for every screened item, including exclusions, using exactly these fields:

```csv
query_id,title,doi_or_url,decision,reason,reviewed_at
```

## 2026-09-05 — Phase 2B-1 full-text acquisition

This was an acquisition pass for the four references selected in `relevance.csv`, not a new literature search. No new candidate was added.

| Citation key | Landing page | Requested/final PDF URL | Outcome |
|---|---|---|---|
| `qomi2022mineralization` | `https://www.osti.gov/biblio/1889464` | `https://www.osti.gov/servlets/purl/1889464` | Valid 47-page PDF, SHA-256 `6974d5c3caaf6719fe63e63ce28b02866c6836d1bf4306f61358e89e3c0daeb2`. OSTI labels the file “Accepted Manuscript (DOE)” and associates the record with DOI `10.1038/s41570-022-00418-1`; the file is not treated as the version of record. |
| `astala2008hapwater` | `https://journals.aps.org/prb/abstract/10.1103/PhysRevB.78.075427` | `https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.78.075427/fulltext` | Valid 11-page APS PDF, SHA-256 `fd022f917fe8043836ff8ebe11c6493d27fc147127e28899d4a78e972684ccb4`; title, authors, and DOI appear in the PDF. |
| `nowicki2024capture` | `https://pubs.rsc.org/en/content/articlelanding/2024/ma/d3ma00909b` | `https://pubs.rsc.org/en/content/articlepdf/2024/ma/d3ma00909b` | `blocked_fulltext`: three fail-closed attempts returned HTTP 403 and `text/html`; no PDF bytes were retained and no locator was created. No mirror was attempted. |
| `kuhne2020cp2k` | `https://arxiv.org/abs/2003.03868` | `https://arxiv.org/pdf/2003.03868v2` | Valid 52-page arXiv v2 preprint, SHA-256 `183c9526232048e4f578a5fe08a5a35cf4fc6eb989b00c3a71239df780d57627`. Official arXiv metadata associates the preprint with DOI `10.1063/5.0007045`; it is not treated as the publisher version. |

All successful files passed `%PDF-`, `pdfinfo`, non-empty `pdftotext`, publication-identity, SHA-256, and locator-page rendering checks. These machine records do not promote an evidence relationship or activate a citation.

## 2026-11-03 — novelty check for candidate contributions

Authorization: Shawn asked for this check in chat on 2026-11-03. It is a new round, separate from the capped round 1. Round cap: 12 new candidates (about 3 per claim). Report: `novelty_check_2026-11.md`.

Databases and services:

- Claude Code WebSearch, standard and extended modes. The index provider is not exposed.
- WebFetch of the DeePMD-kit publication list (`https://blogs.deepmodeling.com/papers/deepmd-kit/`), read in full in three parts.
- WebFetch of the RePEc IDEAS abstract page for `zhu2025portlandite`; nature.com redirected to a sign-in page.
- Crossref REST API and the arXiv API, for metadata preview only. This is not bibliographic verification, and no publisher page was checked.

Exact queries (E = extended mode):

| query_id | Query |
|---|---|
| Q01 | `deep potential machine learning potential hydroxyapatite water interface` (E) |
| Q02 | `machine learning interatomic potential calcium phosphate apatite molecular dynamics DFT accuracy` |
| Q03 | `hydroxyapatite CO2 adsorption carbonate formation ab initio molecular dynamics surface water` |
| Q04 | `"hydroxyapatite" "neural network potential" OR "machine learning potential" OR "deep potential" molecular dynamics` (E) |
| Q05 | `carbonated apatite machine learning interatomic potential carbonate substitution A-type B-type` (E) |
| Q06 | `Ab initio machine-learning simulation of calcium carbonate from aqueous solutions to the solid state PNAS deep potential` |
| Q07 | `"Investigating the Hydroxyl Reorientation in Hydroxyapatite Using Machine Learning Potentials"` |
| Q08 | `machine learning potential hydroxyapatite vacancies sintering H2O OH O2- molecular dynamics` |
| Q09 | `hydroxyapatite CO2 machine learning potential molecular dynamics adsorption surface 2024 2025` (E) |
| Q10 | `"Machine learning insights into calcium phosphate nucleation and aggregation" Acta Biomaterialia` |
| Q11 | `generative model MatterGen candidates DFT validation phonon dynamical stability S.U.N. overestimate stability benchmark` (E) |
| Q12 | `generative crystal structures active learning machine learning potential training set reaction products workflow mineral` |
| Q13 | `CO2 hydroxyapatite surface DFT carbonate bicarbonate formation (001) (010) water periodic` (E) |
| Q14 | `deep potential OR machine learning potential mineral carbonation CO2 portlandite OR brucite OR wollastonite OR CaO water interface reaction` (E) |
| Q15 | `metadynamics active learning neural network potential reactive events urea decomposition water Parrinello` |
| Q16 | `"Establishing baselines for generative discovery of inorganic crystals" journal` |
| Q17 | `hydroxyapatite CO2 adsorption density functional theory surface hydroxyl carbonate formation calculation` (E) |
| Q18 | `Ni-doped hydroxyapatite neuroevolution potential NEP molecular dynamics methane C–H activation Journal of Catalysis` |
| Q19 | `apatite OR hydroxyapatite CO2 capture carbonation atomistic simulation machine-learned potential MACE OR NequIP OR DeePMD 2025 2026` (E) |
| Q20 | `Zhu Pellenq "CO2 mineralization on wetting nanoscale surfaces" portlandite metadynamics ReaxFF thin water film findings` |
| F01 | DeePMD-kit publication list, scanned for apatite, calcium phosphate, carbonate, CO2 and mineral–water titles |

Inclusion rule: a primary study, benchmark or method paper that is the closest prior work for one of the four candidate contributions. Reviews, MOF and molten-salt studies, and items whose identity could not be retrieved were excluded.

Screening record (one row per item screened; the included rows are the 12 new `candidate` registry rows):

```csv
query_id,title,doi_or_url,decision,reason,reviewed_at
Q07|F01,Investigating the Hydroxyl Reorientation in Hydroxyapatite Using Machine Learning Potentials,10.1021/acs.jpcc.3c02426,include,wang2023hapmlp: prior HAP MLIP (claim 1),2026-11-03
F01|Q04,Selective activation of methane on hydroxyapatite surfaces: Insights from machine learning and density functional theory,10.1016/j.nanoen.2024.109762,include,wang2024hapmethane: prior MLIP-assisted HAP surface reaction (claim 1),2026-11-03
Q10,Machine learning insights into calcium phosphate nucleation and aggregation,10.1016/j.actbio.2025.02.036,include,wang2025capnucleation: aqueous calcium phosphate MLIP (claim 1),2026-11-03
Q06,Ab initio machine-learning simulation of calcium carbonate from aqueous solutions to the solid state,10.1073/pnas.2415663122,include,piaggi2025caco3: pipeline accuracy reference and carbonate DP (claims 1 and 4),2026-11-03
Q14|Q20,Molecular mechanisms of CO2 mineralization on wetting nanoscale surfaces using molecular simulations and metadynamics,10.1038/s41467-025-65794-w,include,zhu2025portlandite: closest mechanism framing (claim 4),2026-11-03
Q15,Using metadynamics to build neural network potentials for reactive events: the case of urea decomposition in water,10.1016/j.cattod.2021.03.018,include,yang2022metadynamics: biased active-learning method (claim 4),2026-11-03
Q11|Q16,Establishing baselines for generative discovery of inorganic crystals,10.1039/D5MH00010F,include,szymanski2025baselines: S.U.N. rates versus stability (claim 3),2026-11-03
Q11,PhononBench:A Large-Scale Phonon-Based Benchmark for Dynamical Stability in Crystal Generation,https://arxiv.org/abs/2512.21227,include,han2025phononbench: dynamical stability of generated crystals (claim 3; preprint),2026-11-03
Q17,Acid-base properties of hydroxyapatite(0001) by the adsorption of probe molecules: An ab initio investigation,10.1103/PhysRevMaterials.5.075003,include,bittencourt2021happrobe: static DFT of CO2 on HAP (claim 4),2026-11-03
F01,An ab-initio deep neural network potential for accurate large-scale simulations of the muscovite mica-water interface,10.1080/00268976.2024.2365430,include,raman2025mica: mineral-water DP context (claim 1),2026-11-03
Q03,Structure dynamics and stability of water/scCO2/mineral interfaces from ab initio molecular dynamics simulations,10.1038/srep14857,include,lee2015scco2: AIMD carbonation on a silicate (claim 4),2026-11-03
Q12,Accelerating crystal structure prediction by machine-learning interatomic potentials with active learning,10.1103/PhysRevB.99.064114,include,podryabinkin2019csp: conceptual precedent (claim 2),2026-11-03
F01,Pressure-induced one-dimensional oxygen ion diffusion channel in lead apatite,10.1103/PhysRevB.108.224114,defer,second apatite DP precedent (claim 1); round cap reached; named in the report,2026-11-03
F01,Unlocking the mysterious polytypic features within vaterite CaCO3,10.1038/s41467-023-43625-0,defer,may help identify the s3_055 CaCO3 polymorph; not a novelty precedent,2026-11-03
Q05,Ab initio simulation on the crystal structure and elastic properties of carbonated apatite,https://www.sciencedirect.com/science/article/abs/pii/S1751616113002075,defer,bulk carbonated-apatite DFT; relevant to A/B substitution enumeration later,2026-11-03
Q05,Location of Carbonate Ions in Metal-Doped Carbonated Hydroxylapatites,10.3390/min13101272,defer,same reason,2026-11-03
Q13,The Mutual Incorporation of Mg2+ and CO32- into Hydroxyapatite: A DFT Study,https://pmc.ncbi.nlm.nih.gov/articles/PMC9781354/,defer,same reason,2026-11-03
Q01,Accurate Force Field Parameters and pH Resolved Surface Models for Hydroxyapatite,https://arxiv.org/pdf/1512.00122,exclude,classical force field,2026-11-03
Q01,Computational study of interfaces between hydroxyapatite and water,https://pubs.rsc.org/en/content/articlelanding/2003/cp/b306358e,exclude,classical force field,2026-11-03
Q01,Perspective: Atomistic Simulations of Water and Aqueous Systems with Machine Learning Potentials,https://arxiv.org/pdf/2401.17875,exclude,general review,2026-11-03
Q03,CO2 adsorption mechanisms in hydrated silica nanopores,https://arxiv.org/pdf/2502.07251,exclude,silica; lee2015scco2 kept as the AIMD carbonation example,2026-11-03
Q04,Machine learning on properties of multiscale multisource hydroxyapatite nanoparticles datasets with different morphologies and sizes,https://www.nature.com/articles/s41524-021-00618-1,exclude,property regression; not an interatomic potential,2026-11-03
Q04|Q08,HAP MLIP with H2O/OH/O2- vacancies for sintering (described in a search summary),n/a,exclude,title and DOI not retrieved,2026-11-03
Q04|Q18,Ni-doped HAP NEP molecular dynamics for methane activation (described in a search summary),n/a,exclude,not found by a targeted search; unverified,2026-11-03
Q08,Overcoming Inaccuracies in Machine Learning Interatomic Potential Implementation for Ionic Vacancy Simulations,10.1021/acs.jpclett.4c02934,exclude,vacancy methodology outside scope,2026-11-03
Q09,CO2 in MOFs with machine-learning potentials (ChemRxiv 2023-38kgb; ChemRxiv 2025-c85xt; PMC13071827),n/a,exclude,MOFs,2026-11-03
F01,First-Principles-Based Machine Learning Models for Phase Behavior and Transport Properties of CO2,10.1021/acs.jpcb.3c00610,exclude,bulk CO2 fluid,2026-11-03
F01,Molecular Simulation of Lithium Carbonate Reactive Vapor-Liquid Equilibria Using a Deep Potential Model,10.1021/acs.jced.3c00580,exclude,molten carbonate,2026-11-03
F01,Modeling Chemical Reactions in Alkali Carbonate-Hydroxide Electrolytes with Deep Learning Potentials,10.1021/acs.jctc.2c00816,exclude,molten salt,2026-11-03
F01,TiO2-water Deep Potential studies (e.g. 10.1073/pnas.2212250120),n/a,exclude,raman2025mica kept as the single mineral-water DP example,2026-11-03
Q13,Density Functional Theory Study of Surface Carbonate Formation on BaO(001),10.1021/jp907765e,exclude,other oxide,2026-11-03
Q13,CO2 on anatase TiO2(001),https://arxiv.org/pdf/1811.00341,exclude,other oxide,2026-11-03
Q12,Machine-Learned Potentials by Active Learning from Organic Crystal Structure Prediction Landscapes,https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10860135/,exclude,organic crystals; podryabinkin2019csp kept,2026-11-03
Q11,A generative model for inorganic materials design,10.1038/s41586-025-08628-5,duplicate,already registered as zeni2025mattergen,2026-11-03
Q20,Molecular-scale mechanisms of CO2 mineralization in nanoscale interfacial water films,10.1038/s41570-022-00418-1,duplicate,already registered as qomi2022mineralization,2026-11-03
Q19,TANG-LAB-WHU/HAP_CO2mineralization (GitHub),https://github.com/TANG-LAB-WHU/HAP_CO2mineralization,exclude,this project's own repository,2026-11-03
```

De-duplication: the 12 included DOIs are new to `references.csv`; the checker found no duplicate DOI. Commas inside titles were dropped in the screening CSV and multiple query ids are joined with `|`.

Boundary: all 12 rows are `candidate` with `doi_status=pending`. No evidence relationship was created or changed, and the cap of 12 is reached.
