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
