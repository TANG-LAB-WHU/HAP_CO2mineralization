# Overleaf project: local HAP--CO2 validation draft

Import `HAP_local_validation_overleaf.zip` as a new Overleaf project. Select `main.tex`
as the main document and XeLaTeX as compiler. Use BibTeX (not Biber); Overleaf's normal
recompile resolves it automatically. All referenced figures and bibliography entries
are bundled. No login, upload, external image URL, custom font, or shell escape was used.

This is a scientific working draft about a bounded local workflow and small-model tests.
It does not establish a transferable reactive potential or CO2 mineralization. Read the
failed gates and limitations before using its results. Authorship, affiliations, funding,
and competing interests must be supplied by the research team. Raw data eligibility in
the original repository was not promoted; figure use is restricted to explicitly labeled
technical/negative-finding diagnostics requested for this draft.

`source_to_claim.json` maps numerical/figure claims to local source paths and checksums.
`verification.json` contains actual gate outcomes. `figure_sources.json` identifies the
plot inputs. `literature/` records primary-source access limits; `references.bib` contains
only entries cited by the manuscript. The separately delivered GIF shows actual saved
LAMMPS frames without interpolation. It is not a mineralization movie.

On the configured Shawn-PC checkout, verify the retained run and rebuild the document using:

```powershell
cd D:\project\HAP_CO2mineralization
.\workflow\run-refinement.ps1 -Stage all -Run refine-001
python workflow/build_validation_paper.py --run .research-runs/refine-001 --recovery .research-runs/recovery-001
```

A captured run is immutable: the all-stage command verifies hashes without rewriting it.
Use a new run ID for new calculations, or a separate review copy for figure edits.
These commands depend on the pinned local runtime and retained raw outputs,
documented in `workflow/REFINEMENT.md`; the Overleaf project itself is self-contained for
document compilation, not a redistribution of model weights or all trajectories.
