# Reproduce the plots and inspect the input records

The bundled tables are sufficient to redraw the diagnostic plots without CP2K, LAMMPS, DeepMD, a GPU, or external trajectory files. No simulation or training is part of the plotting workflow.

From this directory, use a Python 3 environment with the plotting dependencies listed in the script:

```sh
python scripts/plot_results.py --verify-only
python scripts/plot_results.py --output-dir figures
```

The first command validates the bundled plotting inputs. The second redraws the figures from the existing compact tables. To preserve the committed figures, choose another output directory. Plotting does not recreate native simulations or independently validate the reported scientific results.

## Provenance and transformations

Read `provenance.json` before using individual files. It links the portable source identifiers to original SHA-256 hashes and records how each public artifact was copied or extracted. Source coordinates and force arrays are retained as recorded; snapshot serialization is not geometry optimization. The force comparison uses one-based atom IDs, forces in eV/Å, and the recorded Hartree/bohr values where supplied. DFT force-unit conversion was performed once in the original comparison; plotting does not apply another conversion.

Machine-specific absolute paths were replaced by portable identifiers or relative paths. The public input files therefore differ from the original inputs byte for byte. Original hashes identify the original records; they must not be interpreted as hashes of the transformed public files. Consult the transformation records and the publication manifest for the corresponding public artifacts.

## Scope of the archived simulation inputs

`inputs/cp2k/` contains `B182-static.inp`, `clean176-static.inp`, `B182-dft4fs.inp`, and `B182-step340.inp`. Paths to method data and any restart assets are represented by relative `data/` and `seed/` locations. Required external files, recorded hashes, seeds, cell settings, method settings, and run lengths are listed in `inputs/run-parameters.json`. Supply and verify those assets before considering a separate calculation. This bundle does not include large wavefunction files or claim bitwise reproduction from inputs alone.

`inputs/lammps/` contains `run0.in`, `md-original-control.in`, `system.data`, and `md.system.data`. These are input records, not a turnkey production workflow. The original control file retains generic runner paths and the known broken stop hook described in the README. The external DeepMD model must match SHA-256 `57c4fc2d7b6890dccdb1346c9d5641f18047feeeb3800df94028c0bbd993fe54`, and its failed qualification remains a blocker to production use. Fixing the stop hook would repair an engineering defect; it would not qualify the model.

Any future simulation is a separate calculation requiring a complete environment and explicit scientific acceptance criteria. In particular, preserve the distinction between the 4 fs direct-DFT trajectory and the 85 fs ML trajectory, and preserve two CO₂ molecules when repeating B182 at fixed coverage. Lateral replication alone does not resolve slab thickness, vacuum, finite-size, sampling, or potential-transferability questions.

All artifacts in this bundle remain manuscript-ineligible. Promotion requires a separate provenance and eligibility decision under the repository's research contract.
