# HAP_CO2mineralization

Combining **AIMD** (ab initio molecular dynamics) and **machine learning** workflows to study artificially accelerated mineralization (AAM–CO₂) of CO₂ on nanocrystalline hydroxyapatite (HAP), amorphous hydroxyapatite, or their hybrids. The project builds a strategic framework of Ca–O–P–H nanomaterial configurations for **CO₂ capture, utilization, and storage (CCUS)**.

---

## Project structure

| Directory | Description |
|-----------|-------------|
| **0.InitialStructureConfig** | Initial HAP structures: PDB (from Materials Studio) and XYZ (from VMD) for multiple Miller indices (e.g. 002, 004, 100–513). |
| **1.GeoOpt** | Geometry optimization with CP2K for each HAP facet (`HAP_xxx_Perfect`), including inputs, Docker setup, and optimized structures. |
| **2.ML_AIMD** | ML–AIMD pipeline: surface/LAMMPS setup, CP2K short MD, DeepMD dataset generation and training, LAMMPS runs with the DeepMD potential. Per-facet folders plus `HAP_hkl_template` with workflow docs and scripts. |
| **3.Data_PostProcess** | Post-processing of LAMMPS trajectories: thermodynamics, chemical/structural analysis, RDF, dynamics, quality control, and figure generation. |
| **4.Figure_Making** | Scripts and assets for final figures and visualizations. |

---

## Workflow overview

1. **Initial structures** — Export HAP slabs (and optional CO₂) from PDB/XYZ for chosen Miller indices.
2. **Geometry optimization** — CP2K optimization of each HAP surface (e.g. `run_cp2k_opt.bat`, Docker).
3. **Surface & LAMMPS prep** — Build HAP slabs (ASE/pymatgen), optional CO₂, and generate LAMMPS data files.
4. **CP2K short MD** — Ab initio MD to produce trajectories, energies, forces, and virials.
5. **DeepMD dataset & training** — Convert CP2K outputs to DeepMD format, train a potential (e.g. DPA-1), freeze to a PyTorch model.
6. **LAMMPS with DeepMD** — Run extended MD using the trained potential.
7. **Post-processing** — Thermodynamics, RDF, chemical/structural and dynamics analysis, and plotting (e.g. `lammps_postprocessing.py` in `3.Data_PostProcess`).

Detailed steps, file names, and renaming rules for different Miller indices are in `2.ML_AIMD/HAP_hkl_template/workflow_documentation.md`.

---

## Requirements and usage

- **CP2K** for DFT and AIMD (e.g. via Docker: `docker-compose-cp2k.yml`).
- **LAMMPS** with DeepMD/DPA support for ML-driven MD.
- **DeepMD-kit** for training and freezing the potential (e.g. `deepmodeling/deepmd-kit` image).
- **Python 3.10+** for structure generation, DeepMD data preparation, and post-processing (ASE, pymatgen, etc.; see per-folder `requirements.txt`).

Entry point for the template workflow: `2.ML_AIMD/HAP_hkl_template/run_workflow.bat`. Configure surface (Miller indices, layers, supercell, vacuum), CO₂ (number, height), and run; logs are written under `logs/`.

---

## License

See [LICENSE](LICENSE) in the repository root.
