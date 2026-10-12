# Figure descriptions and scientific limits

1. **Initial and endpoint local structures.** A91 refines the existing adsorbed
   geometry and reaches the common optimizer criteria. B91 evolves from its
   prescribed approach but remains partial after the last continuation. All atom
   positions are measured; periodic neighbors are shown using the recorded cell.
   Distances and angles describe those snapshots only, not a mechanism or a
   vibrationally verified minimum. Atom IDs are one-based.
2. **Segmented B91 optimization.** Energy is relative to B's own initial energy;
   Cartesian maximum/RMS forces are compared with the unchanged thresholds. The
   original segment, engineering recovery and last continuation remain distinct.
   Restart duplicates and the interrupted SCF-only step are not added as new
   geometries. A/clean/gas endpoint status is contextual, not an energy ranking.
3. **Conditional A adsorption subtraction.** The three matched relaxed reference
   energies give -0.346114 eV/CO2 under the stated periodic model. Displayed
   energies are rounded; the calculation uses the recorded unrounded values.
   B is excluded. No BSSE, isolated-gas or vibrational-minimum claim is made.
4. **Observed DFT-to-ML-to-MD validation route.** The campaign-002 v2 training and
   development diagnostic failed their force criteria. Separately authorized
   exploration passed interface/geometry checks but failed the 10-fs CO2 force
   screen; the 20-fs DFT probe was skipped. Arrows indicate the recorded execution
   route, not a passed qualification chain. No production potential was released.

All four figures are development results with retained evidence-review markers.
