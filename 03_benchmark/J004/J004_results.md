# J004 — SOC Valley-Resolved Electronic Structure of MnTe/ WTe2-Based Heterostructures

## 1. Objective

J004 investigates the spin-orbit-coupled electronic structure and valley-dependent band-edge behaviour of MnTe/WTe2-based heterostructure configurations. The calculations examine opposite magnetic initialisations, a clean magnetic initialisation, changes in Hubbard U treatment, interlayer spacing, and selected composition variants.

The principal quantities extracted include magnetic moments, K- and K'-resolved valence/conduction edges, valley energy differences, interlayer separation, and parser-derived hybridisation indicators.

## 2. Computational setup

The inspected J004 calculations use:

| Parameter | Setting |
|---|---|
| Plane-wave cutoff | 500 eV |
| Electronic convergence, EDIFF | 1 × 10^-8 eV |
| k-point mesh | Gamma-centred 12 × 12 × 1 |
| Bands | 72 |
| Spin treatment | ISPIN = 2 |
| Spin-orbit coupling | LSORBIT = .TRUE. |
| Symmetry | ISYM = -1 |
| Dispersion correction | IVDW = 11 |
| Ionic steps | NSW = 0 |
| Ionic algorithm | IBRION = -1 |
| Hubbard correction | Dudarev, LDAUTYPE = 2, where enabled |

For the standard MnTe/WTe2 setup, the reported Hubbard parameters are U(Mn) = 3.90 eV and U(W) = 2.87 eV. The U0 variant has no LDAU enabled. The UW0 variant retains U(Mn) = 3.90 eV and sets U(W) = 0 eV.

The calculations are static single-point calculations. They do not constitute ionic relaxation or a geometry-convergence study.

## 3. Magnetic initialisation and consistency

Three calculations were compared for the SeH_MnTe configuration:

| Run | Initialisation | Total-cell Mz from OSZICAR (μB) | Total energy (eV) |
|---|---|---:|---:|
| mz | Opposite-z initialisation | -3.1395 | -33.82502404 |
| pz | Positive-z initialisation | +3.1395 | -33.82502416 |
| pz_clean | Mn-only positive-z initialisation | +3.1395 | -33.82502397 |

The three total energies differ by less than 0.0002 meV at the reported precision. The two opposite magnetic initialisations converge to opposite total-cell z magnetisations with essentially identical energies. The clean-start calculation reproduces the positive-z solution.

The OSZICAR total-cell magnetisation and the OUTCAR site-projected magnetisation sum are distinct quantities. For the standard pair, OSZICAR reports approximately ±3.1395 μB, whereas the final OUTCAR site-projected sum is approximately ±3.264 μB. The CSV `mag_cell_z` field corresponds to the OSZICAR value. These values should not be interchanged.

This comparison supports consistency with respect to the tested magnetic starting configurations. It does not by itself establish the physical magnetic ground state.

## 4. Valley-resolved results

The following values are parser-reported. Energies and valley differences are taken from `j004_summary.csv`.

| Configuration | Interlayer distance (Å) | ΔVB (meV) | ΔCB (meV) | ΔK-K' (meV) | Gap K (eV) | Gap K' (eV) | Parser-selected majority valley |
|---|---:|---:|---:|---:|---:|---:|---|
| SeH_MnTe/mz | 2.9748 | +154.97 | -52.86 | +207.82 | 0.4151 | 0.6229 | K |
| SeH_MnTe/pz | 2.9748 | -154.97 | +52.86 | -207.82 | 0.6229 | 0.4151 | K' |
| SeH_MnTe/pz_clean | 2.9748 | -154.97 | +52.86 | -207.83 | 0.6229 | 0.4151 | K' |
| SeH_MnTe_U0/pz_clean | 2.9748 | +131.92 | -14.34 | +146.25 | 0.7257 | 0.8720 | K' |
| SeH_MnTe_UW0/pz_clean | 2.9748 | -151.52 | +22.54 | -174.06 | 0.7298 | 0.5558 | K' |
| SeH_MnW_UW0/pz_clean | 3.2515 | -119.87 | +14.80 | -134.67 | 0.7428 | 0.6082 | K' |
| SeH_d2.900/pz | 2.9000 | -162.96 | +48.22 | -211.19 | 0.6150 | 0.4038 | K' |
| SeH_d3.050/pz | 3.0500 | -147.80 | +56.18 | -203.98 | 0.6316 | 0.4276 | K' |
| SeH_d3.120/pz | 3.1200 | -141.09 | +50.63 | -191.72 | 0.6315 | 0.4398 | K' |
| SeW_MnH_UW0/pz_clean | 3.1915 | -5.22 | +15.98 | -21.20 | 0.7548 | 0.7336 | K' |
| SeW_MnTe_UW0/pz_clean | 3.1296 | -113.76 | +18.42 | -132.18 | 0.7709 | 0.6388 | K' |

Here ΔVB and ΔCB are the parser's signed K-versus-K' edge differences. ΔK-K' is the parser-reported difference between the corresponding valley gaps. The signs depend on the magnetic orientation and the parser's convention.

### 4.1 Opposite magnetic orientations

The mz and pz calculations show reversal of the reported valley preference when the magnetic orientation is reversed. The magnitude of the reported gap difference is approximately 207.8 meV in both cases.

The clean-start pz calculation reproduces the pz result closely: its reported ΔK-K' differs by approximately 0.01 meV. This is a useful internal reproducibility check.

### 4.2 Hubbard-U sensitivity

The U0 calculation gives a smaller reported ΔK-K' than the standard-U calculation. Removing W U while retaining Mn U also changes the reported valley splitting and gap values.

These results indicate sensitivity to the chosen Hubbard treatment. They should not be interpreted as a converged prediction of the physical U dependence without further parameter justification and band-character validation.

### 4.3 Interlayer-distance scan

For the SeH_d2.900, SeH_MnTe, SeH_d3.050, and SeH_d3.120 configurations, the parser reports a decrease in the magnitude of ΔK-K' as the interlayer distance increases over the sampled range.

This is a trend across the sampled static configurations, not yet a fully converged distance-dependent physical law. The structures, band-edge identities, and interlayer geometry should be checked consistently before drawing stronger conclusions.

### 4.4 Composition variants

The SeH_MnW and SeW_MnH variants show different reported valley-edge differences from the SeH_MnTe reference. In particular, SeW_MnH has a much smaller reported ΔK-K' magnitude in this dataset.

These are candidate composition-dependent trends and require consistent structural and electronic validation.

## 5. Early baseline calculations

The separate early-baseline CSV contains two calculations:

| Run | Interlayer distance (Å) | ΔVB (meV) | ΔCB (meV) | ΔK-K' (meV) | Gap K (eV) | Gap K' (eV) |
|---|---:|---:|---:|---:|---:|---:|
| dos_refine | 3.05 | -1.65 | +8.73 | -10.38 | 0.7179 | 0.7075 |
| soc_U0 | — | -1.66 | +8.71 | -10.37 | 0.7178 | 0.7075 |

These calculations are retained as an early baseline. They are not merged into the main SeH_MnTe trend because their calculation provenance and reference setup differ. Total energies from these runs should not be compared directly with the main run family without confirming identical structures and settings.

## 6. Parser methodology and validation limitations

The analysis script `parse_j004.py` reads the calculation outputs, identifies K/K' points, classifies WTe2-like states, and selects candidate valence and conduction edges. It then calculates valley-dependent gaps and related quantities.

Important limitations:

1. The edge-selection routine ranks classified states and can select adjacent classified valence/conduction states while skipping intervening unclassified states. The selected edges therefore require independent inspection of the nearby bands and their orbital/layer character.
2. Some configurations carry warnings about partly WTe2-like states near the selected edges. Such warnings indicate possible hybridisation and ambiguity in assigning the edge states.
3. Parser status alone is not proof that all band-edge assignments are physically valid.
4. The reported quantities `partner_WTe2_share`, `t_eV`, `offset_eV`, and `push_meV` are parser-derived indicators. Their precise definitions and assumptions should be documented and independently checked before physical interpretation.
5. The present comparison is based on selected static calculations. It does not establish convergence with respect to k-mesh, cutoff, vacuum thickness, or all structural degrees of freedom.

Accordingly, the valley splittings in this report should be treated as candidate parser-derived results pending band-identity validation.

## 7. Files and reproducibility

The repository includes:

- `j004_summary.csv` — main J004 parsed results.
- `j004_early_baseline_summary.csv` — separate early-baseline results.
- `parse_j004.py` — analysis script.

Large raw VASP outputs are retained in the calculation area and are not included in the public repository. Licensed POTCAR files must not be committed.

## 8. Conclusions and next steps

J004 establishes a reproducible set of candidate SOC valley-resolved results for several MnTe/WTe2-based configurations. The opposite magnetic initialisations reverse the parser-reported valley preference while retaining essentially identical total energies, and the clean-start calculation reproduces the positive-z result.

The results also suggest sensitivity to Hubbard-U treatment, interlayer distance, and composition. However, the physical interpretation depends on validating the selected band-edge states and their layer/orbital character.

Priority next steps:

1. Inspect the band-resolved PROCAR/EIGENVAL states around both K and K' for the principal SeH_MnTe cases.
2. Verify that the selected valence and conduction edges are the actual adjacent physical band edges.
3. Quantify layer and orbital character of the selected states and nearby competing states.
4. Establish k-mesh and other numerical convergence for the key valley-splitting observable.
5. Only then use the candidate valley splitting for broader physical claims or ML labels.
