# J003 — MnSe2/WTe2 Heterostructure Benchmark

## 1. Objective

J003 investigates the structural energetics and magnetic response of a
MnSe2/WTe2 heterostructure. The benchmark includes:

1. Six interlayer stacking configurations.
2. Isolated-layer reference calculations for binding energy.
3. A MnSe2 lattice-constant scan.
4. A rigid interlayer-distance scan for the lowest-energy relaxed stacking.

These calculations provide a structural and energetic baseline for subsequent
electronic-structure and valleytronics investigations.

## 2. Computational setup

| Parameter | Setting |
|---|---|
| Electronic-structure code | VASP 6.3.2 |
| Exchange-correlation functional | PBE |
| Dispersion correction | DFT-D3(BJ), `IVDW = 12` |
| Plane-wave cutoff | 500 eV |
| Electronic convergence | `EDIFF = 1E-6` eV |
| k-point mesh | Gamma-centered 12 x 12 x 1 |
| Smearing | Gaussian, `ISMEAR = 0`, `SIGMA = 0.05` eV |
| Spin treatment | Spin-polarized, `ISPIN = 2` |
| PAW datasets | W_sv, Te, Mn_pv, Se |
| Dipole correction | Enabled along z |
| Hubbard U | U = 0 eV |
| Ionic relaxation criterion | `EDIFFG = -0.01` eV/Angstrom |
| Heterostructure in-plane lattice constant | 3.5597 Angstrom |
| Out-of-plane cell length | 30 Angstrom |

The heterostructure calculations used the species order W, Te, Mn, Se.
The isolated-layer reference calculations used matching computational settings.

## 3. Stacking comparison

The six stacking configurations were relaxed. Their total energies are
reported relative to the lowest-energy relaxed configuration.

| Stacking | Registry | Energy (eV/cell) | Relative energy (meV/cell) | Interlayer distance (Angstrom) |
|---|---|---:|---:|---:|
| SeW_MnH | Parallel | -39.085657 | 0.0 | 2.978 |
| SeH_MnTe | Parallel | -39.080504 | 5.2 | 3.011 |
| SeW_MnTe | Rotated 60 degrees | -39.075486 | 10.2 | 3.064 |
| SeH_MnW | Rotated 60 degrees | -39.067214 | 18.4 | 3.019 |
| SeTe_MnH | Rotated 60 degrees | -38.964854 | 120.8 | 3.679 |
| SeTe_MnW | Parallel | -38.959254 | 126.4 | 3.699 |

**Result:** SeW_MnH is the lowest-energy configuration among the six
relaxed structures. SeH_MnTe is 5.2 meV/cell higher. The two SeTe-based
configurations have substantially larger separations and higher energies.

The close energy differences among the first four configurations should be
interpreted in the context of numerical convergence and the limitations of
the computational model.

## 4. Binding energy

The isolated-layer reference energies are:

| Reference | Total energy (eV/cell) |
|---|---:|
| MnSe2 | -18.16690140 |
| WTe2 | -20.56735407 |
| SeW_MnH heterostructure | -39.08565685 |

The binding energy is calculated as

    E_bind = E_heterostructure - E_MnSe2 - E_WTe2

    E_bind = -0.35140138 eV/cell

For the six-atom heterostructure cell:

- Binding energy per atom: -58.57 meV/atom.
- Binding energy per in-plane area: approximately -32.02 meV/Angstrom^2.

The negative value indicates that the heterostructure is lower in energy than
the two separately calculated reference layers under this energy convention.
The result depends on the chosen functional, dispersion treatment, geometry,
and consistency of the reference calculations.

## 5. MnSe2 lattice-constant scan

The isolated MnSe2 energy was evaluated at ten in-plane lattice constants.

| a (Angstrom) | Energy (eV/cell) | Relative energy (meV/cell) |
|---:|---:|---:|
| 3.560 | -18.167028 | 6.08 |
| 3.580 | -18.172515 | 0.59 |
| 3.600 | -18.173104 | 0.00 |
| 3.620 | -18.169155 | 3.95 |
| 3.640 | -18.159985 | 13.12 |
| 3.660 | -18.145266 | 27.84 |
| 3.680 | -18.125207 | 47.90 |
| 3.700 | -18.100610 | 72.49 |
| 3.720 | -18.072681 | 100.42 |
| 3.740 | -18.042386 | 130.72 |

A fit to the energy-versus-lattice-constant data gave an estimated
equilibrium lattice constant of 3.5928 Angstrom. The lowest sampled point
was at 3.600 Angstrom.

Using the fitted value, the mismatch relative to the heterostructure lattice
constant of 3.5597 Angstrom is approximately +0.93% for the relaxed MnSe2
lattice relative to WTe2; imposing the WTe2 lattice corresponds to roughly
-0.92% strain in MnSe2.

These mismatch values depend on the definition and reference lattice constants
used. The fitted estimate is a numerical result from this scan, not an
independent experimental measurement.

## 6. Rigid interlayer-distance scan

The lowest-energy relaxed stacking, SeW_MnH, was used as the source geometry.
The MnSe2 layer was shifted rigidly along z while preserving the internal
coordinates of both layers. Static spin-polarized calculations were performed
at nine target separations.

| Target distance (Angstrom) | Total energy (eV/cell) | Relative energy (meV/cell) | Cell magnetization (muB) |
|---:|---:|---:|---:|
| 2.900 | -39.08386796 | 1.637 | 2.8474971 |
| 2.950 | -39.08542619 | 0.079 | 2.8524275 |
| 3.000 | -39.08550538 | 0.000 | 2.8598936 |
| 3.050 | -39.08437739 | 1.128 | 2.8683560 |
| 3.100 | -39.08216852 | 3.337 | 2.8766146 |
| 3.150 | -39.07895566 | 6.550 | 2.8839680 |
| 3.200 | -39.07490892 | 10.596 | 2.8900025 |
| 3.300 | -39.06483307 | 20.672 | 2.8978882 |
| 3.400 | -39.05258157 | 32.924 | 2.9011281 |

The lowest sampled energy occurs at 3.000 Angstrom. The 2.950 Angstrom point
is only 0.079 meV/cell higher, so the data do not establish a sharply resolved
minimum at a precision of 0.05 Angstrom.

The relaxed SeW_MnH geometry has a separation of 2.978 Angstrom and an energy
of -39.08565685 eV/cell, approximately 0.151 meV/cell below the rigid
3.000 Angstrom point. This is consistent with a similar preferred separation,
but the comparison is not a fully relaxed energy curve.

The scan is a rigid-geometry scan: the layers were not re-relaxed at each
distance. The minimum is therefore the minimum among the sampled constrained
structures, not necessarily the exact fully relaxed equilibrium separation.

### Magnetic trend

The calculated total cell magnetization increases from approximately
2.8475 muB at 2.900 Angstrom to 2.9011 muB at 3.400 Angstrom. This indicates
a small change in the magnetic solution with layer separation. Total cell
magnetization alone does not identify which atoms cause the change; site-
projected moments should be analyzed before assigning an atom-specific
mechanism.

### Saved scan outputs

The scan data and plots are saved in the working calculation directory:

- `dscan/dscan_results.csv`
- `dscan/energy_vs_distance.png`
- `dscan/magnetization_vs_distance.png`

## 7. Convergence and quality checks

### Electronic convergence

The final OUTCAR files for all six stacking calculations, two reference
calculations, and ten lattice-scan calculations contained the VASP timing
and accounting marker and the electronic-convergence marker
"aborting loop because EDIFF is reached".

All nine distance-scan calculations also contained these markers and a
final total-energy record.

### Ionic convergence

The six stacking calculations, two reference calculations, and ten
lattice-scan calculations were checked for the VASP ionic-convergence
message "reached required accuracy".

### Maximum final force in stacking calculations

| Stacking | Maximum force (eV/Angstrom) |
|---|---:|
| SeH_MnTe | 0.001532 |
| SeH_MnW | 0.005185 |
| SeTe_MnH | 0.008082 |
| SeTe_MnW | 0.005786 |
| SeW_MnH | 0.003118 |
| SeW_MnTe | 0.008654 |

All six reported maximum forces are below the 0.01 eV/Angstrom relaxation
threshold.

## 8. Main conclusions

1. SeW_MnH is the lowest-energy stacking among the six relaxed configurations.
2. The isolated MnSe2 energy scan gives a fitted lattice constant of
   approximately 3.5928 Angstrom.
3. The heterostructure binding energy is -0.3514 eV/cell under the stated
   reference-energy convention.
4. The rigid distance scan has its lowest sampled energy at 3.000 Angstrom,
   close to the relaxed SeW_MnH separation of 2.978 Angstrom.
5. The total cell magnetization increases slightly as the layers separate.
6. These structural and energetic calculations do not, by themselves,
   establish valley splitting or a valleytronics effect.

## 9. Limitations and next steps

- The very small energy difference between the 2.950 and 3.000 Angstrom
  scan points warrants caution when interpreting the precise minimum.
- The distance scan constrains the layer geometry and does not allow ionic
  relaxation at each point.
- The binding-energy result should be interpreted with attention to reference
  geometry and computational consistency.
- Before making claims about valley physics, calculate and inspect the
  relevant band structures and valley-resolved quantities using a validated
  electronic-structure workflow.
- Preserve raw VASP outputs on Anvil. Track the report, summary CSV files,
  and plots in version control; avoid committing large raw OUTCAR/WAVECAR files
  unless explicitly required.
