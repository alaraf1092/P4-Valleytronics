# J005 — Free-Layer Reference Calculations and Vacuum-Level Alignment

## 1. Objective

J005 calculates isolated MnSe2 and WTe2 reference-layer electronic properties to support the interpretation of MnSe2/WTe2 heterostructure band alignment.

The calculations examine SOC and non-SOC references, Hubbard-U settings, lattice parameters, vacuum levels, work functions, and WTe2 band edges. A refined WTe2 calculation tests the effect of a denser k-point mesh and tighter electronic convergence.

These are free-layer reference calculations. Their vacuum-level alignment is not a direct calculation of the interface dipole or charge transfer in the assembled heterostructure.

## 2. Computational setup

The inspected J005 calculations use:

| Parameter | Setting |
|---|---|
| Plane-wave cutoff | 500 eV |
| Electronic convergence | 1 × 10^-7 eV, or 1 × 10^-8 eV for refined run |
| Standard k-point mesh | Gamma-centred 12 × 12 × 1 |
| Refined k-point mesh | Gamma-centred 18 × 18 × 1 |
| Spin treatment | ISPIN = 2 |
| Spin-orbit coupling | Enabled for SOC runs |
| Symmetry | ISYM = -1 |
| Dispersion correction | IVDW = 11 |
| Ionic steps | NSW = 0 |
| Ionic algorithm | IBRION = -1 |
| Hubbard correction | Dudarev, LDAUTYPE = 2 |

The MnSe2 SOC reference uses U(Mn) = 3.90 eV. The WTe2 SOC reference uses U(W) = 2.87 eV. The inspected SOC runs use `LSORBIT = .TRUE.`.

All listed calculations are static single-point calculations rather than ionic relaxations.

## 3. Vacuum-level and work-function definitions

The parser obtains the vacuum reference from the planar-averaged electrostatic potential in LOCPOT. The work function is calculated as:

WF = V_vac - E_F

where V_vac is the vacuum potential and E_F is the Fermi energy, both expressed in the same energy reference.

For WTe2, the parser also reports K-point valence and conduction energies relative to the vacuum level. These quantities are used to construct a free-layer band-alignment estimate.

The vacuum-level alignment assumes that the isolated-layer vacuum references can be compared. It does not include the interface-induced electrostatic dipole, charge redistribution, or structural changes that can occur after stacking.

## 4. MnSe2 reference results

| Calculation | SOC | U(Mn) (eV) | Lattice parameter (Å) | V_vac (eV) | E_F (eV) | Work function (eV) | Lowest two unoccupied K states (eV) |
|---|---|---:|---:|---:|---:|---:|---:|
| MnSe2 reference | No | 3.90 | 3.5597 | 1.8761 | -3.5510 | 5.4271 | 0.184 |
| MnSe2 SOC | Yes | 3.90 | 3.5597 | 1.8762 | -3.5454 | 5.4216 | 0.143 / 0.201 |
| MnSe2 reference, alternate lattice | No | 3.90 | 3.6400 | 1.7944 | -3.6147 | 5.4091 | 0.250 |
| MnSe2 SOC, alternate lattice | Yes | 3.99 | 3.6400 | 1.7966 | -3.6324 | 5.4290 | 0.216 / 0.278 |

The final column records the parser's lowest two unoccupied K-point states or reported pair of energies. It should not be called a majority-spin doublet for SOC calculations: the parser does not calculate the spin expectation value of each state, and SOC spinors are not generally represented by a simple majority/minority label.

The work functions are close to 5.4 eV across these four reference calculations. The alternate-lattice cases also change the Hubbard parameter slightly, so their differences should not be attributed to lattice parameter alone.

## 5. WTe2 reference results

| Calculation | SOC | U(W) (eV) | Lattice parameter (Å) | V_vac (eV) | E_F (eV) | Work function (eV) | VB(K) vs vacuum (eV) | CB(K) vs vacuum (eV) | Gap at K (eV) | Reported VB pair separation (eV) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| WTe2 reference | No | 0.00 | 3.5597 | 2.8752 | -1.4229 | 4.2981 | -4.672 | -3.601 | 1.071 | — |
| WTe2 reference | No | 2.87 | 3.5597 | 2.8912 | -1.4641 | 4.3553 | -4.624 | -3.580 | 1.044 | — |
| WTe2 SOC | Yes | 0.00 | 3.5597 | 2.8724 | -1.3097 | 4.1821 | -4.395 | -3.648 | 0.747 | 0.486 |
| WTe2 SOC | Yes | 2.87 | 3.5597 | 2.8889 | -1.0621 | 3.9510 | -4.281 | -3.672 | 0.609 | 0.616 |
| WTe2 SOC, alternate lattice | Yes | 2.87 | 3.5400 | 2.9171 | -1.0272 | 3.9443 | -4.266 | -3.598 | 0.668 | 0.610 |
| WTe2 SOC, refined mesh | Yes | 2.87 | 3.5400 | 2.9171 | -1.0645 | 3.9816 | -4.266 | -3.598 | 0.668 | 0.610 |

The reported VB pair separation is the separation of the two highest occupied states selected at the parser's chosen K point. It is not yet independently verified as the desired SOC splitting of a particular orbital- and layer-resolved band.

## 6. SOC and Hubbard-U sensitivity

For WTe2 at a = 3.5597 Å, the parser-reported K-point gap changes from 1.044 eV without SOC at U(W) = 2.87 eV to 0.609 eV with SOC at the same U value.

The reported VB pair separation is 0.616 eV for that SOC calculation. This is a parser-selected state separation and requires inspection of the band identities and orbital character before being interpreted as a specific physical SOC splitting.

The results also show that the calculated work function and band edges depend on the Hubbard-U treatment. These comparisons are useful sensitivity tests, but they do not establish that the chosen U value is uniquely appropriate.

## 7. k-point refinement test

The WTe2 SOC calculation at a = 3.540 Å was repeated with a denser mesh:

| Quantity | 12 × 12 × 1, EDIFF 10^-7 | 18 × 18 × 1, EDIFF 10^-8 |
|---|---:|---:|
| Total energy (eV) | -17.91601938 | -17.91608911 |
| Fermi energy (eV) | -1.0272 | -1.0645 |
| Work function (eV) | 3.9443 | 3.9816 |
| VB(K) vs vacuum (eV) | -4.266 | -4.266 |
| CB(K) vs vacuum (eV) | -3.598 | -3.598 |
| Gap at K (eV) | 0.668 | 0.668 |
| Reported VB pair separation (eV) | 0.610 | 0.610 |

The total-energy difference is approximately 0.0697 meV. The Fermi energy and work function differ by approximately 37 meV, while the printed K-point VB, CB, gap, and pair-separation values remain unchanged to the reported precision.

This supports stability of the listed K-point gap and pair separation at the displayed precision under this particular refinement. It does not establish convergence of every quantity, particularly the Fermi energy and work function. Further convergence tests should focus on the observables used in the final band-alignment argument.

## 8. K-point selection and parser limitations

The parser reads EIGENVAL and identifies a K-corner point using the reciprocal lattice and a geometric criterion. The inspected 12 × 12 × 1 EIGENVAL files contain two opposite corner points, at fractional coordinates approximately (1/3, 1/3, 0) and (-1/3, -1/3, 0).

The current parser selects the first matching corner only. It does not independently calculate and compare both K and K' for every J005 reference.

Consequences:

1. J005 does not establish a complete K-versus-K' comparison for magnetic MnSe2.
2. The reported WTe2 values are selected-K reference quantities, not a two-valley analysis.
3. The SOC state separation requires band- and orbital-character inspection before being assigned to a particular physical splitting.
4. Parser status alone does not guarantee that the EIGENVAL-derived band quantities were successfully and physically interpreted.

The MnSe2 SOC parser output should be described as the lowest two unoccupied K-point states, not as a spin-majority doublet.

## 9. Band alignment interpretation and limitations

The isolated-layer vacuum-level alignment provides a reference estimate for comparing MnSe2 and WTe2 band positions. It is useful for motivating heterostructure calculations and identifying candidate band offsets.

However, this alignment is not a direct prediction of the assembled interface because it does not include:

- Interface dipole formation.
- Interlayer charge transfer.
- Interface-induced hybridisation.
- Structural relaxation after stacking.
- Changes in the electrostatic potential caused by the heterostructure.

The withdrawn +0.66 eV offset should not be used as a J005 result. It combines quantities with different provenance and reference conditions. The free-layer values in this report should be used instead, with the assumptions of vacuum-level alignment stated explicitly.

## 10. Files and reproducibility

The repository includes:

- `j005_summary.csv` — parsed J005 reference results.
- `parse_j005.py` — analysis script.

The original and refined WTe2 calculations are both retained in the calculation area. The refined calculation uses an 18 × 18 × 1 mesh and tighter EDIFF. Large raw VASP outputs are not included in the public repository. Licensed POTCAR files must not be committed.

## 11. Conclusions and next steps

J005 provides free-layer reference calculations for MnSe2 and WTe2, including SOC and non-SOC cases, Hubbard-U variations, and a WTe2 k-point refinement.

The WTe2 K-point gap and reported VB pair separation are unchanged to the displayed precision in the 12 × 12 × 1 versus 18 × 18 × 1 refinement, while the Fermi energy and work function shift by approximately 37 meV. The work-function and band-edge results therefore require observable-specific convergence assessment.

The isolated-layer vacuum alignment is a useful reference for heterostructure design, but it must not be presented as the actual interface band alignment without explicit heterostructure calculations.

Priority next steps:

1. Validate both opposite K corners for the relevant free-layer references.
2. Inspect the orbital and layer character of the selected WTe2 SOC states.
3. Check vacuum plateaus and vacuum-level extraction consistency.
4. Test convergence of work functions and band-edge positions, not only the K-point gap.
5. Compare free-layer alignment with the self-consistent electrostatic potential and band character of the assembled heterostructure.
