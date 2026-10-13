# J002 — MnSe2 Magnetic-Ordering and Hubbard-U Benchmark

**Status:** Six structural relaxations completed; preliminary analysis documented
**Code:** VASP 6.3.2 (`vasp_std`)
**Platform:** Purdue Anvil HPC
**System:** 1T-MnSe2, two Mn atoms and four Se atoms per simulation cell

## 1. Objective

Compare ferromagnetic (FM) and antiferromagnetic (AFM) initial magnetic configurations at three effective Hubbard-U values on Mn-d states:

- U_eff = 0 eV
- U_eff = 2 eV
- U_eff = 4 eV

The purpose is to establish a reproducible baseline before comparing the calculation inputs and magnetic-ordering definitions with the reference study.

## 2. Model and computational settings

| Parameter | Setting |
|---|---|
| Structure label | 1T-MnSe2 |
| Cell composition | Mn2Se4 (6 atoms) |
| Exchange-correlation | PBE, as specified in the project README |
| Dispersion correction | IVDW = 12 (D3-BJ) |
| Plane-wave cutoff | 500 eV |
| Electronic convergence | EDIFF = 1e-6 eV |
| Force convergence | EDIFFG = -0.01 eV/Angstrom |
| Spin treatment | ISPIN = 2, collinear spin |
| Initial Mn moments | FM: +5, +5; AFM: +5, -5 |
| Brillouin-zone sampling | Gamma-centred 6 x 12 x 1 |
| Smearing | ISMEAR = 0, SIGMA = 0.05 eV |
| Ionic relaxation | IBRION = 2, NSW = 60 |
| Cell relaxation | ISIF = 2; cell vectors fixed, ions relaxed |
| Real-space projection | LREAL = .FALSE. |
| PAW datasets | Mn_pv and Se, in POSCAR species order |

For U = 2 and 4 eV, Dudarev DFT+U was enabled using
LDAUTYPE = 2, LDAUL = 2 -1, LDAUU = U 0, and LDAUJ = 0 0.
For U = 0 eV, the LDAU tags were omitted.

The input cell vectors have lengths 7.380, 3.690, and 23.600 Angstrom,
with a non-orthogonal in-plane cell. The in-plane lattice parameter
was seeded at 3.690 Angstrom. Because ISIF = 2 was used, these runs
do not optimize the lattice vectors.

## 3. Relaxation results

Energies below are the final total energies recorded in the VASP outputs.
Define the energy difference consistently as:

Delta E = E_AFM - E_FM

Positive Delta E means the FM solution is lower in energy among the
two converged configurations tested.

| U_eff (eV) | E_FM (eV/cell) | E_AFM (eV/cell) | Delta E (meV/cell) | Lower-energy solution |
|---:|---:|---:|---:|---|
| 0 | -36.22714352 | -36.08664552 | +140.498 | FM |
| 2 | -33.71277818 | -33.51429715 | +198.481 | FM |
| 4 | -31.63802238 | -31.38993675 | +248.086 | FM |

The energy difference increases across the tested U_eff values.
This is an observation for the present computational setup, not a
universal prediction of magnetic ordering.

The energies are total energies of the same six-atom cell within each
FM/AFM pair. They should not be compared directly across different U
values as though they were calculations at the same Hamiltonian.

## 4. Completion checks

All six relaxation outputs were checked for the VASP structural
convergence message and the final timing/accounting section.
The six-run CSV is `j002_summary.csv`.

The parser's `m_Mn_muB` field is NaN for these relaxations because
the runs did not provide the site-resolved magnetization tables needed
by the current parser. Do not interpret NaN as zero magnetic moment.

## 5. Separate U = 0 FM static calculation

A separate static calculation was run in `U0_FM_static` using the
relaxed FM geometry, with LORBIT = 11 to obtain site-projected moments.
It is not included in the six-run relaxation energy comparison.

| Site | Projected total moment (mu_B per atom) |
|---|---:|
| Mn 1 | +3.214 |
| Mn 2 | +3.214 |
| Se 1 | -0.133 |
| Se 2 | -0.133 |
| Se 3 | -0.133 |
| Se 4 | -0.133 |

The OUTCAR projected-table total is +5.895 mu_B. The sum of the
six displayed, rounded site values is +5.896 mu_B, consistent to
the displayed precision. The VASP total cell magnetization is about
+6.000 mu_B; the cell magnetization and atom-projected sum are
different quantities.

These are site-projected moments for this particular U = 0 FM static
calculation, not a U-dependent moment series.

## 6. Reference comparison and limitations

The existing project README records a reference target of approximately
-111 meV for E_AFM - E_FM at zero strain. Our current values are
positive under the explicitly stated convention above. Therefore,
the numerical and sign comparison is not yet resolved.

Before drawing a conclusion, verify against the reference paper and
the forthcoming supervisor inputs:
1. The exact AFM spin pattern and atom mapping.
2. The definition and sign convention used for the reported energy
   difference, including whether energies are normalized per cell or
   per magnetic pair.
3. The structure, lattice constraints, U value, and dispersion method.
4. Whether the reference relaxes the cell differently.

This benchmark uses one chosen AFM pattern and the seeded cell vectors.
It does not establish that all possible magnetic orderings have been
searched, nor does it determine a Curie temperature. No SOC, MAE,
band structure, or valley-splitting calculation is reported here.

## 7. Files and reproducibility

- `j002_summary.csv`: six relaxation energies and metadata.
- `U0_FM_static/site_moments.csv`: site-projected static moments.
- `U0_FM`, `U0_AFM`, `U2_FM`, `U2_AFM`, `U4_FM`, `U4_AFM`:
  relaxation inputs and raw VASP outputs.
- `U0_FM_static`: separate static calculation.
- `../parse_j002.py`: analysis script in the parent J002 directory.
- `../README.md`: original kit description and initial reference target.

Keep the raw OUTCAR files and inputs unchanged. Record any future
changes to the parser or inputs separately so this result remains
reproducible.

## 8. Next steps

1. Compare our inputs against the supervisor's/reference inputs.
2. Resolve the magnetic-ordering pattern and Delta E sign/normalization.
3. Confirm whether lattice relaxation is required.
4. Only after resolving those differences, decide whether additional
   calculations are necessary.
