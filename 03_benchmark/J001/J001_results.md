# J001 — 2H-WTe2 Monolayer Convergence & Structural Benchmark

**Project:** P4 — ML-Guided Systematic Screening of 2D Heterostructures for Room-Temperature Valley Polarization
**Benchmark:** J001 — 2H-WTe2 monolayer
**Date:** 2026-10-05
**Status:** COMPLETED

---

## 1. Objective

Fix the numerical baseline for P4 before any heterostructure calculation: plane-wave cutoff, k-point sampling and the PBE in-plane lattice constant of the 2H-WTe2 host. All J001 runs are non-magnetic and without spin–orbit coupling.

## 2. System

| Item | Value |
|---|---|
| Structure | 2H-WTe2 monolayer, 3 atoms (1 W, 2 Te) |
| Space group | P-6m2 (#187) |
| Starting cell | a = 3.550 Å (seed), c = 23.60 Å (about 20 Å vacuum), Te–Te thickness 3.60 Å |
| Built with | ASE `mx2`, via `setup_j001.py` |

## 3. Method

| Setting | Value |
|---|---|
| Code | VASP 6.3.2 (`vasp_std`) on Purdue Anvil (GCC 11.2.0, OpenMPI 4.1.6, MKL 2019.5) |
| Functional | PBE (from the POTCARs) |
| Spin / SOC | non-magnetic (`ISPIN = 1`), no SOC |
| Electronic | `PREC = Accurate`, `EDIFF = 1E-6`, `ISMEAR = 0`, `SIGMA = 0.05`, `LREAL = .FALSE.` |
| k-points | Γ-centred N×N×1 |
| E–a relaxations | `IBRION = 2`, `ISIF = 2`, `NSW = 60`, `EDIFFG = -0.01` |
| Resources | `shared` partition, 1 node, 8 MPI ranks, `OMP_NUM_THREADS = 1`, account nnt250008 |
| Convergence target | \|E − E_ref\| < 1 meV/atom against the most precise run, held by every denser setting |

PAW datasets (PBE.54):

```text
PAW_PBE W_sv 04Sep2015   (14 valence electrons)
PAW_PBE Te 08Apr2002     (6 valence electrons)
```

## 4. Plane-wave cutoff (k = 12×12×1, a = 3.55 Å; reference 600 eV)

| ENCUT (eV) | E/atom (eV) | \|ΔE\| vs 600 eV (meV/atom) |
|---:|---:|---:|
| 350 | −6.475002 | 0.348 |
| 400 | −6.475194 | 0.156 |
| 450 | −6.475245 | 0.105 |
| 500 | −6.475279 | 0.071 |
| 550 | −6.475308 | 0.042 |
| 600 | −6.475350 | reference |

Energy decreases monotonically (variational basis); every tested cutoff meets the target.
**Adopted: 500 eV.**

## 5. k-point sampling (ENCUT = 500 eV; reference 21×21×1)

| k-mesh | E/atom (eV) | \|ΔE\| vs 21×21×1 (meV/atom) |
|---:|---:|---:|
| 6×6×1 | −6.474058 | 1.229 (fails) |
| 9×9×1 | −6.475212 | 0.075 |
| 12×12×1 | −6.475279 | 0.008 |
| 15×15×1 | −6.475285 | 0.001 |
| 18×18×1 | −6.475287 | ~0 |
| 21×21×1 | −6.475287 | reference |

Converged from 9×9×1; every denser mesh stays converged. The 500 eV / 12×12×1 point appears in both series with the same energy (−6.475279 eV/atom).
**Adopted: Γ-centred 12×12×1.**

## 6. In-plane lattice constant (E–a scan; 500 eV, 12×12×1)

- 11 cells from 3.45 to 3.65 Å in steps of 0.02 Å, ions relaxed in each.
- Quadratic fit: **a0 = 3.5597 Å (≈ 3.56 Å)**, 0.27 % above the seed.
- The minimum lies inside the scanned range.

## 7. Frozen J001 baseline

| Parameter | Value |
|---|---|
| Functional / PAW | PBE, `W_sv` + `Te` (PBE.54) |
| ENCUT | 500 eV |
| k-mesh (1×1 cell) | Γ-centred 12×12×1 |
| a0 | 3.5597 Å ≈ 3.56 Å |
| SOC / magnetism | off / non-magnetic |

## 8. Computational cost and Slurm provenance

| Calculation | Slurm job | Status | Anvil SUs |
|---|---:|---|---:|
| ENCUT convergence | 21082428 | COMPLETED | 0.2200 |
| Initial duplicate k-point submission | 21084888 | CANCELLED | 0.1664 |
| Final k-point convergence submission | 21084915 | COMPLETED | 0.2176 |
| E–a scan | 21084907 | COMPLETED | 0.7488 |
| **Total J001 resource consumption** | | | **1.3528 SU** |

The k-point convergence calculations were completed across two Slurm submissions. Job 21084888 was an accidentally duplicated submission and was cancelled after three VASP subjobs had completed. Job 21084915 completed the remaining k-point calculation work used for the final convergence result.

The `jobsu` values are the authoritative Anvil accounting values. Parser-reported core-hours represent VASP execution time and are retained separately for performance analysis.

For the final scientific J001 dataset, the duplicate cancelled submission is not treated as an additional convergence point; its resource consumption is reported separately for transparent provenance.

## 9. Scope: what J001 does not establish

- Convergence of the valley splitting ΔE itself (to be checked on the SOC runs).
- Magnetic-layer settings (Hubbard U, magnetic ground state), vdW and dipole corrections (J002, J003).
- Vacuum convergence (c = 23.6 Å was not varied).

## 10. Open checks

- [ ] k-series bookkeeping: parser reports 0.353 core-h for `conv_kpts`, but job 21084915 could use at most 0.218 core-h; look for an earlier job.
- [ ] Refit a0 using the 5 points nearest the minimum; quote 3.56 Å if it agrees within about 0.002 Å.
- [ ] Compare a0 with C2DB once the WTe2 row is verified.
- [ ] Run spglib on the relaxed E–a `CONTCAR` files.
- [ ] Optional vacuum test (J001d).


## 11. No-SOC WTe2 Band-Structure Validation

A non-spin-polarized, no-SOC electronic-structure calculation was
performed for the equilibrium WTe2 monolayer using:

- a0 = 3.5597 Å
- ENCUT = 500 eV
- SCF k-mesh = 12×12×1
- Band path = Γ–K–M–K′–Γ
- ISPIN = 1
- SOC = OFF

The SCF calculation (Slurm job 21193201) converged to the requested
EDIFF = 1E-6 and produced the converged charge density `CHGCAR`.
A subsequent non-self-consistent band calculation (Slurm job 21194532)
used `ICHARG = 11`.

### K/K′ Degeneracy Test

The explicit valley points were:

- K  = (1/3, 1/3, 0)
- K′ = (2/3, 2/3, 0)

The extracted energies were:

| Band | K (eV) | K′ (eV) | |ΔE| (eV) |
|---:|---:|---:|---:|
| 13 | -1.040795 | -1.040795 | 0.000000 |
| 14 | 0.007520 | 0.007520 | 0.000000 |

Therefore:

\[
\Delta E_{13}(K,K') =
|E_{13}(K)-E_{13}(K')|
= 0.000000\ \text{eV}
\]

\[
\Delta E_{14}(K,K') =
|E_{14}(K)-E_{14}(K')|
= 0.000000\ \text{eV}
\]

Within the precision of the extracted eigenvalues, the K and K′ states
are degenerate in the no-SOC, nonmagnetic WTe2 baseline.

This provides the required control calculation for the P4 valley-physics
workflow: no artificial K/K′ splitting is observed before introducing
spin–orbit coupling and magnetic proximity effects.

### Resource Usage

| Calculation | Slurm job | Wall time | Cores | Anvil SU |
|---|---:|---:|---:|---:|
| No-SOC SCF | 21193201 | 00:10:23 | 8 | 1.3848 |
| No-SOC BAND | 21194532 | 00:01:25 | 8 | 0.1888 |
| **Subtotal** | | | | **1.5736** |

The cumulative documented J001 Anvil usage is now **2.9264 SU**,
including the earlier convergence, k-point, E–a, SCF, and band
calculations. The cancelled duplicate k-point submission is included
because it consumed 0.1664 SU.

### Baseline Result

$$
\boxed{\Delta E_{K,K'} = 0\ \text{eV}}
$$

within the precision of the extracted no-SOC eigenvalues.