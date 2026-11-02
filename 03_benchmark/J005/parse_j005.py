#!/usr/bin/env python3
"""
parse_j005.py - J005: work functions and band edges of the free layers on a common vacuum scale.

Usage:  python3 parse_j005.py J005 [--no-write-vz]   (run from 03_benchmark; parses both j005 and j005_refined)

Per run: the vacuum level V_vac (planar average of LOCPOT along z, averaged over the middle of the vacuum gap; the
'flat' column is its spread there and should be a few meV), E_F, the work function WF = V_vac - E_F, and the states
at the BZ corner K on the vacuum scale:
  WTe2:  VB(K) and CB(K) (and with SOC the lower VB branch: the SOC splitting at K), IP_K = V_vac - VB(K)
  MnSe2: the lowest two empty majority-spin states at K (the doublet the K' partner band in J004 comes from)
The free-layer offset is WF(MnSe2) - IP_K(WTe2). The paper's approximate SI Fig. S2 values imply only a bound if 4.6 eV is a WTe2 work function: IP_K >= 4.6 eV and offset <= 0.8 eV. The previous +0.66 eV point estimate mixed quantities whose SOC provenance is not established and is withdrawn. The planar-averaged potential is written to <run>/vz.csv unless --no-write-vz is used.
"""
import argparse
import math, os, sys
import numpy as np

# SI Fig. S2 reports approximate work functions. The SOC treatment of the
# WTe2 potential/work-function calculation is not independently confirmed.
# Do not combine WF_W=4.6 eV with the -0.139 eV SOC band position from SI Fig. S1b:
# that mixed-provenance construction produced the withdrawn +0.66 eV estimate.
PAPER = {"WF_W_REPORTED": 4.6, "WF_M_REPORTED": 5.4}


def read_poscar(path):
    L = open(path).read().splitlines()
    s = float(L[1].split()[0])
    A = np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * s
    names, counts = L[5].split(), [int(x) for x in L[6].split()]
    i = 7
    if L[i].strip()[:1].lower() == "s":
        i += 1
    direct = L[i].strip()[:1].lower() == "d"
    pos = np.array([[float(x) for x in L[i + 1 + k].split()[:3]] for k in range(sum(counts))])
    frac = pos if direct else pos @ np.linalg.inv(A)
    sp = [n[:1].upper() + (n[1:2] if n[1:2].islower() else "") for n, c in zip(names, counts) for _ in range(c)]
    return A, frac, sp


def read_outcar(d):
    p = os.path.join(d, "OUTCAR")
    if not os.path.isfile(p):
        return None
    ef, soc, fin, conv = None, False, False, False
    with open(p, errors="replace") as f:
        for line in f:
            if "E-fermi" in line:
                try:
                    ef = float(line.split(":")[1].split()[0])
                except (IndexError, ValueError):
                    pass
            elif "LSORBIT" in line and "=" in line:
                soc = line.split("=")[1].split()[0].upper().startswith("T")
            elif "aborting loop because EDIFF is reached" in line:
                conv = True
            elif "General timing and accounting" in line:
                fin = True
    return {"ef": ef, "soc": soc, "finished": fin, "converged": conv}


def read_incar_u(d, species):
    """U per element from LDAUL / LDAUU (one entry per species, POSCAR order)."""
    tags = {}
    for line in open(os.path.join(d, "INCAR"), errors="replace"):
        line = line.split("#")[0].split("!")[0]
        if "=" in line:
            k, v = line.split("=", 1)
            tags[k.strip().upper()] = v.split()
    uniq = list(dict.fromkeys(species))
    if tags.get("LDAU", [".FALSE."])[0].upper().strip(".").startswith("F"):
        return {}
    L = [int(x) for x in tags.get("LDAUL", [])]
    U = [float(x) for x in tags.get("LDAUU", [])]
    return {s: (U[i] if i < len(U) and i < len(L) and L[i] >= 0 else 0.0) for i, s in enumerate(uniq)}


def read_moment(d):
    p = os.path.join(d, "OSZICAR")
    if not os.path.isfile(p):
        return None
    mag = None
    for line in open(p, errors="replace"):
        if "mag=" in line:
            try:
                vals = [float(x) for x in line.split("mag=")[1].split()]
                mag = vals[-1] if vals else None
            except ValueError:
                pass
    return mag


def read_locpot(path):
    """Grid and the first data block (electrostatic potential with LVHAR), as an (nz, ny, nx) array."""
    with open(path, errors="replace") as f:
        L = f.read().splitlines()
    counts = [int(x) for x in L[6].split()]
    i = 7 + (1 if L[7].strip()[:1].lower() == "s" else 0) + 1 + sum(counts)
    while not L[i].strip():
        i += 1
    nx, ny, nz = [int(x) for x in L[i].split()[:3]]
    n, vals = nx * ny * nz, []
    for line in L[i + 1:]:
        vals.extend(line.split())
        if len(vals) >= n:
            break
    if len(vals) < n:
        raise ValueError(f"LOCPOT has {len(vals)} values, expected {n}")
    return np.array(vals[:n], dtype=float).reshape(nz, ny, nx)


def vacuum_level(A, frac, vz):
    """Mean of the planar average over the middle of the largest gap between atoms along c (cyclic)."""
    c = np.linalg.norm(A[2])
    z = np.sort(frac[:, 2] % 1.0)
    gaps = np.diff(np.concatenate([z, [z[0] + 1.0]]))
    j = int(np.argmax(gaps))
    zc = (z[j] + gaps[j] / 2) % 1.0
    half = min(2.0, gaps[j] * c / 4) / c
    zg = np.arange(len(vz)) / len(vz)
    dist = np.abs(((zg - zc) + 0.5) % 1.0 - 0.5)
    sel = vz[dist <= half]
    return float(sel.mean()), float(sel.max() - sel.min()), gaps[j] * c


def read_eigenval(path):
    L = open(path).read().splitlines()
    ispin = int(L[0].split()[3])
    nk, nb = [int(float(x)) for x in L[5].split()[1:3]]
    i, kpts, E, O = 6, [], [[] for _ in range(ispin)], [[] for _ in range(ispin)]
    for _ in range(nk):
        while not L[i].strip():
            i += 1
        kpts.append([float(x) for x in L[i].split()[:3]])
        i += 1
        rows = [L[i + b].split() for b in range(nb)]
        i += nb
        has_occ = len(rows[0]) >= 1 + 2 * ispin
        for s in range(ispin):
            E[s].append([float(r[1 + s]) for r in rows])
            O[s].append([float(r[1 + ispin + s]) for r in rows] if has_occ else [math.nan] * nb)
    return np.array(kpts), np.array(E), np.array(O), ispin


def corner_index(A, kpts):
    la = np.linalg.norm(A[0])
    B = 2 * math.pi * np.linalg.inv(A).T
    kK = 4 * math.pi / (3 * la)
    for j, k in enumerate(kpts):
        best = min(np.linalg.norm((k + np.array([n1, n2, 0])) @ B) for n1 in (-1, 0, 1) for n2 in (-1, 0, 1))
        if abs(best - kK) < 2e-3 * kK and abs(k[2]) < 1e-6:
            return j
    return None


def filled(e, o, ef):
    """Occupation flags; without printed occupations, split at E_F."""
    return [(y >= 0.5) if not math.isnan(y) else (x <= ef) for x, y in zip(e, o)]


def analyse(d, write_vz=True):
    r = {"run": os.path.basename(d.rstrip("/")), "dir": d}
    oc = read_outcar(d)
    if oc is None or oc["ef"] is None:
        r["status"] = "no OUTCAR / E-fermi (not run yet?)"
        return r
    if not (oc["finished"] and oc["converged"]):
        r["status"] = "not finished or EDIFF not reached: rerun it (sbatch run_j005.slurm on this folder)"
        return r
    pos = os.path.join(d, "POSCAR")
    A, frac, sp = read_poscar(pos)
    r.update(ef=oc["ef"], soc=oc["soc"], a=float(np.linalg.norm(A[0])), U=read_incar_u(d, sp),
             layer="WTe2" if "W" in sp else ("MnSe2" if "Mn" in sp else "?"), mag=read_moment(d))
    lp = os.path.join(d, "LOCPOT")
    if not (os.path.isfile(lp) and os.path.getsize(lp) > 0):
        r["status"] = "no LOCPOT (LVHAR missing?)"
        return r
    vz = read_locpot(lp).mean(axis=(1, 2))
    r["vac"], r["flat"], r["gap"] = vacuum_level(A, frac, vz)
    c = np.linalg.norm(A[2])
    if write_vz:
        with open(os.path.join(d, "vz.csv"), "w") as f:
            f.write("z_A,V_eV\n" + "".join(f"{k * c / len(vz):.4f},{v:.5f}\n" for k, v in enumerate(vz)))
    r["wf"] = r["vac"] - r["ef"]
    ev = os.path.join(d, "EIGENVAL")
    if os.path.isfile(ev) and os.path.getsize(ev) > 0:
        kpts, E, O, ispin = read_eigenval(ev)
        j = corner_index(A, kpts)
        if j is not None:
            if ispin == 2:
                s = 0 if (r["mag"] is None or r["mag"] >= 0) else 1      # majority spin channel
                e, o = E[s][j], O[s][j]
            else:
                e, o = E[0][j], O[0][j]
            fl = filled(e, o, r["ef"])
            occ = sorted(x for x, f_ in zip(e, fl) if f_)
            emp = sorted(x for x, f_ in zip(e, fl) if not f_)
            if r["layer"] == "WTe2" and occ and emp:
                r["vb"], r["cb"] = occ[-1], emp[0]
                if r["soc"] and len(occ) > 1:
                    r["vb2"] = occ[-2]
            if r["layer"] == "MnSe2" and len(emp) > 1:
                r["dbl"] = emp[:2]
    r["status"] = "ok"
    return r


def main():
    ap = argparse.ArgumentParser(
        description="Parse J005 free-layer work functions, K-point edges and vacuum-level alignment."
    )
    ap.add_argument("folder", help="Root folder containing the j005 and/or j005_refined run directories")
    ap.add_argument(
        "--no-write-vz", action="store_true",
        help="Read LOCPOT and report results without rewriting per-run vz.csv files",
    )
    a = ap.parse_args()
    root = os.path.abspath(a.folder)
    runs = sorted({dp for dp, _, fs in os.walk(root) if "prev_" not in dp
                   for f in fs if f == "INCAR"})
    if not runs:
        ap.error(f"no run folders with INCAR found under {a.folder}")

    res = [analyse(d, write_vz=not a.no_write_vz) for d in runs]
    # Use paths relative to the supplied root so same-named refined runs do not
    # overwrite each other in the result dictionary or appear as duplicate labels.
    for r in res:
        rel = os.path.relpath(r["dir"], root).replace(os.sep, "/")
        if rel == "j005_refined/MnSe2_a0":
            r["run"] = "j005_refined/MnSe2_a0_SOC_U3p99_variant"
        elif rel == "j005_refined/WTe2_UW_soc_a0":
            r["run"] = "j005_refined/WTe2_UW_soc_a0_repeat"
        else:
            r["run"] = rel

    print("J005 free layers on a common vacuum scale (energies in eV; vacuum = planar average of LOCPOT)")
    print(f"   {'run':40s} {'layer':5s} {'SOC':3s} {'U (eV)':10s} {'a (A)':>6s} {'V_vac':>8s} {'flat(meV)':>9s} "
          f"{'E_F':>8s} {'WF':>6s}  {'VB(K)-vac':>9s} {'CB(K)-vac':>9s} {'VB(K)-E_F':>9s} {'SOC at K':>8s}")
    ok = {}
    for r in res:
        if r["status"] != "ok":
            print(f"   {r['run']:40s} {r['status']}")
            continue
        ok[r["run"]] = r
        u = " ".join(f"{k} {v:.2f}" for k, v in r["U"].items() if v > 0) or "none"
        line = (f"   {r['run']:40s} {r['layer']:5s} {'yes' if r['soc'] else 'no':3s} {u:10s} {r['a']:6.4f} {r['vac']:8.4f} "
                f"{r['flat'] * 1000:9.1f} {r['ef']:8.4f} {r['wf']:6.3f}")
        if "vb" in r:
            soc_k = f"{r['vb'] - r['vb2']:8.3f}" if "vb2" in r else f"{'—':>8s}"
            line += f"  {r['vb'] - r['vac']:+9.3f} {r['cb'] - r['vac']:+9.3f} {r['vb'] - r['ef']:+9.3f} {soc_k}"
        if "dbl" in r:
            line += (f"  majority doublet at K: {r['dbl'][0] - r['ef']:+.3f} / {r['dbl'][1] - r['ef']:+.3f} above E_F, "
                     f"{r['dbl'][0] - r['vac']:+.3f} / {r['dbl'][1] - r['vac']:+.3f} vs vacuum")
        if r["flat"] > 0.02:
            line += "   ! vacuum not flat (dipole, or too little vacuum)"
        print(line)

    wf_w = PAPER["WF_W_REPORTED"]
    wf_m = PAPER["WF_M_REPORTED"]
    print("\nPUBLISHED ALIGNMENT BOUND (SI Fig. S2; interpretation conditional on the 4.6 eV value being WF)")
    print(f"   reported approximate values: WF(WTe2) {wf_w:.1f} eV; WF(MnSe2) {wf_m:.1f} eV")
    print("   For a semiconductor, WF(WTe2) <= IP_K(WTe2) = V_vac - E_VB(K).")
    print(f"   Therefore IP_K(WTe2) >= {wf_w:.1f} eV and WF(MnSe2) - IP_K(WTe2) <= {wf_m - wf_w:.1f} eV.")
    print("   This is a bound, not a point estimate. The source's SOC treatment / quantity definition is not confirmed.")
    print("   Do not combine the 4.6 eV SI value with the -0.139 eV SOC band position to infer +0.66 eV; that mixed-source estimate is withdrawn.")

    pairs = [
        ("j005/MnSe2", "j005/WTe2_UW0", "common a=3.5597 A; no SOC; U(W)=0"),
        ("j005/MnSe2", "j005/WTe2_UW", "common a=3.5597 A; no SOC; U(W)=2.87"),
        ("j005/MnSe2_soc", "j005/WTe2_UW0_soc", "common a=3.5597 A; both SOC; U(W)=0"),
        ("j005/MnSe2_soc", "j005/WTe2_UW_soc", "common a=3.5597 A; both SOC; U(W)=2.87"),
        ("j005_refined/MnSe2_a0_SOC_U3p99_variant", "j005_refined/WTe2_UW_soc_a0_repeat",
         "near each layer a0; both SOC; U(Mn)=3.99 and U(W)=2.87"),
    ]
    print("\nCOMPUTED FREE-LAYER ALIGNMENT PAIRS")
    for m, w, lab in pairs:
        if m in ok and w in ok and "vb" in ok[w]:
            M, W = ok[m], ok[w]
            ip = W["vac"] - W["vb"]
            line = (f"   {lab}: WF(MnSe2) {M['wf']:.3f}, IP_K(WTe2) {ip:.3f} "
                    f"-> WTe2 VB(K) above MnSe2 E_F {M['wf'] - ip:+.3f} eV")
            if "dbl" in M:
                line += f"; MnSe2 majority doublet minus WTe2 VB(K) {(M['dbl'][0] - M['vac']) - (W['vb'] - W['vac']):+.3f} eV"
            print(line)
    print("   Larger positive offset corresponds to a higher WTe2 VB(K) relative to the MnSe2 Fermi level in this free-layer lineup.")
    print("   The free-layer pair offset is not itself a direct measurement of the heterostructure interface dipole or charge transfer.")
    if a.no_write_vz:
        print("\nPer-run vz.csv files were left unchanged (--no-write-vz).")


if __name__ == "__main__":
    main()
