#!/usr/bin/env python3
"""
parse_j004.py - valley splitting from the J004 SOC runs (Mn moment along +z and along -z).

Usage:  python3 parse_j004.py j004_smoke
        python3 parse_j004.py j004 --show SeW_MnH/pz      (every state near E_F at K and K', with character)
        options: --window 2.0 (eV around E_F), --wl 0.5, --wv 0.3, --wc 0.3 (thresholds, see below)

Reads <folder>/<name>/<tag>/{OUTCAR,OSZICAR,POSCAR,PROCAR,KPOINTS} for every finished run.

K and K' come from the geometry, not from fixed coordinates: they are the two Brillouin-zone corners
on the k-mesh, and K is the corner 30 degrees counter-clockwise (seen from +z) from the in-plane
W->Te bond. For our cells (gamma = 120 deg, setup_j003.py) this is K = (1/3, 1/3), the K of J001.
For another cell setting (e.g. his POSCAR) the same physical valley is found automatically.

WTe2 band edges at K and K' are picked by orbital character, because the heterostructure is metallic
and counting occupied bands does not work. Weights are fractions of the PAW-sphere total (PROCAR):
  WTe2 state : W + Te weight >= wl
  v-type     : W d(xy) + W d(x2-y2) >= wv      (WTe2 valence-band top at K: d orbitals with m = +-2)
  c-type     : W d(z2) >= wc                   (WTe2 conduction-band bottom at K: m = 0)
  VB, CB     : the adjacent (v-type, c-type) pair with the largest gap inside E_F +- window.
Paper's definitions (Fig. S3):  dVB = E_VB(K) - E_VB(K'),  dCB = E_CB(K) - E_CB(K'),
                                dKK' = dVB - dCB  ( = gap(K') - gap(K) ).
Paper, 0 % strain, Fig. S3: dVB = 260, dCB = 22, dKK' = 236 meV. Which corner the paper calls K is
not known yet (it needs his band-path KPOINTS), so compare magnitudes first.

Checks per run: OUTCAR finished, SCF reached EDIFF, SOC on, Mn moment along the intended axis
(tilt below 5 deg), VB(K) and VB(K') carry opposite spin, gaps between 0.2 and 1.6 eV, and no
partly-WTe2 state (hybridised with MnSe2) close to the chosen edges.
Checks per +z/-z pair (time reversal): E(+z) = E(-z) and dKK'(-z) = -dKK'(+z).
Writes <folder>/j004_summary.csv.

v3 (27 Oct): a run whose OUTCAR has the timing block but no E-fermi and no 'EDIFF is reached' is reported
as 'stopped early' (with the SCF step it reached) instead of crashing --show. New PAPER FINGERPRINTS table,
one line per run, for comparing with SI Fig. S4:
  maj./min. valley : the valley whose VB top has Mn's majority / minority spin
  maj VB - min VB  : > 0 means the majority-spin VB is pushed up (our runs with U on Mn), < 0 down (the paper)
  min-spin SOC     : VB(min valley) - lower VB(maj valley), both minority spin; MnSe2 leaves them alone, so this
                     fingerprints U on W (0.49 eV without it, 0.62 eV with 2.87 eV; paper 0.485)
  VB(min) - E_F    : the paper has the WTe2 VB top at E_F (about -0.01 eV)
  partner, t, offset, push : two-level fit of the majority-spin VB with the partly-WTe2 state of the same spin
                     that mixes with it (WTe2 share 0.10-0.50, W d(xy,x2-y2) >= 0.08, within 0.8 eV); offset > 0
                     means the MnSe2 band lies below the VB and pushes it up
"""
import argparse, csv, math, os, re, sys
import numpy as np

FLOAT = re.compile(r"-?\d+\.\d+(?:[Ee][-+]?\d+)?")
PAPER = {"dVB": 260.0, "dCB": 22.0, "dKK": 236.0}
NAN = float("nan")

def element(name):
    return name[:1].upper() + (name[1:2] if name[1:2].islower() else "")

def read_poscar(path):
    L = open(path).read().splitlines()
    s = float(L[1].split()[0])
    cell = np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * s
    names, counts = L[5].split(), [int(x) for x in L[6].split()]
    i = 7
    if L[i].strip()[:1].lower() == "s":
        i += 1
    direct = L[i].strip()[:1].lower() == "d"
    pos = np.array([[float(x) for x in L[i + 1 + k].split()[:3]] for k in range(sum(counts))])
    frac = pos if direct else pos @ np.linalg.inv(cell)
    return cell, frac, [element(n) for n, c in zip(names, counts) for _ in range(c)]

def interlayer_d(cell, frac, sp):
    te_i = [i for i, s in enumerate(sp) if s == "Te"]
    se_i = [i for i, s in enumerate(sp) if s == "Se"]
    z = frac[:, 2] % 1.0
    zs = np.sort(z)
    gaps = np.diff(np.concatenate([zs, [zs[0] + 1.0]]))
    start = zs[(int(np.argmax(gaps)) + 1) % len(zs)]
    zu = ((z - start) % 1.0) * cell[2, 2]
    te, se = zu[te_i], zu[se_i]
    return float(se.min() - te.max()) if se.mean() > te.mean() else float(te.min() - se.max())

# ---------- which mesh points are K and K' ----------
def w_te_bond(cell, frac, sp):
    cart, iw, best = frac @ cell, sp.index("W"), None
    for it in [i for i, s in enumerate(sp) if s == "Te"]:
        for n1 in (-1, 0, 1):
            for n2 in (-1, 0, 1):
                v = (cart[it] + n1 * cell[0] + n2 * cell[1] - cart[iw])[:2]
                if best is None or np.linalg.norm(v) < np.linalg.norm(best) - 1e-9:
                    best = v
    return best

def make_classifier(cell, frac, sp):
    B = 2 * np.pi * np.linalg.inv(cell).T                        # rows: b1, b2, b3
    bond = w_te_bond(cell, frac, sp)
    kK = 4 * np.pi / (3 * np.linalg.norm(cell[0]))
    def classify(kred):
        if len(kred) < 3 or abs(kred[2]) > 1e-6:
            return None
        kc, best = np.asarray(kred) @ B, None
        for n1 in range(-2, 3):
            for n2 in range(-2, 3):
                v = (kc + n1 * B[0] + n2 * B[1])[:2]
                if best is None or np.linalg.norm(v) < np.linalg.norm(best) - 1e-12:
                    best = v
        if abs(np.linalg.norm(best) - kK) > 1e-4 * kK:
            return None
        dphi = math.degrees(math.atan2(best[1], best[0]) - math.atan2(bond[1], bond[0])) % 120.0
        return "K" if abs(dphi - 30) < 5 else ("K'" if abs(dphi - 90) < 5 else None)
    return classify

# ---------- VASP output readers ----------
def read_procar(path, classify):
    """Bands at K and K' only. Each band: E, occ, column names, rows = 4 x nion rows (tot, mx, my, mz)."""
    lines = open(path, errors="replace").read().splitlines()
    m = re.search(r"# of ions:\s*(\d+)", " ".join(lines[:3]))
    if not m:
        raise ValueError("no '# of ions' in the PROCAR header")
    nion, out, label, i, n = int(m.group(1)), {}, None, 0, len(lines)
    while i < n:
        s = lines[i].lstrip()
        if s.startswith("k-point"):
            part = s.split(":", 1)[1].split("weight")[0] if ":" in s else ""
            kred = [float(x) for x in FLOAT.findall(part)][:3]       # handles '0.66666667-0.33333333'
            label = classify(kred)
            if label in out:
                label = None                                        # keep the first copy of each valley
            if label:
                out[label] = {"k": kred, "bands": []}
            i += 1
            continue
        if label and s.startswith("band"):
            mE, mO = re.search(r"energy\s+(\S+)", s), re.search(r"occ\.\s+(\S+)", s)
            band = {"E": float(mE.group(1)), "occ": float(mO.group(1)) if mO else NAN, "cols": None, "rows": []}
            j = i + 1
            while j < n and len(band["rows"]) < 4 * nion:
                t = lines[j].strip()
                if not t or t.startswith("tot"):
                    j += 1
                    continue
                tok = t.split()
                if tok[0] == "ion":
                    band["cols"] = band["cols"] or tok[1:]
                elif tok[0].isdigit():
                    band["rows"].append([float(x) for x in FLOAT.findall(t[len(tok[0]):])])
                else:
                    break                                           # next band or k-point
                j += 1
            out[label]["bands"].append(band)
            i = j
            continue
        i += 1
    return out, nion

def read_outcar(path):
    if not os.path.isfile(path):
        return {"finished": False}
    t = open(path, errors="replace").read()
    r = {"finished": "General timing and accounting" in t,
         "nscf": len(re.findall(r"Iteration\s+\d+\(\s*\d+\)", t)), "has_ef": "E-fermi" in t}
    if not r["finished"]:
        return r
    def last(pat, cast=float, default=NAN):
        v = re.findall(pat, t)
        return cast(v[-1]) if v else default
    r.update(converged="EDIFF is reached" in t and "EDIFF was not reached" not in t,
             E=last(r"free  energy   TOTEN\s*=\s*(\S+)\s*eV"), Ef=last(r"E-fermi :\s*(\S+)"),
             soc=bool(re.search(r"LSORBIT\s*=\s*T", t)), nbands=last(r"NBANDS=\s*(\d+)", int, None),
             nelect=last(r"NELECT\s*=\s*(\S+)"), encut=last(r"ENCUT\s*=\s*(\S+)\s*eV"),
             wall=last(r"Elapsed time \(sec\):\s*(\S+)"))
    lines = t.splitlines()
    for ax in "xyz":
        tab = []
        for k in [i for i, l in enumerate(lines) if l.strip() == f"magnetization ({ax})"]:
            j = k + 1
            while j < len(lines) and not lines[j].strip().startswith("# of ion"):
                j += 1
            j += 2
            cur = []
            while j < len(lines) and lines[j].split() and lines[j].split()[0].isdigit():
                cur.append(float(lines[j].split()[-1]))
                j += 1
            if cur:
                tab = cur
        r["m" + ax] = tab
    return r

def read_oszicar_mag(path):
    if not os.path.isfile(path):
        return [NAN, NAN, NAN]
    mags = [l for l in open(path, errors="replace") if "mag=" in l]
    if not mags:
        return [NAN, NAN, NAN]
    v = [float(x) for x in re.findall(r"-?\d+\.\d*(?:[Ee][-+]?\d+)?", mags[-1].split("mag=")[1])]
    return (v + [NAN, NAN, NAN])[:3] if len(v) >= 3 else [0.0, 0.0, v[0]] if v else [NAN, NAN, NAN]

def read_kmesh(path):
    try:
        L = open(path).read().splitlines()
        return "x".join(L[3].split()[:3])
    except (OSError, IndexError):
        return "?"

# ---------- band-edge picking ----------
def col_index(cols):
    def find(*names):
        for nm in names:
            if nm in cols:
                return cols.index(nm)
        raise ValueError(f"PROCAR columns {cols} lack {names}")
    return find("dxy"), find("dz2"), find("x2-y2", "dx2-y2", "dx2"), find("tot")

def characterise(bands, sp, nion):
    iW = [i for i, s in enumerate(sp) if s == "W"]
    iL = [i for i, s in enumerate(sp) if s in ("W", "Te")]
    out = []
    for b in bands:
        rows = np.array(b["rows"][:nion])
        ixy, iz2, ix2, itot = col_index(b["cols"])
        P = rows[:, itot].sum()
        mz = np.array(b["rows"][3 * nion:4 * nion]) if len(b["rows"]) >= 4 * nion else None
        if P <= 1e-6:
            out.append({"E": b["E"], "occ": b["occ"], "wl": 0.0, "cv": 0.0, "cc": 0.0, "sz": NAN, "typ": None})
            continue
        st = {"E": b["E"], "occ": b["occ"], "wl": rows[iL, itot].sum() / P,
              "cv": (rows[iW, ixy] + rows[iW, ix2]).sum() / P, "cc": rows[iW, iz2].sum() / P,
              "sz": (mz[:, itot].sum() / P) if mz is not None else NAN}
        out.append(st)
    return out

def classify_states(states, a):
    for s in states:
        s["typ"] = None
        if s["wl"] >= a.wl:
            if s["cv"] >= a.wv and s["cv"] > s["cc"]:
                s["typ"] = "v"
            elif s["cc"] >= a.wc and s["cc"] > s["cv"]:
                s["typ"] = "c"

def pick_edges(states, Ef, a):
    typed = sorted([s for s in states if s["typ"] and abs(s["E"] - Ef) <= a.window], key=lambda s: s["E"])
    best = None
    for s1, s2 in zip(typed, typed[1:]):
        if s1["typ"] == "v" and s2["typ"] == "c" and (best is None or s2["E"] - s1["E"] > best[1]["E"] - best[0]["E"]):
            best = (s1, s2)
    return best

def near_misses(states, vb, cb, a):
    """WTe2-like character but below the layer threshold, close to the chosen edges: a hybridisation warning."""
    lo, hi = vb["E"] - 0.10, cb["E"] + 0.10
    return [s for s in states if s["typ"] is None and lo <= s["E"] <= hi and s["wl"] >= 0.25
            and (s["cv"] >= 0.6 * a.wv or s["cc"] >= 0.6 * a.wc)]

# ---------- one run ----------
def analyse(d, a):
    r = {"run": os.path.relpath(d, a.folder), "warn": []}
    o = read_outcar(os.path.join(d, "OUTCAR"))
    if not o["finished"]:
        r["status"] = "not finished"
        return r
    if not o["converged"] and not o["has_ef"]:
        r["status"] = (f"stopped early at SCF step {o['nscf']} (no E-fermi, EDIFF not reached): "
                       "rename its OUTCAR and resubmit")
        return r
    r.update(o)
    r["status"] = "ok"
    if not o["soc"]:
        r["warn"].append("LSORBIT is not on in this run")
    if not o["converged"]:
        r["warn"].append("SCF did not reach EDIFF")
    cell, frac, sp = read_poscar(os.path.join(d, "POSCAR"))
    r["sp"], r["d"], r["kmesh"] = sp, interlayer_d(cell, frac, sp), read_kmesh(os.path.join(d, "KPOINTS"))
    r["mag_cell"] = read_oszicar_mag(os.path.join(d, "OSZICAR"))
    tag = os.path.basename(d)
    r["want"] = +1 if tag.startswith("pz") else (-1 if tag.startswith("mz") else 0)
    iMn = [i for i, s in enumerate(sp) if s == "Mn"]
    iW = [i for i, s in enumerate(sp) if s == "W"]
    if len(o["mz"]) == len(sp):
        mx = o["mx"] if len(o["mx"]) == len(sp) else [0.0] * len(sp)
        my = o["my"] if len(o["my"]) == len(sp) else [0.0] * len(sp)
        m = np.array([sum(mx[i] for i in iMn), sum(my[i] for i in iMn), sum(o["mz"][i] for i in iMn)])
        r["m_Mn"], r["m_W_z"] = m, sum(o["mz"][i] for i in iW)
        r["theta"] = math.degrees(math.acos(min(1.0, abs(m[2]) / max(np.linalg.norm(m), 1e-9))))
        if np.linalg.norm(m) < 2.0:
            r["warn"].append(f"Mn moment only {np.linalg.norm(m):.2f} muB")
        if r["theta"] > 5:
            r["warn"].append(f"Mn moment tilted {r['theta']:.1f} deg from z")
        if r["want"] and np.sign(m[2]) != r["want"]:
            r["warn"].append("Mn moment points the wrong way for this tag")
    else:
        r["m_Mn"], r["m_W_z"], r["theta"] = np.array([NAN] * 3), NAN, NAN
        r["warn"].append("no magnetization (x/y/z) tables in OUTCAR (LORBIT missing?)")
    pro = os.path.join(d, "PROCAR")
    if not os.path.isfile(pro):
        r["warn"].append("no PROCAR")
        return r
    kp, nion = read_procar(pro, make_classifier(cell, frac, sp))
    if set(kp) != {"K", "K'"}:
        r["warn"].append(f"K/K' not both on the k-mesh (found {sorted(kp)}); use a multiple of 3")
        return r
    r["states"] = {}
    for v in ("K", "K'"):
        st = characterise(kp[v]["bands"], sp, nion)
        classify_states(st, a)
        r["states"][v] = st
        r["k_" + v] = kp[v]["k"]
        e = pick_edges(st, o["Ef"], a)
        if e is None:
            r["warn"].append(f"no WTe2 (v, c) edge pair at {v} within E_F +- {a.window} eV")
            return r
        r["VB_" + v], r["CB_" + v] = e
        miss = near_misses(st, e[0], e[1], a)
        if miss:
            r["warn"].append(f"{v}: {len(miss)} partly-WTe2 state(s) near the edges (hybridised; see --show)")
    vk, vkp, ck, ckp = r["VB_K"], r["VB_K'"], r["CB_K"], r["CB_K'"]
    r["dVB"] = (vk["E"] - vkp["E"]) * 1000
    r["dCB"] = (ck["E"] - ckp["E"]) * 1000
    r["dKK"] = r["dVB"] - r["dCB"]
    r["gapK"], r["gapKp"] = ck["E"] - vk["E"], ckp["E"] - vkp["E"]
    for g, v in ((r["gapK"], "K"), (r["gapKp"], "K'")):
        if not 0.2 <= g <= 1.6:
            r["warn"].append(f"gap at {v} = {g:.3f} eV looks wrong for WTe2")
    if not (np.isnan(vk["sz"]) or np.isnan(vkp["sz"])) and vk["sz"] * vkp["sz"] > 0:
        r["warn"].append("VB(K) and VB(K') have the same spin: an edge may be mis-picked")
    r["fp"] = fingerprints(r)
    return r

def fingerprints(r):
    """Quantities that SI Fig. S4 shows directly: which way the majority-spin VB moves, the minority-spin SOC
    splitting (a U-on-W fingerprint), where E_F sits, and a two-level fit of the majority-spin VB with its
    partly-WTe2 partner. Returns {} when the spins are not clear enough."""
    mz = r["m_Mn"][2] if len(r.get("m_Mn", [])) == 3 else NAN
    vb = {v: r["VB_" + v] for v in ("K", "K'")}
    if np.isnan(mz) or any(np.isnan(vb[v]["sz"]) for v in vb):
        return {}
    maj = [v for v in vb if np.sign(vb[v]["sz"]) == np.sign(mz)]
    if len(maj) != 1:
        return {}
    vmaj = maj[0]
    vmin = "K'" if vmaj == "K" else "K"
    fp = {"maj": vmaj, "shift": (vb[vmaj]["E"] - vb[vmin]["E"]) * 1000, "vbmin_ef": vb[vmin]["E"] - r["Ef"]}
    fp["dir"] = "up" if fp["shift"] > 0 else "down"
    st = r["states"][vmaj]
    top = vb[vmaj]
    below = [s for s in st if s["typ"] == "v" and s["E"] < top["E"] - 1e-6 and top["E"] - s["E"] < 1.0
             and not np.isnan(s["sz"]) and np.sign(s["sz"]) != np.sign(top["sz"])]
    fp["soc_min"] = vb[vmin]["E"] - max(below, key=lambda s: s["E"])["E"] if below else NAN
    cands = [s for s in st if s is not top and 0.10 <= s["wl"] < 0.50 and s["cv"] >= 0.08 and abs(s["E"] - top["E"]) <= 0.8
             and not np.isnan(s["sz"]) and np.sign(s["sz"]) == np.sign(top["sz"])]
    if cands:
        p = max(cands, key=lambda s: s["wl"])
        F = top["wl"] + p["wl"]
        x1, x2 = top["wl"] / F, p["wl"] / F
        bare_w = x1 * top["E"] + x2 * p["E"]
        bare_m = x2 * top["E"] + x1 * p["E"]
        fp.update(partner_E=p["E"] - r["Ef"], partner_wl=p["wl"], share_sum=F, offset=bare_w - bare_m,
                  t=abs(top["E"] - p["E"]) * math.sqrt(x1 * x2), push=(top["E"] - bare_w) * 1000)
        fp["rest"] = fp["shift"] - fp["push"]          # what the shift would be without the mixing
    return fp

def show(r):
    for v in ("K", "K'"):
        if v not in r.get("states", {}):
            print(f"\n{r['run']}  {v}: not analysed (the run stopped before this valley could be read)")
            continue
        k = r.get("k_" + v)
        print(f"\n{r['run']}  {v} = ({k[0]:.4f}, {k[1]:.4f})   E_F = {r['Ef']:.4f} eV")
        print(f"   {'E-E_F':>7s} {'occ':>5s} {'W+Te':>5s} {'Wd(xy,x2y2)':>11s} {'Wd(z2)':>7s} {'s_z':>6s} type")
        for s in r["states"][v]:
            if abs(s["E"] - r["Ef"]) > 2.5:
                continue
            mark = "  <- VB" if s is r.get("VB_" + v) else ("  <- CB" if s is r.get("CB_" + v) else "")
            print(f"   {s['E'] - r['Ef']:7.3f} {s['occ']:5.2f} {s['wl']:5.2f} {s['cv']:11.2f} {s['cc']:7.2f} "
                  f"{s['sz']:6.2f} {s['typ'] or '-':>4s}{mark}")

def fp_cells(fp):
    f = lambda k, fmt: format(fp[k], fmt) if k in fp else ""
    return [fp.get("maj", ""), f("shift", ".2f"), f("soc_min", ".4f"), f("vbmin_ef", ".4f"), f("partner_E", ".4f"),
            f("partner_wl", ".3f"), f("t", ".4f"), f("offset", ".4f"), f("push", ".1f")]

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("folder")
    ap.add_argument("--show", default="", help="run to print state by state, e.g. SeW_MnH/pz")
    ap.add_argument("--window", type=float, default=2.0, help="eV around E_F searched for the edges")
    ap.add_argument("--wl", type=float, default=0.5, help="min. W+Te weight for a WTe2 state")
    ap.add_argument("--wv", type=float, default=0.3, help="min. W d(xy)+d(x2-y2) weight for v-type")
    ap.add_argument("--wc", type=float, default=0.3, help="min. W d(z2) weight for c-type")
    a = ap.parse_args()
    runs = sorted(os.path.dirname(os.path.join(dp, f)) for dp, _, fs in os.walk(a.folder) for f in fs if f == "INCAR")
    if not runs:
        sys.exit(f"no run folders under {a.folder}")
    res = [analyse(d, a) for d in runs]
    print(f"J004 valley splitting, SOC. Energies in meV. Paper (0 % strain, Fig. S3): "
          f"dVB {PAPER['dVB']:.0f}, dCB {PAPER['dCB']:.0f}, dKK' {PAPER['dKK']:.0f}")
    hdr = ["run".ljust(22), "k".rjust(7), "E (eV)".rjust(13), "Mn moment (x, y, z)".rjust(21), "tilt".rjust(5),
           "m_W z".rjust(6), "d (A)".rjust(6), "gap K".rjust(6), "gap K'".rjust(6), "dVB".rjust(7),
           "dCB".rjust(7), "dKK'".rjust(7)]
    print(" ".join(hdr))
    for r in res:
        if r["status"] != "ok":
            print(f"{r['run']:22s}   {r['status']}")
            continue
        m = r["m_Mn"]
        line = (f"{r['run']:22s} {r['kmesh']:>7s} {r['E']:13.6f}  ({m[0]:5.2f},{m[1]:5.2f},{m[2]:+6.3f}) "
                f"{r['theta']:5.1f} {r['m_W_z']:+6.3f} {r['d']:6.3f}")
        if "dKK" in r:
            line += f" {r['gapK']:6.3f} {r['gapKp']:6.3f} {r['dVB']:+7.1f} {r['dCB']:+7.1f} {r['dKK']:+7.1f}"
        print(line)
        for w in r["warn"]:
            print(f"{'':24s}! {w}")
    ok = {r["run"]: r for r in res if r["status"] == "ok"}
    printed = False
    for run, r in sorted(ok.items()):
        base, tag = os.path.split(run)
        if not tag.startswith("pz"):
            continue
        partner = ok.get(os.path.join(base, "mz" + tag[2:]))
        if not partner:
            continue
        if not printed:
            print("\nTIME-REVERSAL PAIRS (+z vs -z)")
            printed = True
        dE = (r["E"] - partner["E"]) * 1e6
        line = f"   {run:20s} E(+z) - E(-z) = {dE:+8.1f} ueV"
        verdict = abs(dE) < 100
        if "dKK" in r and "dKK" in partner:
            s = r["dKK"] + partner["dKK"]
            line += f"   dKK'(+z) + dKK'(-z) = {s:+6.2f} meV"
            verdict = verdict and abs(s) < 2.0
        print(line + ("   ok" if verdict else "   CHECK: should both be about 0"))
    printed = False
    for run, r in sorted(ok.items()):
        base, tag = os.path.split(run)
        if tag.endswith("_clean"):
            continue
        other = ok.get(os.path.join(base, tag + "_clean"))
        if not other:
            continue
        if not printed:
            print("\nSTARTING MOMENTS (his: Se, Te at -2 muB; clean: 0)")
            printed = True
        dE = (r["E"] - other["E"]) * 1000
        line = f"   {run:20s} E(his start) - E(clean start) = {dE:+8.3f} meV"
        same = abs(dE) < 0.1
        if "dKK" in r and "dKK" in other:
            line += f"   dKK' differs by {r['dKK'] - other['dKK']:+6.2f} meV"
            same = same and abs(r["dKK"] - other["dKK"]) < 1.0
        print(line + ("   same state" if same else "   DIFFERENT states: send me both OUTCARs"))
    fps = [(run, r["fp"]) for run, r in sorted(ok.items()) if r.get("fp")]
    if fps:
        print("\nPAPER FINGERPRINTS (SI Fig. S4: majority VB pushed down by 256-331 meV, min-spin SOC 0.485 eV, "
              "VB(min) at E_F = -0.01 eV, partner above)")
        print(f"   {'run':22s} maj  maj-min VB    min-spin SOC  VB(min)-E_F   partner (share)   t (eV)  "
              f"offset (eV)  push  rest (meV)")
        for run, fp in fps:
            line = (f"   {run:22s} {fp['maj']:3s} {fp['shift']:+7.1f} {fp['dir']:4s}"
                    f"  {fp['soc_min']:9.3f}   {fp['vbmin_ef']:+10.3f}")
            if "t" in fp:
                line += (f"   {fp['partner_E']:+7.3f} ({fp['partner_wl']:.2f})   {fp['t']:6.3f}   {fp['offset']:+8.3f}"
                         f"  {fp['push']:+5.0f} {fp['rest']:+6.0f}")
            else:
                line += "   no partly-WTe2 partner found"
            print(line)
        print("   offset > 0: the MnSe2 partner lies below the majority-spin VB and pushes it up; < 0: above, pushes it down")
    if a.show:
        r = ok.get(a.show.rstrip("/"))
        if not r or "states" not in r:
            print(f"\n--show {a.show}: no finished run with K/K' states under that name")
        else:
            show(r)
    prod = [r for r in ok.values() if os.path.basename(r["run"]) == "pz" and "dKK" in r]
    if prod:
        print("\nFOR valley_results (Group D), from the +z runs (his start):")
        for r in sorted(prod, key=lambda r: r["run"]):
            print(f"   {os.path.dirname(r['run'])} (k {r['kmesh']}): dE_VB_meV = {r['dVB']:+.1f}, "
                  f"dE_CB_meV = {r['dCB']:+.1f}, dE_KKp_meV = {r['dKK']:+.1f}, m_W_z = {r['m_W_z']:+.3f}"
                  + ("   <- has warnings: do not enter yet" if r["warn"] else ""))
    path = os.path.join(a.folder, "j004_summary.csv")
    cols = ["run", "kmesh", "E_eV", "E_F_eV", "m_Mn_x", "m_Mn_y", "m_Mn_z", "tilt_deg", "m_W_z", "mag_cell_z",
            "d_A", "E_VB_K", "E_VB_Kp", "E_CB_K", "E_CB_Kp", "gap_K", "gap_Kp", "dVB_meV", "dCB_meV",
            "dKKp_meV", "VB_K_minus_EF", "nbands", "nelect", "encut", "wall_s", "maj_valley", "majVB_minus_minVB_meV",
            "soc_minority_eV", "VBmin_minus_EF", "partner_minus_EF", "partner_WTe2_share", "t_eV", "offset_eV", "push_meV",
            "warnings"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in res:
            if r["status"] != "ok":
                continue
            e = lambda key: r[key]["E"] if key in r else NAN
            w.writerow([r["run"], r["kmesh"], f"{r['E']:.8f}", f"{r['Ef']:.4f}",
                        *[f"{x:.4f}" for x in r["m_Mn"]], f"{r['theta']:.2f}", f"{r['m_W_z']:.4f}",
                        f"{r['mag_cell'][2]:.4f}", f"{r['d']:.4f}",
                        *[f"{e(k):.4f}" for k in ("VB_K", "VB_K'", "CB_K", "CB_K'")],
                        f"{r.get('gapK', NAN):.4f}", f"{r.get('gapKp', NAN):.4f}", f"{r.get('dVB', NAN):.2f}",
                        f"{r.get('dCB', NAN):.2f}", f"{r.get('dKK', NAN):.2f}", f"{e('VB_K') - r['Ef']:.4f}",
                        r["nbands"], r["nelect"], r["encut"], f"{r['wall']:.0f}", *fp_cells(r.get("fp", {})),
                        "; ".join(r["warn"])])
    print(f"\nwritten: {path}")

if __name__ == "__main__":
    main()
