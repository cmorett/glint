"""Round 3: headline at 635 nm (measured dR), dR sensitivity at 780/850 nm, random chip
orientation at the organoid centre, flatness geometry, sensitivity ranking, pitch assets.

    python run_v3.py            # ~6 min on 16 cores; writes results/v3/ (rounds 1-2 untouched)
    python run_v3.py --quick    # 10x fewer photons (smoke test)

The decode rule (read time from SNR) is the one calibrated in round 2 on the decoder;
it is loaded from results/v2/decoder_floor.csv and decoder_ber.csv.
"""
import math
import os
import pickle
import sys
import time
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

import config_v3 as C
import run_v2 as R2
import stages
import stages_v2 as S2

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "v3")
QUICK = "--quick" in sys.argv
F = 10 if QUICK else 1
J = R2.J


def sd(*k):
    return np.random.SeedSequence([20261117, 3, *k])


def slab_job(args):
    li, mi, si = args
    lam, mu, seed = C.WAVELENGTHS_V3[li], C.MUSP_780_PER_MM[mi], C.SEEDS_V2[si]
    t0 = time.time()
    rng = np.random.default_rng(sd(10, li, mi, seed))
    s1 = {}
    for gi, (place, P) in enumerate(C.PLACEMENTS.items()):
        slab = stages.tissue_slab(mu, lam, P["n_cam"], P["n_far"])
        depths = [d / 1000 for d in (C.DEPTHS_V2_UM if place == "C_embedded" else [0])]
        s1[place] = S2.stage1_slab(slab, depths, C.N1_SLAB // F, sd(11, li, mi, seed, gi))
    slabC = stages.tissue_slab(mu, lam, C.N_MEDIUM, C.N_DISH)
    rows = []
    cams = list(C.CAMERAS)
    for place in ("A_top", "B_under"):
        side = "above" if place == "A_top" else "below"
        base = dict(geometry="slab", placement=place, D_mm=C.SLAB_THICKNESS_MM, depth_um=0,
                    wavelength_nm=lam, musp_780=mu, seed=seed)
        for j in range(len(C.MUA_PER_MM)):
            r, _ = R2.chip_rows(base, None, None, lambda NA, jj, s=s1[place]: A_slab_bg(s, NA, jj), R2.GLARE[side],
                                cams if j == J else ["scan_NA0.25"], lam, C.SPOT_SIZES_UM if j == J else [100],
                                C.REFL_V3 if j == J else ["paper635"],
                                C.IRRADIANCE_MW_CM2 if j == J else [10.0], rng, j=j)
            rows += r
    for k, d in enumerate(C.DEPTHS_V2_UM):
        X = S2.stage2_slab(slabC, d / 1000, s1["C_embedded"], k, C.N2_SLAB // F, sd(12, li, mi, seed, k))
        base = dict(geometry="slab", placement="C_embedded", D_mm=C.SLAB_THICKNESS_MM, depth_um=d,
                    wavelength_nm=lam, musp_780=mu, seed=seed)
        for j in range(len(C.MUA_PER_MM)):
            r, _ = R2.chip_rows(base, X, X["irr"], lambda NA, jj: A_slab_bg(s1["C_embedded"], NA, jj),
                                R2.GLARE["above"], cams if j == J else ["scan_NA0.25"], lam,
                                C.SPOT_SIZES_UM if j == J else [100], C.REFL_V3 if j == J else ["paper635"],
                                C.IRRADIANCE_MW_CM2 if j == J else [10.0], rng, j=j)
            rows += r
    return dict(kind="slab", rows=rows, sec=time.time() - t0, tag=f"slab lam={lam} mu={mu} seed={seed}")


def A_slab_bg(s, NA, j):
    return R2.A.slab_bg(s, NA, j)


def sphere_job(args):
    Di, mi, li, si = args
    D, mu, lam, seed = C.SPHERE_DIAMETERS_MM[Di], C.MUSP_780_PER_MM[mi], C.WAVELENGTHS_V3[li], C.SEEDS_V2[si]
    t0 = time.time()
    rng = np.random.default_rng(sd(40, Di, mi, li, seed))
    sph = S2.tissue_sphere(mu, lam, D)
    R = D / 2
    s1 = S2.stage1_sphere(sph, [R], C.N1_SPHERE // F, sd(41, Di, mi, li, seed))
    rows = []
    tilts = C.CENTRE_TILTS_DEG if lam in C.CENTRE_TILT_WAVELENGTHS else [0]
    for ti, tilt in enumerate(tilts):
        n2 = (C.N2_SPHERE if tilt == 0 else C.N2_TILT_SPHERE) // F
        X = S2.stage2_sphere(sph, R, s1, 0, n2, sd(42, Di, mi, li, seed, ti), tilt_deg=tilt)
        for side in ("above", "below"):
            base = dict(geometry="sphere", placement="center", side=side, D_mm=D, depth_um=round(R * 1000),
                        wavelength_nm=lam, musp_780=mu, seed=seed)
            if tilt == 0:
                cams, refl = C.HEADLINE_CAMERAS, C.REFL_V3
            else:
                cams, refl = ["scan_NA0.25"], ["paper635"]
            js = range(len(C.MUA_PER_MM)) if (tilt == 0 and D == 1.0) else [J]
            for j in js:
                r, _ = R2.chip_rows(base, X, X["irr"], lambda NA, jj: R2.A.sphere_bg(s1, R, NA, jj), R2.GLARE[side],
                                    cams if j == J else ["scan_NA0.25"], lam, C.SPOT_SIZES_UM if j == J else [100],
                                    refl if j == J else ["paper635"],
                                    C.IRRADIANCE_MW_CM2 if (tilt == 0 and j == J) else [10.0], rng,
                                    tilt=tilt, rho_max_um=R * 1000, j=j)
                rows += r
    return dict(kind="sphere", rows=rows, sec=time.time() - t0, tag=f"sphere D={D} mu={mu} lam={lam} seed={seed}")


KEYS = R2.KEYS + ["mua"]


def aggregate(rows):
    g = defaultdict(list)
    for r in rows:
        g[tuple(r.get(k, "") for k in KEYS)].append(r)
    out = []
    for key, rs in g.items():
        a = dict(zip(KEYS, key)); a["n_seeds"] = len(rs)
        for c in rs[0]:
            if c in KEYS or c == "seed" or isinstance(rs[0][c], str):
                continue
            v = np.array([r[c] for r in rs], float)
            fin = v[np.isfinite(v)]
            a[c] = float(fin.mean()) if fin.size else float("inf")
            a[c + "_err"] = float(fin.std(ddof=1) / math.sqrt(fin.size)) if fin.size > 1 else float("nan")
        out.append(a)
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    sys.path.insert(0, os.path.join(HERE, "tests"))
    import test_validation as TV
    import test_v2 as T2
    import test_v3 as T3
    rep = TV.run_all_checks() + T2.run_all_checks() + T3.run_all_checks()
    print(rep, flush=True)
    open(os.path.join(OUT, "validation_v3.txt"), "w").write(rep)
    if "[FAIL]" in rep:
        sys.exit("validation failed")
    jobs = [("s", (Di, mi, li, si)) for si in range(3) for Di in reversed(range(3)) for mi in reversed(range(3))
            for li in range(3)] + [("l", (li, mi, si)) for si in range(3) for mi in reversed(range(3)) for li in range(3)]
    workers = max(1, (os.cpu_count() or 2) - 1)
    rows_slab, rows_sph = [], []
    print(f"== Monte Carlo: {len(jobs)} jobs ==", flush=True)
    with Pool(workers) as pool:
        asy = [pool.apply_async(sphere_job if k == "s" else slab_job, (a,)) for k, a in jobs]
        for a in asy:
            r = a.get()
            print(f"  {r['tag']}: {r['sec']:.0f}s", flush=True)
            (rows_slab if r["kind"] == "slab" else rows_sph).extend(r["rows"])
    res, sres = aggregate(rows_slab), aggregate(rows_sph)
    pickle.dump(dict(res=res, sres=sres), open(os.path.join(OUT, "mc_v3.pkl"), "wb"))
    R2.write_csv(os.path.join(OUT, "per_chip_v3.csv"), res)
    R2.write_csv(os.path.join(OUT, "sphere_v3.csv"), sres)
    print(f"   MC done at {time.time() - t0:.0f}s", flush=True)
    import make_tables_v3
    make_tables_v3.main(res, sres, OUT)
    import assets_v3
    assets_v3.main(res, OUT)
    with open(os.path.join(OUT, "run_info.txt"), "w") as f:
        f.write(f"total_seconds {time.time() - t0:.0f}\nquick {QUICK}\nN1_slab {C.N1_SLAB // F}\nN2_slab {C.N2_SLAB // F}\n"
                f"N1_sphere {C.N1_SPHERE // F}\nN2_sphere {C.N2_SPHERE // F}\nN2_tilt_sphere {C.N2_TILT_SPHERE // F}\n"
                f"seeds {C.SEEDS_V2}\norientations {C.N_ORIENTATIONS}\nworkers {workers}\n")
    print(f"== done in {time.time() - t0:.0f}s ==")


if __name__ == "__main__":
    main()
