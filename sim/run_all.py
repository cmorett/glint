"""Single entry point: validation, then the full parameter sweep, then tables.

    python run_all.py            # everything (~ a few minutes on a multi-core PC)
    python run_all.py --quick    # 10x fewer photons, for a smoke test
"""
import csv
import math
import os
import sys
import time
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

import config as C
import stages
import analysis as A

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")


def seed_for(*keys):
    return np.random.SeedSequence([20261117, *keys])


def job(args):
    li, mi, si, n1, n2 = args
    lam, musp780, seed = C.WAVELENGTHS_NM[li], C.MUSP_780_PER_MM[mi], C.SEEDS[si]
    t0 = time.time()
    s1 = {}
    for gi, (place, P) in enumerate(C.PLACEMENTS.items()):
        slab = stages.tissue_slab(musp780, lam, P["n_cam"], P["n_far"])
        depths_mm = [d / 1000 for d in P["depths_um"]]
        s1[place] = stages.stage1(slab, depths_mm, n1, seed_for(1, li, mi, seed, gi))
    slabC = stages.tissue_slab(musp780, lam, C.PLACEMENTS["C_embedded"]["n_cam"],
                               C.PLACEMENTS["C_embedded"]["n_far"])
    s2 = {}
    for k, d in enumerate(C.DEPTHS_C_UM):
        s2[d] = stages.stage2(slabC, d / 1000, s1["C_embedded"], k, n2, seed_for(2, li, mi, seed, k))
    rng = np.random.default_rng(seed_for(3, li, mi, seed))
    rows, shift_rows, blur_rows = A.analyse_setting(lam, musp780, s1, s2, rng, seed)
    diag = dict(wavelength_nm=lam, musp_780=musp780, seed=seed, seconds=time.time() - t0,
                lost_s1=max(s["lost"] for s in s1.values()),
                lost_s2=max(v["lost"] for v in s2.values()))
    return rows, shift_rows, blur_rows, diag


def chicken_job(args):
    si, n1, n2 = args
    ck = C.CHICKEN
    seed = C.SEEDS[si]
    slab = stages.tissue_slab(ck["musp_780"], 780, ck["n_cam"], ck["n_far"], mua=ck["mua"],
                              thickness_mm=ck["thickness_mm"])
    d_mm = ck["depth_um"] / 1000
    s1 = stages.stage1(slab, [d_mm], n1, seed_for(4, seed))
    s2 = stages.stage2(slab, d_mm, s1, 0, n2, seed_for(5, seed))
    j = list(C.MUA_PER_MM).index(C.MUA_BASELINE)
    rng = np.random.default_rng(seed_for(6, seed))
    out = []
    rr = np.hypot(s2["x"], s2["y"]) * 1000
    r50p, r90p = A.weighted_quantiles(rr, s2["w"][j].astype(float), [0.5, 0.9])
    for na_i, na in enumerate(C.NA_LIST):
        sel = s2["inv"] <= na
        w = s2["w"][:, sel].astype(float)
        Cna = w[j].sum() / s2["N"]
        F = s1["F_down"][j, 0]
        bg = s1["r_sp"] + s1["Rd_na"][j, na_i]
        for a in ck["spots_um"]:
            img = A.SpotImage(s2["x_bp"][sel] * 1000.0, s2["y_bp"][sel] * 1000.0, w, a, rng)
            st = img.stats(j)
            for E in C.IRRADIANCE_MW_CM2:
                K = A.electrons_per_um2_per_frame(E, 780)
                res = A.snr_models(K, a * a * F * Cna, C.R_MEAN + 0.05, C.R_MEAN - 0.05,
                                   st["EE"], st["peak"], bg)
                out.append(dict(seed=seed, NA=na, spot_um=a, irradiance_mWcm2=E, F_down=F, C_NA=Cna,
                                point_r50_surface_um=r50p, point_r90_surface_um=r90p,
                                img_r50_um=st["r50"], img_r90_um=st["r90"], m=res["m"],
                                snr_cam=res["snr_cam"], snr_ideal=res["snr_ideal"]))
    return out


# ------------------------------------------------------------------ aggregation
def aggregate(rows, key_cols, snr_pairs=()):
    groups = defaultdict(list)
    for r in rows:
        groups[tuple(r[k] for k in key_cols)].append(r)
    out = []
    for key, rs in groups.items():
        agg = dict(zip(key_cols, key))
        agg["n_seeds"] = len(rs)
        for col in rs[0]:
            if col in key_cols or col == "seed" or isinstance(rs[0][col], str):
                continue
            v = np.array([r[col] for r in rs], float)
            fin = v[np.isfinite(v)]
            agg[col] = float(fin.mean()) if fin.size == v.size else (float(np.mean(v)) if fin.size == 0 else float(np.median(v)))
            agg[col + "_err"] = float(fin.std(ddof=1) / math.sqrt(fin.size)) if fin.size > 1 else float("nan")
        for snr_col, prefix in snr_pairs:   # times / verdicts from the seed-mean SNR
            s, e = agg[snr_col], agg[snr_col + "_err"]
            N = A.cycles_needed(s)
            agg["N_" + prefix] = N
            agg["N_" + prefix + "_err"] = 2 * N * e / s if s > 0 and np.isfinite(N) else float("nan")
            agg["t_" + prefix + "_s"] = N / C.BLINK_HZ
            agg["t_" + prefix + "_s_err"] = agg["N_" + prefix + "_err"] / C.BLINK_HZ
            agg["verdict_" + prefix] = A.verdict(N / C.BLINK_HZ)
        out.append(agg)
    return out


def write_csv(path, rows):
    cols = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def main():
    quick = "--quick" in sys.argv
    n1 = C.N_PHOTONS_STAGE1 // (10 if quick else 1)
    n2 = C.N_PHOTONS_STAGE2 // (10 if quick else 1)
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()

    print("== validation ==", flush=True)
    sys.path.insert(0, os.path.join(HERE, "tests"))
    import test_validation as TV
    report = TV.run_all_checks(verbose=True)
    with open(os.path.join(OUT, "validation.txt"), "w") as f:
        f.write(report)
    if "FAIL" in report:
        print("Validation failed; results would not be trustworthy. Stopping.")
        sys.exit(1)

    print(f"== sweep: {n1} stage-1 and {n2} stage-2 photons per run, seeds {C.SEEDS} ==", flush=True)
    jobs = [(li, mi, si, n1, n2) for si in range(len(C.SEEDS))
            for mi in reversed(range(len(C.MUSP_780_PER_MM))) for li in range(len(C.WAVELENGTHS_NM))]
    workers = C.N_WORKERS or max(1, (os.cpu_count() or 2) - 1)
    rows, shift_rows, blur_rows, diags = [], [], [], []
    with Pool(workers) as pool:
        ck_async = pool.map_async(chicken_job, [(si, n1, n2) for si in range(len(C.SEEDS))])
        for r, sr, br, dg in pool.imap_unordered(job, jobs):
            rows += r; shift_rows += sr; blur_rows += br; diags.append(dg)
            print(f"  done lam={dg['wavelength_nm']} musp780={dg['musp_780']} seed={dg['seed']} "
                  f"in {dg['seconds']:.0f}s (unfinished weight s1 {dg['lost_s1']:.1e}, s2 {dg['lost_s2']:.1e})",
                  flush=True)
        ck_rows = [r for rs in ck_async.get() for r in rs]

    keys = ["placement", "depth_um", "spot_um", "wavelength_nm", "musp_780", "mua", "NA",
            "irradiance_mWcm2", "dR"]
    res = aggregate(rows, keys, [("snr_cam", "cam"), ("snr_ideal", "ideal")])
    for r in res:
        r["below_contrast_floor"] = bool(r["m"] < C.CONTRAST_FLOOR)
    write_csv(os.path.join(OUT, "results.csv"), res)
    sres = aggregate(shift_rows, keys + ["filter", "efficiency"], [("snr_cam", "cam"), ("snr_ideal", "ideal")])
    write_csv(os.path.join(OUT, "results_wavelength_shift.csv"), sres)
    bres = aggregate(blur_rows, ["placement", "depth_um", "wavelength_nm", "musp_780", "mua"])
    write_csv(os.path.join(OUT, "blur_point_mirror.csv"), bres)
    cres = aggregate(ck_rows, ["NA", "spot_um", "irradiance_mWcm2"], [("snr_cam", "cam"), ("snr_ideal", "ideal")])
    write_csv(os.path.join(OUT, "sanity_chicken_skin.csv"), cres)

    import make_tables
    make_tables.main(res, sres, bres, cres)
    print(f"== done in {time.time() - t0:.0f}s; photons/run: stage1 {n1}, stage2 {n2} ==")
    with open(os.path.join(OUT, "run_info.txt"), "w") as f:
        f.write(f"total_seconds {time.time() - t0:.0f}\nphotons_stage1 {n1}\nphotons_stage2 {n2}\n"
                f"seeds {C.SEEDS}\nworkers {workers}\n")
        for d in diags:
            f.write(str(d) + "\n")


if __name__ == "__main__":
    main()
