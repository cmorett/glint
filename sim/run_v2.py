"""Round 2: whole-plate optics, spheres, tilt, decoded noise floor, light dose, measured dR.

    python run_v2.py            # everything; writes results/v2/ (round-1 results untouched)
    python run_v2.py --quick    # 10x fewer photons and fewer decoder trials (smoke test)
"""
import csv
import math
import os
import pickle
import sys
import time
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

import config_v2 as C
import analysis_v2 as A
import decode_sim as DS
import stages
import stages_v2 as S2

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "v2")
QUICK = "--quick" in sys.argv
F = 10 if QUICK else 1
J = list(C.MUA_PER_MM).index(C.MUA_BASELINE)
GLARE = {"above": C.fresnel_normal(C.N_AIR, C.N_MEDIUM), "below": C.fresnel_normal(C.N_AIR, C.N_DISH)}


def seed_for(*k):
    return np.random.SeedSequence([20261117, 2, *k])


# ------------------------------------------------------------------ row builder
def chip_rows(base, X, irr, bg_fn, glare, cams, lam, a_list, refl_cases, E_list, rng,
              tilt=0.0, mode="specular", rho_max_um=None, keep_pts=False):
    """Rows for one chip geometry. X = Stage-2 exits (None = no tissue in the path)."""
    To = 1.0 - glare
    rows, pts = [], {}
    for cam_name in cams:
        cam = C.CAMERAS[cam_name]
        NA = cam["NA"]
        if X is None:
            Cna = A.surface_collection(NA, tilt, mode)
            w, xb, yb = np.ones(1), np.zeros(1), np.zeros(1)
            irr_j = 1.0
        else:
            w, xb, yb = A.collected(X, NA, J)
            Cna = w.sum() / X["N"]
            irr_j = float(irr[J])
        bgd = bg_fn(NA) * To * To + glare
        sig = A.diffraction_sigma_um(lam, NA)
        for a in a_list:
            img = A.Image(w, xb * 1000.0, yb * 1000.0, a, tilt, sig, rng)
            EE, peak = img.EE(), img.peak_frac(cam["obj_pixel_um"])
            sig_frac = a * a * irr_j * Cna * To * To
            if keep_pts:
                pts[(cam_name, a)] = DS.image_points(img, rng).astype(np.float32)
            for rc in refl_cases:
                Ron, Roff = C.REFLECTANCE_CASES[rc]
                for E in E_list:
                    K = A.electrons_per_um2_per_frame(E, lam, cam["fps"])
                    r = A.snr_cycle(cam, K, sig_frac, Ron, Roff, EE, peak, bgd, rho_max_um=rho_max_um)
                    cyc = A.cycles_to(r["snr1"])
                    rows.append(dict(base, cam=cam_name, NA=NA, spot_um=a, refl=rc, irradiance=E,
                                     tilt=tilt, mode=mode, irr=irr_j, C_NA=Cna, sig_frac_um2=sig_frac,
                                     bg_dens=bgd, img_r50_um=img.quantile(0.5), img_r90_um=img.quantile(0.9),
                                     K=K, blink_hz=cam["blink_hz"], t_snr5_s=cyc / cam["blink_hz"], **r))
    return rows, pts


# ------------------------------------------------------------------ jobs
def slab_job(args):
    li, mi, si = args
    lam, mu, seed = C.WAVELENGTHS_V2[li], C.MUSP_780_PER_MM[mi], C.SEEDS_V2[si]
    t0 = time.time()
    rng = np.random.default_rng(seed_for(10, li, mi, seed))
    s1 = {}
    for gi, (place, P) in enumerate(C.PLACEMENTS.items()):
        slab = stages.tissue_slab(mu, lam, P["n_cam"], P["n_far"])
        depths = [d / 1000 for d in (C.DEPTHS_V2_UM if place == "C_embedded" else [0])]
        s1[place] = S2.stage1_slab(slab, depths, C.N1_SLAB // F, seed_for(11, li, mi, seed, gi))
    slabC = stages.tissue_slab(mu, lam, C.N_MEDIUM, C.N_DISH)
    keep = (lam == C.HEADLINE_WAVELENGTH and mu == C.MUSP_BASELINE and si == 0)
    rows, pts = [], {}
    cams = list(C.CAMERAS)
    for place in ("A_top", "B_under"):
        side = "above" if place == "A_top" else "below"
        base = dict(geometry="slab", placement=place, D_mm=C.SLAB_THICKNESS_MM, depth_um=0,
                    wavelength_nm=lam, musp_780=mu, seed=seed)
        r, p = chip_rows(base, None, None, lambda NA, s=s1[place]: A.slab_bg(s, NA, J), GLARE[side],
                         cams, lam, C.SPOT_SIZES_UM, list(C.REFLECTANCE_CASES), C.IRRADIANCE_MW_CM2,
                         rng, keep_pts=keep)
        rows += r; pts.update({(place, 0) + k: v for k, v in p.items()})
    for k, d in enumerate(C.DEPTHS_V2_UM):
        X = S2.stage2_slab(slabC, d / 1000, s1["C_embedded"], k, C.N2_SLAB // F, seed_for(12, li, mi, seed, k))
        base = dict(geometry="slab", placement="C_embedded", D_mm=C.SLAB_THICKNESS_MM, depth_um=d,
                    wavelength_nm=lam, musp_780=mu, seed=seed,
                    point_r50_surface_um=_wq(np.hypot(X["x"], X["y"]) * 1000, X["w"][J], 0.5))
        r, p = chip_rows(base, X, X["irr"], lambda NA: A.slab_bg(s1["C_embedded"], NA, J), GLARE["above"],
                         cams, lam, C.SPOT_SIZES_UM, list(C.REFLECTANCE_CASES), C.IRRADIANCE_MW_CM2,
                         rng, keep_pts=keep)
        rows += r; pts.update({("C_embedded", d) + kk: v for kk, v in p.items()})
    return dict(rows=rows, pts=pts, sec=time.time() - t0, tag=f"slab lam={lam} mu={mu} seed={seed}")


def _wq(r, w, q):
    w = np.asarray(w, float)
    if r.size == 0 or w.sum() <= 0:
        return float("nan")
    o = np.argsort(r); cw = np.cumsum(w[o]) / w.sum()
    return float(r[o][min(np.searchsorted(cw, q), r.size - 1)])


def tilt_job(args):
    di, si = args
    d, seed = C.TILT_DEPTHS_UM[di], C.SEEDS_V2[si]
    t0 = time.time()
    lam, mu = C.HEADLINE_WAVELENGTH, C.MUSP_BASELINE
    rng = np.random.default_rng(seed_for(20, di, seed))
    slab = stages.tissue_slab(mu, lam, C.N_MEDIUM, C.N_DISH)
    s1 = S2.stage1_slab(slab, [d / 1000], C.N1_SLAB // (3 * F), seed_for(21, di, seed))
    rows = []
    cams = list(C.CAMERAS)
    jobs = [(t, "specular") for t in C.TILTS_DEG] + [(t, "lambert") for t in C.LAMBERT_TILTS_DEG]
    for ti, (tilt, mode) in enumerate(jobs):
        X = S2.stage2_slab(slab, d / 1000, s1, 0, C.N2_TILT // F, seed_for(22, di, seed, ti), tilt, mode)
        base = dict(geometry="slab", placement="C_embedded", D_mm=C.SLAB_THICKNESS_MM, depth_um=d,
                    wavelength_nm=lam, musp_780=mu, seed=seed)
        r, _ = chip_rows(base, X, X["irr"], lambda NA: A.slab_bg(s1, NA, J), GLARE["above"], cams, lam,
                         [100], [C.HEADLINE_REFLECTANCE], [C.IRRADIANCE_BASELINE], rng, tilt=tilt, mode=mode)
        rows += r
    return dict(rows=rows, sec=time.time() - t0, tag=f"tilt d={d} seed={seed}")


def surface_tilt_rows():
    """No tissue in the path (A on top, B under / SPOT facing the glass): analytic."""
    rows = []
    rng = np.random.default_rng(seed_for(30))
    lam, mu = C.HEADLINE_WAVELENGTH, C.MUSP_BASELINE
    slab = stages.tissue_slab(mu, lam, C.N_DISH, C.N_MEDIUM)
    s1 = S2.stage1_slab(slab, [0.0], C.N1_SLAB // (3 * F), seed_for(31))
    for tilt in sorted(set(C.TILTS_DEG + [0.25, 0.5, 1, 2, 3])):
        for mode in ("specular", "lambert"):
            base = dict(geometry="slab", placement="B_under", D_mm=C.SLAB_THICKNESS_MM, depth_um=0,
                        wavelength_nm=lam, musp_780=mu, seed=1)
            r, _ = chip_rows(base, None, None, lambda NA: A.slab_bg(s1, NA, J), GLARE["below"],
                             list(C.CAMERAS), lam, C.SPOT_SIZES_UM, [C.HEADLINE_REFLECTANCE],
                             [C.IRRADIANCE_BASELINE], rng, tilt=tilt, mode=mode)
            rows += r
    return rows


def sphere_job(args):
    Di, mi, si = args
    D, mu, seed = C.SPHERE_DIAMETERS_MM[Di], C.MUSP_780_PER_MM[mi], C.SEEDS_V2[si]
    t0 = time.time()
    lam = C.HEADLINE_WAVELENGTH
    rng = np.random.default_rng(seed_for(40, Di, mi, seed))
    sph = S2.tissue_sphere(mu, lam, D)
    R = D / 2
    pos = {"halfway": R / 2, "center": R, "far_surface": D - C.FAR_SIDE_GAP_MM}
    s1 = S2.stage1_sphere(sph, list(pos.values()), C.N1_SPHERE // F, seed_for(41, Di, mi, seed))
    rows = []
    cams = C.HEADLINE_CAMERAS
    for side in ("above", "below"):
        base = dict(geometry="sphere", placement="surface", side=side, D_mm=D, depth_um=0,
                    wavelength_nm=lam, musp_780=mu, seed=seed)
        r, _ = chip_rows(base, None, None, lambda NA: A.sphere_bg(s1, 0.0, NA, J), GLARE[side], cams, lam,
                         C.SPOT_SIZES_UM, [C.HEADLINE_REFLECTANCE], [C.IRRADIANCE_BASELINE], rng,
                         rho_max_um=R * 1000)
        rows += r
    for k, (name, d) in enumerate(pos.items()):
        X = S2.stage2_sphere(sph, d, s1, k, C.N2_SPHERE // F, seed_for(42, Di, mi, seed, k))
        for side in ("above", "below"):
            base = dict(geometry="sphere", placement=name, side=side, D_mm=D, depth_um=round(d * 1000),
                        wavelength_nm=lam, musp_780=mu, seed=seed, F_down_sphere=float(s1["F_down"][J, k]))
            r, _ = chip_rows(base, X, X["irr"], lambda NA, d=d: A.sphere_bg(s1, d, NA, J), GLARE[side], cams,
                             lam, C.SPOT_SIZES_UM, [C.HEADLINE_REFLECTANCE], [C.IRRADIANCE_BASELINE], rng,
                             rho_max_um=R * 1000)
            rows += r
    return dict(rows=rows, sec=time.time() - t0, tag=f"sphere D={D} mu={mu} seed={seed}",
                energy=float(s1["r_sp"] + s1["exit_total"][J] + s1["A"] + s1["lost"]))


# ------------------------------------------------------------------ aggregation
KEYS = ["geometry", "placement", "side", "D_mm", "depth_um", "wavelength_nm", "musp_780", "cam", "NA",
        "spot_um", "refl", "irradiance", "tilt", "mode"]


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


def write_csv(path, rows):
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


# ------------------------------------------------------------------ decoder (Task 4)
def floor_search(args):
    """Smallest modulation (as a multiple m_scale of the budget's own) that decodes >= 23/24 chips."""
    budget, pts, n_c, sigma, ref, trials, seed0 = args
    lo, hi = -4.5, 0.5          # log10(m_scale)
    def ok_at(lg, salt):
        good = tot = 0
        for t in range(trials):
            r = DS.trial(dict(budget=budget, pts=pts, n_c=n_c, sigma_src=sigma, ref=ref,
                              seed=seed0 * 1000 + salt * 50 + t, m_scale=10 ** lg))
            good += r["ok"]; tot += r["n"]
        return good >= tot - 1
    if not ok_at(hi, 0):
        return dict(n_c=n_c, sigma=sigma, ref=ref, m_scale_min=float("inf"))
    for i in range(7):
        mid = 0.5 * (lo + hi)
        if ok_at(mid, i + 1):
            hi = mid
        else:
            lo = mid
    return dict(n_c=n_c, sigma=sigma, ref=ref, m_scale_min=10 ** hi)


def ber_job(args):
    budget, pts, n_c, trials, seed0, tag = args
    ok = n = be = bits = false = 0
    T = 0
    for t in range(trials):
        r = DS.trial(dict(budget=budget, pts=pts, n_c=n_c, sigma_src=1e-3, ref=True, seed=seed0 * 100 + t))
        ok += r["ok"]; n += r["n"]; be += r["bit_err"]; bits += r["bits"]; false += r["false"]; T = r["T_s"]
    return dict(tag, n_c=n_c, read_s=T, frames_ok=ok, frames=n, ber=be / bits, false_chips=false)


def budget_for(row):
    cam = C.CAMERAS[row["cam"]]
    Ron, Roff = C.REFLECTANCE_CASES[row["refl"]]
    return dict(e_on=row["sat_k"] * row["K"] * row["sig_frac_um2"] * Ron,
                e_off=row["sat_k"] * row["K"] * row["sig_frac_um2"] * Roff,
                e_bg_um2=row["sat_k"] * row["K"] * row["bg_dens"], roi_um=row["roi_um"], cam=cam)


def main():
    os.makedirs(OUT, exist_ok=True)
    t_start = time.time()
    print("== round-2 validation ==", flush=True)
    sys.path.insert(0, os.path.join(HERE, "tests"))
    import test_validation as TV
    import test_v2 as T2
    rep = TV.run_all_checks() + T2.run_all_checks()
    print(rep, flush=True)
    open(os.path.join(OUT, "validation_v2.txt"), "w").write(rep)
    if "[FAIL]" in rep:
        sys.exit("validation failed")

    workers = max(1, (os.cpu_count() or 2) - 1)
    jobs_slab = [(li, mi, si) for si in range(len(C.SEEDS_V2)) for mi in reversed(range(3))
                 for li in range(len(C.WAVELENGTHS_V2))]
    jobs_tilt = [(di, si) for si in range(len(C.SEEDS_V2)) for di in range(len(C.TILT_DEPTHS_UM))]
    jobs_sph = [(Di, mi, si) for si in range(len(C.SEEDS_V2)) for Di in reversed(range(3)) for mi in reversed(range(3))]
    rows, trows, srows, pts = [], [], [], {}
    print(f"== Monte Carlo: {len(jobs_slab)} slab, {len(jobs_tilt)} tilt, {len(jobs_sph)} sphere jobs ==", flush=True)
    with Pool(workers) as pool:
        asyncs = [pool.apply_async(slab_job, (j,)) for j in jobs_slab] + \
                 [pool.apply_async(tilt_job, (j,)) for j in jobs_tilt] + \
                 [pool.apply_async(sphere_job, (j,)) for j in jobs_sph]
        for a in asyncs:
            r = a.get()
            print(f"  {r['tag']}: {r['sec']:.0f}s" + (f" (sphere energy {r['energy']:.4f})" if "energy" in r else ""),
                  flush=True)
            if r["tag"].startswith("slab"):
                rows += r["rows"]; pts.update(r["pts"])
            elif r["tag"].startswith("tilt"):
                trows += r["rows"]
            else:
                srows += r["rows"]
    trows += surface_tilt_rows()
    res, tres, sres = aggregate(rows), aggregate(trows), aggregate(srows)
    write_csv(os.path.join(OUT, "per_chip.csv"), res)
    write_csv(os.path.join(OUT, "tilt.csv"), tres)
    write_csv(os.path.join(OUT, "sphere.csv"), sres)
    pickle.dump(dict(res=res, tres=tres, sres=sres), open(os.path.join(OUT, "mc_v2.pkl"), "wb"))
    print(f"   MC done at {time.time() - t_start:.0f}s", flush=True)

    # ---------------- Task 4: decoder floors and BER
    def find(**kw):
        for r in res:
            if all(r.get(k) == v for k, v in kw.items()):
                return r
        raise KeyError(kw)
    head = dict(wavelength_nm=C.HEADLINE_WAVELENGTH, musp_780=C.MUSP_BASELINE, refl=C.HEADLINE_REFLECTANCE,
                irradiance=C.IRRADIANCE_BASELINE)
    ref_row = find(placement="C_embedded", depth_um=300, spot_um=100, cam="scan_NA0.25", **head)
    plate_row = find(placement="A_top", depth_um=0, spot_um=100, cam="plate20MP_NA0.02", **head)
    trials = 3 if QUICK else C.DECODE_TRIALS
    fs_jobs = []
    for n_c in ([10, 50] if QUICK else [10, 50, 290]):
        for sg in C.SOURCE_NOISE_RMS:
            for ref in (True, False):
                fs_jobs.append(("scan_NA0.25 a=100 d=300", ref_row, (budget_for(ref_row),
                                pts[("C_embedded", 300, "scan_NA0.25", 100)], n_c, sg, ref, trials, len(fs_jobs) + 1)))
    for n_c in ([2] if QUICK else [2, 9]):
        for sg in (1e-3, 1e-2):
            for ref in (True, False):
                fs_jobs.append(("plate20MP_NA0.02 a=100 surface", plate_row, (budget_for(plate_row),
                                pts[("A_top", 0, "plate20MP_NA0.02", 100)], n_c, sg, ref, trials, len(fs_jobs) + 1)))
    ber_cfgs = []
    for (pl, d) in [("A_top", 0), ("B_under", 0), ("C_embedded", 100), ("C_embedded", 300),
                    ("C_embedded", 500), ("C_embedded", 1000)]:
        for a in (30, 100):
            for cam, ncs in (("scan_NA0.25", [2, 10, 50, 290]), ("plate20MP_NA0.02", [2, 9])):
                r = find(placement=pl, depth_um=d, spot_um=a, cam=cam, **head)
                for n_c in (ncs[:2] if QUICK else ncs):
                    ber_cfgs.append((budget_for(r), pts[(pl, d, cam, a)], n_c, 2 if QUICK else 4,
                                     len(ber_cfgs) + 1, dict(placement=pl, depth_um=d, spot_um=a, cam=cam,
                                                             m=r["m"], snr1=r["snr1"])))
    print(f"== decoder: {len(fs_jobs)} floor searches, {len(ber_cfgs)} BER configs ==", flush=True)
    with Pool(workers) as pool:
        fa = pool.map_async(floor_search, [j[2] for j in fs_jobs])
        ba = pool.map_async(ber_job, ber_cfgs)
        floors, bers = fa.get(), ba.get()
    frows = []
    for (name, row, _), f in zip(fs_jobs, floors):
        cam = C.CAMERAS[row["cam"]]
        snr_sym = row["snr1"] * f["m_scale_min"] * math.sqrt(f["n_c"])
        frows.append(dict(config=name, n_c=f["n_c"], read_s=(C.FRAME_BITS + 2) * f["n_c"] / cam["blink_hz"],
                          source_noise_rms=f["sigma"], reference=f["ref"], m_budget=row["m"],
                          m_min=row["m"] * f["m_scale_min"], snr1_at_min=row["snr1"] * f["m_scale_min"],
                          snr_symbol_ideal_at_min=snr_sym))
    write_csv(os.path.join(OUT, "decoder_floor.csv"), frows)
    write_csv(os.path.join(OUT, "decoder_ber.csv"), bers)
    print(f"   decoder done at {time.time() - t_start:.0f}s", flush=True)

    import make_tables_v2
    make_tables_v2.main(res, tres, sres, frows, bers, OUT)
    with open(os.path.join(OUT, "run_info.txt"), "w") as f:
        f.write(f"total_seconds {time.time() - t_start:.0f}\nquick {QUICK}\nN1_slab {C.N1_SLAB // F}\n"
                f"N2_slab {C.N2_SLAB // F}\nN2_tilt {C.N2_TILT // F}\nN1_sphere {C.N1_SPHERE // F}\n"
                f"N2_sphere {C.N2_SPHERE // F}\nseeds {C.SEEDS_V2}\ndecoder_trials {trials}\nworkers {workers}\n")
    print(f"== done in {time.time() - t_start:.0f}s ==")


if __name__ == "__main__":
    main()
