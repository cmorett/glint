"""Predictions for the proposed bench experiment (results/v2/experiment_plan.md).

Geometry: a SPOT chip lying on a glass slide (n = 1.5) under a cast phantom slab of
thickness d (agarose + Intralipid in PBS, n ~ 1.34), read from above through a thin
PBS layer (n = 1.33). Phantom mu_s' = 1.0 /mm at 780 nm (+/- 20 % band), mu_a 0.003 /mm.
The chip sits AT the bottom face of the slab, so unlike the organoid slab there is no
tissue below it. Outputs modulation m vs d, relative to bare chip, for NA 0.25 and 0.02.

    python bench_prediction.py     (~1 min)
"""
import csv
import os

import numpy as np

import config_v2 as C
import analysis_v2 as A
import mc
import stages
import stages_v2 as S2

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "v2")
THICK_UM = [100, 200, 300, 500, 1000]
MUSP = {"low (0.8)": 0.8, "nominal (1.0)": 1.0, "high (1.2)": 1.2}
SPOTS = [30, 50, 100]
CAMS = ["scan_NA0.25", "plate20MP_NA0.02"]
N = 300_000
MUA = 0.003
N_PH = 1.34


def run():
    rows = []
    rng = np.random.default_rng(7)
    Ron, Roff = C.REFLECTANCE_CASES[C.HEADLINE_REFLECTANCE]
    glare = C.fresnel_normal(C.N_AIR, C.N_MEDIUM)
    To = 1 - glare
    for label, mu in MUSP.items():
        for d_um in THICK_UM:
            d = d_um / 1000
            ms = []
            for seed in C.SEEDS_V2:
                slab = mc.Slab(mus=mu / (1 - C.G), mua=MUA, g=C.G, T=d, n=N_PH, n_top=C.N_MEDIUM, n_bot=C.N_DISH)
                s1 = stages.stage1(slab, [d], N, seed, mua_list=[MUA, MUA, MUA], na_list=C.NA_LIST_V2)
                ms.append((s1, seed))
            for cn in CAMS:
                cam = C.CAMERAS[cn]
                for a in SPOTS:
                    vals = []
                    for s1, seed in ms:
                        Xs = _stage2(s1, d, mu, seed)
                        w, xb, yb = A.collected(Xs, cam["NA"], 1)
                        Cna = w.sum() / Xs["N"]
                        img = A.Image(w, xb * 1000, yb * 1000, a, 0, A.diffraction_sigma_um(780, cam["NA"]), rng)
                        bg = A.slab_bg(s1, cam["NA"], 1) * To * To + glare
                        K = A.electrons_per_um2_per_frame(10.0, 780, cam["fps"])
                        r = A.snr_cycle(cam, K, a * a * Xs["irr"][1] * Cna * To * To, Ron, Roff, img.EE(),
                                        img.peak_frac(cam["obj_pixel_um"]), bg)
                        vals.append((r["m"], img.quantile(0.5)))
                    v = np.array(vals)
                    # bare chip (no phantom, PBS only): no tissue in the path
                    r0 = A.snr_cycle(cam, A.electrons_per_um2_per_frame(10.0, 780, cam["fps"]),
                                     a * a * To * To * A.surface_collection(cam["NA"], 0, "specular"), Ron, Roff,
                                     A.Image(np.ones(1), np.zeros(1), np.zeros(1), a, 0,
                                             A.diffraction_sigma_um(780, cam["NA"]), rng).EE(),
                                     1.0 / (a * a), glare)
                    rows.append(dict(phantom_musp=label, thickness_um=d_um, camera=cn, spot_um=a,
                                     m=float(v[:, 0].mean()), m_err=float(v[:, 0].std(ddof=1) / np.sqrt(len(v))),
                                     m_bare=r0["m"], m_ratio=float(v[:, 0].mean() / r0["m"]),
                                     image_r50_um=float(v[:, 1].mean()), sharp_r50_um=0.399 * a))
            print(f"done {label} d={d_um}", flush=True)
    return rows


_cache = {}


def _stage2(s1, d, mu, seed):
    key = (d, mu, seed)
    if key not in _cache:
        slab = mc.Slab(mus=mu / (1 - C.G), mua=MUA, g=C.G, T=d, n=N_PH, n_top=C.N_MEDIUM, n_bot=C.N_DISH)
        rng = np.random.default_rng(seed + 100)
        rx, ry, rz, w, lb, mean_f = S2.launch(s1["F_ball"][1, 0], s1["theta_hist"][0], s1["theta_hist_up"][0], N, rng)
        start = dict(z=np.full(N, d), ux=rx, uy=ry, uz=rz, w0=w, launch_ballistic=lb)
        o = mc.run(slab, N, rng, mode="mbl", mua_list=[MUA, MUA, MUA], start=start, track_xy=True, record_exits=True)
        X = S2._pack_exits(o["exits"], d, [MUA, MUA, MUA], slab.n, slab.n)
        X.update(N=N, irr=(s1["F_down"][:, 0] + s1["F_up"][:, 0]) * mean_f)
        _cache[key] = X
    return _cache[key]


if __name__ == "__main__":
    rows = run()
    with open(os.path.join(OUT, "bench_prediction.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
    for r in rows:
        if r["phantom_musp"] == "nominal (1.0)":
            print(r)
