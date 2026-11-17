"""Round-2 checks. Run with `pytest tests/` (together with round-1 checks) or via run_v2.py."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config_v2 as C       # noqa: E402
import analysis_v2 as A     # noqa: E402
import decode_sim as DS     # noqa: E402
import decoder              # noqa: E402
import mc_sphere as MS      # noqa: E402
import stages               # noqa: E402
import stages_v2 as S2      # noqa: E402
import timesim              # noqa: E402

J = list(C.MUA_PER_MM).index(C.MUA_BASELINE)


def check_sphere_energy(n=100_000, seed=21):
    sph = S2.tissue_sphere(1.0, 780, 1.0)
    rng = np.random.default_rng(seed)
    start, r_sp = MS.collimated_start(sph, n, rng, r_core=0.15)
    o = MS.run(sph, start, rng, mua_list=C.MUA_PER_MM)
    w_imp_mean = (start["w0"].sum() / n) + r_sp
    tot = (r_sp + o["exit_total"][J] + o["A"] + o["lost"]) / w_imp_mean
    return abs(tot - 1) < 2e-3, f"Sphere energy conservation (D=1 mm): {tot:.5f}"


def check_sphere_vs_slab(seed=22):
    """A large sphere must reproduce the slab near its apex."""
    lam, mu, mua, d = 780, 2.0, 0.1, 0.2
    musp = C.musp_at(mu, lam)
    sph = MS.Sphere(mus=musp / (1 - C.G), mua=mua, g=C.G, R=5.0, n=C.N_TISSUE, n_out=C.N_MEDIUM)
    slab = stages.mc.Slab(mus=musp / (1 - C.G), mua=mua, g=C.G, T=10.0, n=C.N_TISSUE,
                          n_top=C.N_MEDIUM, n_bot=C.N_MEDIUM)
    rng = np.random.default_rng(seed)
    st, _ = MS.collimated_start(sph, 400_000, rng, r_core=0.5)
    o = MS.run(sph, st, rng, mua_list=[mua], planes=[d], r_score=0.1)
    F_sph = o["F_down"][0, 0] * (sph.R / 0.1) ** 2
    s1 = stages.stage1(slab, [d], 200_000, seed, mua_list=[mua], na_list=[0.25])
    F_slab = s1["F_down"][0, 0]
    # Stage 2 from a point mirror: collection into NA 0.25
    sc = (sph.R / 0.1) ** 2
    s1s = dict(F_ball=o["F_ball"] * sc, theta_hist=o["theta_hist"] * sc, theta_hist_up=o["theta_hist_up"] * sc,
               F_down=o["F_down"] * (sph.R / 0.1) ** 2, F_up=o["F_up"] * (sph.R / 0.1) ** 2)
    rng2 = np.random.default_rng(seed + 1)
    rx, ry, rz, w, lb, _ = S2.launch(s1s["F_ball"][0, 0], s1s["theta_hist"][0], s1s["theta_hist_up"][0], 200_000, rng2)
    st2 = dict(x=np.zeros(rx.size), y=np.zeros(rx.size), z=np.full(rx.size, d), ux=rx, uy=ry, uz=rz, w0=w)
    o2 = MS.run(sph, st2, rng2, mua_list=[mua])
    e = o2["exits"]
    inv = C.N_MEDIUM * np.sqrt(np.clip(1 - e["uz"] ** 2, 0, 1))
    C_sph = float((e["w0"] * e["m"] * np.exp(-mua * e["L"]))[inv <= 0.25].sum() / rx.size)
    s2 = stages.stage2(slab, d, s1, 0, 200_000, seed + 2, mua_list=[mua])
    C_slab = float(s2["w"][0][s2["inv"] <= 0.25].sum() / s2["N"])
    eF, eC = abs(F_sph / F_slab - 1), abs(C_sph / C_slab - 1)
    return eF < 0.05 and eC < 0.05, (f"Large sphere (D=10 mm) vs slab at d=0.2 mm: F_down {F_sph:.4f} vs {F_slab:.4f} "
                                     f"({eF:.1%}); C_NA0.25 {C_sph:.4f} vs {C_slab:.4f} ({eC:.1%})")


def check_low_na_scaling(seed=23):
    """Scaled low-NA collection vs direct counting (the scaling assumes flat angular density below NA 0.05)."""
    slab = stages.tissue_slab(1.0, 780, C.N_MEDIUM, C.N_DISH)
    s1 = S2.stage1_slab(slab, [0.1, 0.5], 300_000, seed)
    msgs, ok = [], True
    for k, d in enumerate([0.1, 0.5]):
        X = S2.stage2_slab(slab, d, s1, k, 1_500_000, seed + k + 1)
        for NA in (0.02,):
            w_sc, _, _ = A.collected(X, NA, J)
            direct = X["w"][J][X["inv"] <= NA].sum() / X["N"]
            nd = int(((~X["delta"]) & (X["inv"] <= NA)).sum())
            ratio = (w_sc.sum() / X["N"]) / direct
            ok &= abs(ratio - 1) < 0.3
            msgs.append(f"d={d * 1000:.0f} um NA {NA}: scaled/direct = {ratio:.3f} ({nd} non-ballistic photons counted directly)")
    # background (diffuse reflectance) scaling, slab Stage 1
    nas = list(C.NA_LIST_V2)
    big = stages.stage1(slab, [0.1], 2_000_000, seed + 9, mua_list=C.MUA_PER_MM, na_list=[0.02, 0.05])
    rd_direct = big["Rd_na"][J, 0]
    rd_scaled = big["Rd_na"][J, 1] * (0.02 / 0.05) ** 2
    rb = rd_scaled / rd_direct
    ok &= abs(rb - 1) < 0.15
    msgs.append(f"background Rd at NA 0.02: scaled/direct = {rb:.3f}")
    return ok, "Low-NA scaling: " + "; ".join(msgs)


def check_tilt_zero(seed=24):
    slab = stages.tissue_slab(1.0, 780, C.N_MEDIUM, C.N_DISH)
    s1 = S2.stage1_slab(slab, [0.2], 200_000, seed)
    a = stages.stage2(slab, 0.2, s1, 0, 200_000, seed + 1)
    b = S2.stage2_slab(slab, 0.2, s1, 0, 200_000, seed + 1)
    ca = a["w"][J][a["inv"] <= 0.25].sum() / a["N"]
    cb = b["w"][J][b["inv"] <= 0.25].sum() / b["N"] * b["irr"][J] / s1["F_down"][J, 0]
    e = abs(cb / ca - 1)
    return e < 0.02, f"Tilt=0 launch reproduces round-1 collection: {cb:.4f} vs {ca:.4f} ({e:.2%})"


def check_surface_collection():
    ok = abs(A.surface_collection(0.25, 0, "specular") - 1) < 1e-12
    ok &= A.surface_collection(0.02, 1.0, "specular") == 0.0        # 2*n*sin(2 deg) > 2*NA
    ok &= abs(A.surface_collection(0.25, 0, "lambert") - (0.25 / 1.33) ** 2) < 1e-12
    lim = np.degrees(0.5 * np.arcsin(2 * 0.02 / 1.33))
    return bool(ok), f"Surface-chip collection formulas; specular tilt limit at NA 0.02 = {lim:.2f} deg"


def check_decoder(seed=25):
    cam = C.CAMERAS["scan_NA0.25"]
    rng = np.random.default_rng(seed)
    pts = rng.normal(0, 15, (3000, 2))
    b = dict(e_on=2e5 * 1.05, e_off=2e5 * 0.95, e_bg_um2=30.0, roi_um=30.0, cam=cam)
    good = sum(DS.trial(dict(budget=b, pts=pts, n_c=4, sigma_src=1e-2, ref=True, seed=seed + i))["ok"] for i in range(3))
    # no chips at all: the decoder must not invent any
    b0 = dict(b, e_on=0.0, e_off=0.0)
    false = sum(DS.trial(dict(budget=b0, pts=pts, n_c=4, sigma_src=1e-2, ref=True, seed=seed + 10 + i))["false"]
                for i in range(3))
    return good == 9 and false == 0, f"Decoder: {good}/9 strong chips decoded; {false} false chips on chip-free frames"


def run_all_checks(verbose=False):
    out = []
    for f in (check_sphere_energy, check_sphere_vs_slab, check_low_na_scaling, check_tilt_zero,
              check_surface_collection, check_decoder):
        ok, msg = f()
        out.append(f"[{'PASS' if ok else 'FAIL'}] {msg}")
        if verbose:
            print(out[-1], flush=True)
    return "\n".join(out) + "\n"


def test_sphere_energy():
    assert check_sphere_energy()[0]


def test_sphere_vs_slab():
    assert check_sphere_vs_slab()[0]


def test_low_na_scaling():
    assert check_low_na_scaling()[0]


def test_tilt_zero():
    assert check_tilt_zero()[0]


def test_surface_collection():
    assert check_surface_collection()[0]


def test_decoder():
    assert check_decoder()[0]


if __name__ == "__main__":
    print(run_all_checks(verbose=True))
