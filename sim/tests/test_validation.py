"""Validation checks (prompt section 4). Run with `pytest tests/` or `python tests/test_validation.py`.

run_all.py calls run_all_checks() first and refuses to continue if any check FAILs.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as C   # noqa: E402
import mc            # noqa: E402
import stages        # noqa: E402
import analysis as A  # noqa: E402

# MCML benchmark (Wang, Jacques & Zheng 1995; Prahl et al. 1989; van de Hulst tables)
BENCH_RD, BENCH_TT = 0.09734, 0.66096


def check_mcml_benchmark(mode, n=400_000, seed=11):
    s = mc.Slab(mus=9.0, mua=1.0, g=0.75, T=0.2, n=1.0, n_top=1.0, n_bot=1.0)   # mm units
    o = mc.run(s, n, np.random.default_rng(seed), mode=mode)
    Rd, Tt = float(o["R_top"][0]), float(o["T_bot"][0])
    eR, eT = abs(Rd / BENCH_RD - 1), abs(Tt / BENCH_TT - 1)
    ok = eR < 0.01 and eT < 0.01
    return ok, f"MCML benchmark [{mode}] Rd={Rd:.5f} (ref {BENCH_RD}, err {eR:.2%}), " \
               f"Tt={Tt:.5f} (ref {BENCH_TT}, err {eT:.2%}), N={n}"


def check_zero_scattering(n=200_000, seed=12):
    mua, T, d = 0.5, 2.0, 0.5
    s = mc.Slab(mus=0.0, mua=mua, g=0.9, T=T, n=C.N_TISSUE, n_top=C.N_MEDIUM, n_bot=C.N_DISH)
    s1 = stages.stage1(s, [d], n, seed, mua_list=[mua], na_list=[0.1])
    s2 = stages.stage2(s, d, s1, 0, n, seed + 1, mua_list=[mua])
    rt = C.fresnel_normal(C.N_TISSUE, C.N_MEDIUM)
    rb = C.fresnel_normal(C.N_TISSUE, C.N_DISH)
    loop = 1.0 / (1.0 - rt * rb * np.exp(-2 * mua * T))
    F_an = (1 - rt) * np.exp(-mua * d) * loop
    C_an = (1 - rt) * np.exp(-mua * d) * loop
    F_mc = float(s1["F_down"][0, 0])
    C_mc = float(s2["w"][0][s2["inv"] <= 0.1].sum() / s2["N"])
    sig_mc, sig_an = F_mc * C_mc, F_an * C_an
    err = abs(sig_mc / sig_an - 1)
    # also: index-matched pure Beer-Lambert
    s0 = mc.Slab(mus=0.0, mua=mua, g=0.9, T=T, n=1.37, n_top=1.37, n_bot=1.37)
    F0 = float(stages.stage1(s0, [d], 20_000, seed, mua_list=[mua], na_list=[0.1])["F_down"][0, 0])
    err0 = abs(F0 / np.exp(-mua * d) - 1)
    ok = err < 2e-3 and err0 < 1e-9
    return ok, (f"Zero scattering: signal F*C MC={sig_mc:.6f} analytic={sig_an:.6f} (err {err:.3%}); "
                f"matched-index F_down MC={F0:.8f} vs exp(-mua d)={np.exp(-mua * d):.8f}")


def check_energy(mode, n=100_000, seed=13):
    s = stages.tissue_slab(C.MUSP_BASELINE, 780, C.N_MEDIUM, C.N_DISH)
    o = mc.run(s, n, np.random.default_rng(seed), mode=mode)
    tot = o["r_sp"] + o["R_top"][0] + o["T_bot"][0] + o["A"] + o["lost"]
    ok = abs(tot - 1) < 5e-3
    return ok, (f"Energy [{mode}] r_sp+R+T+A = {o['r_sp']:.5f}+{o['R_top'][0]:.5f}+{o['T_bot'][0]:.5f}"
                f"+{o['A']:.5f} = {tot:.5f}"), o


def check_modes_agree(o1, o2):
    d = abs(o1["R_top"][0] / o2["R_top"][0] - 1)
    ok = d < 0.01
    return ok, f"MCML vs absorption-weighting on mismatched tissue slab: Rd differ by {d:.2%}"


def check_chicken(n=60_000, seed=14):
    """Section 4.4. Reported as CHECK/NOTE, not a hard failure (we report, not tune)."""
    ck = C.CHICKEN
    s = stages.tissue_slab(ck["musp_780"], 780, ck["n_cam"], ck["n_far"], mua=ck["mua"],
                           thickness_mm=ck["thickness_mm"])
    d = ck["depth_um"] / 1000
    s1 = stages.stage1(s, [d], n, seed)
    s2 = stages.stage2(s, d, s1, 0, n, seed + 1)
    j = list(C.MUA_PER_MM).index(C.MUA_BASELINE)
    na_i = C.NA_LIST.index(0.25)
    sel = s2["inv"] <= 0.25
    w = s2["w"][:, sel].astype(float)
    Cna = w[j].sum() / s2["N"]
    F = s1["F_down"][j, 0]
    bg = s1["r_sp"] + s1["Rd_na"][j, na_i]
    rr = np.hypot(s2["x"], s2["y"]) * 1000
    r50p, _ = A.weighted_quantiles(rr, s2["w"][j].astype(float), [0.5, 0.9])
    lines = []
    readable = {}
    for a in ck["spots_um"]:
        img = A.SpotImage(s2["x_bp"][sel] * 1000.0, s2["y_bp"][sel] * 1000.0, w, a,
                          np.random.default_rng(seed))
        st = img.stats(j)
        K = A.electrons_per_um2_per_frame(10.0, 780)
        res = A.snr_models(K, a * a * F * Cna, 0.55, 0.45, st["EE"], st["peak"], bg)
        t = A.cycles_needed(res["snr_cam"]) / C.BLINK_HZ
        readable[a] = t < C.READABLE_S
        lines.append(f"    a={a} um: image r50={st['r50']:.0f} um (sharp square r50 ~{0.399 * a:.0f} um), "
                     f"m={res['m']:.3g}, SNR/cycle (camera)={res['snr_cam']:.3g}, time to SNR 5 = {t:.3g} s")
    ok = readable[300] and r50p > 20
    head = (f"Chicken-skin sanity (250 um, mu_s'={ck['musp_780']:.2f}/mm, 780 nm, NA 0.25, 10 mW/cm2, dR 0.1): "
            f"F_down={F:.3f}, C_NA={Cna:.4f}, point-mirror r50 at surface={r50p:.0f} um")
    return ok, head + "\n" + "\n".join(lines)


def run_all_checks(verbose=False):
    out = []

    def rec(ok, msg, hard=True):
        tag = "PASS" if ok else ("FAIL" if hard else "CHECK-NOT-MET")
        out.append(f"[{tag}] {msg}")
        if verbose:
            print(out[-1], flush=True)

    rec(*check_mcml_benchmark("mcml"))
    rec(*check_mcml_benchmark("mbl"))
    rec(*check_zero_scattering())
    ok1, m1, o1 = check_energy("mcml")
    rec(ok1, m1)
    ok2, m2, o2 = check_energy("mbl")
    rec(ok2, m2)
    rec(*check_modes_agree(o1, o2))
    ok, msg = check_chicken()
    rec(ok, msg, hard=False)
    out.append("[INFO] Statistical errors: every production configuration is run with seeds "
               f"{C.SEEDS}; results.csv reports mean and standard error across seeds.")
    return "\n".join(out) + "\n"


# ---------------- pytest entry points
def test_mcml_benchmark_mcml():
    assert check_mcml_benchmark("mcml")[0]


def test_mcml_benchmark_mbl():
    assert check_mcml_benchmark("mbl")[0]


def test_zero_scattering():
    assert check_zero_scattering()[0]


def test_energy_conservation():
    a, _, o1 = check_energy("mcml")
    b, _, o2 = check_energy("mbl")
    assert a and b and check_modes_agree(o1, o2)[0]


def test_fresnel():
    assert abs(mc.fresnel(1.0, 1.5, np.array([1.0]))[0] - 0.04) < 1e-12
    assert mc.fresnel(1.5, 1.0, np.array([0.5]))[0] == 1.0     # beyond critical angle


if __name__ == "__main__":
    print(run_all_checks(verbose=False))
