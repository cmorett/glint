"""Glue between Monte Carlo photon budgets and the synthetic-frame decoder (SYNTHETIC DATA)."""
import numpy as np

import config_v2 as C
import decoder
import timesim

GRID = (24, 24)
CHIP_SITES = [(6.0, 6.0), (6.0, 18.0), (18.0, 12.0)]


def image_points(img, rng, M=4000):
    """2-D sample offsets (um) drawn from a radially-sorted analysis_v2.Image."""
    if img.tot <= 0:
        return np.zeros((M, 2))
    p = img.w / img.tot
    r = img.r[rng.choice(img.r.size, size=M, p=p)]
    a = rng.uniform(0, 2 * np.pi, M)
    return np.stack([r * np.cos(a), r * np.sin(a)], axis=1)


def budget_from(res, cam, K, sig_frac, R_on, R_off, bg_dens):
    """Whole-image electrons per frame (after the full-well attenuation k in res)."""
    k = res["sat_k"]
    return dict(e_on=k * K * sig_frac * R_on, e_off=k * K * sig_frac * R_off,
                e_bg_um2=k * K * bg_dens, roi_um=res["roi_um"], cam=cam)


def trial(args):
    """One synthetic read of 3 chips. args: dict(budget, pts, T_s or n_c, sigma_src, ref, seed, m_scale)."""
    b = args["budget"]; cam = b["cam"]
    rng = np.random.default_rng(args["seed"])
    p = cam["obj_pixel_um"]
    sp = max(p, b["roi_um"] / 2.5)
    n_bin = (sp / p) ** 2
    fps, fb = cam["fps"], cam["blink_hz"]
    n_c = args["n_c"]
    T = (C.FRAME_BITS + 2) * n_c / fb          # 29 bits need 30 symbols, +1 for unknown timing
    e_mean = 0.5 * (b["e_on"] + b["e_off"])
    e_mod = (b["e_on"] - b["e_off"]) * args.get("m_scale", 1.0)
    chips = []
    for i, (yy, xx) in enumerate(CHIP_SITES):
        chips.append(dict(pos_um=(xx * sp + sp / 2, yy * sp + sp / 2), pts_um=args["pts"],
                          e_mean=e_mean, e_mod=e_mod, id=int(rng.integers(0, 256)),
                          value=int(rng.integers(0, 256))))
    sc = timesim.make_scene(chips, GRID, sp, fps, fb, n_c, b["e_bg_um2"],
                            cam["read_noise"] * np.sqrt(n_bin), rng, sigma_src=args["sigma_src"],
                            motion_um=args.get("motion_um", C.MOTION_UM))
    frames = timesim.render(sc, T, rng)
    radius = int(np.ceil(b["roi_um"] / sp)) + 1
    dec = decoder.decode(frames, fps, fb, n_c, use_reference=args["ref"], radius_px=radius)
    ok = 0; bit_err = 0; bits_tot = 0; false = 0
    used = set()
    for ch in sc["chips"]:
        cy, cx = (ch["pos_um"][1] - sp / 2) / sp, (ch["pos_um"][0] - sp / 2) / sp
        best = None
        for j, d in enumerate(dec):
            dist = np.hypot(d["y"] - cy, d["x"] - cx)
            if dist <= 2.5 and j not in used and (best is None or dist < best[0]):
                best = (dist, j)
        truth = timesim.frame_bits(ch["id"], ch["value"])
        bits_tot += truth.size
        if best is None:
            bit_err += truth.size // 2
            continue
        used.add(best[1]); d = dec[best[1]]
        bit_err += int(np.sum(d["bits"] != truth))
        ok += int(d["id"] == ch["id"] and d["value"] == ch["value"] and d["preamble_ok"])
    false = len(dec) - len(used)
    return dict(ok=ok, n=len(sc["chips"]), bit_err=bit_err, bits=bits_tot, false=false, T_s=T,
                n_c=n_c, seed=args["seed"])
