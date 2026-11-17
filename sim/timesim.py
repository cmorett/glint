"""Synthetic camera frames of a region containing several blinking SPOT chips.

THIS IS SYNTHETIC DATA. Photon budgets come from the Monte Carlo results; the noise is a model.

Signal: each chip's reflectance follows a square-wave carrier at ~f_b (its own clock,
within +/- CLOCK_TOL of nominal, random phase), phase-modulated by DBPSK symbols of
N_c carrier cycles each. Symbols carry a 29-bit frame repeated cyclically:
Barker-13 preamble + 8-bit ID + 8-bit value. Frames integrate the waveform over the
exposure (exposure = frame period).

Noise: shot noise (Poisson), read noise (Gaussian, per (super)pixel), light-source
intensity noise (white frame-to-frame + slow Ornstein-Uhlenbeck drift, common to the
whole field), fixed-pattern noise (static gain PRNU and offset DSNU), static texture in
the tissue background, and slow sample motion (OU random walk, a few um).
"""
import numpy as np

import config_v2 as C


def frame_bits(chip_id, value):
    bits = list(C.PREAMBLE)
    bits += [(chip_id >> (C.ID_BITS - 1 - i)) & 1 for i in range(C.ID_BITS)]
    bits += [(value >> (C.VALUE_BITS - 1 - i)) & 1 for i in range(C.VALUE_BITS)]
    return np.array(bits, int)


def _ou(n, dt, tau, rms, rng):
    if rms == 0 or n == 0:
        return np.zeros(n)
    a = np.exp(-dt / tau)
    e = rng.standard_normal(n) * rms * np.sqrt(1 - a * a)
    out = np.empty(n)
    x = rms * rng.standard_normal()
    for i in range(n):          # n <= ~4e4, fine
        x = a * x + e[i]
        out[i] = x
    return out


def make_scene(chips, grid, sp_um, fps, f_nominal, n_c, e_bg_um2, read_noise_sp, rng,
               sigma_src=1e-3, motion_um=C.MOTION_UM, texture=C.BG_TEXTURE, prnu=C.PRNU,
               dsnu=C.DSNU_E):
    """chips: list of dicts with pos_um (x, y), pts_um (M x 2 image sample offsets, um),
    e_mean (electrons/frame at mean reflectance, whole image), e_mod (e_on - e_off), id, value."""
    ny, nx = grid
    sc = dict(grid=grid, sp_um=sp_um, fps=fps, f_nominal=f_nominal, n_c=n_c, sigma_src=sigma_src,
              motion_um=motion_um, read_noise_sp=read_noise_sp)
    tex = 1 + texture * _smooth_field(ny, nx, 2.0, rng)
    sc["bg"] = e_bg_um2 * sp_um * sp_um * tex
    sc["gain"] = 1 + prnu * rng.standard_normal((ny, nx))
    sc["offset"] = dsnu * rng.standard_normal((ny, nx))
    smax = int(np.ceil(4 * motion_um)) if motion_um > 0 else 0
    sc["shifts"] = np.arange(-smax, smax + 1)
    ex = np.arange(nx + 1) * sp_um; ey = np.arange(ny + 1) * sp_um
    out = []
    for c in chips:
        imgs = np.zeros((len(sc["shifts"]), len(sc["shifts"]), ny, nx), np.float32)
        M = len(c["pts_um"])
        for ix, sx in enumerate(sc["shifts"]):
            for iy, sy in enumerate(sc["shifts"]):
                h, _, _ = np.histogram2d(c["pts_um"][:, 1] + c["pos_um"][1] + sy,
                                         c["pts_um"][:, 0] + c["pos_um"][0] + sx, bins=[ey, ex])
                imgs[ix, iy] = h / M
        cc = dict(c, imgs=imgs,
                  f=f_nominal * (1 + rng.uniform(-C.CLOCK_TOL, C.CLOCK_TOL)),
                  phase=rng.uniform(0, 2 * np.pi), t0=rng.uniform(0, 29 * n_c / f_nominal),
                  bits=frame_bits(c["id"], c["value"]))
        out.append(cc)
    sc["chips"] = out
    return sc


def _smooth_field(ny, nx, corr_px, rng):
    f = rng.standard_normal((ny, nx))
    k = np.exp(-0.5 * (np.arange(-6, 7) / corr_px) ** 2); k /= k.sum()
    f = np.apply_along_axis(lambda v: np.convolve(v, k, "same"), 0, f)
    f = np.apply_along_axis(lambda v: np.convolve(v, k, "same"), 1, f)
    return f / (f.std() + 1e-12)


def carrier(chip, t_start, dt, n_sub=8):
    """Exposure-averaged DBPSK square-wave value in [-1, 1] for frames starting at t_start."""
    ts = t_start[:, None] + (np.arange(n_sub) + 0.5) / n_sub * dt
    sq = np.sign(np.sin(2 * np.pi * chip["f"] * ts + chip["phase"]))
    sym = np.floor((ts + chip["t0"]) * chip["f"] / chip["n_c"]).astype(np.int64)
    nb = len(chip["bits"])
    # DBPSK: phase flips at each '1'; cumulative parity over the cyclic bit stream
    bits = chip["bits"]
    cyc_par = bits.sum() % 2
    cum = np.concatenate([[0], np.cumsum(bits)]) % 2
    par = ((sym // nb) * cyc_par + cum[sym % nb]) % 2
    return (sq * np.where(par == 1, -1.0, 1.0)).mean(axis=1)


def render(sc, T_s, rng, ref_noise=True):
    """All frames for a read of duration T_s. Returns float32 array (n_frames, ny, nx)."""
    fps = sc["fps"]; dt = 1.0 / fps
    n = int(round(T_s * fps))
    t = np.arange(n) * dt
    for ch in sc["chips"]:
        ch.setdefault("n_c", sc["n_c"])
    s = sc["sigma_src"] / np.sqrt(2)
    eps = s * rng.standard_normal(n) + _ou(n, dt, C.DRIFT_TAU_S, s, rng)
    mx = _ou(n, dt, C.MOTION_TAU_S, sc["motion_um"], rng)
    my = _ou(n, dt, C.MOTION_TAU_S, sc["motion_um"], rng)
    sh = sc["shifts"]
    ix = np.clip(np.round(mx).astype(int) - sh[0], 0, len(sh) - 1)
    iy = np.clip(np.round(my).astype(int) - sh[0], 0, len(sh) - 1)
    out = np.empty((n,) + sc["grid"], np.float32)
    for a in range(0, n, 2000):
        b = min(n, a + 2000)
        mean = np.broadcast_to(sc["bg"], (b - a,) + sc["grid"]).copy()
        for ch in sc["chips"]:
            e = ch["e_mean"] + 0.5 * ch["e_mod"] * carrier(ch, t[a:b], dt)
            mean += ch["imgs"][ix[a:b], iy[a:b]] * e[:, None, None]
        mean *= sc["gain"][None] * (1 + eps[a:b])[:, None, None]
        mean = np.maximum(mean, 0)
        frm = rng.poisson(mean).astype(np.float32)
        frm += (sc["read_noise_sp"] * rng.standard_normal(frm.shape)).astype(np.float32)
        frm += sc["offset"][None].astype(np.float32)
        out[a:b] = frm
    return out
