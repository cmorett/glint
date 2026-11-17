"""Round-2 stages: tilted / Lambertian SPOT launch (slab and sphere) and sphere Stage 1/2.

Incoming light at the chip is taken from Stage 1 plane-crossing tallies: the unscattered
downward beam, the scattered downward light (polar histogram) and the upward light
(polar histogram); azimuth is uniform (true for the slab, and on the sphere's axis).
For a mirror whose camera-facing normal is m = (sin t, 0, -cos t), a photon crossing a
horizontal plane with direction u delivers power to the tilted mirror in proportion to
f = max(0, -u.m) / |u_z| (flux through a tilted area vs a horizontal one).
"""
import numpy as np

import config_v2 as C
import mc
import mc_sphere as MS
import stages


def mirror_normal(tilt_deg):
    t = np.radians(tilt_deg)
    return np.array([np.sin(t), 0.0, -np.cos(t)])


def sample_incoming(F_ball, H_down, H_up, n, rng):
    """Directions of light arriving at the chip plane: (ux, uy, uz, is_ballistic)."""
    bins = H_down.size
    w = np.concatenate([[F_ball], H_down, H_up])
    w = np.maximum(w, 0)
    idx = rng.choice(w.size, size=n, p=w / w.sum())
    ball = idx == 0
    down = (idx >= 1) & (idx <= bins)
    b = np.where(down, idx - 1, idx - 1 - bins)
    theta = np.where(ball, 0.0, (b + rng.random(n)) * (np.pi / 2) / bins)
    phi = 2 * np.pi * rng.random(n)
    st = np.sin(theta)
    uz = np.where(down | ball, np.cos(theta), -np.cos(theta))
    return st * np.cos(phi), st * np.sin(phi), uz, ball


def launch(F_ball, H_down, H_up, n, rng, tilt_deg=0.0, mode="specular"):
    """Reflected launch directions and weights for a tilted mirror or Lambertian reflector.

    Returns (ux, uy, uz, w, launch_ballistic, mean_f). The irradiance on the tilted SPOT
    relative to the total crossing flux is mean_f; w is f / mean_f (mean 1).
    """
    ux, uy, uz, ball = sample_incoming(F_ball, H_down, H_up, n, rng)
    m = mirror_normal(tilt_deg)
    udotm = ux * m[0] + uy * m[1] + uz * m[2]
    f = np.maximum(0.0, -udotm) / np.maximum(np.abs(uz), 0.02)
    mean_f = float(f.mean())
    if mode == "specular":
        rx, ry, rz = ux - 2 * udotm * m[0], uy - 2 * udotm * m[1], uz - 2 * udotm * m[2]
        lb = ball
    elif mode == "lambert":
        ct = np.sqrt(rng.random(n)); st = np.sqrt(1 - ct * ct)
        ph = 2 * np.pi * rng.random(n)
        t = np.radians(tilt_deg)
        e1 = np.array([np.cos(t), 0.0, np.sin(t)]); e2 = np.array([0.0, 1.0, 0.0])
        rx = ct * m[0] + st * (np.cos(ph) * e1[0] + np.sin(ph) * e2[0])
        ry = ct * m[1] + st * (np.cos(ph) * e1[1] + np.sin(ph) * e2[1])
        rz = ct * m[2] + st * (np.cos(ph) * e1[2] + np.sin(ph) * e2[2])
        lb = np.zeros(n, bool)
    else:
        raise ValueError(mode)
    keep = f > 0
    w = np.where(keep, f / mean_f, 0.0) if mean_f > 0 else np.zeros(n)
    return rx, ry, rz, w, lb & keep, mean_f


def _pack_exits(e, depth_mm, mua_list, n_out, n_t, lb_key="lb", z_exit=None):
    """Common post-processing: weights per mu_a, NA invariant, back-projection to chip plane."""
    ze = np.zeros_like(e["x"]) if z_exit is None else z_exit
    t = (depth_mm - ze) / e["uz"]              # uz < 0 for light heading to the camera
    x_bp = e["x"] + t * e["ux"]
    y_bp = e["y"] + t * e["uy"]
    w = e["w0"][None, :] * e["m"][None, :] * np.exp(-np.asarray(mua_list)[:, None] * e["L"][None, :])
    inv = n_out * np.sqrt(np.clip(1.0 - e["uz"] ** 2, 0, 1))
    delta = e[lb_key].astype(bool) & (e["nscat"] == 0)
    return dict(x=e["x"].astype(np.float32), y=e["y"].astype(np.float32),
                x_bp=x_bp.astype(np.float32), y_bp=y_bp.astype(np.float32),
                inv=inv.astype(np.float32), w=w.astype(np.float32), delta=delta)


# ------------------------------------------------------------------ slab
def stage1_slab(slab, depths_mm, n, seed):
    return stages.stage1(slab, depths_mm, n, seed, mua_list=C.MUA_PER_MM, na_list=C.NA_LIST_V2)


def stage2_slab(slab, depth_mm, s1, k, n, seed, tilt_deg=0.0, mode="specular"):
    rng = np.random.default_rng(seed)
    ref = list(C.MUA_PER_MM).index(C.MUA_BASELINE)
    rx, ry, rz, w, lb, mean_f = launch(s1["F_ball"][ref, k], s1["theta_hist"][k],
                                       s1["theta_hist_up"][k], n, rng, tilt_deg, mode)
    start = dict(z=np.full(n, depth_mm), ux=rx, uy=ry, uz=rz, w0=w, launch_ballistic=lb)
    o = mc.run(slab, n, rng, mode="mbl", mua_list=C.MUA_PER_MM, start=start, track_xy=True,
               record_exits=True, roulette_threshold=C.ROULETTE_THRESHOLD,
               roulette_survival=C.ROULETTE_SURVIVAL)
    e = o["exits"]
    res = _pack_exits(e, depth_mm, C.MUA_PER_MM, slab.n, slab.n)
    # irradiance factor on the (tilted) mirror for each mu_a, relative to incident
    tot = s1["F_down"][:, k] + s1["F_up"][:, k]
    res.update(N=n, irr=tot * mean_f, depth_mm=depth_mm, tilt=tilt_deg, mode=mode, lost=o["lost"])
    return res


# ------------------------------------------------------------------ sphere
def tissue_sphere(musp_780, wavelength_nm, D_mm, mua=C.MUA_BASELINE):
    musp = C.musp_at(musp_780, wavelength_nm)
    return MS.Sphere(mus=musp / (1 - C.G), mua=mua, g=C.G, R=D_mm / 2, n=C.N_TISSUE, n_out=C.N_MEDIUM)


def stage1_sphere(sph, depths_mm, n, seed):
    rng = np.random.default_rng(seed)
    r_s = C.SPHERE_SCORING_RADIUS_MM
    start, r_sp_mean = MS.collimated_start(sph, n, rng, r_core=3 * r_s)
    o = MS.run(sph, start, rng, mua_list=C.MUA_PER_MM, planes=depths_mm, r_score=r_s,
               theta_bins=C.THETA_BINS, roulette_threshold=C.ROULETTE_THRESHOLD,
               roulette_survival=C.ROULETTE_SURVIVAL)
    scale = (sph.R / r_s) ** 2          # photons represent E * pi R^2 / N each
    e = o["exits"]
    w = e["w0"][None, :] * e["m"][None, :] * np.exp(-np.asarray(C.MUA_PER_MM)[:, None] * e["L"][None, :])
    inv = sph.n_out * np.sqrt(np.clip(1 - e["uz"] ** 2, 0, 1))
    keep = inv <= max(C.NA_LIST_V2)
    return dict(planes=np.asarray(depths_mm), F_down=o["F_down"] * scale, F_ball=o["F_ball"] * scale,
                F_up=o["F_up"] * scale, theta_hist=o["theta_hist"] * scale,
                theta_hist_up=o["theta_hist_up"] * scale,
                r_sp=r_sp_mean, A=o["A"], exit_total=o["exit_total"], lost=o["lost"], N=n,
                R=sph.R, bg=dict(x=e["x"][keep], y=e["y"][keep], z=e["z"][keep], ux=e["ux"][keep],
                                 uy=e["uy"][keep], uz=e["uz"][keep], inv=inv[keep],
                                 w=(w[:, keep] * np.pi * sph.R ** 2 / n).astype(np.float64)))


def stage2_sphere(sph, depth_mm, s1, k, n, seed, tilt_deg=0.0, mode="specular"):
    rng = np.random.default_rng(seed)
    ref = list(C.MUA_PER_MM).index(C.MUA_BASELINE)
    rx, ry, rz, w, lb, mean_f = launch(s1["F_ball"][ref, k], s1["theta_hist"][k],
                                       s1["theta_hist_up"][k], n, rng, tilt_deg, mode)
    start = dict(x=np.zeros(n), y=np.zeros(n), z=np.full(n, depth_mm), ux=rx, uy=ry, uz=rz,
                 w0=w, launch_ballistic=lb)
    o = MS.run(sph, start, rng, mua_list=C.MUA_PER_MM, roulette_threshold=C.ROULETTE_THRESHOLD,
               roulette_survival=C.ROULETTE_SURVIVAL)
    res = _pack_exits(o["exits"], depth_mm, C.MUA_PER_MM, sph.n_out, sph.n, z_exit=o["exits"]["z"])
    tot = s1["F_down"][:, k] + s1["F_up"][:, k]
    res.update(N=n, irr=tot * mean_f, depth_mm=depth_mm, tilt=tilt_deg, mode=mode, lost=o["lost"])
    return res
