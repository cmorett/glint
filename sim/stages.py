"""Stage 1 (wide-field illumination, 1-D) and Stage 2 (light leaving the SPOT) runs.

Stage 1: a collimated beam enters the slab at normal incidence. Because the
  illumination is uniform and the slab is laterally infinite, the problem is 1-D:
  the downward irradiance and its angular distribution at depth d are the same
  everywhere in the plane. We score at every chip depth:
    F_down(d)  downward irradiance / incident irradiance (what lands on the mirror)
    F_ball(d)  the never-scattered part of it
    theta histogram of the scattered downward light
  plus the diffuse reflectance leaving the camera-side surface inside each lens NA
  (this is the background).

Stage 2: photons start at a point on the mirror plane z = d, with directions drawn
  from the Stage 1 downward distribution and reflected specularly (uz -> -uz). We
  track them to the surface and keep exit position, direction and path length.
  The square SPOT is handled later by convolving with its footprint.
"""
import numpy as np

import config as C
import mc


def tissue_slab(musp_780, wavelength_nm, n_cam, n_far, mua=C.MUA_BASELINE,
                thickness_mm=C.SLAB_THICKNESS_MM):
    musp = C.musp_at(musp_780, wavelength_nm)
    return mc.Slab(mus=musp / (1.0 - C.G), mua=mua, g=C.G, T=thickness_mm,
                   n=C.N_TISSUE, n_top=n_cam, n_bot=n_far)


def stage1(slab, depths_mm, n_photons, seed, mua_list=C.MUA_PER_MM, na_list=C.NA_LIST):
    rng = np.random.default_rng(seed)
    planes = [d for d in depths_mm if d > 0]
    o = mc.run(slab, n_photons, rng, mode="mbl", mua_list=mua_list,
               planes=planes if planes else None, theta_bins=C.THETA_BINS,
               roulette_threshold=C.ROULETTE_THRESHOLD,
               roulette_survival=C.ROULETTE_SURVIVAL)
    inv = mc.exit_invariant(slab, o["top_cos"])
    Rd_na = np.array([[o["top_w"][j, inv <= na].sum() / n_photons for na in na_list]
                      for j in range(len(mua_list))])
    res = dict(r_sp=o["r_sp"], Rd=o["R_top"], Rd_na=Rd_na, T=o["T_bot"], A=o["A"],
               lost=o["lost"], mua_list=np.asarray(mua_list), na_list=np.asarray(na_list))
    if planes:
        res.update(planes=np.asarray(planes), F_down=o["F_down"], F_ball=o["F_ball"],
                   F_up=o["F_up"], theta_hist=o["theta_hist"],
                   theta_hist_up=o["theta_hist_up"])
    return res


def sample_launch(F_ball, F_scat_hist, n, rng):
    """Directions of light arriving at the mirror (from Stage 1), reflected upward."""
    tot = F_ball + F_scat_hist.sum()
    p_ball = F_ball / tot
    nb = rng.binomial(n, p_ball)
    bins = len(F_scat_hist)
    pdf = F_scat_hist / F_scat_hist.sum() if F_scat_hist.sum() > 0 else np.ones(bins) / bins
    idx = rng.choice(bins, size=n - nb, p=pdf)
    theta = (idx + rng.random(n - nb)) * (np.pi / 2) / bins
    theta = np.concatenate([np.zeros(nb), theta])
    phi = 2 * np.pi * rng.random(n)
    st = np.sin(theta)
    # incoming direction (ux, uy, +cos); specular mirror in the z-plane flips uz
    return st * np.cos(phi), st * np.sin(phi), -np.cos(theta)


def stage2(slab, depth_mm, s1, plane_index, n_photons, seed, mua_list=C.MUA_PER_MM):
    rng = np.random.default_rng(seed)
    ref = int(np.argmin(np.abs(np.asarray(mua_list) - slab.mua)))
    F_ball = s1["F_ball"][ref, plane_index]
    hist = s1["theta_hist"][plane_index]
    ux, uy, uz = sample_launch(F_ball, hist, n_photons, rng)
    start = dict(z=np.full(n_photons, depth_mm), ux=ux, uy=uy, uz=uz)
    o = mc.run(slab, n_photons, rng, mode="mbl", mua_list=mua_list, start=start,
               track_xy=True, record_exits=True,
               roulette_threshold=C.ROULETTE_THRESHOLD,
               roulette_survival=C.ROULETTE_SURVIVAL)
    e = o["exits"]
    x_bp = e["x"] + depth_mm * e["ux"] / e["uz"]     # ray traced back to the chip plane
    y_bp = e["y"] + depth_mm * e["uy"] / e["uz"]
    w = e["w0"][None, :] * e["m"][None, :] * np.exp(-np.asarray(mua_list)[:, None] * e["L"][None, :])
    inv = mc.exit_invariant(slab, np.abs(e["uz"]))
    return dict(N=n_photons, x=e["x"].astype(np.float32), y=e["y"].astype(np.float32),
                x_bp=x_bp.astype(np.float32), y_bp=y_bp.astype(np.float32),
                inv=inv.astype(np.float32), w=w.astype(np.float32),
                ballistic=(e["nscat"] == 0), R_top=o["R_top"], T_bot=o["T_bot"],
                lost=o["lost"], depth_mm=depth_mm)
