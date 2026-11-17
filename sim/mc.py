"""Vectorized Monte Carlo photon transport in a single turbid slab (NumPy).

Coordinates: z points INTO the tissue; the slab occupies 0 <= z <= T.
z = 0 is the camera-side surface (outer index n_top), z = T the far side (n_bot).

Two absorption schemes, both standard and both validated in tests/:
  mode="mcml": MCML (Wang, Jacques & Zheng 1995). Step length ~ Exp(mu_t); at each
               interaction the weight drops by mu_a/mu_t (implicit capture).
  mode="mbl" : "microscopic Beer-Lambert" (absorption weighting along the path).
               Step length ~ Exp(mu_s); the weight is exp(-mu_a * L) with L the
               total path length. Because L is tracked, any list of mu_a values
               is scored from one run (exact, not an approximation).
Both use Henyey-Greenstein scattering, all-or-none Fresnel reflection/refraction
at the two boundaries, and Russian roulette.
"""
import numpy as np


def fresnel(n1, n2, cos_i):
    """Unpolarized Fresnel reflectance going from n1 into n2. cos_i in [0, 1]."""
    cos_i = np.clip(cos_i, 0.0, 1.0)
    if n1 == n2:
        return np.zeros_like(cos_i)
    sin_i = np.sqrt(1.0 - cos_i * cos_i)
    sin_t = n1 / n2 * sin_i
    tir = sin_t >= 1.0
    cos_t = np.sqrt(np.clip(1.0 - sin_t * sin_t, 0.0, 1.0))
    with np.errstate(invalid="ignore", divide="ignore"):
        rs = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)
        rp = (n1 * cos_t - n2 * cos_i) / (n1 * cos_t + n2 * cos_i)
        r = 0.5 * (rs * rs + rp * rp)
    r = np.where(tir | ~np.isfinite(r), 1.0, r)
    return r


def hg_cos(g, xi):
    if g == 0.0:
        return 2.0 * xi - 1.0
    tmp = (1.0 - g * g) / (1.0 - g + 2.0 * g * xi)
    return np.clip((1.0 + g * g - tmp * tmp) / (2.0 * g), -1.0, 1.0)


class Slab:
    def __init__(self, mus, mua, g, T, n, n_top, n_bot):
        self.mus, self.mua, self.g, self.T = float(mus), float(mua), float(g), float(T)
        self.n, self.n_top, self.n_bot = float(n), float(n_top), float(n_bot)


def run(slab, n_photons, rng, *, mode="mbl", mua_list=None, start=None,
        planes=None, theta_bins=180, track_xy=False, record_exits=False,
        roulette_threshold=1e-4, roulette_survival=0.1, max_iter=200_000):
    """Propagate photons and return tallies.

    start=None -> collimated beam normally incident on z=0 (specular loss r_sp removed
                  deterministically, as in MCML). Otherwise a dict with arrays
                  z, ux, uy, uz (and optional w0) giving launch points inside the slab.
    planes     -> depths (mm) at which DOWNWARD crossings are scored (mbl only):
                  total weight, unscattered weight, and polar-angle histogram.
    mua_list   -> (mbl) absorption coefficients scored simultaneously. Scalar tallies
                  are returned for every entry; the angular histogram uses the entry
                  closest to slab.mua.
    record_exits -> keep per-photon records of everything leaving through z = 0.
    """
    s = slab
    if mua_list is None:
        mua_list = [s.mua]
    mua_arr = np.asarray(mua_list, float)
    nm = len(mua_arr)
    ref = int(np.argmin(np.abs(mua_arr - s.mua)))
    mua_rr = mua_arr.min()          # roulette on the least-absorbing weight (safe for all)

    out = {}
    N = int(n_photons)
    if start is None:
        r_sp = float(fresnel(s.n_top, s.n, np.array([1.0]))[0])
        z = np.zeros(N); ux = np.zeros(N); uy = np.zeros(N); uz = np.ones(N)
        w0 = np.full(N, 1.0 - r_sp)
        out["r_sp"] = r_sp
    else:
        z = np.array(start["z"], float); ux = np.array(start["ux"], float)
        uy = np.array(start["uy"], float); uz = np.array(start["uz"], float)
        w0 = np.array(start.get("w0", np.ones_like(z)), float)
        out["r_sp"] = 0.0
    x = np.zeros(N); y = np.zeros(N)
    L = np.zeros(N)                 # path length (mm)
    m = np.ones(N)                  # roulette multiplier (mbl) or weight (mcml)
    if mode == "mcml":
        m = w0.copy()
    nscat = np.zeros(N, np.int32)

    R_top = np.zeros(nm); T_bot = np.zeros(nm); A = 0.0
    T_unscat = np.zeros(nm)
    lost = 0.0
    if planes is not None:
        planes = np.asarray(planes, float)
        P = len(planes)
        F = np.zeros((nm, P)); Fb = np.zeros((nm, P))
        H = np.zeros((P, theta_bins))
        F_up = np.zeros((nm, P))
    exits_top = {k: [] for k in ("x", "y", "ux", "uy", "uz", "L", "m", "w0", "nscat")}
    top_cos_all = []; top_w_all = []

    mut_step = s.mus if mode == "mbl" else s.mus + s.mua

    def weights(idx_m, idx_L, idx_w0):
        """Weights for every mu_a in the list (mbl), shape (nm, n)."""
        if mode == "mcml":
            return idx_m[None, :] * np.ones((nm, 1))
        return idx_w0[None, :] * idx_m[None, :] * np.exp(-mua_arr[:, None] * idx_L[None, :])

    it = 0
    while z.size and it < max_iter:
        it += 1
        n = z.size
        if mut_step > 0:
            step = -np.log(rng.random(n)) / mut_step
        else:
            step = np.full(n, np.inf)
        with np.errstate(divide="ignore", invalid="ignore"):
            db = np.where(uz > 0, (s.T - z) / uz, np.where(uz < 0, -z / uz, np.inf))
        hit = step >= db
        seg = np.where(hit, db, step)
        z_end = np.where(hit, np.where(uz > 0, s.T, 0.0), z + seg * uz)

        # ---- plane-crossing tallies (downward and upward irradiance)
        if planes is not None:
            for k, p in enumerate(planes):
                c = (uz > 0) & (z < p) & (z_end >= p)
                if c.any():
                    Lc = L[c] + (p - z[c]) / uz[c]
                    W = weights(m[c], Lc, w0[c])
                    F[:, k] += W.sum(1)
                    b = nscat[c] == 0
                    if b.any():
                        Fb[:, k] += W[:, b].sum(1)
                    th = np.arccos(np.clip(uz[c], -1, 1))
                    bins = np.minimum((th / (np.pi / 2) * theta_bins).astype(int), theta_bins - 1)
                    # histogram only scattered light; unscattered is a delta at theta = 0
                    H[k] += np.bincount(bins[~b], weights=W[ref, ~b], minlength=theta_bins)
                cu = (uz < 0) & (z > p) & (z_end <= p)
                if cu.any():
                    Lc = L[cu] + (z[cu] - p) / (-uz[cu])
                    F_up[:, k] += weights(m[cu], Lc, w0[cu]).sum(1)

        # ---- absorption bookkeeping (for the energy-conservation check)
        if mode == "mbl":
            wb = w0 * m * np.exp(-s.mua * L)
            A += float(np.sum(wb * (1.0 - np.exp(-s.mua * seg))))
        if track_xy:
            x = x + seg * ux
            y = y + seg * uy
        z = z_end
        L = L + seg
        alive = np.ones(n, bool)

        # ---- boundaries
        hi = np.nonzero(hit)[0]
        if hi.size:
            up = uz[hi] < 0
            cos_i = np.abs(uz[hi])
            r = np.where(up, fresnel(s.n, s.n_top, cos_i), fresnel(s.n, s.n_bot, cos_i))
            refl = rng.random(hi.size) < r
            uz[hi[refl]] = -uz[hi[refl]]
            esc = hi[~refl]
            esc_up = esc[uz[esc] < 0]
            esc_dn = esc[uz[esc] > 0]
            if esc_up.size:
                W = weights(m[esc_up], L[esc_up], w0[esc_up])
                R_top += W.sum(1)
                top_cos_all.append(np.abs(uz[esc_up]))
                top_w_all.append(W)
                if record_exits:
                    for key, arr in (("x", x), ("y", y), ("ux", ux), ("uy", uy), ("uz", uz),
                                     ("L", L), ("m", m), ("w0", w0), ("nscat", nscat)):
                        exits_top[key].append(arr[esc_up].copy())
            if esc_dn.size:
                W = weights(m[esc_dn], L[esc_dn], w0[esc_dn])
                T_bot += W.sum(1)
                u0 = nscat[esc_dn] == 0
                if u0.any():
                    T_unscat += W[:, u0].sum(1)
            alive[esc] = False

        # ---- scattering events
        si = np.nonzero(~hit)[0]
        if si.size:
            if mode == "mcml":
                dep = m[si] * (s.mua / mut_step)
                A += float(dep.sum())
                m[si] -= dep
            ct = hg_cos(s.g, rng.random(si.size))
            st = np.sqrt(1.0 - ct * ct)
            phi = 2.0 * np.pi * rng.random(si.size)
            cp, sp = np.cos(phi), np.sin(phi)
            uzs = uz[si]
            if track_xy:
                uxs, uys = ux[si], uy[si]
                big = np.abs(uzs) > 0.99999
                tmp = np.sqrt(np.maximum(1.0 - uzs * uzs, 1e-300))
                nux = np.where(big, st * cp, st * (uxs * uzs * cp - uys * sp) / tmp + uxs * ct)
                nuy = np.where(big, st * sp, st * (uys * uzs * cp + uxs * sp) / tmp + uys * ct)
                nuz = np.where(big, ct * np.sign(uzs), -st * cp * tmp + uzs * ct)
                norm = np.sqrt(nux * nux + nuy * nuy + nuz * nuz)
                ux[si], uy[si], uz[si] = nux / norm, nuy / norm, nuz / norm
            else:
                # 1-D problem: only the polar cosine matters (lateral symmetry)
                uz[si] = np.clip(uzs * ct + np.sqrt(np.maximum(1.0 - uzs * uzs, 0.0)) * st * cp, -1, 1)
            nscat[si] += 1

        # ---- Russian roulette
        wref = m * w0 * np.exp(-mua_rr * L) if mode == "mbl" else m
        low = alive & (wref < roulette_threshold)
        if low.any():
            li = np.nonzero(low)[0]
            surv = rng.random(li.size) < roulette_survival
            m[li[surv]] /= roulette_survival
            alive[li[~surv]] = False

        keep = alive
        if not keep.all():
            z, ux, uy, uz, x, y, L, m, w0, nscat = (a[keep] for a in (z, ux, uy, uz, x, y, L, m, w0, nscat))

    if z.size:   # hit iteration cap: report the weight we did not finish
        lost = float(np.sum(w0 * m * np.exp(-s.mua * L)) if mode == "mbl" else np.sum(m))

    out.update(N=N, R_top=R_top / N, T_bot=T_bot / N, T_unscat=T_unscat / N,
               A=A / N, lost=lost / N, iterations=it, mua_list=mua_arr)
    if top_cos_all:
        out["top_cos"] = np.concatenate(top_cos_all)
        out["top_w"] = np.concatenate(top_w_all, axis=1)
    else:
        out["top_cos"] = np.zeros(0); out["top_w"] = np.zeros((nm, 0))
    if planes is not None:
        out.update(planes=planes, F_down=F / N, F_ball=Fb / N, F_up=F_up / N,
                   theta_hist=H / N, theta_bins=theta_bins)
    if record_exits:
        out["exits"] = {k: (np.concatenate(v) if v else np.zeros(0)) for k, v in exits_top.items()}
    return out


def exit_invariant(slab, cos_tissue):
    """n*sin(theta), conserved across parallel interfaces; compare with lens NA (in air)."""
    return slab.n * np.sqrt(np.clip(1.0 - cos_tissue ** 2, 0.0, 1.0))
