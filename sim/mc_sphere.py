"""Monte Carlo transport in a tissue SPHERE surrounded by culture medium (NumPy, vectorized).

Same physics as mc.py (HG scattering, absorption weighting w = exp(-mu_a L), all-or-none
Fresnel, Russian roulette); only the boundary is a sphere instead of two planes.

Coordinates: z points away from the camera. The sphere (radius R) has its centre at
(0, 0, R), so its camera-side apex is at z = 0 and its far apex at z = 2R. Light leaving
the sphere travels straight through the non-scattering medium; we assume it does not
re-enter the sphere (true for a convex body) and ignore light reflected back by the dish.

Collimated illumination from the camera side covers the sphere's whole cross-section.
Launch points are importance-sampled: half uniform over the disk of radius R, half
uniform over a small disk of radius r_c around the axis (where the chips sit); each
photon carries the weight that undoes this bias, so all tallies remain unbiased.
"""
import numpy as np

from mc import fresnel, hg_cos


def _scatter3d(ux, uy, uz, g, rng):
    n = ux.size
    ct = hg_cos(g, rng.random(n))
    st = np.sqrt(1.0 - ct * ct)
    phi = 2.0 * np.pi * rng.random(n)
    cp, sp = np.cos(phi), np.sin(phi)
    big = np.abs(uz) > 0.99999
    tmp = np.sqrt(np.maximum(1.0 - uz * uz, 1e-300))
    nux = np.where(big, st * cp, st * (ux * uz * cp - uy * sp) / tmp + ux * ct)
    nuy = np.where(big, st * sp, st * (uy * uz * cp + ux * sp) / tmp + uy * ct)
    nuz = np.where(big, ct * np.sign(uz), -st * cp * tmp + uz * ct)
    norm = np.sqrt(nux * nux + nuy * nuy + nuz * nuz)
    return nux / norm, nuy / norm, nuz / norm


def refract(u, nrm, cos_i, eta):
    """Refracted direction. nrm: unit normal on the INCIDENT side; cos_i = -u.nrm > 0."""
    sin_t2 = eta * eta * (1.0 - cos_i * cos_i)
    cos_t = np.sqrt(np.clip(1.0 - sin_t2, 0.0, 1.0))
    return eta * u + (eta * cos_i - cos_t)[:, None] * nrm


class Sphere:
    def __init__(self, mus, mua, g, R, n, n_out):
        self.mus, self.mua, self.g, self.R = float(mus), float(mua), float(g), float(R)
        self.n, self.n_out = float(n), float(n_out)


def collimated_start(sph, N, rng, r_core):
    """Launch points on the camera-side hemisphere for a beam travelling along +z."""
    R = sph.R
    r_core = min(r_core, R)
    half = rng.random(N) < 0.5
    rad = np.where(half, R * np.sqrt(rng.random(N)), r_core * np.sqrt(rng.random(N)))
    ang = 2 * np.pi * rng.random(N)
    x, y = rad * np.cos(ang), rad * np.sin(ang)
    # mixture pdf per unit area, relative to uniform over the full disk
    q = 0.5 / (np.pi * R * R) + 0.5 * (rad < r_core) / (np.pi * r_core * r_core)
    w_imp = (1.0 / (np.pi * R * R)) / q
    z = R - np.sqrt(np.maximum(R * R - rad * rad, 0.0))
    nrm = np.stack([x, y, z - R], axis=1) / R          # outward normal (incident side)
    u = np.tile([0.0, 0.0, 1.0], (N, 1))
    cos_i = -nrm[:, 2]
    r = fresnel(sph.n_out, sph.n, cos_i)
    ut = refract(u, nrm, cos_i, sph.n_out / sph.n)
    # nudge inside to avoid starting exactly on the boundary
    pos = np.stack([x, y, z], axis=1) + 1e-9 * ut
    return dict(x=pos[:, 0], y=pos[:, 1], z=pos[:, 2], ux=ut[:, 0], uy=ut[:, 1], uz=ut[:, 2],
                w0=w_imp * (1.0 - r)), float(np.mean(w_imp * r))


def run(sph, start, rng, *, mua_list, planes=None, r_score=0.05, theta_bins=180,
        roulette_threshold=1e-4, roulette_survival=0.1, max_iter=200_000, keep_exit_up=True):
    mua_arr = np.asarray(mua_list, float)
    nm = len(mua_arr)
    ref = int(np.argmin(np.abs(mua_arr - sph.mua)))
    mua_rr = mua_arr.min()
    x, y, z = (np.array(start[k], float) for k in ("x", "y", "z"))
    ux, uy, uz = (np.array(start[k], float) for k in ("ux", "uy", "uz"))
    N = x.size
    w0 = np.array(start.get("w0", np.ones(N)), float)
    lb = np.array(start.get("launch_ballistic", np.zeros(N, bool)))
    L = np.zeros(N); m = np.ones(N); nscat = np.zeros(N, np.int32)
    R = sph.R
    A = 0.0
    exit_tot = np.zeros(nm)
    if planes is not None:
        planes = np.asarray(planes, float); P = planes.size
        F = np.zeros((nm, P)); Fb = np.zeros((nm, P)); Fu = np.zeros((nm, P))
        H = np.zeros((P, theta_bins)); Hu = np.zeros((P, theta_bins))
    ex = {k: [] for k in ("x", "y", "z", "ux", "uy", "uz", "L", "m", "w0", "nscat", "lb")}

    def W(idx_m, idx_L, idx_w0):
        return idx_w0[None, :] * idx_m[None, :] * np.exp(-mua_arr[:, None] * idx_L[None, :])

    it = 0
    while x.size and it < max_iter:
        it += 1
        n = x.size
        step = -np.log(rng.random(n)) / sph.mus if sph.mus > 0 else np.full(n, np.inf)
        px, py, pz = x, y, z - R
        b = px * ux + py * uy + pz * uz
        cc = px * px + py * py + pz * pz - R * R
        tb = -b + np.sqrt(np.maximum(b * b - cc, 0.0))
        hit = step >= tb
        seg = np.where(hit, tb, step)
        xe, ye, ze = x + seg * ux, y + seg * uy, z + seg * uz
        if planes is not None:
            for k, p in enumerate(planes):
                for down in (True, False):
                    c = ((uz > 0) & (z < p) & (ze >= p)) if down else ((uz < 0) & (z > p) & (ze <= p))
                    if not c.any():
                        continue
                    t = (p - z[c]) / uz[c]
                    cx, cy = x[c] + t * ux[c], y[c] + t * uy[c]
                    ins = cx * cx + cy * cy <= r_score * r_score
                    if not ins.any():
                        continue
                    ci = np.nonzero(c)[0][ins]
                    Wc = W(m[ci], L[ci] + t[ins], w0[ci])
                    cosz = np.abs(uz[ci])
                    bins = np.minimum((np.arccos(np.clip(cosz, -1, 1)) / (np.pi / 2) * theta_bins).astype(int),
                                      theta_bins - 1)
                    if down:
                        F[:, k] += Wc.sum(1)
                        bal = (nscat[ci] == 0)
                        Fb[:, k] += Wc[:, bal].sum(1)
                        H[k] += np.bincount(bins[~bal], weights=Wc[ref, ~bal], minlength=theta_bins)
                    else:
                        Fu[:, k] += Wc.sum(1)
                        Hu[k] += np.bincount(bins, weights=Wc[ref], minlength=theta_bins)
        wb = w0 * m * np.exp(-sph.mua * L)
        A += float(np.sum(wb * (1.0 - np.exp(-sph.mua * seg))))
        x, y, z = xe, ye, ze
        L = L + seg
        alive = np.ones(n, bool)

        hi = np.nonzero(hit)[0]
        if hi.size:
            nrm = np.stack([x[hi], y[hi], z[hi] - R], axis=1) / R
            u = np.stack([ux[hi], uy[hi], uz[hi]], axis=1)
            cos_i = np.clip(np.sum(u * nrm, axis=1), 0.0, 1.0)
            r = fresnel(sph.n, sph.n_out, cos_i)
            refl = rng.random(hi.size) < r
            ur = u - 2 * cos_i[:, None] * nrm
            # reflected: step back inside slightly
            ri = hi[refl]
            ux[ri], uy[ri], uz[ri] = ur[refl, 0], ur[refl, 1], ur[refl, 2]
            x[ri] -= 1e-9 * nrm[refl, 0]; y[ri] -= 1e-9 * nrm[refl, 1]; z[ri] -= 1e-9 * nrm[refl, 2]
            esc = ~refl
            if esc.any():
                ei = hi[esc]
                ut = refract(u[esc], -nrm[esc], cos_i[esc], sph.n / sph.n_out)
                Wx = W(m[ei], L[ei], w0[ei])
                exit_tot += Wx.sum(1)
                sel = ut[:, 2] < 0 if keep_exit_up else np.ones(ei.size, bool)
                if sel.any():
                    es = ei[sel]
                    for key, arr in (("x", x), ("y", y), ("z", z), ("L", L), ("m", m), ("w0", w0),
                                     ("nscat", nscat), ("lb", lb)):
                        ex[key].append(arr[es].copy())
                    ex["ux"].append(ut[sel, 0]); ex["uy"].append(ut[sel, 1]); ex["uz"].append(ut[sel, 2])
                alive[ei] = False

        si = np.nonzero(~hit)[0]
        if si.size:
            ux[si], uy[si], uz[si] = _scatter3d(ux[si], uy[si], uz[si], sph.g, rng)
            nscat[si] += 1

        wref = m * w0 * np.exp(-mua_rr * L)
        low = alive & (wref < roulette_threshold * np.maximum(w0, 1e-300))
        if low.any():
            li = np.nonzero(low)[0]
            surv = rng.random(li.size) < roulette_survival
            m[li[surv]] /= roulette_survival
            alive[li[~surv]] = False
        if not alive.all():
            x, y, z, ux, uy, uz, L, m, w0, nscat, lb = (a[alive] for a in (x, y, z, ux, uy, uz, L, m, w0, nscat, lb))

    lost = float(np.sum(w0 * m * np.exp(-sph.mua * L))) if x.size else 0.0
    out = dict(N=N, A=A / N, exit_total=exit_tot / N, lost=lost / N, iterations=it, mua_list=mua_arr)
    out["exits"] = {k: (np.concatenate(v) if v else np.zeros(0)) for k, v in ex.items()}
    if planes is not None:
        out.update(planes=planes, F_down=F / N, F_ball=Fb / N, F_up=Fu / N, theta_hist=H / N,
                   theta_hist_up=Hu / N)
    return out
