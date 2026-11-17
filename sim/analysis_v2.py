"""Round-2 analysis: any NA (incl. whole-plate NA 0.01-0.05), any camera, tilted chips, spheres.

Low-NA collection (NA < NA_EFF = 0.05): light that was never scattered on either leg and
reflected specularly forms an exact delta in angle, so it is counted directly. Everything
else is counted inside NA_EFF and scaled by (NA/NA_EFF)^2, i.e. its angular density is
taken as constant inside a cone of NA 0.05. tests/test_v2.py checks this against direct
counting; for shallow chips it slightly UNDER-estimates signal (forward-peaked light).
"""
import numpy as np

import config_v2 as C

H = 6.62607015e-34
CLIGHT = 2.99792458e8


def electrons_per_um2_per_frame(E_mW_cm2, lam_nm, fps):
    flux = E_mW_cm2 * 10.0 / (H * CLIGHT / (lam_nm * 1e-9)) / 1e12      # photons/um^2/s
    return C.QE * flux / fps


def diffraction_sigma_um(lam_nm, NA):
    return 0.21 * lam_nm * 1e-3 / NA      # Gaussian approximation of the Airy core


def collected(X, NA, j):
    """(weights, x_bp, y_bp) of light collected inside NA, with low-NA scaling."""
    w = X["w"][j].astype(float)
    if NA >= C.NA_EFF_SCALING:
        sel = X["inv"] <= NA
        return w[sel], X["x_bp"][sel].astype(float), X["y_bp"][sel].astype(float)
    d = X["delta"] & (X["inv"] <= NA + 1e-12)
    nd = (~X["delta"]) & (X["inv"] <= C.NA_EFF_SCALING)
    sc = (NA / C.NA_EFF_SCALING) ** 2
    ww = np.concatenate([w[d], w[nd] * sc])
    return ww, np.concatenate([X["x_bp"][d], X["x_bp"][nd]]).astype(float), \
        np.concatenate([X["y_bp"][d], X["y_bp"][nd]]).astype(float)


def slab_bg(s1, NA, j):
    """Tissue background (fraction of incident per unit area) inside NA for the slab."""
    nas = list(C.NA_LIST_V2)
    if NA in nas and NA >= C.NA_EFF_SCALING:
        return s1["r_sp"] + s1["Rd_na"][j, nas.index(NA)]
    ref = nas.index(C.NA_EFF_SCALING)
    return s1["r_sp"] + s1["Rd_na"][j, ref] * (NA / C.NA_EFF_SCALING) ** 2


def sphere_bg(s1, depth_mm, NA, j):
    """Background density near the axis of a sphere image, in the chip plane."""
    b = s1["bg"]
    t = (depth_mm - b["z"]) / b["uz"]
    xb, yb = b["x"] + t * b["ux"], b["y"] + t * b["uy"]
    rho_min = max(0.25 * s1["R"], 0.05)
    near = xb * xb + yb * yb <= rho_min ** 2
    na_eff = max(NA, C.NA_EFF_SCALING)
    sel = near & (b["inv"] <= na_eff)
    dens = b["w"][j][sel].sum() / (np.pi * rho_min ** 2)
    if NA < na_eff:
        dens *= (NA / na_eff) ** 2
    return dens + C.fresnel_normal(C.N_MEDIUM, C.N_TISSUE)


class Image:
    """Chip-plane image of a uniformly lit a x (a cos tilt) SPOT, plus optical blur."""
    RHO_UM = np.geomspace(1.0, 4000.0, 160)

    def __init__(self, w, xb_um, yb_um, a_um, tilt_deg, sigma_um, rng, min_samples=150_000):
        n = max(w.size, 1)
        if w.size == 0:
            w, xb_um, yb_um = np.zeros(1), np.zeros(1), np.zeros(1)
        k = max(1, int(np.ceil(min_samples / n)))
        b = a_um * abs(np.cos(np.radians(min(tilt_deg, 180 - tilt_deg) if tilt_deg <= 180 else 0)))
        N = w.size * k
        xs = np.repeat(xb_um, k) + (rng.random(N) - 0.5) * a_um + sigma_um * rng.standard_normal(N)
        ys = np.repeat(yb_um, k) + (rng.random(N) - 0.5) * b + sigma_um * rng.standard_normal(N)
        r = np.hypot(xs, ys)
        o = np.argsort(r)
        self.r = r[o]
        self.w = np.repeat(w, k)[o] / k
        self.tot = self.w.sum()

    def EE(self):
        if self.tot <= 0:
            return np.zeros_like(self.RHO_UM)
        cw = np.cumsum(self.w) / self.tot
        return np.interp(self.RHO_UM, self.r, cw, left=0.0, right=1.0)

    def quantile(self, q):
        if self.tot <= 0:
            return np.nan
        cw = np.cumsum(self.w) / self.tot
        return float(self.r[min(np.searchsorted(cw, q), cw.size - 1)])

    def peak_frac(self, pixel_um):
        """Signal fraction per um^2 in the central pixel-sized disk."""
        if self.tot <= 0:
            return 0.0
        rc = pixel_um / np.sqrt(np.pi)
        cw = np.cumsum(self.w) / self.tot
        return float(np.interp(rc, self.r, cw) / (np.pi * rc * rc))


def snr_cycle(cam, K, sig_frac, R_on, R_off, EE, peak, bg_dens, eta=1.0, filt=1.0, rho_max_um=None):
    """Per-blink-cycle SNR for a camera with SAMPLES_PER_CYCLE frames per cycle.

    Ideal synchronous detection: the 'on' half and 'off' half of each cycle are each
    SAMPLES_PER_CYCLE/2 frames. The brightest pixel is kept below full well by
    attenuating the light (factor k). Read noise is added for every frame and pixel.
    Returns the best ROI and the electron budget (used by the time-domain decoder sim).
    """
    p = cam["obj_pixel_um"]
    rho = Image.RHO_UM
    ok = rho >= p / np.sqrt(np.pi)
    if rho_max_um is not None:
        ok &= rho <= max(rho_max_um, p)
    area = np.pi * rho ** 2
    s_on = K * sig_frac * R_on * EE * eta
    s_off = K * sig_frac * R_off * EE * eta
    B = K * bg_dens * area / filt
    peak_pix = K * (sig_frac * R_on * eta * peak + bg_dens / filt) * p * p
    k = min(1.0, cam["full_well"] / peak_pix) if peak_pix > 0 else 1.0
    half = C.SAMPLES_PER_CYCLE / 2
    npix = np.maximum(1.0, area / (p * p))
    diff = half * k * (s_on - s_off)
    var = half * k * (s_on + s_off + 2 * B) + C.SAMPLES_PER_CYCLE * npix * cam["read_noise"] ** 2
    snr = np.where(ok, diff / np.sqrt(var), -np.inf)
    i = int(np.argmax(snr))
    m = (s_on[i] - s_off[i]) / (0.5 * (s_on[i] + s_off[i]) + B[i]) if (s_on[i] + B[i]) > 0 else 0.0
    return dict(snr1=float(max(snr[i], 0.0)), roi_um=float(rho[i]), m=float(m), sat_k=float(k),
                npix=float(npix[i]), e_on_frame=float(k * s_on[i]), e_off_frame=float(k * s_off[i]),
                e_bg_frame=float(k * B[i]))


def cycles_to(snr1, target=C.SNR_TARGET):
    if not np.isfinite(snr1) or snr1 <= 0:
        return np.inf
    return max(1.0, np.ceil((target / snr1) ** 2))


def surface_collection(NA, tilt_deg, mode, n_layer=C.N_MEDIUM):
    """Fraction of E*a^2*R collected from a SPOT with no tissue in the path.

    Specular: illumination fills the lens NA (Koehler/coaxial); a tilt shifts the reflected
    cone by n*sin(2*tilt) in NA space; collected fraction = cos(tilt) * overlap of two disks
    of radius NA. Lambertian: cos(tilt)^2 * (NA/n)^2 (small-NA limit).
    """
    t = np.radians(tilt_deg)
    if tilt_deg >= 90:
        return 0.0
    if mode == "lambert":
        return float(np.cos(t) ** 2 * (NA / n_layer) ** 2)
    s = n_layer * np.sin(min(2 * t, np.pi / 2))
    if 2 * t >= np.pi / 2 or s >= 2 * NA:
        return 0.0
    r = NA
    ov = 2 * r * r * np.arccos(s / (2 * r)) - 0.5 * s * np.sqrt(4 * r * r - s * s)
    return float(np.cos(t) * ov / (np.pi * r * r))
