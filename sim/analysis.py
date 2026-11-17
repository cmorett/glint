"""Turn Stage 1 / Stage 2 tallies into signal, blur, background, SNR and verdicts."""
import numpy as np

import config as C

H = 6.62607015e-34
CLIGHT = 2.99792458e8
RHO_GRID_UM = np.geomspace(2.0, 4000.0, 140)    # candidate ROI radii (image of chip plane)


def electrons_per_um2_per_frame(E_mW_cm2, lam_nm):
    """QE * photon fluence per frame per um^2 for irradiance E (incident, in air)."""
    flux = E_mW_cm2 * 10.0 / (H * CLIGHT / (lam_nm * 1e-9)) / 1e12   # photons / um^2 / s
    return C.QE * flux * C.EXPOSURE_S


def weighted_quantiles(r, w, qs):
    if r.size == 0 or w.sum() <= 0:
        return [np.nan] * len(qs)
    o = np.argsort(r)
    cw = np.cumsum(w[o]); cw /= cw[-1]
    return [float(r[o][min(np.searchsorted(cw, q), r.size - 1)]) for q in qs]


class SpotImage:
    """Image (in the chip plane) of a uniformly lit square mirror of side a.

    Built by Monte Carlo convolution: each collected photon from the point-mirror
    run is shifted by a uniform random offset inside the square.
    """

    def __init__(self, xb_um, yb_um, w, a_um, rng, min_samples=200_000):
        n = xb_um.size
        k = max(1, int(np.ceil(min_samples / max(n, 1))))
        xs = np.repeat(xb_um, k) + (rng.random(n * k) - 0.5) * a_um
        ys = np.repeat(yb_um, k) + (rng.random(n * k) - 0.5) * a_um
        r = np.hypot(xs, ys)
        self.order = np.argsort(r)
        self.r = r[self.order]
        self.k = k
        self.w_all = np.repeat(w, k, axis=-1)[..., self.order] / k   # (n_mua, n*k)
        self.a = a_um

    def stats(self, j):
        w = self.w_all[j]
        tot = w.sum()
        if tot <= 0:
            return dict(EE=np.zeros_like(RHO_GRID_UM), r50=np.nan, r90=np.nan, fwhm=np.nan, peak=0.0)
        cw = np.cumsum(w) / tot
        EE = np.interp(RHO_GRID_UM, self.r, cw, left=0.0, right=1.0)
        r50 = float(self.r[min(np.searchsorted(cw, 0.5), cw.size - 1)])
        r90 = float(self.r[min(np.searchsorted(cw, 0.9), cw.size - 1)])
        # radially averaged profile -> FWHM (approximate, MC noise)
        edges = np.linspace(0, max(1.5 * r90, 2.0), 76)
        h, _ = np.histogram(self.r, bins=edges, weights=w)
        dens = h / (np.pi * (edges[1:] ** 2 - edges[:-1] ** 2))
        c0 = dens[0]
        below = np.nonzero(dens < 0.5 * c0)[0]
        fwhm = float(2 * 0.5 * (edges[below[0]] + edges[below[0] + 1])) if below.size else np.nan
        # central density (fraction of signal per um^2) over one pixel area
        rc = C.PIXEL_OBJ_UM / np.sqrt(np.pi)
        peak = float(np.interp(rc, self.r, cw) / (np.pi * rc * rc))
        return dict(EE=EE, r50=r50, r90=r90, fwhm=fwhm, peak=peak)


def snr_models(K, sig_frac, R_on, R_off, EE, peak_frac, bg_dens, eta=1.0, filt=1.0):
    """Shot-noise SNR per blink cycle on the ROI grid; returns best ROI for two detectors.

    K        electrons per um^2 per frame per unit reflectance (incident irradiance)
    sig_frac a^2 * F_down * C_NA * T_outer^2  (um^2): signal per unit R before ROI
    EE       encircled fraction of signal for each ROI radius
    peak_frac signal fraction per um^2 at the image centre
    bg_dens  background per um^2 as a fraction of incident (tissue + glare)
    """
    area = np.pi * RHO_GRID_UM ** 2
    s_on = K * sig_frac * R_on * EE * eta
    s_off = K * sig_frac * R_off * EE * eta
    B = K * bg_dens * area / filt
    dS = s_on - s_off
    var = s_on + s_off + 2 * B
    with np.errstate(divide="ignore", invalid="ignore"):
        snr_ideal = np.where(var > 0, dS / np.sqrt(var), 0.0)
    i_id = int(np.argmax(snr_ideal))
    # camera: brightest pixel must not exceed full well -> attenuate everything by k
    peak_pix = K * (sig_frac * R_on * eta * peak_frac + bg_dens / filt) * C.PIXEL_OBJ_UM ** 2
    k = min(1.0, C.FULL_WELL_E / peak_pix) if peak_pix > 0 else 1.0
    npix = np.maximum(1.0, area / C.PIXEL_OBJ_UM ** 2)
    var_c = k * var + 2 * npix * C.READ_NOISE_E ** 2
    snr_cam = k * dS / np.sqrt(var_c)
    i_c = int(np.argmax(snr_cam))
    m = dS[i_c] / (0.5 * (s_on[i_c] + s_off[i_c]) + B[i_c]) if (s_on[i_c] + B[i_c]) > 0 else 0.0
    return dict(snr_cam=float(snr_cam[i_c]), roi_cam_um=float(RHO_GRID_UM[i_c]),
                snr_ideal=float(snr_ideal[i_id]), roi_ideal_um=float(RHO_GRID_UM[i_id]),
                m=float(m), S_on=float(s_on[i_c]), S_off=float(s_off[i_c]), B=float(B[i_c]),
                sat_atten=float(k))


def cycles_needed(snr1):
    if not np.isfinite(snr1) or snr1 <= 0:
        return np.inf
    return max(1.0, np.ceil((C.SNR_TARGET / snr1) ** 2))


def verdict(t):
    if t < C.READABLE_S:
        return "readable"
    if t < C.MARGINAL_S:
        return "marginal"
    return "not_readable"


def analyse_setting(lam, musp780, s1_by_geom, s2_by_depth, seed_rng, seed):
    """All rows for one (wavelength, mu_s') tissue setting and one seed."""
    rows, shift_rows, blur_rows = [], [], []
    mua_list = list(C.MUA_PER_MM)
    for place, P in C.PLACEMENTS.items():
        s1 = s1_by_geom[place]
        glare = C.GLARE[place]
        To = 1.0 - glare
        for d_um in P["depths_um"]:
            d_mm = d_um / 1000.0
            if d_um == 0:
                F = np.ones(len(mua_list)); Fb = np.ones(len(mua_list))
                n_ex = 1
                exits = None
            else:
                k_pl = int(np.argmin(np.abs(s1["planes"] - d_mm)))
                F = s1["F_down"][:, k_pl]; Fb = s1["F_ball"][:, k_pl]
                exits = s2_by_depth[d_um]
            # point-mirror exit-surface blur (all exiting light)
            for j, mua in enumerate(mua_list):
                if exits is None:
                    r50p = r90p = 0.0; Rexit = 1.0
                else:
                    rr = np.hypot(exits["x"], exits["y"]) * 1000.0
                    r50p, r90p = weighted_quantiles(rr, exits["w"][j].astype(float), [0.5, 0.9])
                    Rexit = float(exits["w"][j].sum() / exits["N"])
                blur_rows.append(dict(placement=place, depth_um=d_um, wavelength_nm=lam,
                                      musp_780=musp780, mua=mua, seed=seed,
                                      point_r50_surface_um=r50p, point_r90_surface_um=r90p,
                                      exit_fraction_all_angles=Rexit))
            for na_i, na in enumerate(C.NA_LIST):
                if exits is None:
                    Cna = np.ones(len(mua_list))
                    xb = np.zeros(1); yb = np.zeros(1); w = np.ones((len(mua_list), 1))
                else:
                    sel = exits["inv"] <= na
                    w = exits["w"][:, sel].astype(float)
                    Cna = w.sum(1) / exits["N"]
                    xb = exits["x_bp"][sel].astype(float) * 1000.0
                    yb = exits["y_bp"][sel].astype(float) * 1000.0
                bg_tissue = s1["r_sp"] + s1["Rd_na"][:, na_i]
                for a in C.SPOT_SIZES_UM:
                    img = SpotImage(xb, yb, w, a, seed_rng)
                    for j, mua in enumerate(mua_list):
                        st = img.stats(j)
                        sig_frac = a * a * F[j] * Cna[j] * To * To
                        bg_dens = bg_tissue[j] * To * To + glare
                        base = dict(placement=place, depth_um=d_um, spot_um=a, wavelength_nm=lam,
                                    musp_780=musp780, musp=C.musp_at(musp780, lam), mua=mua, NA=na,
                                    seed=seed, F_down=F[j], F_ballistic=Fb[j], C_NA=Cna[j],
                                    signal_frac=F[j] * Cna[j],
                                    img_r50_um=st["r50"], img_r90_um=st["r90"], img_fwhm_um=st["fwhm"],
                                    bg_tissue_NA=bg_tissue[j], glare=glare, bg_per_area=bg_dens,
                                    power_on_cell_nW_per_mWcm2=100.0 * F[j] * To,
                                    E_needed_mWcm2=(C.CHIP_NEED_NW_ELECTRICAL / C.PV_EFFICIENCY)
                                    / (100.0 * F[j] * To) if F[j] > 0 else np.inf)
                        for E in C.IRRADIANCE_MW_CM2:
                            K = electrons_per_um2_per_frame(E, lam)
                            for dR in C.DELTA_R:
                                Ron, Roff = C.R_MEAN + dR / 2, C.R_MEAN - dR / 2
                                res = snr_models(K, sig_frac, Ron, Roff, st["EE"], st["peak"], bg_dens)
                                Nc, Ni = cycles_needed(res["snr_cam"]), cycles_needed(res["snr_ideal"])
                                row = dict(base, irradiance_mWcm2=E, dR=dR, **res,
                                           N_cam=Nc, t_cam_s=Nc / C.BLINK_HZ,
                                           N_ideal=Ni, t_ideal_s=Ni / C.BLINK_HZ,
                                           power_on_cell_nW=100.0 * F[j] * To * E)
                                rows.append(row)
                                if dR == C.DELTA_R_BASELINE and mua == C.MUA_BASELINE:
                                    for filt in C.SHIFT_FILTER:
                                        for eta in C.SHIFT_EFFICIENCY:
                                            rs = snr_models(K, sig_frac, Ron, Roff, st["EE"], st["peak"],
                                                            bg_dens, eta=eta, filt=filt)
                                            Nc2, Ni2 = cycles_needed(rs["snr_cam"]), cycles_needed(rs["snr_ideal"])
                                            shift_rows.append(dict(
                                                placement=place, depth_um=d_um, spot_um=a,
                                                wavelength_nm=lam, musp_780=musp780, mua=mua, NA=na,
                                                irradiance_mWcm2=E, dR=dR, filter=filt, efficiency=eta,
                                                seed=seed, m=rs["m"], snr_cam=rs["snr_cam"],
                                                snr_ideal=rs["snr_ideal"], N_cam=Nc2,
                                                t_cam_s=Nc2 / C.BLINK_HZ, N_ideal=Ni2,
                                                t_ideal_s=Ni2 / C.BLINK_HZ))
    return rows, shift_rows, blur_rows
