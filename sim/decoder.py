"""Reusable decoder for blinking SPOT chips in a camera frame sequence.

Pipeline:
  1. locate   - per-pixel lock-in power at the nominal blink frequency (short windows,
                summed incoherently), robust threshold, local maxima -> chip candidates.
  2. combine  - pixels around each candidate weighted by estimated amplitude/noise (matched filter).
  3. reference- pixels with no carrier power, far from chips; the chip signal is divided by
                the reference to cancel common-mode light-source fluctuations (optional).
  4. lock-in  - search carrier frequency (+/- tolerance) and symbol timing that maximise
                the symbol energy; complex demodulation per symbol.
  5. DBPSK    - differential detection between consecutive symbols; soft bits folded
                modulo the frame length; cyclic correlation with the Barker-13 preamble
                gives frame alignment; ID and value are read out.
Works on any array frames[t, y, x]; nothing in it knows the ground truth.
"""
import numpy as np

import config_v2 as C


def _windows_power(frames, fps, f0, win_frames):
    n = frames.shape[0]
    t = np.arange(n) / fps
    ph = np.exp(-2j * np.pi * f0 * t).astype(np.complex64)
    P = np.zeros(frames.shape[1:], np.float64)
    nw = 0
    for a in range(0, n - win_frames + 1, win_frames):
        blk = frames[a:a + win_frames]
        blk = blk - blk.mean(axis=0, keepdims=True)
        z = np.tensordot(ph[a:a + win_frames], blk, axes=(0, 0))
        P += np.abs(z) ** 2
        nw += 1
    return P / max(nw, 1)


def _box(P, r=1):
    out = np.zeros_like(P)
    ny, nx = P.shape
    cnt = np.zeros_like(P)
    for dy in range(-r, r + 1):
        for dx in range(-r, r + 1):
            ys = slice(max(0, dy), ny + min(0, dy)); yd = slice(max(0, -dy), ny + min(0, -dy))
            xs = slice(max(0, dx), nx + min(0, dx)); xd = slice(max(0, -dx), nx + min(0, -dx))
            out[yd, xd] += P[ys, xs]; cnt[yd, xd] += 1
    return out / cnt


def locate(frames, fps, f0, cycles_per_window=16, spc=C.SAMPLES_PER_CYCLE, k_sigma=8.0,
           radius_px=6, max_chips=20):
    P = _windows_power(frames, fps, f0, max(spc, int(cycles_per_window * spc)))
    Ps = _box(P, 2)
    med = np.median(Ps)
    mad = np.median(np.abs(Ps - med)) * 1.4826 + 1e-30
    thr = med + k_sigma * mad
    found = []
    work = Ps.copy()
    ny, nx = P.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    for _ in range(max_chips):
        i = int(np.argmax(work))
        if work.flat[i] <= thr:
            break
        y0, x0 = divmod(i, nx)
        found.append((y0, x0))
        work[(yy - y0) ** 2 + (xx - x0) ** 2 <= (2 * radius_px) ** 2] = -np.inf
    noise_med = np.median(P)
    chips = []
    for (y0, x0) in found:
        near = (yy - y0) ** 2 + (xx - x0) ** 2 <= radius_px ** 2
        amp = np.sqrt(np.maximum(Ps - med, 0)) * near      # smoothed map: less noisy weights
        w = amp / max(noise_med, 1e-30)
        if w.sum() <= 0:
            continue
        cy = float((w * yy).sum() / w.sum()); cx = float((w * xx).sum() / w.sum())
        chips.append(dict(y=cy, x=cx, weights=w / w.sum()))
    far = np.ones((ny, nx), bool)
    for (y0, x0) in found:
        far &= (yy - y0) ** 2 + (xx - x0) ** 2 > (2 * radius_px) ** 2
    ref_mask = far & (Ps < med + 2 * mad)
    return chips, ref_mask, Ps


def demodulate(y, fps, f0, n_c, tol=C.CLOCK_TOL, n_offsets=8):
    """Search carrier frequency and symbol timing; return complex symbol amplitudes."""
    n = y.size
    t = np.arange(n) / fps
    y = y - y.mean()
    best = (-1.0, None)
    T_nom = n_c / f0
    df = 0.05 / T_nom
    nf = int(np.ceil(1.2 * tol * f0 / df))
    for f in f0 + df * np.arange(-nf, nf + 1):
        base = y * np.exp(-2j * np.pi * f * t)
        cs = np.concatenate([[0], np.cumsum(base)])
        Ts = n_c / f
        for o in range(n_offsets):
            t0 = o / n_offsets * Ts
            k = np.arange(int((t[-1] - t0) / Ts))
            if k.size < 3:
                continue
            idx = np.clip(np.round((t0 + np.append(k, k[-1] + 1) * Ts) * fps).astype(int), 0, n)
            z = cs[idx[1:]] - cs[idx[:-1]]
            e = float(np.sum(np.abs(z) ** 2))
            if e > best[0]:
                best = (e, z, f, t0)
    return best[1], best[2]


def dbpsk_bits(z, n_bits=C.FRAME_BITS, preamble=C.PREAMBLE):
    d = -np.real(z[1:] * np.conj(z[:-1]))            # > 0 means a phase flip (bit 1)
    d = d / (np.abs(z[1:]) * np.abs(z[:-1]) + 1e-30)
    soft = np.zeros(n_bits); cnt = np.zeros(n_bits)
    for i, v in enumerate(d):
        soft[i % n_bits] += v; cnt[i % n_bits] += 1
    soft = soft / np.maximum(cnt, 1)
    pre = 2 * np.asarray(preamble) - 1
    corr = [np.dot(pre, soft[(r + np.arange(len(pre))) % n_bits]) for r in range(n_bits)]
    r = int(np.argmax(corr))
    srt = np.sort(corr)
    bits = (soft[(r + np.arange(n_bits)) % n_bits] > 0).astype(int)
    return bits, float(srt[-1] / len(pre)), float((srt[-1] - srt[-2]) / len(pre))


def bits_to_fields(bits):
    p = len(C.PREAMBLE)
    idv = int("".join(map(str, bits[p:p + C.ID_BITS])), 2)
    val = int("".join(map(str, bits[p + C.ID_BITS:p + C.ID_BITS + C.VALUE_BITS])), 2)
    pre_ok = list(bits[:p]) == list(C.PREAMBLE)
    return idv, val, pre_ok


def decode(frames, fps, f0, n_c, use_reference=True, radius_px=6):
    if use_reference:
        # common-mode correction before anything else: divide each frame by its median
        # (chips cover a small part of the field, so the median tracks the light source)
        med = np.median(frames.reshape(frames.shape[0], -1), axis=1)
        frames = frames / (med / med.mean())[:, None, None]
    chips, ref_mask, _ = locate(frames, fps, f0, radius_px=radius_px)
    ref = frames[:, ref_mask].sum(axis=1) if (use_reference and ref_mask.sum() > 0) else None
    out = []
    for c in chips:
        y = np.tensordot(frames, c["weights"], axes=([1, 2], [0, 1])).astype(np.float64)
        if ref is not None:
            y = y / ref * ref.mean()
        z, f = demodulate(y, fps, f0, n_c)
        if z is None:
            continue
        bits, score, margin = dbpsk_bits(z)
        idv, val, pre_ok = bits_to_fields(bits)
        out.append(dict(x=c["x"], y=c["y"], id=idv, value=val, preamble_ok=pre_ok, score=score,
                        margin=margin, f=f, bits=bits))
    return out
