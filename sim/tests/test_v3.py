"""Round-3 checks (post-processing only; rounds 1-2 physics checks still run first)."""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config_v3 as C        # noqa: E402
import make_tables_v3 as M3  # noqa: E402


class _StepRule:
    """Readable (1 s) below a cut-off SNR, never above; stands in for the calibrated rule."""
    def t_read(self, row):
        return 1.0 if row["snr1"] >= 1.0 else math.inf


def check_reflectance_scaling():
    on, off = C.REFLECTANCE_CASES["paper635"]
    ok = True
    for f in (0.5, 0.25):
        a, b = C.REFLECTANCE_CASES[f"paper635_x{f}"]
        ok &= abs((a + b) - (on + off)) < 1e-12 and abs((a - b) - f * (on - off)) < 1e-12
    return ok, f"dR scaling keeps mean R = {(on + off) / 2:.4f} and scales the swing (x0.5, x0.25)"


def check_orientation_sampler():
    # SNR high below 60 deg tilt, zero above: readable fraction must be (1 - cos 60)/2 = 0.25
    rows = [dict(tilt=t, snr1=(10.0 if t <= 60 else 1e-6), blink_hz=150.0, depth_um=0)
            for t in (0, 30, 59.999, 60.001, 90, 180)]
    s = M3.orientation_stats(_StepRule(), rows, n=200_000)
    ok = abs(s["f60"] - 0.25) < 0.005
    return ok, f"Random-orientation sampler: readable fraction {s['f60']:.4f} vs exact 0.25"


def check_blink_noise_zero():
    row = dict(snr1=3.7, e_on_frame=1.1e5, e_off_frame=0.9e5, e_bg_frame=5e4, blink_hz=150.0)
    a = M3.snr_with_blink_noise(row, 0.0)
    b = M3.snr_with_blink_noise(row, 1e-3)
    return abs(a - 3.7) < 1e-12 and b < a, f"Blink-band noise term: 0 noise reproduces SNR ({a:.3f}); 1e-3 lowers it ({b:.3f})"


def check_flatness_geometry():
    h = M3.h_max_um(150, 0.5)
    return abs(h - 150 * math.tan(math.radians(0.5))) < 1e-12 and 1.3 < h < 1.32, f"h_max(150 um, 0.5 deg) = {h:.3f} um"


def run_all_checks(verbose=False):
    out = []
    for f in (check_reflectance_scaling, check_orientation_sampler, check_blink_noise_zero, check_flatness_geometry):
        ok, msg = f()
        out.append(f"[{'PASS' if ok else 'FAIL'}] {msg}")
        if verbose:
            print(out[-1])
    return "\n".join(out) + "\n"


def test_reflectance_scaling():
    assert check_reflectance_scaling()[0]


def test_orientation_sampler():
    assert check_orientation_sampler()[0]


def test_blink_noise_zero():
    assert check_blink_noise_zero()[0]


def test_flatness_geometry():
    assert check_flatness_geometry()[0]


if __name__ == "__main__":
    print(run_all_checks(verbose=True))
