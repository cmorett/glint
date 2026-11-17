"""Round-2 assumptions. Imports round 1 (config.py) unchanged and adds to it.

Tags: SOURCE (reference), ASSUMPTION (our choice), DERIVED (computed from others).
"""
import math

from config import *  # noqa: F401,F403  (round-1 values stay the single source of truth)
import config as C1

# ---------------------------------------------------------------- Task 6: reflectance
# SOURCE: Campbell et al., arXiv:2606.23916 (file Y-Prize_Technology_Article_SPOTs.pdf),
# Supplementary Section S4, Fig. S11 and Eq. S40: sigmoid fit of reflectance vs potential
# for a 50x50 um SPOT, ~45 nm PEDOT:PSS on Pt, 1x PBS, 635 nm laser:
#   rho(V) = chi1 + (chi2 - chi1) / (1 + exp(chi3 (V - chi4)))
#   chi1 = 47.0 %, chi2 = 32.5 %, chi3 = 0.005 /mV, chi4 = -400 mV.
# Evaluated at the paper's standard drive limits (+200 mV / -800 mV, main text Fig. 1e):
PAPER_FIT = dict(chi1=0.470, chi2=0.325, chi3=0.005, chi4=-400.0)


def paper_reflectance(v_mV, f=PAPER_FIT):
    return f["chi1"] + (f["chi2"] - f["chi1"]) / (1.0 + math.exp(f["chi3"] * (v_mV - f["chi4"])))


R_ON_PAPER = paper_reflectance(200.0)    # 0.463
R_OFF_PAPER = paper_reflectance(-800.0)  # 0.342
# The paper gives no number at 780/850 nm; it says 635 nm is near the wavelength where
# the +200/-800 mV reflectance ratio is largest (Section S2.1, Fig. S12). Applying the
# 635 nm value at 780 nm is therefore an OPTIMISTIC ASSUMPTION, labelled wherever used.
REFLECTANCE_CASES = {
    "paper635": (R_ON_PAPER, R_OFF_PAPER),       # measured (635 nm)
    "sweep_dR0.05": (0.525, 0.475),              # round-1 sweep, kept for comparison
    "sweep_dR0.2": (0.60, 0.40),
}
HEADLINE_REFLECTANCE = "paper635"

# ---------------------------------------------------------------- wavelengths in round 2
# 635 nm: where the paper measured dR. 780 nm: round-1 baseline (paper's tissue imaging).
WAVELENGTHS_V2 = [635, 780]
HEADLINE_WAVELENGTH = 780

# ---------------------------------------------------------------- Task 1: cameras / optics
# Plate field (SOURCE: ANSI/SLAS 96 and 384 plates share the 127.76 x 85.48 mm footprint;
# the well area is about 110 x 75 mm, as given in the round-2 prompt).
PLATE_FIELD_MM = (110.0, 75.0)
WELLS = {"96": 96, "384": 384}


def obj_pixel_um(sensor_px):
    """Object-space pixel so the whole field fits the sensor (DERIVED)."""
    return max(PLATE_FIELD_MM[0] * 1000 / sensor_px[0], PLATE_FIELD_MM[1] * 1000 / sensor_px[1])


# ASSUMPTIONS: typical datasheet-class sensors, not specific products.
SENSORS = {
    "20MP": dict(px=(5472, 3648), pixel_um=2.4, full_well=15000.0, read_noise=2.0, fps=20.0),
    "60MP": dict(px=(9576, 6388), pixel_um=3.76, full_well=50000.0, read_noise=3.0, fps=5.0),
}
for _k, _s in SENSORS.items():
    _s["obj_pixel_um"] = obj_pixel_um(_s["px"])
    _s["magnification"] = _s["pixel_um"] / _s["obj_pixel_um"]
# Lens: object-side TELECENTRIC lens with coaxial illumination (needed so that every well
# sees normal incidence; with an ordinary macro lens a flat mirror at the field edge
# reflects the illumination outside the lens). Object NA = magnification / (2 * f-number):
#   20MP (m ~ 0.12): f/1.2..f/6  ->  NA 0.05..0.01;   60MP (m ~ 0.32): f/3.2..f/16 -> NA 0.05..0.01
LOW_NA = [0.01, 0.02, 0.05]
NA_LIST_V2 = [0.01, 0.02, 0.05, 0.1, 0.25, 0.5]
# Scanning head (ASSUMPTION): NA 0.25 objective, sCMOS in ROI mode.
SCAN_HEAD = dict(NA=0.25, obj_pixel_um=5.0, full_well=30000.0, read_noise=1.6, fps=600.0,
                 move_settle_s=0.5)   # ASSUMPTION: well-to-well move + settle + autofocus
SAMPLES_PER_CYCLE = 4       # asynchronous lock-in needs >= 4 frames per blink cycle
NA_EFF_SCALING = 0.05       # below this NA, non-ballistic light is scaled as NA^2 (checked in tests)

CAMERAS = {"scan_NA0.25": dict(SCAN_HEAD, kind="scan")}
for _k, _s in SENSORS.items():
    for _na in LOW_NA:
        CAMERAS[f"plate{_k}_NA{_na}"] = dict(kind="plate", NA=_na, obj_pixel_um=_s["obj_pixel_um"],
                                            full_well=_s["full_well"], read_noise=_s["read_noise"],
                                            fps=_s["fps"])
for _c in CAMERAS.values():
    _c["blink_hz"] = _c["fps"] / SAMPLES_PER_CYCLE

# ---------------------------------------------------------------- Task 2: spheres
SPHERE_DIAMETERS_MM = [0.5, 1.0, 2.0]
SPHERE_SCORING_RADIUS_MM = 0.05   # irradiance at the chip averaged over this disk
FAR_SIDE_GAP_MM = 0.01            # "chip at the far surface" is placed this far inside

# ---------------------------------------------------------------- Task 3: tilt
TILTS_DEG = [0, 5, 10, 20, 45, 90, 135, 180]
LAMBERT_TILTS_DEG = [0, 20, 45, 90, 135, 180]
TILT_DEPTHS_UM = [50, 100, 300]
SETTLE_TILT_DEG = 2.0   # ASSUMPTION: residual tilt of a chip settled on glass (debris, roughness)

# ---------------------------------------------------------------- Task 4: noise model
SOURCE_NOISE_RMS = [1e-4, 1e-3, 1e-2]   # prompt: 0.01 %, 0.1 %, 1 %
DRIFT_TAU_S = 2.0                       # ASSUMPTION: correlation time of slow drift
PRNU = 0.01                             # ASSUMPTION: 1 % pixel gain non-uniformity
DSNU_E = 1.0                            # ASSUMPTION: static offset non-uniformity (e-)
MOTION_UM = 3.0                         # prompt: slow sample motion of a few um
MOTION_TAU_S = 5.0
BG_TEXTURE = 0.2                        # ASSUMPTION: +/-20 % spatial texture in tissue background
CLOCK_TOL = 0.003                       # ASSUMPTION: chip oscillator within +/-0.3 %
PREAMBLE = [1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1]   # Barker-13
ID_BITS, VALUE_BITS = 8, 8
FRAME_BITS = len(PREAMBLE) + ID_BITS + VALUE_BITS     # 29 bits, DBPSK
DECODE_TRIALS = 8
READ_LIMIT_S = 60.0                     # a chip "decodes reliably" if a frame decodes in <= 60 s

# ---------------------------------------------------------------- Task 5: dose
READ_INTERVALS_MIN = [1, 5, 15]
CHIP_STARTUP_S = 0.05    # ASSUMPTION: PV charge-up + oscillator start before the first bit

# ---------------------------------------------------------------- Monte Carlo sizes
SEEDS_V2 = [1, 2, 3]
N1_SLAB = 1_000_000
N2_SLAB = 600_000
N2_TILT = 300_000
N1_SPHERE = 600_000
N2_SPHERE = 300_000

# Extra optics for comparison (same sensor as the scan head) and multi-camera plate readers:
# k cameras each imaging 1/k of the plate -> magnification and NA x sqrt(k), pixel / sqrt(k).
for _na in (0.1, 0.5):
    CAMERAS[f"scan_NA{_na}"] = dict(SCAN_HEAD, NA=_na, kind="scan")
for _k in (4, 16):
    _s = SENSORS["20MP"]
    CAMERAS[f"{_k}x20MP_NA{0.02 * math.sqrt(_k):.2f}"] = dict(
        kind="multi", NA=0.02 * math.sqrt(_k), obj_pixel_um=_s["obj_pixel_um"] / math.sqrt(_k),
        full_well=_s["full_well"], read_noise=_s["read_noise"], fps=_s["fps"], n_cameras=_k)
for _c in CAMERAS.values():
    _c["blink_hz"] = _c["fps"] / SAMPLES_PER_CYCLE

DEPTHS_V2_UM = [50, 100, 150, 200, 250, 300, 500, 1000]
HEADLINE_CAMERAS = ["scan_NA0.25", "plate20MP_NA0.02"]
