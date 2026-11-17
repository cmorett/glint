"""Round-3 additions. No new physics: only new wavelengths, reflectance scalings, orientations
and post-processing. Everything else comes from config_v2 / config (rounds 1-2)."""
from config_v2 import *  # noqa: F401,F403
import config_v2 as C2

# Task 1: the measured reflectance change is at 635 nm. At 780/850 nm the paper gives no
# number, so we test the true dR being 100 %, 50 % or 25 % of the 635 nm value, keeping
# the mean reflectance unchanged (ASSUMPTION: only the swing is unknown).
_mean = 0.5 * (C2.R_ON_PAPER + C2.R_OFF_PAPER)
_half = 0.5 * (C2.R_ON_PAPER - C2.R_OFF_PAPER)
for _f in (0.5, 0.25):
    C2.REFLECTANCE_CASES[f"paper635_x{_f}"] = (_mean + _f * _half, _mean - _f * _half)
REFL_V3 = ["paper635", "paper635_x0.5", "paper635_x0.25"]
WAVELENGTHS_V3 = [635, 780, 850]
HEADLINE_WAVELENGTH_V3 = 635

# Task 2: random orientation of a chip at the centre of a sphere.
CENTRE_TILTS_DEG = [0, 10, 20, 30, 45, 60, 90, 120, 150, 180]
CENTRE_TILT_WAVELENGTHS = [635, 780]
N_ORIENTATIONS = 20000
N2_TILT_SPHERE = 200_000

# Task 3: flatness geometry.
CHIP_SIDES_UM = [150, 200, 300]
FLAT_NA = {"plate20MP_NA0.02": 0.02, "plate20MP_NA0.05": 0.05, "16x20MP_NA0.08": 0.08, "scan_NA0.25": 0.25}
DEBRIS_UM = 5.0            # "a few um of debris" (ASSUMPTION for the robustness statement)

# Task 4: added noise at the blink frequency, expressed as the apparent modulation it
# produces in a 1-second average (1/sqrt(Hz)). 0 = round-2 model. ASSUMPTION: range.
BLINK_NOISE = [0.0, 1e-4, 1e-3, 1e-2]
