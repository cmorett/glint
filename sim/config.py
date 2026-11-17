"""All assumptions for the Glint light-transport simulation live here.

Every number is tagged SOURCE (from a reference) or ASSUMPTION (our choice).
Change values here; nothing else in the code hard-codes physics parameters.
Lengths are in mm inside the Monte Carlo; reported in um in results.
"""

# ---------------------------------------------------------------- geometry
# ASSUMPTION: organoid modelled as a flat slab, infinite laterally, of
# thickness T = organoid diameter. 2 mm lets us place a chip at 1000 um
# (the centre). A real organoid is a sphere: light also leaks out its sides.
SLAB_THICKNESS_MM = 2.0

N_TISSUE = 1.37   # SOURCE: prompt; typical soft-tissue index (Jacques 2013 lists 1.33-1.45)
N_MEDIUM = 1.33   # culture medium ~ water
N_DISH = 1.50     # ASSUMPTION: clear glass-bottom dish / plate
N_AIR = 1.00

# Placements: (name, depth list in um, index on camera side, index on far side)
#   A: chip on top of tissue, camera above.         d = 0
#   B: chip under tissue on the dish, camera below.  d = 0 (from dish side)
#   C: chip embedded, camera above, d below the top surface.
DEPTHS_C_UM = [50, 100, 150, 200, 300, 500, 1000]
PLACEMENTS = {
    "A_top":       dict(depths_um=[0],          n_cam=N_MEDIUM, n_far=N_DISH),
    "B_under":     dict(depths_um=[0],          n_cam=N_DISH,   n_far=N_MEDIUM),
    "C_embedded":  dict(depths_um=DEPTHS_C_UM,  n_cam=N_MEDIUM, n_far=N_DISH),
}
# Specular "glare" from the outermost flat interface on the camera side,
# uniform over the field for collimated epi-illumination (Fresnel, normal incidence).
#   A/C: air -> medium meniscus (treated as flat). B: air -> dish bottom.
def fresnel_normal(n1, n2):
    return ((n1 - n2) / (n1 + n2)) ** 2
GLARE = {
    "A_top": fresnel_normal(N_AIR, N_MEDIUM),
    "C_embedded": fresnel_normal(N_AIR, N_MEDIUM),
    "B_under": fresnel_normal(N_AIR, N_DISH),
}

# ---------------------------------------------------------------- tissue optics
WAVELENGTHS_NM = [635, 780, 850]
G = 0.9                    # SOURCE: prompt; typical soft tissue anisotropy
# Reduced scattering is specified AT 780 nm (the paper's tissue wavelength) and
# scaled to other wavelengths with mu_s'(lam) = mu_s'(780) * (lam/780)^-b.
# SOURCE for b: Jacques 2013 (Phys Med Biol 58 R37), Table 2, "other soft
# tissues" mean b = 1.286 (SD 0.521). For reference that row's a = 18.9 cm^-1
# gives mu_s'(780 nm) = 1.07 mm^-1, consistent with our 1.0 baseline.
MUSP_780_PER_MM = [0.3, 1.0, 2.0]   # ASSUMPTION (prompt): loose / baseline / dense
MUSP_BASELINE = 1.0
SCATTER_POWER_B = 1.286
# ASSUMPTION: absorption is low (no blood) and taken wavelength-independent.
MUA_PER_MM = [0.01, 0.02, 0.05]
MUA_BASELINE = 0.02


def musp_at(musp_780, wavelength_nm):
    return musp_780 * (wavelength_nm / 780.0) ** (-SCATTER_POWER_B)


# Sanity case: the SPOT paper imaged a 50x50 um SPOT through ~250 um of
# chicken skin at 780 nm (arXiv:2606.23916, Section S9). Skin optics:
# Jacques 2013 Table 2 "skin" a = 46.0 cm^-1, b = 1.421 -> mu_s'(780) = 2.45 mm^-1.
# Top surface wet skin in air (n_cam = 1.0); chip directly under the skin.
CHICKEN = dict(musp_780=4.60 * (780 / 500.0) ** (-1.421), mua=0.02, thickness_mm=0.25,
               depth_um=250, n_cam=N_AIR, n_far=N_TISSUE, spots_um=[50, 300])

# ---------------------------------------------------------------- SPOT
SPOT_SIZES_UM = [30, 100, 300]
R_MEAN = 0.5                     # ASSUMPTION (prompt): absolute paper values not extracted
DELTA_R = [0.05, 0.1, 0.2]       # ASSUMPTION (prompt): swept R_on - R_off
DELTA_R_BASELINE = 0.1
BLINK_HZ = 150.0                 # SOURCE: paper, efficient around 150-300 Hz

# ---------------------------------------------------------------- illumination
# ASSUMPTION: collimated, normally incident, uniform wide-field epi-illumination.
IRRADIANCE_MW_CM2 = [1.0, 10.0, 100.0]   # ASSUMPTION: cell-safety-limited sweep
IRRADIANCE_BASELINE = 10.0

# ---------------------------------------------------------------- detection
NA_LIST = [0.1, 0.25, 0.5]
NA_BASELINE = 0.25
QE = 0.6                          # SOURCE: prompt
EXPOSURE_S = 1.0 / (2 * BLINK_HZ)  # one half-cycle per frame (camera at 300 fps)
# Camera model (ASSUMPTIONS, typical scientific sCMOS):
PIXEL_OBJ_UM = 5.0     # object-space pixel size (e.g. 6.5 um pixel at 1.3x)
FULL_WELL_E = 30000.0  # electrons per pixel per frame
READ_NOISE_E = 1.6     # electrons rms per pixel per frame
SNR_TARGET = 5.0
READABLE_S = 60.0
MARGINAL_S = 600.0
# ASSUMPTION: below this modulation depth, real-world drifts (LED noise, motion,
# fixed-pattern) likely dominate shot noise. Reported as a flag only.
CONTRAST_FLOOR = 1e-3   # tables also report the looser 1e-4

# ---------------------------------------------------------------- wavelength shift (hypothetical)
SHIFT_FILTER = [1e2, 1e3, 1e4]     # background suppression factor
SHIFT_EFFICIENCY = [0.01, 0.1, 1.0]  # fraction of reflected signal converted

# ---------------------------------------------------------------- power check
SOLAR_CELL_UM = 100.0
CHIP_NEED_NW_ELECTRICAL = 100.0   # ASSUMPTION: our own rough estimate
PV_EFFICIENCY = 0.10              # ASSUMPTION: our own rough estimate

# ---------------------------------------------------------------- Monte Carlo
SEEDS = [1, 2, 3]
N_PHOTONS_STAGE1 = 300_000
N_PHOTONS_STAGE2 = 300_000
ROULETTE_THRESHOLD = 1e-4
ROULETTE_SURVIVAL = 0.1
THETA_BINS = 180                  # 0.5 deg bins for the downward angular distribution
N_WORKERS = None                  # None -> os.cpu_count() - 1
