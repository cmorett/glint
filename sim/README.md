# Glint light-transport simulation (`glint/sim`)

A Monte Carlo light-transport model that asks whether a camera can read a blinking SPOT chip on, under or inside an organoid, and how deep the chip can be. The verdict is in `results/SUMMARY.md`.

## Run

```
python run_all.py            # validation, then the full sweep (142 s on 16 cores; ~30 min on 1 core)
python run_all.py --quick    # 10x fewer photons, smoke test (~50 s)
python make_tables.py        # rebuild results/tables.md from the CSVs
pytest tests/                # validation checks alone
```

You need Python 3 and NumPy. Numba is not used. Work runs in parallel over CPU cores (`multiprocessing`). All assumptions live in `config.py`.

Outputs in `results/`:
- `results.csv`: one row per configuration (placement × depth × SPOT size × wavelength × μs′ × μa × NA × irradiance × ΔR). Each row has the mean over 3 seeds and a `_err` column (standard error).
- `results_wavelength_shift.csv`: the hypothetical wavelength-shift case.
- `blur_point_mirror.csv`, `sanity_chicken_skin.csv`, `validation.txt`, `tables.md`, `run_info.txt`.

## Method

1. **Core (`mc.py`).** This is MCML-style transport in one slab. It uses exponential free paths, Henyey–Greenstein scattering (g = 0.9), all-or-none Fresnel reflection and refraction at both faces, and Russian roulette. There are two absorption schemes. Classic MCML weight-drop is used only for validation. Production runs use absorption weighting along the path, w = exp(−μa·L): steps are sampled with μs only and the path length L is kept, so one run scores μa = 0.01, 0.02 and 0.05 exactly.
2. **Stage 1 (`stages.stage1`).** A collimated beam hits the slab at normal incidence. Wide-field illumination of an infinite slab is translation-invariant, so this is a 1-D problem. At every chip depth d it scores the downward irradiance F_down(d), its unscattered part, and the polar-angle histogram of the scattered part. It also scores diffuse reflectance leaving the camera side inside each NA, which is the background. The NA test uses the invariant n·sinθ, which is conserved across flat interfaces.
3. **Stage 2 (`stages.stage2`).** Photons start at a point mirror at depth d. Their directions come from the Stage 1 distribution, reflected specularly. Each photon is tracked until it leaves. The output is the collection fraction C_NA and the exit positions. Each exit ray is also traced back to the chip plane, which is where a camera focused on the chip forms its image.
4. **SPOT size.** Each collected photon gets a uniform random offset inside the a × a square. This is a Monte Carlo convolution of the point-mirror PSF with the SPOT footprint.
5. **Detection (`analysis.py`).** Absolute electron counts per frame come from irradiance, wavelength, QE = 0.6 and a 3.33 ms exposure (one half-cycle at 150 Hz):
   - signal: S = E·a²·F_down·R·C_NA·EE(ρ)·T_outer²
   - background: B = E·πρ²·[(r_sp + R_d,NA)·T_outer² + glare]

   The ROI radius ρ is chosen to maximise SNR, which is standard aperture photometry. There are two detector models:
   - **camera** (sCMOS): 5 µm object-space pixels, 30 000 e⁻ full well, 1.6 e⁻ read noise. The scene is dimmed until the brightest pixel stops saturating.
   - **ideal**: pure shot noise, for example a photodiode with lock-in and a matched field stop.

   SNR per cycle = ΔS/√(S_on + S_off + 2B). The number of cycles needed is N = (5/SNR)², and the time is N/150 s. A configuration is readable under 60 s, marginal under 10 min, and not readable beyond that.
6. **Validation (`tests/`).** The checks are:
   - the MCML benchmark in both modes
   - the zero-scattering analytic limit, including multiple Fresnel bounces
   - energy conservation
   - agreement between the two modes
   - the paper's chicken-skin observation

   `run_all.py` stops if a hard check fails.

Simplifications: flat laterally infinite slab instead of a sphere; chip does not shadow or perturb the light field; ideal aberration-free optics focused on the chip plane; collimated normal illumination; shot noise (+ camera full-well/read noise) only; shifted wavelength assumed to see the same tissue optics.
