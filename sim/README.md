# Glint light-transport simulation (`sim/`)

A Monte Carlo light-transport model that asks whether a camera can read a blinking SPOT chip on, under or inside an organoid, and how deep the chip can be.

- Round 1 verdict: `results/SUMMARY.md` (flat slab, microscope-like optics).
- Round 2 verdict: `results/v2/SUMMARY_V2.md` (whole-plate optics, spheres, chip tilt, a working decoder, light dose, the paper's measured ΔR).
- Round 3 verdict (final, supersedes round 2's claims): `results/v3/SUMMARY_V3.md`, covering the headline at 635 nm, random orientation, flatness and a sensitivity ranking. Rounds 2 and 3 are described at the end of this file.

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

Round-1 simplifications: flat laterally infinite slab instead of a sphere; chip does not shadow or perturb the light field; ideal aberration-free optics focused on the chip plane; collimated normal illumination; shot noise (+ camera full-well/read noise) only; shifted wavelength assumed to see the same tissue optics.

## Round 2 (`run_v2.py`, outputs in `results/v2/`)

```
python run_v2.py             # round-1 + round-2 validation, Monte Carlo, decoder sims, tables (~8 min on 16 cores)
python bench_prediction.py   # predictions for the bench experiment (~1 min)
python run_v2.py --quick     # 10x fewer photons, fewer decoder trials (~2 min)
python make_tables_v2.py     # rebuild results/v2/tables_v2.md from the saved results
python decode_demo.py        # decoder demo on SYNTHETIC frames (see --help)
pytest tests/                # all round-1 and round-2 checks
```

Round 2 reuses `mc.py`, `stages.py`, `analysis.py` and `config.py` unchanged, apart from two additions to `mc.py`: it records the up-going angle histogram and a "launched unscattered" flag. New assumptions live in `config_v2.py`.

- **Reflectance (`config_v2.py`).** R_on/R_off = 0.463/0.342 come from the paper's sigmoid fit (Supplementary Fig. S11, Eq. S40; 45 nm PEDOT:PSS on Pt, 635 nm), evaluated at its +200/−800 mV drive. The paper reports no value at 780 nm. Using the 635 nm value at 780 nm is labelled an optimistic assumption. The round-1 ΔR sweep is kept for comparison.
- **Low NA and cameras (`analysis_v2.py`).**
  - Whole-plate optics are an object-side telecentric lens with coaxial illumination, NA 0.01–0.05. The object pixel comes from a 110 × 75 mm field on a 20 MP or 60 MP sensor.
  - Each camera has its own frame rate, full well and read noise. The blink rate is fps/4, because asynchronous lock-in needs 4 frames per cycle.
  - Diffraction blur is a Gaussian approximation of the Airy core.
  - Below NA 0.05, light that never scattered and reflects straight back is counted exactly. All other light is counted inside NA 0.05 and scaled as NA². A test checks this against direct counting.
- **Tilt and Lambertian SPOTs (`stages_v2.launch`).**
  - Incoming light at the chip is rebuilt from the Stage 1 histograms in both hemispheres.
  - A tilted mirror receives flux in proportion to max(0, −u·m)/|u_z| and reflects about its tilted normal m. A Lambertian SPOT re-emits in a cosine distribution about m.
  - With no tissue in the path, an analytic two-disk overlap model is used: the illumination fills the lens NA.
- **Spheres (`mc_sphere.py`).** A tissue sphere sits in medium (n = 1.33), lit by a collimated beam over its whole cross-section. Launch points are importance-sampled near the axis, with weights that remove the bias. Chips sit on the axis: at the camera-side surface, halfway, at the centre, and just inside the far surface. Background is the light leaving the sphere toward the camera, traced back to the chip plane near the axis. Readout from above and from below differ only in the outer-interface glare (2% vs 4%). Light reflected by the dish back into the sphere is ignored.
- **Decoder (`timesim.py`, `decoder.py`, `decode_sim.py`).** The decoder runs on synthetic data only.
  - *Frames:* `timesim.py` renders camera frames from the Monte Carlo photon budgets. It adds shot noise, read noise, light-source noise (white plus slow drift), pixel gain and offset non-uniformity, background texture and slow sample motion.
  - *Chip signal:* each chip blinks a 29-bit frame (Barker-13 preamble + 8-bit ID + 8-bit value). It uses differential phase-shift keying (DBPSK) of a square-wave carrier near fps/4, with a clock tolerance of ±0.3%.
  - *Decoding:* `decoder.py` does not see the truth. It removes common-mode light changes, locates chips with a lock-in power map, combines pixels with matched-filter weights, and searches carrier frequency and symbol timing. It then detects bits differentially and aligns the frame with the preamble.
  - *Read rule:* bisection finds the lowest modulation that decodes reliably. That fixes the rule "SNR per cycle × √(cycles per bit) ≥ threshold", which then gives the read time for every Monte Carlo configuration. `decoder_ber.csv` checks the rule against actual decoder runs.
- **Validation (`tests/test_v2.py`).**
  - sphere energy conservation
  - large sphere vs slab near the apex
  - low-NA scaling vs direct counting
  - the tilt-0 launch reproduces round 1
  - the surface-collection formulas
  - the decoder decodes strong chips and invents none on chip-free frames

  All round-1 checks still run first.

Round-2 simplifications:
- no re-entry of light into the sphere after it leaves, and no dish back-reflection
- chips are point sources at the chip centre for the transport, though their image has the full size
- tilt is about one axis only
- ideal telecentric optics
- the time-domain model covers only the listed noise sources (no biological or chip-electrical noise near the blink frequency)

## Round 3 (`run_v3.py`, outputs in `results/v3/`)

```
python run_v3.py             # all checks, Monte Carlo at 635/780/850 nm, tables, pitch assets (~10 min on 16 cores)
python run_v3.py --quick     # smoke test (~2 min)
python make_tables_v3.py     # rebuild tables_v3.md from results/v3/mc_v3.pkl
python assets_v3.py          # rebuild the chart and decoder GIF
```

Round 3 adds no new physics. It reruns the round-2 slab and sphere models with these additions:
- 635 nm, where the paper measured ΔR, becomes the headline wavelength.
- 850 nm is added.
- Two reflectance cases model a smaller swing at 780 and 850 nm: 50% and 25% of the 635 nm swing, at the same mean R.
- The centre chip is simulated at tilts 0–180° (`stages_v2.launch`).
- μa variants are added for the sensitivity ranking.

The read rule is the one calibrated in round 2 on the decoder (`results/v2/decoder_*.csv`).

Round-3 post-processing (`make_tables_v3.py`):
- **Random orientation.** cos(tilt) is uniform. log(SNR) is interpolated between the simulated tilts, and 20,000 sampled orientations give the read-time distribution.
- **Flatness.** h_max = side × tan(tilt limit). The tilt limits come from round 2's no-tissue tilt runs.
- **Landing fixes.** These are face and tilt bookkeeping only.
- **Sensitivity ranking.** Each uncertain input is varied over its range. Extra noise at the blink frequency is added to the SNR as an apparent modulation of RMS m_n/√T.

`assets_v3.py` uses matplotlib and Pillow, only for the two pitch assets. Both assets are labelled as simulation or synthetic data.
